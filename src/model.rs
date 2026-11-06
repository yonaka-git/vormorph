//! 時間ループと質量収支。
//!
//! ## 1ステップの流れ
//!
//! 1. [`HydroSolver::step`] で η と u_n を時刻 n+1 へ進める。時刻は開始時刻 + ステップ数 × Δt
//!    で計算し、Δt の足し算による誤差の累積を避ける。
//! 2. t > `time.morphology_start` かつ移動床なら、
//!    [`SedimentSolver::compute_bedload`]、[`SedimentSolver::step_suspended`]、
//!    [`SedimentSolver::update_bed`] の順に流砂と河床を更新する（流砂は未実装で、
//!    これらの関数は `todo!()` のまま）。
//! 3. 水の収支（[`MassBalance`]）を更新する。上流端・下流端の流量は、連続式と同じ θ 平均流量
//!    （[`HydroSolver::theta_flux`]）から求める。
//! 4. 出力時刻なら [`VtkWriter::write_step`] で書き出す。[`Model::run`] は最初に初期状態も書く。
//!
//! ## 質量収支
//!
//! 水の体積は、流入・流出の累積と、河床更新のときに湛水セルの水位を保ったことによる
//! 体積の変化で説明できなければならない。
//!
//! ```text
//! V(t) − V(0) = Σ (Q_in − Q_out) Δt + ΔV_bed
//! ```
//!
//! 半陰解法では、Newton 反復と PCG の収束判定値の範囲でこれが成り立つ。収支の誤差は
//! 検証ケース（閉鎖水域の振動など）で確かめる。

use crate::config::{Config, ConfigError, DownstreamBoundary, UpstreamBoundary};
use crate::hydro::{self, HydroError, HydroParams, HydroSolver, StepReport};
use crate::io::VtkWriter;
use crate::mesh::{EdgeKind, Mesh};
use crate::sediment::{SedimentParams, SedimentSolver};
use crate::state::State;
use std::error::Error;
use std::fmt;

/// 水の収支（累積）\[m³\]。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct MassBalance {
    /// 初期の水の体積。
    pub initial_volume: f64,
    /// 現在の水の体積。
    pub volume: f64,
    /// 上流端からの累積流入量。
    pub inflow: f64,
    /// 下流端からの累積流出量。
    pub outflow: f64,
    /// 河床更新のときに湛水セルの水位を保ったことによる、体積の累積変化。
    pub bed_exchange: f64,
}

impl MassBalance {
    /// 収支の誤差 V − V₀ − (流入 − 流出) − ΔV_bed \[m³\]。
    pub fn error(&self) -> f64 {
        self.volume - self.initial_volume - (self.inflow - self.outflow) - self.bed_exchange
    }
}

/// 土砂の収支（累積、空隙を含まない実体積）\[m³\]。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SedimentBalance {
    /// 初期からの河床の土砂量の変化（(1 − λ) Σ A_i Δz_{b,i}）。
    pub bed_change: f64,
    /// 現在の浮遊砂の量（Σ A_i (hC)_i）。
    pub suspended: f64,
    /// 初期の浮遊砂の量。
    pub initial_suspended: f64,
    /// 上流端からの累積給砂量。
    pub inflow: f64,
    /// 下流端からの累積流出量。
    pub outflow: f64,
}

impl SedimentBalance {
    /// 収支の誤差 \[m³\]。
    pub fn error(&self) -> f64 {
        self.bed_change + (self.suspended - self.initial_suspended) - (self.inflow - self.outflow)
    }
}

/// モデルの組み立てと実行の失敗。
#[derive(Debug)]
pub enum ModelError {
    /// 設定が正しくない。
    Config(ConfigError),
    /// 格子・状態・設定の組み合わせが正しくない（配列の長さ、境界条件と境界辺の対応など）。
    Setup(String),
    /// 流れの計算が失敗した。
    Hydro(HydroError),
    /// 入出力の失敗。
    Io(std::io::Error),
}

impl fmt::Display for ModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Config(e) => write!(f, "configuration error: {e}"),
            Self::Setup(message) => write!(f, "invalid model setup: {message}"),
            Self::Hydro(e) => write!(f, "hydrodynamic solver error: {e}"),
            Self::Io(e) => write!(f, "I/O error: {e}"),
        }
    }
}

impl Error for ModelError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Config(e) => Some(e),
            Self::Setup(_) => None,
            Self::Hydro(e) => Some(e),
            Self::Io(e) => Some(e),
        }
    }
}

impl From<ConfigError> for ModelError {
    fn from(e: ConfigError) -> Self {
        Self::Config(e)
    }
}

impl From<HydroError> for ModelError {
    fn from(e: HydroError) -> Self {
        Self::Hydro(e)
    }
}

impl From<std::io::Error> for ModelError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

/// 河床変動モデル。格子・状態・各ソルバー・出力をまとめ、時間ループを回す。
#[derive(Debug)]
pub struct Model {
    config: Config,
    mesh: Mesh,
    state: State,
    /// 初期の河床高（出力の dzb に使う）
    zb_initial: Vec<f64>,
    hydro: HydroSolver,
    sediment: Option<SedimentSolver>,
    writer: VtkWriter,
    water: MassBalance,
    sediment_balance: SedimentBalance,
    /// 開始時刻 \[s\]
    start_time: f64,
    /// 進めたステップ数
    steps: u64,
    /// 次に出力する時刻 \[s\]
    next_output: f64,
}

impl Model {
    /// 設定から格子と初期条件を作り、モデルを組み立てる。
    pub fn new(config: Config) -> Result<Self, ModelError> {
        let _ = config;
        todo!("格子の生成・初期条件の設定")
    }

    /// 格子と初期状態を外から与えて組み立てる（検証ケースで任意の地形を与えるため）。
    ///
    /// 設定の `[mesh]`、`[bed]`、`[[structures]]` は使わない。状態は次のように整える。
    ///
    /// - η < z_b のセルは η = z_b にする（乾いたセル。水の体積は変わらない）。
    /// - 壁の辺（上流端・下流端が壁の設定なら、その辺も）の u_n を 0 にする。
    /// - `u_cell` を `u_n` から Perot 法で再構成する。
    ///
    /// 配列の長さが格子と合わないとき、流量境界なのに上流端の辺がないとき、水位・水深境界
    /// なのに下流端の辺がないときは [`ModelError::Setup`] を返す。
    pub fn with_state(config: Config, mesh: Mesh, mut state: State) -> Result<Self, ModelError> {
        config.validate()?;
        check_setup(&config, &mesh, &state)?;

        for (eta, &zb) in state.eta.iter_mut().zip(&state.zb) {
            *eta = eta.max(zb);
        }
        let closed = |kind: EdgeKind| match kind {
            EdgeKind::Interior => false,
            EdgeKind::Wall => true,
            EdgeKind::Upstream => matches!(config.boundary.upstream, UpstreamBoundary::Wall),
            EdgeKind::Downstream => {
                matches!(config.boundary.downstream, DownstreamBoundary::Wall)
            }
        };
        for (u, &kind) in state.u_n.iter_mut().zip(&mesh.edge_kind) {
            if closed(kind) {
                *u = 0.0;
            }
        }
        hydro::reconstruct_cell_velocity(&mesh, &state.u_n, &mut state.u_cell);

        let hydro = HydroSolver::new(&mesh, HydroParams::from_config(&config));
        let sediment = SedimentParams::from_config(&config).map(|p| SedimentSolver::new(&mesh, p));
        let writer = VtkWriter::new(
            config.output.dir.clone(),
            config.title.clone(),
            config.hydro.dry_depth,
        );
        let volume = water_volume(&mesh, &state);
        let suspended = suspended_volume(&mesh, &state);
        Ok(Self {
            zb_initial: state.zb.clone(),
            start_time: state.time,
            next_output: state.time + config.output.interval,
            water: MassBalance {
                initial_volume: volume,
                volume,
                ..MassBalance::default()
            },
            sediment_balance: SedimentBalance {
                suspended,
                initial_suspended: suspended,
                ..SedimentBalance::default()
            },
            steps: 0,
            config,
            mesh,
            state,
            hydro,
            sediment,
            writer,
        })
    }

    /// 設定。
    pub fn config(&self) -> &Config {
        &self.config
    }

    /// 格子。
    pub fn mesh(&self) -> &Mesh {
        &self.mesh
    }

    /// 現在の状態。
    pub fn state(&self) -> &State {
        &self.state
    }

    /// 流れのソルバー（直前のステップの θ 平均流量などを読むため）。
    pub fn hydro(&self) -> &HydroSolver {
        &self.hydro
    }

    /// 進めたステップ数。
    pub fn steps(&self) -> u64 {
        self.steps
    }

    /// 水の収支。
    pub fn water_balance(&self) -> &MassBalance {
        &self.water
    }

    /// 土砂の収支。
    pub fn sediment_balance(&self) -> &SedimentBalance {
        &self.sediment_balance
    }

    /// 1ステップ進め、流れの計算の記録を返す。出力時刻を過ぎたら VTK を書く。
    pub fn step(&mut self) -> Result<StepReport, ModelError> {
        let dt = self.config.time.dt;
        let report = self
            .hydro
            .step(&self.mesh, &mut self.state, &self.config.boundary, dt)?;
        self.steps += 1;
        self.state.time = self.start_time + self.steps as f64 * dt;

        let (inflow, outflow) = boundary_discharge(&self.mesh, self.hydro.theta_flux());
        self.water.inflow += dt * inflow;
        self.water.outflow += dt * outflow;

        if let Some(sediment) = self.sediment.as_mut() {
            if self.state.time > self.config.time.morphology_start {
                // 河床変動の段階で、河床の更新による水の体積の変化（bed_exchange）と
                // 土砂の収支を加える。
                sediment.compute_bedload(&self.mesh, &self.state, dt);
                sediment.step_suspended(&self.mesh, &mut self.state, self.hydro.theta_flux(), dt);
                sediment.update_bed(&self.mesh, &mut self.state, dt);
            }
        }
        self.water.volume = water_volume(&self.mesh, &self.state);

        if self.state.time >= self.next_output - 1.0e-9 * dt {
            self.write_output()?;
            let interval = self.config.output.interval;
            while self.next_output <= self.state.time + 1.0e-9 * dt {
                self.next_output += interval;
            }
        }
        Ok(report)
    }

    /// 初期状態を書き出し（まだ何も書いていなければ）、設定の終了時刻まで進める。
    pub fn run(&mut self) -> Result<(), ModelError> {
        if self.writer.entries().is_empty() {
            self.write_output()?;
        }
        let (dt, end) = (self.config.time.dt, self.config.time.end);
        while self.state.time < end - 1.0e-9 * dt {
            self.step()?;
        }
        Ok(())
    }

    /// 現在の状態を `.vtu` に書き、`.pvd` を書き直す（出力時刻とは別に書きたいとき）。
    pub fn write_output(&mut self) -> Result<(), ModelError> {
        self.writer
            .write_step(&self.mesh, &self.state, &self.zb_initial)?;
        self.writer.write_collection()?;
        Ok(())
    }
}

/// 配列の長さと、境界条件・境界辺の対応を確かめる。
fn check_setup(config: &Config, mesh: &Mesh, state: &State) -> Result<(), ModelError> {
    let (n_cells, n_edges) = (mesh.n_cells(), mesh.n_edges());
    let mut problems = Vec::new();
    for (name, len, expected) in [
        ("eta", state.eta.len(), n_cells),
        ("zb", state.zb.len(), n_cells),
        ("z_ne", state.z_ne.len(), n_cells),
        ("hc", state.hc.len(), n_cells),
        ("u_cell", state.u_cell.len(), n_cells),
        ("u_n", state.u_n.len(), n_edges),
    ] {
        if len != expected {
            problems.push(format!(
                "state.{name} has {len} values, expected {expected}"
            ));
        }
    }
    let has = |kind: EdgeKind| mesh.edge_kind.contains(&kind);
    if matches!(config.boundary.upstream, UpstreamBoundary::Discharge { .. })
        && !has(EdgeKind::Upstream)
    {
        problems.push("a discharge boundary needs upstream edges in the mesh".to_owned());
    }
    if !matches!(config.boundary.downstream, DownstreamBoundary::Wall) && !has(EdgeKind::Downstream)
    {
        problems
            .push("a water level or depth boundary needs downstream edges in the mesh".to_owned());
    }
    if problems.is_empty() {
        Ok(())
    } else {
        Err(ModelError::Setup(problems.join("; ")))
    }
}

/// 水の体積 Σ A_i max(η_i − z_{b,i}, 0) \[m³\]。
fn water_volume(mesh: &Mesh, state: &State) -> f64 {
    (0..mesh.n_cells())
        .map(|c| mesh.area_cell[c] * state.depth(c))
        .sum()
}

/// 浮遊砂の量 Σ A_i (hC)_i \[m³\]。
fn suspended_volume(mesh: &Mesh, state: &State) -> f64 {
    mesh.area_cell
        .iter()
        .zip(&state.hc)
        .map(|(a, hc)| a * hc)
        .sum()
}

/// 上流端からの流入量と下流端からの流出量 \[m³/s\]（θ 平均流量から）。
fn boundary_discharge(mesh: &Mesh, theta_flux: &[f64]) -> (f64, f64) {
    let (mut inflow, mut outflow) = (0.0, 0.0);
    for (e, &kind) in mesh.edge_kind.iter().enumerate() {
        // 境界辺の法線は領域の外を向くので、L q は流出を正とする流量
        let flux = mesh.dv_edge[e] * theta_flux[e];
        match kind {
            EdgeKind::Upstream => inflow -= flux,
            EdgeKind::Downstream => outflow += flux,
            EdgeKind::Interior | EdgeKind::Wall => {}
        }
    }
    (inflow, outflow)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::ChannelLayout;
    use crate::mesh::channel::{self, ChannelSpec};

    fn config(upstream: &str, downstream: &str) -> Config {
        format!(
            r#"
title = "test"
[mesh]
kind = "file"
path = "unused"
[time]
dt = 0.5
end = 2.0
[hydro]
theta = 0.55
manning_n = 0.02
[boundary.upstream]
{upstream}
[boundary.downstream]
{downstream}
[output]
dir = "unused"
interval = 100.0
"#
        )
        .parse()
        .expect("valid config")
    }

    fn basin() -> Mesh {
        channel::build(&ChannelSpec {
            length: 4.0,
            width: 2.0,
            cell_size: 1.0,
            layout: ChannelLayout::Rectangular,
        })
        .expect("valid channel")
    }

    const WALL: &str = "kind = \"wall\"";

    #[test]
    fn state_is_normalized_and_time_advances_by_whole_steps() {
        let mesh = basin();
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        state.zb = (0..mesh.n_cells()).map(|c| 0.1 * c as f64).collect();
        state.eta = vec![0.3; mesh.n_cells()];
        state.u_n = vec![0.2; mesh.n_edges()];
        let mut model = Model::with_state(config(WALL, WALL), mesh, state).expect("valid setup");

        let (mesh, state) = (model.mesh(), model.state());
        assert!((0..mesh.n_cells()).all(|c| state.eta[c] >= state.zb[c]));
        for e in 0..mesh.n_edges() {
            let expected = if mesh.is_boundary_edge(e) { 0.0 } else { 0.2 };
            assert_eq!(state.u_n[e], expected);
        }
        let v0 = model.water_balance().initial_volume;
        assert!(v0 > 0.0);

        // run() は初期状態を書き出すので、ここでは step() で進める（出力間隔は終了時刻より長い）
        while model.state().time < model.config().time.end {
            model.step().expect("the step must succeed");
        }
        assert_eq!(model.steps(), 4);
        assert_eq!(model.state().time, 2.0);
        assert!(model.water_balance().error().abs() <= 1e-12 * v0);
    }

    #[test]
    fn mismatched_setups_are_rejected() {
        let mesh = basin();
        let state = State::new(mesh.n_cells() + 1, mesh.n_edges());
        let err = Model::with_state(config(WALL, WALL), mesh.clone(), state).unwrap_err();
        assert!(
            matches!(&err, ModelError::Setup(m) if m.contains("state.eta")),
            "{err}"
        );

        let mut mesh_without_ends = mesh;
        for kind in &mut mesh_without_ends.edge_kind {
            if *kind != EdgeKind::Interior {
                *kind = EdgeKind::Wall;
            }
        }
        let state = State::new(mesh_without_ends.n_cells(), mesh_without_ends.n_edges());
        let discharge = "kind = \"discharge\"\ndischarge = 1.0";
        let level = "kind = \"water_level\"\nlevel = 1.0";
        let err =
            Model::with_state(config(discharge, level), mesh_without_ends, state).unwrap_err();
        let message = err.to_string();
        assert!(message.contains("upstream edges"), "{message}");
        assert!(message.contains("downstream edges"), "{message}");
    }

    #[test]
    fn water_balance_error_accounts_for_all_terms() {
        let b = MassBalance {
            initial_volume: 100.0,
            volume: 112.0,
            inflow: 20.0,
            outflow: 5.0,
            bed_exchange: -3.0,
        };
        assert_eq!(b.error(), 0.0);
    }
}
