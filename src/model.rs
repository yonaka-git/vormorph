//! 時間ループと質量収支。
//!
//! ## モデルの組み立て
//!
//! - [`Model::new`]：設定から直線水路の格子（[`crate::mesh::channel`]）、初期河床、床止めを作る。
//!   河床高は母点 x_i で評価し（z_b = `z_downstream` + `slope` × (L − x_i)）、床止めごとに
//!   x_i ≤ `x_end` のセルを `height` だけかさ上げし、`x_start` ≤ x_i ≤ `x_end` のセルを
//!   非侵食層（z_ne = 初期の z_b）とする。母点が `x_start`・`x_end` にちょうど乗る場合に
//!   丸め誤差で分かれないよう、比較には格子幅の 10⁻⁹ 倍の許容差をつける。初期状態は
//!   乾いた河床（η = z_b、u = 0）で、水は境界から入る。
//! - [`Model::with_state`]：格子と初期状態を外から与える（検証ケースで任意の地形を使うため）。
//!
//! ## 1ステップの流れ
//!
//! 1. 時間刻み Δt を決める（[`TimeStep`]）。
//!    - 固定：`time.dt`。時刻は開始時刻 + ステップ数 × Δt で計算し、足し算による誤差の累積を
//!      避ける。
//!    - 可変：時刻 n の状態でセルの流出量の比 r = max_i Σ_流出 L q / V
//!      （[`HydroSolver::max_outflow_rate`]）を求め、Δt = min(C / r, `dt_max`) とする
//!      （流出がなければ `dt_max`）。次の区切りの時刻（出力時刻、`morphology_start`、終了時刻）
//!      までの残りが Δt 以下ならちょうどそこまで進め、時刻をその値にそろえる。残りが 2Δt 未満
//!      なら、半分ずつ2ステップで進める（極端に短いステップを作らないため）。ADR-0010。
//! 2. [`HydroSolver::step`] で η と u_n を時刻 n+1 へ進める。
//! 3. t > `time.morphology_start` かつ移動床なら、
//!    [`SedimentSolver::compute_bedload`]、[`SedimentSolver::step_suspended`]、
//!    [`SedimentSolver::update_bed`] の順に流砂と河床を更新する（流砂は未実装なので、
//!    `[sediment]` のある設定は組み立ての時点でエラーにする）。
//! 4. 水の収支（[`MassBalance`]）を更新する。上流端・下流端の流量は、連続式と同じ θ 平均流量
//!    （[`HydroSolver::theta_flux`]）から求める。
//! 5. 出力時刻なら [`VtkWriter::write_step`] で書き出す。[`Model::run`] は最初に初期状態も書く。
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

use crate::config::{
    BedConfig, Config, ConfigError, DownstreamBoundary, MeshConfig, StructureConfig, TimeStep,
    UpstreamBoundary,
};
use crate::hydro::{self, HydroError, HydroParams, HydroSolver, StepReport};
use crate::io::VtkWriter;
use crate::mesh::channel::{self, ChannelSpec};
use crate::mesh::{EdgeKind, Mesh, MeshError};
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
    /// 格子を作れない。
    Mesh(MeshError),
    /// 流れの計算が失敗した。
    Hydro(HydroError),
    /// 可変の時間刻みが正の有限の値にならない（時刻と、求めた値 \[s\]）。
    TimeStep {
        /// 時刻 \[s\]。
        time: f64,
        /// 求めた時間刻み \[s\]。
        dt: f64,
    },
    /// 入出力の失敗。
    Io(std::io::Error),
}

impl fmt::Display for ModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Config(e) => write!(f, "configuration error: {e}"),
            Self::Setup(message) => write!(f, "invalid model setup: {message}"),
            Self::Mesh(e) => write!(f, "cannot build the mesh: {e}"),
            Self::Hydro(e) => write!(f, "hydrodynamic solver error: {e}"),
            Self::TimeStep { time, dt } => write!(
                f,
                "the adaptive time step became {dt:e} s at t = {time} s (expected a positive finite value)"
            ),
            Self::Io(e) => write!(f, "I/O error: {e}"),
        }
    }
}

impl Error for ModelError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Config(e) => Some(e),
            Self::Setup(_) | Self::TimeStep { .. } => None,
            Self::Mesh(e) => Some(e),
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

impl From<MeshError> for ModelError {
    fn from(e: MeshError) -> Self {
        Self::Mesh(e)
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
    /// 時間刻みの決め方
    time_step: TimeStep,
    /// 開始時刻 \[s\]
    start_time: f64,
    /// 進めたステップ数
    steps: u64,
    /// 次に出力する時刻 \[s\]
    next_output: f64,
}

/// [`Model::run_with`] が出力のたびに渡す、直前の出力からの計算の記録。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Progress {
    /// 直前の出力から進めたステップ数。
    pub steps: u64,
    /// 時間刻みの最小値 \[s\]。
    pub dt_min: f64,
    /// 時間刻みの最大値 \[s\]。
    pub dt_max: f64,
    /// セルの流出量の比（クーラン数）の最大値 \[-\]（[`StepReport::outflow_courant`]）。
    pub max_outflow_courant: f64,
    /// 1ステップの Newton 法の反復回数の最大値。
    pub max_newton_iterations: usize,
    /// PCG の反復回数の合計。
    pub pcg_iterations: usize,
    /// 移流の重みを頭打ちにした辺の数の合計（[`StepReport::limited_edges`]）。
    pub limited_edges: usize,
}

impl Default for Progress {
    fn default() -> Self {
        Self {
            steps: 0,
            dt_min: f64::INFINITY,
            dt_max: 0.0,
            max_outflow_courant: 0.0,
            max_newton_iterations: 0,
            pcg_iterations: 0,
            limited_edges: 0,
        }
    }
}

impl Progress {
    /// 1ステップの記録を加える。
    pub fn add(&mut self, report: &StepReport) {
        self.steps += 1;
        self.dt_min = self.dt_min.min(report.dt);
        self.dt_max = self.dt_max.max(report.dt);
        self.max_outflow_courant = self.max_outflow_courant.max(report.outflow_courant);
        self.max_newton_iterations = self.max_newton_iterations.max(report.newton_iterations);
        self.pcg_iterations += report.pcg_iterations;
        self.limited_edges += report.limited_edges;
    }
}

impl Model {
    /// 設定から格子と初期条件を作り、モデルを組み立てる。
    ///
    /// `mesh.kind = "channel"` のとき、直線水路の格子を作り、`[bed]` と `[[structures]]` から
    /// 初期河床と非侵食層を決める（モジュールの説明の「モデルの組み立て」）。初期状態は乾いた
    /// 河床（η = z_b、u = 0）。
    ///
    /// `mesh.kind = "file"` は、格子ファイルにセルの河床高がないので、まだ使えない
    /// （[`ModelError::Setup`]。[`Model::with_state`] で状態を与える）。
    pub fn new(config: Config) -> Result<Self, ModelError> {
        config.validate()?;
        let (mesh, state) = match &config.mesh {
            MeshConfig::Channel {
                length,
                width,
                cell_size,
                layout,
            } => {
                let mesh = channel::build(&ChannelSpec {
                    length: *length,
                    width: *width,
                    cell_size: *cell_size,
                    layout: *layout,
                })?;
                let bed = config.bed.as_ref().ok_or_else(|| {
                    ModelError::Setup("[bed] is required when mesh.kind = \"channel\"".to_owned())
                })?;
                let state = channel_state(&mesh, *length, *cell_size, bed, &config.structures);
                (mesh, state)
            }
            MeshConfig::File { path } => {
                return Err(ModelError::Setup(format!(
                    "mesh.kind = \"file\" ({}) is not supported by Model::new yet: mesh files have no bed elevations (build the state yourself and use Model::with_state)",
                    path.display()
                )));
            }
        };
        Self::with_state(config, mesh, state)
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
    /// なのに下流端の辺がないとき、まだ実装していない流砂（`[sediment]`）を含むときは
    /// [`ModelError::Setup`] を返す。
    pub fn with_state(config: Config, mesh: Mesh, mut state: State) -> Result<Self, ModelError> {
        config.validate()?;
        check_setup(&config, &mesh, &state)?;
        let time_step = config.time.step_control().ok_or_else(|| {
            ModelError::Setup("time: give either dt or courant and dt_max".to_owned())
        })?;

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
            time_step,
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

    /// 時間刻みの決め方。
    pub fn time_step(&self) -> TimeStep {
        self.time_step
    }

    /// 計算が設定の終了時刻に達したか。
    pub fn finished(&self) -> bool {
        self.state.time >= self.config.time.end - 1.0e-9 * self.time_step.max_dt()
    }

    /// 1ステップ進め、流れの計算の記録を返す。出力時刻を過ぎたら VTK を書く。
    pub fn step(&mut self) -> Result<StepReport, ModelError> {
        let (dt, landing) = match self.time_step {
            TimeStep::Fixed { dt } => (dt, None),
            TimeStep::Adaptive { courant, dt_max } => {
                let rate =
                    self.hydro
                        .max_outflow_rate(&self.mesh, &self.state, &self.config.boundary);
                let dt = if rate > 0.0 {
                    (courant / rate).min(dt_max)
                } else {
                    dt_max
                };
                if !(dt.is_finite() && dt > 0.0) {
                    return Err(ModelError::TimeStep {
                        time: self.state.time,
                        dt,
                    });
                }
                self.fit_to_next_event(dt)
            }
        };
        let report = self
            .hydro
            .step(&self.mesh, &mut self.state, &self.config.boundary, dt)?;
        self.steps += 1;
        self.state.time = match (self.time_step, landing) {
            (TimeStep::Fixed { dt }, _) => self.start_time + self.steps as f64 * dt,
            (TimeStep::Adaptive { .. }, Some(event)) => event,
            (TimeStep::Adaptive { .. }, None) => self.state.time,
        };

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

    /// 可変の時間刻み `dt` を、次の区切りの時刻（出力時刻、河床変動の開始時刻、終了時刻）に
    /// 合わせる。区切りまでの残りが `dt` 以下ならちょうどそこまで進め（その時刻を返す）、
    /// 2`dt` 未満なら残りの半分にする。
    fn fit_to_next_event(&self, dt: f64) -> (f64, Option<f64>) {
        let time = self.state.time;
        let time_config = &self.config.time;
        let next_event = [
            self.next_output,
            time_config.morphology_start,
            time_config.end,
        ]
        .into_iter()
        .filter(|&t| t > time)
        .fold(f64::INFINITY, f64::min);
        let remaining = next_event - time;
        if remaining <= dt {
            (remaining, Some(next_event))
        } else if remaining < 2.0 * dt {
            (0.5 * remaining, None)
        } else {
            (dt, None)
        }
    }

    /// 初期状態を書き出し（まだ何も書いていなければ）、設定の終了時刻まで進める。
    pub fn run(&mut self) -> Result<(), ModelError> {
        self.run_with(|_, _| {})
    }

    /// [`Model::run`] と同じく終了時刻まで進め、出力を書くたびに `on_output` を呼ぶ。
    ///
    /// `on_output` には、モデル（書き出した時刻の状態）と、直前の出力からの計算の記録
    /// （[`Progress`]）を渡す。初期状態を書いたときは呼ばない。
    pub fn run_with(
        &mut self,
        mut on_output: impl FnMut(&Self, &Progress),
    ) -> Result<(), ModelError> {
        if self.writer.entries().is_empty() {
            self.write_output()?;
        }
        let mut progress = Progress::default();
        while !self.finished() {
            let written = self.writer.entries().len();
            let report = self.step()?;
            progress.add(&report);
            if self.writer.entries().len() > written {
                on_output(self, &progress);
                progress = Progress::default();
            }
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
    if config.sediment.is_some() {
        problems.push(
            "sediment transport is not implemented yet: remove [sediment] to run a fixed-bed simulation"
                .to_owned(),
        );
    }
    if problems.is_empty() {
        Ok(())
    } else {
        Err(ModelError::Setup(problems.join("; ")))
    }
}

/// 直線水路の初期状態：河床高と非侵食層を母点で評価し、乾いた河床（η = z_b、u = 0）にする。
///
/// z_b = `z_downstream` + `slope` × (L − x_i)。床止めごとに、x_i ≤ `x_end` のセルを `height`
/// だけかさ上げし、`x_start` ≤ x_i ≤ `x_end` のセルの非侵食層の高さを、すべての床止めを
/// 加えた後の z_b にする（天端の上に土砂がない状態）。比較には `cell_size` の 10⁻⁹ 倍の
/// 許容差をつける。
fn channel_state(
    mesh: &Mesh,
    length: f64,
    cell_size: f64,
    bed: &BedConfig,
    structures: &[StructureConfig],
) -> State {
    let tol = 1.0e-9 * cell_size;
    let mut state = State::new(mesh.n_cells(), mesh.n_edges());
    for (c, &x) in mesh.x_cell.iter().enumerate() {
        let mut zb = bed.z_downstream + bed.slope * (length - x);
        for structure in structures {
            match structure {
                StructureConfig::Sill { x_end, height, .. } => {
                    if x <= x_end + tol {
                        zb += height;
                    }
                }
            }
        }
        state.zb[c] = zb;
        state.eta[c] = zb;
        let fixed = structures.iter().any(|structure| match structure {
            StructureConfig::Sill { x_start, x_end, .. } => x_start - tol <= x && x <= x_end + tol,
        });
        if fixed {
            state.z_ne[c] = zb;
        }
    }
    state
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

    /// 直線水路の設定（正方形配置、格子幅 0.5 m で母点は x = 0.25, 0.75, …, 9.75 m）。
    fn channel_config(name: &str, structures: &str, time: &str) -> Config {
        let dir = std::env::temp_dir().join(format!("vormorph-{name}-{}", std::process::id()));
        format!(
            r#"
title = "{name}"
[mesh]
kind = "channel"
length = 10.0
width = 2.0
cell_size = 0.5
layout = "rectangular"
[bed]
slope = 0.01
z_downstream = 0.5
{structures}
[time]
{time}
end = 3.0
morphology_start = 1.5
[hydro]
theta = 0.55
manning_n = 0.03
[boundary.upstream]
kind = "discharge"
discharge = 0.5
[boundary.downstream]
kind = "depth"
depth = 0.2
[output]
dir = {dir:?}
interval = 1.0
"#
        )
        .parse()
        .expect("valid channel config")
    }

    const SILL: &str = r#"
[[structures]]
kind = "sill"
x_start = 4.0
x_end = 5.0
height = 1.0
"#;

    #[test]
    fn new_evaluates_the_bed_and_the_sill_at_the_generators() {
        let model =
            Model::new(channel_config("bed", SILL, "dt = 0.1")).expect("the channel must build");
        let (mesh, state) = (model.mesh(), model.state());
        assert_eq!(mesh.n_cells(), 80);
        let mut fixed = 0;
        for c in 0..mesh.n_cells() {
            let x = mesh.x_cell[c];
            let raised = if x <= 5.0 { 1.0 } else { 0.0 };
            let zb = 0.5 + 0.01 * (10.0 - x) + raised;
            assert!((state.zb[c] - zb).abs() < 1e-14, "zb at x = {x}");
            assert_eq!(state.eta[c], state.zb[c], "the initial bed must be dry");
            if (4.0..=5.0).contains(&x) {
                assert_eq!(state.z_ne[c], state.zb[c], "non-erodible at x = {x}");
                fixed += 1;
            } else {
                assert_eq!(state.z_ne[c], f64::NEG_INFINITY, "erodible at x = {x}");
            }
        }
        // x = 4.25, 4.75 m の2列（各 4 セル）
        assert_eq!(fixed, 8);
        assert!(state.u_n.iter().all(|&u| u == 0.0));
        assert_eq!(model.water_balance().initial_volume, 0.0);
    }

    /// 母点が床止めの端にちょうど乗る場合を、丸め誤差で取りこぼさない。
    #[test]
    fn generators_on_the_sill_ends_are_included_despite_rounding() {
        let mut config = channel_config("sill-ends", SILL, "dt = 0.1");
        let StructureConfig::Sill { x_start, x_end, .. } = &mut config.structures[0];
        *x_start = 4.25 + 1e-12;
        *x_end = 5.25 - 1e-12;
        let model = Model::new(config).expect("the channel must build");
        let (mesh, state) = (model.mesh(), model.state());
        for c in 0..mesh.n_cells() {
            let x = mesh.x_cell[c];
            let in_band = (4.2..5.3).contains(&x);
            assert_eq!(state.z_ne[c].is_finite(), in_band, "x = {x}");
            let raised = state.zb[c] - (0.5 + 0.01 * (10.0 - x));
            assert!((raised - if x < 5.3 { 1.0 } else { 0.0 }).abs() < 1e-14);
        }
    }

    #[test]
    fn new_rejects_file_meshes_and_sediment_for_now() {
        let mut config = channel_config("reject", "", "dt = 0.1");
        config.mesh = MeshConfig::File {
            path: "mesh.toml".into(),
        };
        let err = Model::new(config).unwrap_err();
        assert!(
            matches!(&err, ModelError::Setup(m) if m.contains("not supported")),
            "{err}"
        );

        let mut config = channel_config("reject", "", "dt = 0.1");
        config.sediment =
            crate::config::Config::load(concat!(env!("CARGO_MANIFEST_DIR"), "/examples/sill.toml"))
                .expect("examples/sill.toml must be valid")
                .sediment;
        let err = Model::new(config).unwrap_err();
        assert!(
            matches!(&err, ModelError::Setup(m) if m.contains("sediment transport is not implemented")),
            "{err}"
        );
    }

    /// 可変の時間刻みは上限とクーラン数を守り、出力時刻・河床変動の開始時刻・終了時刻に
    /// ちょうど止まる。
    #[test]
    fn adaptive_steps_respect_the_limits_and_land_on_events() {
        let (courant, dt_max) = (0.5, 0.2);
        let config = channel_config(
            "adaptive",
            SILL,
            &format!("courant = {courant:?}\ndt_max = {dt_max:?}"),
        );
        let mut model = Model::new(config).expect("the channel must build");
        assert_eq!(model.time_step(), TimeStep::Adaptive { courant, dt_max });
        let mut times = Vec::new();
        let mut limited_by_courant = 0;
        while !model.finished() {
            let report = model.step().expect("the step must succeed");
            assert!(report.dt > 0.0 && report.dt <= dt_max, "dt = {}", report.dt);
            assert!(
                report.outflow_courant <= courant * (1.0 + 1e-12),
                "Courant {} at t = {}",
                report.outflow_courant,
                model.state().time
            );
            if report.outflow_courant > 0.99 * courant {
                limited_by_courant += 1;
            }
            times.push(model.state().time);
        }
        assert!(
            limited_by_courant > 0,
            "the Courant condition must be active"
        );
        for event in [1.0, 1.5, 2.0, 3.0] {
            assert!(times.contains(&event), "t = {event} must be hit exactly");
        }
        assert_eq!(model.state().time, 3.0);
        assert!(times.windows(2).all(|w| w[0] < w[1]));
        let water = model.water_balance();
        assert!(water.inflow > 0.0);
        assert!(water.error().abs() <= 1e-10 * water.inflow, "{water:?}");
    }

    /// `run_with` は出力のたびに呼ばれ、直前の出力からのステップ数を渡す。
    #[test]
    fn run_with_reports_progress_at_each_output() {
        let mut model =
            Model::new(channel_config("progress", "", "dt = 0.1")).expect("the channel must build");
        let mut outputs = Vec::new();
        model
            .run_with(|model, progress| {
                outputs.push((model.state().time, progress.steps, progress.dt_min));
            })
            .expect("the run must succeed");
        assert_eq!(outputs.len(), 3);
        for (k, &(time, steps, dt)) in outputs.iter().enumerate() {
            assert!((time - (k + 1) as f64).abs() < 1e-12);
            assert_eq!(steps, 10);
            assert_eq!(dt, 0.1);
        }
        // 初期状態を含めて 4 回書く
        assert_eq!(model.writer.entries().len(), 4);
    }
}
