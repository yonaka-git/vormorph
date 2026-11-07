//! 流砂と河床変動：掃流砂、浮遊砂、Exner 式（均一粒径）。
//!
//! 離散式の詳細は `docs/design/overview.md` の第4.6〜4.8節と第4.10節、判断の理由は ADR-0011
//! （掃流砂）、ADR-0012（河床の更新）、ADR-0019（二次流）。浮遊砂
//! （[`SedimentSolver::step_suspended`]）は未実装。
//!
//! ## 掃流砂（辺の法線成分）
//!
//! 掃流砂量は、辺 j の法線成分 q_{b,j} を直接評価する（[`SedimentSolver::compute_bedload`]）。
//!
//! ```text
//! q_{b,j} = √(s g d³) Φ(τ*_j) [ u_j / |u|_j − γ_j ∂z_b/∂n ]
//! τ*_j    = g n² |u|_j² / (H_j^{1/3} s g d),   γ_j = √(τ*_c / (μ_s μ_k τ*_j))
//! ```
//!
//! - 流れの値は、流れを進めた後（時刻 n+1）の u_j、η を使う。辺の水深 H_j は流れの計算と
//!   同じ規則（[`hydro::edge_depths`]：風上側の水位 − 両側の高い方の河床高）、|u|_j は摩擦と
//!   同じ値（[`hydro::edge_speed`]）。したがって掃流力は、運動量式の摩擦 ρ g n² |u|² / H^{1/3} と
//!   一致する。H_j が `dry_depth` 未満、または τ* ≤ τ*_c の辺では q_{b,j} = 0。
//! - γ_j は長谷川の式による河床勾配の補正係数。内部の辺では、母点を結ぶ線分が辺に直交するので、
//!   法線方向の勾配は2点差分 (z_{b,2(j)} − z_{b,1(j)}) / δ_j になる（全方向にかける）。
//!   周期境界では河床高を縦断勾配 I を除いた値で持つので、勾配は
//!   (z_{b,2(j)} − z_{b,1(j)}) / δ_j − I n_{x,j} になる（[`SedimentParams::background_slope`]、
//!   ADR-0014）。安息角補正の高低差にも、除いた勾配の分 I n_{x,j} δ_j を戻す。
//! - Φ は [`formulas`] の MPM 式または芦田・道上式（設定で切り替える）で、有効掃流力
//!   τ*_e = μ τ*_j（μ は砂漣係数 `ripple_factor`、既定 1。[`formulas::bedload_intensity`]、
//!   ADR-0021）を使う。MPM 式は τ*_e で限界と比べ、芦田・道上式は τ*_e を先頭の係数だけに使う。
//!   勾配補正の γ_j は全掃流力 τ*_j のまま。
//! - 二次流（`[sediment.secondary_flow]`）を考えるときは、向きの項 u_j / |u|_j を、河床付近の
//!   流速の向き（水深平均流速から内岸側へ tan δ = N* I / |u| だけ回した向き）の法線成分に置き換える
//!   （[`secondary_flow::bedload_direction`]）。二次流の強さ I は、流れを進めた後に
//!   [`SedimentSolver::step_secondary_flow`] で求める（[`secondary_flow`]）。
//! - 境界辺：壁では 0。上流端（平衡給砂）と下流端（自由流出）では、内側のセルの水深と辺の流速から
//!   同じ式で求め、勾配 ∂z_b/∂n には内側のセルの河床勾配（隣のセルとの最小二乗。[`bed_gradient`]）
//!   の法線成分を使う。上流端は流入だけ、下流端は流出だけを許す。勾配を内側のセルから
//!   外挿するので、一様な勾配の等流では境界辺の流砂量も内部の辺と一致し、平衡河床が保たれる
//!   （検証ケース V12）。
//! - 上流端の給砂：平衡給砂（`equilibrium`）は、各ステップの上流端の流れから求めた上の量。
//!   一定の給砂量（`constant`）は、河床変動を始めた最初のステップの平衡給砂量の合計
//!   Σ L_j q_{b,j} を保ち、各ステップの平衡給砂量に比例して辺に配る（ADR-0013）。給砂しない
//!   （`none`）なら 0。
//!
//! **非侵食層（床止め）**：セル i から1ステップで流出する土砂量
//! O_i = M Δt Σ_{流出} s L q_b と、時刻 n に非侵食層より上にある土砂量
//! S_i = (1 − λ) A_i max(z_{b,i} − z_{ne,i}, 0) を比べ、O_i > S_i なら、セル i から流出する辺の
//! 流砂量をすべて S_i / O_i 倍に縮める。同じステップの流入は数えないので、反復なしで
//! z_b ≥ z_ne が保証される（床止めの上に来た土砂は1ステップ遅れで通過する）。
//!
//! ## 河床（Exner 式）
//!
//! ```text
//! (1 − λ) A_i Δz_{b,i} = −M Δt Σ_j s_ij L_j q_{b,j}
//! ```
//!
//! M は MORFAC（流れの時間 Δt を M 倍した地形の時間で河床を進める）。続いて
//! [`SedimentSolver::update_bed`] が次の処理を行う。
//!
//! 1. **安息角補正**：勾配が tan φ を超える内部の辺で、その辺だけならちょうど tan φ になる
//!    体積の半分（緩和係数 0.5）を高い側から低い側へ移す。全辺を同時に評価する Jacobi 型で、
//!    辺のループで移す量を辺の配列に書き、セルのループで集める（両側の値が丸めまで一致し、
//!    体積が厳密に保存される）。高い側のセルから出る量の合計は、非侵食層より上の量で頭打ちに
//!    する。1回の反復で動いた厚さの最大値が粒径の 10⁻³ 倍未満になるか、100 回で止める。
//! 2. **水位の調整**：水深が `dry_depth` 以上のセルでは水位 η を保つ（河床が水面を超えたら
//!    η = z_b）。それ未満のセルでは水深を保つ（ほぼ乾いたセルで水を作らない）。水の体積の
//!    変化は [`BedUpdateReport::water_volume_change`] で返し、モデルが水の収支に記録する。
//!
//! ## 浮遊砂
//!
//! hC の移流拡散式を、連続式と同じ θ 平均流量 q_j^θ（[`crate::hydro::HydroSolver::theta_flux`]）
//! で解く。沈降 w_f c_b は陰的に扱い、水深が `min_depth` 未満のセルでは浮上量を 0 にする。
//! 未実装（`[sediment.suspended]` のある設定は、モデルの組み立ての時点でエラーにする）。
//!
//! ## 参考文献
//!
//! - Meyer-Peter & Müller (1948), 芦田・道上 (1972): 掃流砂量式
//! - Struiksma (1983): 砂漣係数 μ = (C/C_90)^{3/2} を掛けた MPM 式
//! - 長谷川 (1984)（式の形は長谷川, 1983 と Nays2DH マニュアルで確認）: 河床勾配の補正
//! - Itakura & Kishi (1980), Garcia & Parker (1991): 浮遊砂の浮上量
//! - Rubey (1933): 沈降速度

pub mod formulas;
pub mod secondary_flow;

use crate::config::{Config, DownstreamBoundary, SedimentConfig, SedimentSupply, UpstreamBoundary};
use crate::hydro;
use crate::mesh::{EdgeKind, Mesh, NONE, gradient};
use crate::state::State;
use secondary_flow::{SecondaryFlow, SecondaryFlowParams};

/// 安息角補正の緩和係数（その辺だけなら勾配がちょうど tan φ になる体積に掛ける）。
const AVALANCHE_RELAXATION: f64 = 0.5;
/// 安息角補正の反復の上限。
const AVALANCHE_MAX_ITERATIONS: usize = 100;
/// 安息角補正の収束判定：1回の反復で動いた厚さの最大値 / 粒径。
const AVALANCHE_TOLERANCE_PER_GRAIN: f64 = 1.0e-3;

/// 流砂の計算に使うパラメータ。物理定数は設定から受け取る。
#[derive(Debug, Clone, PartialEq)]
pub struct SedimentParams {
    /// 重力加速度 g \[m/s²\]。
    pub gravity: f64,
    /// 水の密度 ρ \[kg/m³\]。
    pub water_density: f64,
    /// 水の動粘性係数 ν \[m²/s\]。
    pub kinematic_viscosity: f64,
    /// カルマン定数 κ（二次流の適応長に使う）。
    pub von_karman: f64,
    /// マニングの粗度係数 n（掃流力の評価に使う）。
    pub manning_n: f64,
    /// wet/dry を判定する水深の閾値 \[m\]。
    pub dry_depth: f64,
    /// 流砂の設定。
    pub sediment: SedimentConfig,
    /// 上流端の給砂の条件（上流端が壁なら `None`）。
    pub supply: Option<SedimentSupply>,
    /// 下流端から土砂が流出できるか（下流端が壁なら `false`）。
    pub downstream_open: bool,
    /// 周期境界で河床高から除いた縦断勾配 I \[-\]（x の正の向きに下がる向きを正）。
    /// 周期境界でなければ 0（ADR-0014）。
    pub background_slope: f64,
}

impl SedimentParams {
    /// 設定から作る。`[sediment]` が省略されている（固定床）なら `None`。
    pub fn from_config(config: &Config) -> Option<Self> {
        let sediment = config.sediment.clone()?;
        Some(Self {
            gravity: config.physics.gravity,
            water_density: config.physics.water_density,
            kinematic_viscosity: config.physics.kinematic_viscosity,
            von_karman: config.physics.von_karman,
            manning_n: config.hydro.manning_n,
            dry_depth: config.hydro.dry_depth,
            sediment,
            supply: match &config.boundary.upstream {
                UpstreamBoundary::Discharge {
                    sediment_supply, ..
                } => Some(*sediment_supply),
                UpstreamBoundary::Wall | UpstreamBoundary::Periodic { .. } => None,
            },
            downstream_open: matches!(
                config.boundary.downstream,
                DownstreamBoundary::WaterLevel { .. } | DownstreamBoundary::Depth { .. }
            ),
            background_slope: hydro::background_slope(config),
        })
    }

    /// 水中比重 s = ρ_s / ρ − 1 \[-\]。
    pub fn relative_density(&self) -> f64 {
        self.sediment.density / self.water_density - 1.0
    }

    /// 無次元掃流砂量を掃流砂量に直す係数 √(s g d³) \[m²/s\]。
    pub fn bedload_scale(&self) -> f64 {
        let d = self.sediment.grain_size;
        (self.relative_density() * self.gravity * d * d * d).sqrt()
    }

    /// 流速の大きさ `speed` \[m/s\]、水深 `depth` \[m\] での無次元掃流力
    /// τ* = u*² / (s g d)、u*² = g n² |u|² / h^{1/3}。
    pub fn shields(&self, speed: f64, depth: f64) -> f64 {
        let g = self.gravity;
        let u_star2 = g * self.manning_n * self.manning_n * speed * speed / depth.cbrt();
        u_star2 / (self.relative_density() * g * self.sediment.grain_size)
    }
}

/// 掃流砂量の計算の記録（[`SedimentSolver::compute_bedload`]）。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct BedloadReport {
    /// 非侵食層で流出量を頭打ちにしたセルの数。
    pub capped_cells: usize,
    /// 上流端からの給砂量（地形の時間 M Δt での体積、空隙を含まない）\[m³\]。
    pub inflow: f64,
    /// 下流端からの流出量（同上）\[m³\]。
    pub outflow: f64,
}

/// 河床更新の記録（[`SedimentSolver::update_bed`]）。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct BedUpdateReport {
    /// 河床高の変化量の最大値 |Δz_b|（Exner 式と安息角補正の合計）\[m\]。
    pub max_change: f64,
    /// 安息角補正で土砂を動かした反復の回数。
    pub avalanche_iterations: usize,
    /// 安息角補正が反復の上限までに収束しなかったか。
    pub avalanche_unconverged: bool,
    /// 河床の土砂の体積の変化 (1 − λ) Σ A_i Δz_{b,i}（空隙を含まない）\[m³\]。
    pub bed_volume_change: f64,
    /// 水位の調整による水の体積の変化 Σ A_i (h_i^新 − h_i^旧) \[m³\]。
    pub water_volume_change: f64,
}

/// 流砂と河床変動のソルバー。
#[derive(Debug, Clone)]
pub struct SedimentSolver {
    params: SedimentParams,
    /// 辺の掃流砂量 q_{b,j}（法線方向、単位幅あたり、頭打ちの後）\[m²/s\]
    bedload_flux: Vec<f64>,
    /// 辺の水深 H_j \[m\]（作業配列）
    edge_depth: Vec<f64>,
    /// 境界辺の外側の水位（すべて −∞。[`hydro::edge_depths`] に渡す）
    outside_level: Vec<f64>,
    /// セルの流出量を縮める倍率（非侵食層の頭打ち、安息角補正）
    cap_factor: Vec<f64>,
    /// 安息角補正で辺を通して移す河床の体積（法線の向きを正）\[m³\]
    avalanche_flux: Vec<f64>,
    /// 更新前の河床高 \[m\]
    zb_old: Vec<f64>,
    /// 二次流の強さのソルバー（`[sediment.secondary_flow]` がなければ `None`）
    secondary_flow: Option<SecondaryFlow>,
    /// 一定の給砂量（`constant`）で保つ、上流端からの給砂量の合計 Σ L_j q_{b,j}（流入を負）
    /// \[m³/s\]（最初の [`SedimentSolver::compute_bedload`] で決める。それまでは `None`）
    held_supply: Option<f64>,
    /// セルの浮上量 E_i \[m/s\]
    #[expect(dead_code, reason = "浮遊砂（V08）の実装で使う")]
    entrainment: Vec<f64>,
    /// セルの沈降量 w_f c_{b,i} \[m/s\]
    #[expect(dead_code, reason = "浮遊砂（V08）の実装で使う")]
    deposition: Vec<f64>,
}

impl SedimentSolver {
    /// 格子の大きさに合わせて作業配列を確保する。
    pub fn new(mesh: &Mesh, params: SedimentParams) -> Self {
        let (n_cells, n_edges) = (mesh.n_cells(), mesh.n_edges());
        let secondary_flow = params.sediment.secondary_flow.as_ref().map(|config| {
            SecondaryFlow::new(
                mesh,
                SecondaryFlowParams {
                    config: config.clone(),
                    gravity: params.gravity,
                    manning_n: params.manning_n,
                    von_karman: params.von_karman,
                    dry_depth: params.dry_depth,
                    background_slope: params.background_slope,
                },
            )
        });
        Self {
            params,
            secondary_flow,
            bedload_flux: vec![0.0; n_edges],
            edge_depth: vec![0.0; n_edges],
            outside_level: vec![f64::NEG_INFINITY; n_edges],
            cap_factor: vec![1.0; n_cells],
            avalanche_flux: vec![0.0; n_edges],
            zb_old: vec![0.0; n_cells],
            held_supply: None,
            entrainment: vec![0.0; n_cells],
            deposition: vec![0.0; n_cells],
        }
    }

    /// パラメータ。
    pub fn params(&self) -> &SedimentParams {
        &self.params
    }

    /// 二次流の強さのソルバー（`[sediment.secondary_flow]` がなければ `None`）。
    pub fn secondary_flow(&self) -> Option<&SecondaryFlow> {
        self.secondary_flow.as_ref()
    }

    /// 直前に計算した辺の掃流砂量 q_{b,j}（頭打ちの後）\[m²/s\]。
    pub fn bedload_flux(&self) -> &[f64] {
        &self.bedload_flux
    }

    /// 辺の掃流砂量 q_{b,j} を計算する（勾配補正、境界の給砂と流出、非侵食層の頭打ちを含む）。
    ///
    /// `state` は流れを進めた後の状態（時刻 n+1 の η、u_n、u_cell）、`dt` は流れの時間刻み
    /// （頭打ちの判定に M Δt の流出量を使う）。
    pub fn compute_bedload(&mut self, mesh: &Mesh, state: &State, dt: f64) -> BedloadReport {
        let Self {
            params,
            bedload_flux,
            edge_depth,
            outside_level,
            cap_factor,
            held_supply,
            ..
        } = self;
        hydro::edge_depths(
            mesh,
            state,
            outside_level,
            params.dry_depth,
            params.background_slope,
            edge_depth,
        );
        let scale = params.bedload_scale();
        for (e, q) in bedload_flux.iter_mut().enumerate() {
            *q = edge_capacity(params, mesh, state, edge_depth, scale, e);
        }
        // 一定の給砂量：最初のステップの合計を保ち、今の平衡給砂量に比例して辺に配る。辺ごとに
        // 値を保つと、河床が上がって流量の配分が減った辺にも同じ量を入れ続けるので、上流端の
        // セルの河床の偏りが増え続ける（ADR-0013）
        if params.supply == Some(SedimentSupply::Constant) {
            let inlet = || (0..mesh.n_edges()).filter(|&e| mesh.edge_kind[e] == EdgeKind::Upstream);
            let total: f64 = inlet().map(|e| mesh.dv_edge[e] * bedload_flux[e]).sum();
            let held = *held_supply.get_or_insert(total);
            if total < 0.0 {
                let ratio = held / total;
                for e in inlet() {
                    bedload_flux[e] *= ratio;
                }
            } else if held < 0.0 {
                // 平衡給砂量が 0（流れが止まった）なら、辺の長さに比例して配る
                let length: f64 = inlet().map(|e| mesh.dv_edge[e]).sum();
                for e in inlet() {
                    bedload_flux[e] = held / length;
                }
            }
        }

        // 非侵食層の頭打ち：セルごとに倍率を求め、辺は流出側のセルの倍率を掛ける
        let factor = params.sediment.morfac * dt;
        let solid = 1.0 - params.sediment.porosity;
        let mut report = BedloadReport::default();
        for (c, alpha) in cap_factor.iter_mut().enumerate() {
            let outflow = factor * cell_outflow(mesh, c, bedload_flux, true);
            let available = solid * mesh.area_cell[c] * (state.zb[c] - state.z_ne[c]).max(0.0);
            *alpha = if outflow > available {
                report.capped_cells += 1;
                available / outflow
            } else {
                1.0
            };
        }
        for (e, q) in bedload_flux.iter_mut().enumerate() {
            let [c1, c2] = mesh.cells_on_edge[e];
            if *q > 0.0 {
                *q *= cap_factor[c1];
            } else if *q < 0.0 && c2 != NONE {
                *q *= cap_factor[c2];
            }
        }

        for (e, &q) in bedload_flux.iter().enumerate() {
            let flux = factor * mesh.dv_edge[e] * q;
            match mesh.edge_kind[e] {
                EdgeKind::Upstream => report.inflow -= flux,
                EdgeKind::Downstream => report.outflow += flux,
                EdgeKind::Interior | EdgeKind::Wall => {}
            }
        }
        report
    }

    /// 二次流の強さ I（`state.spiral`）を `dt` だけ進める（[`secondary_flow`]）。
    /// `[sediment.secondary_flow]` がなければ何もしない。
    ///
    /// `state` は流れを進めた後（時刻 n+1）の状態、`theta_flux` はその流れの計算の θ 平均の
    /// 単位幅流量。河床変動を始める前（`time.morphology_start` まで）も毎ステップ呼ぶ。
    pub fn step_secondary_flow(
        &mut self,
        mesh: &Mesh,
        state: &mut State,
        theta_flux: &[f64],
        dt: f64,
    ) {
        if let Some(solver) = self.secondary_flow.as_mut() {
            solver.step(mesh, state, theta_flux, dt);
        }
    }

    /// 浮遊砂 hC を `dt` だけ進める。`theta_flux` は連続式と同じ θ 平均の単位幅流量。
    pub fn step_suspended(&mut self, mesh: &Mesh, state: &mut State, theta_flux: &[f64], dt: f64) {
        let _ = (mesh, state, theta_flux, dt);
        todo!("hC の移流拡散（沈降は陰的）")
    }

    /// 直前の [`SedimentSolver::compute_bedload`] の掃流砂量で、Exner 式（MORFAC）により河床を
    /// 更新し、安息角補正と水位の調整を行う。`dt` は流れの時間刻み。
    pub fn update_bed(&mut self, mesh: &Mesh, state: &mut State, dt: f64) -> BedUpdateReport {
        let Self {
            params,
            bedload_flux,
            cap_factor,
            avalanche_flux,
            zb_old,
            ..
        } = self;
        let sediment = &params.sediment;
        let factor = sediment.morfac * dt;
        let solid = 1.0 - sediment.porosity;
        zb_old.copy_from_slice(&state.zb);

        // 1. Exner 式
        for (c, zb) in state.zb.iter_mut().enumerate() {
            let net = cell_net_outflow(mesh, c, bedload_flux);
            let z = zb_old[c] - factor * net / (solid * mesh.area_cell[c]);
            // 頭打ちで z ≥ z_ne は保証されている。丸め誤差の分だけそろえる
            *zb = z.max(state.z_ne[c]);
        }

        // 2. 安息角補正
        let tan_phi = sediment.angle_of_repose.to_radians().tan();
        let tolerance = AVALANCHE_TOLERANCE_PER_GRAIN * sediment.grain_size;
        let (iterations, unconverged) = avalanche(
            mesh,
            &mut state.zb,
            &state.z_ne,
            Avalanche {
                tan_phi,
                tolerance,
                background_slope: params.background_slope,
            },
            avalanche_flux,
            cap_factor,
        );

        // 3. 水位の調整と、体積の変化
        let mut report = BedUpdateReport {
            avalanche_iterations: iterations,
            avalanche_unconverged: unconverged,
            ..BedUpdateReport::default()
        };
        for (c, (eta, &zb)) in state.eta.iter_mut().zip(&state.zb).enumerate() {
            let h_old = (*eta - zb_old[c]).max(0.0);
            if h_old >= params.dry_depth {
                *eta = eta.max(zb);
            } else {
                *eta = zb + h_old;
            }
            let h_new = (*eta - zb).max(0.0);
            let area = mesh.area_cell[c];
            report.water_volume_change += area * (h_new - h_old);
            report.bed_volume_change += solid * area * (zb - zb_old[c]);
            report.max_change = report.max_change.max((zb - zb_old[c]).abs());
        }
        report
    }
}

/// 辺 `e` の掃流砂量 q_{b,j}（頭打ちの前）\[m²/s\]。
fn edge_capacity(
    params: &SedimentParams,
    mesh: &Mesh,
    state: &State,
    edge_depth: &[f64],
    scale: f64,
    e: usize,
) -> f64 {
    let [c1, c2] = mesh.cells_on_edge[e];
    let u = state.u_n[e];
    let kind = mesh.edge_kind[e];
    let (depth, slope) = match kind {
        EdgeKind::Wall => return 0.0,
        EdgeKind::Interior => (
            edge_depth[e],
            (state.zb[c2] - state.zb[c1]) / mesh.dc_edge[e]
                - params.background_slope * mesh.normal_edge[e][0],
        ),
        EdgeKind::Upstream => {
            // 法線は領域の外を向くので、流入は u < 0
            let supplied = matches!(
                params.supply,
                Some(SedimentSupply::Equilibrium | SedimentSupply::Constant)
            );
            if !supplied || u >= 0.0 {
                return 0.0;
            }
            (state.depth(c1), boundary_slope(mesh, &state.zb, e))
        }
        EdgeKind::Downstream => {
            if !params.downstream_open || u <= 0.0 {
                return 0.0;
            }
            (state.depth(c1), boundary_slope(mesh, &state.zb, e))
        }
    };
    if depth < params.dry_depth {
        return 0.0;
    }
    let speed = hydro::edge_speed(mesh, state, e);
    if speed <= 0.0 {
        return 0.0;
    }
    let bedload = &params.sediment.bedload;
    let tau = params.shields(speed, depth);
    let phi = formulas::bedload_intensity(
        bedload.formula,
        tau,
        bedload.ripple_factor,
        bedload.critical_shields,
    );
    if phi <= 0.0 {
        return 0.0;
    }
    let gamma = formulas::slope_correction_coefficient(
        tau,
        bedload.critical_shields,
        bedload.static_friction,
        bedload.kinetic_friction,
    );
    let direction = match &params.sediment.secondary_flow {
        Some(secondary) => {
            secondary_flow::bedload_direction(mesh, state, e, speed, secondary.n_star)
        }
        None => u / speed,
    };
    let q = scale * phi * (direction - gamma * slope);
    match kind {
        EdgeKind::Upstream => q.min(0.0),
        EdgeKind::Downstream => q.max(0.0),
        EdgeKind::Interior | EdgeKind::Wall => q,
    }
}

/// 境界辺 `e` の法線方向の河床勾配：内側のセルの河床勾配ベクトル（[`bed_gradient`]）と、
/// 外向きの法線との内積。
fn boundary_slope(mesh: &Mesh, zb: &[f64], e: usize) -> f64 {
    let [gx, gy] = bed_gradient(mesh, zb, mesh.cells_on_edge[e][0]);
    let [nx, ny] = mesh.normal_edge[e];
    gx * nx + gy * ny
}

/// セル `c` の河床勾配ベクトル \[-\]：隣のセル（境界の外は除く）との差の最小二乗
/// （[`crate::mesh::gradient`]）。
///
/// 河床が平面なら厳密。隣のセルが一直線に並ぶ（幅が1セルの水路など）ときは、その方向の
/// 成分だけを求め、直交する成分は 0 とする。周期境界の継ぎ目の向こうのセルは、平行移動した
/// 位置（[`Mesh::neighbor_position_from`]）で数える（河床高は勾配を除いた値のまま）。
/// 開境界の辺でしか使わないので、重みはその都度計算する。
pub fn bed_gradient(mesh: &Mesh, zb: &[f64], c: usize) -> [f64; 2] {
    let weights = gradient::least_squares_weights(mesh, c);
    let mut g = [0.0, 0.0];
    for (&k, [wx, wy]) in mesh.cells_on_cell[c].iter().zip(weights) {
        if k != NONE {
            g[0] += wx * (zb[k] - zb[c]);
            g[1] += wy * (zb[k] - zb[c]);
        }
    }
    g
}

/// セル `c` から辺を通って出る量の合計 Σ_k max(s L f, 0)（`with_length` が偽なら L を掛けない）。
fn cell_outflow(mesh: &Mesh, c: usize, flux: &[f64], with_length: bool) -> f64 {
    mesh.edges_on_cell[c]
        .iter()
        .zip(&mesh.edge_sign_on_cell[c])
        .map(|(&e, &s)| {
            let length = if with_length { mesh.dv_edge[e] } else { 1.0 };
            (s * length * flux[e]).max(0.0)
        })
        .sum()
}

/// セル `c` から辺を通って出る量の正味 Σ_k s L f。
fn cell_net_outflow(mesh: &Mesh, c: usize, flux: &[f64]) -> f64 {
    mesh.edges_on_cell[c]
        .iter()
        .zip(&mesh.edge_sign_on_cell[c])
        .map(|(&e, &s)| s * mesh.dv_edge[e] * flux[e])
        .sum()
}

/// 安息角補正の条件。
#[derive(Debug, Clone, Copy)]
struct Avalanche {
    /// 水中安息角の正接 tan φ。
    tan_phi: f64,
    /// 収束判定：1回の反復で動いた厚さの最大値 \[m\]。
    tolerance: f64,
    /// 周期境界で河床高から除いた縦断勾配 I（高低差に I n_x δ を戻す）。
    background_slope: f64,
}

/// 安息角補正（Jacobi 型の反復）。土砂を動かした反復の回数と、上限までに収束しなかったかを返す。
///
/// `flux` は辺を通して移す河床の体積（法線の向きを正）、`scale` はセルの流出量を縮める倍率の
/// 作業配列。
fn avalanche(
    mesh: &Mesh,
    zb: &mut [f64],
    z_ne: &[f64],
    condition: Avalanche,
    flux: &mut [f64],
    scale: &mut [f64],
) -> (usize, bool) {
    let Avalanche {
        tan_phi,
        tolerance,
        background_slope,
    } = condition;
    for iteration in 0..AVALANCHE_MAX_ITERATIONS {
        // 辺：勾配の超過から移す体積（高い側から低い側へ）
        let mut steep = false;
        for (e, v) in flux.iter_mut().enumerate() {
            let [c1, c2] = mesh.cells_on_edge[e];
            *v = 0.0;
            if c2 == NONE {
                continue;
            }
            // 実際の高低差（周期境界では、除いた勾配の分を戻す）
            let drop =
                zb[c1] - zb[c2] + background_slope * mesh.normal_edge[e][0] * mesh.dc_edge[e];
            let excess = drop.abs() - mesh.dc_edge[e] * tan_phi;
            if excess > 0.0 {
                let inverse_areas = 1.0 / mesh.area_cell[c1] + 1.0 / mesh.area_cell[c2];
                *v = (AVALANCHE_RELAXATION * excess / inverse_areas).copysign(drop);
                steep = true;
            }
        }
        if !steep {
            return (iteration, false);
        }
        // セル：出る量の合計を、非侵食層より上の量で頭打ちにする
        for (c, alpha) in scale.iter_mut().enumerate() {
            let outflow = cell_outflow(mesh, c, flux, false);
            let available = mesh.area_cell[c] * (zb[c] - z_ne[c]).max(0.0);
            *alpha = if outflow > available {
                available / outflow
            } else {
                1.0
            };
        }
        for (e, v) in flux.iter_mut().enumerate() {
            let [c1, c2] = mesh.cells_on_edge[e];
            if *v > 0.0 {
                *v *= scale[c1];
            } else if *v < 0.0 {
                *v *= scale[c2];
            }
        }
        // セル：移した体積を集める
        let mut moved = 0.0_f64;
        for (c, z) in zb.iter_mut().enumerate() {
            let net: f64 = mesh.edges_on_cell[c]
                .iter()
                .zip(&mesh.edge_sign_on_cell[c])
                .map(|(&e, &s)| s * flux[e])
                .sum();
            let dz = -net / mesh.area_cell[c];
            *z = (*z + dz).max(z_ne[c]);
            moved = moved.max(dz.abs());
        }
        if moved < tolerance {
            // 頭打ちで何も動かなかった（床止めの下流面だけが急な）ときは、反復に数えない
            return (iteration + usize::from(moved > 0.0), false);
        }
    }
    (AVALANCHE_MAX_ITERATIONS, true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{BedloadConfig, BedloadFormula, ChannelLayout};
    use crate::mesh::channel::{self, ChannelSpec};

    const MORFAC: f64 = 10.0;

    fn mesh(length: f64, width: f64, cell_size: f64, layout: ChannelLayout) -> Mesh {
        channel::build(&ChannelSpec {
            length,
            width,
            cell_size,
            layout,
        })
        .expect("valid channel")
    }

    /// `examples/sill.toml` と同じ流砂の条件（d = 2 mm、MPM、MORFAC 10）。
    fn params(supply: Option<SedimentSupply>, downstream_open: bool) -> SedimentParams {
        SedimentParams {
            gravity: 9.81,
            water_density: 1000.0,
            kinematic_viscosity: 1.0e-6,
            von_karman: 0.41,
            manning_n: 0.025,
            dry_depth: 1.0e-3,
            sediment: SedimentConfig {
                grain_size: 0.002,
                density: 2650.0,
                porosity: 0.4,
                morfac: MORFAC,
                angle_of_repose: 32.0,
                bedload: BedloadConfig {
                    formula: BedloadFormula::Mpm,
                    critical_shields: 0.047,
                    static_friction: 1.0,
                    kinetic_friction: 0.45,
                    ripple_factor: 1.0,
                },
                suspended: None,
                secondary_flow: None,
            },
            supply,
            downstream_open,
            background_slope: 0.0,
        }
    }

    /// 勾配 `slope` の河床の上の、水深 `h`・流速 `u`（x 方向）の一様流。
    fn uniform_flow(mesh: &Mesh, slope: f64, h: f64, u: f64) -> State {
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        for c in 0..mesh.n_cells() {
            state.zb[c] = -slope * mesh.x_cell[c];
            state.eta[c] = state.zb[c] + h;
        }
        for e in 0..mesh.n_edges() {
            if mesh.edge_kind[e] != EdgeKind::Wall {
                state.u_n[e] = u * mesh.normal_edge[e][0];
            }
        }
        hydro::reconstruct_cell_velocity(mesh, &state.u_n, &mut state.u_cell);
        state
    }

    fn total_bed(mesh: &Mesh, zb: &[f64]) -> f64 {
        (0..mesh.n_cells()).map(|c| mesh.area_cell[c] * zb[c]).sum()
    }

    /// 決まった乱数列（線形合同法）で [0, 1) の値を返す。
    fn pseudo_random(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (*seed >> 11) as f64 / (1u64 << 53) as f64
    }

    /// 一様な勾配の等流では、辺の流砂量は q_b0 (1 + γ I) n_x で、境界辺（平衡給砂・自由流出）も
    /// 内部の辺と一致する。q_b0 は手計算の値。
    #[test]
    fn uniform_flow_gives_the_same_bedload_on_every_edge() {
        let (slope, h, u) = (0.002, 1.069235, 1.870496);
        for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
            let mesh = mesh(20.0, 6.0, 1.0, layout);
            let state = uniform_flow(&mesh, slope, h, u);
            let p = params(Some(SedimentSupply::Equilibrium), true);
            let mut solver = SedimentSolver::new(&mesh, p.clone());
            let report = solver.compute_bedload(&mesh, &state, 0.1);
            assert_eq!(report.capped_cells, 0);

            // τ* = n² u² / (h^{1/3} s d)、q_b0 = √(s g d³) 8 (τ* − τ*_c)^{3/2}
            let s = 1.65;
            let tau = 0.025_f64.powi(2) * u * u / (h.cbrt() * s * 0.002);
            let qb0 = (s * 9.81 * 0.002_f64.powi(3)).sqrt() * 8.0 * (tau - 0.047).powf(1.5);
            let gamma = (0.047 / (0.45 * tau)).sqrt();
            assert!((tau - 0.648).abs() < 1e-3, "tau* = {tau}");
            for e in 0..mesh.n_edges() {
                let expected = match mesh.edge_kind[e] {
                    EdgeKind::Wall => 0.0,
                    _ => qb0 * (1.0 + gamma * slope) * mesh.normal_edge[e][0],
                };
                let q = solver.bedload_flux()[e];
                assert!(
                    (q - expected).abs() <= 1e-12 * qb0,
                    "{layout:?} edge {e} ({:?}): {q} vs {expected}",
                    mesh.edge_kind[e]
                );
            }
            // 給砂量と流出量は等しく、内部のセルの河床は変わらない
            assert!((report.inflow - report.outflow).abs() <= 1e-12 * report.inflow);
            let mut state = state;
            let bed = solver.update_bed(&mesh, &mut state, 0.1);
            assert!(bed.max_change < 1e-15, "{layout:?}: {}", bed.max_change);
            assert_eq!(bed.avalanche_iterations, 0);
        }
    }

    /// 砂漣係数 μ：等流の辺の流砂量は √(s g d³) 8 (μτ* − τ*_c)^{3/2} (1 + γ I) n_x で、
    /// γ は全掃流力 τ* のまま（ADR-0021）。給砂量と流出量も同じ式で、河床は変わらない。
    #[test]
    fn ripple_factor_reduces_the_bedload_but_not_the_slope_correction() {
        let (slope, h, u, mu) = (0.002, 1.069235, 1.870496, 0.35);
        let mesh = mesh(20.0, 6.0, 1.0, ChannelLayout::Hexagonal);
        let mut state = uniform_flow(&mesh, slope, h, u);
        let mut p = params(Some(SedimentSupply::Equilibrium), true);
        p.sediment.bedload.ripple_factor = mu;
        let mut solver = SedimentSolver::new(&mesh, p.clone());
        let report = solver.compute_bedload(&mesh, &state, 0.1);

        let tau = p.shields(u, h);
        let qb0 = p.bedload_scale() * 8.0 * (mu * tau - 0.047).powf(1.5);
        let gamma = (0.047 / (0.45 * tau)).sqrt();
        for e in 0..mesh.n_edges() {
            let expected = match mesh.edge_kind[e] {
                EdgeKind::Wall => 0.0,
                _ => qb0 * (1.0 + gamma * slope) * mesh.normal_edge[e][0],
            };
            let q = solver.bedload_flux()[e];
            assert!(
                (q - expected).abs() <= 1e-12 * qb0,
                "edge {e} ({:?}): {q} vs {expected}",
                mesh.edge_kind[e]
            );
        }
        assert!((report.inflow - report.outflow).abs() <= 1e-12 * report.inflow);
        let bed = solver.update_bed(&mesh, &mut state, 0.1);
        assert!(bed.max_change < 1e-15, "{}", bed.max_change);

        // μ τ* ≤ τ*_c なら掃流砂は動かない
        p.sediment.bedload.ripple_factor = 0.04 / tau;
        let mut still = SedimentSolver::new(&mesh, p);
        let report = still.compute_bedload(&mesh, &state, 0.1);
        assert!(still.bedload_flux().iter().all(|&q| q == 0.0));
        assert_eq!(report.inflow, 0.0);
    }

    /// 二次流の設定（輸送方程式型、N* = 7）。
    fn with_secondary_flow(p: SedimentParams) -> SedimentParams {
        let mut p = p;
        p.sediment.secondary_flow = Some(crate::config::SecondaryFlowConfig {
            kind: crate::config::SecondaryFlowKind::Transport,
            n_star: 7.0,
            diffusivity_coef: 0.15,
        });
        p
    }

    /// 二次流：一様な勾配の等流に一様な I を与えると、辺の掃流砂量は
    /// q_b0 (cos δ n_x + sin δ n_y + γ I_b n_x)（tan δ = N* I / |u|。内岸側＝左へ δ 回る）。
    /// 境界辺も同じ式で、I = 0 なら二次流を考えないときと一致する。
    #[test]
    fn secondary_flow_turns_the_bedload_towards_the_inner_bank() {
        let (slope, h, u, intensity) = (0.002, 1.069235, 1.870496, 0.05);
        for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
            let mesh = mesh(20.0, 6.0, 1.0, layout);
            let mut state = uniform_flow(&mesh, slope, h, u);
            let p = with_secondary_flow(params(Some(SedimentSupply::Equilibrium), true));
            let mut solver = SedimentSolver::new(&mesh, p.clone());
            let mut plain =
                SedimentSolver::new(&mesh, params(Some(SedimentSupply::Equilibrium), true));
            solver.compute_bedload(&mesh, &state, 0.1);
            plain.compute_bedload(&mesh, &state, 0.1);
            assert_eq!(solver.bedload_flux(), plain.bedload_flux());

            state.spiral.fill(intensity);
            solver.compute_bedload(&mesh, &state, 0.1);
            let tau = p.shields(u, h);
            let qb0 = p.bedload_scale() * 8.0 * (tau - 0.047).powf(1.5);
            let gamma = (0.047 / (0.45 * tau)).sqrt();
            let delta = (7.0 * intensity / u).atan();
            for e in 0..mesh.n_edges() {
                let [nx, ny] = mesh.normal_edge[e];
                let expected = match mesh.edge_kind[e] {
                    EdgeKind::Wall => 0.0,
                    _ => qb0 * (delta.cos() * nx + delta.sin() * ny + gamma * slope * nx),
                };
                let q = solver.bedload_flux()[e];
                assert!(
                    (q - expected).abs() <= 1e-12 * qb0,
                    "{layout:?} edge {e} ({:?}): {q} vs {expected}",
                    mesh.edge_kind[e]
                );
            }
        }
    }

    /// 十分発達した湾曲部の横断勾配（合成の状態）：x 方向の一様流に一様な I を与え、河床と水面を
    /// 横断方向に勾配 S = sin δ / γ で傾ける（水深は一様）と、横断方向の掃流砂が 0 になり、
    /// 河床は変わらない（周期境界の正方形配置）。I = 0 なら勾配を下る掃流砂で側壁沿いが変わる。
    #[test]
    fn transverse_bed_slope_balances_the_secondary_flow() {
        let (h, u, intensity) = (1.0, 1.8, 0.02);
        let mesh = channel::build_periodic(&ChannelSpec {
            length: 12.0,
            width: 6.0,
            cell_size: 1.0,
            layout: ChannelLayout::Rectangular,
        })
        .expect("valid periodic channel");
        let p = with_secondary_flow(params(None, false));
        let tau = p.shields(u, h);
        let gamma = (0.047 / (0.45 * tau)).sqrt();
        let delta = (7.0 * intensity / u).atan();
        let transverse = delta.sin() / gamma;
        // 内岸（左、y の大きい側）が高い
        let tilted = |spiral: f64| {
            let mut state = uniform_flow(&mesh, 0.0, h, u);
            for c in 0..mesh.n_cells() {
                state.zb[c] = transverse * mesh.y_cell[c];
                state.eta[c] = state.zb[c] + h;
            }
            state.spiral.fill(spiral);
            state
        };
        let mut state = tilted(intensity);
        let mut solver = SedimentSolver::new(&mesh, p.clone());
        solver.compute_bedload(&mesh, &state, 1.0);
        let qb0 = p.bedload_scale() * 8.0 * (tau - 0.047).powf(1.5);
        for e in (0..mesh.n_edges()).filter(|&e| mesh.normal_edge[e][1].abs() > 0.5) {
            let q = solver.bedload_flux()[e];
            assert!(q.abs() <= 1e-14 * qb0, "edge {e}: {q}");
        }
        let report = solver.update_bed(&mesh, &mut state, 1.0);
        assert!(report.max_change < 1e-15, "{report:?}");
        assert_eq!(report.avalanche_iterations, 0);

        let mut state = tilted(0.0);
        solver.compute_bedload(&mesh, &state, 1.0);
        let report = solver.update_bed(&mesh, &mut state, 1.0);
        assert!(report.max_change > 1.0e-4, "{report:?}");
    }

    /// 周期境界：勾配を除いた一様な河床の等流でも、辺の掃流砂量は除いた勾配の分の補正を含む
    /// q_b0 (1 + γ I) n_x になる（勾配のある水路の等流と同じ。ADR-0014）。
    #[test]
    fn periodic_uniform_flow_includes_the_removed_slope() {
        let (slope, h, u) = (0.002, 1.069235, 1.870496);
        for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
            let mesh = channel::build_periodic(&ChannelSpec {
                length: 20.0,
                width: 6.0,
                cell_size: 1.0,
                layout,
            })
            .expect("valid periodic channel");
            // 勾配を除いた河床は一様
            let state = uniform_flow(&mesh, 0.0, h, u);
            let p = SedimentParams {
                background_slope: slope,
                ..params(None, false)
            };
            let mut solver = SedimentSolver::new(&mesh, p);
            solver.compute_bedload(&mesh, &state, 0.1);
            let s = 1.65;
            let tau = 0.025_f64.powi(2) * u * u / (h.cbrt() * s * 0.002);
            let qb0 = (s * 9.81 * 0.002_f64.powi(3)).sqrt() * 8.0 * (tau - 0.047).powf(1.5);
            let gamma = (0.047 / (0.45 * tau)).sqrt();
            for e in 0..mesh.n_edges() {
                let expected = match mesh.edge_kind[e] {
                    EdgeKind::Wall => 0.0,
                    _ => qb0 * (1.0 + gamma * slope) * mesh.normal_edge[e][0],
                };
                let q = solver.bedload_flux()[e];
                assert!(
                    (q - expected).abs() <= 1e-12 * qb0,
                    "{layout:?} edge {e}: {q} vs {expected}"
                );
            }
            let mut state = state;
            let bed = solver.update_bed(&mesh, &mut state, 0.1);
            assert!(bed.max_change < 1e-15, "{layout:?}: {}", bed.max_change);
        }
    }

    /// 周期境界の安息角補正は、勾配を除いた河床に除いた勾配を戻した実際の高低差で判定する。
    /// 勾配のある水路（勾配を除かない）と同じ量を動かす。
    #[test]
    fn periodic_avalanching_uses_the_actual_drop() {
        let (slope, raise) = (0.5, 0.2);
        let spec = ChannelSpec {
            length: 9.0,
            width: 3.0,
            cell_size: 1.0,
            layout: ChannelLayout::Rectangular,
        };
        let periodic = channel::build_periodic(&spec).expect("valid periodic channel");
        let open = channel::build(&spec).expect("valid channel");
        let peak = 13; // 中央の行の x = 4.5 m
        let run = |mesh: &Mesh, background: f64| {
            let mut state = State::new(mesh.n_cells(), mesh.n_edges());
            for c in 0..mesh.n_cells() {
                // 勾配を除いた河床（周期境界）、または勾配のある実際の河床
                state.zb[c] = if background > 0.0 {
                    0.0
                } else {
                    -slope * mesh.x_cell[c]
                };
            }
            state.zb[peak] += raise;
            let zb0 = state.zb.clone();
            let p = SedimentParams {
                background_slope: background,
                ..params(None, false)
            };
            let mut solver = SedimentSolver::new(mesh, p);
            solver.compute_bedload(mesh, &state, 1.0);
            let report = solver.update_bed(mesh, &mut state, 1.0);
            let change: Vec<f64> = (0..mesh.n_cells()).map(|c| state.zb[c] - zb0[c]).collect();
            (report, change)
        };
        let (with_trend, change_periodic) = run(&periodic, slope);
        let (_, change_open) = run(&open, 0.0);
        // 下流側の高低差 0.5 + 0.2 = 0.7 > tan 32° なので崩れる（勾配を除いた差 0.2 だけなら崩れない）
        assert!(with_trend.avalanche_iterations > 0, "{with_trend:?}");
        assert!(change_periodic[peak] < -1e-3, "{}", change_periodic[peak]);
        for c in 0..periodic.n_cells() {
            assert!(
                (change_periodic[c] - change_open[c]).abs() < 1e-14,
                "cell {c}: {} vs {}",
                change_periodic[c],
                change_open[c]
            );
        }
    }

    /// 一定の給砂量は、最初のステップの平衡給砂量の合計を保ち、上流端の流れが変わっても
    /// 変えない（辺への配り方は今の平衡給砂量に比例）。平衡給砂は流れに追従する。
    #[test]
    fn constant_supply_keeps_the_first_equilibrium_rate() {
        let mesh = mesh(10.0, 2.0, 1.0, ChannelLayout::Rectangular);
        let fast = uniform_flow(&mesh, 0.002, 1.0, 2.0);
        let slow = uniform_flow(&mesh, 0.002, 1.0, 1.6);
        let mut constant = SedimentSolver::new(&mesh, params(Some(SedimentSupply::Constant), true));
        let mut equilibrium =
            SedimentSolver::new(&mesh, params(Some(SedimentSupply::Equilibrium), true));
        let first = constant.compute_bedload(&mesh, &fast, 0.1);
        let reference = equilibrium.compute_bedload(&mesh, &fast, 0.1);
        assert!(first.inflow > 0.0);
        assert_eq!(first.inflow, reference.inflow);
        let later = constant.compute_bedload(&mesh, &slow, 0.1);
        assert!((later.inflow / first.inflow - 1.0).abs() < 1e-14);
        let following = equilibrium.compute_bedload(&mesh, &slow, 0.1);
        assert!(following.inflow < 0.5 * first.inflow, "{following:?}");
        // 辺への配り方は、今の平衡給砂量に比例する
        let ratio = later.inflow / following.inflow;
        for e in (0..mesh.n_edges()).filter(|&e| mesh.edge_kind[e] == EdgeKind::Upstream) {
            let (q, q_eq) = (constant.bedload_flux()[e], equilibrium.bedload_flux()[e]);
            assert!((q - ratio * q_eq).abs() <= 1e-14 * q.abs(), "edge {e}");
        }
    }

    /// 給砂しない設定と壁の境界では、境界辺の流砂量は 0。流れの向きが逆なら流出もしない。
    #[test]
    fn closed_and_unsupplied_boundaries_carry_no_bedload() {
        let mesh = mesh(10.0, 2.0, 1.0, ChannelLayout::Rectangular);
        let state = uniform_flow(&mesh, 0.002, 1.0, 2.0);
        let mut solver = SedimentSolver::new(&mesh, params(Some(SedimentSupply::None), false));
        let report = solver.compute_bedload(&mesh, &state, 0.1);
        assert_eq!((report.inflow, report.outflow), (0.0, 0.0));
        for e in (0..mesh.n_edges()).filter(|&e| mesh.is_boundary_edge(e)) {
            assert_eq!(solver.bedload_flux()[e], 0.0);
        }
        // 逆流（下流端から流入）では、下流端の辺から土砂は入らない
        let reverse = uniform_flow(&mesh, 0.0, 1.0, -2.0);
        let mut solver = SedimentSolver::new(&mesh, params(None, true));
        let report = solver.compute_bedload(&mesh, &reverse, 0.1);
        assert_eq!((report.inflow, report.outflow), (0.0, 0.0));
    }

    /// 壁で囲んだ水域では、頭打ちと安息角補正を含めても、河床の土砂の総量が丸め誤差の範囲で
    /// 保存される（S2）。非侵食層より下は掘らない。
    #[test]
    fn exner_with_caps_and_avalanching_conserves_sediment() {
        let mesh = mesh(12.0, 6.0, 0.5, ChannelLayout::Hexagonal);
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        let mut seed = 12345;
        for c in 0..mesh.n_cells() {
            state.zb[c] = 0.3 * pseudo_random(&mut seed);
            state.eta[c] = 1.5;
            // 中央の帯は、薄い土砂が非侵食層に乗っている
            if (5.0..7.0).contains(&mesh.x_cell[c]) {
                state.z_ne[c] = state.zb[c] - 0.002 * pseudo_random(&mut seed);
            }
        }
        // 高い山（安息角を超える勾配）
        let peak = (0..mesh.n_cells())
            .min_by(|&a, &b| {
                let d = |c: usize| (mesh.x_cell[c] - 3.0).hypot(mesh.y_cell[c] - 3.0);
                d(a).total_cmp(&d(b))
            })
            .expect("cells");
        state.zb[peak] += 1.0;
        for e in 0..mesh.n_edges() {
            if !mesh.is_boundary_edge(e) {
                state.u_n[e] = 3.0 * (pseudo_random(&mut seed) - 0.3);
            }
        }
        hydro::reconstruct_cell_velocity(&mesh, &state.u_n, &mut state.u_cell);

        let mut solver = SedimentSolver::new(&mesh, params(None, false));
        let before = total_bed(&mesh, &state.zb);
        let bedload = solver.compute_bedload(&mesh, &state, 0.5);
        let bed = solver.update_bed(&mesh, &mut state, 0.5);
        assert!(bedload.capped_cells > 0, "the cap must be active");
        assert!(bed.avalanche_iterations > 0, "avalanching must be active");
        assert!(!bed.avalanche_unconverged);
        let after = total_bed(&mesh, &state.zb);
        assert!(
            (after - before).abs() <= 1e-12 * before.abs(),
            "{before} -> {after}"
        );
        assert!(bed.bed_volume_change.abs() <= 1e-12 * before.abs());
        assert!((0..mesh.n_cells()).all(|c| state.zb[c] >= state.z_ne[c]));
    }

    /// 薄い土砂のセルから強く流出させると、河床はちょうど非侵食層の高さで止まる（S3）。
    /// 上流側は土砂のない非侵食層（床止め）で、流入はない。
    #[test]
    fn non_erodible_layer_stops_erosion_exactly() {
        let mesh = mesh(10.0, 1.0, 1.0, ChannelLayout::Rectangular);
        let mut state = uniform_flow(&mesh, 0.0, 0.5, 3.0);
        let c = (0..mesh.n_cells())
            .find(|&c| (mesh.x_cell[c] - 4.5).abs() < 1e-12)
            .expect("a cell at x = 4.5 m");
        for k in 0..mesh.n_cells() {
            if mesh.x_cell[k] < mesh.x_cell[c] {
                state.z_ne[k] = state.zb[k];
            }
        }
        state.z_ne[c] = state.zb[c] - 1.0e-4;
        let mut solver = SedimentSolver::new(&mesh, params(Some(SedimentSupply::None), true));
        let bedload = solver.compute_bedload(&mesh, &state, 1.0);
        // 床止めの 4 セルと、薄い土砂のセル
        assert_eq!(bedload.capped_cells, 5);
        let before = total_bed(&mesh, &state.zb);
        let report = solver.update_bed(&mesh, &mut state, 1.0);
        assert_eq!(state.zb[c], state.z_ne[c]);
        // 失われた量は下流端からの流出だけ
        let lost = 0.6 * (before - total_bed(&mesh, &state.zb));
        assert!((lost - bedload.outflow).abs() <= 1e-12 * bedload.outflow);
        assert!((report.bed_volume_change + bedload.outflow).abs() <= 1e-12 * bedload.outflow);
    }

    /// 1セルだけ高い山は、勾配が安息角以下になるまで崩れ、体積は保存される。高い側が
    /// 非侵食層の段差（床止めの下流面）では何も動かない（S4）。
    #[test]
    fn avalanching_relaxes_slopes_and_respects_the_non_erodible_layer() {
        let mesh = mesh(9.0, 9.0, 1.0, ChannelLayout::Hexagonal);
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        let peak = 40;
        state.zb[peak] = 3.0;
        let mut solver = SedimentSolver::new(&mesh, params(None, false));
        let before = total_bed(&mesh, &state.zb);
        solver.compute_bedload(&mesh, &state, 1.0);
        let report = solver.update_bed(&mesh, &mut state, 1.0);
        assert!(report.avalanche_iterations > 1 && !report.avalanche_unconverged);
        assert!((total_bed(&mesh, &state.zb) - before).abs() <= 1e-12 * before);
        let tan_phi = 32.0_f64.to_radians().tan();
        let tolerance = 3.0 * AVALANCHE_TOLERANCE_PER_GRAIN * 0.002;
        for e in (0..mesh.n_edges()).filter(|&e| !mesh.is_boundary_edge(e)) {
            let [c1, c2] = mesh.cells_on_edge[e];
            let slope = (state.zb[c2] - state.zb[c1]).abs() / mesh.dc_edge[e];
            assert!(
                slope <= tan_phi + tolerance / mesh.dc_edge[e],
                "edge {e}: {slope}"
            );
        }

        // 高い側に土砂がない段差は崩れない
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        state.zb[peak] = 3.0;
        state.z_ne[peak] = 3.0;
        let report = solver.update_bed(&mesh, &mut state, 1.0);
        assert_eq!(report.avalanche_iterations, 0);
        assert_eq!(report.max_change, 0.0);
    }

    /// 河床の変化に対する水位の調整：湛水セルは水位を保ち、水面を超えたら干上がり、ほぼ乾いた
    /// セルは水深を保つ。水の体積の変化を返す（S5）。
    #[test]
    fn water_level_is_kept_in_wet_cells() {
        let mesh = mesh(4.0, 1.0, 1.0, ChannelLayout::Rectangular);
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        // セル 0：湛水（水深 1 m）、1：湛水（水深 0.05 m）、2：ほぼ乾いた（0.5 mm）、3：乾いた
        let depth = [1.0, 0.05, 0.0005, 0.0];
        state.eta.copy_from_slice(&depth);
        let mut solver = SedimentSolver::new(&mesh, params(None, false));
        // 辺の流砂量を直接与える：セル 0 → 1 → 2 → 3 へ、厚さで 0.1, 0.2, 0.1 m 分
        let factor = MORFAC * 1.0 / 0.6;
        for e in 0..mesh.n_edges() {
            let [c1, c2] = mesh.cells_on_edge[e];
            if c2 != NONE {
                let (lo, sign) = if c1 < c2 { (c1, 1.0) } else { (c2, -1.0) };
                let thickness = [0.1, 0.2, 0.1][lo];
                solver.bedload_flux[e] = sign * thickness / factor;
            }
        }
        let report = solver.update_bed(&mesh, &mut state, 1.0);
        let expected_zb = [-0.1, -0.1, 0.1, 0.1];
        for (c, (&zb, expected)) in state.zb.iter().zip(expected_zb).enumerate() {
            assert!((zb - expected).abs() < 1e-14, "zb[{c}]");
        }
        assert!((state.eta[0] - 1.0).abs() < 1e-14, "wet: the level is kept");
        assert!(
            (state.eta[1] - 0.05).abs() < 1e-14,
            "wet: the level is kept"
        );
        assert!(
            (state.eta[2] - 0.1005).abs() < 1e-14,
            "nearly dry: the depth is kept"
        );
        assert!((state.eta[3] - 0.1).abs() < 1e-14, "dry: eta = zb");
        // 水の体積の変化：セル 0 と 1 が 0.1 m ずつ深くなった
        assert!((report.water_volume_change - 0.2).abs() < 1e-14);

        // 湛水セルで、河床が水面を超えると干上がる（セル 2 は 0.1 m 上がる）
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        state.eta[2] = 0.05;
        let report = solver.update_bed(&mesh, &mut state, 1.0);
        assert!((state.zb[2] - 0.1).abs() < 1e-14);
        assert_eq!(state.eta[2], state.zb[2]);
        assert!((report.water_volume_change + 0.05).abs() < 1e-14);
    }

    /// 河床勾配の最小二乗は、平面なら厳密で、隣が一直線に並ぶ（幅1セル）ときもその方向の
    /// 成分を返す。
    #[test]
    fn bed_gradient_is_exact_for_planes() {
        for (width, layout) in [
            (5.0, ChannelLayout::Hexagonal),
            (5.0, ChannelLayout::Rectangular),
            (1.0, ChannelLayout::Rectangular),
        ] {
            let mesh = mesh(8.0, width, 1.0, layout);
            let zb: Vec<f64> = (0..mesh.n_cells())
                .map(|c| 0.3 - 0.02 * mesh.x_cell[c] + 0.05 * mesh.y_cell[c])
                .collect();
            for c in 0..mesh.n_cells() {
                let [gx, gy] = bed_gradient(&mesh, &zb, c);
                let expected_gy = if width > 1.0 { 0.05 } else { 0.0 };
                assert!((gx + 0.02).abs() < 1e-12, "{layout:?} cell {c}: {gx}");
                assert!(
                    (gy - expected_gy).abs() < 1e-12,
                    "{layout:?} cell {c}: {gy}"
                );
            }
        }
    }

    /// 周期的な水路では、継ぎ目の向こうのセルを平行移動した位置で数える。正方形配置の
    /// 最小二乗は x 方向の中心差分 (z_E − z_W) / 2Δx になり、継ぎ目のセルでも同じ。
    #[test]
    fn bed_gradient_wraps_around_the_periodic_seam() {
        let length = 8.0;
        let mesh = channel::build_periodic(&ChannelSpec {
            length,
            width: 3.0,
            cell_size: 1.0,
            layout: ChannelLayout::Rectangular,
        })
        .expect("valid periodic channel");
        let wave = |x: f64| (2.0 * std::f64::consts::PI * x / length).cos();
        let zb: Vec<f64> = mesh.x_cell.iter().map(|&x| wave(x)).collect();
        for c in 0..mesh.n_cells() {
            let x = mesh.x_cell[c];
            let [gx, gy] = bed_gradient(&mesh, &zb, c);
            let expected = (wave(x + 1.0) - wave(x - 1.0)) / 2.0;
            assert!(
                (gx - expected).abs() < 1e-14,
                "cell {c}: {gx} vs {expected}"
            );
            assert!(gy.abs() < 1e-14, "cell {c}: {gy}");
        }
    }
}
