//! 流れの計算：Casulli 型 θ 法の半陰解法（C 格子）。
//!
//! ## 離散式
//!
//! 詳細と導出は `docs/design/overview.md` の第4節。記号は i がセル、j が辺、
//! 1(j), 2(j) が `cells_on_edge[j][0]`, `[1]`、s_ij が `edge_sign_on_cell`、
//! k(i, j) が辺 j を挟んだ隣のセル。
//!
//! 辺 j の法線流速は、移流を陽的に、水位勾配を θ 法で、摩擦を陰的に扱う。
//!
//! ```text
//! u_j^{n+1} = G_j − θ g Δt / (δ_j (1 + Δt γ_j)) · (η_{2(j)}^{n+1} − η_{1(j)}^{n+1})
//! G_j = [F u_j^n − (1 − θ) g Δt / δ_j · (η_{2(j)}^n − η_{1(j)}^n)] / (1 + Δt γ_j)
//! γ_j = g n² |u|_j / H_j^{4/3}
//! ```
//!
//! |u|_j は、法線成分 u_j^n と、両側のセルの Perot 流速の平均の接線成分から求める
//! （一様流なら厳密に元の速さになる）。これを連続式に代入すると、水位について次の
//! 連立方程式が得られる。
//!
//! ```text
//! V_i(η_i) + Σ_j w_j (η_i − η_{k(i,j)}) = b_i
//! w_j = g θ² Δt² L_j H_j / (δ_j (1 + Δt γ_j))
//! b_i = V_i(η_i^n) − Δt Σ_j s_ij L_j H_j (θ G_j + (1 − θ) u_j^n)
//! ```
//!
//! V_i(η) = A_i max(η − z_{b,i}, 0) は wet/dry の境目（η = z_b）で折れる区分線形関数で、
//! これを Newton 法（Casulli, 2009; Casulli & Zanolli, 2012）で解く。
//!
//! ## wet/dry と Newton 法
//!
//! - 辺の水深 H_j が `dry_depth` 未満の辺は乾いているとみなし、u_j^{n+1} = 0、w_j = 0 とする
//!   （[`edge_depths`]）。体積関数 V_i には閾値を入れない（水の体積を正確に保つため）。
//! - Newton 法は修正量の形 (P + T) δ = r で解く。r = b − V(η^m) − T η^m は残差、P は
//!   η_i^m ≥ z_{b,i} のセルで A_i、それ以外で 0 の対角行列、T は重み w_j のグラフラプラシアン。
//!   PCG（[`crate::linsolve::pcg`]）の相対判定値はその時点の残差に対してかかるので、反復が
//!   進むほど解の精度が上がる。収束判定は max_i |r_i| / A_i ≤ `newton_tol`。
//! - 初期値は η^n とする。ステップの最後に乾いたセルの水位を河床高に揃える
//!   （η_i = max(η_i, z_{b,i})。V_i は変わらない）ので、全セルで η^n ≥ z_b となり、最初の段の
//!   係数行列は diag(A) + T（対称正定値）になる。V が凸なので、2段目以降の反復値は解以上に
//!   あり、単調に減少して有限回で収束する（`docs/design/overview.md` の第4.3節）。
//!
//! ## 1ステップの手順
//!
//! 1. 辺の役割（内部・壁）と辺の水深 H_j^n を決める（[`edge_depths`]）。H_j は風上側の水位から
//!    辺の河床高 max(z_{b,1}, z_{b,2}) を引いた値（Stelling & Duinmeijer, 2003）。
//! 2. 移流を評価して F u_j^n を求める。**移流は次の段階で実装する**（[`advection`]）。
//!    それまでは F u_j^n = u_j^n とする。
//! 3. 摩擦係数 γ_j、G_j、重み w_j、右辺 b_i を組み立てる。
//! 4. Newton 法＋行列フリー PCG で η^{n+1} を求める（[`WaterLevelOperator`]）。
//! 5. u_j^{n+1} を更新し、θ 平均の単位幅流量 q_j^θ = H_j^n (θ u_j^{n+1} + (1 − θ) u_j^n) を
//!    保存する。浮遊砂の移流に同じ流量を使い、濃度が一様な状態が保たれるようにする。
//! 6. 乾いたセルの水位を河床高に揃え、Perot 法でセル中心の流速を再構成する
//!    （[`reconstruct_cell_velocity`]）。
//!
//! ## 参考文献
//!
//! - Casulli & Walters (2000): 非構造格子上の半陰解法
//! - Casulli (2009), Casulli & Zanolli (2012): wet/dry と Newton 法
//! - Stelling & Duinmeijer (2003), Kramer & Stelling (2008): 保存型風上差分
//! - Kernkamp et al. (2011): 非構造格子への適用（D-Flow FM）
//! - Perot (2000): セル中心流速の再構成

use crate::config::{BoundaryConfig, Config, DownstreamBoundary, SolverConfig, UpstreamBoundary};
use crate::linsolve::{self, LinearOperator, PcgError, PcgSettings, PcgWorkspace};
use crate::mesh::{EdgeKind, Mesh, NONE};
use crate::state::State;
use std::error::Error;
use std::fmt;

/// 流れの計算に使うパラメータ。物理定数は設定から受け取る。
#[derive(Debug, Clone, PartialEq)]
pub struct HydroParams {
    /// 重力加速度 g \[m/s²\]。
    pub gravity: f64,
    /// θ 法の重み θ。
    pub theta: f64,
    /// マニングの粗度係数 n [s/m^(1/3)]。
    pub manning_n: f64,
    /// wet/dry を判定する水深の閾値 \[m\]。
    pub dry_depth: f64,
    /// 反復解法の設定。
    pub solver: SolverConfig,
}

impl HydroParams {
    /// 設定から作る。
    pub fn from_config(config: &Config) -> Self {
        Self {
            gravity: config.physics.gravity,
            theta: config.hydro.theta,
            manning_n: config.hydro.manning_n,
            dry_depth: config.hydro.dry_depth,
            solver: config.hydro.solver.clone(),
        }
    }
}

/// 1ステップの計算の記録。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct StepReport {
    /// Newton 法の反復回数（解いた線形系の数）。
    pub newton_iterations: usize,
    /// PCG の反復回数の合計。
    pub pcg_iterations: usize,
    /// 最終的な残差（体積残差をセル面積で割った値の最大値）\[m\]。乾いたセルを含む。
    pub residual: f64,
}

/// 流れの計算の失敗。
#[derive(Debug, Clone, PartialEq)]
pub enum HydroError {
    /// Newton 法が収束しなかった。
    NewtonNotConverged {
        /// 反復回数。
        iterations: usize,
        /// 最終的な残差 \[m\]。
        residual: f64,
    },
    /// PCG が失敗した。
    Pcg(PcgError),
}

impl fmt::Display for HydroError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NewtonNotConverged {
                iterations,
                residual,
            } => write!(
                f,
                "Newton iteration did not converge in {iterations} iterations (residual {residual:e} m)"
            ),
            Self::Pcg(e) => write!(f, "linear solver failed: {e}"),
        }
    }
}

impl Error for HydroError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::NewtonNotConverged { .. } => None,
            Self::Pcg(e) => Some(e),
        }
    }
}

impl From<PcgError> for HydroError {
    fn from(e: PcgError) -> Self {
        Self::Pcg(e)
    }
}

/// 辺の役割（辺の種類に境界条件を当てはめたもの）。
#[derive(Debug, Clone, Copy, PartialEq)]
enum EdgeRole {
    /// 2つのセルに挟まれた内部の辺。
    Interior,
    /// 壁（法線流速 0）。
    Wall,
}

/// 辺 `e` の役割を決める。
fn edge_role(mesh: &Mesh, boundary: &BoundaryConfig, e: usize) -> EdgeRole {
    match mesh.edge_kind[e] {
        EdgeKind::Interior => EdgeRole::Interior,
        EdgeKind::Wall => EdgeRole::Wall,
        EdgeKind::Upstream => match boundary.upstream {
            UpstreamBoundary::Wall => EdgeRole::Wall,
            UpstreamBoundary::Discharge { .. } => todo!("上流端の流量境界"),
        },
        EdgeKind::Downstream => match boundary.downstream {
            DownstreamBoundary::Wall => EdgeRole::Wall,
            DownstreamBoundary::WaterLevel { .. } | DownstreamBoundary::Depth { .. } => {
                todo!("下流端の水位・水深境界")
            }
        },
    }
}

/// 水位勾配を評価する辺なら (η_{1(j)}, η_{2(j)}, δ_j) を返す。壁なら `None`。
fn gradient_points(mesh: &Mesh, role: EdgeRole, eta: &[f64], e: usize) -> Option<(f64, f64, f64)> {
    let [c1, c2] = mesh.cells_on_edge[e];
    match role {
        EdgeRole::Interior => Some((eta[c1], eta[c2], mesh.dc_edge[e])),
        EdgeRole::Wall => None,
    }
}

/// 半陰解法のソルバー。作業配列を持ち回して、ステップごとの確保を避ける。
#[derive(Debug, Clone)]
pub struct HydroSolver {
    params: HydroParams,
    /// 辺の役割
    role: Vec<EdgeRole>,
    /// 境界辺の外側の水位（外に水がなければ −∞）\[m\]
    outside_level: Vec<f64>,
    /// 辺の水深 H_j^n \[m\]
    edge_depth: Vec<f64>,
    /// 移流を評価した F u_j^n \[m/s\]
    advected: Vec<f64>,
    /// 摩擦係数 γ_j \[1/s\]
    friction: Vec<f64>,
    /// 陽的な部分 G_j \[m/s\]
    explicit_part: Vec<f64>,
    /// 水位方程式の辺の重み w_j \[m³/m\]
    weight: Vec<f64>,
    /// 右辺に入る辺の流量 L_j H_j (θ G_j + (1 − θ) u_j^n) \[m³/s\]
    explicit_flux: Vec<f64>,
    /// 水位方程式の右辺 b_i \[m³\]
    rhs: Vec<f64>,
    /// Newton 法の係数行列の対角 P_i \[m²\]
    wet_area: Vec<f64>,
    /// Newton 法の残差 r_i = b_i − V_i − (Tη)_i \[m³\]
    residual: Vec<f64>,
    /// Newton 法の修正量 δ_i \[m\]
    correction: Vec<f64>,
    /// θ 平均の単位幅流量 q_j^θ \[m²/s\]
    theta_flux: Vec<f64>,
    /// PCG の作業配列
    pcg: PcgWorkspace,
}

impl HydroSolver {
    /// 格子の大きさに合わせて作業配列を確保する。
    pub fn new(mesh: &Mesh, params: HydroParams) -> Self {
        let (n_cells, n_edges) = (mesh.n_cells(), mesh.n_edges());
        Self {
            params,
            role: vec![EdgeRole::Wall; n_edges],
            outside_level: vec![f64::NEG_INFINITY; n_edges],
            edge_depth: vec![0.0; n_edges],
            advected: vec![0.0; n_edges],
            friction: vec![0.0; n_edges],
            explicit_part: vec![0.0; n_edges],
            weight: vec![0.0; n_edges],
            explicit_flux: vec![0.0; n_edges],
            rhs: vec![0.0; n_cells],
            wet_area: vec![0.0; n_cells],
            residual: vec![0.0; n_cells],
            correction: vec![0.0; n_cells],
            theta_flux: vec![0.0; n_edges],
            pcg: PcgWorkspace::new(n_cells),
        }
    }

    /// パラメータ。
    pub fn params(&self) -> &HydroParams {
        &self.params
    }

    /// 直前のステップの θ 平均の単位幅流量 q_j^θ \[m²/s\]（浮遊砂の移流に使う）。
    ///
    /// 向きは u_n と同じ（`cells_on_edge[e][0]` から外向きを正）。連続式は
    /// V_i^{n+1} = V_i^n − Δt Σ_j s_ij L_j q_j^θ を、Newton 法の残差の範囲で満たす。
    pub fn theta_flux(&self) -> &[f64] {
        &self.theta_flux
    }

    /// 時刻を `dt` だけ進め、`state` の η, u_n, u_cell を更新する。
    ///
    /// `state` は、全セルで η ≥ z_b、壁の辺で u_n = 0、`u_cell` が `u_n` を再構成した値で
    /// あることを前提にする（ステップの最後にこの状態に揃えて返す）。失敗したときの
    /// `state` は途中の値になる。
    pub fn step(
        &mut self,
        mesh: &Mesh,
        state: &mut State,
        boundary: &BoundaryConfig,
        dt: f64,
    ) -> Result<StepReport, HydroError> {
        let Self {
            params,
            role,
            outside_level,
            edge_depth,
            advected,
            friction,
            explicit_part,
            weight,
            explicit_flux,
            rhs,
            wet_area,
            residual,
            correction,
            theta_flux,
            pcg,
        } = self;
        let (g, theta) = (params.gravity, params.theta);
        let n2 = params.manning_n * params.manning_n;

        // 1. 辺の役割と辺の水深
        for (e, (r, level)) in role.iter_mut().zip(outside_level.iter_mut()).enumerate() {
            *r = edge_role(mesh, boundary, e);
            *level = f64::NEG_INFINITY;
        }
        edge_depths(mesh, state, outside_level, params.dry_depth, edge_depth);

        // 2. 移流（次の段階で実装する）
        advected.copy_from_slice(&state.u_n);

        // 3. 摩擦係数、G_j、重み w_j、右辺に入る流量
        for e in 0..mesh.n_edges() {
            let h = edge_depth[e];
            let (mut gamma, mut big_g, mut w, mut flux) = (0.0, 0.0, 0.0, 0.0);
            if h > 0.0 {
                if let Some((eta1, eta2, delta)) = gradient_points(mesh, role[e], &state.eta, e) {
                    gamma = g * n2 * edge_speed(mesh, state, e) / h.powf(4.0 / 3.0);
                    let denom = 1.0 + dt * gamma;
                    big_g = (advected[e] - (1.0 - theta) * g * dt * (eta2 - eta1) / delta) / denom;
                    w = g * theta * theta * dt * dt * mesh.dv_edge[e] * h / (delta * denom);
                    flux = mesh.dv_edge[e] * h * (theta * big_g + (1.0 - theta) * state.u_n[e]);
                }
            }
            friction[e] = gamma;
            explicit_part[e] = big_g;
            weight[e] = w;
            explicit_flux[e] = flux;
        }
        for (c, b) in rhs.iter_mut().enumerate() {
            let outflow: f64 = mesh.edges_on_cell[c]
                .iter()
                .zip(&mesh.edge_sign_on_cell[c])
                .map(|(&e, &s)| s * explicit_flux[e])
                .sum();
            *b = volume(mesh.area_cell[c], state.eta[c], state.zb[c]) - dt * outflow;
        }

        // 4. Newton 法
        let solver = &params.solver;
        let settings = PcgSettings::from_config(solver);
        let mut report = StepReport::default();
        loop {
            report.residual = water_level_residual(
                mesh,
                &state.zb,
                &state.eta,
                weight,
                outside_level,
                rhs,
                residual,
            );
            if report.residual <= solver.newton_tol {
                break;
            }
            if report.newton_iterations >= solver.newton_max_iter {
                return Err(HydroError::NewtonNotConverged {
                    iterations: report.newton_iterations,
                    residual: report.residual,
                });
            }
            for (c, p) in wet_area.iter_mut().enumerate() {
                *p = if state.eta[c] >= state.zb[c] {
                    mesh.area_cell[c]
                } else {
                    0.0
                };
            }
            correction.fill(0.0);
            let operator = WaterLevelOperator {
                mesh,
                wet_area,
                weight,
            };
            let stats = linsolve::pcg(&operator, residual, correction, &settings, pcg)?;
            report.newton_iterations += 1;
            report.pcg_iterations += stats.iterations;
            for (eta, d) in state.eta.iter_mut().zip(correction.iter()) {
                *eta += d;
            }
        }

        // 5. 流速と θ 平均の流量
        for e in 0..mesh.n_edges() {
            let h = edge_depth[e];
            let (mut u_new, mut q) = (0.0, 0.0);
            if h > 0.0 {
                if let Some((eta1, eta2, delta)) = gradient_points(mesh, role[e], &state.eta, e) {
                    u_new = explicit_part[e]
                        - theta * g * dt * (eta2 - eta1) / (delta * (1.0 + dt * friction[e]));
                    q = h * (theta * u_new + (1.0 - theta) * state.u_n[e]);
                }
            }
            state.u_n[e] = u_new;
            theta_flux[e] = q;
        }

        // 6. 乾いたセルの水位と、セル中心の流速
        for (eta, &zb) in state.eta.iter_mut().zip(&state.zb) {
            *eta = eta.max(zb);
        }
        reconstruct_cell_velocity(mesh, &state.u_n, &mut state.u_cell);
        state.time += dt;
        Ok(report)
    }
}

/// セルの水の体積 V = A max(η − z_b, 0) \[m³\]。
fn volume(area: f64, eta: f64, zb: f64) -> f64 {
    area * (eta - zb).max(0.0)
}

/// 辺 `e` での流速の大きさ |u|_j \[m/s\]。
///
/// 法線成分は u_j、接線成分は両側のセルの Perot 流速の平均（境界辺では内側のセルの値）の
/// 接線方向の成分とする。
fn edge_speed(mesh: &Mesh, state: &State, e: usize) -> f64 {
    let [c1, c2] = mesh.cells_on_edge[e];
    let [mut vx, mut vy] = state.u_cell[c1];
    if c2 != NONE {
        let [wx, wy] = state.u_cell[c2];
        vx = 0.5 * (vx + wx);
        vy = 0.5 * (vy + wy);
    }
    let [nx, ny] = mesh.normal_edge[e];
    let tangential = vy * nx - vx * ny;
    state.u_n[e].hypot(tangential)
}

/// 水位方程式の残差 r_i = b_i − V_i(η_i) − Σ_j w_j (η_i − η_{k(i,j)}) を `out` に書き、
/// max_i |r_i| / A_i を返す（NaN があれば NaN）。境界辺の η_k は `outside_level`。
fn water_level_residual(
    mesh: &Mesh,
    zb: &[f64],
    eta: &[f64],
    weight: &[f64],
    outside_level: &[f64],
    rhs: &[f64],
    out: &mut [f64],
) -> f64 {
    let mut max = 0.0_f64;
    for (c, r) in out.iter_mut().enumerate() {
        let mut lhs = volume(mesh.area_cell[c], eta[c], zb[c]);
        for (&e, &k) in mesh.edges_on_cell[c].iter().zip(&mesh.cells_on_cell[c]) {
            let w = weight[e];
            // 重み 0 の辺は飛ばす（壁の外側の水位 −∞ を掛けて NaN にしないため）
            if w != 0.0 {
                let eta_k = if k == NONE { outside_level[e] } else { eta[k] };
                lhs += w * (eta[c] - eta_k);
            }
        }
        *r = rhs[c] - lhs;
        let scaled = r.abs() / mesh.area_cell[c];
        max = if max.is_nan() || scaled.is_nan() {
            f64::NAN
        } else {
            max.max(scaled)
        };
    }
    max
}

/// 辺の水深 H_j \[m\] を計算する。
///
/// H_j = max(η_up − max(z_{b,1}, z_{b,2}), 0)。η_up は u_j の向きで決まる風上側のセルの水位
/// （u_j = 0 のときは高い方）。`dry_depth` 未満の辺は乾いているとみなして 0 にする。
///
/// 境界辺では、外側の水位 `outside_level[e]`（水位境界の水位。外に水がなければ −∞）を
/// 隣のセルの水位とし、辺の河床高には内側のセルの値を使う。`outside_level` と `out` の長さは
/// 辺の数でなければならない。
///
/// 静水（η が一定、u = 0）では、乾いたセル（η = z_b > η_0）と接する辺の H_j は
/// z_b − z_b = 0 になるので、wet/dry の境界で偽の流れは生じない（検証ケース V03）。
pub fn edge_depths(
    mesh: &Mesh,
    state: &State,
    outside_level: &[f64],
    dry_depth: f64,
    out: &mut [f64],
) {
    debug_assert_eq!(
        outside_level.len(),
        mesh.n_edges(),
        "outside_level must have one value per edge"
    );
    debug_assert_eq!(
        out.len(),
        mesh.n_edges(),
        "out must have one value per edge"
    );
    for (e, h) in out.iter_mut().enumerate() {
        let [c1, c2] = mesh.cells_on_edge[e];
        let (eta1, zb1) = (state.eta[c1], state.zb[c1]);
        let (eta2, zb_edge) = if c2 == NONE {
            (outside_level[e], zb1)
        } else {
            (state.eta[c2], zb1.max(state.zb[c2]))
        };
        let u = state.u_n[e];
        let eta_up = if u > 0.0 {
            eta1
        } else if u < 0.0 {
            eta2
        } else {
            eta1.max(eta2)
        };
        let depth = eta_up - zb_edge;
        *h = if depth >= dry_depth { depth } else { 0.0 };
    }
}

/// Perot 法でセル中心の流速ベクトルを再構成する。
///
/// ```text
/// u_i = (1 / A_i) Σ_{j ∈ ∂i} s_ij L_j u_j (x_j − x_i)
/// ```
///
/// x_j は辺の中点、x_i は母点。発散定理から、一様流なら丸め誤差の範囲で厳密に元の
/// 流速に一致する（検証ケース V02）。境界辺も含めてセルを囲むすべての辺を使うので、
/// 境界辺の u_j（壁なら 0、流入・流出端なら境界条件で決まる値）を入れてから呼ぶ。
///
/// `u_n` の長さは辺の数、`out` の長さはセル数でなければならない。
pub fn reconstruct_cell_velocity(mesh: &Mesh, u_n: &[f64], out: &mut [[f64; 2]]) {
    debug_assert_eq!(
        u_n.len(),
        mesh.n_edges(),
        "u_n must have one value per edge"
    );
    debug_assert_eq!(
        out.len(),
        mesh.n_cells(),
        "out must have one value per cell"
    );
    for (c, u) in out.iter_mut().enumerate() {
        let (xc, yc) = (mesh.x_cell[c], mesh.y_cell[c]);
        let (mut sx, mut sy) = (0.0, 0.0);
        for (&e, &s) in mesh.edges_on_cell[c].iter().zip(&mesh.edge_sign_on_cell[c]) {
            let flux = s * mesh.dv_edge[e] * u_n[e];
            sx += flux * (mesh.x_edge[e] - xc);
            sy += flux * (mesh.y_edge[e] - yc);
        }
        *u = [sx / mesh.area_cell[c], sy / mesh.area_cell[c]];
    }
}

/// 移流項を保存型風上差分で陽的に評価し、F u_j = u_j − Δt (移流加速度)_j を `out` に書く。
///
/// セル i の移流加速度を、流入側の辺だけの寄与
/// a_i = (1 / V_i) Σ_{j: s_ij q_j < 0} s_ij L_j q_j (u_{k(i,j)} − u_i) で求め、
/// 辺の両側のセルの値を重み付き平均して法線成分をとる（Kramer & Stelling, 2008）。
/// 減速流では運動量保存形、加速流ではエネルギー水頭保存形に切り替える拡張を予定している。
///
/// 次の段階（検証ケース V05, V07）で実装する。それまで [`HydroSolver::step`] はこの関数を
/// 呼ばず、F u_j = u_j とする。
pub fn advection(mesh: &Mesh, state: &State, edge_depth: &[f64], dt: f64, out: &mut [f64]) {
    let _ = (mesh, state, edge_depth, dt, out);
    todo!("保存型風上差分による移流")
}

/// Newton 法の各段の係数行列 P + T（行列フリー）。
///
/// (P + T) x の第 i 成分は P_i x_i + Σ_{j ∈ ∂i} w_j (x_i − x_{k(i,j)})。境界辺（k = [`NONE`]）は、
/// 水位境界なら外側の水位を既知として右辺へ移すので、ここでは w_j x_i だけを残す
/// （壁と流量境界の辺は w_j = 0）。T は対称で、P + T は対称正定値になる。
#[derive(Debug, Clone, Copy)]
pub struct WaterLevelOperator<'a> {
    /// 格子。
    pub mesh: &'a Mesh,
    /// P の対角成分：wet セルでは A_i、dry セルでは 0 \[m²\]。
    pub wet_area: &'a [f64],
    /// 辺の重み w_j。
    pub weight: &'a [f64],
}

impl LinearOperator for WaterLevelOperator<'_> {
    fn dim(&self) -> usize {
        self.mesh.n_cells()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let mesh = self.mesh;
        for (c, yc) in y.iter_mut().enumerate() {
            let mut sum = self.wet_area[c] * x[c];
            for (&e, &k) in mesh.edges_on_cell[c].iter().zip(&mesh.cells_on_cell[c]) {
                let w = self.weight[e];
                sum += if k == NONE {
                    w * x[c]
                } else {
                    w * (x[c] - x[k])
                };
            }
            *yc = sum;
        }
    }

    fn diagonal(&self, d: &mut [f64]) {
        let mesh = self.mesh;
        for (c, dc) in d.iter_mut().enumerate() {
            let edges: f64 = mesh.edges_on_cell[c].iter().map(|&e| self.weight[e]).sum();
            *dc = self.wet_area[c] + edges;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::ChannelLayout;
    use crate::mesh::channel::{self, ChannelSpec};

    fn mesh(length: f64, width: f64, cell_size: f64, layout: ChannelLayout) -> Mesh {
        channel::build(&ChannelSpec {
            length,
            width,
            cell_size,
            layout,
        })
        .expect("valid channel")
    }

    fn params(manning_n: f64) -> HydroParams {
        HydroParams {
            gravity: 9.81,
            theta: 0.55,
            manning_n,
            dry_depth: 1e-3,
            solver: SolverConfig::default(),
        }
    }

    fn closed() -> BoundaryConfig {
        BoundaryConfig {
            upstream: UpstreamBoundary::Wall,
            downstream: DownstreamBoundary::Wall,
        }
    }

    fn total_volume(mesh: &Mesh, state: &State) -> f64 {
        (0..mesh.n_cells())
            .map(|c| volume(mesh.area_cell[c], state.eta[c], state.zb[c]))
            .sum()
    }

    /// 2セルの水路（正方形配置、セル 0 が x の小さい側）の内部の辺。
    fn two_cells() -> (Mesh, usize) {
        let mesh = mesh(2.0, 1.0, 1.0, ChannelLayout::Rectangular);
        assert_eq!(mesh.n_cells(), 2);
        let e = (0..mesh.n_edges())
            .find(|&e| !mesh.is_boundary_edge(e))
            .expect("one interior edge");
        (mesh, e)
    }

    #[test]
    fn edge_depth_uses_the_upwind_level_and_the_higher_bed() {
        let (mesh, e) = two_cells();
        let [c1, c2] = mesh.cells_on_edge[e];
        let outside = vec![f64::NEG_INFINITY; mesh.n_edges()];
        let mut state = State::new(2, mesh.n_edges());
        state.eta[c1] = 1.0;
        state.eta[c2] = 0.8;
        state.zb[c1] = 0.2;
        state.zb[c2] = 0.5;
        let mut h = vec![0.0; mesh.n_edges()];

        let mut depth_for = |u: f64| {
            state.u_n[e] = u;
            edge_depths(&mesh, &state, &outside, 1e-3, &mut h);
            h[e]
        };
        // 河床高は高い方（0.5）、水位は風上側
        assert!((depth_for(0.3) - 0.5).abs() < 1e-15);
        assert!((depth_for(-0.3) - 0.3).abs() < 1e-15);
        assert!((depth_for(0.0) - 0.5).abs() < 1e-15);
    }

    #[test]
    fn shallow_and_dry_edges_have_zero_depth() {
        let (mesh, e) = two_cells();
        let [c1, c2] = mesh.cells_on_edge[e];
        let outside = vec![f64::NEG_INFINITY; mesh.n_edges()];
        let mut state = State::new(2, mesh.n_edges());
        let mut h = vec![0.0; mesh.n_edges()];

        // 湛水セル（η = 1）と、水面より高い乾いたセル（η = z_b = 1.2）
        state.eta = vec![1.0, 1.2];
        state.zb = vec![0.0, 1.2];
        edge_depths(&mesh, &state, &outside, 1e-3, &mut h);
        assert_eq!(h[e], 0.0);

        // 閾値より浅い
        state.eta[c1] = 0.5005;
        state.eta[c2] = 0.5;
        state.zb[c1] = 0.0;
        state.zb[c2] = 0.5;
        edge_depths(&mesh, &state, &outside, 1e-3, &mut h);
        assert_eq!(h[e], 0.0);
        edge_depths(&mesh, &state, &outside, 1e-4, &mut h);
        assert!((h[e] - 5e-4).abs() < 1e-12);
    }

    #[test]
    fn water_level_operator_is_symmetric_and_consistent_with_its_diagonal() {
        let mesh = mesh(10.0, 4.0, 1.0, ChannelLayout::Hexagonal);
        let n = mesh.n_cells();
        let wet_area: Vec<f64> = (0..n)
            .map(|c| if c % 5 == 0 { 0.0 } else { mesh.area_cell[c] })
            .collect();
        // 境界辺の一部にも重みを持たせる（水位境界の形）
        let weight: Vec<f64> = (0..mesh.n_edges())
            .map(|e| 0.1 + (e % 4) as f64 * 0.3)
            .collect();
        let op = WaterLevelOperator {
            mesh: &mesh,
            wet_area: &wet_area,
            weight: &weight,
        };
        let x: Vec<f64> = (0..n).map(|i| (0.7 * i as f64).sin()).collect();
        let y: Vec<f64> = (0..n).map(|i| (1.3 * i as f64).cos()).collect();
        let (mut ax, mut ay) = (vec![0.0; n], vec![0.0; n]);
        op.apply(&x, &mut ax);
        op.apply(&y, &mut ay);
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(p, q)| p * q).sum::<f64>();
        let (xay, yax) = (dot(&x, &ay), dot(&y, &ax));
        assert!(
            (xay - yax).abs() <= 1e-12 * xay.abs().max(1.0),
            "{xay} != {yax}"
        );

        let mut d = vec![0.0; n];
        op.diagonal(&mut d);
        for c in [0, 1, n / 2, n - 1] {
            let mut unit = vec![0.0; n];
            unit[c] = 1.0;
            let mut column = vec![0.0; n];
            op.apply(&unit, &mut column);
            assert!((column[c] - d[c]).abs() <= 1e-14 * d[c]);
        }
    }

    #[test]
    fn lake_at_rest_stays_at_rest() {
        let mesh = mesh(20.0, 10.0, 1.0, ChannelLayout::Hexagonal);
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        for c in 0..mesh.n_cells() {
            let (x, y) = (mesh.x_cell[c] - 10.0, mesh.y_cell[c] - 5.0);
            state.zb[c] = 2.0 * (-(x * x + y * y) / 8.0).exp();
            state.eta[c] = state.zb[c].max(1.0);
        }
        let initial = state.clone();
        let mut solver = HydroSolver::new(&mesh, params(0.03));
        for _ in 0..10 {
            let report = solver.step(&mesh, &mut state, &closed(), 0.5).unwrap();
            assert_eq!(report.newton_iterations, 0);
            assert_eq!(report.residual, 0.0);
        }
        assert_eq!(state.eta, initial.eta);
        assert!(state.u_n.iter().all(|&u| u == 0.0));
    }

    /// 1ステップの計算を `steps` 回繰り返し、体積の保存、Newton 法の収束、η ≥ z_b を確かめる。
    /// 戻り値は Newton 法の反復回数の最大値。
    fn run_closed(mesh: &Mesh, state: &mut State, manning_n: f64, dt: f64, steps: usize) -> usize {
        let v0 = total_volume(mesh, state);
        let mut solver = HydroSolver::new(mesh, params(manning_n));
        let tol = solver.params().solver.newton_tol;
        let mut max_newton = 0;
        for _ in 0..steps {
            let report = solver.step(mesh, state, &closed(), dt).unwrap();
            assert!(report.residual <= tol, "residual {:e}", report.residual);
            max_newton = max_newton.max(report.newton_iterations);
        }
        let v = total_volume(mesh, state);
        assert!(((v - v0) / v0).abs() <= 1e-12, "volume {v} != {v0}");
        assert!((0..mesh.n_cells()).all(|c| state.eta[c] >= state.zb[c]));
        assert!(state.u_n.iter().any(|&u| u != 0.0));
        max_newton
    }

    /// 傾いた水面から振動させ、斜面の冠水と干出を含めて体積が保存されることを確かめる。
    #[test]
    fn sloshing_with_wetting_and_drying_conserves_volume() {
        let mesh = mesh(20.0, 4.0, 1.0, ChannelLayout::Hexagonal);
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        for c in 0..mesh.n_cells() {
            let x = mesh.x_cell[c];
            // 両端が高い谷形の河床。左側だけ水位を上げる。
            state.zb[c] = 0.02 * (x - 10.0).powi(2);
            let level = if x < 10.0 { 1.0 } else { 0.5 };
            state.eta[c] = state.zb[c].max(level);
        }
        run_closed(&mesh, &mut state, 0.02, 0.2, 200);
    }

    /// 浅い棚の水が1ステップで流れ落ちる（解でセルが干上がる）ときは、Newton 法の反復が
    /// 2回以上必要になる。その場合も体積が保存されることを確かめる。
    #[test]
    fn draining_a_shelf_needs_newton_iterations_and_conserves_volume() {
        let mesh = mesh(20.0, 2.0, 1.0, ChannelLayout::Rectangular);
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        for c in 0..mesh.n_cells() {
            let shelf = mesh.x_cell[c] < 10.0;
            state.zb[c] = if shelf { 0.5 } else { 0.0 };
            state.eta[c] = if shelf { 0.52 } else { 0.2 };
        }
        let max_newton = run_closed(&mesh, &mut state, 0.02, 2.0, 20);
        assert!(max_newton >= 2, "max Newton iterations {max_newton}");
    }
}
