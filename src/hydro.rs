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
//! γ_j = g n² |u|_j / H_j^{4/3}
//! ```
//!
//! これを連続式に代入すると、水位について次の連立方程式が得られる。
//!
//! ```text
//! V_i(η_i) + Σ_j w_j (η_i − η_{k(i,j)}) = b_i
//! w_j = g θ² Δt² L_j H_j / (δ_j (1 + Δt γ_j))
//! ```
//!
//! V_i(η) = A_i max(η − z_{b,i}, 0) は wet/dry の境目（η = z_b）で折れる区分線形関数で、
//! 係数行列は対称正定値になる。これを Newton 法（Casulli, 2009; Casulli & Zanolli, 2012）で
//! 解き、各段の線形系を行列フリーの PCG（[`crate::linsolve::pcg`]）で解く。
//!
//! ## 1ステップの手順
//!
//! 1. 辺の水深 H_j^n と wet/dry の判定（[`edge_depths`]）。H_j は風上側の水位から辺の河床高
//!    max(z_{b,1}, z_{b,2}) を引いた値（Stelling & Duinmeijer, 2003）。
//! 2. Perot 法でセル中心の流速を再構成する（[`reconstruct_cell_velocity`]）。
//! 3. 移流を保存型風上差分で陽的に評価し、F u_j^n を求める（[`advection`]）。
//! 4. 摩擦係数 γ_j、G_j、重み w_j、右辺 b_i を組み立てる。
//! 5. Newton 法＋行列フリー PCG で η^{n+1} を求める（[`WaterLevelOperator`]）。
//! 6. u_j^{n+1} を更新し、θ 平均の単位幅流量 q_j^θ = H_j^n (θ u_j^{n+1} + (1 − θ) u_j^n) を
//!    保存する。浮遊砂の移流に同じ流量を使い、濃度が一様な状態が保たれるようにする。
//!
//! ## 参考文献
//!
//! - Casulli & Walters (2000): 非構造格子上の半陰解法
//! - Casulli (2009), Casulli & Zanolli (2012): wet/dry と Newton 法
//! - Stelling & Duinmeijer (2003), Kramer & Stelling (2008): 保存型風上差分
//! - Kernkamp et al. (2011): 非構造格子への適用（D-Flow FM）
//! - Perot (2000): セル中心流速の再構成

use crate::config::{BoundaryConfig, Config, SolverConfig};
use crate::linsolve::{LinearOperator, PcgError, PcgWorkspace};
use crate::mesh::Mesh;
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
    /// Newton 法の反復回数。
    pub newton_iterations: usize,
    /// PCG の反復回数の合計。
    pub pcg_iterations: usize,
    /// 最終的な残差（体積残差をセル面積で割った値の最大値）\[m\]。
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

/// 半陰解法のソルバー。作業配列を持ち回して、ステップごとの確保を避ける。
#[derive(Debug, Clone)]
#[expect(dead_code, reason = "v0.1.0 の雛形。step の実装で使う")]
pub struct HydroSolver {
    params: HydroParams,
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
    /// 水位方程式の右辺 b_i \[m³\]
    rhs: Vec<f64>,
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
            edge_depth: vec![0.0; n_edges],
            advected: vec![0.0; n_edges],
            friction: vec![0.0; n_edges],
            explicit_part: vec![0.0; n_edges],
            weight: vec![0.0; n_edges],
            rhs: vec![0.0; n_cells],
            theta_flux: vec![0.0; n_edges],
            pcg: PcgWorkspace::new(n_cells),
        }
    }

    /// パラメータ。
    pub fn params(&self) -> &HydroParams {
        &self.params
    }

    /// 直前のステップの θ 平均の単位幅流量 q_j^θ \[m²/s\]（浮遊砂の移流に使う）。
    pub fn theta_flux(&self) -> &[f64] {
        &self.theta_flux
    }

    /// 時刻を `dt` だけ進め、`state` の η, u_n, u_cell を更新する。
    pub fn step(
        &mut self,
        mesh: &Mesh,
        state: &mut State,
        boundary: &BoundaryConfig,
        dt: f64,
    ) -> Result<StepReport, HydroError> {
        let _ = (mesh, state, boundary, dt);
        todo!("Casulli 型 θ 法の1ステップ")
    }
}

/// 辺の水深 H_j \[m\] を計算する。
///
/// H_j = max(η_up − max(z_{b,1}, z_{b,2}), 0)。η_up は u_j の向きで決まる風上側のセルの水位
/// （u_j = 0 のときは高い方）。`dry_depth` 未満の辺は乾いているとみなして 0 にする。
pub fn edge_depths(mesh: &Mesh, state: &State, dry_depth: f64, out: &mut [f64]) {
    let _ = (mesh, state, dry_depth, out);
    todo!("辺の水深（風上水位 − 辺の河床高）")
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
pub fn advection(mesh: &Mesh, state: &State, edge_depth: &[f64], dt: f64, out: &mut [f64]) {
    let _ = (mesh, state, edge_depth, dt, out);
    todo!("保存型風上差分による移流")
}

/// Newton 法の各段の係数行列 P + T（行列フリー）。
///
/// (P + T) x の第 i 成分は P_i x_i + Σ_{j ∈ ∂i} w_j (x_i − x_{k(i,j)})。境界辺は、水位境界なら
/// 隣を既知の水位として右辺へ移し、ここでは w_j x_i だけを残す。
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
        let _ = (x, y);
        todo!("(P + T) x のセル単位のループ")
    }

    fn diagonal(&self, d: &mut [f64]) {
        let _ = d;
        todo!("P + T の対角成分")
    }
}
