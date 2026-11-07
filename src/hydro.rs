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
//! G_j = [F u_j^n − (1 − θ) g Δt / δ_j · (η_{2(j)}^n − η_{1(j)}^n) + Δt g I n_{x,j}] / (1 + Δt γ_j)
//! γ_j = g n² |u|_j / H_j^{4/3}
//! ```
//!
//! g I n_{x,j} は周期境界の体積力で、I は河床高と水位から除いた縦断勾配
//! （[`HydroParams::background_slope`]。周期境界でなければ 0。ADR-0014）。
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
//! ## 境界条件
//!
//! - 壁：u_j = 0、w_j = 0。
//! - 上流端（流量）：Q を通水能 K_j = L_j max(η_s − z_{b,i}, `dry_depth`)^{5/3} / n に比例して
//!   配分し、既知の流入 Q_j として右辺 b_i に入れる（w_j = 0）。η_s は上流端の湿ったセルの
//!   水位の平均（断面で共通の水位。ADR-0007）、z_{b,i} は内側のセルの河床高。境界辺の u_j は、
//!   単位幅流量 Q_j / L_j を新しい水深で割った値にする（Perot 再構成と出力に使う）。
//! - 周期境界：上流端と下流端を継いだ格子（[`crate::mesh::channel::build_periodic`]）では、
//!   継ぎ目の辺は内部の辺として扱い、開境界の辺はない。河床高と水位は縦断勾配 I を除いた値で
//!   持ち、勾配は体積力 g I n_x として運動量式に入れる。辺の水深の風上と高い方の河床は、
//!   除いた勾配を戻した実際の高さで選ぶ（[`edge_depths`]、ADR-0015）。Perot 法では、継ぎ目の
//!   辺の中点をセルの側の座標（[`Mesh::edge_midpoint_from`]）で使う。
//! - 下流端（水位・水深）：境界辺の上の水位 η_b を既知とし、母点から辺までの距離
//!   δ_j = `dc_edge` / 2 で水位勾配を評価する（ADR-0007）。w_j η_b を右辺へ移す。水深一定の
//!   境界では η_b = z_{b,i} + 水深（z_{b,i} は内側のセルの河床高）。辺の水深の風上は、流出なら
//!   内側のセル、流入なら η_b、u_j = 0 なら高い方の水位で決め、河床高は内側のセルの値を使う。
//!
//! ## 移流
//!
//! 各辺の流速を、両側のセルに流れ込む水の流速（法線成分）へ、局所クーラン数 σ_j の割合で
//! 近づける（ADR-0008, ADR-0009）。
//!
//! ```text
//! S_i   = Σ_{l ∈ ∂i, s_il q_l < 0} |L_l q_l|,  P_i = (1 / A_i) Σ_{同じ辺} |L_l q_l| u_{k(i,l)}
//! h*_i  = h_i − (Δt / A_i) Σ_l s_il L_l q_l
//! σ_j   = Δt (S_1 / A_1 + S_2 / A_2) / (h*_1 + h*_2),  ū_j = (P_1 + P_2)·n_j / (S_1 / A_1 + S_2 / A_2)
//! F u_j = u_j + min(σ_j, 1) (ū_j − u_j)
//! ```
//!
//! q_l = H_l^n u_l^n、u_k は Perot 法で再構成したセル中心の流速、h*_i は時刻 n の流量で連続式を
//! 陽的に進めたセルの水深。保存形の移流 Σ s L q u_up − u_j Σ s L q を流入側の辺で評価したもので、
//! 辺は母点を結ぶ線分の垂直二等分線なので、両側のセルの重みは 1/2 ずつになる。
//!
//! - 1次元の定常流では、水位勾配の項と合わせて q u + g h²/2 が離散的に保存され、跳水の
//!   共役水深が Bélanger の式を満たす（Stelling & Duinmeijer, 2003）。
//! - 辺の流速そのものを上流の流速へ近づけるので、六角形配置で生じる、セル中心の流速に
//!   現れない格子規模の辺の流速の成分も下流へ運び去られる。
//! - σ_j ≤ 1 なら F u_j は u_j と ū_j の凸結合で、汀線の浅いセルや段波の先端でも発散しない。
//! - 境界辺から流れ込む水の流速は、その辺の法線流速 u_j n_j とする（流量境界の流入）。
//!   水位境界の辺はセル 1(j) の値だけで評価する。
//!
//! `hydro.advection = "second_order"` では、セル中心の流速の勾配（最小二乗。
//! [`crate::mesh::gradient`]）で面の値と辺へ外挿したセルの流速を作り、空間について2次精度に
//! する（[`advection_second_order`]、ADR-0016）。勾配を 0 にすると上の形に一致する。汀線の
//! 近くと、開境界の辺を持つセルでは勾配を 0 にする（ADR-0018）。
//!
//! ## 時間積分
//!
//! 下の「1ステップの手順」の 1〜2 を陽的な部分（辺の水深 H_j と流量、移流の増分
//! Δ_j = F u_j − u_j、摩擦係数に使う |u|_j）、3〜6 を陰的な部分とする。
//!
//! - `hydro.time_integration = "euler"`：時刻 n の陽的な部分で1回解く（移流は前進 Euler）。
//! - `hydro.time_integration = "predictor_corrector"`：予測子・修正子（ADR-0017）。時刻 n の
//!   陽的な部分で1回解いて予測子 U* を求め、U* の陽的な部分と平均して（H̄ = (H^n + H*)/2、
//!   移流は辺の運動量と体積を2つの状態の流束の平均で進めた値、流量と |u| も平均）、時刻 n から
//!   解き直す。陽的な部分は Heun 法、
//!   水位勾配と摩擦は θ 法の1回なので、平均流のない重力波の減衰は θ 法と同じ。計算量は約2倍だが、
//!   1次元の線形解析ではどのフルード数でもクーラン数 1 まで安定で（`tools/advection_stability.py`）、
//!   定常解は1段の定常解と同じ（Δt によらない）。2次の移流は1段では不安定なので、この形と
//!   組み合わせる。
//!
//! ## 1ステップの手順
//!
//! 1. 辺の役割（内部・壁・流量境界・水位境界）と辺の水深 H_j^n を決める（[`edge_depths`]）。
//!    H_j は風上側の水位から辺の河床高 max(z_{b,1}, z_{b,2}) を引いた値（Stelling & Duinmeijer,
//!    2003）。上流端の流量 Q を、通水能に比例して境界辺に配分する（上の「境界条件」）。
//!    辺の単位幅流量 q_j^n = H_j^n u_j^n（流量境界では配分した流入量）を求める。
//!    この流量で、セルの流出量の比（クーラン数）Δt max_i Σ_流出 L q / V を記録する
//!    （[`HydroSolver::max_outflow_rate`]。可変の時間刻みはこの比から決める。ADR-0010）。
//! 2. 移流を保存型の風上差分で陽的に評価し、F u_j^n を求める（[`advection`]、
//!    [`advection_second_order`]、上の「移流」）。`hydro.advection = "none"` のときは
//!    F u_j^n = u_j^n とする。
//! 3. 摩擦係数 γ_j、G_j、重み w_j、右辺 b_i を組み立てる。
//! 4. Newton 法＋行列フリー PCG で η^{n+1} を求める（[`WaterLevelOperator`]）。
//! 5. u_j^{n+1} を更新し、θ 平均の単位幅流量 q_j^θ = H_j^n (θ u_j^{n+1} + (1 − θ) u_j^n) を
//!    保存する。浮遊砂の移流に同じ流量を使い、濃度が一様な状態が保たれるようにする。
//! 6. 乾いたセルの水位を河床高に揃え、Perot 法でセル中心の流速を再構成する
//!    （[`reconstruct_cell_velocity`]）。
//!
//! 予測子・修正子では、1〜6 で予測子を求め、予測子の状態で 1〜2 を行って時刻 n の値と平均し、
//! 時刻 n の状態から 3〜6 を行う（上の「時間積分」）。
//!
//! ## 参考文献
//!
//! - Casulli & Walters (2000): 非構造格子上の半陰解法
//! - Casulli (2009), Casulli & Zanolli (2012): wet/dry と Newton 法
//! - Stelling & Duinmeijer (2003), Kramer & Stelling (2008): 保存型風上差分
//! - Kernkamp et al. (2011): 非構造格子への適用（D-Flow FM）
//! - Deltares (2026) D-Flow FM Technical Reference Manual §6.2: Perot 型の移流の離散化
//! - Perot (2000): セル中心流速の再構成
//! - Gottlieb, Shu & Tadmor (2001): SSP（強安定性保存）の Runge–Kutta 法

use crate::config::{
    AdvectionScheme, BoundaryConfig, Config, DownstreamBoundary, SolverConfig, TimeIntegration,
    UpstreamBoundary,
};
use crate::linsolve::{self, LinearOperator, PcgError, PcgSettings, PcgWorkspace};
use crate::mesh::gradient::CellGradient;
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
    /// 運動量の移流の扱い。
    pub advection: AdvectionScheme,
    /// 時間積分の方法。
    pub time_integration: TimeIntegration,
    /// 反復解法の設定。
    pub solver: SolverConfig,
    /// 周期境界で河床高と水位から除いた縦断勾配 I \[-\]（x の正の向きに下がる向きを正）。
    /// 運動量式に体積力 g I n_x として入れる。周期境界でなければ 0（ADR-0014）。
    pub background_slope: f64,
}

impl HydroParams {
    /// 設定から作る。
    pub fn from_config(config: &Config) -> Self {
        Self {
            gravity: config.physics.gravity,
            theta: config.hydro.theta,
            manning_n: config.hydro.manning_n,
            dry_depth: config.hydro.dry_depth,
            advection: config.hydro.advection,
            time_integration: config.hydro.time_integration,
            solver: config.hydro.solver.clone(),
            background_slope: background_slope(config),
        }
    }
}

/// 周期境界で河床高と水位から除いた縦断勾配（`[bed] slope`）。周期境界でなければ 0。
pub fn background_slope(config: &Config) -> f64 {
    match (&config.bed, config.boundary.is_periodic()) {
        (Some(bed), true) => bed.slope,
        _ => 0.0,
    }
}

/// 1ステップの計算の記録。
///
/// 予測子・修正子（`hydro.time_integration = "predictor_corrector"`）では、反復回数は2回の合計、
/// セルの流出量の比 `outflow_courant` は時刻 n の値（可変の時間刻みを決めた値）、移流の局所
/// クーラン数と残差は2回の大きい方を記録する。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct StepReport {
    /// 時間刻み Δt \[s\]。
    pub dt: f64,
    /// セルの流出量の比（クーラン数）Δt Σ_{流出する辺} L_j q_j^n / V_i の最大値 \[-\]。
    /// 1 以下なら、陽的に扱う部分で水深が負にならない（[`HydroSolver::max_outflow_rate`]）。
    pub outflow_courant: f64,
    /// Newton 法の反復回数（解いた線形系の数）。
    pub newton_iterations: usize,
    /// PCG の反復回数の合計。
    pub pcg_iterations: usize,
    /// 最終的な残差（体積残差をセル面積で割った値の最大値）\[m\]。乾いたセルを含む。
    pub residual: f64,
    /// 移流の局所クーラン数 σ_j の最大値（重みを 1 で頭打ちにする前）\[-\]。移流を省くときは 0。
    pub max_courant: f64,
    /// 移流の重みを σ_j > 1 のために 1 で頭打ちにした辺の数（時間刻みが大きすぎる兆候）。
    pub limited_edges: usize,
}

/// 移流の関数（[`advection`] と [`advection_second_order`]）の型。
type AdvectionFn =
    fn(&Mesh, &State, &[f64], &[bool], f64, &mut AdvectionWorkspace, &mut [f64]) -> AdvectionStats;

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
    /// 流量境界（流入量を与える）。
    Discharge,
    /// 水位境界。値は境界辺の上の水位 η_b \[m\]。
    Level(f64),
}

/// 辺 `e` の役割を決める。水深一定の境界では、内側のセルの河床高に水深を加えた値を η_b とする。
fn edge_role(mesh: &Mesh, state: &State, boundary: &BoundaryConfig, e: usize) -> EdgeRole {
    match mesh.edge_kind[e] {
        EdgeKind::Interior => EdgeRole::Interior,
        EdgeKind::Wall => EdgeRole::Wall,
        // 周期境界の格子には開境界の辺がない（継ぎ目の辺は内部の辺）
        EdgeKind::Upstream => match boundary.upstream {
            UpstreamBoundary::Wall | UpstreamBoundary::Periodic { .. } => EdgeRole::Wall,
            UpstreamBoundary::Discharge { .. } => EdgeRole::Discharge,
        },
        EdgeKind::Downstream => match boundary.downstream {
            DownstreamBoundary::Wall | DownstreamBoundary::Periodic => EdgeRole::Wall,
            DownstreamBoundary::WaterLevel { level } => EdgeRole::Level(level),
            DownstreamBoundary::Depth { depth } => {
                EdgeRole::Level(state.zb[mesh.cells_on_edge[e][0]] + depth)
            }
        },
    }
}

/// 水位勾配を評価する辺なら (η_{1(j)}, η_{2(j)}, δ_j) を返す。壁と流量境界なら `None`。
///
/// 水位境界では、η_2 を境界辺の上の水位 η_b、δ_j を母点から辺までの距離（`dc_edge` の半分）
/// とする（ADR-0007）。
fn gradient_points(mesh: &Mesh, role: EdgeRole, eta: &[f64], e: usize) -> Option<(f64, f64, f64)> {
    let [c1, c2] = mesh.cells_on_edge[e];
    match role {
        EdgeRole::Interior => Some((eta[c1], eta[c2], mesh.dc_edge[e])),
        EdgeRole::Level(level) => Some((eta[c1], level, 0.5 * mesh.dc_edge[e])),
        EdgeRole::Wall | EdgeRole::Discharge => None,
    }
}

/// 流量境界の辺 `e` の通水能の重み L_j max(η_s − z_{b,i}, `dry_depth`)^{5/3}。
///
/// η_s は上流端の断面の水位（[`inflow_section_level`]）、z_{b,i} は内側のセルの河床高。
/// 粗度係数は辺によらないので、配分の比 K_j / Σ K_j では打ち消し合う。水深の下限を
/// `dry_depth` とするので、上流端のセルがすべて乾いていれば辺の長さに比例して配分される。
fn conveyance(mesh: &Mesh, state: &State, section_level: f64, dry_depth: f64, e: usize) -> f64 {
    let h = (section_level - state.zb[mesh.cells_on_edge[e][0]]).max(dry_depth);
    mesh.dv_edge[e] * h.powf(5.0 / 3.0)
}

/// 上流端の断面の水位 η_s \[m\]：流量境界の辺の内側の湿ったセル（水深 ≥ `dry_depth`）の水位を、
/// 辺の長さで重み付けして平均した値。湿ったセルがなければ −∞。
///
/// 各セルの水深で通水能を計算すると、深くなったセルほど流入が増えてさらに深くなる正の
/// フィードバックが生じ、横断方向に不安定になる（ADR-0007）。断面で共通の水位を使うと、
/// 配分は河床高の違いだけで決まる。
fn inflow_section_level(mesh: &Mesh, state: &State, role: &[EdgeRole], dry_depth: f64) -> f64 {
    let (mut sum, mut length) = (0.0, 0.0);
    for (e, _) in role
        .iter()
        .enumerate()
        .filter(|&(_, &r)| r == EdgeRole::Discharge)
    {
        let c = mesh.cells_on_edge[e][0];
        if state.depth(c) >= dry_depth {
            sum += mesh.dv_edge[e] * state.eta[c];
            length += mesh.dv_edge[e];
        }
    }
    if length > 0.0 {
        sum / length
    } else {
        f64::NEG_INFINITY
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
    /// 辺の単位幅流量 q_j^n（流量境界では配分した流入量）\[m²/s\]
    edge_flux: Vec<f64>,
    /// 運動量式を解く辺（内部の辺と水位境界の辺のうち、乾いていない辺）
    momentum_edge: Vec<bool>,
    /// 移流の作業配列
    advection: AdvectionWorkspace,
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
    /// 摩擦係数に使う辺の流速の大きさ |u|_j \[m/s\]
    speed: Vec<f64>,
    /// 時刻 n の流量で連続式を陽的に進めたセルの水深 h*_i \[m\]（移流と同じ値）
    advanced_depth: Vec<f64>,
    /// 辺の両側のセルの水深の和 h_1 + h_2 \[m\]（境界辺は内側のセルだけ）
    side_depth: Vec<f64>,
    /// 辺の両側のセルの h* の和 h*_1 + h*_2 \[m\]（境界辺は内側のセルだけ）
    advanced_side_depth: Vec<f64>,
    /// 予測子・修正子の、時刻 n の水位 \[m\]
    initial_level: Vec<f64>,
    /// 予測子・修正子の、時刻 n の辺の流速 u_j^n \[m/s\]
    initial_velocity: Vec<f64>,
    /// 予測子・修正子の、時刻 n の辺の水深 H_j^n \[m\]
    initial_edge_depth: Vec<f64>,
    /// 予測子・修正子の、時刻 n の辺の単位幅流量 q_j^n \[m²/s\]（流量境界では配分した流入量）
    initial_edge_flux: Vec<f64>,
    /// 予測子・修正子の、時刻 n の |u|_j \[m/s\]
    initial_speed: Vec<f64>,
    /// 予測子・修正子の、時刻 n の移流の増分 F u_j^n − u_j^n \[m/s\]
    initial_increment: Vec<f64>,
    /// 予測子・修正子の、時刻 n の h_1 + h_2 \[m\]
    initial_side_depth: Vec<f64>,
    /// 予測子・修正子の、時刻 n の h*_1 + h*_2 \[m\]
    initial_advanced_side_depth: Vec<f64>,
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
            edge_flux: vec![0.0; n_edges],
            momentum_edge: vec![false; n_edges],
            advection: AdvectionWorkspace::new(mesh),
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
            speed: vec![0.0; n_edges],
            advanced_depth: vec![0.0; n_cells],
            side_depth: vec![0.0; n_edges],
            advanced_side_depth: vec![0.0; n_edges],
            initial_level: vec![0.0; n_cells],
            initial_velocity: vec![0.0; n_edges],
            initial_edge_depth: vec![0.0; n_edges],
            initial_edge_flux: vec![0.0; n_edges],
            initial_speed: vec![0.0; n_edges],
            initial_increment: vec![0.0; n_edges],
            initial_side_depth: vec![0.0; n_edges],
            initial_advanced_side_depth: vec![0.0; n_edges],
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
    /// 予測子・修正子では、修正子の流量（辺の水深は時刻 n と予測子の平均）。
    pub fn theta_flux(&self) -> &[f64] {
        &self.theta_flux
    }

    /// 時刻 n の状態で、セルの流出量の比 max_i Σ_{流出する辺} L_j q_j^n / V_i \[1/s\] を返す。
    ///
    /// q_j^n = H_j^n u_j^n は [`HydroSolver::step`] の最初の段と同じ値（流量境界の辺は流入なので
    /// 数えない）。これに Δt を掛けた値が、そのステップのセルの流出量の比（クーラン数）になり、
    /// 1 以下なら陽的に扱う部分（連続式の辺の水深と移流）で水深が負にならない
    /// （`docs/design/overview.md` の第7章）。可変の時間刻みは Δt = C / この値で決める
    /// （ADR-0010）。どのセルからも流出がなければ 0。
    pub fn max_outflow_rate(
        &mut self,
        mesh: &Mesh,
        state: &State,
        boundary: &BoundaryConfig,
    ) -> f64 {
        self.prepare_edges(mesh, state, boundary);
        outflow_rate(mesh, state, &self.edge_flux)
    }

    /// ステップの最初の段：辺の役割、辺の水深 H_j^n、辺の単位幅流量 q_j^n（流量境界では
    /// 配分した流入量）、運動量式を解く辺を決める。
    fn prepare_edges(&mut self, mesh: &Mesh, state: &State, boundary: &BoundaryConfig) {
        let Self {
            params,
            role,
            outside_level,
            edge_depth,
            edge_flux,
            momentum_edge,
            ..
        } = self;
        for (e, (r, level)) in role.iter_mut().zip(outside_level.iter_mut()).enumerate() {
            *r = edge_role(mesh, state, boundary, e);
            *level = match *r {
                EdgeRole::Level(eta_b) => eta_b,
                _ => f64::NEG_INFINITY,
            };
        }
        edge_depths(
            mesh,
            state,
            outside_level,
            params.dry_depth,
            params.background_slope,
            edge_depth,
        );
        let discharge = match boundary.upstream {
            UpstreamBoundary::Discharge { discharge, .. } => discharge,
            UpstreamBoundary::Wall | UpstreamBoundary::Periodic { .. } => 0.0,
        };
        let section_level = inflow_section_level(mesh, state, role, params.dry_depth);
        let total_conveyance: f64 = (0..mesh.n_edges())
            .filter(|&e| role[e] == EdgeRole::Discharge)
            .map(|e| conveyance(mesh, state, section_level, params.dry_depth, e))
            .sum();

        for (e, (q, active)) in edge_flux
            .iter_mut()
            .zip(momentum_edge.iter_mut())
            .enumerate()
        {
            let h = edge_depth[e];
            (*q, *active) = match role[e] {
                // 法線は領域の外を向くので、流入は負の流量
                EdgeRole::Discharge => (
                    -discharge * conveyance(mesh, state, section_level, params.dry_depth, e)
                        / (total_conveyance * mesh.dv_edge[e]),
                    false,
                ),
                EdgeRole::Wall => (0.0, false),
                EdgeRole::Interior | EdgeRole::Level(_) => (h * state.u_n[e], h > 0.0),
            };
        }
    }

    /// 時刻を `dt` だけ進め、`state` の η, u_n, u_cell を更新する。
    ///
    /// `hydro.time_integration` が `"euler"` なら半陰解法の1段、`"predictor_corrector"` なら
    /// 予測子・修正子（モジュールの説明の「時間積分」）で進める。
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
        match self.params.time_integration {
            TimeIntegration::Euler => {
                let report = self.explicit_terms(mesh, state, boundary, dt);
                self.implicit_solve(mesh, state, dt, report)
            }
            TimeIntegration::PredictorCorrector => {
                self.step_predictor_corrector(mesh, state, boundary, dt)
            }
        }
    }

    /// 予測子・修正子（ADR-0017）。
    ///
    /// 1. 時刻 n の陽的な部分（辺の水深 H^n、流量、移流の増分 Δ^n = F u^n − u^n、|u|^n）で
    ///    1段解き、予測子 U* を得る。
    /// 2. U* の陽的な部分（H*、流量、Δ* = F* u* − u*、|u|*）を求める。
    /// 3. 陽的な部分を平均し、時刻 n の状態から水位と流速を解き直す（修正子）。辺の水深は
    ///    H̄ = (H^n + H*)/2、流量と |u| も平均する。移流は、辺の運動量 (h_1 + h_2) u と体積
    ///    h_1 + h_2 を、2つの状態の流束の平均で陽的に進めた値とする（ADR-0008 の h* で割る形を
    ///    2つの状態で平均したもの）。
    ///
    /// ```text
    /// F u = [h^n u^n + (P^n + P*)/2] / [h^n − (D^n + D*)/2]
    /// P_s = h*_s F_s u_s − h_s u_s（運動量の流束 × Δt）,  D_s = h_s − h*_s（体積の流出 × Δt）
    /// ```
    ///
    /// h_s、h*_s は状態 s の両側のセルの水深と h* の和。流速の増分 F_s u_s − u_s をそのまま
    /// 平均すると、動く段波で運動量が保存されず、段波の速さがずれる（V05）。どちらかの状態で
    /// h* が 0 の辺（干出しかけた辺）では、流速の増分の平均にする。
    ///
    /// 陽的な部分については Heun 法（SSP-RK2 と同じ）、水位勾配と摩擦については θ 法の1回なので、
    /// 平均流のない重力波の減衰は θ 法と同じになる。連続式は修正子の1回の流量で閉じるので、
    /// V^{n+1} = V^n − Δt Σ_j s_ij L_j q_j^θ が Newton 法の残差の範囲で成り立つ。
    /// 定常解は1段の定常解と同じ（U* = U^n なら修正子は予測子と同じ式になる）。
    fn step_predictor_corrector(
        &mut self,
        mesh: &Mesh,
        state: &mut State,
        boundary: &BoundaryConfig,
        dt: f64,
    ) -> Result<StepReport, HydroError> {
        let time = state.time;
        // 1. 時刻 n の陽的な部分と予測子
        let predictor = self.explicit_terms(mesh, state, boundary, dt);
        self.initial_level.copy_from_slice(&state.eta);
        self.initial_velocity.copy_from_slice(&state.u_n);
        self.initial_edge_depth.copy_from_slice(&self.edge_depth);
        self.initial_edge_flux.copy_from_slice(&self.edge_flux);
        self.initial_speed.copy_from_slice(&self.speed);
        self.initial_side_depth.copy_from_slice(&self.side_depth);
        self.initial_advanced_side_depth
            .copy_from_slice(&self.advanced_side_depth);
        for (d, (&f, &u)) in self
            .initial_increment
            .iter_mut()
            .zip(self.advected.iter().zip(&state.u_n))
        {
            *d = f - u;
        }
        let first = self.implicit_solve(mesh, state, dt, predictor)?;

        // 2. 予測子 U* の陽的な部分
        let explicit = self.explicit_terms(mesh, state, boundary, dt);

        // 3. 陽的な部分を平均して、時刻 n から解き直す
        let dry_depth = self.params.dry_depth;
        for e in 0..mesh.n_edges() {
            let depth = 0.5 * (self.initial_edge_depth[e] + self.edge_depth[e]);
            self.edge_depth[e] = if depth >= dry_depth { depth } else { 0.0 };
            self.edge_flux[e] = 0.5 * (self.initial_edge_flux[e] + self.edge_flux[e]);
            self.speed[e] = 0.5 * (self.initial_speed[e] + self.speed[e]);
            // 移流：辺の運動量 (h_1 + h_2) u と体積 h_1 + h_2 を、2つの状態の流束の平均で陽的に
            // 進める（ADR-0008 の h* で割る形の平均）。流速の増分をそのまま平均すると、動く段波で
            // 運動量が保存されない
            let (u0, d0) = (self.initial_velocity[e], self.initial_increment[e]);
            let (h0, s0) = (
                self.initial_side_depth[e],
                self.initial_advanced_side_depth[e],
            );
            let (u1, f1) = (state.u_n[e], self.advected[e]);
            let (h1, s1) = (self.side_depth[e], self.advanced_side_depth[e]);
            let volume = 0.5 * (h0 + s0 + s1 - h1);
            self.advected[e] = if s0 > 0.0 && s1 > 0.0 && volume > 0.0 {
                0.5 * (h0 * u0 + s0 * (u0 + d0) + s1 * f1 - h1 * u1) / volume
            } else {
                u0 + 0.5 * (d0 + f1 - u1)
            };
        }
        state.eta.copy_from_slice(&self.initial_level);
        state.u_n.copy_from_slice(&self.initial_velocity);
        reconstruct_cell_velocity(mesh, &state.u_n, &mut state.u_cell);
        state.time = time;
        let corrector = StepReport {
            dt,
            max_courant: first.max_courant.max(explicit.max_courant),
            limited_edges: first.limited_edges.max(explicit.limited_edges),
            ..StepReport::default()
        };
        let second = self.implicit_solve(mesh, state, dt, corrector)?;
        Ok(StepReport {
            outflow_courant: first.outflow_courant,
            newton_iterations: first.newton_iterations + second.newton_iterations,
            pcg_iterations: first.pcg_iterations + second.pcg_iterations,
            residual: first.residual.max(second.residual),
            ..second
        })
    }

    /// 1ステップの手順の 1〜2（陽的な部分）：辺の役割、辺の水深 H_j、辺の流量 q_j、移流 F u_j
    /// （`advected`）、摩擦係数に使う |u|_j（`speed`）を `state` から求める。
    ///
    /// 返す記録には、Δt、セルの流出量の比、移流の局所クーラン数を入れる。
    fn explicit_terms(
        &mut self,
        mesh: &Mesh,
        state: &State,
        boundary: &BoundaryConfig,
        dt: f64,
    ) -> StepReport {
        self.prepare_edges(mesh, state, boundary);
        let Self {
            params,
            edge_flux,
            momentum_edge,
            advection: advection_workspace,
            advected,
            speed,
            ..
        } = self;
        let mut report = StepReport {
            dt,
            outflow_courant: dt * outflow_rate(mesh, state, edge_flux),
            ..StepReport::default()
        };
        let scheme = match params.advection {
            AdvectionScheme::Upwind => Some(advection as AdvectionFn),
            AdvectionScheme::SecondOrder => Some(advection_second_order as AdvectionFn),
            AdvectionScheme::Disabled => None,
        };
        match scheme {
            Some(scheme) => {
                let stats = scheme(
                    mesh,
                    state,
                    edge_flux,
                    momentum_edge,
                    dt,
                    advection_workspace,
                    advected,
                );
                report.max_courant = stats.max_courant;
                report.limited_edges = stats.limited_edges;
            }
            None => advected.copy_from_slice(&state.u_n),
        }
        for (e, s) in speed.iter_mut().enumerate() {
            *s = edge_speed(mesh, state, e);
        }
        // 移流と同じ h*_i と、辺の両側の和（予測子・修正子で運動量を平均するため）
        let Self {
            edge_flux,
            advanced_depth,
            side_depth,
            advanced_side_depth,
            ..
        } = self;
        for (c, h) in advanced_depth.iter_mut().enumerate() {
            let outflow: f64 = mesh.edges_on_cell[c]
                .iter()
                .zip(&mesh.edge_sign_on_cell[c])
                .map(|(&e, &sign)| sign * mesh.dv_edge[e] * edge_flux[e])
                .sum();
            *h = (state.depth(c) - dt * outflow / mesh.area_cell[c]).max(0.0);
        }
        for (e, (h, h_star)) in side_depth
            .iter_mut()
            .zip(advanced_side_depth.iter_mut())
            .enumerate()
        {
            (*h, *h_star) = (0.0, 0.0);
            for c in mesh.cells_on_edge[e] {
                if c != NONE {
                    *h += state.depth(c);
                    *h_star += advanced_depth[c];
                }
            }
        }
        report
    }

    /// 1ステップの手順の 3〜6（陰的な部分）：[`HydroSolver::explicit_terms`] が求めた辺の水深・
    /// 流量・F u・|u| を使って、水位の連立方程式を解き、`state` を `dt` だけ進める。
    ///
    /// `report` の Newton 法・PCG の反復回数と残差を書き足して返す。
    fn implicit_solve(
        &mut self,
        mesh: &Mesh,
        state: &mut State,
        dt: f64,
        mut report: StepReport,
    ) -> Result<StepReport, HydroError> {
        let Self {
            params,
            role,
            outside_level,
            edge_depth,
            edge_flux,
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
            speed,
            ..
        } = self;
        let (g, theta) = (params.gravity, params.theta);
        let n2 = params.manning_n * params.manning_n;
        // 周期境界の体積力 g I（x 方向）
        let body_force = g * params.background_slope;

        // 3. 摩擦係数、G_j、重み w_j、右辺に入る流量（流量境界では配分した流入量）
        for e in 0..mesh.n_edges() {
            let h = edge_depth[e];
            let (mut gamma, mut big_g, mut w, mut flux) = (0.0, 0.0, 0.0, 0.0);
            if role[e] == EdgeRole::Discharge {
                flux = mesh.dv_edge[e] * edge_flux[e];
            } else if h > 0.0 {
                if let Some((eta1, eta2, delta)) = gradient_points(mesh, role[e], &state.eta, e) {
                    gamma = g * n2 * speed[e] / h.powf(4.0 / 3.0);
                    let denom = 1.0 + dt * gamma;
                    big_g = (advected[e] - (1.0 - theta) * g * dt * (eta2 - eta1) / delta
                        + dt * body_force * mesh.normal_edge[e][0])
                        / denom;
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
            if role[e] == EdgeRole::Discharge {
                // 流入の単位幅流量を与え、流速は新しい水深から求める（Perot 再構成と出力用）
                q = explicit_flux[e] / mesh.dv_edge[e];
                let depth = state.depth(mesh.cells_on_edge[e][0]).max(params.dry_depth);
                u_new = q / depth;
            } else if h > 0.0 {
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

/// セルの流出量の比 max_i Σ_{流出する辺} L_j q_j / V_i \[1/s\]（流出がなければ 0）。
///
/// 流出する辺の水深は風上側（このセル）の水位から求めるので、流出のあるセルは水深が
/// `dry_depth` 以上で、V_i > 0 になっている。
fn outflow_rate(mesh: &Mesh, state: &State, edge_flux: &[f64]) -> f64 {
    let mut max = 0.0_f64;
    for c in 0..mesh.n_cells() {
        let outflow: f64 = mesh.edges_on_cell[c]
            .iter()
            .zip(&mesh.edge_sign_on_cell[c])
            .map(|(&e, &s)| (s * mesh.dv_edge[e] * edge_flux[e]).max(0.0))
            .sum();
        if outflow > 0.0 {
            let v = volume(mesh.area_cell[c], state.eta[c], state.zb[c]);
            max = max.max(outflow / v);
        }
    }
    max
}

/// 辺 `e` での流速の大きさ |u|_j \[m/s\]。
///
/// 法線成分は u_j、接線成分は [`edge_tangential_velocity`] とする。摩擦係数と、掃流砂の掃流力
/// （[`crate::sediment`]）で同じ値を使う。一様流では、Perot 法の再構成が厳密なので、元の速さに
/// 丸め誤差の範囲で一致する。
pub fn edge_speed(mesh: &Mesh, state: &State, e: usize) -> f64 {
    state.u_n[e].hypot(edge_tangential_velocity(mesh, state, e))
}

/// 辺 `e` での流速の接線成分 u_t = **u**·(**k** × **n**_j) \[m/s\]（**k** × **n**_j は法線を
/// 反時計回りに 90° 回した向き）。
///
/// 両側のセルの Perot 流速の平均（境界辺では内側のセルの値）の接線方向の成分。
pub fn edge_tangential_velocity(mesh: &Mesh, state: &State, e: usize) -> f64 {
    let [c1, c2] = mesh.cells_on_edge[e];
    let [mut vx, mut vy] = state.u_cell[c1];
    if c2 != NONE {
        let [wx, wy] = state.u_cell[c2];
        vx = 0.5 * (vx + wx);
        vy = 0.5 * (vy + wy);
    }
    let [nx, ny] = mesh.normal_edge[e];
    vy * nx - vx * ny
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
/// `background_slope` は、周期境界で河床高と水位から除いた縦断勾配 I（周期境界でなければ 0）。
/// 隣のセルの水位と河床高に、除いた勾配の分 −I (x_2 − x_1)（x_2 はセル 1 の側から見た母点の位置）
/// を戻してから、風上と高い方の河床を選ぶ。勾配を除いたまま max をとると、河床が下流へ上がる
/// 所で下流側のセルが選ばれ（風下）、河床に 2Δx の振動が育つ（ADR-0015）。
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
    background_slope: f64,
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
            // セル 2 の実際の高さを、セル 1 の基準に直す（周期境界でなければ 0）
            let trend = if background_slope == 0.0 {
                0.0
            } else {
                background_slope * (mesh.neighbor_position_from(e, c1)[0] - mesh.x_cell[c1])
            };
            (state.eta[c2] - trend, zb1.max(state.zb[c2] - trend))
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
/// 周期境界の継ぎ目の辺では、中点をセルの側の座標（[`Mesh::edge_midpoint_from`]）で使う。
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
            let [xe, ye] = mesh.edge_midpoint_from(e, c);
            sx += flux * (xe - xc);
            sy += flux * (ye - yc);
        }
        *u = [sx / mesh.area_cell[c], sy / mesh.area_cell[c]];
    }
}

/// 移流の作業配列（セルごと）と、2次の再構成に使う最小二乗の重み。
#[derive(Debug, Clone)]
pub struct AdvectionWorkspace {
    /// 流れ込む運動量 P_i = (1 / A_i) Σ_{流入側の辺} |L_l q_l| u_{k(i,l)} \[m²/s²\]
    upwind: Vec<[f64; 2]>,
    /// 流入量 S_i / A_i \[m/s\]
    inflow: Vec<f64>,
    /// 時刻 n の流量で連続式を陽的に進めた水深 h*_i \[m\]
    depth: Vec<f64>,
    /// セル中心の流速の最小二乗の重み（2次の形だけで使う）
    weights: CellGradient,
    /// 内部の辺がすべて濡れているセル（2次の形だけで使う）
    wet: Vec<bool>,
    /// セル中心の流速の勾配 ∂U_a/∂x_b（`[a][b]`）\[1/s\]（2次の形だけで使う）
    gradient: Vec<[[f64; 2]; 2]>,
    /// 面の値の運動量の収支 B_i = −(1 / A_i) Σ_l s_il L_l q_l û_l \[m²/s²\]（2次の形だけで使う）
    balance: Vec<[f64; 2]>,
    /// 正味の流出量 N_i = (1 / A_i) Σ_l s_il L_l q_l \[m/s\]（2次の形だけで使う）
    net: Vec<f64>,
}

impl AdvectionWorkspace {
    /// 格子 `mesh` 用に確保し、最小二乗の重みを計算する。
    pub fn new(mesh: &Mesh) -> Self {
        let n_cells = mesh.n_cells();
        Self {
            upwind: vec![[0.0; 2]; n_cells],
            inflow: vec![0.0; n_cells],
            depth: vec![0.0; n_cells],
            weights: CellGradient::new(mesh),
            wet: vec![false; n_cells],
            gradient: vec![[[0.0; 2]; 2]; n_cells],
            balance: vec![[0.0; 2]; n_cells],
            net: vec![0.0; n_cells],
        }
    }
}

/// [`advection`] の記録。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct AdvectionStats {
    /// 局所クーラン数 σ_j の最大値（重みを 1 で頭打ちにする前）\[-\]。
    pub max_courant: f64,
    /// σ_j > 1 のために重みを 1 で頭打ちにした辺の数（時間刻みが大きすぎる兆候）。
    pub limited_edges: usize,
}

/// 移流項を保存型の1次風上差分で陽的に評価し、F u_j を `out` に書く（ADR-0008, ADR-0009）。
///
/// 各辺の流速を、両側のセルに流れ込む水の流速（法線成分）へ、局所クーラン数 σ_j の割合で
/// 近づける。
///
/// ```text
/// S_i   = Σ_{l ∈ ∂i, s_il q_l < 0} |L_l q_l|                     （セル i への流入量）
/// P_i   = (1 / A_i) Σ_{l ∈ ∂i, s_il q_l < 0} |L_l q_l| u_{k(i,l)}
/// h*_i  = h_i − (Δt / A_i) Σ_{l ∈ ∂i} s_il L_l q_l
/// σ_j   = Δt (S_1 / A_1 + S_2 / A_2) / (h*_1 + h*_2)
/// ū_j   = (P_1 + P_2)·n_j / (S_1 / A_1 + S_2 / A_2)
/// F u_j = u_j + min(σ_j, 1) (ū_j − u_j)
/// ```
///
/// σ_j ≤ 1 のとき、これは u_j − Δt [Σ_i (S_i / A_i) u_j − (P_1 + P_2)·n_j] / (h*_1 + h*_2) と同じで、
/// 運動量の流束 Σ s L q u_up から連続式の項 u_j Σ s L q を引いた保存形の移流を、流入側の辺で
/// 評価したものにあたる。h* は、セルの運動量と体積を時刻 n の同じ流量で陽的に進めたときの
/// 体積から求めた水深で、定常流では h に一致する。
///
/// - 辺の流速そのものを上流の流速へ近づけるので、セル中心の流速に現れない格子規模の辺の
///   流速の成分（六角形配置で生じる）も下流へ運び去られる。
/// - 一様流では ū_j = u_j なので、F u_j = u_j（Perot 再構成が一様流で厳密なため）。
/// - 1次元の定常流では、水位勾配の項と合わせて q u + g h²/2 が離散的に保存され、跳水の
///   共役水深が Bélanger の式を満たす（`docs/design/adr/0009-advection-face-velocity.md`）。
/// - セルからの流出が体積を超えない（時間刻みの条件。`docs/design/overview.md` の第7章）限り
///   σ_j ≤ 1 で、F u_j は u_j と ū_j の凸結合になる。σ_j > 1 の辺では重みを 1 にする（安全策）。
/// - `flux` は辺の単位幅流量 q_j^n \[m²/s\]（向きは u_n と同じ）。壁と乾いた辺では 0。
/// - u_k は `state.u_cell`（Perot 法で再構成したセル中心の流速）。境界辺（隣のセルが [`NONE`]）
///   から流れ込む水の流速は、その辺の法線流速 u_l n_l とする。
/// - `active` が真の辺だけを計算し、それ以外の辺は F u_j = u_j とする。境界辺はセル 1(j) の
///   値だけで評価する。h*_1 + h*_2 = 0 の辺と、どちらのセルにも流入がない辺も F u_j = u_j。
///
/// 計算は2つのループに分け、セルのループでは各セルが自分の値だけを、辺のループでは
/// 各辺が自分の F u_j だけを書く。`flux`, `active`, `out` の長さは辺の数、`workspace` は
/// セル数に合わせて確保したものでなければならない。
pub fn advection(
    mesh: &Mesh,
    state: &State,
    flux: &[f64],
    active: &[bool],
    dt: f64,
    workspace: &mut AdvectionWorkspace,
    out: &mut [f64],
) -> AdvectionStats {
    let n_edges = mesh.n_edges();
    debug_assert_eq!(flux.len(), n_edges, "flux must have one value per edge");
    debug_assert_eq!(active.len(), n_edges, "active must have one value per edge");
    debug_assert_eq!(out.len(), n_edges, "out must have one value per edge");
    debug_assert_eq!(
        workspace.inflow.len(),
        mesh.n_cells(),
        "the workspace must have one value per cell"
    );

    // セルに流れ込む運動量 P_i、流入量 S_i / A_i、陽的に進めた水深 h*_i
    for (c, ((p, s), hs)) in workspace
        .upwind
        .iter_mut()
        .zip(workspace.inflow.iter_mut())
        .zip(workspace.depth.iter_mut())
        .enumerate()
    {
        let (mut px, mut py, mut inflow, mut net) = (0.0, 0.0, 0.0, 0.0);
        let neighbors = mesh.edges_on_cell[c]
            .iter()
            .zip(&mesh.edge_sign_on_cell[c])
            .zip(&mesh.cells_on_cell[c]);
        for ((&e, &sign), &k) in neighbors {
            // セル c から流れ出る向きを正とする流量 [m³/s]
            let q = sign * mesh.dv_edge[e] * flux[e];
            net += q;
            if q < 0.0 {
                let [vx, vy] = if k == NONE {
                    let [nx, ny] = mesh.normal_edge[e];
                    [state.u_n[e] * nx, state.u_n[e] * ny]
                } else {
                    state.u_cell[k]
                };
                px -= q * vx;
                py -= q * vy;
                inflow -= q;
            }
        }
        let area = mesh.area_cell[c];
        *p = [px / area, py / area];
        *s = inflow / area;
        *hs = (state.depth(c) - dt * net / area).max(0.0);
    }

    // 辺の F u_j
    let mut stats = AdvectionStats::default();
    for (e, f) in out.iter_mut().enumerate() {
        let u = state.u_n[e];
        *f = u;
        if !active[e] {
            continue;
        }
        let (mut px, mut py, mut inflow, mut depth) = (0.0, 0.0, 0.0, 0.0);
        for c in mesh.cells_on_edge[e] {
            if c != NONE {
                let [x, y] = workspace.upwind[c];
                px += x;
                py += y;
                inflow += workspace.inflow[c];
                depth += workspace.depth[c];
            }
        }
        if depth <= 0.0 || inflow <= 0.0 {
            continue;
        }
        let [nx, ny] = mesh.normal_edge[e];
        let courant = dt * inflow / depth;
        let target = (px * nx + py * ny) / inflow;
        stats.max_courant = stats.max_courant.max(courant);
        if courant > 1.0 {
            stats.limited_edges += 1;
        }
        *f = u + courant.min(1.0) * (target - u);
    }
    stats
}

/// 移流項を保存型の2次精度の差分で陽的に評価し、F u_j を `out` に書く（ADR-0016）。
///
/// セル中心の流速 U の勾配 ∇U（最小二乗。[`CellGradient`]）で、面の値と、辺へ外挿した
/// セルの流速を作る。辺 j の両側のセル c ごとに次の R_{c,j} を求め、h* の和で割る。
///
/// ```text
/// û_l     = U_k + ∇U_k (x_l − x_k)             （面 l の値。k は面 l の上流のセル）
/// ũ_{c,j} = U_c + ∇U_c (x_j − x_c)             （セル c の流速を辺 j まで外挿）
/// R_{c,j} = −(1 / A_c) Σ_{l ∈ ∂c} s_cl L_l q_l (û_l − ũ_{c,j})·n_j − (S_c / A_c) (u_j − ũ_{c,j}·n_j)
/// F u_j   = u_j + Δt (R_{1,j} + R_{2,j}) / (h*_1 + h*_2)
/// ```
///
/// - 第1項は、セルの運動量の収支を全部の面で保存形に評価したもの。
/// - 第2項は、辺の流速を、両側のセルの流速を外挿した値へ近づける。滑らかな流れではこの差は
///   格子幅の2乗で小さく、セル中心の流速に現れない格子規模の辺の流速の成分（六角形配置で
///   生じる）だけに強く効く（ADR-0009 の性質）。
/// - 勾配を 0 にすると、[`advection`] と丸め誤差の範囲で一致する。
/// - 1次元の定常流では、面の値が面ごとに1つ、外挿の勾配がセルごとに1つなので、各辺の式の和が
///   打ち消し合い、跳水で Bélanger の式が成り立つ。
/// - 時間については前進 Euler（1段）。θ 法の水位と連続式の風上の水深と組み合わせると、
///   前進 Euler のままでは不安定なので、予測子・修正子（`time_integration = "predictor_corrector"`、
///   [`HydroSolver::step`]）と組み合わせて使う（ADR-0016、ADR-0017）。
/// - 勾配は、セルとその隣のセルの内部の辺がすべて濡れている（`active`）ときだけ求め、
///   それ以外は 0 とする（汀線と、乾いた辺のある段差の近くでは1次に戻る）。開境界の辺を持つセルも
///   0 とする（境界から流れ込む流れで発散しないため）。
/// - 境界辺から流れ込む水の流速は、その辺の法線流速 u_l n_l とする（再構成しない）。
/// - σ_j > 1 の辺（時間刻みが大きすぎるとき）は、1次の形と同じく F u_j = ū_j とする（安全策）。
///
/// 引数と戻り値は [`advection`] と同じ。
pub fn advection_second_order(
    mesh: &Mesh,
    state: &State,
    flux: &[f64],
    active: &[bool],
    dt: f64,
    workspace: &mut AdvectionWorkspace,
    out: &mut [f64],
) -> AdvectionStats {
    velocity_gradients(mesh, state, active, workspace);
    second_order_update(mesh, state, flux, active, dt, workspace, out)
}

/// セル中心の流速の勾配 ∇U を `workspace.gradient` に書く。
///
/// セル c と、境界の外を除くすべての隣のセルで、内部の辺がすべて `active` のときだけ求め、
/// それ以外は 0 とする。汀線のセルの Perot 流速は乾いた辺を流速 0 として数えるので小さく出る
/// （ADR-0008）。その値を使った勾配で、隣のセルの面の値を作らないようにするため。
///
/// 開境界（上流端・下流端）の辺を持つセルも 0 とする。勾配は内側だけの片側の差で、境界から
/// 水が流れ込むとき、境界へ向かう流速の傾きを外挿で強めて発散するため（乾いた河床に下流端の
/// 水位境界から水が入る `examples/sill.toml` の初期で確認。ADR-0018）。壁の辺は数えない。
fn velocity_gradients(
    mesh: &Mesh,
    state: &State,
    active: &[bool],
    workspace: &mut AdvectionWorkspace,
) {
    let AdvectionWorkspace {
        weights,
        wet,
        gradient,
        ..
    } = workspace;
    for (c, w) in wet.iter_mut().enumerate() {
        *w = mesh.edges_on_cell[c]
            .iter()
            .zip(&mesh.cells_on_cell[c])
            .all(|(&e, &k)| k == NONE || active[e]);
    }
    let open = |c: usize| {
        mesh.edges_on_cell[c]
            .iter()
            .any(|&e| matches!(mesh.edge_kind[e], EdgeKind::Upstream | EdgeKind::Downstream))
    };
    for (c, g) in gradient.iter_mut().enumerate() {
        let neighbors = &mesh.cells_on_cell[c];
        *g = if !wet[c] || open(c) || neighbors.iter().any(|&k| k != NONE && !wet[k]) {
            [[0.0; 2]; 2]
        } else {
            weights.of_vector(mesh, &state.u_cell, c)
        };
    }
}

/// セル `c` の流速を、勾配で位置 `c` からのずれ `(dx, dy)` だけ動かした値。
fn extrapolate(state: &State, gradient: &[[[f64; 2]; 2]], c: usize, dx: f64, dy: f64) -> [f64; 2] {
    let [ux, uy] = state.u_cell[c];
    let g = gradient[c];
    [
        ux + g[0][0] * dx + g[0][1] * dy,
        uy + g[1][0] * dx + g[1][1] * dy,
    ]
}

/// 面 `e` の値 û = U_k + ∇U_k (x_e − x_k)（`k` は面の上流のセル）。
///
/// 面の両側のセルから呼んでも、同じ上流のセルの値から同じ順で計算するので、ビット単位で一致する
/// （1次元の定常流で各辺の式の和が打ち消し合う条件）。
fn face_value(
    mesh: &Mesh,
    state: &State,
    gradient: &[[[f64; 2]; 2]],
    e: usize,
    k: usize,
) -> [f64; 2] {
    let [xe, ye] = mesh.edge_midpoint_from(e, k);
    extrapolate(state, gradient, k, xe - mesh.x_cell[k], ye - mesh.y_cell[k])
}

/// [`advection_second_order`] の本体（勾配 `workspace.gradient` は求めてあるものとする）。
///
/// 単体テストで、勾配を 0 にして1次の形と比べるために分けている。
fn second_order_update(
    mesh: &Mesh,
    state: &State,
    flux: &[f64],
    active: &[bool],
    dt: f64,
    workspace: &mut AdvectionWorkspace,
    out: &mut [f64],
) -> AdvectionStats {
    let n_edges = mesh.n_edges();
    debug_assert_eq!(flux.len(), n_edges, "flux must have one value per edge");
    debug_assert_eq!(active.len(), n_edges, "active must have one value per edge");
    debug_assert_eq!(out.len(), n_edges, "out must have one value per edge");
    let AdvectionWorkspace {
        upwind,
        inflow,
        depth,
        gradient,
        balance,
        net,
        ..
    } = workspace;
    debug_assert_eq!(
        gradient.len(),
        mesh.n_cells(),
        "the workspace must match the mesh"
    );

    // セルの値：1次の形と同じ P_i, S_i / A_i, h*_i と、面の値の収支 B_i、正味の流出量 N_i
    for c in 0..mesh.n_cells() {
        let (mut px, mut py, mut s, mut outflow) = (0.0, 0.0, 0.0, 0.0);
        let (mut bx, mut by) = (0.0, 0.0);
        let neighbors = mesh.edges_on_cell[c]
            .iter()
            .zip(&mesh.edge_sign_on_cell[c])
            .zip(&mesh.cells_on_cell[c]);
        for ((&e, &sign), &k) in neighbors {
            // セル c から流れ出る向きを正とする流量 [m³/s]
            let q = sign * mesh.dv_edge[e] * flux[e];
            outflow += q;
            if q < 0.0 {
                let (velocity, face) = if k == NONE {
                    let [nx, ny] = mesh.normal_edge[e];
                    let v = [state.u_n[e] * nx, state.u_n[e] * ny];
                    (v, v)
                } else {
                    (state.u_cell[k], face_value(mesh, state, gradient, e, k))
                };
                px -= q * velocity[0];
                py -= q * velocity[1];
                s -= q;
                bx -= q * face[0];
                by -= q * face[1];
            } else if q > 0.0 {
                let face = face_value(mesh, state, gradient, e, c);
                bx -= q * face[0];
                by -= q * face[1];
            }
        }
        let area = mesh.area_cell[c];
        upwind[c] = [px / area, py / area];
        inflow[c] = s / area;
        depth[c] = (state.depth(c) - dt * outflow / area).max(0.0);
        balance[c] = [bx / area, by / area];
        net[c] = outflow / area;
    }

    // 辺の F u_j = u_j + Δt (R_{1,j} + R_{2,j}) / (h*_1 + h*_2)
    let mut stats = AdvectionStats::default();
    for (e, f) in out.iter_mut().enumerate() {
        let u = state.u_n[e];
        *f = u;
        if !active[e] {
            continue;
        }
        let [nx, ny] = mesh.normal_edge[e];
        let (mut px, mut py, mut s, mut h, mut rate) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for c in mesh.cells_on_edge[e] {
            if c == NONE {
                continue;
            }
            let [x, y] = upwind[c];
            px += x;
            py += y;
            s += inflow[c];
            h += depth[c];
            // R_{c,j} = B_c·n_j + N_c ũ_{c,j}·n_j − (S_c / A_c)(u_j − ũ_{c,j}·n_j)
            let [xe, ye] = mesh.edge_midpoint_from(e, c);
            let [vx, vy] =
                extrapolate(state, gradient, c, xe - mesh.x_cell[c], ye - mesh.y_cell[c]);
            let extrapolated = vx * nx + vy * ny;
            let [bx, by] = balance[c];
            rate += bx * nx + by * ny + net[c] * extrapolated - inflow[c] * (u - extrapolated);
        }
        if h <= 0.0 || s <= 0.0 {
            continue;
        }
        let courant = dt * s / h;
        stats.max_courant = stats.max_courant.max(courant);
        if courant > 1.0 {
            // 安全策：1次の形と同じく、重みを 1 にした値
            stats.limited_edges += 1;
            *f = (px * nx + py * ny) / s;
            continue;
        }
        *f = u + dt * rate / h;
    }
    stats
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
    use std::f64::consts::PI;

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
            background_slope: 0.0,
            gravity: 9.81,
            theta: 0.55,
            manning_n,
            dry_depth: 1e-3,
            advection: AdvectionScheme::Upwind,
            time_integration: TimeIntegration::Euler,
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

    /// 直線水路の側壁ぎわのセルで、一次の流速場の Perot 法の再構成の誤差が小さいこと。
    ///
    /// 側壁は最も外側の行の母点と鏡映点の二等分線なので、側壁ぎわのセルでも母点は重心の近くに
    /// ある。正方形配置は、母点が長方形のセルの中心なので厳密になる。六角形配置の上下流端の近くの
    /// セルは、端まで広がった多角形になるので除く（`mesh::channel`）。
    #[test]
    fn perot_reconstruction_of_linear_fields_is_accurate_at_side_walls() {
        let fields: [fn(f64, f64) -> [f64; 2]; 4] = [
            |x, _| [0.2 * x, 0.0],
            |_, y| [0.2 * y, 0.0],
            |x, y| [-0.2 * y, 0.2 * x],
            |x, y| [0.2 * x, -0.2 * y],
        ];
        let (length, h) = (10.0, 0.25);
        for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
            let mesh = mesh(length, 2.0, h, layout);
            let away_from_ends = |c: usize| (2.0 * h..=length - 2.0 * h).contains(&mesh.x_cell[c]);
            let at_wall = |c: usize| {
                mesh.edges_on_cell[c]
                    .iter()
                    .any(|&e| mesh.is_boundary_edge(e))
            };
            for field in fields {
                let u_n: Vec<f64> = (0..mesh.n_edges())
                    .map(|e| {
                        let [u, v] = field(mesh.x_edge[e], mesh.y_edge[e]);
                        u * mesh.normal_edge[e][0] + v * mesh.normal_edge[e][1]
                    })
                    .collect();
                let mut u_cell = vec![[0.0; 2]; mesh.n_cells()];
                reconstruct_cell_velocity(&mesh, &u_n, &mut u_cell);
                // 母点での厳密な値との差 / 1セル分の流速の変化 0.2h
                let error = |c: usize| {
                    let [u, v] = field(mesh.x_cell[c], mesh.y_cell[c]);
                    (u_cell[c][0] - u).hypot(u_cell[c][1] - v) / (0.2 * h)
                };
                let max_of = |wall: bool| {
                    (0..mesh.n_cells())
                        .filter(|&c| away_from_ends(c) && at_wall(c) == wall)
                        .map(error)
                        .fold(0.0, f64::max)
                };
                let (interior, wall) = (max_of(false), max_of(true));
                // 内部の六角形は母点について対称なので厳密。側壁ぎわの五角形は母点が重心から
                // わずかにずれ（約 0.005h）、誤差は 0.06 程度
                let wall_limit = match layout {
                    ChannelLayout::Rectangular => 1e-12,
                    ChannelLayout::Hexagonal => 0.1,
                };
                assert!(interior < 1e-12, "{layout:?}: interior error {interior}");
                assert!(wall < wall_limit, "{layout:?}: wall error {wall}");
            }
        }
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
            edge_depths(&mesh, &state, &outside, 1e-3, 0.0, &mut h);
            h[e]
        };
        // 河床高は高い方（0.5）、水位は風上側
        assert!((depth_for(0.3) - 0.5).abs() < 1e-15);
        assert!((depth_for(-0.3) - 0.3).abs() < 1e-15);
        assert!((depth_for(0.0) - 0.5).abs() < 1e-15);
    }

    /// 周期境界の辺の水深は、勾配を除いた河床高と水位から、勾配のある水路の実際の値で求めた
    /// ものと一致する（ADR-0015）。勾配を戻さずに max をとると、勾配を除いた河床が下流へ上がる
    /// 辺で下流側のセルが選ばれ、一致しない。
    #[test]
    fn periodic_edge_depths_match_the_sloped_channel() {
        let (slope, length) = (0.01, 20.0);
        let spec = ChannelSpec {
            length,
            width: 4.0,
            cell_size: 1.0,
            layout: ChannelLayout::Rectangular,
        };
        let periodic = channel::build_periodic(&spec).expect("valid periodic channel");
        let open = channel::build(&spec).expect("valid channel");
        let detrended = |m: &Mesh| {
            let mut state = State::new(m.n_cells(), m.n_edges());
            for c in 0..m.n_cells() {
                let (x, y) = (m.x_cell[c], m.y_cell[c]);
                // 勾配を除いた河床は下流へ上がる所と下がる所がある（勾配 I より緩い）
                state.zb[c] = 0.02 * (2.0 * PI * x / length).sin() * (PI * y / 4.0).cos();
                state.eta[c] = 0.8 + 0.01 * (2.0 * PI * x / length).cos();
            }
            for e in 0..m.n_edges() {
                state.u_n[e] = 0.5 * m.normal_edge[e][0] + 0.05 * m.normal_edge[e][1];
            }
            state
        };
        let mut actual = detrended(&open);
        for c in 0..open.n_cells() {
            let trend = slope * (length - open.x_cell[c]);
            actual.zb[c] += trend;
            actual.eta[c] += trend;
        }
        let depths = |m: &Mesh, state: &State, background: f64| {
            let mut h = vec![0.0; m.n_edges()];
            let outside = vec![f64::NEG_INFINITY; m.n_edges()];
            edge_depths(m, state, &outside, 1e-3, background, &mut h);
            h
        };
        let h_open = depths(&open, &actual, 0.0);
        let state = detrended(&periodic);
        let h_periodic = depths(&periodic, &state, slope);
        let h_unadjusted = depths(&periodic, &state, 0.0);
        let mut compared = 0;
        let mut differs = 0;
        for e in (0..open.n_edges()).filter(|&e| !open.is_boundary_edge(e)) {
            let cells = open.cells_on_edge[e];
            let p = (0..periodic.n_edges())
                .find(|&p| periodic.cells_on_edge[p] == cells)
                .expect("the interior edges are the same");
            assert!(
                (h_periodic[p] - h_open[e]).abs() < 1e-12,
                "edge {cells:?}: {} vs {}",
                h_periodic[p],
                h_open[e]
            );
            compared += 1;
            if (h_unadjusted[p] - h_open[e]).abs() > 1e-6 {
                differs += 1;
            }
        }
        assert!(compared > 50);
        assert!(differs > 0, "the unadjusted maximum must differ somewhere");
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
        edge_depths(&mesh, &state, &outside, 1e-3, 0.0, &mut h);
        assert_eq!(h[e], 0.0);

        // 閾値より浅い
        state.eta[c1] = 0.5005;
        state.eta[c2] = 0.5;
        state.zb[c1] = 0.0;
        state.zb[c2] = 0.5;
        edge_depths(&mesh, &state, &outside, 1e-3, 0.0, &mut h);
        assert_eq!(h[e], 0.0);
        edge_depths(&mesh, &state, &outside, 1e-4, 0.0, &mut h);
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
        for time_integration in [TimeIntegration::Euler, TimeIntegration::PredictorCorrector] {
            let mut state = initial.clone();
            let params = HydroParams {
                time_integration,
                ..params(0.03)
            };
            let mut solver = HydroSolver::new(&mesh, params);
            for _ in 0..10 {
                let report = solver.step(&mesh, &mut state, &closed(), 0.5).unwrap();
                assert_eq!(report.newton_iterations, 0);
                assert_eq!(report.residual, 0.0);
            }
            assert_eq!(state.eta, initial.eta, "{time_integration:?}");
            assert!(state.u_n.iter().all(|&u| u == 0.0));
            assert_eq!(state.time, 5.0);
        }
    }

    /// 周期境界の等流は1段の定常解なので、予測子・修正子でも（1次・2次の移流とも、正方形・
    /// 六角形配置とも）変わらない。予測子・修正子の定常解は Δt によらない。
    #[test]
    fn predictor_corrector_keeps_the_periodic_uniform_flow() {
        let (slope, manning_n, q): (f64, f64, f64) = (0.002, 0.025, 2.0);
        let depth = (manning_n * q / slope.sqrt()).powf(0.6);
        for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
            let mesh = channel::build_periodic(&ChannelSpec {
                length: 40.0,
                width: 10.0,
                cell_size: 2.0,
                layout,
            })
            .expect("valid periodic channel");
            let mut initial = State::new(mesh.n_cells(), mesh.n_edges());
            initial.eta.fill(depth);
            initial.u_n = mesh.normal_edge.iter().map(|n| q / depth * n[0]).collect();
            reconstruct_cell_velocity(&mesh, &initial.u_n, &mut initial.u_cell);
            for advection in [AdvectionScheme::Upwind, AdvectionScheme::SecondOrder] {
                let params = HydroParams {
                    advection,
                    time_integration: TimeIntegration::PredictorCorrector,
                    background_slope: slope,
                    ..params(manning_n)
                };
                let mut solver = HydroSolver::new(&mesh, params);
                let mut state = initial.clone();
                for _ in 0..10 {
                    solver.step(&mesh, &mut state, &closed(), 1.0).unwrap();
                }
                let label = format!("{layout:?}, {advection:?}");
                for (&eta, &eta0) in state.eta.iter().zip(&initial.eta) {
                    assert!((eta - eta0).abs() <= 1e-12, "{label}: {eta} != {eta0}");
                }
                for (&u, &u0) in state.u_n.iter().zip(&initial.u_n) {
                    assert!((u - u0).abs() <= 1e-12, "{label}: {u} != {u0}");
                }
            }
        }
    }

    /// セルの流出量の比は、正方形配置の一様流では u/Δx になり（流量境界の流入は数えない）、
    /// 流れがなければ 0 になる。`step` の記録の `outflow_courant` は、その Δt 倍。
    #[test]
    fn outflow_rate_is_velocity_over_cell_size_in_uniform_flow() {
        let mesh = mesh(10.0, 2.0, 0.5, ChannelLayout::Rectangular);
        let (u, h) = (0.8, 1.0);
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        state.eta.fill(h);
        for e in 0..mesh.n_edges() {
            state.u_n[e] = u * mesh.normal_edge[e][0];
        }
        reconstruct_cell_velocity(&mesh, &state.u_n, &mut state.u_cell);
        let boundary = BoundaryConfig {
            upstream: UpstreamBoundary::Discharge {
                discharge: u * h * 2.0,
                sediment_supply: crate::config::SedimentSupply::default(),
            },
            downstream: DownstreamBoundary::WaterLevel { level: h },
        };
        let mut solver = HydroSolver::new(&mesh, params(0.03));
        let rate = solver.max_outflow_rate(&mesh, &state, &boundary);
        assert!((rate - u / 0.5).abs() < 1e-12, "rate = {rate}");

        let dt = 0.1;
        let report = solver
            .step(&mesh, &mut state.clone(), &boundary, dt)
            .expect("the step must succeed");
        assert_eq!(report.dt, dt);
        assert!((report.outflow_courant - dt * rate).abs() < 1e-15);

        state.u_n.fill(0.0);
        assert_eq!(solver.max_outflow_rate(&mesh, &state, &closed()), 0.0);
    }

    /// 上流端の流量は、断面で共通の水位から求めた通水能に比例して配分される。
    #[test]
    fn discharge_is_distributed_by_conveyance_with_a_common_level() {
        let mesh = mesh(4.0, 4.0, 1.0, ChannelLayout::Rectangular);
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        for c in 0..mesh.n_cells() {
            // 横断方向に河床高が変わる断面。水位は一様。
            state.zb[c] = 0.1 * mesh.y_cell[c];
            state.eta[c] = 1.0;
        }
        let boundary = BoundaryConfig {
            upstream: UpstreamBoundary::Discharge {
                discharge: 3.0,
                sediment_supply: crate::config::SedimentSupply::None,
            },
            downstream: DownstreamBoundary::Depth { depth: 0.5 },
        };
        let mut solver = HydroSolver::new(&mesh, params(0.03));
        let zb_before = state.zb.clone();
        solver.step(&mesh, &mut state, &boundary, 0.1).unwrap();

        let upstream: Vec<usize> = (0..mesh.n_edges())
            .filter(|&e| mesh.edge_kind[e] == EdgeKind::Upstream)
            .collect();
        assert_eq!(upstream.len(), 4);
        let inflow: f64 = upstream
            .iter()
            .map(|&e| -mesh.dv_edge[e] * solver.theta_flux()[e])
            .sum();
        assert!((inflow - 3.0).abs() < 1e-12, "inflow {inflow}");
        for &e in &upstream {
            let c = mesh.cells_on_edge[e][0];
            let expected = 3.0 * (1.0 - zb_before[c]).powf(5.0 / 3.0)
                / upstream
                    .iter()
                    .map(|&f| (1.0 - zb_before[mesh.cells_on_edge[f][0]]).powf(5.0 / 3.0))
                    .sum::<f64>();
            let share = -mesh.dv_edge[e] * solver.theta_flux()[e];
            assert!((share - expected).abs() < 1e-12, "{share} != {expected}");
            assert!(state.u_n[e] < 0.0, "inflow must point into the domain");
        }

        // 水深一定の境界の水位は、内側のセルの河床高 + 水深
        for e in (0..mesh.n_edges()).filter(|&e| mesh.edge_kind[e] == EdgeKind::Downstream) {
            let c = mesh.cells_on_edge[e][0];
            assert_eq!(
                edge_role(&mesh, &state, &boundary, e),
                EdgeRole::Level(state.zb[c] + 0.5)
            );
        }
    }

    /// 状態の u_n から u_cell を再構成し、辺の単位幅流量 q_j = H_j u_j と、水深のある辺を
    /// 計算の対象とする `active` を返す。境界辺の H_j は内側のセルの水深とする（流入も表す）。
    fn advection_inputs(mesh: &Mesh, state: &mut State) -> (Vec<f64>, Vec<bool>) {
        reconstruct_cell_velocity(mesh, &state.u_n, &mut state.u_cell);
        let outside = vec![f64::NEG_INFINITY; mesh.n_edges()];
        let mut h = vec![0.0; mesh.n_edges()];
        edge_depths(mesh, state, &outside, 1e-3, 0.0, &mut h);
        for (e, depth) in h.iter_mut().enumerate() {
            if mesh.is_boundary_edge(e) {
                *depth = state.depth(mesh.cells_on_edge[e][0]);
            }
        }
        let flux = (0..mesh.n_edges()).map(|e| h[e] * state.u_n[e]).collect();
        let active = h.iter().map(|&d| d > 0.0).collect();
        (flux, active)
    }

    /// 1次と2次の移流の関数（テストで両方を確かめるため）。
    const SCHEMES: [(&str, AdvectionFn); 2] = [
        ("upwind", advection),
        ("second_order", advection_second_order),
    ];

    /// 一様流では、Perot 再構成が厳密なので移流の項が 0 になる（境界の流入を含む）。
    /// 2次の形でも、勾配が 0 になるので同じ。周期境界の継ぎ目でも同じ。
    #[test]
    fn advection_vanishes_in_uniform_flow() {
        let spec = |layout| ChannelSpec {
            length: 20.0,
            width: 6.0,
            cell_size: 1.0,
            layout,
        };
        let meshes = [ChannelLayout::Rectangular, ChannelLayout::Hexagonal].map(|layout| {
            [
                channel::build(&spec(layout)).expect("valid channel"),
                channel::build_periodic(&spec(layout)).expect("valid channel"),
            ]
        });
        for (name, scheme) in SCHEMES {
            for mesh in meshes.iter().flatten() {
                let periodic = mesh.is_periodic();
                let mut state = State::new(mesh.n_cells(), mesh.n_edges());
                state.eta.fill(1.3);
                // 側壁に平行な一様流（上流端から流入し、下流端から流出する）
                state.u_n = mesh.normal_edge.iter().map(|n| 0.8 * n[0]).collect();
                let (flux, active) = advection_inputs(mesh, &mut state);
                let mut out = vec![0.0; mesh.n_edges()];
                let mut work = AdvectionWorkspace::new(mesh);
                let stats = scheme(mesh, &state, &flux, &active, 0.5, &mut work, &mut out);
                assert_eq!(stats.limited_edges, 0);
                for (e, (&f, &u)) in out.iter().zip(&state.u_n).enumerate() {
                    assert!(
                        (f - u).abs() <= 1e-13,
                        "{name}, periodic {periodic}: edge {e}: {f} != {u}"
                    );
                }

                // 静水では流量が 0 なので、何も変わらない
                state.u_n.fill(0.0);
                let (flux, active) = advection_inputs(mesh, &mut state);
                let stats = scheme(mesh, &state, &flux, &active, 0.5, &mut work, &mut out);
                assert!(out.iter().all(|&u| u == 0.0));
                assert_eq!(stats.max_courant, 0.0);
            }
        }
    }

    /// 六角形配置の水路に、x 方向に流れる滑らかな流速場を与える（側壁で法線成分が 0）。
    fn smooth_flow(mesh: &Mesh, state: &mut State, field: impl Fn(f64, f64) -> [f64; 2]) {
        state.eta.fill(1.0);
        for e in 0..mesh.n_edges() {
            let [ux, uy] = field(mesh.x_edge[e], mesh.y_edge[e]);
            let [nx, ny] = mesh.normal_edge[e];
            state.u_n[e] = ux * nx + uy * ny;
        }
    }

    /// 2次の形で勾配を 0 にすると、1次の形と丸め誤差の範囲で一致する。
    #[test]
    fn second_order_without_gradients_matches_upwind() {
        let mesh = mesh(20.0, 6.0, 1.0, ChannelLayout::Hexagonal);
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        smooth_flow(&mesh, &mut state, |x, y| {
            [
                1.0 + 0.3 * (0.4 * x).sin() * (0.5 * y).cos(),
                0.2 * (0.3 * x).cos() * (PI * y / 6.0).sin(),
            ]
        });
        for c in 0..mesh.n_cells() {
            state.eta[c] = 1.0 + 0.1 * (0.2 * mesh.x_cell[c]).cos();
        }
        let (flux, active) = advection_inputs(&mesh, &mut state);
        let mut work = AdvectionWorkspace::new(&mesh);
        let (mut first, mut second) = (vec![0.0; mesh.n_edges()], vec![0.0; mesh.n_edges()]);
        let reference = advection(&mesh, &state, &flux, &active, 0.3, &mut work, &mut first);
        work.gradient.fill([[0.0; 2]; 2]);
        let stats = second_order_update(&mesh, &state, &flux, &active, 0.3, &mut work, &mut second);
        assert_eq!(stats.max_courant, reference.max_courant);
        assert!(
            first
                .iter()
                .zip(&state.u_n)
                .any(|(f, u)| (f - u).abs() > 1e-3)
        );
        for (e, (&a, &b)) in first.iter().zip(&second).enumerate() {
            assert!((a - b).abs() <= 1e-14, "edge {e}: {a} != {b}");
        }
    }

    /// 滑らかな流速場の移流項 −(u·∇)u·n の誤差が、2次の形では格子幅の約2乗で、1次の形では
    /// 約1乗で減る（周期的な水路、正方形配置と六角形配置。横断方向に 10・20・40 行）。
    /// 正方形配置は 1.97・1.97、六角形配置は 1.97・1.72（側壁ぎわのセルを除くと 2.0）。
    ///
    /// 誤差は、辺ごとの誤差を Perot 法でセル中心へ再構成した値（セル中心の流速に現れる成分）の
    /// RMS で測る。六角形配置では、辺ごとの誤差に格子幅の1乗の成分が残るが、それはセル中心の
    /// 流速に現れない格子規模の成分で（辺の流速を外挿した値へ近づける項の、辺の向きによる
    /// 2次の誤差が両側で打ち消し合わないため）、移流で減衰する。運動量式を解かない側壁の辺の
    /// 誤差は 0 として数える。
    #[test]
    fn second_order_advection_converges_faster() {
        let (length, width) = (40.0, 20.0);
        let (k, m) = (2.0 * PI / length, PI / width);
        let field = |x: f64, y: f64| {
            [
                1.0 + 0.2 * (k * x).sin() * (m * y).cos(),
                0.1 * (k * x).cos() * (m * y).sin(),
            ]
        };
        // (u·∇)u
        let exact = |x: f64, y: f64| {
            let [ux, uy] = field(x, y);
            let dux = [
                0.2 * k * (k * x).cos() * (m * y).cos(),
                -0.2 * m * (k * x).sin() * (m * y).sin(),
            ];
            let duy = [
                -0.1 * k * (k * x).sin() * (m * y).sin(),
                0.1 * m * (k * x).cos() * (m * y).cos(),
            ];
            [ux * dux[0] + uy * dux[1], ux * duy[0] + uy * duy[1]]
        };
        let dt = 1e-4;
        for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
            let mut errors = Vec::new();
            for rows in [10, 20, 40] {
                let cell_size = match layout {
                    ChannelLayout::Rectangular => width / rows as f64,
                    ChannelLayout::Hexagonal => width / (rows as f64 * 0.5 * 3f64.sqrt()),
                };
                let mesh = channel::build_periodic(&ChannelSpec {
                    length,
                    width,
                    cell_size,
                    layout,
                })
                .expect("valid channel");
                let mut state = State::new(mesh.n_cells(), mesh.n_edges());
                smooth_flow(&mesh, &mut state, field);
                let (flux, mut active) = advection_inputs(&mesh, &mut state);
                for (e, a) in active.iter_mut().enumerate() {
                    *a = *a && !mesh.is_boundary_edge(e);
                }
                let mut work = AdvectionWorkspace::new(&mesh);
                let mut out = vec![0.0; mesh.n_edges()];
                let mut rms = [0.0; 2];
                for ((_, scheme), r) in SCHEMES.iter().zip(rms.iter_mut()) {
                    scheme(&mesh, &state, &flux, &active, dt, &mut work, &mut out);
                    let error: Vec<f64> = (0..mesh.n_edges())
                        .map(|e| {
                            if !active[e] {
                                return 0.0;
                            }
                            let [ax, ay] = exact(mesh.x_edge[e], mesh.y_edge[e]);
                            let [nx, ny] = mesh.normal_edge[e];
                            (state.u_n[e] - out[e]) / dt - (ax * nx + ay * ny)
                        })
                        .collect();
                    let mut projected = vec![[0.0; 2]; mesh.n_cells()];
                    reconstruct_cell_velocity(&mesh, &error, &mut projected);
                    let sum: f64 = projected.iter().map(|[x, y]| x * x + y * y).sum();
                    *r = (sum / mesh.n_cells() as f64).sqrt();
                }
                errors.push(rms);
            }
            let order = |i: usize, a: usize| (errors[a][i] / errors[a + 1][i]).log2();
            for a in 0..2 {
                assert!(
                    (0.7..=1.2).contains(&order(0, a)),
                    "{layout:?}: upwind {errors:?}"
                );
                // 六角形配置の側壁ぎわのセル（五角形、片側だけの勾配）で、20〜40 行は約 1.7
                assert!(order(1, a) >= 1.6, "{layout:?}: second order {errors:?}");
            }
            // 20 行で、誤差が1次の形の 1/10 以下
            assert!(errors[1][1] <= 0.1 * errors[1][0], "{layout:?}: {errors:?}");
        }
    }

    /// 1次元の定常な跳水：両側の一様な区間の間に任意の遷移区間をはさんでも、移流と水位勾配の
    /// 項を辺の運動量 (h_1 + h_2)/2 Δx で重み付けして足すと [q u + g h²/2] の差になる。両側の
    /// 水深が Bélanger の式を満たせば和は 0 になり、離散式の上で跳水が釣り合う。2次の移流でも
    /// 同じ（面の値が面ごとに1つ、外挿の勾配がセルごとに1つなので、各辺の式の和が打ち消し合う）。
    #[test]
    fn steady_jump_balances_momentum_for_conjugate_depths() {
        const G: f64 = 9.81;
        let (q, h1) = (1.0, 0.2);
        let froude = q / (h1 * (G * h1).sqrt());
        let h2 = h1 * ((1.0 + 8.0 * froude * froude).sqrt() - 1.0) / 2.0;
        let mesh = mesh(20.0, 1.0, 1.0, ChannelLayout::Rectangular);
        assert_eq!(mesh.n_cells(), 20);
        // セル番号は x の順。遷移区間（x = 8〜12）の水深は単調でなくてよい
        let transition = [0.45, 0.31, 0.9, 0.6];
        let balance = |h2: f64, name: &str, scheme: AdvectionFn| {
            let depth: Vec<f64> = (0..20)
                .map(|c| match c {
                    0..=7 => h1,
                    8..=11 => transition[c - 8],
                    _ => h2,
                })
                .collect();
            let mut state = State::new(mesh.n_cells(), mesh.n_edges());
            state.eta.clone_from(&depth);
            // 単位幅流量 q が一様な定常流：辺の水深は風上（左）のセルの水深、u_j = q / H_j
            let mut flux = vec![0.0; mesh.n_edges()];
            for (e, f) in flux.iter_mut().enumerate() {
                let [c1, c2] = mesh.cells_on_edge[e];
                let [nx, _] = mesh.normal_edge[e];
                if nx.abs() < 0.5 {
                    continue; // 側壁
                }
                let upwind = if c2 != NONE && mesh.x_cell[c2] < mesh.x_cell[c1] {
                    c2
                } else {
                    c1
                };
                *f = q * nx;
                state.u_n[e] = *f / depth[upwind];
            }
            reconstruct_cell_velocity(&mesh, &state.u_n, &mut state.u_cell);
            let active: Vec<bool> = (0..mesh.n_edges())
                .map(|e| !mesh.is_boundary_edge(e))
                .collect();
            let dt = 0.01;
            let mut out = vec![0.0; mesh.n_edges()];
            let mut work = AdvectionWorkspace::new(&mesh);
            let stats = scheme(&mesh, &state, &flux, &active, dt, &mut work, &mut out);
            assert_eq!(stats.limited_edges, 0, "{name}");
            // 両端から3セル以上内側の内部の辺について、(h_1 + h_2)/2 Δx の重みで足す
            let mut sum = 0.0;
            for e in (0..mesh.n_edges()).filter(|&e| active[e]) {
                let [c1, c2] = mesh.cells_on_edge[e];
                let x = mesh.x_edge[e];
                if !(3.0..=17.0).contains(&x) {
                    continue;
                }
                let h_bar = 0.5 * (depth[c1] + depth[c2]);
                let advective = (state.u_n[e] - out[e]) / dt;
                let gradient = G * (depth[c2] - depth[c1]) / mesh.dc_edge[e];
                sum += h_bar * mesh.dc_edge[e] * (advective + gradient) * mesh.normal_edge[e][0];
            }
            sum
        };
        let momentum = |h: f64| q * q / h + 0.5 * G * h * h;
        // Bélanger の式を満たす水深では釣り合い、満たさなければ運動量関数の差が残る
        for (name, scheme) in SCHEMES {
            let imbalance = balance(h2, name, scheme);
            assert!(imbalance.abs() <= 1e-12, "{name}: imbalance {imbalance}");
            let other = 1.1 * h2;
            let expected = momentum(other) - momentum(h1);
            let actual = balance(other, name, scheme);
            assert!(
                (actual - expected).abs() <= 1e-12 * expected.abs(),
                "{name}: {actual} != {expected}"
            );
        }
    }

    /// 六角形配置で一様流に乗せた格子規模の辺の流速の乱れは、移流を繰り返すと減衰する（辺の
    /// 流速そのものを上流の流速へ近づけるため。5回で約 0.23 倍）。セル中心の流速だけで増分を
    /// 作る形（ADR-0008）では、乱れのうちセル中心の流速に現れない成分が減らずに残り、5回で
    /// 約 0.59 倍にしかならない。
    #[test]
    fn advection_damps_grid_scale_noise_on_hexagons() {
        let mesh = mesh(20.0, 6.0, 1.0, ChannelLayout::Hexagonal);
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        state.eta.fill(1.0);
        let exact: Vec<f64> = mesh.normal_edge.iter().map(|n| 0.8 * n[0]).collect();
        state.u_n = exact.clone();
        let (flux, active) = advection_inputs(&mesh, &mut state);
        for e in (0..mesh.n_edges()).filter(|&e| !mesh.is_boundary_edge(e)) {
            state.u_n[e] += 0.05 * (12.9898 * e as f64).sin();
        }
        let noise = |u: &[f64]| {
            (0..mesh.n_edges())
                .map(|e| (u[e] - exact[e]).powi(2))
                .sum::<f64>()
                .sqrt()
        };
        let before = noise(&state.u_n);
        let mut out = vec![0.0; mesh.n_edges()];
        let mut work = AdvectionWorkspace::new(&mesh);
        // 流量は一様流のまま固定し、辺の流速だけを5回進める
        for _ in 0..5 {
            reconstruct_cell_velocity(&mesh, &state.u_n, &mut state.u_cell);
            let stats = advection(&mesh, &state, &flux, &active, 0.75, &mut work, &mut out);
            assert_eq!(stats.limited_edges, 0, "{stats:?}");
            state.u_n.copy_from_slice(&out);
        }
        let after = noise(&state.u_n);
        assert!(after <= 0.35 * before, "noise {before} -> {after}");
    }

    /// 汀線：深いセルから、濡れたばかりの浅いセルへ流れ込む。h* で割るので σ_j ≤ 1 で、F u_j は
    /// 辺の流速と上流の流速の凸結合に収まり、時間刻みが大きくても発散しない。流出が体積を
    /// 超える時間刻み（σ_j > 1）では、安全策として重みを 1 にする。
    #[test]
    fn advection_stays_bounded_at_a_wetting_front() {
        let mesh = mesh(8.0, 1.0, 1.0, ChannelLayout::Rectangular);
        let cell_at = |x: f64| {
            (0..mesh.n_cells())
                .find(|&c| (mesh.x_cell[c] - x).abs() < 0.1)
                .expect("a cell at x")
        };
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        for c in 0..mesh.n_cells() {
            let x = mesh.x_cell[c];
            // x < 3 は深く、x = 3.5 のセルは濡れたばかり、その先は乾いている
            state.eta[c] = if x < 3.0 {
                1.0
            } else if x < 4.0 {
                0.01
            } else {
                0.0
            };
        }
        for e in 0..mesh.n_edges() {
            let [c1, c2] = mesh.cells_on_edge[e];
            if c2 != NONE && mesh.x_cell[c1].max(mesh.x_cell[c2]) < 4.0 {
                state.u_n[e] = mesh.normal_edge[e][0];
            }
        }
        let (flux, active) = advection_inputs(&mesh, &mut state);
        let front = cell_at(3.5);
        let ahead = (0..mesh.n_edges())
            .find(|&e| {
                let [c1, c2] = mesh.cells_on_edge[e];
                c1 == front && c2 == cell_at(4.5) || c2 == front && c1 == cell_at(4.5)
            })
            .expect("the edge ahead of the front");
        assert!(active[ahead] && state.u_n[ahead] == 0.0);

        let mut out = vec![0.0; mesh.n_edges()];
        let mut work = AdvectionWorkspace::new(&mesh);
        for (dt, limited) in [(0.5, false), (5.0, true)] {
            let stats = advection(&mesh, &state, &flux, &active, dt, &mut work, &mut out);
            assert_eq!(stats.limited_edges > 0, limited, "dt = {dt}: {stats:?}");
            assert!(out[ahead] > 0.0, "the front must receive momentum");
            // F u_j は u_j と上流のセルの流速（法線成分）の凸結合
            for (e, (&f, &u)) in out.iter().zip(&state.u_n).enumerate() {
                let [nx, ny] = mesh.normal_edge[e];
                let projected = state.u_cell.iter().map(|v| v[0] * nx + v[1] * ny);
                let (lo, hi) = projected.fold((u, u), |(lo, hi), p| (lo.min(p), hi.max(p)));
                assert!(
                    lo - 1e-12 <= f && f <= hi + 1e-12,
                    "dt = {dt}: edge {e}: {f} outside [{lo}, {hi}]"
                );
            }
        }
    }

    /// `advection = "none"` では移流項を省く。
    #[test]
    fn advection_can_be_disabled() {
        let mesh = mesh(20.0, 4.0, 1.0, ChannelLayout::Hexagonal);
        let mut initial = State::new(mesh.n_cells(), mesh.n_edges());
        for c in 0..mesh.n_cells() {
            initial.eta[c] = 1.0 + 0.1 * (0.3 * mesh.x_cell[c]).sin();
        }
        let run = |scheme: AdvectionScheme| {
            let mut state = initial.clone();
            let mut solver = HydroSolver::new(
                &mesh,
                HydroParams {
                    advection: scheme,
                    ..params(0.0)
                },
            );
            let mut report = StepReport::default();
            for _ in 0..5 {
                report = solver.step(&mesh, &mut state, &closed(), 0.2).unwrap();
            }
            (state, report)
        };
        let (with, report_with) = run(AdvectionScheme::Upwind);
        let (without, report_without) = run(AdvectionScheme::Disabled);
        assert!(report_with.max_courant > 0.0);
        assert_eq!(report_without.max_courant, 0.0);
        assert_ne!(with.u_n, without.u_n);
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

    /// 閉じた水路の静振（摩擦・移流なし、1周期 50 ステップ）：予測子・修正子の1周期後の振幅は、
    /// θ 法の1段とほぼ同じ（重力波を余分に減衰させない）。1ステップ全体を SSP-RK2 にした形
    /// （ADR-0016）では、1次元の線形解析で、この条件の振幅が約 0.65 倍になる。
    #[test]
    fn predictor_corrector_does_not_damp_seiches() {
        let (length, depth, amplitude): (f64, f64, f64) = (20.0, 1.0, 0.01);
        let mesh = mesh(length, 1.0, 0.5, ChannelLayout::Rectangular);
        let period = 2.0 * length / (9.81 * depth).sqrt();
        let steps = 50;
        let dt = period / steps as f64;
        let mut initial = State::new(mesh.n_cells(), mesh.n_edges());
        for c in 0..mesh.n_cells() {
            initial.eta[c] = depth + amplitude * (PI * mesh.x_cell[c] / length).cos();
        }
        let first = (0..mesh.n_cells())
            .min_by(|&a, &b| mesh.x_cell[a].total_cmp(&mesh.x_cell[b]))
            .expect("cells");
        let run = |time_integration| {
            let params = HydroParams {
                advection: AdvectionScheme::Disabled,
                time_integration,
                ..params(0.0)
            };
            let mut solver = HydroSolver::new(&mesh, params);
            let mut state = initial.clone();
            for _ in 0..steps {
                solver.step(&mesh, &mut state, &closed(), dt).unwrap();
            }
            (state.eta[first] - depth) / (initial.eta[first] - depth)
        };
        let euler = run(TimeIntegration::Euler);
        let predictor_corrector = run(TimeIntegration::PredictorCorrector);
        assert!(euler > 0.9, "euler {euler}");
        assert!(
            (predictor_corrector - euler).abs() <= 0.02,
            "predictor-corrector {predictor_corrector}, euler {euler}"
        );
    }

    /// 予測子・修正子でも、セルごとに V^{n+1} − V^n = −Δt Σ_j s_ij L_j q_j^θ が Newton 法の判定値の
    /// 範囲で成り立つ（連続式は修正子の1回の流量で閉じる）。冠水と干出を含む振動で、1次と2次の
    /// 移流について確かめる。
    #[test]
    fn predictor_corrector_volume_change_matches_the_theta_flux() {
        let mesh = mesh(20.0, 4.0, 1.0, ChannelLayout::Hexagonal);
        let mut initial = State::new(mesh.n_cells(), mesh.n_edges());
        for c in 0..mesh.n_cells() {
            let x = mesh.x_cell[c];
            initial.zb[c] = 0.02 * (x - 10.0).powi(2);
            let level = if x < 10.0 { 1.0 } else { 0.5 };
            initial.eta[c] = initial.zb[c].max(level);
        }
        let dt = 0.2;
        for advection in [AdvectionScheme::Upwind, AdvectionScheme::SecondOrder] {
            let params = HydroParams {
                advection,
                time_integration: TimeIntegration::PredictorCorrector,
                ..params(0.02)
            };
            let (tol, dry_depth) = (params.solver.newton_tol, params.dry_depth);
            let mut solver = HydroSolver::new(&mesh, params);
            let mut state = initial.clone();
            let cell_volume =
                |state: &State, c: usize| volume(mesh.area_cell[c], state.eta[c], state.zb[c]);
            let (mut wetted, mut dried) = (false, false);
            for _ in 0..200 {
                let before: Vec<f64> = (0..mesh.n_cells())
                    .map(|c| cell_volume(&state, c))
                    .collect();
                let report = solver.step(&mesh, &mut state, &closed(), dt).unwrap();
                assert!(report.residual <= tol, "residual {:e}", report.residual);
                let flux = solver.theta_flux();
                for (c, &v0) in before.iter().enumerate() {
                    let outflow: f64 = mesh.edges_on_cell[c]
                        .iter()
                        .zip(&mesh.edge_sign_on_cell[c])
                        .map(|(&e, &s)| s * mesh.dv_edge[e] * flux[e])
                        .sum();
                    let v = cell_volume(&state, c);
                    let error = (v - v0 + dt * outflow).abs();
                    assert!(
                        error <= 1.01 * tol * mesh.area_cell[c] + 1e-13,
                        "{advection:?}: cell {c}: {error:e}"
                    );
                    let (h0, h) = (v0 / mesh.area_cell[c], state.depth(c));
                    wetted |= h0 == 0.0 && h > 0.0;
                    dried |= h0 >= dry_depth && h < dry_depth;
                }
                assert!((0..mesh.n_cells()).all(|c| state.eta[c] >= state.zb[c]));
            }
            assert!(
                wetted && dried,
                "{advection:?}: wetting {wetted}, drying {dried}"
            );
            let (v, v0) = (total_volume(&mesh, &state), total_volume(&mesh, &initial));
            assert!(
                ((v - v0) / v0).abs() <= 1e-12,
                "{advection:?}: volume {v} != {v0}"
            );
        }
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
