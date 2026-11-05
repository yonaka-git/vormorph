//! 流砂と河床変動：掃流砂、浮遊砂、Exner 式（均一粒径）。
//!
//! 離散式の詳細は `docs/design/overview.md` の第4.6〜4.8節。
//!
//! ## 掃流砂（辺の法線成分）
//!
//! 掃流砂量は、辺 j の法線成分 q_{b,j} を直接評価する。
//!
//! ```text
//! q_{b,j} = √(s g d³) Φ(τ*_j) [ u_j / |u|_j − γ_j (z_{b,2(j)} − z_{b,1(j)}) / δ_j ]
//! γ_j = √(τ*_c / (μ_s μ_k τ*_j))
//! ```
//!
//! γ_j は長谷川の式による河床勾配の補正係数。ボロノイ格子では母点を結ぶ線分が辺に直交
//! するので、法線方向の河床勾配は2点差分で求まる。Φ は [`formulas`] の MPM 式または
//! 芦田・道上式（設定で切り替える）。
//!
//! **非侵食層（床止め）**：あるセルから1ステップで流出する土砂量が、非侵食層より上にある
//! 土砂量 (1 − λ) A_i (z_{b,i} − z_{ne,i}) を超えるときは、そのセルから流出する辺の流砂量を
//! 同じ比率で縮めて頭打ちにする。
//!
//! ## 浮遊砂
//!
//! hC の移流拡散式を、連続式と同じ θ 平均流量 q_j^θ（[`crate::hydro::HydroSolver::theta_flux`]）
//! で解く。沈降 w_f c_b は陰的に扱い、水深が `min_depth` 未満のセルでは浮上量を 0 にする。
//!
//! ## 河床（Exner 式）
//!
//! ```text
//! (1 − λ) A_i Δz_{b,i} = −M Δt [ Σ_j s_ij L_j q_{b,j} + A_i (E_i − w_f c_{b,i}) ]
//! ```
//!
//! M は MORFAC。更新後、河床勾配が水中安息角を超える辺では土砂を再配分する（安息角補正）。
//! 湛水しているセルでは水位 η を保ち（水深が変わる）、乾いたセルでは η = z_b とする。
//!
//! ## 参考文献
//!
//! - Meyer-Peter & Müller (1948), 芦田・道上 (1972): 掃流砂量式
//! - 長谷川 (1984)（式の形は長谷川, 1983 と Nays2DH マニュアルで確認）: 河床勾配の補正
//! - Itakura & Kishi (1980), Garcia & Parker (1991): 浮遊砂の浮上量
//! - Rubey (1933): 沈降速度

pub mod formulas;

use crate::config::{Config, SedimentConfig};
use crate::mesh::Mesh;
use crate::state::State;

/// 流砂の計算に使うパラメータ。物理定数は設定から受け取る。
#[derive(Debug, Clone, PartialEq)]
pub struct SedimentParams {
    /// 重力加速度 g \[m/s²\]。
    pub gravity: f64,
    /// 水の密度 ρ \[kg/m³\]。
    pub water_density: f64,
    /// 水の動粘性係数 ν \[m²/s\]。
    pub kinematic_viscosity: f64,
    /// マニングの粗度係数 n（掃流力の評価に使う）。
    pub manning_n: f64,
    /// wet/dry を判定する水深の閾値 \[m\]。
    pub dry_depth: f64,
    /// 流砂の設定。
    pub sediment: SedimentConfig,
}

impl SedimentParams {
    /// 設定から作る。`[sediment]` が省略されている（固定床）なら `None`。
    pub fn from_config(config: &Config) -> Option<Self> {
        let sediment = config.sediment.clone()?;
        Some(Self {
            gravity: config.physics.gravity,
            water_density: config.physics.water_density,
            kinematic_viscosity: config.physics.kinematic_viscosity,
            manning_n: config.hydro.manning_n,
            dry_depth: config.hydro.dry_depth,
            sediment,
        })
    }

    /// 水中比重 s = ρ_s / ρ − 1 \[-\]。
    pub fn relative_density(&self) -> f64 {
        self.sediment.density / self.water_density - 1.0
    }
}

/// 河床更新の記録。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct BedUpdateReport {
    /// 河床高の変化量の最大値 |Δz_b| \[m\]。
    pub max_change: f64,
    /// 非侵食層で流出量を頭打ちにしたセルの数。
    pub capped_cells: usize,
    /// 安息角補正の反復回数。
    pub avalanche_iterations: usize,
}

/// 流砂と河床変動のソルバー。
#[derive(Debug, Clone)]
#[expect(dead_code, reason = "v0.1.0 の雛形。浮遊砂と河床更新の実装で使う")]
pub struct SedimentSolver {
    params: SedimentParams,
    /// 辺の掃流砂量 q_{b,j}（法線方向、単位幅あたり）\[m²/s\]
    bedload_flux: Vec<f64>,
    /// セルの浮上量 E_i \[m/s\]
    entrainment: Vec<f64>,
    /// セルの沈降量 w_f c_{b,i} \[m/s\]
    deposition: Vec<f64>,
}

impl SedimentSolver {
    /// 格子の大きさに合わせて作業配列を確保する。
    pub fn new(mesh: &Mesh, params: SedimentParams) -> Self {
        Self {
            params,
            bedload_flux: vec![0.0; mesh.n_edges()],
            entrainment: vec![0.0; mesh.n_cells()],
            deposition: vec![0.0; mesh.n_cells()],
        }
    }

    /// パラメータ。
    pub fn params(&self) -> &SedimentParams {
        &self.params
    }

    /// 直前に計算した辺の掃流砂量 q_{b,j} \[m²/s\]。
    pub fn bedload_flux(&self) -> &[f64] {
        &self.bedload_flux
    }

    /// 辺の掃流砂量 q_{b,j} を計算する（勾配補正と非侵食層の頭打ちを含む）。
    pub fn compute_bedload(&mut self, mesh: &Mesh, state: &State, dt: f64) {
        let _ = (mesh, state, dt);
        todo!("辺の法線方向の掃流砂量")
    }

    /// 浮遊砂 hC を `dt` だけ進める。`theta_flux` は連続式と同じ θ 平均の単位幅流量。
    pub fn step_suspended(&mut self, mesh: &Mesh, state: &mut State, theta_flux: &[f64], dt: f64) {
        let _ = (mesh, state, theta_flux, dt);
        todo!("hC の移流拡散（沈降は陰的）")
    }

    /// Exner 式で河床を更新し、安息角補正と水位の調整を行う。
    pub fn update_bed(&mut self, mesh: &Mesh, state: &mut State, dt: f64) -> BedUpdateReport {
        let _ = (mesh, state, dt);
        todo!("Exner 式・MORFAC・安息角補正")
    }
}
