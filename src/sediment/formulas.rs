//! 流砂量式などの経験式。
//!
//! 無次元量の定義：
//!
//! - 水中比重 s = ρ_s / ρ − 1
//! - 無次元掃流力 τ* = u*² / (s g d)
//! - 無次元掃流砂量 Φ = q_b / √(s g d³)
//!
//! 関数はすべて純粋関数とし、物理定数は引数で受け取る。式の係数は原著
//! （`docs/references.md`）で確かめてから実装する。

use crate::config::{BedloadFormula, EntrainmentFormula};

/// Meyer-Peter & Müller (1948) の掃流砂量式。
///
/// Φ = 8 (τ* − τ*_c)^{3/2}（τ* ≤ τ*_c では 0）。原著の τ*_c は 0.047。
pub fn meyer_peter_muller(tau_star: f64, tau_star_c: f64) -> f64 {
    todo!("MPM 式 (tau* = {tau_star}, tau*c = {tau_star_c})")
}

/// 芦田・道上 (1972) の掃流砂量式。
///
/// Φ = 17 τ*_e^{3/2} (1 − τ*_c / τ*) (1 − √(τ*_c / τ*))（τ* ≤ τ*_c では 0）。
/// τ*_e は有効掃流力。当面は τ*_e = τ* として扱う予定で、扱いは実装時に決める。
pub fn ashida_michiue(tau_star: f64, tau_star_c: f64) -> f64 {
    todo!("芦田・道上式 (tau* = {tau_star}, tau*c = {tau_star_c})")
}

/// 設定で選んだ掃流砂量式で、無次元掃流砂量 Φ を返す。
pub fn bedload_intensity(formula: BedloadFormula, tau_star: f64, tau_star_c: f64) -> f64 {
    match formula {
        BedloadFormula::Mpm => meyer_peter_muller(tau_star, tau_star_c),
        BedloadFormula::AshidaMichiue => ashida_michiue(tau_star, tau_star_c),
    }
}

/// 長谷川の式による河床勾配の補正係数 γ = √(τ*_c / (μ_s μ_k τ*))。
///
/// μ_s は静止摩擦係数、μ_k は動摩擦係数。辺の法線方向の掃流砂量に
/// −γ ∂z_b/∂n の形で加える。
pub fn slope_correction_coefficient(tau_star: f64, tau_star_c: f64, mu_s: f64, mu_k: f64) -> f64 {
    todo!(
        "長谷川の式の係数 (tau* = {tau_star}, tau*c = {tau_star_c}, mu_s = {mu_s}, mu_k = {mu_k})"
    )
}

/// Rubey (1933) の沈降速度式 \[m/s\]。
///
/// w_f = √(s g d) [ √(2/3 + 36 ν² / (s g d³)) − √(36 ν² / (s g d³)) ]
pub fn rubey_settling_velocity(
    grain_size: f64,
    relative_density: f64,
    gravity: f64,
    nu: f64,
) -> f64 {
    todo!("Rubey 式 (d = {grain_size}, s = {relative_density}, g = {gravity}, nu = {nu})")
}

/// 浮上量式に渡す量。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EntrainmentInput {
    /// 摩擦速度 u* \[m/s\]。
    pub friction_velocity: f64,
    /// 沈降速度 w_f \[m/s\]。
    pub settling_velocity: f64,
    /// 無次元掃流力 τ* \[-\]。
    pub shields: f64,
    /// 水中比重 s \[-\]。
    pub relative_density: f64,
    /// 粒径 d \[m\]。
    pub grain_size: f64,
    /// 重力加速度 g \[m/s²\]。
    pub gravity: f64,
    /// 動粘性係数 ν \[m²/s\]。
    pub kinematic_viscosity: f64,
}

/// 板倉・岸の浮上量式 \[m/s\]（Itakura & Kishi, 1980）。
///
/// 係数と関数 Ω の定義は原著で確かめてから実装する。
pub fn itakura_kishi(input: &EntrainmentInput) -> f64 {
    todo!("板倉・岸の式 ({input:?})")
}

/// Garcia & Parker (1991) の浮上量式 \[m/s\]（均一粒径の場合）。
///
/// E = w_f A Z⁵ / (1 + A Z⁵ / 0.3)、Z = (u* / w_f) Re_p^{0.6}、Re_p = √(s g d) d / ν、
/// A = 1.3 × 10⁻⁷（無次元の浮上率に w_f を掛けて体積フラックスにする）。
pub fn garcia_parker(input: &EntrainmentInput) -> f64 {
    todo!("Garcia & Parker の式 ({input:?})")
}

/// 設定で選んだ浮上量式で、浮上量 E \[m/s\] を返す。
pub fn entrainment_rate(formula: EntrainmentFormula, input: &EntrainmentInput) -> f64 {
    match formula {
        EntrainmentFormula::ItakuraKishi => itakura_kishi(input),
        EntrainmentFormula::GarciaParker => garcia_parker(input),
    }
}
