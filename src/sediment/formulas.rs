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
    if tau_star > tau_star_c {
        8.0 * (tau_star - tau_star_c).powf(1.5)
    } else {
        0.0
    }
}

/// 芦田・道上 (1972) の掃流砂量式。
///
/// Φ = 17 τ*_e^{3/2} (1 − τ*_c / τ*) (1 − √(τ*_c / τ*))（τ* ≤ τ*_c では 0）。
/// 係数 17 は原著の β/μ_f（実験資料から定めた値）。τ*_e は有効掃流力で、河床波による形状抵抗を
/// 除いた掃流力。この実装では河床波を考えないので τ*_e = τ* とする（`docs/design/overview.md`
/// の第2.2節）。
pub fn ashida_michiue(tau_star: f64, tau_star_c: f64) -> f64 {
    if tau_star > tau_star_c {
        let ratio = tau_star_c / tau_star;
        17.0 * tau_star.powf(1.5) * (1.0 - ratio) * (1.0 - ratio.sqrt())
    } else {
        0.0
    }
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
/// μ_s は静止摩擦係数、μ_k は動摩擦係数（長谷川, 1983 の式 (17) の Γ）。辺の法線方向の
/// 掃流砂量に −γ ∂z_b/∂n の形で加える。掃流砂が動く τ* > τ*_c でだけ使う（τ* = 0 では
/// 無限大になる）。
pub fn slope_correction_coefficient(tau_star: f64, tau_star_c: f64, mu_s: f64, mu_k: f64) -> f64 {
    (tau_star_c / (mu_s * mu_k * tau_star)).sqrt()
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

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * b.abs().max(1.0)
    }

    #[test]
    fn mpm_matches_hand_calculation_and_vanishes_below_the_threshold() {
        // 8 × (0.647 − 0.047)^{3/2} = 8 × 0.6^{3/2}
        assert!(close(
            meyer_peter_muller(0.647, 0.047),
            8.0 * 0.6_f64.powf(1.5)
        ));
        assert!(close(meyer_peter_muller(1.047, 0.047), 8.0));
        assert_eq!(meyer_peter_muller(0.047, 0.047), 0.0);
        assert_eq!(meyer_peter_muller(0.01, 0.047), 0.0);
        assert_eq!(meyer_peter_muller(0.0, 0.047), 0.0);
        // τ*_c の直上で連続（Φ → 0）
        assert!(meyer_peter_muller(0.047 + 1e-9, 0.047) < 1e-12);
    }

    #[test]
    fn ashida_michiue_matches_hand_calculation_and_vanishes_below_the_threshold() {
        // τ* = 4 τ*_c：17 τ*^{3/2} (1 − 1/4)(1 − 1/2)
        let (tau_c, tau) = (0.05, 0.2);
        let expected = 17.0 * 0.2_f64.powf(1.5) * 0.75 * 0.5;
        assert!(close(ashida_michiue(tau, tau_c), expected));
        assert_eq!(ashida_michiue(tau_c, tau_c), 0.0);
        assert_eq!(ashida_michiue(0.0, tau_c), 0.0);
        assert!(ashida_michiue(tau_c * (1.0 + 1e-9), tau_c) < 1e-12);
        // 掃流力が大きいと 17 τ*^{3/2} に近づく
        let tau = 100.0;
        assert!((ashida_michiue(tau, tau_c) / (17.0 * tau.powf(1.5)) - 1.0).abs() < 0.03);
    }

    #[test]
    fn bedload_intensity_dispatches_on_the_formula() {
        assert_eq!(
            bedload_intensity(BedloadFormula::Mpm, 0.5, 0.047),
            meyer_peter_muller(0.5, 0.047)
        );
        assert_eq!(
            bedload_intensity(BedloadFormula::AshidaMichiue, 0.5, 0.05),
            ashida_michiue(0.5, 0.05)
        );
    }

    #[test]
    fn slope_coefficient_follows_hasegawa() {
        // τ* = τ*_c では 1/√(μ_s μ_k)
        assert!(close(
            slope_correction_coefficient(0.047, 0.047, 1.0, 0.45),
            1.0 / 0.45_f64.sqrt()
        ));
        // τ* を 4 倍にすると半分
        assert!(close(
            slope_correction_coefficient(0.188, 0.047, 1.0, 0.45),
            0.5 / 0.45_f64.sqrt()
        ));
    }
}
