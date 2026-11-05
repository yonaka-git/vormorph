//! 行列フリーの前処理付き共役勾配法（PCG）。
//!
//! 水位の Newton 反復の各段で現れる線形系 (P + T) x = r を解く。P は wet セルの面積を
//! 並べた対角行列、T は辺の重み w_j ≥ 0 を持つグラフラプラシアン型の対称行列で、
//! P + T は、各連結成分に wet セルか水位境界が1つでもあれば対称正定値になる
//! （`docs/design/overview.md` の第4.3節）。
//!
//! 行列は組み立てず、[`LinearOperator::apply`] で y = A x を計算する。`apply` はセル単位の
//! ループで書き、各セルが自分の y_i だけを書き込むようにする（将来 rayon で並列化するとき、
//! 書き込みの競合が起きないようにするため）。前処理は対角スケーリング（Jacobi）から始める。
//!
//! 収束判定は相対残差 ‖r_k‖₂ / ‖b‖₂ ≤ `rel_tol`。
//!
//! ## 参考文献
//!
//! - Casulli & Zanolli (2012): 区分線形な非線形系に対する Newton 法と、内側の線形反復

use crate::config::SolverConfig;
use std::error::Error;
use std::fmt;

/// 行列を組み立てずに y = A x を計算できる線形作用素。A は対称正定値を想定する。
pub trait LinearOperator {
    /// 未知数の数。
    fn dim(&self) -> usize;

    /// y = A x を計算する。`x` と `y` の長さは [`LinearOperator::dim`]。
    fn apply(&self, x: &[f64], y: &mut [f64]);

    /// A の対角成分を `d` に書く（Jacobi 前処理に使う）。
    fn diagonal(&self, d: &mut [f64]);
}

/// PCG の設定。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcgSettings {
    /// 相対残差の収束判定値。
    pub rel_tol: f64,
    /// 最大反復回数。
    pub max_iter: usize,
}

impl PcgSettings {
    /// 設定ファイルの `[hydro.solver]` から作る。
    pub fn from_config(solver: &SolverConfig) -> Self {
        Self {
            rel_tol: solver.pcg_rel_tol,
            max_iter: solver.pcg_max_iter,
        }
    }
}

/// PCG の結果。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PcgStats {
    /// 反復回数。
    pub iterations: usize,
    /// 最終的な相対残差。
    pub rel_residual: f64,
}

/// PCG の失敗。
#[derive(Debug, Clone, PartialEq)]
pub enum PcgError {
    /// 最大反復回数までに収束しなかった。
    NotConverged(PcgStats),
    /// pᵀAp ≤ 0 となった（係数行列が正定値でない。wet/dry の判定か境界条件の誤りを疑う）。
    NotPositiveDefinite {
        /// 発生した反復の番号。
        iteration: usize,
    },
}

impl fmt::Display for PcgError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotConverged(s) => write!(
                f,
                "PCG did not converge in {} iterations (relative residual {:e})",
                s.iterations, s.rel_residual
            ),
            Self::NotPositiveDefinite { iteration } => {
                write!(f, "PCG breakdown: p^T A p <= 0 at iteration {iteration}")
            }
        }
    }
}

impl Error for PcgError {}

/// PCG の作業配列。ステップごとの確保を避けるため、ソルバーが持ち回す。
#[derive(Debug, Clone, Default)]
#[expect(dead_code, reason = "v0.1.0 の雛形。pcg の実装で使う")]
pub struct PcgWorkspace {
    /// 残差 r
    r: Vec<f64>,
    /// 前処理後の残差 z = M⁻¹ r
    z: Vec<f64>,
    /// 探索方向 p
    p: Vec<f64>,
    /// A p
    q: Vec<f64>,
    /// 対角の逆数 M⁻¹
    inv_diag: Vec<f64>,
}

impl PcgWorkspace {
    /// 未知数 `n` 個分の作業配列を確保する。
    pub fn new(n: usize) -> Self {
        Self {
            r: vec![0.0; n],
            z: vec![0.0; n],
            p: vec![0.0; n],
            q: vec![0.0; n],
            inv_diag: vec![0.0; n],
        }
    }
}

/// Jacobi 前処理付き共役勾配法で A x = b を解く。
///
/// `x` には初期値を入れて渡す（Newton 反復では前の反復の解を使う）。解で上書きされる。
pub fn pcg<A: LinearOperator + ?Sized>(
    a: &A,
    b: &[f64],
    x: &mut [f64],
    settings: &PcgSettings,
    work: &mut PcgWorkspace,
) -> Result<PcgStats, PcgError> {
    let _ = (a, b, x, settings, work);
    todo!("Jacobi 前処理付き共役勾配法")
}
