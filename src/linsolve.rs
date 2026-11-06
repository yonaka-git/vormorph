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

    /// 長さが `n` でなければ確保し直す。
    fn ensure_len(&mut self, n: usize) {
        if self.r.len() != n {
            *self = Self::new(n);
        }
    }
}

/// Jacobi 前処理付き共役勾配法で A x = b を解く。
///
/// `x` には初期値を入れて渡す。解で上書きされる。
///
/// - `b` が 0 なら、反復せずに x = 0 を返す（相対残差が定義できないため）。
/// - 前処理は対角の逆数。対角成分が 0 の行は前処理で 0 を掛ける（その行が他と結合して
///   いなければ、x の成分は初期値のまま残る）。対角成分が負か有限でなければ、A は正定値で
///   ないので [`PcgError::NotPositiveDefinite`]（`iteration = 0`）を返す。
/// - 収束判定には漸化式で更新した残差を使う。
pub fn pcg<A: LinearOperator + ?Sized>(
    a: &A,
    b: &[f64],
    x: &mut [f64],
    settings: &PcgSettings,
    work: &mut PcgWorkspace,
) -> Result<PcgStats, PcgError> {
    let n = a.dim();
    debug_assert_eq!(b.len(), n, "b must have one value per unknown");
    debug_assert_eq!(x.len(), n, "x must have one value per unknown");
    work.ensure_len(n);
    let PcgWorkspace {
        r,
        z,
        p,
        q,
        inv_diag,
    } = work;

    let b_norm = norm(b);
    if b_norm == 0.0 {
        x.fill(0.0);
        return Ok(PcgStats::default());
    }

    a.diagonal(inv_diag);
    for d in inv_diag.iter_mut() {
        if !(d.is_finite() && *d >= 0.0) {
            return Err(PcgError::NotPositiveDefinite { iteration: 0 });
        }
        *d = if *d > 0.0 { 1.0 / *d } else { 0.0 };
    }

    // r = b − A x
    a.apply(x, q);
    for ((ri, bi), qi) in r.iter_mut().zip(b).zip(q.iter()) {
        *ri = bi - qi;
    }
    let mut rel_residual = norm(r) / b_norm;
    if rel_residual <= settings.rel_tol {
        return Ok(PcgStats {
            iterations: 0,
            rel_residual,
        });
    }

    precondition(inv_diag, r, z);
    p.copy_from_slice(z);
    let mut rz = dot(r, z);

    for iteration in 1..=settings.max_iter {
        a.apply(p, q);
        let pq = dot(p, q);
        // NaN も正定値でないとみなす
        if pq.is_nan() || pq <= 0.0 {
            return Err(PcgError::NotPositiveDefinite { iteration });
        }
        let alpha = rz / pq;
        for ((xi, ri), (pi, qi)) in x.iter_mut().zip(r.iter_mut()).zip(p.iter().zip(q.iter())) {
            *xi += alpha * pi;
            *ri -= alpha * qi;
        }
        rel_residual = norm(r) / b_norm;
        if rel_residual <= settings.rel_tol {
            return Ok(PcgStats {
                iterations: iteration,
                rel_residual,
            });
        }
        precondition(inv_diag, r, z);
        let rz_next = dot(r, z);
        let beta = rz_next / rz;
        rz = rz_next;
        for (pi, zi) in p.iter_mut().zip(z.iter()) {
            *pi = zi + beta * *pi;
        }
    }
    Err(PcgError::NotConverged(PcgStats {
        iterations: settings.max_iter,
        rel_residual,
    }))
}

/// z = M⁻¹ r（Jacobi 前処理）。
fn precondition(inv_diag: &[f64], r: &[f64], z: &mut [f64]) {
    for ((zi, ri), di) in z.iter_mut().zip(r).zip(inv_diag) {
        *zi = di * ri;
    }
}

/// 内積。
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// 2 ノルム。
fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 密行列で表した線形作用素（テスト用）。
    struct Dense(Vec<Vec<f64>>);

    impl LinearOperator for Dense {
        fn dim(&self) -> usize {
            self.0.len()
        }

        fn apply(&self, x: &[f64], y: &mut [f64]) {
            for (yi, row) in y.iter_mut().zip(&self.0) {
                *yi = dot(row, x);
            }
        }

        fn diagonal(&self, d: &mut [f64]) {
            for (i, di) in d.iter_mut().enumerate() {
                *di = self.0[i][i];
            }
        }
    }

    /// 対角が場所によって異なる、1次元のラプラシアン＋対角行列（水位方程式と同じ形）。
    fn laplacian_plus_diagonal(n: usize) -> Dense {
        let mut a = vec![vec![0.0; n]; n];
        for i in 0..n {
            let w = 1.0 + (i % 7) as f64; // 辺の重み
            a[i][i] += 0.1 + (i % 3) as f64; // P
            if i + 1 < n {
                a[i][i] += w;
                a[i + 1][i + 1] += w;
                a[i][i + 1] -= w;
                a[i + 1][i] -= w;
            }
        }
        Dense(a)
    }

    fn settings(rel_tol: f64, max_iter: usize) -> PcgSettings {
        PcgSettings { rel_tol, max_iter }
    }

    #[test]
    fn solves_a_symmetric_positive_definite_system() {
        let n = 60;
        let a = laplacian_plus_diagonal(n);
        let x_true: Vec<f64> = (0..n).map(|i| (0.3 * i as f64).sin() + 2.0).collect();
        let mut b = vec![0.0; n];
        a.apply(&x_true, &mut b);

        let mut x = vec![0.0; n];
        let mut work = PcgWorkspace::new(n);
        let stats = pcg(&a, &b, &mut x, &settings(1e-14, 500), &mut work).expect("must converge");
        assert!(stats.iterations > 0 && stats.iterations <= n);
        assert!(stats.rel_residual <= 1e-14);
        for (xi, ti) in x.iter().zip(&x_true) {
            assert!((xi - ti).abs() <= 1e-11, "{xi} != {ti}");
        }
    }

    #[test]
    fn zero_right_hand_side_gives_zero() {
        let a = laplacian_plus_diagonal(5);
        let mut x = vec![1.0; 5];
        let mut work = PcgWorkspace::default();
        let stats = pcg(&a, &[0.0; 5], &mut x, &settings(1e-12, 10), &mut work).unwrap();
        assert_eq!(stats, PcgStats::default());
        assert_eq!(x, vec![0.0; 5]);
    }

    #[test]
    fn exact_initial_guess_needs_no_iteration() {
        let a = laplacian_plus_diagonal(4);
        let x0 = vec![1.0, -2.0, 0.5, 3.0];
        let mut b = vec![0.0; 4];
        a.apply(&x0, &mut b);
        let mut x = x0.clone();
        let stats = pcg(
            &a,
            &b,
            &mut x,
            &settings(1e-12, 10),
            &mut PcgWorkspace::new(4),
        )
        .unwrap();
        assert_eq!(stats.iterations, 0);
        assert_eq!(x, x0);
    }

    #[test]
    fn indefinite_matrix_is_detected() {
        // 対角は正だが固有値 3 と −1 を持つ。
        let a = Dense(vec![vec![1.0, 2.0], vec![2.0, 1.0]]);
        let mut x = vec![0.0; 2];
        let err = pcg(
            &a,
            &[1.0, -1.0],
            &mut x,
            &settings(1e-12, 10),
            &mut PcgWorkspace::new(2),
        );
        assert_eq!(err, Err(PcgError::NotPositiveDefinite { iteration: 1 }));

        let a = Dense(vec![vec![1.0, 0.0], vec![0.0, -1.0]]);
        let err = pcg(
            &a,
            &[1.0, 1.0],
            &mut x,
            &settings(1e-12, 10),
            &mut PcgWorkspace::new(2),
        );
        assert_eq!(err, Err(PcgError::NotPositiveDefinite { iteration: 0 }));
    }

    #[test]
    fn iteration_limit_is_reported() {
        let n = 30;
        let a = laplacian_plus_diagonal(n);
        let b: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let mut x = vec![0.0; n];
        match pcg(
            &a,
            &b,
            &mut x,
            &settings(1e-14, 2),
            &mut PcgWorkspace::new(n),
        ) {
            Err(PcgError::NotConverged(stats)) => {
                assert_eq!(stats.iterations, 2);
                assert!(stats.rel_residual > 1e-14);
            }
            other => panic!("expected NotConverged, got {other:?}"),
        }
    }

    #[test]
    fn workspace_of_another_size_is_reallocated() {
        let a = laplacian_plus_diagonal(8);
        let b = vec![1.0; 8];
        let mut x = vec![0.0; 8];
        let mut work = PcgWorkspace::new(3);
        pcg(&a, &b, &mut x, &settings(1e-12, 100), &mut work).expect("must converge");
        let mut ax = vec![0.0; 8];
        a.apply(&x, &mut ax);
        assert!(ax.iter().zip(&b).all(|(u, v)| (u - v).abs() <= 1e-10));
    }
}
