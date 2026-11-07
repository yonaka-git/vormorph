//! セル中心の値の勾配（最小二乗）。
//!
//! セル c の勾配 ∇φ_c を、隣のセルとの差 φ_k − φ_c の最小二乗で求める。
//!
//! ```text
//! M = Σ_k d_k d_kᵀ,   ∇φ_c = M⁻¹ Σ_k d_k (φ_k − φ_c) = Σ_k w_ck (φ_k − φ_c),   w_ck = M⁻¹ d_k
//! ```
//!
//! - d_k = x_k − x_c は母点の差。周期境界の継ぎ目の向こうのセルは、平行移動した位置
//!   （[`Mesh::neighbor_position_from`]）で数える。境界の外（[`NONE`]）は数えない。
//! - 1次式の場なら厳密。重み w_ck は格子だけで決まるので、[`CellGradient`] は最初に1回だけ
//!   計算して持ち回す。
//! - 隣のセルが一直線に並ぶ（幅が1セルの水路など）ときは、M の階数が 1 になる。その向きの
//!   成分だけを求め、直交する成分は 0 とする。
//!
//! 使う所：開境界の辺の河床勾配（[`crate::sediment::bed_gradient`]）、2次の移流の流速の勾配
//! （[`crate::hydro::advection_second_order`]）、二次流の遠心加速度
//! （[`crate::sediment::secondary_flow`]）。

use super::{Jagged, Mesh, NONE};

/// セル `c` の最小二乗の重み w_ck を、`edges_on_cell[c]`（`cells_on_cell[c]`）の並びで返す。
///
/// 境界辺（隣が [`NONE`]）の重みは 0。隣のセルがなければすべて 0。
pub fn least_squares_weights(mesh: &Mesh, c: usize) -> Vec<[f64; 2]> {
    let (xc, yc) = (mesh.x_cell[c], mesh.y_cell[c]);
    let offsets: Vec<Option<[f64; 2]>> = mesh.edges_on_cell[c]
        .iter()
        .zip(&mesh.cells_on_cell[c])
        .map(|(&e, &k)| {
            (k != NONE).then(|| {
                let [xk, yk] = mesh.neighbor_position_from(e, c);
                [xk - xc, yk - yc]
            })
        })
        .collect();
    let (mut mxx, mut mxy, mut myy) = (0.0, 0.0, 0.0);
    for [dx, dy] in offsets.iter().flatten() {
        mxx += dx * dx;
        mxy += dx * dy;
        myy += dy * dy;
    }
    let trace = mxx + myy;
    let det = mxx * myy - mxy * mxy;
    offsets
        .iter()
        .map(|offset| match *offset {
            None => [0.0, 0.0],
            Some(_) if trace <= 0.0 => [0.0, 0.0],
            Some([dx, dy]) if det > 1.0e-10 * trace * trace => {
                [(myy * dx - mxy * dy) / det, (mxx * dy - mxy * dx) / det]
            }
            Some([dx, dy]) => {
                // 階数 1：M = trace v vᵀ（v は隣のセルの並ぶ向き）なので、勾配は v 方向の成分だけ
                let (vx, vy) = if mxx >= myy { (mxx, mxy) } else { (mxy, myy) };
                let norm = vx.hypot(vy);
                let (vx, vy) = (vx / norm, vy / norm);
                let along = (vx * dx + vy * dy) / trace;
                [along * vx, along * vy]
            }
        })
        .collect()
}

/// 全セルの最小二乗の重み（[`least_squares_weights`]）。
#[derive(Debug, Clone, PartialEq)]
pub struct CellGradient {
    weights: Jagged<[f64; 2]>,
}

impl CellGradient {
    /// 格子の全セルについて重みを計算する。
    pub fn new(mesh: &Mesh) -> Self {
        Self {
            weights: Jagged::from_rows((0..mesh.n_cells()).map(|c| least_squares_weights(mesh, c))),
        }
    }

    /// セル `c` の重み w_ck \[1/m\]（`edges_on_cell[c]` の並び。境界辺は 0）。
    pub fn weights(&self, c: usize) -> &[[f64; 2]] {
        self.weights.row(c)
    }

    /// セル中心の値 `values` の、セル `c` での勾配 ∇φ_c。
    pub fn of(&self, mesh: &Mesh, values: &[f64], c: usize) -> [f64; 2] {
        let mut g = [0.0, 0.0];
        for (&k, &[wx, wy]) in mesh.cells_on_cell[c].iter().zip(self.weights(c)) {
            if k != NONE {
                let d = values[k] - values[c];
                g[0] += wx * d;
                g[1] += wy * d;
            }
        }
        g
    }

    /// セル中心のベクトル `values` の、セル `c` での勾配 ∂v_a/∂x_b（`[a][b]`）。
    ///
    /// 成分ごとに [`CellGradient::of`] と同じ和をとる。
    pub fn of_vector(&self, mesh: &Mesh, values: &[[f64; 2]], c: usize) -> [[f64; 2]; 2] {
        let mut g = [[0.0; 2]; 2];
        let [ux, uy] = values[c];
        for (&k, &[wx, wy]) in mesh.cells_on_cell[c].iter().zip(self.weights(c)) {
            if k != NONE {
                let [vx, vy] = values[k];
                g[0][0] += wx * (vx - ux);
                g[0][1] += wy * (vx - ux);
                g[1][0] += wx * (vy - uy);
                g[1][1] += wy * (vy - uy);
            }
        }
        g
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::ChannelLayout;
    use crate::mesh::channel::{self, ChannelSpec};

    fn spec(width: f64, layout: ChannelLayout) -> ChannelSpec {
        ChannelSpec {
            length: 8.0,
            width,
            cell_size: 1.0,
            layout,
        }
    }

    /// 1次式の場では厳密。隣が一直線に並ぶ（幅1セル）ときは、その方向の成分だけを返す。
    /// 周期境界の継ぎ目でも、平行移動した位置で数えるので厳密。
    #[test]
    fn gradient_is_exact_for_linear_fields() {
        let cases = [
            (5.0, ChannelLayout::Hexagonal, false),
            (5.0, ChannelLayout::Rectangular, false),
            (1.0, ChannelLayout::Rectangular, false),
            (5.0, ChannelLayout::Hexagonal, true),
        ];
        for (width, layout, periodic) in cases {
            let spec = spec(width, layout);
            let mesh = if periodic {
                channel::build_periodic(&spec)
            } else {
                channel::build(&spec)
            }
            .expect("valid channel");
            let gradient = CellGradient::new(&mesh);
            // 周期境界では x 方向に一様（継ぎ目で値がつながる場）
            let gx = if periodic { 0.0 } else { -0.02 };
            let values: Vec<f64> = (0..mesh.n_cells())
                .map(|c| 0.3 + gx * mesh.x_cell[c] + 0.05 * mesh.y_cell[c])
                .collect();
            // ベクトル：x 成分は上の場、y 成分は −2 倍の場
            let vectors: Vec<[f64; 2]> = values.iter().map(|&v| [v, -2.0 * v]).collect();
            for c in 0..mesh.n_cells() {
                let [x, y] = gradient.of(&mesh, &values, c);
                let expected_y = if width > 1.0 { 0.05 } else { 0.0 };
                assert!((x - gx).abs() < 1e-12, "{layout:?} cell {c}: {x}");
                assert!((y - expected_y).abs() < 1e-12, "{layout:?} cell {c}: {y}");
                let g = gradient.of_vector(&mesh, &vectors, c);
                assert_eq!(g[0], [x, y]);
                assert!((g[1][0] + 2.0 * x).abs() < 1e-15 && (g[1][1] + 2.0 * y).abs() < 1e-15);
            }
        }
    }

    /// 境界辺の重みは 0 で、重みの和 Σ_k w_ck d_kᵀ は単位行列（1次式で厳密になる条件）。
    #[test]
    fn weights_reproduce_the_identity() {
        let mesh = channel::build(&spec(5.0, ChannelLayout::Hexagonal)).expect("valid channel");
        let gradient = CellGradient::new(&mesh);
        for c in 0..mesh.n_cells() {
            let mut m = [[0.0; 2]; 2];
            let rows = mesh.edges_on_cell[c].iter().zip(&mesh.cells_on_cell[c]);
            for ((&e, &k), &[wx, wy]) in rows.zip(gradient.weights(c)) {
                if k == NONE {
                    assert_eq!([wx, wy], [0.0, 0.0]);
                    continue;
                }
                let [xk, yk] = mesh.neighbor_position_from(e, c);
                let (dx, dy) = (xk - mesh.x_cell[c], yk - mesh.y_cell[c]);
                m[0][0] += wx * dx;
                m[0][1] += wx * dy;
                m[1][0] += wy * dx;
                m[1][1] += wy * dy;
            }
            for (i, row) in m.iter().enumerate() {
                for (j, &v) in row.iter().enumerate() {
                    let expected = if i == j { 1.0 } else { 0.0 };
                    assert!((v - expected).abs() < 1e-12, "cell {c}: {m:?}");
                }
            }
        }
    }
}
