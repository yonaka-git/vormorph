//! ボロノイ格子のデータ構造。
//!
//! 接続配列と幾何量は、MPAS のメッシュ仕様（MPAS Mesh Specification v1.0, 2015;
//! Ringler et al., 2010）の名前と意味に合わせ、Rust の命名規則に従って snake_case にしている。
//!
//! ## 配列の一覧
//!
//! | フィールド | MPAS 名 | 位置 | 意味 |
//! |---|---|---|---|
//! | `x_cell`, `y_cell` | `xCell`, `yCell` | セル | 母点の座標 \[m\] |
//! | `area_cell` | `areaCell` | セル | セル面積 A_i \[m²\] |
//! | `edges_on_cell` | `edgesOnCell` | セル | セルを囲む辺（反時計回り） |
//! | `cells_on_cell` | `cellsOnCell` | セル | 辺を挟んだ隣のセル（境界辺では [`NONE`]） |
//! | `vertices_on_cell` | `verticesOnCell` | セル | セルの頂点（反時計回り） |
//! | `edge_sign_on_cell` | `edgeSignOnCell`（符号が逆、下記） | セル | 辺の法線がセルから外向きなら +1 |
//! | `cells_on_edge` | `cellsOnEdge` | 辺 | 辺を挟む2セル。法線は `[0]` から `[1]` へ向く |
//! | `vertices_on_edge` | `verticesOnEdge` | 辺 | 辺の両端の頂点 |
//! | `dc_edge` | `dcEdge` | 辺 | 母点間の距離 δ_j \[m\] |
//! | `dv_edge` | `dvEdge` | 辺 | 辺（ボロノイ辺）の長さ L_j \[m\] |
//! | `x_edge`, `y_edge` | （定義が異なる、下記） | 辺 | 辺の中点 \[m\] |
//! | `normal_edge` | — | 辺 | 単位法線ベクトル n_j |
//! | `edge_kind` | — | 辺 | 内部・壁・上流端・下流端の別 |
//! | `x_vertex`, `y_vertex` | `xVertex`, `yVertex` | 頂点 | ボロノイ頂点の座標 \[m\] |
//!
//! ## 規約
//!
//! - 辺 e の法線 n_e は `cells_on_edge[e][0]` から外向き（`[0]` → `[1]` の向き）。
//!   法線流速 u_n もこの向きを正とする。
//! - `edge_sign_on_cell[c][k]` は、辺 `edges_on_cell[c][k]` の法線がセル `c` から見て外向きなら
//!   `+1.0`、内向きなら `-1.0`。したがってセル c から外へ出る流量は
//!   Σ_k `edge_sign_on_cell[c][k]` × L × H × u_n。
//!   **MPAS-Ocean の `edgeSignOnCell` は `cellsOnEdge(1)` 側を −1 とする逆符号の定義**なので、
//!   MPAS 形式のデータを読み込むときは符号を反転する。
//! - 境界辺では `cells_on_edge[e] = [c, NONE]`（領域の内側が `[0]`）で、法線は領域の外を向く。
//! - `edges_on_cell[c][k]` は頂点 `vertices_on_cell[c][k]` と `vertices_on_cell[c][(k + 1) % n]`
//!   を結ぶ（この並びはこのクレートの約束。外部データを読むときに合わせる）。
//! - `dc_edge` は、境界辺では母点を辺に関して鏡映した点までの距離（母点から辺までの距離の
//!   2倍）とする。下流端の水位境界で水位勾配を評価するために使う。
//! - `x_edge`, `y_edge` は辺の両端の頂点の中点。Perot 法の再構成が一様流で厳密になるのは
//!   この点を使ったときなので、母点間の線分と辺の交点（MPAS の `xEdge` の定義に近い）とは
//!   区別する。
//!
//! 接続配列は行ごとに長さが異なるので、CSR 形式の [`Jagged`] で持つ。

pub mod channel;
pub mod voronoi;

use std::ops::{Index, IndexMut};

/// 領域の外側を表すセル番号。
pub const NONE: usize = usize::MAX;

/// 行ごとに長さが異なる2次元配列（CSR 形式）。
///
/// `mesh.edges_on_cell[c]` のように行番号で添字アクセスすると、その行のスライスが返る。
/// 全行を1本の `Vec` に詰めて持つので、セル単位のループでメモリアクセスが連続する。
#[derive(Debug, Clone, PartialEq)]
pub struct Jagged<T> {
    offsets: Vec<usize>,
    data: Vec<T>,
}

impl<T> Jagged<T> {
    /// 0 行の配列を作る。
    pub fn new() -> Self {
        Self {
            offsets: vec![0],
            data: Vec::new(),
        }
    }

    /// 行の並びから作る。
    pub fn from_rows<R, I>(rows: R) -> Self
    where
        R: IntoIterator<Item = I>,
        I: IntoIterator<Item = T>,
    {
        let mut out = Self::new();
        for row in rows {
            out.push_row(row);
        }
        out
    }

    /// 末尾に1行追加する。
    pub fn push_row<I: IntoIterator<Item = T>>(&mut self, row: I) {
        self.data.extend(row);
        self.offsets.push(self.data.len());
    }

    /// 行数。
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    /// 行が1つもなければ `true`。
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// 全要素数。
    pub fn total_len(&self) -> usize {
        self.data.len()
    }

    /// `i` 行目の要素。
    pub fn row(&self, i: usize) -> &[T] {
        &self.data[self.offsets[i]..self.offsets[i + 1]]
    }

    /// `i` 行目の要素（書き換え可能）。
    pub fn row_mut(&mut self, i: usize) -> &mut [T] {
        let (start, end) = (self.offsets[i], self.offsets[i + 1]);
        &mut self.data[start..end]
    }

    /// 行を先頭から順に返す。
    pub fn rows(&self) -> impl ExactSizeIterator<Item = &[T]> {
        self.offsets.windows(2).map(|w| &self.data[w[0]..w[1]])
    }
}

impl<T> Default for Jagged<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Index<usize> for Jagged<T> {
    type Output = [T];

    fn index(&self, i: usize) -> &[T] {
        self.row(i)
    }
}

impl<T> IndexMut<usize> for Jagged<T> {
    fn index_mut(&mut self, i: usize) -> &mut [T] {
        self.row_mut(i)
    }
}

/// 辺の種類。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EdgeKind {
    /// 2つのセルに挟まれた内部の辺。
    #[default]
    Interior,
    /// 壁。法線流速 0、流砂の出入りなし。
    Wall,
    /// 上流端（流量境界）。
    Upstream,
    /// 下流端（水位・水深境界）。
    Downstream,
}

/// ボロノイ格子。
///
/// フィールドの意味と規約は[モジュールの説明](self)を参照。
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Mesh {
    /// 母点の x 座標 \[m\]（MPAS: `xCell`）。
    pub x_cell: Vec<f64>,
    /// 母点の y 座標 \[m\]（MPAS: `yCell`）。
    pub y_cell: Vec<f64>,
    /// セル面積 A_i \[m²\]（MPAS: `areaCell`）。
    pub area_cell: Vec<f64>,
    /// セルを囲む辺、反時計回り（MPAS: `edgesOnCell`）。
    pub edges_on_cell: Jagged<usize>,
    /// 辺を挟んだ隣のセル、境界辺では [`NONE`]（MPAS: `cellsOnCell`）。
    pub cells_on_cell: Jagged<usize>,
    /// セルの頂点、反時計回り（MPAS: `verticesOnCell`）。
    pub vertices_on_cell: Jagged<usize>,
    /// 辺の法線がセルから外向きなら +1、内向きなら −1（MPAS-Ocean の `edgeSignOnCell` と逆符号）。
    pub edge_sign_on_cell: Jagged<f64>,
    /// 辺を挟む2セル。法線は `[0]` から `[1]` へ向く（MPAS: `cellsOnEdge`）。
    pub cells_on_edge: Vec<[usize; 2]>,
    /// 辺の両端の頂点（MPAS: `verticesOnEdge`）。
    pub vertices_on_edge: Vec<[usize; 2]>,
    /// 母点間の距離 δ_j \[m\]。境界辺では鏡映点までの距離（MPAS: `dcEdge`）。
    pub dc_edge: Vec<f64>,
    /// 辺の長さ L_j \[m\]（MPAS: `dvEdge`）。
    pub dv_edge: Vec<f64>,
    /// 辺の中点の x 座標 \[m\]。
    pub x_edge: Vec<f64>,
    /// 辺の中点の y 座標 \[m\]。
    pub y_edge: Vec<f64>,
    /// 辺の単位法線ベクトル n_j（`cells_on_edge[e][0]` から外向き）。
    pub normal_edge: Vec<[f64; 2]>,
    /// 辺の種類。
    pub edge_kind: Vec<EdgeKind>,
    /// ボロノイ頂点の x 座標 \[m\]（MPAS: `xVertex`）。
    pub x_vertex: Vec<f64>,
    /// ボロノイ頂点の y 座標 \[m\]（MPAS: `yVertex`）。
    pub y_vertex: Vec<f64>,
}

impl Mesh {
    /// セル数。
    pub fn n_cells(&self) -> usize {
        self.x_cell.len()
    }

    /// 辺の数。
    pub fn n_edges(&self) -> usize {
        self.cells_on_edge.len()
    }

    /// 頂点の数。
    pub fn n_vertices(&self) -> usize {
        self.x_vertex.len()
    }

    /// 辺 `e` が境界辺なら `true`。
    pub fn is_boundary_edge(&self, e: usize) -> bool {
        self.cells_on_edge[e][1] == NONE
    }

    /// 辺 `e` を挟んでセル `c` の反対側にあるセル。境界辺なら [`NONE`]。
    pub fn other_cell(&self, e: usize, c: usize) -> usize {
        let [c0, c1] = self.cells_on_edge[e];
        debug_assert!(c == c0 || c == c1, "cell {c} is not on edge {e}");
        if c == c0 { c1 } else { c0 }
    }

    /// 辺 `e` の法線がセル `c` から見て外向きなら `+1.0`、内向きなら `-1.0`。
    pub fn edge_sign(&self, e: usize, c: usize) -> f64 {
        if self.cells_on_edge[e][0] == c {
            1.0
        } else {
            -1.0
        }
    }

    /// `cells_on_edge` と `edges_on_cell` から、`cells_on_cell` と `edge_sign_on_cell` を作り直す。
    ///
    /// 符号の規約をコードで定めている箇所なので、格子の生成・読み込みでは必ずこれを通す。
    pub fn derive_cell_connectivity(&mut self) {
        let mut cells_on_cell = Jagged::new();
        let mut edge_sign_on_cell = Jagged::new();
        for (c, edges) in self.edges_on_cell.rows().enumerate() {
            cells_on_cell.push_row(edges.iter().map(|&e| self.other_cell(e, c)));
            edge_sign_on_cell.push_row(edges.iter().map(|&e| self.edge_sign(e, c)));
        }
        self.cells_on_cell = cells_on_cell;
        self.edge_sign_on_cell = edge_sign_on_cell;
    }

    /// 格子の幾何を検査する（検証ケース V01）。
    pub fn check_geometry(&self) -> GeometryReport {
        todo!("総面積・直交誤差・閉合誤差の計算")
    }
}

/// [`Mesh::check_geometry`] の結果。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct GeometryReport {
    /// 全セル面積の和 \[m²\]。
    pub total_area: f64,
    /// 内部の辺の直交誤差の最大値 max_j |cos(母点間ベクトル, 辺の向き)| \[-\]。
    pub max_orthogonality_error: f64,
    /// セルの閉合誤差の最大値 max_i |Σ_k s_ik L_k n_k| / Σ_k L_k \[-\]。
    pub max_closure_error: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jagged_rows_and_indexing() {
        let mut j = Jagged::from_rows([vec![1, 2, 3], vec![], vec![4]]);
        assert_eq!(j.len(), 3);
        assert_eq!(j.total_len(), 4);
        assert_eq!(&j[0], &[1, 2, 3]);
        assert!(j[1].is_empty());
        j[2][0] = 5;
        assert_eq!(
            j.rows().collect::<Vec<_>>(),
            vec![&[1, 2, 3][..], &[], &[5]]
        );
        assert!(Jagged::<i32>::new().is_empty());
    }

    /// 2つのセルが辺 0 を共有し、セル 0 だけが境界辺 1 を持つ最小の接続で、符号の規約を確かめる。
    #[test]
    fn edge_sign_is_positive_outward() {
        let mut mesh = Mesh {
            cells_on_edge: vec![[0, 1], [0, NONE]],
            edges_on_cell: Jagged::from_rows([vec![0, 1], vec![0]]),
            ..Mesh::default()
        };
        mesh.derive_cell_connectivity();
        assert_eq!(&mesh.edge_sign_on_cell[0], &[1.0, 1.0]);
        assert_eq!(&mesh.edge_sign_on_cell[1], &[-1.0]);
        assert_eq!(&mesh.cells_on_cell[0], &[1, NONE]);
        assert_eq!(&mesh.cells_on_cell[1], &[0]);
        assert!(mesh.is_boundary_edge(1));
        assert!(!mesh.is_boundary_edge(0));
    }
}
