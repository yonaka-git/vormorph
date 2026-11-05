//! 外部で生成した格子の読み込みと、ボロノイ格子の幾何量の計算。
//!
//! ## 直交性
//!
//! デローニー三角形分割の外心を頂点とするボロノイ図では、隣り合う母点を結ぶ線分と、
//! 両者を隔てるボロノイ辺が直交する。C 格子の水位勾配 (η_2 − η_1)/δ_j や、河床の法線方向
//! 勾配が2点差分で厳密になるのはこの性質による。鈍角三角形では外心が三角形の外に出て、
//! 母点の順序と辺の向きが食い違うことがあるため、格子生成の側で局所的に直交性のよい分割
//! （JIGSAW の locally orthogonal grid; Engwirda, 2017）を作る。
//!
//! ## 入力形式
//!
//! 当面は JIGSAW（jigsawpy）が出力する三角形分割を読み、その双対としてボロノイ格子を
//! 組み立てる想定。具体的なファイル形式、境界の種類（壁・上流端・下流端）の指定方法、
//! 河床高の与え方は未定で、実装時に ADR に記録する。
//!
//! ## 参考文献
//!
//! - Engwirda (2017): JIGSAW-GEO による局所直交格子の生成
//! - MPAS Mesh Specification v1.0 (2015): 接続配列の定義

use super::Mesh;
use std::io;
use std::path::Path;

/// 三角形分割からその双対のボロノイ格子を組み立てる。
///
/// `points` は三角形の頂点（ボロノイ格子の母点）、`triangles` は反時計回りに並べた
/// 頂点番号。ボロノイ頂点は各三角形の外心になる。領域の境界上の母点のセルは、
/// 境界線で閉じる。
pub fn dual_of_triangulation(points: &[[f64; 2]], triangles: &[[usize; 3]]) -> Mesh {
    todo!(
        "デローニー三角形分割の双対 ({} points, {} triangles)",
        points.len(),
        triangles.len()
    )
}

/// 座標と接続から幾何量（`area_cell`, `dc_edge`, `dv_edge`, `x_edge`, `y_edge`,
/// `normal_edge`）を計算する。
pub fn compute_geometry(mesh: &mut Mesh) {
    todo!("幾何量の計算 ({} cells)", mesh.n_cells())
}

/// JIGSAW の格子ファイルを読み込む（形式は未定）。
pub fn read_jigsaw_msh(path: &Path) -> io::Result<Mesh> {
    todo!("JIGSAW 格子の読み込み ({})", path.display())
}
