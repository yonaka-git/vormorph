//! 多角形のセルの並びから、MPAS 形式の接続配列と幾何量を組み立てる。
//!
//! 格子の作り方によらず（直線水路の [`channel`](super::channel)、格子生成ツールの出力を読む
//! [`file`](super::file)）、[`PolygonMesh`] を作ってから [`build`] を通す。辺の向きと符号の
//! 規約を決めるのはこのモジュールと [`Mesh::derive_cell_connectivity`] だけにする。
//!
//! ## 直交性
//!
//! ボロノイ図では、隣り合う母点を結ぶ線分と、両者を隔てるボロノイ辺が直交する。C 格子の
//! 水位勾配 (η_2 − η_1)/δ_j や、河床の法線方向勾配が2点差分で厳密になるのはこの性質による。
//! 境界辺は、母点をその辺の直線に関して鏡映した点との二等分線になるように格子を作る
//! （`tools/vormesh` の鏡映対の母点、`docs/design/mesh-generation.md`）。
//!
//! ## 番号のつけ方
//!
//! - セルと頂点の番号は [`PolygonMesh`] のまま使う。格子生成ツールはセルをヒルベルト曲線の順に
//!   並べ、頂点をセルの順に初めて現れた順に並べて書き出す。
//! - 辺の番号は、セルを番号順にたどり、各セルの辺を反時計回りにたどって初めて現れた順に
//!   つける。辺の `cells_on_edge[e][0]` は先に現れたセル（番号の小さいセル）になり、
//!   `vertices_on_edge[e]` はそのセルがたどる向きになる。
//!   これにより、セル単位のループで参照する辺と頂点も、番号の近い範囲に集まる。
//!
//! ## 幾何量
//!
//! | 量 | 定義 |
//! |---|---|
//! | `area_cell` | 頂点を結ぶ多角形の面積（母点を原点にした shoelace 公式） |
//! | `dv_edge` | 辺の両端の頂点の距離 L |
//! | `x_edge`, `y_edge` | 辺の両端の頂点の中点 |
//! | `normal_edge` | 辺の接線 t = (x_b − x_a)/L を時計回りに 90° 回したもの (t_y, −t_x) |
//! | `dc_edge` | 内部の辺は母点間の距離。境界辺は母点から辺の直線までの距離の2倍 |
//!
//! ## 参考文献
//!
//! - MPAS Mesh Specification v1.0 (2015): 接続配列の定義
//! - Du, Faber & Gunzburger (1999): 重心ボロノイ分割

use super::{EdgeKind, Jagged, Mesh, MeshError, NONE};
use std::collections::HashMap;

/// 境界辺の両端の頂点と種類。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoundaryEdge {
    /// 両端の頂点番号。セルから見て反時計回りの向き。
    pub vertices: [usize; 2],
    /// 辺の種類（[`EdgeKind::Interior`] 以外）。
    pub kind: EdgeKind,
}

/// 多角形のセルの並びで表したボロノイ格子（接続配列を組み立てる前の形）。
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PolygonMesh {
    /// 母点の x 座標 \[m\]。
    pub x_cell: Vec<f64>,
    /// 母点の y 座標 \[m\]。
    pub y_cell: Vec<f64>,
    /// ボロノイ頂点の x 座標 \[m\]。
    pub x_vertex: Vec<f64>,
    /// ボロノイ頂点の y 座標 \[m\]。
    pub y_vertex: Vec<f64>,
    /// セルの頂点番号（反時計回り）。セル c の辺 k は `vertices_on_cell[c][k]` →
    /// `vertices_on_cell[c][(k + 1) % n]`。
    pub vertices_on_cell: Jagged<usize>,
    /// 境界辺。相手のセルがない辺は、すべてここに種類つきで挙げる。
    pub boundary_edges: Vec<BoundaryEdge>,
}

impl PolygonMesh {
    /// セル数。
    pub fn n_cells(&self) -> usize {
        self.x_cell.len()
    }

    /// 頂点の数。
    pub fn n_vertices(&self) -> usize {
        self.x_vertex.len()
    }
}

/// 多角形のセルの並びから、MPAS 形式の接続配列と幾何量を持つ [`Mesh`] を組み立てる。
///
/// 辺の番号は[モジュールの説明](self)のとおりにつけ、[`Mesh::derive_cell_connectivity`] と
/// [`compute_geometry`] を通す。次の場合は、見つかった問題をすべて挙げて失敗する。
///
/// - 配列の長さが合わない、頂点番号が範囲外、頂点が3つ未満、同じ頂点を2度含むセルがある
/// - 内部の辺を2つのセルが逆向きに共有していない（同じ向きに共有する、3つ以上のセルが共有する）
/// - 相手のいない辺が `boundary_edges` にない（T 字の接合など）、または `boundary_edges` の辺が
///   セルの辺でない・内部の辺である・向きが逆である・重複している
/// - 面積が正でない（時計回りの）セル、長さ 0 の辺、母点が一致する内部の辺、母点が境界辺の
///   外側にある境界辺がある
pub fn build(polygons: &PolygonMesh) -> Result<Mesh, MeshError> {
    MeshError::check(check_input(polygons))?;
    let edges = number_edges(polygons);
    let (edge_kind, boundary_problems) = classify_edges(polygons, &edges);
    let mut problems = edges.problems;
    problems.extend(boundary_problems);
    MeshError::check(problems)?;

    let mut mesh = Mesh {
        x_cell: polygons.x_cell.clone(),
        y_cell: polygons.y_cell.clone(),
        edges_on_cell: edges.edges_on_cell,
        vertices_on_cell: polygons.vertices_on_cell.clone(),
        cells_on_edge: edges.cells_on_edge,
        vertices_on_edge: edges.vertices_on_edge,
        edge_kind,
        x_vertex: polygons.x_vertex.clone(),
        y_vertex: polygons.y_vertex.clone(),
        ..Mesh::default()
    };
    mesh.derive_cell_connectivity();
    compute_geometry(&mut mesh);
    MeshError::check(check_built(&mesh))?;
    Ok(mesh)
}

/// 座標と接続から幾何量（`area_cell`, `dc_edge`, `dv_edge`, `x_edge`, `y_edge`,
/// `normal_edge`）を計算する。定義は[モジュールの説明](self)の表のとおり。
///
/// `x_cell`, `y_cell`, `x_vertex`, `y_vertex`, `vertices_on_cell`, `cells_on_edge`,
/// `vertices_on_edge` が規約どおりに入っていることを前提にする。長さ 0 の辺の法線は NaN になる
/// （[`build`] はそのような格子を拒否する）。
pub fn compute_geometry(mesh: &mut Mesh) {
    let area_cell = (0..mesh.n_cells())
        .map(|c| {
            let (xc, yc) = (mesh.x_cell[c], mesh.y_cell[c]);
            let verts = &mesh.vertices_on_cell[c];
            let n = verts.len();
            0.5 * (0..n)
                .map(|k| {
                    let (a, b) = (verts[k], verts[(k + 1) % n]);
                    let (xa, ya) = (mesh.x_vertex[a] - xc, mesh.y_vertex[a] - yc);
                    let (xb, yb) = (mesh.x_vertex[b] - xc, mesh.y_vertex[b] - yc);
                    xa * yb - xb * ya
                })
                .sum::<f64>()
        })
        .collect();

    let n_edges = mesh.n_edges();
    let mut dv_edge = Vec::with_capacity(n_edges);
    let mut dc_edge = Vec::with_capacity(n_edges);
    let mut x_edge = Vec::with_capacity(n_edges);
    let mut y_edge = Vec::with_capacity(n_edges);
    let mut normal_edge = Vec::with_capacity(n_edges);
    for e in 0..n_edges {
        let [a, b] = mesh.vertices_on_edge[e];
        let (xa, ya) = (mesh.x_vertex[a], mesh.y_vertex[a]);
        let (xb, yb) = (mesh.x_vertex[b], mesh.y_vertex[b]);
        let length = (xb - xa).hypot(yb - ya);
        let normal = [(yb - ya) / length, -(xb - xa) / length];
        let [c0, c1] = mesh.cells_on_edge[e];
        let (x0, y0) = (mesh.x_cell[c0], mesh.y_cell[c0]);
        let dc = if c1 == NONE {
            // 母点から辺の直線までの距離（外向きの法線の向きに測る）の2倍 = 鏡映点までの距離
            2.0 * ((xa - x0) * normal[0] + (ya - y0) * normal[1])
        } else {
            (mesh.x_cell[c1] - x0).hypot(mesh.y_cell[c1] - y0)
        };
        dv_edge.push(length);
        dc_edge.push(dc);
        x_edge.push(0.5 * (xa + xb));
        y_edge.push(0.5 * (ya + yb));
        normal_edge.push(normal);
    }

    mesh.area_cell = area_cell;
    mesh.dv_edge = dv_edge;
    mesh.dc_edge = dc_edge;
    mesh.x_edge = x_edge;
    mesh.y_edge = y_edge;
    mesh.normal_edge = normal_edge;
}

/// 配列の長さ、番号の範囲、セルの頂点の数と重複、境界辺の種類を検査する。
///
/// ここで問題が見つかったら、接続の組み立てに進まない（範囲外の番号で添字を引かないため）。
fn check_input(p: &PolygonMesh) -> Vec<String> {
    let mut problems = Vec::new();
    let (n_cells, n_vertices) = (p.n_cells(), p.n_vertices());
    if n_cells == 0 {
        problems.push("the mesh has no cells".to_string());
    }
    if p.y_cell.len() != n_cells || p.vertices_on_cell.len() != n_cells {
        problems.push(format!(
            "x_cell, y_cell and vertices_on_cell must have the same length ({}, {}, {})",
            n_cells,
            p.y_cell.len(),
            p.vertices_on_cell.len()
        ));
    }
    if p.y_vertex.len() != n_vertices {
        problems.push(format!(
            "x_vertex and y_vertex must have the same length ({}, {})",
            n_vertices,
            p.y_vertex.len()
        ));
    }
    let coords = p.x_cell.iter().chain(&p.y_cell);
    if !coords
        .chain(&p.x_vertex)
        .chain(&p.y_vertex)
        .all(|v| v.is_finite())
    {
        problems.push("coordinates must be finite".to_string());
    }
    for (c, verts) in p.vertices_on_cell.rows().enumerate() {
        if verts.len() < 3 {
            problems.push(format!(
                "cell {c} has {} vertices (at least 3 are needed)",
                verts.len()
            ));
        }
        if let Some(&v) = verts.iter().find(|&&v| v >= n_vertices) {
            problems.push(format!(
                "cell {c} refers to vertex {v}, which is out of range (n_vertices = {n_vertices})"
            ));
        }
        let mut sorted = verts.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        if sorted.len() != verts.len() {
            problems.push(format!("cell {c} lists a vertex twice"));
        }
    }
    for (k, edge) in p.boundary_edges.iter().enumerate() {
        let [a, b] = edge.vertices;
        if a >= n_vertices || b >= n_vertices || a == b {
            problems.push(format!(
                "boundary edge {k} ({a}, {b}) is not a pair of distinct vertices in range"
            ));
        }
        if edge.kind == EdgeKind::Interior {
            problems.push(format!("boundary edge {k} ({a}, {b}) has kind Interior"));
        }
    }
    problems
}

/// [`number_edges`] の結果。
struct NumberedEdges {
    /// 向きつきの辺 (a, b) → 辺番号。辺を初めてたどったセルの向きだけを登録する。
    edge_of: HashMap<(usize, usize), usize>,
    edges_on_cell: Jagged<usize>,
    cells_on_edge: Vec<[usize; 2]>,
    vertices_on_edge: Vec<[usize; 2]>,
    problems: Vec<String>,
}

/// セルを番号順に、各セルの辺を反時計回りにたどり、初めて現れた順に辺の番号をつける。
///
/// 逆向きの辺がすでにあれば、2つのセルが共有する内部の辺とする。
fn number_edges(p: &PolygonMesh) -> NumberedEdges {
    let n_directed = p.vertices_on_cell.total_len();
    let mut out = NumberedEdges {
        edge_of: HashMap::with_capacity(n_directed),
        edges_on_cell: Jagged::new(),
        cells_on_edge: Vec::with_capacity(n_directed / 2 + 1),
        vertices_on_edge: Vec::with_capacity(n_directed / 2 + 1),
        problems: Vec::new(),
    };
    for (c, verts) in p.vertices_on_cell.rows().enumerate() {
        let n = verts.len();
        let mut row = Vec::with_capacity(n);
        for k in 0..n {
            let (a, b) = (verts[k], verts[(k + 1) % n]);
            let e = if let Some(&e) = out.edge_of.get(&(b, a)) {
                let [c0, c1] = out.cells_on_edge[e];
                if c1 == NONE {
                    out.cells_on_edge[e][1] = c;
                } else {
                    out.problems.push(format!(
                        "edge ({a}, {b}) is shared by more than two cells ({c0}, {c1} and {c})"
                    ));
                }
                e
            } else if let Some(&e) = out.edge_of.get(&(a, b)) {
                out.problems.push(format!(
                    "cells {} and {c} traverse edge ({a}, {b}) in the same direction \
                     (overlapping cells or inconsistent orientation)",
                    out.cells_on_edge[e][0]
                ));
                e
            } else {
                let e = out.cells_on_edge.len();
                out.edge_of.insert((a, b), e);
                out.cells_on_edge.push([c, NONE]);
                out.vertices_on_edge.push([a, b]);
                e
            };
            row.push(e);
        }
        out.edges_on_cell.push_row(row);
    }
    out
}

/// `boundary_edges` を相手のいない辺と照合し、辺の種類と見つかった問題を返す。
fn classify_edges(p: &PolygonMesh, edges: &NumberedEdges) -> (Vec<EdgeKind>, Vec<String>) {
    let (edge_of, cells_on_edge) = (&edges.edge_of, &edges.cells_on_edge);
    let mut problems = Vec::new();
    let mut edge_kind = vec![EdgeKind::Interior; cells_on_edge.len()];
    let mut listed = vec![false; cells_on_edge.len()];
    for edge in &p.boundary_edges {
        let [a, b] = edge.vertices;
        match (edge_of.get(&(a, b)), edge_of.get(&(b, a))) {
            (Some(&e), _) | (_, Some(&e)) if cells_on_edge[e][1] != NONE => {
                let [c0, c1] = cells_on_edge[e];
                problems.push(format!(
                    "boundary edge ({a}, {b}) is an interior edge shared by cells {c0} and {c1}"
                ));
            }
            (Some(&e), _) => {
                if listed[e] {
                    problems.push(format!("boundary edge ({a}, {b}) is listed twice"));
                }
                listed[e] = true;
                edge_kind[e] = edge.kind;
            }
            (None, Some(_)) => problems.push(format!(
                "boundary edge ({a}, {b}) runs clockwise around its cell; expected ({b}, {a})"
            )),
            (None, None) => {
                problems.push(format!(
                    "boundary edge ({a}, {b}) is not an edge of any cell"
                ));
            }
        }
    }
    for (e, &[c0, c1]) in cells_on_edge.iter().enumerate() {
        if c1 == NONE && !listed[e] {
            let [a, b] = edges.vertices_on_edge[e];
            problems.push(format!(
                "edge ({a}, {b}) of cell {c0} has no neighbouring cell and is not listed as \
                 a boundary edge (a T-junction or a gap between cells?)"
            ));
        }
    }
    (edge_kind, problems)
}

/// 幾何量を計算したあとで、面積・辺の長さ・母点間の距離が正であることを検査する。
fn check_built(mesh: &Mesh) -> Vec<String> {
    let mut problems = Vec::new();
    for (c, &area) in mesh.area_cell.iter().enumerate() {
        if !positive(area) {
            problems.push(format!(
                "cell {c} is not counterclockwise (signed area {area:e})"
            ));
        }
    }
    for e in 0..mesh.n_edges() {
        let [a, b] = mesh.vertices_on_edge[e];
        let [c0, c1] = mesh.cells_on_edge[e];
        if !positive(mesh.dv_edge[e]) {
            problems.push(format!("edge ({a}, {b}) has zero length"));
        } else if !positive(mesh.dc_edge[e]) {
            problems.push(if c1 == NONE {
                format!("the generator of cell {c0} is not inside its boundary edge ({a}, {b})")
            } else {
                format!("cells {c0} and {c1} have the same generator")
            });
        }
    }
    problems
}

/// `x > 0`（NaN は偽）。
fn positive(x: f64) -> bool {
    x > 0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 2つの正方形のセル（x = 0..1 と 1..2、y = 0..1）。左端が上流、右端が下流。
    ///
    /// セル 0 の母点は中心から x 方向にずらしてあり、`dc_edge` の内部と境界での定義の違いが
    /// 値に表れる。
    ///
    /// ```text
    /// 3 ─── 2 ─── 5
    /// │  0  │  1  │
    /// 0 ─── 1 ─── 4
    /// ```
    fn two_squares() -> PolygonMesh {
        let wall = |a, b| BoundaryEdge {
            vertices: [a, b],
            kind: EdgeKind::Wall,
        };
        PolygonMesh {
            x_cell: vec![0.4, 1.5],
            y_cell: vec![0.5, 0.5],
            x_vertex: vec![0.0, 1.0, 1.0, 0.0, 2.0, 2.0],
            y_vertex: vec![0.0, 0.0, 1.0, 1.0, 0.0, 1.0],
            vertices_on_cell: Jagged::from_rows([vec![0, 1, 2, 3], vec![1, 4, 5, 2]]),
            boundary_edges: vec![
                wall(0, 1),
                wall(2, 3),
                wall(1, 4),
                wall(5, 2),
                BoundaryEdge {
                    vertices: [3, 0],
                    kind: EdgeKind::Upstream,
                },
                BoundaryEdge {
                    vertices: [4, 5],
                    kind: EdgeKind::Downstream,
                },
            ],
        }
    }

    fn rejected(p: &PolygonMesh) -> Vec<String> {
        match build(p) {
            Ok(_) => panic!("the mesh must be rejected"),
            Err(err) => err.problems().to_vec(),
        }
    }

    fn assert_reported(problems: &[String], needle: &str) {
        assert!(
            problems.iter().any(|p| p.contains(needle)),
            "expected a problem containing {needle:?}, got {problems:?}"
        );
    }

    #[test]
    fn two_squares_connectivity() {
        let mesh = build(&two_squares()).expect("valid mesh");
        assert_eq!(mesh.n_cells(), 2);
        assert_eq!(mesh.n_edges(), 7);
        // セル 0 の辺が 0..4、セル 1 の新しい辺が 4..7。共有辺はセル 0 の向きで登録される
        assert_eq!(&mesh.edges_on_cell[0], &[0, 1, 2, 3]);
        assert_eq!(&mesh.edges_on_cell[1], &[4, 5, 6, 1]);
        assert_eq!(mesh.cells_on_edge[1], [0, 1]);
        assert_eq!(mesh.vertices_on_edge[1], [1, 2]);
        assert_eq!(mesh.cells_on_edge[5], [1, NONE]);
        assert_eq!(&mesh.edge_sign_on_cell[1], &[1.0, 1.0, 1.0, -1.0]);
        assert_eq!(&mesh.cells_on_cell[0], &[NONE, 1, NONE, NONE]);
        let kinds = [
            EdgeKind::Wall,
            EdgeKind::Interior,
            EdgeKind::Wall,
            EdgeKind::Upstream,
            EdgeKind::Wall,
            EdgeKind::Downstream,
            EdgeKind::Wall,
        ];
        assert_eq!(mesh.edge_kind, kinds);
    }

    #[test]
    fn two_squares_geometry() {
        let mesh = build(&two_squares()).expect("valid mesh");
        let close = |a: f64, b: f64| (a - b).abs() < 1e-14;
        assert!(mesh.area_cell.iter().all(|&a| close(a, 1.0)));
        assert!(mesh.dv_edge.iter().all(|&l| close(l, 1.0)));
        // 共有辺：法線はセル 0 → 1、母点間の距離
        assert_eq!(mesh.normal_edge[1], [1.0, 0.0]);
        assert!(close(mesh.dc_edge[1], 1.1));
        assert!(close(mesh.x_edge[1], 1.0) && close(mesh.y_edge[1], 0.5));
        // 上流端：法線は −x（領域の外）、母点から辺までの距離 0.4 の2倍
        assert_eq!(mesh.normal_edge[3], [-1.0, 0.0]);
        assert!(close(mesh.dc_edge[3], 0.8));
        // 下流端：母点から辺までの距離 0.5 の2倍
        assert_eq!(mesh.normal_edge[5], [1.0, 0.0]);
        assert!(close(mesh.dc_edge[5], 1.0));
        // 底面の壁：法線は −y
        assert_eq!(mesh.normal_edge[0], [0.0, -1.0]);

        let report = mesh.check_geometry();
        assert!(close(report.total_area, 2.0));
        assert_eq!(report.max_orthogonality_error, 0.0);
        assert!(report.max_closure_error < 1e-15);
    }

    #[test]
    fn missing_boundary_edge_is_reported() {
        let mut p = two_squares();
        p.boundary_edges.retain(|e| e.kind != EdgeKind::Upstream);
        assert_reported(
            &rejected(&p),
            "edge (3, 0) of cell 0 has no neighbouring cell",
        );
    }

    #[test]
    fn inconsistent_boundary_edges_are_reported() {
        let mut p = two_squares();
        let upstream = p.boundary_edges[4];
        p.boundary_edges.push(upstream);
        p.boundary_edges[0].vertices = [1, 0];
        p.boundary_edges.push(BoundaryEdge {
            vertices: [1, 2],
            kind: EdgeKind::Wall,
        });
        p.boundary_edges.push(BoundaryEdge {
            vertices: [0, 2],
            kind: EdgeKind::Wall,
        });
        let problems = rejected(&p);
        assert_reported(&problems, "(3, 0) is listed twice");
        assert_reported(&problems, "(1, 0) runs clockwise");
        assert_reported(&problems, "(1, 2) is an interior edge");
        assert_reported(&problems, "(0, 2) is not an edge of any cell");
    }

    #[test]
    fn t_junction_is_reported() {
        // 左の縦長のセル 0 の右辺 (1, 2) に、右の2つのセルの辺 (6, 1), (2, 6) が T 字に接する
        let p = PolygonMesh {
            x_cell: vec![0.5, 1.5, 1.5],
            y_cell: vec![1.0, 0.5, 1.5],
            x_vertex: vec![0.0, 1.0, 1.0, 0.0, 2.0, 2.0, 1.0, 2.0],
            y_vertex: vec![0.0, 0.0, 2.0, 2.0, 0.0, 1.0, 1.0, 2.0],
            vertices_on_cell: Jagged::from_rows([
                vec![0, 1, 2, 3],
                vec![1, 4, 5, 6],
                vec![6, 5, 7, 2],
            ]),
            boundary_edges: [[0, 1], [2, 3], [3, 0], [1, 4], [4, 5], [5, 7], [7, 2]]
                .into_iter()
                .map(|vertices| BoundaryEdge {
                    vertices,
                    kind: EdgeKind::Wall,
                })
                .collect(),
        };
        let problems = rejected(&p);
        assert_reported(&problems, "edge (1, 2) of cell 0 has no");
        assert_reported(&problems, "edge (6, 1) of cell 1 has no");
        assert_reported(&problems, "edge (2, 6) of cell 2 has no");
    }

    #[test]
    fn overlapping_cells_are_reported() {
        let mut p = two_squares();
        p.vertices_on_cell = Jagged::from_rows([vec![0, 1, 2, 3], vec![0, 1, 2, 3]]);
        assert_reported(&rejected(&p), "in the same direction");
    }

    #[test]
    fn clockwise_cell_is_reported() {
        let p = PolygonMesh {
            x_cell: vec![0.5],
            y_cell: vec![0.5],
            x_vertex: vec![0.0, 1.0, 1.0, 0.0],
            y_vertex: vec![0.0, 0.0, 1.0, 1.0],
            vertices_on_cell: Jagged::from_rows([vec![3, 2, 1, 0]]),
            boundary_edges: [[3, 2], [2, 1], [1, 0], [0, 3]]
                .into_iter()
                .map(|vertices| BoundaryEdge {
                    vertices,
                    kind: EdgeKind::Wall,
                })
                .collect(),
        };
        let problems = rejected(&p);
        assert_reported(&problems, "cell 0 is not counterclockwise");
        // 時計回りのセルでは、内側にある母点が境界辺の外側にあるように見える
        assert_reported(&problems, "not inside its boundary edge");
    }

    #[test]
    fn generator_outside_boundary_edge_is_reported() {
        let mut p = two_squares();
        p.x_cell[0] = -0.1;
        assert_reported(&rejected(&p), "generator of cell 0 is not inside");
    }

    #[test]
    fn invalid_input_is_reported_without_panicking() {
        let mut p = two_squares();
        p.vertices_on_cell = Jagged::from_rows([vec![0, 1, 2, 3], vec![1, 4, 6, 2]]);
        p.y_cell.pop();
        p.boundary_edges[0].kind = EdgeKind::Interior;
        let problems = rejected(&p);
        assert_reported(&problems, "must have the same length");
        assert_reported(&problems, "cell 1 refers to vertex 6");
        assert_reported(&problems, "has kind Interior");
        assert_eq!(
            rejected(&PolygonMesh::default()),
            vec!["the mesh has no cells".to_string()]
        );
    }
}
