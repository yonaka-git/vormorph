//! 格子生成ツール（Python）が書き出した格子ファイルを Rust 側で読み、両者の規約が
//! 一致していることを確かめる。
//!
//! `tests/data/island_coarse.mesh.toml` は次のコマンドで作ったもの（gmsh 4.15.2）。
//!
//! ```text
//! uv run --project tools python tools/make_mesh.py examples/domains/island.toml \
//!     --cell-size 4.0 -o tests/data/island_coarse.mesh.toml
//! ```
//!
//! 格子ファイルの形式を変えたら、このファイルも作り直す。

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use vormorph::mesh::file::{self, CellOrder, MeshFile};
use vormorph::mesh::voronoi::PolygonMesh;
use vormorph::mesh::{EdgeKind, NONE};

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/island_coarse.mesh.toml")
}

fn load() -> MeshFile {
    MeshFile::load(&fixture()).expect("the fixture must be a valid mesh file")
}

fn meta_f64(file: &MeshFile, key: &str) -> f64 {
    file.meta
        .get(key)
        .and_then(|v| v.as_float())
        .unwrap_or_else(|| panic!("meta.{key} must be a float"))
}

/// セル c の頂点座標（反時計回り）。
fn polygon(p: &PolygonMesh, c: usize) -> Vec<[f64; 2]> {
    p.vertices_on_cell[c]
        .iter()
        .map(|&v| [p.x_vertex[v], p.y_vertex[v]])
        .collect()
}

fn signed_area(poly: &[[f64; 2]]) -> f64 {
    let n = poly.len();
    0.5 * (0..n)
        .map(|k| {
            let [x0, y0] = poly[k];
            let [x1, y1] = poly[(k + 1) % n];
            x0 * y1 - x1 * y0
        })
        .sum::<f64>()
}

/// 向きつきの辺 (a, b) → それを反時計回りに持つセル。
fn directed_edges(p: &PolygonMesh) -> HashMap<(usize, usize), usize> {
    let mut out = HashMap::new();
    for c in 0..p.n_cells() {
        let verts = &p.vertices_on_cell[c];
        for k in 0..verts.len() {
            let edge = (verts[k], verts[(k + 1) % verts.len()]);
            assert!(out.insert(edge, c).is_none(), "edge {edge:?} appears twice");
        }
    }
    out
}

#[test]
fn header_and_counts() {
    let file = load();
    assert_eq!(file.title, "island");
    assert_eq!(file.cell_order, CellOrder::Hilbert);
    // 周囲も島もすべて壁
    assert!(file.boundary.upstream.is_empty() && file.boundary.downstream.is_empty());
    assert!(!file.boundary.wall.is_empty());
}

#[test]
fn cells_tile_the_domain() {
    let file = load();
    let p = file.to_polygons();
    let mut total = 0.0;
    for c in 0..p.n_cells() {
        let a = signed_area(&polygon(&p, c));
        assert!(a > 0.0, "cell {c} is not counterclockwise");
        total += a;
    }
    let domain_area = meta_f64(&file, "domain_area");
    assert!(
        ((total - domain_area) / domain_area).abs() < 1e-12,
        "total area {total} != domain area {domain_area}"
    );
}

#[test]
fn edges_are_shared_or_listed_as_boundary() {
    let p = load().to_polygons();
    let edges = directed_edges(&p);
    let boundary: HashMap<(usize, usize), EdgeKind> = p
        .boundary_edges
        .iter()
        .map(|e| ((e.vertices[0], e.vertices[1]), e.kind))
        .collect();
    assert_eq!(
        boundary.len(),
        p.boundary_edges.len(),
        "duplicate boundary edges"
    );
    for &(a, b) in edges.keys() {
        let shared = edges.contains_key(&(b, a));
        let listed = boundary.contains_key(&(a, b));
        assert!(
            shared != listed,
            "edge ({a}, {b}) must be either shared or a boundary edge"
        );
    }
    for &(a, b) in boundary.keys() {
        assert!(
            edges.contains_key(&(a, b)),
            "boundary edge ({a}, {b}) is not a cell edge"
        );
    }
}

#[test]
fn interior_edges_are_orthogonal_and_boundary_generators_are_inside() {
    let p = load().to_polygons();
    let edges = directed_edges(&p);
    let g = |c: usize| [p.x_cell[c], p.y_cell[c]];
    let v = |k: usize| [p.x_vertex[k], p.y_vertex[k]];
    for (&(a, b), &c1) in &edges {
        let [ax, ay] = v(a);
        let [bx, by] = v(b);
        let (tx, ty) = (bx - ax, by - ay);
        let length = tx.hypot(ty);
        match edges.get(&(b, a)) {
            Some(&c2) if length > 1e-6 => {
                let ([x1, y1], [x2, y2]) = (g(c1), g(c2));
                let (dx, dy) = (x2 - x1, y2 - y1);
                let cos = (dx * tx + dy * ty).abs() / (dx.hypot(dy) * length);
                assert!(
                    cos < 1e-10,
                    "edge ({a}, {b}) is not orthogonal: |cos| = {cos}"
                );
            }
            Some(_) => {}
            None => {
                // 境界辺：母点が辺の内側（左側）にあれば、鏡映点までの距離が正になる
                let [x, y] = g(c1);
                let d = (tx * (y - ay) - ty * (x - ax)) / length;
                assert!(
                    d > 0.0,
                    "generator of cell {c1} is not inside boundary edge ({a}, {b})"
                );
            }
        }
    }
}

#[test]
fn consecutive_cells_are_close_in_hilbert_order() {
    let file = load();
    let p = file.to_polygons();
    let h = meta_f64(&file, "cell_size");
    let mut steps: Vec<f64> = (1..p.n_cells())
        .map(|c| (p.x_cell[c] - p.x_cell[c - 1]).hypot(p.y_cell[c] - p.y_cell[c - 1]))
        .collect();
    steps.sort_by(f64::total_cmp);
    let median = steps[steps.len() / 2];
    assert!(
        median < 1.5 * h,
        "median step {median} is not local (h = {h})"
    );
}

/// Rust 側で接続配列と幾何量を組み立て、格子生成ツールの検査と同じ結果になることを確かめる。
#[test]
fn rust_builds_connectivity_and_geometry() {
    let file = load();
    let mesh = file::read(&fixture()).expect("the fixture must build into a mesh");
    assert_eq!(mesh.n_cells(), file.cells.x.len());
    assert_eq!(mesh.n_vertices(), file.vertices.x.len());

    let n_boundary = (0..mesh.n_edges())
        .filter(|&e| mesh.is_boundary_edge(e))
        .count();
    assert_eq!(n_boundary, file.n_boundary_edges());
    assert!(
        (0..mesh.n_edges())
            .all(|e| mesh.is_boundary_edge(e) == (mesh.edge_kind[e] == EdgeKind::Wall))
    );
    // 内部の辺の [0] は番号の小さいセル
    assert!(
        mesh.cells_on_edge
            .iter()
            .all(|&[c0, c1]| c1 == NONE || c0 < c1)
    );
    // オイラー標数：島が1つ（穴が1つ）の領域では V − E + C = 0
    let euler = mesh.n_vertices() as i64 - mesh.n_edges() as i64 + mesh.n_cells() as i64;
    assert_eq!(euler, 0);

    let report = mesh.check_geometry();
    let domain_area = meta_f64(&file, "domain_area");
    assert!(
        ((report.total_area - domain_area) / domain_area).abs() < 1e-12,
        "total area {} != domain area {domain_area}",
        report.total_area
    );
    assert!(
        report.max_orthogonality_error < 1e-10,
        "orthogonality error {:e}",
        report.max_orthogonality_error
    );
    assert!(
        report.max_closure_error < 1e-12,
        "closure error {:e}",
        report.max_closure_error
    );
}
