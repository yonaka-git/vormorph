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
use vormorph::hydro;
use vormorph::mesh::file::{self, CellOrder, MeshFile};
use vormorph::mesh::voronoi::{self, PolygonMesh};
use vormorph::mesh::{EdgeKind, Mesh, NONE};

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

/// 格子を平面直角座標の桁（10⁵ m）だけ平行移動して組み立てても、幾何量が座標の丸めの範囲で
/// 変わらないこと。面積は母点を原点にし、ほかの量は座標の差から計算しているため
/// （`mesh::voronoi::compute_geometry`）。
#[test]
fn geometry_does_not_depend_on_the_coordinate_origin() {
    const SHIFT: f64 = 1.0e5;
    // 平行移動した座標の丸め（10⁵ m で約 1.5e-11 m）の数十倍 [m]
    const TOL: f64 = 1e-9;
    let polygons = load().to_polygons();
    let mut moved = polygons.clone();
    for coords in [
        &mut moved.x_cell,
        &mut moved.y_cell,
        &mut moved.x_vertex,
        &mut moved.y_vertex,
    ] {
        coords.iter_mut().for_each(|x| *x += SHIFT);
    }
    let near = voronoi::build(&polygons).expect("the fixture must build");
    let far = voronoi::build(&moved).expect("the translated fixture must build");

    for c in 0..near.n_cells() {
        let (a, b) = (far.area_cell[c], near.area_cell[c]);
        assert!((a - b).abs() <= TOL * b.sqrt(), "cell {c}: area {a} != {b}");
    }
    for e in 0..near.n_edges() {
        let close = |a: f64, b: f64| (a - b).abs() <= TOL;
        assert!(close(far.dv_edge[e], near.dv_edge[e]), "edge {e}: dv_edge");
        assert!(close(far.dc_edge[e], near.dc_edge[e]), "edge {e}: dc_edge");
        assert!(
            close(far.x_edge[e] - SHIFT, near.x_edge[e]),
            "edge {e}: x_edge"
        );
        assert!(
            close(far.y_edge[e] - SHIFT, near.y_edge[e]),
            "edge {e}: y_edge"
        );
        // 法線は短い辺ほど丸めの影響を受けるので、長さを掛けた L n で比べる
        for k in 0..2 {
            let (a, b) = (
                far.dv_edge[e] * far.normal_edge[e][k],
                near.dv_edge[e] * near.normal_edge[e][k],
            );
            assert!(close(a, b), "edge {e}: normal_edge");
        }
    }
    let report = far.check_geometry();
    assert!(
        report.max_orthogonality_error < 1e-9,
        "orthogonality error {:e}",
        report.max_orthogonality_error
    );
    assert!(report.max_closure_error < 1e-12);
}

/// 位置 (x, y) \[m\] での流速 (u, v) \[m/s\]。
type Field = fn(f64, f64) -> [f64; 2];

/// 一次の流速場（勾配の成分は ±0.2 1/s または 0）。
const LINEAR_FIELDS: [(&str, Field); 4] = [
    ("u = 0.2x", |x, _| [0.2 * x, 0.0]),
    ("u = 0.2y", |_, y| [0.2 * y, 0.0]),
    ("rotation", |x, y| [-0.2 * y, 0.2 * x]),
    ("acceleration", |x, y| [0.2 * x, -0.2 * y]),
];

/// 一次の流速場 `field` の法線成分を辺の中点で与え、Perot 法で再構成したセルの流速と、母点での
/// 厳密な値との差を、1セル分の流速の変化 0.2h で割った値（セルごと）。
fn perot_errors(mesh: &Mesh, field: Field, h: f64) -> Vec<f64> {
    let u_n: Vec<f64> = (0..mesh.n_edges())
        .map(|e| {
            let [u, v] = field(mesh.x_edge[e], mesh.y_edge[e]);
            let [nx, ny] = mesh.normal_edge[e];
            u * nx + v * ny
        })
        .collect();
    let mut u_cell = vec![[0.0; 2]; mesh.n_cells()];
    hydro::reconstruct_cell_velocity(mesh, &u_n, &mut u_cell);
    (0..mesh.n_cells())
        .map(|c| {
            let [u, v] = field(mesh.x_cell[c], mesh.y_cell[c]);
            (u_cell[c][0] - u).hypot(u_cell[c][1] - v) / (0.2 * h)
        })
        .collect()
}

/// 格子生成ツールの格子で、壁ぎわのセルの Perot 法の精度が内部のセルと大きく変わらないこと。
///
/// 境界の母点は境界から約 0.43h の位置にあり、壁はその母点と鏡映点の二等分線なので、壁ぎわの
/// セルでも母点は重心の近くにある（ずれは格子幅 0.5 m の同じ領域で最大 0.06h）。例外は外周の
/// 凸な角に接するセルで、ずれが約 0.2h あり、誤差が大きいので除く
/// （`docs/design/mesh-generation.md` の第8節）。
#[test]
fn perot_reconstruction_is_accurate_at_walls() {
    let file = load();
    let mesh = file::read(&fixture()).expect("the fixture must build into a mesh");
    let h = meta_f64(&file, "cell_size");
    // examples/domains/island.toml の外周の角
    let corners = [[0.0, 0.0], [60.0, 0.0], [60.0, 30.0], [0.0, 30.0]];
    let at_corner = |c: usize| {
        mesh.vertices_on_cell[c].iter().any(|&v| {
            corners
                .iter()
                .any(|p| (mesh.x_vertex[v] - p[0]).hypot(mesh.y_vertex[v] - p[1]) < 1e-9)
        })
    };
    let at_wall = |c: usize| {
        mesh.edges_on_cell[c]
            .iter()
            .any(|&e| mesh.is_boundary_edge(e))
    };
    for (name, field) in LINEAR_FIELDS {
        let errors = perot_errors(&mesh, field, h);
        let max_of = |keep: &dyn Fn(usize) -> bool| {
            (0..mesh.n_cells())
                .filter(|&c| keep(c))
                .map(|c| errors[c])
                .fold(0.0, f64::max)
        };
        let interior = max_of(&|c| !at_wall(c));
        let wall = max_of(&|c| at_wall(c) && !at_corner(c));
        let corner = max_of(&|c| at_corner(c));
        println!(
            "{name}: max error / (0.2 h) {interior:.3} (interior), {wall:.3} (walls), {corner:.3} (corners)"
        );
        assert!(interior < 0.1, "{name}: interior error {interior}");
        assert!(wall < 0.15, "{name}: wall error {wall}");
    }
}
