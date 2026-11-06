//! 直線矩形水路の格子生成（検証ケースとデモ用）。
//!
//! x 軸を流下方向（上流端 x = 0、下流端 x = L）、y 軸を横断方向（側壁 y = 0 と y = W）にとる。
//!
//! - [`ChannelLayout::Rectangular`]：母点を正方格子に置く。セルは長方形（ほぼ正方形）。
//! - [`ChannelLayout::Hexagonal`]：母点を三角格子に置く。内部のセルは六角形で、側壁と
//!   上下流端に接するセルは、ボロノイ図を水路の矩形で切り取った多角形になる。
//!
//! 境界辺の種類は、x = 0 が [`EdgeKind::Upstream`]、x = L が [`EdgeKind::Downstream`]、
//! y = 0 と y = W が [`EdgeKind::Wall`]。
//!
//! ## 母点の配置
//!
//! 縦断方向の間隔 Δx = L / round(L/h)、横断方向の間隔 Δy = W / n_y とする（h は `cell_size`）。
//! 正方形配置では n_y = round(W/h)、六角形配置では n_y = round(W / (√3 h/2)) とし、
//! 水路の寸法に合わせて間隔をわずかに調整する（六角形は正六角形からわずかにつぶれる）。
//!
//! - 行 j（j = 0, 1, …, n_y − 1）の母点は y = (j + ½)Δy に置く。最も外側の行が側壁から Δy/2 の
//!   位置に来るので、側壁はその行の母点と鏡映点の二等分線に一致する
//!   （`docs/design/mesh-generation.md` の鏡映対の母点と同じ性質）。
//! - 正方形配置と、六角形配置の偶数行では、母点を x = (i + ½)Δx（i = 0, …, n_x − 1）に置く。
//! - 六角形配置の奇数行は半間隔ずらして x = iΔx（i = 1, …, n_x − 1）に置く。上下流端に母点が
//!   乗らないように、行の端の母点を1つ減らしている。上下流端に接する偶数行のセルは、その分
//!   端まで広がった多角形になる。
//!
//! ## セルの作り方
//!
//! 各セルは、水路の矩形を、近くの母点との垂直二等分線で順に切り取って求める（凸多角形と
//! 半平面の交差）。母点 g のセルの頂点のうち g から最も遠いものまでの距離を r とすると、
//! g から 2r より遠い母点の二等分線はセルを切らない。探索した範囲がこの条件を満たすことを、
//! セルごとに確かめる。
//!
//! 隣り合うセルで別々に計算した頂点は、母点の間隔の [`TOLERANCE`] 倍以内なら同じ頂点とみなす。
//! 4つ以上の母点が同一円周上にある縮退した配置（正方形配置のすべての頂点や、六角形配置で
//! Δy = (√3/2)Δx のときの上下流端）でも、長さが丸め誤差程度の辺は残らない。
//!
//! セルは行ごと（y の小さい行から、各行を x の順）に並べる。小さな検証用の格子なので、
//! ヒルベルト曲線の順への並べ替えはしない。生成した多角形の並び（[`PolygonMesh`]）は
//! [`voronoi::build`] を通して、符号と幾何量をほかの格子と同じ手順で計算する。

use super::voronoi::{self, BoundaryEdge, PolygonMesh};
use super::{EdgeKind, Jagged, Mesh, MeshError};
use crate::config::{ChannelLayout, MeshConfig};
use std::collections::HashMap;

/// 頂点の同一視と、頂点が二等分線の上にあるかの判定に使う許容誤差
/// （母点の間隔 min(Δx, Δy) に対する比）。
pub const TOLERANCE: f64 = 1.0e-9;

/// 直線水路の諸元。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelSpec {
    /// 水路延長 L \[m\]。
    pub length: f64,
    /// 水路幅 W \[m\]。
    pub width: f64,
    /// 代表格子幅（母点の間隔）\[m\]。
    pub cell_size: f64,
    /// 母点の配置。
    pub layout: ChannelLayout,
}

impl ChannelSpec {
    /// 設定の `[mesh]` から作る。`kind = "channel"` 以外なら `None`。
    pub fn from_config(mesh: &MeshConfig) -> Option<Self> {
        match *mesh {
            MeshConfig::Channel {
                length,
                width,
                cell_size,
                layout,
            } => Some(Self {
                length,
                width,
                cell_size,
                layout,
            }),
            MeshConfig::File { .. } => None,
        }
    }

    /// 寸法が有限の正の値であることを検査する。
    fn problems(&self) -> Vec<String> {
        [
            ("length", self.length),
            ("width", self.width),
            ("cell_size", self.cell_size),
        ]
        .into_iter()
        .filter(|&(_, v)| !(v.is_finite() && v > 0.0))
        .map(|(name, v)| format!("channel {name} must be finite and > 0, got {v}"))
        .collect()
    }
}

/// 直線水路の格子を生成する。
///
/// 横断方向のセル数は W / `cell_size`（六角形配置では W / (√3/2 × `cell_size`)）を丸めた
/// 整数とし、母点の間隔はそれに合わせてわずかに調整する（縦断方向も同様）。
/// 寸法が有限の正の値でなければ失敗する。
pub fn build(spec: &ChannelSpec) -> Result<Mesh, MeshError> {
    MeshError::check(spec.problems())?;
    let lattice = Lattice::new(spec);
    let polygons = lattice.polygons()?;
    voronoi::build(&polygons)
}

/// 母点の配置。
struct Lattice {
    length: f64,
    width: f64,
    /// 縦断方向の母点の間隔 Δx \[m\]。
    dx: f64,
    /// 行の間隔 Δy \[m\]。
    dy: f64,
    /// 母点の座標（行ごと、各行を x の順）。
    generators: Vec<[f64; 2]>,
}

impl Lattice {
    fn new(spec: &ChannelSpec) -> Self {
        let count = |extent: f64, spacing: f64| (extent / spacing).round().max(1.0) as usize;
        let nx = count(spec.length, spec.cell_size);
        let row_spacing = match spec.layout {
            ChannelLayout::Rectangular => spec.cell_size,
            ChannelLayout::Hexagonal => 0.5 * 3f64.sqrt() * spec.cell_size,
        };
        let ny = count(spec.width, row_spacing);
        let (dx, dy) = (spec.length / nx as f64, spec.width / ny as f64);

        let mut generators = Vec::new();
        for j in 0..ny {
            let y = (j as f64 + 0.5) * dy;
            let staggered = spec.layout == ChannelLayout::Hexagonal && j % 2 == 1;
            if staggered {
                generators.extend((1..nx).map(|i| [i as f64 * dx, y]));
            } else {
                generators.extend((0..nx).map(|i| [(i as f64 + 0.5) * dx, y]));
            }
        }
        Self {
            length: spec.length,
            width: spec.width,
            dx,
            dy,
            generators,
        }
    }

    /// 各母点のボロノイセルを水路の矩形で切り取り、多角形の並びにする。
    fn polygons(&self) -> Result<PolygonMesh, MeshError> {
        let tol = TOLERANCE * self.dx.min(self.dy);
        let radius = 3.0 * self.dx.max(self.dy);
        let buckets = Buckets::new(&self.generators, radius, self.length, self.width);
        let rectangle = [
            [0.0, 0.0],
            [self.length, 0.0],
            [self.length, self.width],
            [0.0, self.width],
        ];

        let mut vertices = VertexIndex::new(tol);
        let mut vertices_on_cell = Jagged::new();
        let mut problems = Vec::new();
        for (c, &g) in self.generators.iter().enumerate() {
            let mut candidates: Vec<(f64, usize)> = buckets
                .near(g)
                .filter(|&q| q != c)
                .map(|q| (distance(g, self.generators[q]), q))
                .collect();
            candidates.sort_by(|a, b| a.0.total_cmp(&b.0));

            let mut polygon = rectangle.to_vec();
            for &(d, q) in &candidates {
                if d > 2.0 * max_distance(g, &polygon) {
                    break;
                }
                polygon = clip(&polygon, g, self.generators[q], tol);
            }
            if 2.0 * max_distance(g, &polygon) > radius {
                problems.push(format!(
                    "internal error: the neighbour search around cell {c} is too narrow"
                ));
            }

            let mut row: Vec<usize> = polygon.iter().map(|&p| vertices.insert(p)).collect();
            row.dedup();
            while row.len() > 1 && row.first() == row.last() {
                row.pop();
            }
            if row.len() < 3 {
                problems.push(format!("cell {c} degenerated to {} vertices", row.len()));
            }
            vertices_on_cell.push_row(row);
        }
        MeshError::check(problems)?;

        // 境界の近くで同一視した頂点を、矩形の辺の上に正確に乗せる
        for x in &mut vertices.x {
            snap(x, 0.0, tol);
            snap(x, self.length, tol);
        }
        for y in &mut vertices.y {
            snap(y, 0.0, tol);
            snap(y, self.width, tol);
        }

        let boundary_edges = self.boundary_edges(&vertices, &vertices_on_cell);
        Ok(PolygonMesh {
            x_cell: self.generators.iter().map(|g| g[0]).collect(),
            y_cell: self.generators.iter().map(|g| g[1]).collect(),
            x_vertex: vertices.x,
            y_vertex: vertices.y,
            vertices_on_cell,
            boundary_edges,
        })
    }

    /// 両端の頂点が矩形の同じ辺の上にあるセルの辺を、境界辺として種類を決める。
    fn boundary_edges(
        &self,
        vertices: &VertexIndex,
        vertices_on_cell: &Jagged<usize>,
    ) -> Vec<BoundaryEdge> {
        // 矩形の辺（x = 0, x = L, y = 0, y = W）とその種類
        const KINDS: [EdgeKind; 4] = [
            EdgeKind::Upstream,
            EdgeKind::Downstream,
            EdgeKind::Wall,
            EdgeKind::Wall,
        ];
        let on_sides = |v: usize| {
            let (x, y) = (vertices.x[v], vertices.y[v]);
            [x == 0.0, x == self.length, y == 0.0, y == self.width]
        };
        let mut out = Vec::new();
        for verts in vertices_on_cell.rows() {
            let n = verts.len();
            for k in 0..n {
                let (a, b) = (verts[k], verts[(k + 1) % n]);
                let (sa, sb) = (on_sides(a), on_sides(b));
                if let Some(kind) = (0..4).find(|&s| sa[s] && sb[s]).map(|s| KINDS[s]) {
                    out.push(BoundaryEdge {
                        vertices: [a, b],
                        kind,
                    });
                }
            }
        }
        out
    }
}

/// 凸多角形 `polygon`（反時計回り）のうち、母点 `g` と `q` の垂直二等分線より `g` の側を残す
/// （Sutherland–Hodgman 法の1段）。
///
/// 二等分線からの距離が `tol` 以内の頂点は線の上にあるとみなして残し、新しい交点を作らない。
/// 縮退した配置で、丸め誤差程度の長さの辺ができるのを防ぐためである。
fn clip(polygon: &[[f64; 2]], g: [f64; 2], q: [f64; 2], tol: f64) -> Vec<[f64; 2]> {
    let d = distance(g, q);
    let (nx, ny) = ((q[0] - g[0]) / d, (q[1] - g[1]) / d);
    let (mx, my) = (0.5 * (g[0] + q[0]), 0.5 * (g[1] + q[1]));
    // 二等分線からの符号つき距離（q の側が正）
    let side = |p: [f64; 2]| (p[0] - mx) * nx + (p[1] - my) * ny;

    let n = polygon.len();
    let mut out = Vec::with_capacity(n + 1);
    for k in 0..n {
        let (p, s) = (polygon[k], polygon[(k + 1) % n]);
        let (fp, fs) = (side(p), side(s));
        if fp <= tol {
            out.push(p);
        }
        if (fp < -tol && fs > tol) || (fp > tol && fs < -tol) {
            let t = fp / (fp - fs);
            out.push([p[0] + t * (s[0] - p[0]), p[1] + t * (s[1] - p[1])]);
        }
    }
    out
}

fn distance(a: [f64; 2], b: [f64; 2]) -> f64 {
    (b[0] - a[0]).hypot(b[1] - a[1])
}

/// `g` から多角形の頂点までの距離の最大値。
fn max_distance(g: [f64; 2], polygon: &[[f64; 2]]) -> f64 {
    polygon.iter().map(|&p| distance(g, p)).fold(0.0, f64::max)
}

/// `value` が `target` から `tol` 以内なら `target` にする。
fn snap(value: &mut f64, target: f64, tol: f64) {
    if (*value - target).abs() <= tol {
        *value = target;
    }
}

/// 母点を一辺 `size` の正方形のます目に振り分けた索引。
///
/// [`Buckets::near`] は、点を含むます目とその周りの 8 つのます目の母点を返すので、
/// 点から距離 `size` 以内の母点はすべて含まれる。
struct Buckets {
    size: f64,
    nx: usize,
    ny: usize,
    members: Vec<Vec<usize>>,
}

impl Buckets {
    fn new(points: &[[f64; 2]], size: f64, length: f64, width: f64) -> Self {
        let nx = ((length / size).ceil() as usize).max(1);
        let ny = ((width / size).ceil() as usize).max(1);
        let mut out = Self {
            size,
            nx,
            ny,
            members: vec![Vec::new(); nx * ny],
        };
        for (k, &p) in points.iter().enumerate() {
            let (i, j) = out.index(p);
            out.members[j * nx + i].push(k);
        }
        out
    }

    fn index(&self, p: [f64; 2]) -> (usize, usize) {
        let i = ((p[0] / self.size) as usize).min(self.nx - 1);
        let j = ((p[1] / self.size) as usize).min(self.ny - 1);
        (i, j)
    }

    /// 点 `p` の近くの母点の番号。
    fn near(&self, p: [f64; 2]) -> impl Iterator<Item = usize> + '_ {
        let (i, j) = self.index(p);
        let columns = i.saturating_sub(1)..=(i + 1).min(self.nx - 1);
        let rows = j.saturating_sub(1)..=(j + 1).min(self.ny - 1);
        rows.flat_map(move |jj| {
            columns
                .clone()
                .flat_map(move |ii| self.members[jj * self.nx + ii].iter().copied())
        })
    }
}

/// 座標の近い点を同じ頂点とみなして番号をつける（初めて現れた順）。
struct VertexIndex {
    tol: f64,
    x: Vec<f64>,
    y: Vec<f64>,
    /// 一辺 `tol` のます目 → そこに入る頂点。
    grid: HashMap<(i64, i64), Vec<usize>>,
}

impl VertexIndex {
    fn new(tol: f64) -> Self {
        Self {
            tol,
            x: Vec::new(),
            y: Vec::new(),
            grid: HashMap::new(),
        }
    }

    /// 点 `p` から `tol` 以内の頂点があればその番号を、なければ新しい頂点の番号を返す。
    fn insert(&mut self, p: [f64; 2]) -> usize {
        let key = |v: f64| (v / self.tol).floor() as i64;
        let (kx, ky) = (key(p[0]), key(p[1]));
        for jx in kx - 1..=kx + 1 {
            for jy in ky - 1..=ky + 1 {
                if let Some(list) = self.grid.get(&(jx, jy)) {
                    let found = list
                        .iter()
                        .copied()
                        .find(|&v| distance([self.x[v], self.y[v]], p) <= self.tol);
                    if let Some(v) = found {
                        return v;
                    }
                }
            }
        }
        let v = self.x.len();
        self.x.push(p[0]);
        self.y.push(p[1]);
        self.grid.entry((kx, ky)).or_default().push(v);
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(length: f64, width: f64, layout: ChannelLayout) -> ChannelSpec {
        ChannelSpec {
            length,
            width,
            cell_size: 1.0,
            layout,
        }
    }

    /// 境界の種類ごとの辺の長さの合計（上流端、下流端、壁）。
    fn boundary_lengths(mesh: &Mesh) -> [f64; 3] {
        let total = |kind| {
            (0..mesh.n_edges())
                .filter(|&e| mesh.edge_kind[e] == kind)
                .map(|e| mesh.dv_edge[e])
                .sum::<f64>()
        };
        [
            total(EdgeKind::Upstream),
            total(EdgeKind::Downstream),
            total(EdgeKind::Wall),
        ]
    }

    fn assert_close(actual: f64, expected: f64, what: &str) {
        assert!(
            (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
            "{what}: {actual} != {expected}"
        );
    }

    /// V − E + C = 1（穴のない領域）と、幾何の検査。
    fn assert_valid(mesh: &Mesh, length: f64, width: f64) {
        let euler = mesh.n_vertices() as i64 - mesh.n_edges() as i64 + mesh.n_cells() as i64;
        assert_eq!(euler, 1, "Euler characteristic");
        let [up, down, wall] = boundary_lengths(mesh);
        assert_close(up, width, "upstream length");
        assert_close(down, width, "downstream length");
        assert_close(wall, 2.0 * length, "wall length");
        let report = mesh.check_geometry();
        assert_close(report.total_area, length * width, "total area");
        assert!(report.max_orthogonality_error <= 1e-12, "{report:?}");
        assert!(report.max_closure_error <= 1e-12, "{report:?}");
    }

    #[test]
    fn rectangular_cells_are_unit_squares_in_row_order() {
        let mesh = build(&spec(10.0, 4.0, ChannelLayout::Rectangular)).expect("valid channel");
        assert_eq!(mesh.n_cells(), 40);
        assert_eq!(mesh.n_vertices(), 55);
        assert_eq!(mesh.n_edges(), 94);
        assert_valid(&mesh, 10.0, 4.0);
        for c in 0..mesh.n_cells() {
            let (i, j) = (c % 10, c / 10);
            assert_eq!(mesh.x_cell[c], i as f64 + 0.5);
            assert_eq!(mesh.y_cell[c], j as f64 + 0.5);
            assert_eq!(mesh.vertices_on_cell[c].len(), 4);
            assert_close(mesh.area_cell[c], 1.0, "cell area");
        }
        // 境界辺の dc_edge は鏡映点までの距離 = 母点の間隔
        assert!(mesh.dc_edge.iter().all(|&d| (d - 1.0).abs() < 1e-12));
    }

    /// 六角形配置。幅を変えて、上下流端の縮退の3つの場合を確かめる
    /// （Δx = 1 m で、x* = 3Δx/4 − Δy²/Δx の符号が変わる）。
    ///
    /// - W = 4.0（Δy = 0.8）：x* > 0。上下流端で、偶数行のセルどうしが短い辺を共有する
    /// - W = 4.5（Δy = 0.9）：x* < 0。奇数行の端のセルも上下流端に短い境界辺を持つ
    /// - W = 5 × √3/2（Δy = √3/2）：x* = 0。4つのセルの頂点が上下流端の上に乗る
    #[test]
    fn hexagonal_layout_is_valid_in_all_end_configurations() {
        for (width, odd_cells_touch_ends) in [(4.0, false), (4.5, true), (2.5 * 3f64.sqrt(), false)]
        {
            let mesh = build(&spec(10.0, width, ChannelLayout::Hexagonal)).expect("valid channel");
            // 5 行：偶数行 3 × 10 セル、奇数行 2 × 9 セル
            assert_eq!(mesh.n_cells(), 48, "W = {width}");
            assert_valid(&mesh, 10.0, width);

            // 6つの隣を持つセル（最初と最後の行を除き、上下流端から 1.5Δx 以上内側）は六角形。
            // 奇数行の端のセルは左（右）隣の母点がないので五角形になる
            let dy = width / 5.0;
            let interior: Vec<usize> = (0..mesh.n_cells())
                .filter(|&c| {
                    let (x, y) = (mesh.x_cell[c], mesh.y_cell[c]);
                    (dy..=width - dy).contains(&y) && (1.5..=8.5).contains(&x)
                })
                .collect();
            assert_eq!(interior.len(), 8 + 2 * 7, "W = {width}");
            assert!(
                interior
                    .iter()
                    .all(|&c| mesh.vertices_on_cell[c].len() == 6),
                "W = {width}"
            );

            // 側壁の境界辺の dc_edge は行の間隔 Δy（母点は側壁から Δy/2）
            for e in 0..mesh.n_edges() {
                if mesh.edge_kind[e] == EdgeKind::Wall {
                    assert_close(mesh.dc_edge[e], dy, "wall dc_edge");
                }
            }

            // 奇数行の上流端のセル（行 1 の先頭 = セル 10）が上流端に接するかどうか
            let touches = mesh.edges_on_cell[10]
                .iter()
                .any(|&e| mesh.edge_kind[e] == EdgeKind::Upstream);
            assert_eq!(touches, odd_cells_touch_ends, "W = {width}");
        }
    }

    #[test]
    fn single_row_and_single_column_are_valid() {
        for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
            let mesh = build(&spec(10.0, 1.0, layout)).expect("one row");
            assert_eq!(mesh.n_cells(), 10, "{layout:?}");
            assert_valid(&mesh, 10.0, 1.0);
            let mesh = build(&spec(1.0, 4.0, layout)).expect("one column");
            assert_valid(&mesh, 1.0, 4.0);
        }
    }

    #[test]
    fn invalid_dimensions_are_rejected() {
        let mut bad = spec(-1.0, 4.0, ChannelLayout::Hexagonal);
        bad.cell_size = f64::NAN;
        let err = build(&bad).expect_err("invalid dimensions");
        assert_eq!(err.problems().len(), 2, "{err}");
        assert!(err.problems()[0].contains("length"), "{err}");
        assert!(err.problems()[1].contains("cell_size"), "{err}");
    }
}
