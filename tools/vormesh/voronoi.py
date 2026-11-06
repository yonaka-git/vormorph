"""境界に沿うボロノイ格子：境界の母点の対、Lloyd 法、多角形格子の組み立てと検査。

## 境界の扱い（鏡映対の母点）

境界の折れ線の節点 b_k（gmsh が格子幅 h で分割したもの）を中心に半径 r_k の円を描く。
隣り合う円 (b_k, r_k) と (b_{k+1}, r_{k+1}) の2つの交点は、線分 b_k b_{k+1} の直線に関して
対称なので、領域の内側の交点を境界辺 k の母点、外側の交点をその鏡映点とする
（Abdelkader et al., 2018 の VoroCrust と同じ考え方の2次元版）。

- 母点と鏡映点の二等分線は境界の直線そのものなので、線分 b_k b_{k+1} がそのままボロノイ辺に
  なり、節点 b_k はボロノイ頂点になる。境界辺でも直交性が厳密に成り立ち、Rust 側の規約
  （境界辺の `dc_edge` は鏡映点までの距離）とそのまま整合する。
- これが成り立つ条件は、どの円の内側にも、その円の上の4点（隣り合う2つの対）以外の母点が
  ないこと。内部の母点は円の外に置き、境界の母点についても検査する。
- 半径は r_k = 0.66 × (隣り合う2辺の長さの平均)。一様な分割では母点が境界から約 0.43h の
  位置に来る（境界に平行な六角形の列と同じ間隔）。凸な角は 90° まで、凹な角は制限なし。

## 重心ボロノイ分割

領域を境界から (√3/4) h だけ縮めた多角形の、gmsh の三角形分割の内部の節点を初期の母点とする
（`vormesh.domain`）。gmsh の分割は正六角形の格子より数％密なので、母点の総数が正六角形の
格子の個数 A / (√3/2 h²) になるまで、混み合った内部の母点から間引く。そのうえで、内部の母点
だけを各セルの重心へ移す操作（Lloyd 法）を繰り返す。境界の母点は動かさない。収束すると重心
ボロノイ分割（CVT; Du et al., 1999）に近づく。
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import shapely
from scipy.spatial import KDTree, Voronoi

from .domain import KINDS, Discretization

BALL_RATIO = 0.66
"""境界の節点の円の半径と、隣り合う2辺の長さの平均との比。"""


class MeshError(RuntimeError):
    """格子を組み立てられない、または検査に通らない。"""


# --- 境界と境界の母点 -----------------------------------------------------------------------------


@dataclass(frozen=True)
class Boundary:
    """計算領域の境界。`rings[r]` は閉曲線 r の節点、`kinds[r][k]` は辺 b_k → b_{k+1} の種類。"""

    polygon: shapely.Polygon
    rings: list[np.ndarray]
    kinds: list[list[str]]

    @classmethod
    def from_discretization(cls, disc: Discretization) -> Boundary:
        loops = [disc.outer, *disc.holes]
        rings = [np.vstack([p.points[:-1] for p in loop]) for loop in loops]
        kinds = [[p.kind for p in loop for _ in range(len(p.points) - 1)] for loop in loops]
        polygon = shapely.Polygon(rings[0], holes=rings[1:])
        if not polygon.is_valid:
            reason = shapely.is_valid_reason(polygon)
            raise MeshError(f"the domain boundary is not a valid polygon: {reason}")
        shapely.prepare(polygon)
        return cls(polygon, rings, kinds)

    @property
    def nodes(self) -> np.ndarray:
        return np.vstack(self.rings)


@dataclass(frozen=True)
class BoundarySeeds:
    """境界辺ごとの母点（領域の内側）と鏡映点（外側）、節点の円の中心と半径。"""

    inner: np.ndarray
    outer: np.ndarray
    centers: np.ndarray
    radii: np.ndarray
    owner: np.ndarray
    """境界辺 k の母点が乗っている2つの円の番号 (k, k+1)。形は (辺の数, 2)。"""

    def outside_balls(self, points: np.ndarray, margin: float = 1.0e-9) -> np.ndarray:
        """どの円の内側にもない点なら True。"""
        tree = KDTree(self.centers)
        ok = np.ones(len(points), dtype=bool)
        for i, near in enumerate(tree.query_ball_point(points, float(self.radii.max()))):
            if near:
                d = np.hypot(*(self.centers[near] - points[i]).T)
                ok[i] = bool((d >= self.radii[near] * (1.0 - margin)).all())
        return ok

    def push_out(self, points: np.ndarray) -> np.ndarray:
        """円の内側に入った点を、円の外へ半径方向に押し出す。"""
        pts = points.copy()
        tree = KDTree(self.centers)
        for _ in range(3):
            moved = False
            for i, near in enumerate(tree.query_ball_point(pts, float(self.radii.max()))):
                for c in near:
                    v = pts[i] - self.centers[c]
                    d = float(np.hypot(*v))
                    if d < self.radii[c]:
                        if d == 0.0:
                            v, d = np.array([1.0, 0.0]), 1.0
                        pts[i] = self.centers[c] + v / d * self.radii[c] * (1.0 + 1.0e-3)
                        moved = True
            if not moved:
                break
        return pts


def boundary_seeds(boundary: Boundary, ratio: float = BALL_RATIO) -> BoundarySeeds:
    """境界の節点の円の交点から、境界辺ごとの母点と鏡映点を作る。"""
    inner, outer, centers, radii, owner = [], [], [], [], []
    offset = 0
    for ring in boundary.rings:
        n = len(ring)
        nxt = np.roll(ring, -1, axis=0)
        length = np.hypot(*(nxt - ring).T)
        r = ratio * 0.5 * (np.roll(length, 1) + length)
        r_next = np.roll(r, -1)
        a = (r**2 - r_next**2 + length**2) / (2.0 * length)
        h2 = r**2 - a**2
        if (h2 <= 0.0).any():
            k = int(np.argmin(h2))
            raise MeshError(f"boundary circles do not intersect near {ring[k].tolist()}")
        u = (nxt - ring) / length[:, None]
        normal = np.column_stack([-u[:, 1], u[:, 0]])
        base = ring + a[:, None] * u
        left = base + np.sqrt(h2)[:, None] * normal
        right = base - np.sqrt(h2)[:, None] * normal
        left_in = shapely.contains_xy(boundary.polygon, left[:, 0], left[:, 1])
        right_in = shapely.contains_xy(boundary.polygon, right[:, 0], right[:, 1])
        if not (left_in ^ right_in).all():
            k = int(np.argmin(left_in ^ right_in))
            raise MeshError(f"the domain is too narrow or too sharp near {ring[k].tolist()}")
        inner.append(np.where(left_in[:, None], left, right))
        outer.append(np.where(left_in[:, None], right, left))
        centers.append(ring)
        radii.append(r)
        idx = np.arange(n)
        owner.append(np.column_stack([idx, (idx + 1) % n]) + offset)
        offset += n
    seeds = BoundarySeeds(
        np.vstack(inner), np.vstack(outer), np.vstack(centers), np.concatenate(radii), np.vstack(owner)
    )
    _check_empty_balls(seeds)
    return seeds


def _check_empty_balls(seeds: BoundarySeeds) -> None:
    """境界の母点と鏡映点が、自分の乗る2つの円以外の円の内側にないことを確かめる。"""
    tree = KDTree(seeds.centers)
    for points in (seeds.inner, seeds.outer):
        for k, near in enumerate(tree.query_ball_point(points, float(seeds.radii.max()))):
            for c in near:
                if c in seeds.owner[k]:
                    continue
                if np.hypot(*(points[k] - seeds.centers[c])) < seeds.radii[c] * (1.0 - 1.0e-9):
                    raise MeshError(
                        f"the boundary is too sharply bent near {seeds.centers[c].tolist()} "
                        "(convex corners must be 90 degrees or more); use a smaller cell_size"
                    )


# --- ボロノイ図と Lloyd 法 ------------------------------------------------------------------------


HEXAGON_AREA = np.sqrt(3.0) / 2.0
"""母点の間隔 h の正六角形のセルの面積 / h²。"""


def thin_out(points: np.ndarray, fixed: np.ndarray, n_remove: int, cell_size: float) -> np.ndarray:
    """最も近い点との距離が小さい（混み合った）点から順に `n_remove` 個を取り除く。

    `fixed` は取り除かない点（境界の母点）。1回の走査では、互いに近い点を同時に取り除かない。
    """
    pts = np.asarray(points, dtype=float)
    remaining = n_remove
    while remaining > 0 and len(pts) > 0:
        allpts = np.vstack([fixed, pts])
        dist, _ = KDTree(allpts).query(pts, k=2)
        crowd = dist[:, 1]
        tree = KDTree(pts)
        removed = np.zeros(len(pts), dtype=bool)
        blocked = np.zeros(len(pts), dtype=bool)
        for i in np.argsort(crowd, kind="stable"):
            if remaining == 0:
                break
            if blocked[i]:
                continue
            removed[i] = True
            remaining -= 1
            blocked[tree.query_ball_point(pts[i], 1.5 * cell_size)] = True
        pts = pts[~removed]
    return pts


def voronoi_cells(real: np.ndarray, ghosts: np.ndarray) -> list[np.ndarray]:
    """`real` の各点のボロノイ領域（反時計回りの頂点列）を返す。`ghosts` は鏡映点。"""
    pts = np.vstack([real, ghosts])
    lower, upper = pts.min(axis=0), pts.max(axis=0)
    center = 0.5 * (lower + upper)
    radius = 10.0 * float(np.hypot(*(upper - lower))) + 1.0
    angles = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
    frame = center + radius * np.column_stack([np.cos(angles), np.sin(angles)])
    vor = Voronoi(np.vstack([pts, frame]))
    cells = []
    for i in range(len(real)):
        region = vor.regions[vor.point_region[i]]
        if -1 in region or len(region) < 3:
            raise MeshError(f"Voronoi region of generator {i} is unbounded")
        verts = vor.vertices[region]
        angle = np.arctan2(verts[:, 1] - real[i, 1], verts[:, 0] - real[i, 0])
        cells.append(verts[np.argsort(angle)])
    return cells


def polygon_area_centroid(poly: np.ndarray) -> tuple[float, np.ndarray]:
    """反時計回りの多角形の面積と重心。

    頂点の平均を原点にした座標で shoelace の式を計算し、重心に原点を足し戻す。座標のまま
    計算すると、原点から X 離れた面積 A の多角形で、重心におよそ εX³/A（ε は倍精度の丸めの
    単位）の桁落ちが生じる（X = 10⁵ m、A = 1 m² で約 0.1 m）。
    """
    origin = poly.mean(axis=0)
    x, y = poly[:, 0] - origin[0], poly[:, 1] - origin[1]
    xn, yn = np.roll(x, -1), np.roll(y, -1)
    cross = x * yn - xn * y
    area = 0.5 * float(cross.sum())
    cx = float(((x + xn) * cross).sum()) / (6.0 * area)
    cy = float(((y + yn) * cross).sum()) / (6.0 * area)
    return area, origin + np.array([cx, cy])


@dataclass
class LloydResult:
    """Lloyd 法の結果。

    `displacement` と `max_displacement` は、各反復での内部の母点の移動量の二乗平均平方根と
    最大値 / h。収束の判定には二乗平均平方根を使う（円に接して行き来する少数の母点で、
    最大値は下がりきらないことがあるため）。
    """

    generators: np.ndarray
    cells: list[np.ndarray]
    n_boundary: int
    n_removed: int = 0
    displacement: list[float] = field(default_factory=list)
    max_displacement: list[float] = field(default_factory=list)

    @property
    def iterations(self) -> int:
        return len(self.displacement)

    def translated(self, offset: np.ndarray) -> LloydResult:
        """母点とセルの頂点を `offset` だけ平行移動した結果。"""
        shift = np.asarray(offset, dtype=float)
        return replace(self, generators=self.generators + shift, cells=[c + shift for c in self.cells])


def lloyd(
    interior: np.ndarray,
    boundary: Boundary,
    cell_size: float,
    iterations: int,
    tolerance: float,
) -> LloydResult:
    """境界の母点を固定し、内部の母点をセルの重心へ移す操作を、移動量の二乗平均平方根が
    `tolerance * cell_size` を下回るまで（最大 `iterations` 回）繰り返す。"""
    seeds = boundary_seeds(boundary)
    pts = np.asarray(interior, dtype=float)
    keep = shapely.contains_xy(boundary.polygon, pts[:, 0], pts[:, 1]) & seeds.outside_balls(pts)
    pts = pts[keep].copy()
    m = len(seeds.inner)
    target = round(boundary.polygon.area / (HEXAGON_AREA * cell_size**2))
    removed = max(0, m + len(pts) - target)
    if removed:
        pts = thin_out(pts, seeds.inner, removed, cell_size)
    rms_history: list[float] = []
    max_history: list[float] = []
    for _ in range(iterations):
        cells = voronoi_cells(np.vstack([seeds.inner, pts]), seeds.outer)
        centroids = np.array([polygon_area_centroid(c)[1] for c in cells[m:]]).reshape(-1, 2)
        centroids = seeds.push_out(centroids)
        ok = shapely.contains_xy(boundary.polygon, centroids[:, 0], centroids[:, 1])
        new = np.where(ok[:, None], centroids, pts)
        moved = np.hypot(*(new - pts).T) / cell_size if len(pts) else np.zeros(1)
        pts = new
        rms_history.append(float(np.sqrt(np.mean(moved**2))))
        max_history.append(float(moved.max()))
        if rms_history[-1] < tolerance:
            break
    generators = np.vstack([seeds.inner, pts])
    cells = voronoi_cells(generators, seeds.outer)
    return LloydResult(generators, cells, m, removed, rms_history, max_history)


# --- 多角形格子 ---------------------------------------------------------------------------------


@dataclass
class PolygonMesh:
    """多角形のセルで表したボロノイ格子（Rust 側の `mesh::voronoi::PolygonMesh` に対応）。

    - `cells[i]` はセル i の頂点番号（反時計回り）。辺 k は `cells[i][k]` → `cells[i][k + 1]`。
    - `boundary[kind]` は境界辺の両端の頂点番号の組（セルから見て反時計回りの向き）。
    """

    x_cell: np.ndarray
    y_cell: np.ndarray
    x_vertex: np.ndarray
    y_vertex: np.ndarray
    cells: list[list[int]]
    boundary: dict[str, list[tuple[int, int]]]

    @property
    def n_cells(self) -> int:
        return len(self.x_cell)

    @property
    def n_vertices(self) -> int:
        return len(self.x_vertex)

    def cell_polygon(self, i: int) -> np.ndarray:
        idx = self.cells[i]
        return np.column_stack([self.x_vertex[idx], self.y_vertex[idx]])

    def translated(self, offset: np.ndarray) -> PolygonMesh:
        """母点と頂点を `offset` だけ平行移動した格子（接続はそのまま）。"""
        dx, dy = float(offset[0]), float(offset[1])
        return replace(
            self,
            x_cell=self.x_cell + dx,
            y_cell=self.y_cell + dy,
            x_vertex=self.x_vertex + dx,
            y_vertex=self.y_vertex + dy,
        )

    def directed_edges(self):
        """(頂点 a, 頂点 b, セル) を、セルの番号順・反時計回りに返す。"""
        for i, cell in enumerate(self.cells):
            for a, b in zip(cell, cell[1:] + cell[:1], strict=True):
                yield a, b, i

    def interior_edges(self) -> list[tuple[int, int, int, int]]:
        """内部の辺を (頂点 a, 頂点 b, a→b を持つセル, b→a を持つセル) で、出現順に返す。"""
        owner = {(a, b): i for a, b, i in self.directed_edges()}
        out = []
        for (a, b), i in owner.items():
            j = owner.get((b, a))
            if j is not None and i < j:
                out.append((a, b, i, j))
        return out


def build_polygon_mesh(
    generators: np.ndarray,
    cells: list[np.ndarray],
    boundary: Boundary,
    cell_size: float,
) -> PolygonMesh:
    """ボロノイ領域の頂点を共有する多角形格子を作る。

    頂点の番号は、セルを番号順にたどって最初に現れた順につける。セルをヒルベルト曲線の順に
    並べておけば、頂点（と Rust 側で作る辺）の番号もその順に近くなる。
    """
    coords = np.vstack(cells)
    offsets = np.cumsum([0] + [len(c) for c in cells])

    # 縮退（同じ円の上の4点）で重なったボロノイ頂点を1つにまとめる
    tol = 1.0e-8 * cell_size
    parent = np.arange(len(coords))

    def find(k: int) -> int:
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    for a, b in KDTree(coords).query_pairs(tol):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    roots = np.array([find(k) for k in range(len(coords))])

    # 境界の節点に一致する頂点は、節点の座標そのものにそろえる
    exact = coords.copy()
    nodes = boundary.nodes
    dist, hit = KDTree(coords).query(nodes)
    if (dist > 1.0e3 * tol).any():
        k = int(np.argmax(dist))
        raise MeshError(f"boundary node {nodes[k].tolist()} is not a Voronoi vertex")
    node_root = roots[hit]
    exact[node_root] = nodes

    number: dict[int, int] = {}
    vertex_xy: list[np.ndarray] = []
    polys: list[list[int]] = []
    for i in range(len(cells)):
        poly: list[int] = []
        for k in range(offsets[i], offsets[i + 1]):
            r = int(roots[k])
            if r not in number:
                number[r] = len(vertex_xy)
                vertex_xy.append(exact[r])
            v = number[r]
            if not poly or poly[-1] != v:
                poly.append(v)
        while len(poly) > 1 and poly[0] == poly[-1]:
            poly.pop()
        if len(poly) < 3:
            raise MeshError(f"cell {i} degenerated to fewer than three vertices")
        polys.append(poly)
    vxy = np.array(vertex_xy)

    # 境界の辺（節点の組）と種類
    ring_edges: dict[tuple[int, int], str] = {}
    start = 0
    for ring, kinds in zip(boundary.rings, boundary.kinds, strict=True):
        ids = [number[int(r)] for r in node_root[start : start + len(ring)]]
        for k, kind in enumerate(kinds):
            a, b = ids[k], ids[(k + 1) % len(ids)]
            ring_edges[(min(a, b), max(a, b))] = kind
        start += len(ring)

    mesh = PolygonMesh(
        x_cell=generators[:, 0].copy(),
        y_cell=generators[:, 1].copy(),
        x_vertex=vxy[:, 0].copy(),
        y_vertex=vxy[:, 1].copy(),
        cells=polys,
        boundary={k: [] for k in KINDS},
    )
    directed = {(a, b) for a, b, _ in mesh.directed_edges()}
    for a, b, _ in mesh.directed_edges():
        if (b, a) in directed:
            continue
        kind = ring_edges.pop((min(a, b), max(a, b)), None)
        if kind is None:
            raise MeshError(
                f"edge {vxy[a].tolist()} -> {vxy[b].tolist()} has no neighbor but is not a "
                "boundary segment (the cells do not conform to the boundary)"
            )
        mesh.boundary[kind].append((a, b))
    if ring_edges:
        (a, b), _ = next(iter(ring_edges.items()))
        raise MeshError(f"boundary segment {vxy[a].tolist()} - {vxy[b].tolist()} is not a cell edge")
    return mesh


# --- 検査 -------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Quality:
    """格子の検査結果。長さは格子幅 h で、面積は h² で無次元化している。"""

    n_cells: int
    n_vertices: int
    n_interior_edges: int
    n_boundary_edges: dict[str, int]
    total_area: float
    domain_area: float
    area_min: float
    area_max: float
    edge_min: float
    max_orthogonality_error: float
    min_boundary_distance: float
    max_centroid_offset: float
    """内部のセル（境界辺を持たないセル）での、重心と母点の距離の最大値 / h。"""
    max_centroid_offset_boundary: float
    """境界のセルでの同じ量。境界の母点は動かさないので、内部より大きい。"""

    @property
    def area_error(self) -> float:
        """セル面積の和と領域の面積の相対誤差。"""
        return abs(self.total_area - self.domain_area) / self.domain_area

    def lines(self) -> list[str]:
        b = ", ".join(f"{k} {v}" for k, v in self.n_boundary_edges.items())
        return [
            f"cells          : {self.n_cells}",
            f"vertices       : {self.n_vertices}",
            f"edges          : {self.n_interior_edges} interior, {b}",
            f"area           : {self.total_area:.6g} m^2 (relative error {self.area_error:.1e})",
            f"cell area / h^2: {self.area_min:.3f} .. {self.area_max:.3f}",
            f"min edge / h   : {self.edge_min:.3g}",
            f"orthogonality  : max |cos| = {self.max_orthogonality_error:.1e}",
            f"boundary       : min generator-to-edge distance / h = {self.min_boundary_distance:.3f}",
            f"centroid offset: max |x_c - x_g| / h = {self.max_centroid_offset:.3g} (interior), "
            f"{self.max_centroid_offset_boundary:.3g} (boundary cells)",
        ]


SHORT_EDGE = 1.0e-6
"""直交誤差の評価から除く、ごく短い辺の長さ / h（縮退に近い4点の円で生じる）。"""


def check(mesh: PolygonMesh, domain_area: float, cell_size: float) -> Quality:
    """格子を検査し、規約を満たさなければ `MeshError` を送出する。"""
    h = cell_size
    vx, vy = mesh.x_vertex, mesh.y_vertex
    areas = np.empty(mesh.n_cells)
    offsets = np.empty(mesh.n_cells)
    for i in range(mesh.n_cells):
        poly = mesh.cell_polygon(i)
        a, c = polygon_area_centroid(poly)
        if not a > 0.0:
            raise MeshError(f"cell {i} is not counterclockwise (area {a:.3g})")
        areas[i] = a
        offsets[i] = np.hypot(c[0] - mesh.x_cell[i], c[1] - mesh.y_cell[i])
        if not shapely.Polygon(poly).contains(shapely.Point(mesh.x_cell[i], mesh.y_cell[i])):
            raise MeshError(f"generator {i} is not inside its cell")

    edges = mesh.interior_edges()
    lengths = [float(np.hypot(vx[b] - vx[a], vy[b] - vy[a])) for a, b, _, _ in edges]
    cos_max = 0.0
    for (a, b, i, j), length in zip(edges, lengths, strict=True):
        if length < SHORT_EDGE * h:
            continue
        gx, gy = mesh.x_cell[j] - mesh.x_cell[i], mesh.y_cell[j] - mesh.y_cell[i]
        tx, ty = vx[b] - vx[a], vy[b] - vy[a]
        cos_max = max(cos_max, abs(gx * tx + gy * ty) / (np.hypot(gx, gy) * length))

    # 境界辺：母点が辺の内側にあること（鏡映点までの距離 dc = 2d > 0）
    owner = {(a, b): i for a, b, i in mesh.directed_edges()}
    d_min = np.inf
    for pairs in mesh.boundary.values():
        for a, b in pairs:
            i = owner[(a, b)]
            tx, ty = vx[b] - vx[a], vy[b] - vy[a]
            length = float(np.hypot(tx, ty))
            lengths.append(length)
            d = (tx * (mesh.y_cell[i] - vy[a]) - ty * (mesh.x_cell[i] - vx[a])) / length
            d_min = min(d_min, d)
    if not d_min > 0.0:
        raise MeshError("a generator lies on or outside one of its boundary edges")
    on_boundary = np.zeros(mesh.n_cells, dtype=bool)
    for pairs in mesh.boundary.values():
        for a, b in pairs:
            on_boundary[owner[(a, b)]] = True

    n_directed = sum(len(c) for c in mesh.cells)
    n_boundary = sum(len(v) for v in mesh.boundary.values())
    if 2 * len(edges) + n_boundary != n_directed:
        raise MeshError("cells do not share edges consistently (T-junctions)")

    quality = Quality(
        n_cells=mesh.n_cells,
        n_vertices=mesh.n_vertices,
        n_interior_edges=len(edges),
        n_boundary_edges={k: len(v) for k, v in mesh.boundary.items()},
        total_area=float(areas.sum()),
        domain_area=float(domain_area),
        area_min=float(areas.min()) / h**2,
        area_max=float(areas.max()) / h**2,
        edge_min=float(min(lengths)) / h,
        max_orthogonality_error=float(cos_max),
        min_boundary_distance=float(d_min) / h,
        max_centroid_offset=float(offsets[~on_boundary].max(initial=0.0)) / h,
        max_centroid_offset_boundary=float(offsets[on_boundary].max(initial=0.0)) / h,
    )
    if quality.area_error > 1.0e-9:
        raise MeshError(f"cell areas do not add up to the domain area ({quality.area_error:.1e})")
    if quality.max_orthogonality_error > 1.0e-8:
        raise MeshError(f"interior edges are not orthogonal (max |cos| {cos_max:.1e})")
    return quality
