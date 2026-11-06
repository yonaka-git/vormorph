"""gmsh を使わずに、境界の折れ線と初期母点を直接与えて格子生成の後半を確かめる。"""

import itertools

import numpy as np
import pytest
import shapely

from vormesh import generate_from_discretization, meshfile
from vormesh.domain import BoundaryPiece, Discretization, LloydSettings
from vormesh.hilbert import hilbert_order
from vormesh.voronoi import polygon_area_centroid

H = 1.0


def polyline(points, h=H):
    """折れ線の各辺を間隔 h 以下に分割する。"""
    out = [np.asarray(points[0], dtype=float)]
    for p, q in itertools.pairwise(points):
        p, q = np.asarray(p, dtype=float), np.asarray(q, dtype=float)
        n = max(1, int(np.ceil(np.hypot(*(q - p)) / h - 1e-9)))
        out.extend(p + (q - p) * t for t in np.linspace(0.0, 1.0, n + 1)[1:])
    return np.array(out)


def interior_points(polygon, h=H, seed=0):
    """三角格子状の点に小さな乱れを加え、領域の内側（境界から h/2 以上）に残す。"""
    rng = np.random.default_rng(seed)
    x0, y0, x1, y1 = polygon.bounds
    dy = h * np.sqrt(3.0) / 2.0
    pts = [
        (x + (0.5 * h if k % 2 else 0.0), y)
        for k, y in enumerate(np.arange(y0, y1 + dy, dy))
        for x in np.arange(x0, x1 + h, h)
    ]
    pts = np.array(pts) + rng.uniform(-0.15 * h, 0.15 * h, size=(len(pts), 2))
    keep = shapely.contains_xy(polygon.buffer(-0.5 * h), pts[:, 0], pts[:, 1])
    return pts[keep]


def rectangle_disc(length=10.0, width=6.0):
    pieces = [
        BoundaryPiece("wall", polyline([(0, 0), (length, 0)])),
        BoundaryPiece("downstream", polyline([(length, 0), (length, width)])),
        BoundaryPiece("wall", polyline([(length, width), (0, width)])),
        BoundaryPiece("upstream", polyline([(0, width), (0, 0)])),
    ]
    poly = shapely.Polygon([(0, 0), (length, 0), (length, width), (0, width)])
    return Discretization(pieces, [], interior_points(poly))


def l_shape_disc():
    corners = [(0, 0), (12, 0), (12, 5), (5, 5), (5, 12), (0, 12)]
    pieces = [BoundaryPiece("wall", polyline([corners[k], corners[(k + 1) % 6]])) for k in range(6)]
    pieces[0] = BoundaryPiece("upstream", pieces[0].points)
    pieces[1] = BoundaryPiece("downstream", pieces[1].points)
    return Discretization(pieces, [], interior_points(shapely.Polygon(corners)))


def hole_disc():
    outer = [(0, 0), (14, 0), (14, 10), (0, 10)]
    hole = [(5, 3), (5, 7), (9, 7), (9, 3)]  # 時計回り
    pieces = [BoundaryPiece("wall", polyline([outer[k], outer[(k + 1) % 4]])) for k in range(4)]
    hole_pieces = [BoundaryPiece("wall", polyline([hole[k], hole[(k + 1) % 4]])) for k in range(4)]
    poly = shapely.Polygon(outer, holes=[hole])
    return Discretization(pieces, [hole_pieces], interior_points(poly))


SETTINGS = LloydSettings(iterations=60, tolerance=1e-3)


@pytest.mark.parametrize("make", [rectangle_disc, l_shape_disc, hole_disc])
def test_mesh_satisfies_the_conventions(make):
    result = generate_from_discretization(make(), H, SETTINGS)
    q = result.quality
    assert q.area_error < 1e-12
    assert q.max_orthogonality_error < 1e-10
    assert q.min_boundary_distance > 0.2
    assert q.area_min > 0.3 and q.area_max < 1.6
    # Lloyd 法で母点がセルの重心に近づいている
    assert q.max_centroid_offset < 0.05
    assert result.lloyd.displacement[-1] < result.lloyd.displacement[0]


def test_boundary_kinds_cover_the_whole_boundary():
    result = generate_from_discretization(rectangle_disc(), H, SETTINGS)
    mesh = result.mesh
    vx, vy = mesh.x_vertex, mesh.y_vertex

    def total_length(kind):
        return sum(np.hypot(vx[b] - vx[a], vy[b] - vy[a]) for a, b in mesh.boundary[kind])

    assert total_length("upstream") == pytest.approx(6.0, rel=1e-12)
    assert total_length("downstream") == pytest.approx(6.0, rel=1e-12)
    assert total_length("wall") == pytest.approx(20.0, rel=1e-12)
    # 上流端の辺はすべて x = 0 の上にある
    for a, b in mesh.boundary["upstream"]:
        assert abs(vx[a]) < 1e-9 and abs(vx[b]) < 1e-9


def test_cells_and_vertices_follow_the_hilbert_curve():
    mesh = generate_from_discretization(rectangle_disc(), H, SETTINGS).mesh
    gens = np.column_stack([mesh.x_cell, mesh.y_cell])
    # 書き出したセルの並びは、すでにヒルベルト曲線の順
    assert hilbert_order(gens).tolist() == list(range(mesh.n_cells))
    # 番号の隣り合うセルは近い（多くは隣接セル）
    steps = np.hypot(*np.diff(gens, axis=0).T)
    assert np.median(steps) < 1.2 * H
    # 頂点はセルの順に初めて現れた順に番号がつく
    seen = []
    for cell in mesh.cells:
        seen.extend(v for v in cell if v not in seen)
    assert seen == list(range(mesh.n_vertices))


def test_mesh_file_round_trip(tmp_path):
    mesh = generate_from_discretization(hole_disc(), H, SETTINGS).mesh
    path = tmp_path / "hole.mesh.toml"
    meshfile.write(path, mesh, "hole", "hilbert", {"cell_size": H, "note": 'a "quoted" text'})
    back, header = meshfile.read(path)
    assert header["title"] == "hole" and header["cell_order"] == "hilbert"
    assert header["meta"]["note"] == 'a "quoted" text'
    np.testing.assert_array_equal(back.x_cell, mesh.x_cell)
    np.testing.assert_array_equal(back.y_vertex, mesh.y_vertex)
    assert back.cells == mesh.cells
    assert back.boundary == {k: [tuple(p) for p in v] for k, v in mesh.boundary.items()}


# --- 平面直角座標のような、原点から遠い座標 --------------------------------------------------------

FAR = np.array([1.0e5, 1.0e5])
"""平面直角座標の桁の平行移動 [m]。"""


def dyadic(points, bits=20):
    """座標を 2^-bits の倍数に丸める（FAR だけ平行移動しても丸めが起きない値にする）。"""
    return np.round(np.asarray(points, dtype=float) * 2.0**bits) / 2.0**bits


def dyadic_disc(disc):
    """境界の折れ線と初期母点の座標を `dyadic` で丸めた Discretization。"""

    def move(pieces):
        return [BoundaryPiece(p.kind, dyadic(p.points)) for p in pieces]

    return Discretization(move(disc.outer), [move(h) for h in disc.holes], dyadic(disc.interior_points))


def test_centroid_does_not_lose_digits_far_from_the_origin():
    # 格子幅 0.5 m 程度の六角形（面積 0.22 m²）。座標のまま shoelace の式を使うと、FAR だけ
    # 離れた位置で重心が 0.5 m 程度ずれる
    angles = np.arange(6) * np.pi / 3.0 + 0.1
    hexagon = dyadic(np.column_stack([0.29 * np.cos(angles) + 3.7, 0.29 * np.sin(angles) + 1.3]))
    area, centroid = polygon_area_centroid(hexagon)
    far_area, far_centroid = polygon_area_centroid(hexagon + FAR)
    assert far_area == pytest.approx(area, rel=1e-12)
    np.testing.assert_allclose(far_centroid - FAR, centroid, rtol=0.0, atol=1e-10)
    # 頂点が中心について対称なので、重心は中心に一致する
    np.testing.assert_allclose(centroid, [3.7, 1.3], rtol=0.0, atol=1e-6)


@pytest.mark.parametrize("make", [rectangle_disc, l_shape_disc, hole_disc])
def test_translated_input_gives_the_same_mesh(make):
    # 座標を丸めておき、FAR の平行移動を厳密にする。局所座標での計算がまったく同じになるので、
    # 格子は足し戻しの丸め（10⁻¹¹ m 程度）の範囲で一致する
    disc = dyadic_disc(make())
    near = generate_from_discretization(disc, H, SETTINGS)
    far = generate_from_discretization(disc.translated(FAR), H, SETTINGS)
    assert far.quality == near.quality
    assert far.mesh.cells == near.mesh.cells
    assert far.mesh.boundary == near.mesh.boundary
    for axis, names in enumerate((("x_cell", "x_vertex"), ("y_cell", "y_vertex"))):
        for name in names:
            moved = getattr(far.mesh, name) - FAR[axis]
            np.testing.assert_allclose(moved, getattr(near.mesh, name), rtol=0.0, atol=1e-10)


@pytest.mark.parametrize("make", [rectangle_disc, l_shape_disc, hole_disc])
def test_mesh_far_from_the_origin_satisfies_the_conventions(make):
    # 丸めの起きる平行移動（負の座標を含む）。格子は近くの値と完全には一致しないが、規約と
    # 品質は原点の近くと同じ基準を満たす
    shift = np.array([123456.789, -98765.4321])
    result = generate_from_discretization(make().translated(shift), H, SETTINGS)
    q = result.quality
    assert q.area_error < 1e-12
    assert q.max_orthogonality_error < 1e-10
    assert q.min_boundary_distance > 0.2
    assert q.max_centroid_offset < 0.05
    assert result.lloyd.displacement[-1] < result.lloyd.displacement[0]
    mesh = result.mesh
    assert mesh.x_vertex.min() == pytest.approx(shift[0], abs=1e-9)
    assert mesh.y_vertex.min() == pytest.approx(shift[1], abs=1e-9)
    # どのセルでも、重心がセルの内側にある
    for i in range(mesh.n_cells):
        poly = mesh.cell_polygon(i)
        _, c = polygon_area_centroid(poly)
        assert shapely.Polygon(poly).contains(shapely.Point(*c)), f"centroid of cell {i} is outside"
