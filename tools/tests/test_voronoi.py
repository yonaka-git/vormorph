"""gmsh を使わずに、境界の折れ線と初期母点を直接与えて格子生成の後半を確かめる。"""

import itertools

import numpy as np
import pytest
import shapely

from vormesh import generate_from_discretization, meshfile
from vormesh.domain import BoundaryPiece, Discretization, LloydSettings
from vormesh.hilbert import hilbert_order

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
