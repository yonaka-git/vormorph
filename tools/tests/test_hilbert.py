import numpy as np
import pytest

from vormesh.hilbert import hilbert_index, hilbert_order


def test_order_one_visits_quadrants_in_u_shape():
    # (0,0) → (0,1) → (1,1) → (1,0)
    d = hilbert_index(np.array([0, 0, 1, 1]), np.array([0, 1, 1, 0]), order=1)
    assert d.tolist() == [0, 1, 2, 3]


@pytest.mark.parametrize("order", [2, 3, 5])
def test_index_is_a_bijection_and_consecutive_cells_are_adjacent(order):
    n = 1 << order
    ix, iy = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    d = hilbert_index(ix.ravel(), iy.ravel(), order)
    assert sorted(d.tolist()) == list(range(n * n))
    perm = np.argsort(d)
    steps = np.abs(np.diff(ix.ravel()[perm])) + np.abs(np.diff(iy.ravel()[perm]))
    assert (steps == 1).all()
    assert (ix.ravel()[perm[-1]], iy.ravel()[perm[-1]]) == (n - 1, 0)


def test_out_of_range_coordinates_are_rejected():
    with pytest.raises(ValueError):
        hilbert_index(np.array([4]), np.array([0]), order=2)


def test_order_of_points_follows_the_curve_and_is_stable():
    # 4×4 の格子の中心に置いた点（行優先の順）を並べ替えると、曲線の順になる
    n = 4
    xs, ys = np.meshgrid(np.arange(n) + 0.5, np.arange(n) + 0.5)
    pts = np.column_stack([xs.ravel(), ys.ravel()])
    perm = hilbert_order(pts, order=2)
    walk = pts[perm]
    assert (np.abs(np.diff(walk, axis=0)).sum(axis=1) == 1.0).all()
    # 同じ点が重なっていても元の順を保つ
    same = np.zeros((3, 2))
    assert hilbert_order(same).tolist() == [0, 1, 2]
