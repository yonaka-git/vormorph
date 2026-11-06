"""gmsh を使う格子生成の全体の流れ。

gmsh を読み込めない環境では飛ばす。ただし環境変数 `VORMORPH_REQUIRE_GMSH=1` のとき
（CI と Codespaces）は、共有ライブラリの不足に気づけるように、飛ばさずに失敗させる。
"""

import dataclasses
import os
from pathlib import Path

import numpy as np
import pytest

try:
    import gmsh  # noqa: F401
except (ImportError, OSError) as err:  # 共有ライブラリ（libGLU など）がないと OSError になる
    if os.environ.get("VORMORPH_REQUIRE_GMSH") == "1":
        pytest.fail(f"gmsh is required (VORMORPH_REQUIRE_GMSH=1) but not available: {err}", pytrace=False)
    pytest.skip(f"gmsh is not available: {err}", allow_module_level=True)

from vormesh import Domain, generate, meshfile

EXAMPLES = sorted((Path(__file__).parents[2] / "examples" / "domains").glob("*.toml"))


def boundary_length(mesh, kind):
    vx, vy = mesh.x_vertex, mesh.y_vertex
    return sum(np.hypot(vx[b] - vx[a], vy[b] - vy[a]) for a, b in mesh.boundary[kind])


def channel_with_island(x0=0.0, y0=0.0):
    """長さ 30 m・幅 12 m の水路の中央に半径 2.5 m の円形の島がある領域（左下が (x0, y0)）。"""

    def p(x, y):
        return [x0 + x, y0 + y]

    return Domain.from_dict(
        {
            "title": "channel-island",
            "cell_size": 1.0,
            "outer": [
                {"kind": "wall", "line": [p(0.0, 0.0), p(30.0, 0.0)]},
                {"kind": "downstream", "line": [p(30.0, 0.0), p(30.0, 12.0)]},
                {"kind": "wall", "line": [p(30.0, 12.0), p(0.0, 12.0)]},
                {"kind": "upstream", "line": [p(0.0, 12.0), p(0.0, 0.0)]},
            ],
            "holes": [{"circle": {"center": p(15.0, 6.0), "radius": 2.5}}],
        }
    )


def test_channel_with_a_circular_island():
    q = generate(channel_with_island()).quality
    assert q.max_orthogonality_error < 1e-10
    assert q.area_error < 1e-12
    # 円の島は折れ線で近似されるので、面積は少しだけ大きい
    assert q.domain_area == pytest.approx(30.0 * 12.0 - np.pi * 2.5**2, rel=5e-3)
    assert q.n_cells == pytest.approx(q.domain_area / (np.sqrt(3) / 2), rel=0.02)


def test_domain_far_from_the_origin():
    # 平面直角座標の桁（10⁵ m）に置いた領域でも、原点の近くと同じ品質の格子ができる。gmsh の
    # 分割は平行移動でわずかに変わるので、格子そのものは比べない
    near = generate(channel_with_island()).quality
    far_result = generate(channel_with_island(1.0e5, 1.0e5))
    far = far_result.quality
    assert far.max_orthogonality_error < 1e-10
    assert far.area_error < 1e-12
    assert far.domain_area == pytest.approx(near.domain_area, rel=1e-9)
    assert far.n_cells == pytest.approx(near.n_cells, rel=0.02)
    assert far.max_centroid_offset < 0.05 and near.max_centroid_offset < 0.05
    assert far_result.mesh.x_vertex.min() == pytest.approx(1.0e5, abs=1e-9)


def test_open_boundaries_have_the_given_lengths():
    domain = Domain.load(next(p for p in EXAMPLES if p.stem == "bend"))
    mesh = generate(dataclasses.replace(domain, cell_size=2.0)).mesh
    assert boundary_length(mesh, "upstream") == pytest.approx(10.0, rel=1e-12)
    assert boundary_length(mesh, "downstream") == pytest.approx(10.0, rel=1e-12)


@pytest.mark.parametrize("path", EXAMPLES, ids=lambda p: p.stem)
def test_examples_generate_valid_meshes(path, tmp_path):
    domain = Domain.load(path)
    # テストの時間を抑えるため、格子幅を2倍にする
    result = generate(dataclasses.replace(domain, cell_size=2.0 * domain.cell_size))
    out = tmp_path / f"{path.stem}.mesh.toml"
    meshfile.write(out, result.mesh, domain.title, "hilbert", {"cell_size": 2.0 * domain.cell_size})
    back, header = meshfile.read(out)
    assert back.cells == result.mesh.cells
    assert header["title"] == domain.title
