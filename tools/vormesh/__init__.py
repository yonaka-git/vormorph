"""vormorph の格子生成ツール。

計算領域の定義（TOML）から、領域の境界に沿ったボロノイ格子を作り、セルをヒルベルト曲線の
順に並べて、vormorph の格子ファイル（TOML）に書き出す。

1. gmsh で境界を格子幅 h の折れ線に分け、内部に初期の母点を置く（`domain`）。
2. 母点のボロノイ領域を計算領域で切り取ってセルとし、Lloyd 法で重心ボロノイ分割に近づける
   （`voronoi`）。
3. 母点をヒルベルト曲線の順に並べ替え、頂点もその順に番号をつける（`hilbert`）。
4. 面積・直交性・境界の規約を検査して書き出す（`voronoi.check`, `meshfile`）。

手順の説明と格子ファイルの仕様は `docs/design/mesh-generation.md` にある。
"""

from __future__ import annotations

from dataclasses import dataclass

from .domain import Discretization, Domain, LloydSettings, discretize
from .hilbert import hilbert_order
from .voronoi import Boundary, LloydResult, PolygonMesh, Quality, build_polygon_mesh, check, lloyd

__all__ = ["Domain", "MeshResult", "generate", "generate_from_discretization"]


@dataclass(frozen=True)
class MeshResult:
    """格子生成の結果。"""

    mesh: PolygonMesh
    quality: Quality
    lloyd: LloydResult
    discretization: Discretization


def generate(domain: Domain, verbose: bool = False) -> MeshResult:
    """領域の定義から格子を作る（gmsh を使う）。"""
    disc = discretize(domain, verbose=verbose)
    return generate_from_discretization(disc, domain.cell_size, domain.lloyd)


def generate_from_discretization(
    disc: Discretization, cell_size: float, settings: LloydSettings
) -> MeshResult:
    """境界の折れ線と初期母点から格子を作る（gmsh を使わない部分）。"""
    boundary = Boundary.from_discretization(disc)
    result = lloyd(disc.interior_points, boundary, cell_size, settings.iterations, settings.tolerance)
    perm = hilbert_order(result.generators)
    generators = result.generators[perm]
    cells = [result.cells[k] for k in perm]
    mesh = build_polygon_mesh(generators, cells, boundary, cell_size)
    quality = check(mesh, boundary.polygon.area, cell_size)
    return MeshResult(mesh, quality, result, disc)
