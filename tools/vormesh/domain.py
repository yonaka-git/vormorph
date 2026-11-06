"""計算領域の定義（TOML）の読み込みと、gmsh による境界の分割・初期母点の生成。

領域は、外周を反時計回りにたどる区間の列と、島などの穴で与える。

```toml
title = "bend180"
cell_size = 1.0                       # 目標の格子幅（母点の間隔）[m]

[[outer]]                             # 外周の区間（反時計回りの順）
kind = "upstream"                     # wall / upstream / downstream
line = [[25.0, -40.0], [35.0, -40.0]] # 折れ線

[[outer]]
kind = "wall"
arc = { center = [0.0, 0.0], radius = 35.0, start = 0.0, end = 180.0 }  # 円弧（度、start → end）

[[holes]]                             # 穴（境界の種類は wall）
circle = { center = [30.0, 15.0], radius = 5.0 }   # または polygon = [[x, y], ...]

[lloyd]                               # 省略可
iterations = 100                      # Lloyd 法の反復回数の上限
tolerance = 1e-3                      # 収束判定：母点の移動量の二乗平均平方根 / cell_size
```

gmsh は次の2つに使う。ボロノイ格子そのものは `vormesh.voronoi` で作る。

1. 領域の境界を、格子幅 h に合わせて折れ線に分割する（境界の節点）。
2. 領域を境界から (√3/4) h だけ内側に縮めた多角形を Frontal-Delaunay 法で三角形分割し、
   その内部の節点を初期の母点にする。境界の母点の列は境界から約 (√3/4) h にあるので、
   内部の母点はその次の列（正六角形の格子の列の間隔 (√3/2) h だけ内側）から始まる。
"""

from __future__ import annotations

import math
import tomllib
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

KINDS = ("wall", "upstream", "downstream")
"""境界の種類。Rust 側の `mesh::EdgeKind` の Wall / Upstream / Downstream に対応する。"""

INTERIOR_OFFSET = math.sqrt(3.0) / 4.0
"""初期の母点を置く領域を、境界から内側へ縮める幅 / h（境界の母点の列の位置）。"""

MAX_ARC_DEGREES = 90.0
"""gmsh の組み込みカーネルの円弧は π 未満に限られるので、円弧をこの角度以下に分ける。"""


class DomainError(ValueError):
    """領域の定義の誤り。"""


@dataclass(frozen=True)
class Segment:
    """外周または穴の境界の1区間。`points` は始点から終点までの折れ線（円弧は中心も持つ）。"""

    kind: str
    points: np.ndarray
    arc_center: np.ndarray | None = None

    @property
    def start(self) -> np.ndarray:
        return self.points[0]

    @property
    def end(self) -> np.ndarray:
        return self.points[-1]


@dataclass(frozen=True)
class LloydSettings:
    """Lloyd 法（重心ボロノイ分割への反復）の設定。"""

    iterations: int = 100
    tolerance: float = 1.0e-3


@dataclass(frozen=True)
class Domain:
    """計算領域。`outer` は反時計回りの閉じた区間の列、`holes` は穴ごとの区間の列。"""

    title: str
    cell_size: float
    outer: list[Segment]
    holes: list[list[Segment]] = field(default_factory=list)
    lloyd: LloydSettings = LloydSettings()

    @classmethod
    def load(cls, path: Path) -> Domain:
        with open(path, "rb") as f:
            data = tomllib.load(f)
        return cls.from_dict(data, source=str(path))

    @classmethod
    def from_dict(cls, data: dict, source: str = "<domain>") -> Domain:
        known = {"title", "cell_size", "outer", "holes", "lloyd"}
        unknown = sorted(set(data) - known)
        if unknown:
            raise DomainError(f"{source}: unknown keys: {', '.join(unknown)}")
        try:
            title = str(data["title"])
            cell_size = float(data["cell_size"])
            outer_spec = data["outer"]
        except KeyError as err:
            raise DomainError(f"{source}: missing key {err}") from None
        if not cell_size > 0.0:
            raise DomainError(f"{source}: cell_size must be positive")

        outer = [_segment(s, f"{source}: outer[{k}]", cell_size) for k, s in enumerate(outer_spec)]
        holes = [_hole(h, f"{source}: holes[{k}]") for k, h in enumerate(data.get("holes", []))]
        lloyd_spec = data.get("lloyd", {})
        unknown = sorted(set(lloyd_spec) - {"iterations", "tolerance"})
        if unknown:
            raise DomainError(f"{source}: unknown keys in [lloyd]: {', '.join(unknown)}")
        lloyd = LloydSettings(**lloyd_spec)
        domain = cls(title, cell_size, outer, holes, lloyd)
        domain.validate(source)
        return domain

    def validate(self, source: str = "<domain>") -> None:
        """区間がつながって閉じていること、外周が反時計回りであることを確かめる。"""
        if self.lloyd.iterations < 0 or not self.lloyd.tolerance > 0.0:
            raise DomainError(f"{source}: lloyd.iterations must be >= 0 and tolerance > 0")
        loops = [("outer", self.outer)] + [(f"holes[{k}]", h) for k, h in enumerate(self.holes)]
        for name, loop in loops:
            if not loop:
                raise DomainError(f"{source}: {name} is empty")
            tol = 1.0e-9 * max(1.0, _extent(loop))
            for k, seg in enumerate(loop):
                nxt = loop[(k + 1) % len(loop)]
                gap = float(np.hypot(*(seg.end - nxt.start)))
                if gap > tol:
                    raise DomainError(
                        f"{source}: {name}[{k}] ends at {seg.end.tolist()} but the next segment "
                        f"starts at {nxt.start.tolist()} (gap {gap:.3g})"
                    )
        if _signed_area(self.outer) <= 0.0:
            raise DomainError(f"{source}: outer must be counterclockwise (interior on the left)")


def _segment(spec: dict, where: str, cell_size: float) -> Segment:
    kind = spec.get("kind", "wall")
    if kind not in KINDS:
        raise DomainError(f"{where}: kind must be one of {KINDS}, got {kind!r}")
    shapes = [k for k in ("line", "arc") if k in spec]
    extra = sorted(set(spec) - {"kind", "line", "arc"})
    if len(shapes) != 1 or extra:
        raise DomainError(f"{where}: give exactly one of 'line' or 'arc' (and 'kind')")
    if "line" in spec:
        points = np.asarray(spec["line"], dtype=float)
        if points.ndim != 2 or points.shape[1] != 2 or len(points) < 2:
            raise DomainError(f"{where}: line must be a list of at least two [x, y] points")
        return Segment(kind, points)
    arc = spec["arc"]
    try:
        center = np.asarray(arc["center"], dtype=float)
        radius = float(arc["radius"])
        start, end = math.radians(float(arc["start"])), math.radians(float(arc["end"]))
    except KeyError as err:
        raise DomainError(f"{where}: arc needs center, radius, start and end ({err})") from None
    if not radius > 0.0 or start == end:
        raise DomainError(f"{where}: arc needs radius > 0 and start != end")
    pieces = max(1, math.ceil(abs(end - start) / math.radians(MAX_ARC_DEGREES) - 1.0e-12))
    angles = np.linspace(start, end, pieces + 1)
    points = center + radius * np.column_stack([np.cos(angles), np.sin(angles)])
    return Segment(kind, points, arc_center=center)


def _hole(spec: dict, where: str) -> list[Segment]:
    extra = sorted(set(spec) - {"circle", "polygon", "kind"})
    if extra or ("circle" in spec) == ("polygon" in spec):
        raise DomainError(f"{where}: give exactly one of 'circle' or 'polygon'")
    if spec.get("kind", "wall") != "wall":
        raise DomainError(f"{where}: holes must be walls")
    if "circle" in spec:
        c = spec["circle"]
        arc = {"center": c["center"], "radius": c["radius"], "start": 0.0, "end": 360.0}
        seg = _segment({"kind": "wall", "arc": arc}, where, 1.0)
        return [seg]
    points = np.asarray(spec["polygon"], dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 3:
        raise DomainError(f"{where}: polygon must have at least three [x, y] points")
    closed = np.vstack([points, points[:1]])
    return [Segment("wall", closed)]


def _extent(loop: list[Segment]) -> float:
    pts = np.vstack([s.points for s in loop])
    return float((pts.max(axis=0) - pts.min(axis=0)).max())


def _signed_area(loop: list[Segment]) -> float:
    pts = np.vstack([s.points[:-1] for s in loop])
    x, y = pts[:, 0], pts[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


# --- gmsh ------------------------------------------------------------------------------------


@dataclass(frozen=True)
class BoundaryPiece:
    """境界の折れ線の1区間（gmsh の曲線1本ぶん）。`points` は始点から終点の順。"""

    kind: str
    points: np.ndarray


@dataclass(frozen=True)
class Discretization:
    """gmsh の結果：境界の折れ線（閉曲線ごと）と、領域の内部の初期母点。"""

    outer: list[BoundaryPiece]
    holes: list[list[BoundaryPiece]]
    interior_points: np.ndarray

    def rings(self) -> list[np.ndarray]:
        """外周と穴の閉じた折れ線（始点を繰り返さない頂点列）。先頭が外周。"""
        return [_ring(self.outer)] + [_ring(h) for h in self.holes]

    def local_origin(self) -> np.ndarray:
        """局所座標の原点 [m]：外周の外接矩形の左下の点を、1 m 単位に切り下げたもの。

        原点を整数にしておくと、座標の大きさに比べて領域が小さい場合（平面直角座標など）は、
        原点を引く計算と足し戻す計算が丸めなしで済み、境界の節点の座標が元の値に戻る。
        """
        return np.floor(_ring(self.outer).min(axis=0))

    def translated(self, offset: np.ndarray) -> Discretization:
        """境界の折れ線と初期母点を `offset` だけ平行移動したもの。"""
        shift = np.asarray(offset, dtype=float)

        def move(pieces: list[BoundaryPiece]) -> list[BoundaryPiece]:
            return [replace(p, points=p.points + shift) for p in pieces]

        return replace(
            self,
            outer=move(self.outer),
            holes=[move(h) for h in self.holes],
            interior_points=self.interior_points + shift,
        )


def _ring(pieces: list[BoundaryPiece]) -> np.ndarray:
    return np.vstack([p.points[:-1] for p in pieces])


def discretize(domain: Domain, verbose: bool = False) -> Discretization:
    """gmsh で境界を格子幅に合わせて分割し、内部の初期母点を作る。"""
    import gmsh  # 実行時にだけ読み込む（共有ライブラリが必要なため）
    import shapely

    h = domain.cell_size
    gmsh.initialize(readConfigFiles=False, interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay：正三角形に近い分割
        gmsh.option.setNumber("Mesh.MeshSizeMin", h)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.model.add(domain.title)
        geo = gmsh.model.geo

        def add_loop(segments: list[Segment]) -> tuple[int, list[tuple[int, Segment, int]]]:
            tags: list[int] = []
            curves: list[tuple[int, Segment, int]] = []  # (曲線の番号, 区間, 区間内の何本目)
            first = geo.addPoint(*segments[0].start, 0.0, h)
            current = first
            for k, seg in enumerate(segments):
                last_segment = k == len(segments) - 1
                center = geo.addPoint(*seg.arc_center, 0.0, h) if seg.arc_center is not None else 0
                for m in range(1, len(seg.points)):
                    closing = last_segment and m == len(seg.points) - 1
                    nxt = first if closing else geo.addPoint(*seg.points[m], 0.0, h)
                    if seg.arc_center is None:
                        tag = geo.addLine(current, nxt)
                    else:
                        tag = geo.addCircleArc(current, center, nxt)
                    tags.append(tag)
                    curves.append((tag, seg, m - 1))
                    current = nxt
            return geo.addCurveLoop(tags), curves

        outer_loop, outer_curves = add_loop(domain.outer)
        hole_loops = [add_loop(h_) for h_ in domain.holes]
        geo.addPlaneSurface([outer_loop] + [loop for loop, _ in hole_loops])
        geo.synchronize()
        gmsh.model.mesh.generate(1)

        def pieces(curves: list[tuple[int, Segment, int]]) -> list[BoundaryPiece]:
            out = []
            for tag, seg, m in curves:
                start, end = seg.points[m], seg.points[m + 1]
                _, coords, u = gmsh.model.mesh.getNodes(
                    1, tag, includeBoundary=False, returnParametricCoord=True
                )
                inner = np.asarray(coords, dtype=float).reshape(-1, 3)[:, :2]
                inner = inner[np.argsort(np.asarray(u))]
                if len(inner) > 1 and np.hypot(*(inner[0] - start)) > np.hypot(*(inner[-1] - start)):
                    inner = inner[::-1]
                out.append(BoundaryPiece(seg.kind, np.vstack([start, inner, end])))
            return out

        outer = pieces(outer_curves)
        holes = [pieces(curves) for _, curves in hole_loops]

        # 境界から内側に縮めた領域を三角形分割し、その内部の節点を初期の母点にする
        rings = [_ring(outer)] + [_ring(p) for p in holes]
        shrunk = shapely.Polygon(rings[0], holes=rings[1:]).buffer(
            -INTERIOR_OFFSET * h, join_style="mitre", mitre_limit=2.0
        )
        shrunk = shrunk.simplify(0.05 * h)
        gmsh.model.add(f"{domain.title}-interior")
        surfaces = []
        for part in shapely.get_parts(shrunk):
            if not isinstance(part, shapely.Polygon) or part.is_empty:
                continue
            loops = []
            for ring in [part.exterior, *part.interiors]:
                xy = np.asarray(ring.coords)[:-1]
                tags = [geo.addPoint(x, y, 0.0, h) for x, y in xy]
                lines = [geo.addLine(a, b) for a, b in zip(tags, tags[1:] + tags[:1], strict=True)]
                loops.append(geo.addCurveLoop(lines))
            surfaces.append(geo.addPlaneSurface(loops))
        geo.synchronize()
        gmsh.model.mesh.generate(2)
        interior = [np.zeros((0, 2))]
        for surface in surfaces:
            _, coords, _ = gmsh.model.mesh.getNodes(2, surface, includeBoundary=False)
            interior.append(np.asarray(coords, dtype=float).reshape(-1, 3)[:, :2])
        return Discretization(outer=outer, holes=holes, interior_points=np.vstack(interior))
    finally:
        gmsh.finalize()
