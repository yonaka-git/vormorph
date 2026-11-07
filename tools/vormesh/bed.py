"""河床の定義（領域の定義の `[bed]`）と、母点での河床高の計算。

vormorph の格子ファイルにはセルの河床高がないので、`tools/make_bed.py` がこの定義から母点の河床高を
求め、セルの値のファイル（`cellsfile`）に書く（ADR-0020）。

```toml
[bed]                                  # 中心線に沿った一様な勾配（湾曲水路など）
kind = "centerline"
slope = 0.002                          # 中心線に沿った勾配（下流へ下がる向きを正）[-]
z_downstream = 0.0                     # 中心線の終点（下流端）の河床高 [m]
centerline = [                         # 上流端から下流端へたどる区間の列（外周と同じ書き方）
  { line = [[30.0, -40.0], [30.0, 0.0]] },
  { arc = { center = [0.0, 0.0], radius = 30.0, start = 0.0, end = 180.0 } },
  { line = [[-30.0, 0.0], [-30.0, -40.0]] },
]

[bed]                                  # 平面
kind = "plane"
z0 = 0.5                               # 原点の河床高 [m]
gradient = [-0.002, 0.0]               # ∂z/∂x, ∂z/∂y [-]
```

中心線の型では、母点を中心線へ射影した点までの、中心線に沿った距離 s を求め、
z_b = z_downstream + slope (s_end − s) とする（s_end は中心線の全長）。射影は区間ごとに最も近い点を
求め、最も近い区間を選ぶ（円弧は厳密に扱う）。中心線の両端より外の母点は、端の点に射影する
（中心線は上流端から下流端までをおおうように与える）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

KINDS = ("centerline", "plane")
"""河床の定義の種類。"""


class BedError(ValueError):
    """河床の定義の誤り。"""


@dataclass(frozen=True)
class LinePiece:
    """中心線の線分。"""

    start: np.ndarray
    end: np.ndarray

    @property
    def length(self) -> float:
        return float(np.hypot(*(self.end - self.start)))

    def project(self, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """各点の最も近い点までの距離と、線分の始点からの距離（線分に沿った長さ）。"""
        d = self.end - self.start
        t = np.clip(((points - self.start) @ d) / float(d @ d), 0.0, 1.0)
        nearest = self.start + t[:, None] * d
        return np.hypot(*(points - nearest).T), t * self.length


@dataclass(frozen=True)
class ArcPiece:
    """中心線の円弧（`start` から `end` へ。角度はラジアン、反時計回りが正）。"""

    center: np.ndarray
    radius: float
    start: float
    end: float

    @property
    def length(self) -> float:
        return self.radius * abs(self.end - self.start)

    def project(self, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """各点の最も近い点までの距離と、円弧の始点からの円弧に沿った長さ。"""
        rel = points - self.center
        sweep = self.end - self.start
        sign = 1.0 if sweep > 0.0 else -1.0
        # 始点から円弧の向きに測った点の角度を [0, 2π) にし、円弧の範囲に収める
        angle = np.mod(sign * (np.arctan2(rel[:, 1], rel[:, 0]) - self.start), 2.0 * math.pi)
        span = abs(sweep)
        # 範囲の外の点は、近い方の端に寄せる
        beyond = angle > span
        to_start = (2.0 * math.pi - angle) < (angle - span)
        angle = np.where(beyond, np.where(to_start, 0.0, span), angle)
        theta = self.start + sign * angle
        nearest = self.center + self.radius * np.column_stack([np.cos(theta), np.sin(theta)])
        return np.hypot(*(points - nearest).T), self.radius * angle


@dataclass(frozen=True)
class CenterlineBed:
    """中心線に沿った一様な勾配の河床。"""

    pieces: list[LinePiece | ArcPiece]
    slope: float
    z_downstream: float

    @property
    def length(self) -> float:
        return float(sum(p.length for p in self.pieces))

    def distance(self, points: np.ndarray) -> np.ndarray:
        """各点を中心線へ射影した点までの、中心線に沿った距離 s [m]（上流端が 0）。"""
        points = np.asarray(points, dtype=float).reshape(-1, 2)
        best = np.full(len(points), np.inf)
        s = np.zeros(len(points))
        offset = 0.0
        for piece in self.pieces:
            gap, along = piece.project(points)
            closer = gap < best
            best = np.where(closer, gap, best)
            s = np.where(closer, offset + along, s)
            offset += piece.length
        return s

    def elevation(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        s = self.distance(np.column_stack([x, y]))
        return self.z_downstream + self.slope * (self.length - s)


@dataclass(frozen=True)
class PlaneBed:
    """平面の河床 z_b = z0 + a x + b y。"""

    z0: float
    gradient: tuple[float, float]

    def elevation(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        a, b = self.gradient
        return self.z0 + a * np.asarray(x, dtype=float) + b * np.asarray(y, dtype=float)


Bed = CenterlineBed | PlaneBed


def parse_bed(spec: dict, source: str = "<bed>") -> Bed:
    """領域の定義の `[bed]` を読む。"""
    kind = spec.get("kind")
    if kind == "centerline":
        _check_keys(spec, {"kind", "slope", "z_downstream", "centerline"}, source)
        try:
            slope, z_downstream = float(spec["slope"]), float(spec["z_downstream"])
            pieces_spec = spec["centerline"]
        except KeyError as err:
            raise BedError(f"{source}: kind = 'centerline' needs {err}") from None
        if not (math.isfinite(slope) and math.isfinite(z_downstream)):
            raise BedError(f"{source}: slope and z_downstream must be finite")
        pieces: list[LinePiece | ArcPiece] = []
        for k, piece in enumerate(pieces_spec):
            pieces.extend(_pieces(piece, f"{source}: centerline[{k}]"))
        if not pieces:
            raise BedError(f"{source}: centerline is empty")
        _check_connected(pieces, source)
        return CenterlineBed(pieces, slope, z_downstream)
    if kind == "plane":
        _check_keys(spec, {"kind", "z0", "gradient"}, source)
        try:
            z0 = float(spec["z0"])
            a, b = (float(g) for g in spec["gradient"])
        except (KeyError, ValueError, TypeError) as err:
            raise BedError(f"{source}: kind = 'plane' needs z0 and gradient = [a, b] ({err})") from None
        if not all(math.isfinite(v) for v in (z0, a, b)):
            raise BedError(f"{source}: z0 and gradient must be finite")
        return PlaneBed(z0, (a, b))
    raise BedError(f"{source}: kind must be one of {KINDS}, got {kind!r}")


def _check_keys(spec: dict, known: set[str], source: str) -> None:
    unknown = sorted(set(spec) - known)
    if unknown:
        raise BedError(f"{source}: unknown keys: {', '.join(unknown)}")


def _pieces(spec: dict, where: str) -> list[LinePiece | ArcPiece]:
    if set(spec) not in ({"line"}, {"arc"}):
        raise BedError(f"{where}: give exactly one of 'line' or 'arc'")
    if "line" in spec:
        points = np.asarray(spec["line"], dtype=float)
        if points.ndim != 2 or points.shape[1] != 2 or len(points) < 2:
            raise BedError(f"{where}: line must be a list of at least two [x, y] points")
        return [LinePiece(a, b) for a, b in zip(points[:-1], points[1:], strict=True)]
    arc = spec["arc"]
    try:
        center = np.asarray(arc["center"], dtype=float)
        radius = float(arc["radius"])
        start, end = math.radians(float(arc["start"])), math.radians(float(arc["end"]))
    except KeyError as err:
        raise BedError(f"{where}: arc needs center, radius, start and end ({err})") from None
    if not radius > 0.0 or start == end:
        raise BedError(f"{where}: arc needs radius > 0 and start != end")
    if abs(end - start) > 2.0 * math.pi:
        raise BedError(f"{where}: an arc of the centerline must not exceed 360 degrees")
    return [ArcPiece(center, radius, start, end)]


def _endpoints(piece: LinePiece | ArcPiece) -> tuple[np.ndarray, np.ndarray]:
    if isinstance(piece, LinePiece):
        return piece.start, piece.end

    def point(angle: float) -> np.ndarray:
        return piece.center + piece.radius * np.array([math.cos(angle), math.sin(angle)])

    return point(piece.start), point(piece.end)


def _check_connected(pieces: list[LinePiece | ArcPiece], source: str) -> None:
    ends = [_endpoints(p) for p in pieces]
    scale = max(1.0, max(float(np.abs(np.vstack(e)).max()) for e in ends))
    for k in range(len(pieces) - 1):
        gap = float(np.hypot(*(ends[k][1] - ends[k + 1][0])))
        if gap > 1.0e-9 * scale:
            raise BedError(
                f"{source}: centerline piece {k} ends at {ends[k][1].tolist()} but piece {k + 1} "
                f"starts at {ends[k + 1][0].tolist()} (gap {gap:.3g})"
            )
