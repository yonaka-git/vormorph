"""床止めの例題（examples/sill.toml）の計算結果から、河床変動量の平面図と縦断図を描く。

README の冒頭の図（docs/images/hero.png）を、この図で作る。

使い方::

    cargo run --release -- run examples/sill.toml
    uv run --project tools python tools/plot_sill.py output/sill/sill.pvd -o docs/images/hero.png

前提（src/io.rs の出力の規約）:
    `.pvd` の最初の時刻が初期状態、最後の時刻が描く状態。セルデータに x_cell, y_cell, zb, eta, h,
    dzb がある。床止めの位置（非侵食層の帯）は fixed_bed から求める。
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import meshio
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap

from plot_profile import cell_field, longitudinal_profile
from plot_verification import read_collection

EROSION = "#e34948"
"""洗掘（河床が下がる）の色（発散配色の赤の極）。"""
DEPOSITION = "#2a78d6"
"""堆積（河床が上がる）の色（発散配色の青の極）。"""
NEUTRAL = "#f0efec"
"""変化のない所の色（発散配色の中点の灰色）。"""
WATER = "#2a78d6"
"""水面の線の色。"""
INITIAL_BED = "0.65"
"""初期の河床の線の色（無彩色）。"""
BED = "0.2"
"""計算後の河床の線の色。"""


def bed_change_colormap() -> LinearSegmentedColormap:
    """河床変動量の発散配色（洗掘 = 赤、0 = 灰色、堆積 = 青）。0 を中央に置いて使う。"""
    return LinearSegmentedColormap.from_list("bed_change", [EROSION, NEUTRAL, DEPOSITION])


def symmetric_limit(values: np.ndarray) -> float:
    """0 を中央に置く色の範囲の上限（絶対値の最大値。すべて 0 なら 1）。"""
    limit = float(np.nanmax(np.abs(values))) if values.size else 0.0
    return limit if limit > 0.0 else 1.0


def cell_polygons(mesh: meshio.Mesh) -> list[np.ndarray]:
    """セルの多角形（頂点の xy 座標）を、セルデータと同じ順に返す。"""
    points = mesh.points[:, :2]
    return [points[cell] for block in mesh.cells for cell in block.data]


def sill_extent(mesh: meshio.Mesh) -> tuple[float, float] | None:
    """非侵食層のあるセル（fixed_bed = 1）の x の範囲（なければ None）。"""
    if "fixed_bed" not in mesh.cell_data:
        return None
    fixed = cell_field(mesh, "fixed_bed") > 0.5
    if not fixed.any():
        return None
    polygons = cell_polygons(mesh)
    xs = np.concatenate([polygons[c][:, 0] for c in np.flatnonzero(fixed)])
    return float(xs.min()), float(xs.max())


def plot_sill(pvd: Path, output: Path, x_range: tuple[float, float]) -> None:
    """平面図（河床変動量）と、中心線の縦断図（初期の河床、計算後の河床、水面）を描く。"""
    entries = read_collection(pvd)
    initial = meshio.read(entries[0][1])
    time, path = entries[-1]
    final = meshio.read(path)

    x = cell_field(final, "x_cell")
    y = cell_field(final, "y_cell")
    dzb = cell_field(final, "dzb")
    polygons = cell_polygons(final)
    inside = (x >= x_range[0]) & (x <= x_range[1])
    limit = symmetric_limit(dzb[inside])

    width = float(y.max() + y.min())
    fig, (plan, profile) = plt.subplots(
        2, 1, figsize=(11, 6.2), constrained_layout=True, height_ratios=(1.0, 1.15)
    )
    cells = PolyCollection(
        [polygons[c] for c in np.flatnonzero(inside)],
        array=dzb[inside],
        cmap=bed_change_colormap(),
        edgecolors="face",
        linewidths=0.2,
    )
    cells.set_clim(-limit, limit)
    plan.add_collection(cells)
    plan.set_xlim(*x_range)
    plan.set_ylim(0.0, width)
    plan.set_aspect("equal")
    plan.set_xlabel("x [m]")
    plan.set_ylabel("y [m]")
    colorbar = fig.colorbar(cells, ax=plan, shrink=0.9, pad=0.01)
    colorbar.set_label("bed change [m]\n(− scour, + deposition)")
    sill = sill_extent(final)
    if sill is not None:
        plan.text(
            0.5 * (sill[0] + sill[1]),
            0.93 * width,
            "sill",
            ha="center",
            va="top",
            fontsize="small",
            color="0.2",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
        )
    plan.set_title(
        f"(a) bed change after t = {time:g} s (plan view)",
        loc="left",
        fontsize="medium",
    )

    center = 0.5 * width
    x0, zb0, _ = longitudinal_profile(initial, y_center=center, half_band=1.0, bin_width=0.5)
    x1, zb1, eta1 = longitudinal_profile(final, y_center=center, half_band=1.0, bin_width=0.5)
    keep0 = (x0 >= x_range[0]) & (x0 <= x_range[1])
    keep1 = (x1 >= x_range[0]) & (x1 <= x_range[1])
    zb0_on_1 = np.interp(x1, x0, zb0)
    profile.fill_between(
        x1[keep1],
        zb0_on_1[keep1],
        zb1[keep1],
        where=zb1[keep1] < zb0_on_1[keep1],
        color=EROSION,
        alpha=0.35,
        lw=0.0,
        label="scour",
    )
    profile.fill_between(
        x1[keep1],
        zb0_on_1[keep1],
        zb1[keep1],
        where=zb1[keep1] > zb0_on_1[keep1],
        color=DEPOSITION,
        alpha=0.35,
        lw=0.0,
        label="deposition",
    )
    profile.plot(x0[keep0], zb0[keep0], color=INITIAL_BED, lw=1.5, label="initial bed")
    profile.plot(x1[keep1], zb1[keep1], color=BED, lw=1.5, label=f"bed at t = {time:g} s")
    profile.plot(x1[keep1], eta1[keep1], color=WATER, lw=1.2, label="water surface")
    profile.set_xlim(*x_range)
    profile.set_xlabel("x [m]")
    profile.set_ylabel("elevation [m]")
    profile.set_title("(b) longitudinal profile along the centerline", loc="left", fontsize="medium")
    profile.legend(loc="upper right", fontsize="small", ncol=3)
    profile.grid(alpha=0.3)
    fig.savefig(output, dpi=150)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("pvd", type=Path, help="the .pvd written by `vormorph run examples/sill.toml`")
    parser.add_argument("-o", "--output", type=Path, required=True, help="output image")
    parser.add_argument(
        "--x-range", type=float, nargs=2, default=(50.0, 150.0), metavar=("X0", "X1"), help="x range [m]"
    )
    args = parser.parse_args(argv)
    plot_sill(args.pvd, args.output, tuple(args.x_range))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
