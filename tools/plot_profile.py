"""VTK（.vtu）の計算結果から縦断図（河床高と水位）を描く。

水路の中心線の近くにあるセルを切り出し、x 方向のビンごとに平均して描く。

使い方::

    uv run --project tools python tools/plot_profile.py output/sill/sill_*.vtu -o profile.png

前提（src/io.rs の出力の規約）:
    セルデータに x_cell, y_cell, zb, eta がある。h があれば、乾いたセルの水位は描かない。
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import meshio  # noqa: E402
import numpy as np  # noqa: E402

DRY_DEPTH = 1.0e-3
"""これより浅いセルは乾いているとみなし、水位を描かない [m]。"""


def cell_field(mesh: meshio.Mesh, name: str) -> np.ndarray:
    """セルデータ `name` を、全セルブロックをつないだ1次元配列として返す。"""
    if name not in mesh.cell_data:
        available = ", ".join(sorted(mesh.cell_data)) or "(none)"
        raise KeyError(f"cell data '{name}' not found (available: {available})")
    return np.concatenate([np.asarray(block).reshape(len(block), -1)[:, 0] for block in mesh.cell_data[name]])


def longitudinal_profile(
    mesh: meshio.Mesh,
    y_center: float | None,
    half_band: float | None,
    bin_width: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """中心線の近くのセルを x 方向のビンで平均し、(x, 河床高, 水位) を返す。"""
    x = cell_field(mesh, "x_cell")
    y = cell_field(mesh, "y_cell")
    zb = cell_field(mesh, "zb")
    eta = cell_field(mesh, "eta")
    if "h" in mesh.cell_data:
        eta = np.where(cell_field(mesh, "h") > DRY_DEPTH, eta, np.nan)

    if y_center is None:
        y_center = 0.5 * (y.min() + y.max())
    if half_band is None:
        half_band = 0.05 * (y.max() - y.min()) + 0.5 * bin_width
    selected = np.abs(y - y_center) <= half_band
    if not selected.any():
        raise ValueError(f"no cells within {half_band} m of y = {y_center} m")

    x, zb, eta = x[selected], zb[selected], eta[selected]
    edges = np.arange(x.min(), x.max() + bin_width, bin_width)
    index = np.clip(np.digitize(x, edges) - 1, 0, len(edges) - 2)
    centers, zb_mean, eta_mean = [], [], []
    for k in np.unique(index):
        in_bin = index == k
        centers.append(x[in_bin].mean())
        zb_mean.append(zb[in_bin].mean())
        wet = eta[in_bin][~np.isnan(eta[in_bin])]
        eta_mean.append(wet.mean() if wet.size else np.nan)
    return np.array(centers), np.array(zb_mean), np.array(eta_mean)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("files", nargs="+", type=Path, help=".vtu files (one line per file)")
    parser.add_argument("-o", "--output", type=Path, default=Path("profile.png"), help="output image")
    parser.add_argument("--y-center", type=float, default=None, help="centerline y [m] (default: mid-width)")
    parser.add_argument("--band", type=float, default=None, help="half width of the sampling band [m]")
    parser.add_argument("--bin", type=float, default=1.0, help="bin width along x [m]")
    args = parser.parse_args(argv)

    fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.0, 0.9, len(args.files)))
    for path, color in zip(sorted(args.files), colors, strict=True):
        x, zb, eta = longitudinal_profile(meshio.read(path), args.y_center, args.band, args.bin)
        ax.plot(x, zb, color=color, lw=1.6, label=f"bed ({path.stem})")
        ax.plot(x, eta, color=color, lw=1.0, ls="--", label=f"water surface ({path.stem})")

    ax.set_xlabel("x [m]")
    ax.set_ylabel("elevation [m]")
    ax.grid(alpha=0.3)
    if len(args.files) <= 4:
        ax.legend(loc="best", fontsize="small")
    fig.savefig(args.output, dpi=150)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
