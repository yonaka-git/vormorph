"""格子の図：セルをヒルベルト曲線の順（セル番号）で塗り、母点を番号順に結ぶ。"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import LinearSegmentedColormap

from .voronoi import PolygonMesh

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
SEQUENTIAL = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7"]
"""セル番号の塗り分け（単一色相の連続的なランプ）。番号順の線が読めるよう、暗い段は使わない。"""
BOUNDARY_STYLE = {
    "wall": {"color": TEXT_SECONDARY, "label": "wall"},
    "upstream": {"color": "#eb6834", "label": "upstream"},
    "downstream": {"color": "#1baf7a", "label": "downstream"},
}


def plot_mesh(mesh: PolygonMesh, path: Path, title: str = "", show_path: bool = True) -> None:
    """格子を PNG などに描く。"""
    polys = [mesh.cell_polygon(i) for i in range(mesh.n_cells)]
    xy = np.column_stack([mesh.x_vertex, mesh.y_vertex])
    lower, upper = xy.min(axis=0), xy.max(axis=0)
    span = upper - lower
    width = 10.0
    height = float(np.clip(width * span[1] / span[0], 2.5, 12.0)) + 1.2

    fig, ax = plt.subplots(figsize=(width, height), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    cmap = LinearSegmentedColormap.from_list("cell_index", SEQUENTIAL)
    cells = PolyCollection(polys, array=np.arange(mesh.n_cells), cmap=cmap, edgecolor=SURFACE, linewidth=0.4)
    ax.add_collection(cells)

    if show_path:
        ax.plot(mesh.x_cell, mesh.y_cell, color=TEXT, lw=0.6, alpha=0.85, label="cell order (Hilbert curve)")
    for kind, style in BOUNDARY_STYLE.items():
        pairs = mesh.boundary.get(kind, [])
        if not pairs:
            continue
        segs = [[xy[a], xy[b]] for a, b in pairs]
        ax.add_collection(
            LineCollection(
                segs, colors=style["color"], linewidths=2.0, capstyle="round", label=style["label"]
            )
        )

    pad = 0.03 * float(span.max())
    ax.set_xlim(lower[0] - pad, upper[0] + pad)
    ax.set_ylim(lower[1] - pad, upper[1] + pad)
    ax.set_aspect("equal")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(TEXT_SECONDARY)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)
    ax.set_xlabel("x [m]", color=TEXT_SECONDARY, fontsize=9)
    ax.set_ylabel("y [m]", color=TEXT_SECONDARY, fontsize=9)
    if title:
        ax.set_title(title, color=TEXT, fontsize=11, loc="left")
    legend = ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=8)
    for text in legend.get_texts():
        text.set_color(TEXT)
    bar = fig.colorbar(cells, ax=ax, fraction=0.025, pad=0.02)
    bar.set_label("cell index", color=TEXT_SECONDARY, fontsize=8)
    bar.ax.tick_params(colors=TEXT_SECONDARY, labelsize=7)
    bar.outline.set_visible(False)
    fig.savefig(path, dpi=150, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
