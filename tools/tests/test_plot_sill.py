"""plot_sill.py のテスト。"""

from __future__ import annotations

from pathlib import Path

import meshio
import numpy as np
import pytest

from plot_sill import bed_change_colormap, main, sill_extent, symmetric_limit


def write_case(directory: Path) -> Path:
    """4 × 2 個の正方形セルの水路（x = 0〜4、y = 0〜2）の2時刻を書き、.pvd を返す。

    x = 1〜2 の列を床止め（fixed_bed = 1）とし、最後の時刻に下流側を掘る。
    """
    nx, ny = 4, 2
    points = np.array([[i, j, 0.0] for j in range(ny + 1) for i in range(nx + 1)], dtype=float)
    cells = []
    for j in range(ny):
        for i in range(nx):
            a = j * (nx + 1) + i
            cells.append([a, a + 1, a + nx + 2, a + nx + 1])
    xc = np.array([i + 0.5 for j in range(ny) for i in range(nx)])
    yc = np.array([j + 0.5 for j in range(ny) for i in range(nx)])
    zb0 = 0.01 * (nx - xc)
    entries = []
    for k, dz in enumerate([np.zeros_like(xc), np.where(xc > 2.0, -0.2, 0.0)]):
        zb = zb0 + dz
        data = {
            "x_cell": xc,
            "y_cell": yc,
            "zb": zb,
            "eta": zb0 + 1.0,
            "h": zb0 + 1.0 - zb,
            "dzb": dz,
            "fixed_bed": np.where((xc > 1.0) & (xc < 2.0), 1.0, 0.0),
        }
        name = f"case_{k:06d}.vtu"
        meshio.write(
            directory / name,
            meshio.Mesh(points, [("quad", np.array(cells))], cell_data={key: [v] for key, v in data.items()}),
        )
        entries.append((60.0 * k, name))
    pvd = directory / "case.pvd"
    body = "\n".join(f'    <DataSet timestep="{t}" group="" part="0" file="{n}"/>' for t, n in entries)
    pvd.write_text(
        f'<?xml version="1.0"?>\n<VTKFile type="Collection"><Collection>\n{body}\n</Collection></VTKFile>\n',
        encoding="utf-8",
    )
    return pvd


def test_symmetric_limit_centers_the_colors_on_zero() -> None:
    assert symmetric_limit(np.array([-0.3, 0.1, 0.2])) == pytest.approx(0.3)
    assert symmetric_limit(np.zeros(3)) == 1.0


def test_colormap_runs_from_scour_through_neutral_to_deposition() -> None:
    cmap = bed_change_colormap()
    red, middle, blue = cmap(0.0), cmap(0.5), cmap(1.0)
    assert red[0] > red[2] and blue[2] > blue[0]
    assert max(middle[:3]) - min(middle[:3]) < 0.03


def test_sill_extent_and_figure(tmp_path: Path) -> None:
    pvd = write_case(tmp_path)
    assert sill_extent(meshio.read(tmp_path / "case_000001.vtu")) == (1.0, 2.0)
    output = tmp_path / "hero.png"
    assert main([str(pvd), "-o", str(output), "--x-range", "0", "4"]) == 0
    assert output.stat().st_size > 10_000
