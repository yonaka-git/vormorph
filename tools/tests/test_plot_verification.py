"""plot_verification.py の補助関数のテスト。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from plot_verification import (
    V04_A,
    V04_H0,
    V04_L,
    normal_depth,
    read_collection,
    read_history,
    thacker_centroid,
    thacker_depth,
    thacker_omega,
)


def test_read_collection_sorts_entries_by_time(tmp_path: Path) -> None:
    pvd = tmp_path / "case.pvd"
    pvd.write_text(
        """<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
    <DataSet timestep="600" group="" part="0" file="case_000002.vtu"/>
    <DataSet timestep="0" group="" part="0" file="case_000000.vtu"/>
    <DataSet timestep="300" group="" part="0" file="case_000001.vtu"/>
  </Collection>
</VTKFile>
""",
        encoding="utf-8",
    )
    entries = read_collection(pvd)
    assert [t for t, _ in entries] == [0.0, 300.0, 600.0]
    assert entries[-1][1] == tmp_path / "case_000002.vtu"


def test_read_collection_rejects_an_empty_collection(tmp_path: Path) -> None:
    pvd = tmp_path / "empty.pvd"
    pvd.write_text('<VTKFile type="Collection"><Collection/></VTKFile>', encoding="utf-8")
    with pytest.raises(ValueError, match="no DataSet"):
        read_collection(pvd)


def test_normal_depth_matches_the_sill_example() -> None:
    # examples/sill.toml のコメントの目安（h ≈ 1.07 m）
    assert normal_depth(40.0, 20.0, 0.002, 0.025) == pytest.approx(1.069235, abs=1e-6)


@pytest.mark.parametrize("fraction", [0.0, 0.3, 0.75])
def test_thacker_centroid_and_volume_follow_the_depth_field(fraction: float) -> None:
    # 解析解の水深を細かい格子で積分すると、重心は式の値に一致し、体積は時刻によらない
    t = fraction * 2.0 * np.pi / thacker_omega()
    n = 801
    x, y = np.meshgrid(np.linspace(0.0, V04_L, n), np.linspace(0.0, V04_L, n))
    h = thacker_depth(x, y, t)
    volume = h.sum()
    cx, cy = (h * x).sum() / volume, (h * y).sum() / volume
    ex, ey = thacker_centroid(t)
    assert cx == pytest.approx(ex, abs=1e-4)
    assert cy == pytest.approx(ey, abs=1e-4)
    # 放物面の容器の静水の体積 π a² h₀ / 2
    cell = (V04_L / (n - 1)) ** 2
    assert volume * cell == pytest.approx(np.pi * V04_A**2 * V04_H0 / 2.0, rel=1e-3)


def test_read_history_reads_named_columns(tmp_path: Path) -> None:
    path = tmp_path / "centroid.csv"
    path.write_text("time,centroid_x\n0.5,2.25\n1.0,2.5\n", encoding="utf-8")
    history = read_history(path)
    assert history["time"].tolist() == [0.5, 1.0]
    assert history["centroid_x"].tolist() == [2.25, 2.5]
    one_row = tmp_path / "one.csv"
    one_row.write_text("time,centroid_x\n0.5,2.25\n", encoding="utf-8")
    assert read_history(one_row)["time"].tolist() == [0.5]
