"""plot_verification.py の補助関数のテスト。"""

from __future__ import annotations

from pathlib import Path

import pytest

from plot_verification import normal_depth, read_collection


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
