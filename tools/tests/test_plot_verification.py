"""plot_verification.py の補助関数のテスト。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from plot_verification import (
    V04_A,
    V04_H0,
    V04_L,
    V05_G,
    V05_HL,
    V05_L,
    V05_T,
    V05_X0,
    V07_G,
    V07_Q,
    bar_beta,
    bar_growth_curve,
    bar_refinement_pairs,
    bump_depth,
    bump_jump_solution,
    bump_linear_celerity,
    crest_position,
    dam_break_depth,
    normal_depth,
    read_bar_summary,
    read_collection,
    read_history,
    stoker_cm,
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


@pytest.mark.parametrize("hr", [0.0, 0.001])
def test_dam_break_solutions_conserve_the_water_volume(hr: float) -> None:
    # 壁に届く前なので、水の体積は初期の h_l x₀ + h_r (L − x₀) のまま
    x = np.linspace(0.0, V05_L, 200_001)
    h = dam_break_depth(x, V05_T, hr)
    volume = np.trapezoid(h, x)
    assert volume == pytest.approx(V05_HL * V05_X0 + hr * (V05_L - V05_X0), rel=1e-4)
    assert h[0] == V05_HL and h[-1] == hr


def test_stoker_state_satisfies_the_shock_relations() -> None:
    # 中間の状態と下流の状態は、衝撃波の速さ s で質量と運動量の Rankine–Hugoniot の関係を満たす
    hr = 0.001
    cm = stoker_cm(hr)
    hm, um = cm * cm / V05_G, 2.0 * (np.sqrt(V05_G * V05_HL) - cm)
    s = hm * um / (hm - hr)
    momentum = (hm * um * um + V05_G * hm * hm / 2.0) - V05_G * hr * hr / 2.0
    assert momentum == pytest.approx(s * hm * um, rel=1e-9)
    # SWASHES の値（h_m ≈ 2.54 mm、x_C(6 s) ≈ 6.26 m）
    assert hm == pytest.approx(2.539e-3, abs=1e-6)
    assert V05_X0 + V05_T * s == pytest.approx(6.2598, abs=1e-3)


def test_bump_jump_satisfies_belanger_and_matches_swashes() -> None:
    x_s, h1, h2, upstream, downstream = bump_jump_solution()
    froude = V07_Q / (h1 * np.sqrt(V07_G * h1))
    assert h2 / h1 == pytest.approx((np.sqrt(1.0 + 8.0 * froude**2) - 1.0) / 2.0, rel=1e-9)
    # SWASHES の図の跳水の位置（約 11.7 m）と、上流の水深（約 0.41 m）
    assert x_s == pytest.approx(11.6656, abs=1e-3)
    assert bump_depth(np.array([2.0]))[0] == pytest.approx(0.4137, abs=1e-3)
    # 跳水で全水頭が失われ、下流の水深は境界条件の 0.33 m
    assert upstream > downstream
    assert bump_depth(np.array([20.0]))[0] == pytest.approx(0.33, abs=1e-9)


def test_bump_linear_celerity_matches_the_rust_value() -> None:
    # tests/verification.rs の V11 が表示する線形理論の値（0.1184 m/s）
    assert bump_linear_celerity() == pytest.approx(0.1184, abs=5e-5)


@pytest.mark.parametrize("peak", [10.3, 10.5, 10.74])
def test_crest_position_is_exact_for_a_parabola(peak: float) -> None:
    x = np.arange(0.0, 20.0, 0.5)
    b = 1.0 - (x - peak) ** 2
    assert crest_position(x, b) == pytest.approx(peak, abs=1e-12)


def test_read_bar_summary_and_refinement_pairs(tmp_path: Path) -> None:
    summary = tmp_path / "v09_summary.csv"
    summary.write_text(
        "case,layout,rows,width,wavelength,advection,cells,steps,growth,celerity,theory_growth,theory_celerity\n"
        "a,Rectangular,10,40,211.4,Upwind,530,1459,1.0e-3,0.3,2.0e-3,0.34\n"
        "b,Rectangular,20,40,211.4,Upwind,2120,2740,1.5e-3,0.33,2.0e-3,0.34\n"
        "c,Hexagonal,10,10,142.5,Upwind,1230,3852,-1.0e-2,0.28,-5.0e-3,0.31\n"
        "d,Rectangular,5,10,142.5,Upwind,355,1666,-6.0e-3,0.31,-5.0e-3,0.31\n"
        "e,Rectangular,10,10,142.5,Upwind,1430,3354,-5.5e-3,0.31,-5.0e-3,0.31\n"
        "f,Rectangular,20,40,211.4,SecondOrder,2120,1714,1.8e-3,0.34,2.0e-3,0.34\n"
        "g,Rectangular,40,40,211.4,SecondOrder,8480,3400,1.9e-3,0.34,2.0e-3,0.34\n",
        encoding="utf-8",
    )
    rows = read_bar_summary(summary)
    assert rows[0]["case"] == "a" and rows[0]["layout"] == "Rectangular"
    assert rows[5]["advection"] == "SecondOrder"
    assert rows[1]["growth"] == pytest.approx(1.5e-3)
    pairs = bar_refinement_pairs(rows)
    # 六角形の幅 10 m は1つしかないので組にならない。移流の違う 20 行どうしも組にしない
    assert [(layout, coarse["case"], fine["case"]) for layout, coarse, fine in pairs] == [
        ("Rectangular", "a", "b"),
        ("Rectangular", "d", "e"),
        ("Rectangular", "f", "g"),
    ]


def test_bar_growth_curve_matches_the_v09_reference_value() -> None:
    # tests/verification.rs の v09_alternate_bars（MORFAC 100）の理論値
    assert bar_growth_curve(40.0, np.array([211.4]))[0] == pytest.approx(2.5895e-3, rel=2e-4)
    assert bar_beta(40.0) == pytest.approx(18.705, abs=1e-3)
