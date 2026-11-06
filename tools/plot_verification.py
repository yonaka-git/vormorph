"""検証ケースの計算結果（VTK）から、docs/verification.md に載せる図を描く。

検証テストに VTK を書き出させてから、ケースごとのサブコマンドで描く。

使い方::

    VORMORPH_VERIFICATION_OUTPUT=output/verification \\
        cargo test --release --test verification v06 -- --nocapture
    uv run --project tools python tools/plot_verification.py v06 output/verification \\
        -o docs/images/v06_uniform_flow.png

前提:
    `<dir>/<ケース名>/<ケース名>.pvd` と、そこから参照する `.vtu`（src/io.rs の出力の規約）。
    V04 は、検証テストが書く `<dir>/<ケース名>/centroid.csv`（水塊の重心の時系列）も読む。
    ケースの条件（寸法・流量など）は tests/verification.rs の定数と一致させる。
"""

from __future__ import annotations

import argparse
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import meshio
import numpy as np

from plot_profile import cell_field, longitudinal_profile

# 図の色。系列は検証済みのカテゴリ配色の最初の3色を、この順に使う（3色までなら、どの2色の
# 組み合わせも色覚の多様性を考慮した色差の基準を満たす）。解析解と河床は無彩色にする。
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")
"""系列の色（青、橙、水色）。"""
REFERENCE = "0.65"
"""解析解・理論値の色（太い灰色の線）。"""
BED = "0.25"
"""河床の色。"""

# V06 等流（tests/verification.rs の v06_uniform_flow と同じ条件）
V06_LENGTH = 200.0
"""水路延長 [m]。"""
V06_WIDTH = 20.0
"""水路幅 [m]。"""
V06_SLOPE = 0.002
"""河床勾配 [-]。"""
V06_DISCHARGE = 40.0
"""流量 [m³/s]。"""
V06_MANNING_N = 0.025
"""マニングの粗度係数 [s/m^(1/3)]。"""
V06_CENTRAL = (50.0, 150.0)
"""合格基準を確かめる中央部の範囲 [m]。"""
V06_TOLERANCE = 1.0e-4
"""中央部の水深の相対誤差の合格基準 [-]。"""
V06_CASES = (
    ("v06_rectangular_2m", "rectangular, 2 m"),
    ("v06_hexagonal_2m", "hexagonal, 2 m"),
    ("v06_hexagonal_1m", "hexagonal, 1 m"),
)
"""（ケース名, 凡例）。ケース名は出力のディレクトリ名。"""

# V04 閉鎖水域の振動（tests/verification.rs の thacker モジュールと同じ条件）
V04_A = 1.0
"""容器の半径（静水の汀線の半径）a [m]。"""
V04_H0 = 0.1
"""中心の静水深 h₀ [m]。"""
V04_SHIFT = 0.5
"""水塊の中心の振れ幅 η [m]。"""
V04_L = 4.0
"""計算領域の一辺 [m]（中心は (L/2, L/2)）。"""
V04_G = 9.81
"""重力加速度 [m/s²]。"""
V04_CASES = (
    ("v04_rectangular", "upwind advection, rectangular", "-"),
    ("v04_hexagonal", "upwind advection, hexagonal", (0, (4, 2))),
    ("v04_rectangular_none", "no advection, rectangular", "-"),
)
"""（ケース名, 凡例, 線種）。正方形配置と六角形配置の線は重なるので、六角形配置を破線にする。"""


def thacker_omega() -> float:
    """V04 の角振動数 ω = √(2 g h₀) / a [1/s]。"""
    return float(np.sqrt(2.0 * V04_G * V04_H0) / V04_A)


def thacker_bed(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """V04 の河床高 z_b = −h₀ (1 − r²/a²) [m]。"""
    r2 = (x - 0.5 * V04_L) ** 2 + (y - 0.5 * V04_L) ** 2
    return -V04_H0 * (1.0 - r2 / V04_A**2)


def thacker_depth(x: np.ndarray, y: np.ndarray, t: float) -> np.ndarray:
    """V04 の解析解の水深 h = max(ζ − z_b, 0) [m]（Thacker, 1981）。"""
    wt = thacker_omega() * t
    xs, ys = x - 0.5 * V04_L, y - 0.5 * V04_L
    surface = V04_SHIFT * V04_H0 / V04_A**2 * (2.0 * xs * np.cos(wt) + 2.0 * ys * np.sin(wt) - V04_SHIFT)
    return np.maximum(surface - thacker_bed(x, y), 0.0)


def thacker_centroid(t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """V04 の解析解の水塊の重心 (L/2 + η cos ωt, L/2 + η sin ωt) [m]。"""
    wt = thacker_omega() * np.asarray(t)
    return 0.5 * V04_L + V04_SHIFT * np.cos(wt), 0.5 * V04_L + V04_SHIFT * np.sin(wt)


def read_history(path: Path) -> dict[str, np.ndarray]:
    """検証テストが書いた CSV（1行目が列名）を、列名から配列への辞書として読む。"""
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.size == 0:
        raise ValueError(f"{path}: no rows")
    data = np.atleast_1d(data)
    return {name: np.asarray(data[name], dtype=float) for name in data.dtype.names}


def read_collection(path: Path) -> list[tuple[float, Path]]:
    """`.pvd` を読み、（時刻, `.vtu` のパス）を時刻の順に返す。"""
    root = ET.parse(path).getroot()
    entries = [
        (float(item.attrib["timestep"]), path.parent / item.attrib["file"]) for item in root.iter("DataSet")
    ]
    if not entries:
        raise ValueError(f"{path}: no DataSet entries")
    return sorted(entries, key=lambda entry: entry[0])


def last_step(case_dir: Path) -> tuple[float, meshio.Mesh]:
    """ケースのディレクトリから、最後の時刻の結果を読む。"""
    time, vtu = read_collection(case_dir / f"{case_dir.name}.pvd")[-1]
    return time, meshio.read(vtu)


def normal_depth(discharge: float, width: float, slope: float, manning_n: float) -> float:
    """マニング式の等流水深 h_n = (n q / √I)^(3/5) [m]（広幅矩形水路）。"""
    return float((manning_n * discharge / width / np.sqrt(slope)) ** 0.6)


def plot_v06(root: Path, output: Path) -> None:
    """V06 の図：縦断図と、水深の相対誤差の分布。"""
    h_n = normal_depth(V06_DISCHARGE, V06_WIDTH, V06_SLOPE, V06_MANNING_N)
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True, constrained_layout=True)

    time, mesh = last_step(root / V06_CASES[-1][0])
    x, zb, eta = longitudinal_profile(mesh, y_center=None, half_band=None, bin_width=1.0)
    line = np.linspace(0.0, V06_LENGTH, 3)
    bed_line = V06_SLOPE * (V06_LENGTH - line)
    top.plot(line, bed_line + h_n, color=REFERENCE, lw=3.0, label=f"uniform flow (h_n = {h_n:.4f} m)")
    top.plot(x, eta, color=SERIES[0], lw=1.2, label=f"water surface ({V06_CASES[-1][1]}, t = {time:.0f} s)")
    top.plot(x, zb, color=BED, lw=1.2, label="bed")
    top.set_ylabel("elevation [m]")
    top.legend(loc="upper right", fontsize="small")
    top.grid(alpha=0.3)

    bottom.axvspan(*V06_CENTRAL, color="0.92", label="central part")
    bottom.axhline(V06_TOLERANCE, color="k", lw=0.8, ls="--", label="criterion")
    for (case, label), marker, color in zip(V06_CASES, ("s", "o", "."), SERIES, strict=True):
        _, result = last_step(root / case)
        error = np.abs(cell_field(result, "h") - h_n) / h_n
        bottom.semilogy(
            cell_field(result, "x_cell"),
            np.maximum(error, 1e-16),
            marker,
            color=color,
            ms=2.5,
            ls="none",
            label=label,
        )
    bottom.set_xlabel("x [m]")
    bottom.set_ylabel("|h − h_n| / h_n")
    bottom.set_ylim(1e-13, 1e-2)
    bottom.legend(loc="upper right", fontsize="small", ncol=2)
    bottom.grid(alpha=0.3, which="both")
    fig.savefig(output, dpi=150)


def plot_v04(root: Path, output: Path) -> None:
    """V04 の図：水塊の重心の軌跡、重心の x 座標の時系列、3周期後の水面形。"""
    period = 2.0 * np.pi / thacker_omega()
    fig = plt.figure(figsize=(10, 8), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=(1.0, 1.6))
    orbit, series = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    profile = fig.add_subplot(grid[1, :])

    t = np.linspace(0.0, 3.0 * period, 601)
    ex, ey = thacker_centroid(t)
    orbit.plot(ex, ey, color=REFERENCE, lw=3.0)
    series.plot(t / period, ex - 0.5 * V04_L, color=REFERENCE, lw=3.0)
    for (case, _, style), color in zip(V04_CASES, SERIES, strict=True):
        history = read_history(root / case / "centroid.csv")
        orbit.plot(history["centroid_x"], history["centroid_y"], color=color, ls=style, lw=1.2)
        x = history["centroid_x"] - 0.5 * V04_L
        series.plot(history["time"] / period, x, color=color, ls=style, lw=1.2)
    orbit.set_aspect("equal")
    orbit.set_xlabel("centroid x [m]")
    orbit.set_ylabel("centroid y [m]")
    orbit.set_title("(a) path of the centroid (3 periods)", loc="left", fontsize="medium")
    orbit.grid(alpha=0.3)
    series.set_xlabel("t / T")
    series.set_ylabel("centroid x − L/2 [m]")
    series.set_title("(b) centroid x", loc="left", fontsize="medium")
    series.grid(alpha=0.3)

    y_center = 0.5 * V04_L
    time, _ = last_step(root / V04_CASES[0][0])
    line = np.linspace(0.5 * V04_L - 1.6 * V04_A, 0.5 * V04_L + 1.6 * V04_A, 801)
    exact = thacker_depth(line, np.full_like(line, y_center), time)
    bed = thacker_bed(line, np.full_like(line, y_center))
    profile.plot(line, np.where(exact > 0.0, bed + exact, np.nan), color=REFERENCE, lw=3.0, label="analytic")
    for (case, label, style), color in zip(V04_CASES, SERIES, strict=True):
        _, result = last_step(root / case)
        x, _, eta = longitudinal_profile(result, y_center=y_center, half_band=0.03, bin_width=0.05)
        profile.plot(x, eta, color=color, ls=style, lw=1.2, label=label)
    profile.plot(line, bed, color=BED, lw=1.2, label="bed")
    profile.set_xlim(line[0], line[-1])
    profile.set_ylim(-1.1 * V04_H0, 1.5 * V04_H0)
    profile.set_xlabel(f"x [m] (y = {y_center:g} m)")
    profile.set_ylabel("elevation [m]")
    profile.set_title(f"(c) water surface at t = {time / period:.0f} T", loc="left", fontsize="medium")
    profile.grid(alpha=0.3)
    fig.legend(*profile.get_legend_handles_labels(), loc="outside lower center", ncol=5, fontsize="small")
    fig.savefig(output, dpi=150)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case", choices=["v04", "v06"], help="verification case")
    parser.add_argument("root", type=Path, help="VORMORPH_VERIFICATION_OUTPUT of the test run")
    parser.add_argument("-o", "--output", type=Path, required=True, help="output image")
    args = parser.parse_args(argv)
    plots = {"v04": plot_v04, "v06": plot_v06}
    plots[args.case](args.root, args.output)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
