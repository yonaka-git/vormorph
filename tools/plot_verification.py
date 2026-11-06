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
    V05 は、`<dir>/v05_convergence.csv`（格子幅ごとの L1 誤差）も読む。
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


# V05 ダム破壊（tests/verification.rs の dam_break モジュールと同じ条件）
V05_L = 10.0
"""水路の延長 [m]。"""
V05_X0 = 5.0
"""ダムの位置 [m]。"""
V05_HL = 0.005
"""上流側の水深 [m]。"""
V05_T = 6.0
"""比べる時刻 [s]。"""
V05_G = 9.81
"""重力加速度 [m/s²]。"""
V05_CASES = (("stoker", 0.001, "Stoker (h_r = 1 mm)"), ("ritter", 0.0, "Ritter (dry bed)"))
"""（ケース名, 下流側の水深 [m], 題名）。"""
V05_PROFILES = (
    ("rectangular", 0.1, "rectangular, 0.1 m", "-"),
    ("rectangular", 0.0125, "rectangular, 0.0125 m", "-"),
    ("hexagonal", 0.0125, "hexagonal, 0.0125 m", (0, (4, 2))),
)
"""水深の分布を描く計算（配置, 格子幅 [m], 凡例, 線種）。"""


def stoker_cm(hr: float) -> float:
    """Stoker の解の中間の状態の波速 c_m = √(g h_m) [m/s]（二分法）。"""
    c0 = np.sqrt(V05_G * V05_HL)

    def f(c: float) -> float:
        return -8.0 * V05_G * hr * c * c * (c0 - c) ** 2 + (c * c - V05_G * hr) ** 2 * (c * c + V05_G * hr)

    lo, hi = np.sqrt(V05_G * hr), c0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) < 0.0 else (lo, mid)
    return float(0.5 * (lo + hi))


def dam_break_depth(x: np.ndarray, t: float, hr: float) -> np.ndarray:
    """ダム破壊の解析解の水深 [m]。`hr` = 0 なら Ritter、`hr` > 0 なら Stoker の解。"""
    c0 = np.sqrt(V05_G * V05_HL)
    xi = (np.asarray(x, dtype=float) - V05_X0) / t
    rarefaction = 4.0 / (9.0 * V05_G) * (c0 - 0.5 * xi) ** 2
    if hr == 0.0:
        return np.where(xi <= -c0, V05_HL, np.where(xi <= 2.0 * c0, rarefaction, 0.0))
    cm = stoker_cm(hr)
    shock = 2.0 * cm * cm * (c0 - cm) / (cm * cm - V05_G * hr)
    return np.select(
        [xi <= -c0, xi <= 2.0 * c0 - 3.0 * cm, xi <= shock],
        [V05_HL, rarefaction, cm * cm / V05_G],
        default=hr,
    )


def binned_depth(mesh: meshio.Mesh, bin_width: float) -> tuple[np.ndarray, np.ndarray]:
    """セルの水深を x 方向のビンごとに平均し、(x, h) を返す（乾いたセルも 0 として含める）。"""
    x, h = cell_field(mesh, "x_cell"), cell_field(mesh, "h")
    index = np.floor(x / bin_width).astype(int)
    bins = np.unique(index)
    return (
        np.array([x[index == k].mean() for k in bins]),
        np.array([h[index == k].mean() for k in bins]),
    )


# V07 こぶを越える跳水（tests/verification.rs の bump_jump モジュールと同じ条件）
V07_Q = 0.18
"""単位幅流量 [m²/s]。"""
V07_H_DOWNSTREAM = 0.33
"""下流端の水深 [m]。"""
V07_Z_MAX = 0.2
"""こぶの頂部の高さ [m]。"""
V07_WIDTH = 0.2
"""水路の幅 [m]。"""
V07_G = 9.81
"""重力加速度 [m/s²]。"""
V07_CASES = (
    ("v07_rectangular_0.1m", 0.1, "rectangular, 0.1 m", "-"),
    ("v07_rectangular_0.05m", 0.05, "rectangular, 0.05 m", "-"),
    ("v07_hexagonal_0.05m", 0.05, "hexagonal, 0.05 m", (0, (4, 2))),
)
"""（ケース名, 格子幅 [m], 凡例, 線種）。"""


def bump_bed(x: np.ndarray) -> np.ndarray:
    """V07 の河床高 z = 0.2 − 0.05 (x − 10)²（8 m < x < 12 m）[m]。"""
    x = np.asarray(x, dtype=float)
    return np.where((x > 8.0) & (x < 12.0), V07_Z_MAX - 0.05 * (x - 10.0) ** 2, 0.0)


def bump_depth_for_head(head: float, z: float, supercritical: bool) -> float:
    """全水頭 `head` の流れの、河床高 `z` での射流側または常流側の水深 [m]（なければ NaN）。"""
    hc = (V07_Q**2 / V07_G) ** (1.0 / 3.0)
    energy = head - z

    def specific_energy(h: float) -> float:
        return h + V07_Q**2 / (2.0 * V07_G * h * h)

    if energy < specific_energy(hc):
        return float("nan")
    lo, hi = (1e-6, hc) if supercritical else (hc, energy)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (specific_energy(mid) > energy) == supercritical:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def bump_jump_solution() -> tuple[float, float, float, float, float]:
    """V07 の解析解の（跳水の位置 [m], 共役水深 h₁, h₂ [m], 上流側・下流側の全水頭 [m]）。"""
    hc = (V07_Q**2 / V07_G) ** (1.0 / 3.0)
    upstream = V07_Z_MAX + 1.5 * hc
    downstream = V07_H_DOWNSTREAM + V07_Q**2 / (2.0 * V07_G * V07_H_DOWNSTREAM**2)

    def balance(x: float) -> float:
        z = float(bump_bed(x))
        h1 = bump_depth_for_head(upstream, z, True)
        h2 = bump_depth_for_head(downstream, z, False)
        return V07_Q**2 * (1.0 / h1 - 1.0 / h2) + 0.5 * V07_G * (h1 * h1 - h2 * h2)

    lo, hi = 10.0, 12.0
    while np.isnan(balance(lo)):
        lo += 1e-3
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if balance(mid) > 0.0 else (lo, mid)
    x_s = 0.5 * (lo + hi)
    z = float(bump_bed(x_s))
    h1 = bump_depth_for_head(upstream, z, True)
    h2 = bump_depth_for_head(downstream, z, False)
    return x_s, h1, h2, upstream, downstream


def bump_depth(x: np.ndarray) -> np.ndarray:
    """V07 の解析解の水深 [m]。"""
    x_s, _, _, upstream, downstream = bump_jump_solution()
    return np.array(
        [
            bump_depth_for_head(upstream, float(bump_bed(xi)), xi > 10.0)
            if xi < x_s
            else bump_depth_for_head(downstream, float(bump_bed(xi)), False)
            for xi in np.atleast_1d(x)
        ]
    )


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


def plot_v05(root: Path, output: Path) -> None:
    """V05 の図：Stoker と Ritter の水深の分布（t = 6 s）と、L1 誤差の収束。"""
    fig = plt.figure(figsize=(10, 8), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)
    profiles = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    convergence = fig.add_subplot(grid[1, :])

    line = np.linspace(0.0, V05_L, 2001)
    for ax, (case, hr, title), letter in zip(profiles, V05_CASES, "ab", strict=True):
        ax.plot(line, 1e3 * dam_break_depth(line, V05_T, hr), color=REFERENCE, lw=3.0, label="analytic")
        for (layout, dx, label, style), color in zip(V05_PROFILES, SERIES, strict=True):
            time, mesh = last_step(root / f"v05_{case}_{layout}_{dx}m")
            x, h = binned_depth(mesh, dx)
            ax.plot(x, 1e3 * h, color=color, ls=style, lw=1.2, label=label)
        ax.set_xlim(2.5, 8.5)
        ax.set_ylim(-0.1, 5.4)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("depth [mm]")
        ax.set_title(f"({letter}) {title}, t = {time:g} s", loc="left", fontsize="medium")
        ax.grid(alpha=0.3)
    profiles[0].legend(loc="lower left", fontsize="small")

    table = np.genfromtxt(
        root / "v05_convergence.csv", delimiter=",", names=True, dtype=None, encoding="utf-8"
    )
    for (case, _, title), color in zip(V05_CASES, SERIES[:2], strict=True):
        for layout, style, marker in (("rectangular", "-", "s"), ("hexagonal", (0, (4, 2)), "o")):
            rows = table[(table["case"] == case) & (table["layout"] == layout)]
            convergence.loglog(
                rows["cell_size"],
                100.0 * rows["l1_error"],
                color=color,
                ls=style,
                marker=marker,
                ms=6,
                lw=1.2,
                label=f"{title}, {layout}",
            )
    sizes = np.array([0.1, 0.0125])
    for order, style in ((1.0, ":"), (0.5, "--")):
        convergence.loglog(
            sizes, 0.25 * (sizes / sizes[-1]) ** order, color=BED, ls=style, lw=0.8, label=f"order {order:g}"
        )
    convergence.set_xlabel("cell size [m]")
    convergence.set_ylabel("L1 error of depth [%]")
    convergence.set_title(
        "(c) convergence (L1 error relative to the water volume)", loc="left", fontsize="medium"
    )
    convergence.legend(loc="upper left", fontsize="small", ncol=2)
    convergence.grid(alpha=0.3, which="both")
    fig.savefig(output, dpi=150)


def plot_v07(root: Path, output: Path) -> None:
    """V07 の図：水面の縦断形、跳水の近くの水深、横断方向の流速。"""
    x_s, h1, h2, _, _ = bump_jump_solution()
    fig, (top, middle, bottom) = plt.subplots(3, 1, figsize=(10, 9.5), constrained_layout=True)
    line = np.linspace(0.0, 25.0, 2501)
    exact = bump_depth(line)
    bed = bump_bed(line)
    zoom = (line > 10.8) & (line < 12.6)
    top.plot(line, bed + exact, color=REFERENCE, lw=3.0, label="analytic")
    middle.plot(line[zoom], exact[zoom], color=REFERENCE, lw=3.0, label="analytic")
    for (case, dx, label, style), color in zip(V07_CASES, SERIES, strict=True):
        time, mesh = last_step(root / case)
        x, _, eta = longitudinal_profile(mesh, y_center=0.5 * V07_WIDTH, half_band=V07_WIDTH, bin_width=dx)
        top.plot(x, eta, color=color, ls=style, lw=1.2, label=label)
        x, h = binned_depth(mesh, dx)
        near = (x > 10.8) & (x < 12.6)
        middle.plot(x[near], h[near], color=color, ls=style, lw=1.2, marker="o", ms=3, label=label)
        blocks = mesh.cell_data["velocity"]
        v = np.abs(np.concatenate([np.asarray(b).reshape(len(b), -1) for b in blocks])[:, 1])
        xc = cell_field(mesh, "x_cell")
        index = np.floor(xc / dx).astype(int)
        bins = np.unique(index)
        bottom.semilogy(
            [xc[index == k].mean() for k in bins],
            np.maximum([v[index == k].max() for k in bins], 1e-16),
            color=color,
            ls=style,
            lw=1.2,
            label=label,
        )
    top.plot(line, bed, color=BED, lw=1.2, label="bed")
    top.set_xlim(5.0, 16.0)
    top.set_ylim(-0.02, 0.48)
    top.set_xlabel("x [m]")
    top.set_ylabel("elevation [m]")
    top.set_title(f"(a) water surface, t = {time:g} s", loc="left", fontsize="medium")
    top.legend(loc="upper right", fontsize="small", ncol=2)
    top.grid(alpha=0.3)

    middle.plot([x_s, x_s], [h1, h2], color=BED, lw=0.8, ls=":")
    for depth, name, offset in ((h1, "h₁", (8, -4)), (h2, "h₂", (-80, 0))):
        middle.annotate(
            f"{name} = {depth:.3f} m",
            (x_s, depth),
            xytext=offset,
            textcoords="offset points",
            fontsize="small",
            va="center",
        )
    middle.set_ylim(0.06, 0.345)
    middle.set_xlabel("x [m]")
    middle.set_ylabel("depth [m]")
    middle.set_title(
        f"(b) depth near the jump (analytic jump at x = {x_s:.2f} m, h₂/h₁ = {h2 / h1:.2f})",
        loc="left",
        fontsize="medium",
    )
    middle.legend(loc="lower right", fontsize="small")
    middle.grid(alpha=0.3)

    bottom.set_xlim(5.0, 16.0)
    bottom.set_ylim(1e-16, 1.0)
    bottom.set_xlabel("x [m]")
    bottom.set_ylabel("max |v| in the cross-section [m/s]")
    bottom.set_title("(c) transverse velocity", loc="left", fontsize="medium")
    bottom.legend(loc="center right", fontsize="small")
    bottom.grid(alpha=0.3)
    fig.savefig(output, dpi=150)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case", choices=["v04", "v05", "v06", "v07"], help="verification case")
    parser.add_argument("root", type=Path, help="VORMORPH_VERIFICATION_OUTPUT of the test run")
    parser.add_argument("-o", "--output", type=Path, required=True, help="output image")
    args = parser.parse_args(argv)
    plots = {"v04": plot_v04, "v05": plot_v05, "v06": plot_v06, "v07": plot_v07}
    plots[args.case](args.root, args.output)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
