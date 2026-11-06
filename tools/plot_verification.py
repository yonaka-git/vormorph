"""検証ケースの計算結果（VTK）から、docs/verification.md に載せる図を描く。

検証テストに VTK を書き出させてから、ケースごとのサブコマンドで描く。

使い方::

    VORMORPH_VERIFICATION_OUTPUT=output/verification \\
        cargo test --release --test verification v06 -- --nocapture
    uv run --project tools python tools/plot_verification.py v06 output/verification \\
        -o docs/images/v06_uniform_flow.png

前提:
    `<dir>/<ケース名>/<ケース名>.pvd` と、そこから参照する `.vtu`（src/io.rs の出力の規約）。
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
    top.plot(line, bed_line + h_n, color="0.6", lw=3.0, label=f"uniform flow (h_n = {h_n:.4f} m)")
    top.plot(x, eta, color="C0", lw=1.2, label=f"water surface ({V06_CASES[-1][1]}, t = {time:.0f} s)")
    top.plot(x, zb, color="C3", lw=1.2, label="bed")
    top.set_ylabel("elevation [m]")
    top.legend(loc="upper right", fontsize="small")
    top.grid(alpha=0.3)

    bottom.axvspan(*V06_CENTRAL, color="0.92", label="central part")
    bottom.axhline(V06_TOLERANCE, color="k", lw=0.8, ls="--", label="criterion")
    for (case, label), marker in zip(V06_CASES, ("s", "o", "."), strict=True):
        _, result = last_step(root / case)
        error = np.abs(cell_field(result, "h") - h_n) / h_n
        bottom.semilogy(
            cell_field(result, "x_cell"), np.maximum(error, 1e-16), marker, ms=2.5, ls="none", label=label
        )
    bottom.set_xlabel("x [m]")
    bottom.set_ylabel("|h − h_n| / h_n")
    bottom.set_ylim(1e-13, 1e-2)
    bottom.legend(loc="upper right", fontsize="small", ncol=2)
    bottom.grid(alpha=0.3, which="both")
    fig.savefig(output, dpi=150)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case", choices=["v06"], help="verification case")
    parser.add_argument("root", type=Path, help="VORMORPH_VERIFICATION_OUTPUT of the test run")
    parser.add_argument("-o", "--output", type=Path, required=True, help="output image")
    args = parser.parse_args(argv)
    if args.case == "v06":
        plot_v06(args.root, args.output)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
