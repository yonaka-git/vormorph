"""交互砂州の線形安定解析（vormorph と同じ支配方程式）。

検証ケース V09 の比較対象を計算する（`docs/verification.md`）。支配方程式は vormorph と同じにする
（`docs/design/overview.md`）。

- 水深平均の浅水流方程式（移流あり、水平渦粘性なし）、マニング式の摩擦 g n² |u| u / h^{4/3}
- 掃流砂：MPM 式または芦田・道上式。長谷川の式による全方向の勾配補正
  q_b = q̂(τ*) [u/|u| − γ ∇z_b]、γ = √(τ*_c / (μ_s μ_k τ*))、τ* = n² |u|² / (h^{1/3} s d)
- Exner 式 (1 − λ) ∂z_b/∂t + M ∇·q_b = 0（MORFAC M は河床の式だけに掛ける）
- 二次流（任意。ADR-0019）：強さ I の輸送方程式（遠心加速度を生成項とし、適応長 L_I で平衡値へ
  緩和する）か平衡型。掃流砂の向きを tan δ = N* I / |u| だけ内岸側へ回す（`Closure.secondary_flow`）
- 側壁は滑り条件。勾配 I の等流（水深 H、流速 U）のまわりで線形化し、横断方向のモード m の
  u, h, z_b ∝ cos(mπy/B) e^{ikx}、v ∝ sin(mπy/B) e^{ikx} の振幅 X = (u, v, h, z) について
  dX/dt = J X を立てる（`jacobian`）。二次流の輸送方程式を解くときは、I ∝ sin(mπy/B) e^{ikx} を
  5番目の状態に加える。

J の固有値のうち、MORFAC が小さい極限で 0 に近づく根が河床のモードで、実部が成長率、
−Im ω / k が移動速度になる（`bed_mode`）。どちらも流れの時間あたりの値（MORFAC を掛けた後）。
MORFAC が大きいと、流れの調整の遅れのために準定常の理論からずれる（B = 40 m、MORFAC 100 で
成長率が +14%）ので、vormorph と比べるときは同じ MORFAC で計算する。

使い方：

    uv run --project tools python tools/bar_theory.py --width 40 --morfac 100 --wavelength 211.4
    uv run --project tools python tools/bar_theory.py --width 40 --scan        # 最も速く育つ波長
    uv run --project tools python tools/bar_theory.py --critical-width         # 発生の境の幅
    uv run --project tools python tools/bar_theory.py --width 40 --morfac 100 --wavelength 211.4 \
        --secondary-flow transport                                               # 二次流あり

参考：線形安定解析の考え方は Colombini, Seminara & Tubino (1987) などの交互砂州の理論と同じで、
Federici & Seminara (2003) は、その不安定が移流型であることを示した。
"""

from __future__ import annotations

import argparse
import itertools
from dataclasses import dataclass, replace

import numpy as np

GRAVITY = 9.81


@dataclass(frozen=True)
class Closure:
    """流れと流砂の条件（SI 単位）。既定値は `examples/sill.toml` と同じ。"""

    unit_discharge: float = 2.0  # 単位幅流量 q [m²/s]
    slope: float = 0.002  # 河床勾配 I [-]
    manning_n: float = 0.025  # マニングの粗度係数
    grain_size: float = 0.002  # 粒径 d [m]
    relative_density: float = 1.65  # 水中比重 s = ρ_s/ρ − 1
    critical_shields: float = 0.047  # τ*_c
    static_friction: float = 1.0  # μ_s
    kinetic_friction: float = 0.45  # μ_k
    porosity: float = 0.4  # λ
    morfac: float = 1.0  # M
    gravity: float = GRAVITY
    formula: str = "mpm"  # "mpm" または "ashida_michiue"
    advection: bool = True  # 運動量の移流を含めるか（vormorph の advection = "none" に対応）
    secondary_flow: str = "none"  # 二次流："none"、"transport"（輸送方程式型）、"equilibrium"（平衡型）
    n_star: float = 7.0  # 掃流砂の偏角の係数 N*
    diffusivity_coef: float = 0.15  # 二次流の強さの拡散係数 D = α_D u* h の α_D
    von_karman: float = 0.41  # κ


@dataclass(frozen=True)
class Base:
    """等流の基本状態。"""

    depth: float  # H [m]
    velocity: float  # U [m/s]
    shields: float  # τ*0
    bedload: float  # q̂0 = √(s g d³) Φ(τ*0) [m²/s]
    phi_t: float  # d ln Φ / d ln τ*（τ*0 での値）
    gamma: float  # 長谷川の γ0


def intensity(formula: str, tau: float, tau_c: float) -> float:
    """無次元掃流砂量 Φ(τ*)（τ* ≤ τ*_c では 0）。"""
    if tau <= tau_c:
        return 0.0
    if formula == "mpm":
        return 8.0 * (tau - tau_c) ** 1.5
    if formula == "ashida_michiue":
        r = tau_c / tau
        return 17.0 * tau**1.5 * (1.0 - r) * (1.0 - np.sqrt(r))
    raise ValueError(f"unknown formula {formula!r}")


def shields(c: Closure, speed: float, depth: float) -> float:
    """無次元掃流力 τ* = n² |u|² / (h^{1/3} s d)。"""
    return c.manning_n**2 * speed**2 / (depth ** (1.0 / 3.0) * c.relative_density * c.grain_size)


def base_state(c: Closure) -> Base:
    """等流の水深 H = (n q / √I)^{3/5} と、その状態の掃流砂の諸量。"""
    h = (c.manning_n * c.unit_discharge / np.sqrt(c.slope)) ** 0.6
    u = c.unit_discharge / h
    tau = shields(c, u, h)
    phi = intensity(c.formula, tau, c.critical_shields)
    if phi <= 0.0:
        raise ValueError(f"no bedload in the uniform flow (tau* = {tau:.4f})")
    eps = 1e-6
    phi_t = (
        np.log(intensity(c.formula, tau * (1 + eps), c.critical_shields))
        - np.log(intensity(c.formula, tau * (1 - eps), c.critical_shields))
    ) / (2 * eps)
    scale = np.sqrt(c.relative_density * c.gravity * c.grain_size**3)
    gamma = np.sqrt(c.critical_shields / (c.static_friction * c.kinetic_friction * tau))
    return Base(h, u, tau, scale * phi, phi_t, gamma)


def adaptation_length(c: Closure, depth: float) -> float:
    """二次流の適応長 L_I = (1 − 2α) h / (2κ² α)、α = √g n / (κ h^{1/6})（α ≥ 1/2 なら 0）。"""
    alpha = np.sqrt(c.gravity) * c.manning_n / (c.von_karman * depth ** (1.0 / 6.0))
    return 0.0 if alpha >= 0.5 else (1.0 - 2.0 * alpha) * depth / (2.0 * c.von_karman**2 * alpha)


def secondary_flow_diffusivity(c: Closure, speed: float, depth: float) -> float:
    """二次流の強さの拡散係数 D = α_D u* h、u* = √g n |u| / h^{1/6} [m²/s]。"""
    return c.diffusivity_coef * np.sqrt(c.gravity) * c.manning_n * speed / depth ** (1.0 / 6.0) * depth


def jacobian(c: Closure, k: float, width: float, mode: int = 1) -> np.ndarray:
    """状態 X = (u, v, h, z[, I]) の線形系 dX/dt = J X の係数行列（4×4 か 5×5、複素数）。

    z は勾配を除いた河床高の振れ。x 方向の掃流砂量は q̂ (1 + γI − γ ∂z/∂x) なので、τ* の振れ
    τ*'/τ*0 = 2u/U − h/(3H) と、γ ∝ τ*^{-1/2} の振れの両方を含める。連続式は
    ∂h/∂t = −∂z/∂t − ∇·(h u)（河床の変化の分だけ水深が変わる。水位を保つ vormorph の河床の
    更新と同じ）。

    二次流（`secondary_flow`。vormorph の `[sediment.secondary_flow]`、ADR-0019）：等流では I = 0
    で、振れの平衡値は I_e = h a_n / |u| = H ∂v/∂x（a_n ≈ U ∂v/∂x）。横断方向の掃流砂量に
    q̂ N* I / U が加わる（tan δ = N* I / |u|。向きの正規化 √(1 + tan²δ) は2次の微小量）。

    - `"transport"`：I ∝ sin(μy) を5番目の状態にし、
      ∂I/∂t = −U ∂I/∂x + D ∇²I + (U / L_I)(I_e − I)（hI の式から連続式を引いたもの）。
      側壁で ∂I/∂y = 0 の条件は満たさない（D μ² の項は側壁の境界層を無視した近似）。
    - `"equilibrium"`：I = I_e を河床の式に代入する（4×4 のまま）。
    """
    b = base_state(c)
    g, i, h, u = c.gravity, c.slope, b.depth, b.velocity
    mu = mode * np.pi / width
    a = b.phi_t * (1.0 + b.gamma * i) - 0.5 * b.gamma * i
    coef = -c.morfac * b.bedload / (1.0 - c.porosity)
    bed = np.array(
        [
            coef * 1j * k * a * 2.0 / u,
            coef * mu / u,
            -coef * 1j * k * a / (3.0 * h),
            coef * b.gamma * (k * k + mu * mu),
        ],
        dtype=complex,
    )
    # 横断方向の掃流砂量の二次流の項 q̂ N* I / U の、河床の式への寄与（I の係数）
    spiral = coef * mu * c.n_star / u
    if c.secondary_flow == "equilibrium":
        bed[1] += spiral * h * 1j * k
    elif c.secondary_flow not in ("none", "transport"):
        raise ValueError(f"unknown secondary flow {c.secondary_flow!r}")
    size = 5 if c.secondary_flow == "transport" else 4
    adv = 1j * k * u if c.advection else 0.0
    j = np.zeros((size, size), dtype=complex)
    j[0, :4] = [-adv - 2 * g * i / u, 0, -1j * k * g + 4 / 3 * g * i / h, -1j * k * g]
    j[1, :4] = [0, -adv - g * i / u, g * mu, g * mu]
    j[3, :4] = bed
    j[2, :4] = -bed + np.array([-1j * k * h, -h * mu, -1j * k * u, 0])
    if size == 5:
        length = adaptation_length(c, h)
        rate = u / length
        diffusion = secondary_flow_diffusivity(c, u, h)
        j[3, 4] = spiral
        j[2, 4] = -spiral
        j[4, 1] = rate * h * 1j * k
        j[4, 4] = -1j * k * u - diffusion * (k * k + mu * mu) - rate
    return j


def bed_mode(c: Closure, k: float, width: float, mode: int = 1) -> complex:
    """河床のモードの固有値 ω \\[1/s\\]（流れの時間あたり）。

    MORFAC が小さい極限では、河床のモードは M に比例して 0 に近づく。M を 10⁻⁶ 倍から少しずつ
    上げながら、前の段の値を M の比で延ばした予測に最も近い根を追う。
    """
    factors = np.geomspace(1e-6, 1.0, 49)
    ev = np.linalg.eigvals(jacobian(replace(c, morfac=factors[0] * c.morfac), k, width, mode))
    omega = ev[np.argmin(abs(ev))]
    for previous, f in itertools.pairwise(factors):
        ev = np.linalg.eigvals(jacobian(replace(c, morfac=f * c.morfac), k, width, mode))
        omega = ev[np.argmin(abs(ev - omega * f / previous))]
    return complex(omega)


def quasi_steady_mode(c: Closure, k: float, width: float, mode: int = 1) -> complex:
    """流れを準定常とした河床のモードの ω（M → 0 の極限。検算用）。"""
    j = jacobian(replace(c, morfac=1.0), k, width, mode)
    # 流れの式（z 以外）で ∂/∂t を省き、z = 1 に対する (u, v, h[, I]) を解いて、河床の式に代入する。
    # 水深の式の ∂h/∂t = −∂z/∂t − ∇·(hu) の −∂z/∂t（河床の式）を、h の行に戻してから解く
    flow = [r for r in range(j.shape[0]) if r != 3]
    a = j[np.ix_(flow, flow)].copy()
    rhs = -j[flow, 3].copy()
    a[2] += j[3, flow]
    rhs[2] -= j[3, 3]
    x = np.linalg.solve(a, rhs)
    return complex(c.morfac * (j[3, flow] @ x + j[3, 3]))


def growth_and_celerity(c: Closure, wavelength: float, width: float, mode: int = 1) -> tuple[float, float]:
    """成長率 [1/s] と移動速度 [m/s]（流れの時間あたり）。"""
    k = 2 * np.pi / wavelength
    omega = bed_mode(c, k, width, mode)
    return omega.real, -omega.imag / k


def most_unstable(c: Closure, width: float, wavelengths: np.ndarray) -> tuple[float, float]:
    """`wavelengths` のうち最も速く育つ波長 [m] と、その成長率 [1/s]。"""
    growth = np.array([growth_and_celerity(c, wl, width)[0] for wl in wavelengths])
    best = int(np.argmax(growth))
    return float(wavelengths[best]), float(growth[best])


def critical_width(c: Closure, wavelengths: np.ndarray, low: float = 1.0, high: float = 200.0) -> float:
    """交互砂州が発生する（どれかの波長で成長率が正になる）最小の水路幅 [m]（二分法）。"""

    def unstable(width: float) -> bool:
        return most_unstable(c, width, wavelengths)[1] > 0.0

    if unstable(low) or not unstable(high):
        raise ValueError("the critical width is not between the bounds")
    for _ in range(50):
        mid = 0.5 * (low + high)
        if unstable(mid):
            high = mid
        else:
            low = mid
    return high


def main() -> None:
    parser = argparse.ArgumentParser(description="交互砂州の線形安定解析（vormorph と同じ式）")
    parser.add_argument("--q", type=float, default=2.0, help="単位幅流量 [m²/s]")
    parser.add_argument("--slope", type=float, default=0.002, help="河床勾配")
    parser.add_argument("--n", type=float, default=0.025, help="マニングの粗度係数")
    parser.add_argument("--d", type=float, default=0.002, help="粒径 [m]")
    parser.add_argument("--formula", default="mpm", choices=["mpm", "ashida_michiue"])
    parser.add_argument("--width", type=float, default=40.0, help="水路幅 [m]")
    parser.add_argument("--morfac", type=float, default=1.0)
    parser.add_argument("--wavelength", type=float, nargs="*", default=[], help="波長 [m]")
    parser.add_argument("--scan", action="store_true", help="最も速く育つ波長を探す")
    parser.add_argument("--critical-width", action="store_true", help="発生の境の幅を求める")
    parser.add_argument(
        "--secondary-flow", default="none", choices=["none", "transport", "equilibrium"], help="二次流"
    )
    parser.add_argument("--n-star", type=float, default=7.0, help="掃流砂の偏角の係数 N*")
    parser.add_argument("--diffusivity", type=float, default=0.15, help="二次流の拡散係数の係数 α_D")
    args = parser.parse_args()

    c = Closure(
        args.q,
        args.slope,
        args.n,
        args.d,
        formula=args.formula,
        morfac=args.morfac,
        secondary_flow=args.secondary_flow,
        n_star=args.n_star,
        diffusivity_coef=args.diffusivity,
    )
    b = base_state(c)
    froude = b.velocity / np.sqrt(c.gravity * b.depth)
    print(
        f"H = {b.depth:.6f} m, U = {b.velocity:.6f} m/s, Fr = {froude:.4f}, tau* = {b.shields:.4f}, "
        f"q_b0 = {b.bedload:.4e} m²/s, Phi_T = {b.phi_t:.4f}, gamma = {b.gamma:.4f}"
    )
    print(f"B = {args.width} m, beta = B/(2H) = {args.width / (2 * b.depth):.3f}, MORFAC {args.morfac}")
    for wl in args.wavelength:
        s, v = growth_and_celerity(c, wl, args.width)
        print(f"L = {wl:9.3f} m: growth {s:+.6e} 1/s, celerity {v:+.6e} m/s")
    wavelengths = np.geomspace(2.0 * b.depth, 400.0 * b.depth, 400)
    if args.scan:
        wl, s = most_unstable(c, args.width, wavelengths)
        print(f"most unstable: L = {wl:.2f} m ({wl / args.width:.2f} B), growth {s:+.6e} 1/s")
    if args.critical_width:
        w = critical_width(c, wavelengths)
        print(f"critical width: B = {w:.3f} m, beta_c = {w / (2 * b.depth):.3f}")


if __name__ == "__main__":
    main()
