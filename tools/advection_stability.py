"""運動量の移流と時間積分の、1次元の線形安定解析（ADR-0016）。

vormorph の半陰解法（`docs/design/overview.md` の4.1〜4.4節）を、1次元・周期・等間隔の C 格子で
一様流（流速 U > 0、水深 H）のまわりに線形化し、1ステップの増幅行列のスペクトル半径 max|λ| を
求める。辺 j（x = jΔx）に流速 u_j、セル j（辺 j−1 と j の間）に水位 η_j を置く。

    運動量：(1 + Δt γ) u^{n+1} + θ (gΔt/Δx) D η^{n+1} = A u^n − (1 − θ)(gΔt/Δx) D η^n
    連続式：η^{n+1} + (Δt/Δx) D' [θ H u^{n+1}] = η^n − (Δt/Δx) D' [(1 − θ) H u^n + U η_up^n]

- (D η)_j = η_{j+1} − η_j、(D' q)_i = q_i − q_{i−1}。連続式の辺の水深は風上なので、辺 j の流量の
  振れに U η_j が入る（水深を陽的な風上差分で移流しているのと同じ。overview の7章）。
- γ は線形化した摩擦係数 g n² U / H^{4/3}（`--friction`。既定は摩擦なし）。
- A は移流の行列（u^n → F u^n）。流量は一様流の値に固定して線形化する（局所クーラン数 σ_j は
  C = UΔt/Δx）。セル中心の流速は Perot 法で U_j = (u_{j−1} + u_j)/2。
  - upwind（ADR-0009）：F u_j = u_j + C (ū_j − u_j)、ū_j = (U_{j−1} + U_j)/2
  - second_order（ADR-0016）：最小二乗の勾配は G_j = (U_{j+1} − U_{j−1})/2（Δx あたり）で、
    F u_j = u_j + (C/2)(û_{j−1} − û_{j+1} + ũ_{j,j} + ũ_{j+1,j} − 2u_j)。
    面の値 û_l = U_l + G_l/2、辺へ外挿したセルの流速 ũ_{j,j} = U_j + G_j/2、ũ_{j+1,j} = U_{j+1} − G_{j+1}/2

時間積分（S は上の1段。2次の移流の比較のため、採らなかった案も計算できる）：

| 名前 | 内容 |
|---|---|
| `euler` | S を1回（移流は前進 Euler） |
| `ssp_rk2` | 1ステップ全体の SSP-RK2：(I + S²)/2（vormorph の `time_integration = "ssp_rk2"`） |
| `lax_wendroff` | 面の値を Δt/2 だけ上流で評価する（û_l = U_l + (1 − C) G_l/2）。2次の移流だけ |
| `advection_rk2` | 移流だけ SSP-RK2：A = I + ΔtL + (ΔtL)²/2（水位は1回） |
| `ab3` | 移流を3段の Adams–Bashforth で積分する（前の2ステップの移流の項を使う） |

使い方：

    uv run --project tools python tools/advection_stability.py              # 限界クーラン数の表
    uv run --project tools python tools/advection_stability.py --theta 0.6 --friction

限界クーラン数は、C を 1 まで二分法で探した値（その C 以下で max|λ| ≤ 1）。この1次元の解析は、
2次元の V09 で見えた横断方向の蛇行（前進 Euler と2次の移流）は表さないが、同じ組み合わせが
1次元でも不安定になることを示す。
"""

from __future__ import annotations

import argparse

import numpy as np

GRAVITY = 9.81
VELOCITY = 1.87  # 一様流の流速 [m/s]（V06・V09 の等流と同じ）
CELL_SIZE = 2.0  # 格子幅 [m]（V09 の 20 行の正方形配置と同じ）
MANNING_N = 0.025
ADVECTION = ("upwind", "second_order")
TIME_INTEGRATION = ("euler", "ssp_rk2", "lax_wendroff", "advection_rk2", "ab3")
FROUDE = (0.2, 0.3, 0.58, 0.95)
METHODS = (
    ("upwind", "euler"),
    ("upwind", "ssp_rk2"),
    ("second_order", "euler"),
    ("second_order", "lax_wendroff"),
    ("second_order", "advection_rk2"),
    ("second_order", "ab3"),
    ("second_order", "ssp_rk2"),
)


def _shift(u: np.ndarray, k: int) -> np.ndarray:
    """(shift u)_j = u_{j+k}（周期）。"""
    return np.roll(u, -k)


def advected(u: np.ndarray, courant: float, scheme: str, lax_wendroff: bool = False) -> np.ndarray:
    """辺の流速 `u` に移流を前進 Euler で1回施した F u（流量は一様流の値に固定）。"""
    cell = 0.5 * (_shift(u, -1) + u)  # セル j の Perot 流速 U_j
    if scheme == "upwind":
        mean = 0.5 * (_shift(cell, -1) + cell)
        return u + courant * (mean - u)
    if scheme != "second_order":
        raise ValueError(f"unknown advection scheme: {scheme}")
    slope = 0.5 * (_shift(cell, 1) - _shift(cell, -1))  # G_j（Δx あたり）
    back = courant if lax_wendroff else 0.0
    face = cell + 0.5 * (1.0 - back) * slope  # 面 l（セル l の下流の面）の値
    own = cell + 0.5 * slope  # セル j を辺 j まで外挿
    next_ = _shift(cell - 0.5 * slope, 1)  # セル j+1 を辺 j まで外挿
    return u + 0.5 * courant * (_shift(face, -1) - _shift(face, 1) + own + next_ - 2.0 * u)


def advection_matrix(n: int, courant: float, scheme: str, lax_wendroff: bool = False) -> np.ndarray:
    """`advected` を行列にしたもの（n 本の辺）。"""
    columns = [advected(e, courant, scheme, lax_wendroff) for e in np.eye(n)]
    return np.array(columns).T


def _blocks(n: int, depth: float, dt: float, theta: float, gamma: float):
    """半陰解法の1段の左辺 L と、右辺のうち移流によらない部分（状態 (u, η)）。"""
    eye = np.eye(n)
    grad = _shift_matrix(n, 1) - eye  # (D η)_j = η_{j+1} − η_j
    div = eye - _shift_matrix(n, -1)  # (D' q)_i = q_i − q_{i−1}
    a = GRAVITY * dt / CELL_SIZE
    b = dt / CELL_SIZE
    lhs = np.block([[(1.0 + dt * gamma) * eye, theta * a * grad], [b * theta * depth * div, eye]])
    eta_part = np.block(
        [[-(1.0 - theta) * a * grad], [eye - b * VELOCITY * div]]
    )  # η^n の列
    u_part = -b * (1.0 - theta) * depth * div  # 連続式の u^n の項
    return lhs, eta_part, u_part


def _shift_matrix(n: int, k: int) -> np.ndarray:
    """(M u)_j = u_{j+k} の行列。"""
    return np.roll(np.eye(n), k, axis=1)


def depth_for_froude(froude: float) -> float:
    """流速 `VELOCITY` でフルード数が `froude` になる水深 [m]。"""
    return VELOCITY**2 / (GRAVITY * froude**2)


def friction_coefficient(depth: float) -> float:
    """線形化した摩擦係数 γ = g n² U / H^{4/3} [1/s]。"""
    return GRAVITY * MANNING_N**2 * VELOCITY / depth ** (4.0 / 3.0)


def amplification_matrix(
    advection: str,
    time_integration: str,
    froude: float,
    courant: float,
    theta: float = 0.55,
    friction: bool = False,
    n: int = 48,
) -> np.ndarray:
    """1ステップの増幅行列（状態は (u, η)、ab3 では前の2ステップの u も含む）。"""
    if advection not in ADVECTION:
        raise ValueError(f"unknown advection scheme: {advection}")
    if time_integration not in TIME_INTEGRATION:
        raise ValueError(f"unknown time integration: {time_integration}")
    depth = depth_for_froude(froude)
    dt = courant * CELL_SIZE / VELOCITY
    gamma = friction_coefficient(depth) if friction else 0.0
    lhs, eta_part, u_part = _blocks(n, depth, dt, theta, gamma)
    eye = np.eye(n)

    def stage(a: np.ndarray) -> np.ndarray:
        rhs = np.hstack([np.vstack([a, u_part]), eta_part])
        return np.linalg.solve(lhs, rhs)

    euler = advection_matrix(n, courant, advection)
    if time_integration == "euler":
        return stage(euler)
    if time_integration == "ssp_rk2":
        s = stage(euler)
        return 0.5 * (np.eye(2 * n) + s @ s)
    if advection == "upwind":
        raise ValueError(f"{time_integration} is defined only for second_order advection")
    if time_integration == "lax_wendroff":
        return stage(advection_matrix(n, courant, advection, lax_wendroff=True))
    increment = euler - eye  # Δt L
    if time_integration == "advection_rk2":
        return stage(eye + increment + 0.5 * increment @ increment)
    # ab3：F u = u^n + Δt L (23 u^n − 16 u^{n−1} + 5 u^{n−2}) / 12。状態 (u^n, η^n, u^{n−1}, u^{n−2})
    weights = (23.0 / 12.0, -16.0 / 12.0, 5.0 / 12.0)
    size = 4 * n
    big_lhs = np.eye(size)
    big_lhs[: 2 * n, : 2 * n] = lhs
    big_rhs = np.zeros((size, size))
    big_rhs[:n, :n] = eye + weights[0] * increment
    big_rhs[:n, 2 * n : 3 * n] = weights[1] * increment
    big_rhs[:n, 3 * n :] = weights[2] * increment
    big_rhs[n : 2 * n, :n] = u_part
    big_rhs[: 2 * n, n : 2 * n] = eta_part
    big_rhs[2 * n : 3 * n, :n] = eye  # u^{n−1} ← u^n
    big_rhs[3 * n :, 2 * n : 3 * n] = eye  # u^{n−2} ← u^{n−1}
    return np.linalg.solve(big_lhs, big_rhs)


def spectral_radius(*args, **kwargs) -> float:
    """`amplification_matrix` の固有値の絶対値の最大値。"""
    return float(np.abs(np.linalg.eigvals(amplification_matrix(*args, **kwargs))).max())


def critical_courant(
    advection: str,
    time_integration: str,
    froude: float,
    theta: float = 0.55,
    friction: bool = False,
    tol: float = 1e-9,
) -> float:
    """max|λ| ≤ 1 + `tol` となるクーラン数の上限（1 まで。二分法）。"""

    def stable(courant: float) -> bool:
        radius = spectral_radius(advection, time_integration, froude, courant, theta, friction)
        return radius <= 1.0 + tol

    if stable(1.0):
        return 1.0
    low, high = 0.0, 1.0
    for _ in range(20):
        mid = 0.5 * (low + high)
        if stable(mid):
            low = mid
        else:
            high = mid
    return low


def main() -> None:
    parser = argparse.ArgumentParser(description="移流と時間積分の1次元の線形安定解析（ADR-0016）")
    parser.add_argument("--theta", type=float, default=0.55, help="θ 法の重み（既定 0.55）")
    parser.add_argument("--friction", action="store_true", help="マニング式の摩擦を入れる（n = 0.025）")
    args = parser.parse_args()
    label = "摩擦あり" if args.friction else "摩擦なし"
    print(f"限界クーラン数（θ = {args.theta}、{label}、U = {VELOCITY} m/s、Δx = {CELL_SIZE} m）")
    print("| 移流 | 時間積分 | " + " | ".join(f"Fr {fr}" for fr in FROUDE) + " |")
    print("|---|---|" + "---|" * len(FROUDE))
    for advection, time_integration in METHODS:
        values = [
            critical_courant(advection, time_integration, fr, args.theta, args.friction)
            for fr in FROUDE
        ]
        cells = " | ".join("≥ 1" if v >= 1.0 else f"{v:.2f}" for v in values)
        print(f"| {advection} | {time_integration} | {cells} |")


if __name__ == "__main__":
    main()
