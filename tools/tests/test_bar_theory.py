"""交互砂州の線形安定解析（tools/bar_theory.py）のテスト。

手で導いた係数行列 J を、非線形の支配方程式に小さな振れを与えて数値的に線形化したものと
照合する。V09（tests/verification.rs）が比べる理論の値も、ここで確かめる。
"""

from dataclasses import replace

import numpy as np
import pytest

import bar_theory as bt

SILL = bt.Closure()  # examples/sill.toml の流れと流砂の条件


def test_uniform_flow_matches_the_v06_depth_and_shields_stress():
    b = bt.base_state(SILL)
    assert b.depth == pytest.approx(1.069235, abs=1e-6)
    assert b.velocity == pytest.approx(2.0 / 1.069235, rel=1e-6)
    assert b.shields == pytest.approx(0.648, abs=1e-3)
    # MPM：d ln Φ / d ln τ* = 1.5 τ* / (τ* − τ*_c)
    assert b.phi_t == pytest.approx(1.5 * b.shields / (b.shields - 0.047), rel=1e-8)
    assert b.gamma == pytest.approx(np.sqrt(0.047 / (0.45 * b.shields)), rel=1e-12)


def _tendencies(c: bt.Closure, fields, x: float, y: float) -> np.ndarray:
    """非線形の支配方程式の時間微分 (∂u/∂t, ∂v/∂t, ∂h/∂t, ∂z/∂t[, ∂I/∂t]) を点 (x, y) で評価する。

    `fields(x, y)` は (u, v, h, z_b, I) を返す（z_b は勾配を含む実際の河床高）。空間微分は中心差分で
    求める（場は滑らかなので、刻み 1e-3 m の打ち切り誤差は十分小さい）。二次流は vormorph と同じ
    式（ADR-0019）：掃流砂の向きは河床付近の流速 u + tan δ (k × u)（tan δ = N* I / |u|）、
    平衡型なら I = h a_n / |u|、a_n = (u × (u·∇)u)_z / |u|。
    """
    d = 1e-3
    g, n2 = c.gravity, c.manning_n**2

    def ddx(f, x, y):
        return (f(x + d, y) - f(x - d, y)) / (2 * d)

    def ddy(f, x, y):
        return (f(x, y + d) - f(x, y - d)) / (2 * d)

    def component(i):
        return lambda x, y: fields(x, y)[i]

    u_f, v_f, h_f, z_f, i_f = (component(i) for i in range(5))

    def eta(x, y):
        return h_f(x, y) + z_f(x, y)

    def equilibrium_intensity(x, y):
        u, v, h, _, _ = fields(x, y)
        ax = u * ddx(u_f, x, y) + v * ddy(u_f, x, y)
        ay = u * ddx(v_f, x, y) + v * ddy(v_f, x, y)
        return h * (u * ay - v * ax) / (u * u + v * v)

    def intensity(x, y):
        if c.secondary_flow == "equilibrium":
            return equilibrium_intensity(x, y)
        if c.secondary_flow == "transport":
            return i_f(x, y)
        return 0.0

    def bedload(x, y):
        u, v, h, _, _ = fields(x, y)
        speed = np.hypot(u, v)
        tau = bt.shields(c, speed, h)
        scale = np.sqrt(c.relative_density * c.gravity * c.grain_size**3)
        q = scale * bt.intensity(c.formula, tau, c.critical_shields)
        gamma = np.sqrt(c.critical_shields / (c.static_friction * c.kinetic_friction * tau))
        tan = c.n_star * intensity(x, y) / speed
        # 河床付近の流速 u + tan δ (k × u) の向き
        bx, by = (u - tan * v), (v + tan * u)
        norm = np.hypot(bx, by)
        return q * (bx / norm - gamma * ddx(z_f, x, y)), q * (by / norm - gamma * ddy(z_f, x, y))

    u, v, h, _, spiral = fields(x, y)
    speed = np.hypot(u, v)
    friction = g * n2 * speed / h ** (4.0 / 3.0)
    adv = 1.0 if c.advection else 0.0
    du = -adv * (u * ddx(u_f, x, y) + v * ddy(u_f, x, y)) - g * ddx(eta, x, y) - friction * u
    dv = -adv * (u * ddx(v_f, x, y) + v * ddy(v_f, x, y)) - g * ddy(eta, x, y) - friction * v
    div_qb = ddx(lambda x, y: bedload(x, y)[0], x, y) + ddy(lambda x, y: bedload(x, y)[1], x, y)
    dz = -c.morfac * div_qb / (1.0 - c.porosity)
    div_q = ddx(lambda x, y: h_f(x, y) * u_f(x, y), x, y) + ddy(lambda x, y: h_f(x, y) * v_f(x, y), x, y)
    dh = -dz - div_q
    if c.secondary_flow != "transport":
        return np.array([du, dv, dh, dz])
    # ∂(hI)/∂t + ∇·(qI) = ∇·(hD∇I) + h|u|(I_e − I)/L_I から連続式を引いた I の式
    length = bt.adaptation_length(c, h)

    def hd_ddx(x, y):
        uu, vv, hh, _, _ = fields(x, y)
        return hh * bt.secondary_flow_diffusivity(c, np.hypot(uu, vv), hh) * ddx(i_f, x, y)

    def hd_ddy(x, y):
        uu, vv, hh, _, _ = fields(x, y)
        return hh * bt.secondary_flow_diffusivity(c, np.hypot(uu, vv), hh) * ddy(i_f, x, y)

    diffusion = (ddx(hd_ddx, x, y) + ddy(hd_ddy, x, y)) / h
    di = (
        -(u * ddx(i_f, x, y) + v * ddy(i_f, x, y))
        + diffusion
        + speed * (equilibrium_intensity(x, y) - spiral) / length
    )
    return np.array([du, dv, dh, dz, di])


@pytest.mark.parametrize(
    "closure",
    [
        replace(SILL, morfac=100.0),
        replace(SILL, morfac=3.0, formula="ashida_michiue", critical_shields=0.05),
        replace(SILL, morfac=10.0, advection=False),
        replace(SILL, morfac=100.0, secondary_flow="transport"),
        replace(SILL, morfac=100.0, secondary_flow="equilibrium", n_star=11.0),
    ],
    ids=["mpm", "ashida_michiue", "no_advection", "secondary_flow", "equilibrium_secondary_flow"],
)
def test_jacobian_matches_a_numerical_linearization(closure):
    width, wavelength, mode = 40.0, 150.0, 1
    k, mu = 2 * np.pi / wavelength, mode * np.pi / width
    b = bt.base_state(closure)
    j = bt.jacobian(closure, k, width, mode)
    size = j.shape[0]
    rng = np.random.default_rng(1)
    amplitude = rng.normal(size=size) + 1j * rng.normal(size=size)
    eps = 1e-6

    def shapes_at(y):
        # u, h, z ∝ cos(μy)、v と I ∝ sin(μy)
        return np.array([np.cos(mu * y), np.sin(mu * y), np.cos(mu * y), np.cos(mu * y), np.sin(mu * y)])[:size]

    def fields(x, y, scale):
        p = scale * (amplitude * np.exp(1j * k * x) * shapes_at(y)).real
        spiral = p[4] if size == 5 else 0.0
        return (b.velocity + p[0], p[1], b.depth + p[2], -closure.slope * x + p[3], spiral)

    for x, y in [(13.0, 7.0), (61.0, 29.0)]:
        base = _tendencies(closure, lambda x, y: fields(x, y, 0.0), x, y)
        perturbed = _tendencies(closure, lambda x, y: fields(x, y, eps), x, y)
        numerical = (perturbed - base) / eps
        shapes = shapes_at(y)
        analytic = ((j @ amplitude) * np.exp(1j * k * x) * shapes).real
        # 等流はつり合っている（流速と河床の時間微分は 0）
        assert np.allclose(base, 0.0, atol=1e-10)
        scale = np.abs(j @ amplitude) * np.abs(shapes)
        assert np.allclose(numerical, analytic, rtol=0.0, atol=2e-4 * scale + 1e-12), (
            numerical,
            analytic,
        )


def test_bed_mode_approaches_the_quasi_steady_limit():
    k, width = 2 * np.pi / 211.4, 40.0
    small = replace(SILL, morfac=1e-3)
    omega = bt.bed_mode(small, k, width)
    limit = bt.quasi_steady_mode(small, k, width)
    assert omega == pytest.approx(limit, rel=1e-4)
    # MORFAC 100 では、流れの調整の遅れで成長率が準定常の理論より 14% 大きい
    large = replace(SILL, morfac=100.0)
    ratio = bt.bed_mode(large, k, width).real / bt.quasi_steady_mode(large, k, width).real
    assert ratio == pytest.approx(1.1437, abs=1e-3)


# V09（tests/verification.rs の v09_alternate_bars）が比べる理論の値（MORFAC 100）
V09_CASES = [
    # 幅 [m], 波長 [m], 成長率 [1/s], 移動速度 [m/s]
    (40.0, 211.4, 2.5895e-3, 0.33858),
    (40.0, 400.0, 1.6758e-3, 0.22031),
    (10.0, 142.5, -5.0585e-3, 0.31199),
]


@pytest.mark.parametrize(("width", "wavelength", "growth", "celerity"), V09_CASES)
def test_v09_reference_values(width, wavelength, growth, celerity):
    s, v = bt.growth_and_celerity(replace(SILL, morfac=100.0), wavelength, width)
    assert s == pytest.approx(growth, rel=2e-4)
    assert v == pytest.approx(celerity, rel=2e-4)


def test_alternate_bars_form_only_in_wide_channels():
    c = SILL
    wavelengths = np.geomspace(5.0, 500.0, 120)
    # 幅 10 m（β = 4.7）では全波長で減衰、40 m（β = 18.7）では育つ
    assert bt.most_unstable(c, 10.0, wavelengths)[1] < 0.0
    wavelength, growth = bt.most_unstable(c, 40.0, wavelengths)
    assert growth > 0.0
    assert 150.0 < wavelength < 300.0
    b = bt.base_state(c)
    beta_c = bt.critical_width(c, wavelengths) / (2 * b.depth)
    assert 7.5 < beta_c < 8.5


def test_moving_bars_travel_downstream():
    for width, wavelength in [(40.0, 100.0), (40.0, 211.4), (10.0, 142.5)]:
        assert bt.growth_and_celerity(SILL, wavelength, width)[1] > 0.0


def test_adaptation_length_matches_vormorph():
    # src/sediment/secondary_flow.rs の単体テストと同じ値（h = 1 m、n = 0.025、κ = 0.41）
    assert bt.adaptation_length(SILL, 1.0) == pytest.approx(9.625571, abs=1e-6)
    assert bt.adaptation_length(SILL, 3.0e-3) == 0.0
    b = bt.base_state(SILL)
    assert bt.adaptation_length(SILL, b.depth) == pytest.approx(10.48, abs=0.01)


@pytest.mark.parametrize("kind", ["transport", "equilibrium"])
def test_secondary_flow_without_deflection_changes_nothing(kind):
    k, width = 2 * np.pi / 211.4, 40.0
    reference = bt.bed_mode(replace(SILL, morfac=100.0), k, width)
    omega = bt.bed_mode(replace(SILL, morfac=100.0, secondary_flow=kind, n_star=0.0), k, width)
    assert omega == pytest.approx(reference, rel=1e-10)


def test_secondary_flow_approaches_the_quasi_steady_limit():
    k, width = 2 * np.pi / 211.4, 40.0
    small = replace(SILL, morfac=1e-3, secondary_flow="transport")
    assert bt.bed_mode(small, k, width) == pytest.approx(bt.quasi_steady_mode(small, k, width), rel=1e-4)


# V09 の二次流ありのケース（輸送方程式型、N* = 7、α_D = 0.15、MORFAC 100）
V09_SECONDARY_FLOW_CASES = [
    # 種類, 幅 [m], 波長 [m], 成長率 [1/s], 移動速度 [m/s]
    ("transport", 40.0, 211.4, 3.7512e-3, 0.35747),
    ("equilibrium", 40.0, 211.4, 3.9178e-3, 0.34800),
]


@pytest.mark.parametrize(("kind", "width", "wavelength", "growth", "celerity"), V09_SECONDARY_FLOW_CASES)
def test_v09_secondary_flow_reference_values(kind, width, wavelength, growth, celerity):
    c = replace(SILL, morfac=100.0, secondary_flow=kind)
    s, v = bt.growth_and_celerity(c, wavelength, width)
    assert s == pytest.approx(growth, rel=2e-4)
    assert v == pytest.approx(celerity, rel=2e-4)


def test_secondary_flow_enhances_alternate_bars():
    # この条件（β = 18.7）では、二次流で成長率が大きくなり、最も速く育つ波長が短くなる
    wavelengths = np.geomspace(20.0, 500.0, 120)
    plain = bt.most_unstable(SILL, 40.0, wavelengths)
    secondary = bt.most_unstable(replace(SILL, secondary_flow="transport"), 40.0, wavelengths)
    assert secondary[1] > plain[1]
    assert secondary[0] < plain[0]
