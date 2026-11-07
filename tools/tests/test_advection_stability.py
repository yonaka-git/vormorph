"""移流と時間積分の1次元の線形安定解析（tools/advection_stability.py）のテスト。

移流の行列が移流項 −U ∂u/∂x の1次・2次の近似になっていることと、ADR-0016・ADR-0017 が根拠に
する性質（予測子・修正子と1ステップ全体の SSP-RK2 はクーラン数 1 まで安定、2次の移流を
前進 Euler で進めると不安定、SSP-RK2 は重力波を余分に減衰させ、予測子・修正子は θ 法と同じ）を
確かめる。
"""

import numpy as np
import pytest

import advection_stability as stab


@pytest.mark.parametrize("scheme, order", [("upwind", 1.0), ("second_order", 2.0)])
def test_advection_is_consistent_with_the_expected_order(scheme, order):
    # 辺の流速 sin(2πx/λ) に移流を施した (F u − u)/C は、−Δx ∂u/∂x を近似する
    errors = []
    for n in (32, 64, 128):
        x = np.arange(n) / n  # 周期 1、Δx = 1/n（辺 j は x = j/n）
        u = np.sin(2 * np.pi * x)
        courant = 1e-6
        rate = (stab.advected(u, courant, scheme) - u) / courant
        exact = -(1.0 / n) * 2 * np.pi * np.cos(2 * np.pi * x)
        errors.append(np.abs(rate - exact).max() / np.abs(exact).max())
    observed = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert observed == pytest.approx([order, order], abs=0.1)


def test_uniform_flow_is_unchanged_by_advection():
    u = np.full(16, 0.7)
    for scheme in stab.ADVECTION:
        assert stab.advected(u, 0.5, scheme) == pytest.approx(u, abs=1e-15)


def test_advection_matrix_matches_the_function():
    rng = np.random.default_rng(1)
    u = rng.normal(size=12)
    for scheme in stab.ADVECTION:
        matrix = stab.advection_matrix(12, 0.4, scheme)
        assert matrix @ u == pytest.approx(stab.advected(u, 0.4, scheme), abs=1e-14)


@pytest.mark.parametrize("method", ["predictor_corrector", "ssp_rk2"])
@pytest.mark.parametrize("scheme", stab.ADVECTION)
@pytest.mark.parametrize("froude", [0.2, 0.58, 0.95])
def test_two_stage_methods_are_stable_up_to_courant_one(method, scheme, froude):
    for courant in (0.25, 0.5, 0.8, 1.0):
        for friction in (False, True):
            radius = stab.spectral_radius(scheme, method, froude, courant, friction=friction, n=24)
            assert radius <= 1.0 + 1e-9, (courant, friction, radius)


def test_predictor_corrector_damps_gravity_waves_like_the_theta_method():
    # 1周期 200 ステップ（V04 と同じ）
    euler = stab.gravity_wave_amplitude("euler", 200)
    assert euler == pytest.approx(0.990, abs=1e-3)
    assert stab.gravity_wave_amplitude("predictor_corrector", 200) == pytest.approx(euler, abs=1e-12)
    # SSP-RK2 は時間について1次で余分に減衰する：ステップ数を2倍にすると減衰の不足分が約半分
    ssp = [1.0 - stab.gravity_wave_amplitude("ssp_rk2", steps) for steps in (200, 400)]
    assert ssp[0] == pytest.approx(0.103, abs=2e-3)
    assert ssp[1] / ssp[0] == pytest.approx(0.5, abs=0.03)
    # θ = 0.5 なら θ 法は減衰しないが、SSP-RK2 は減衰する
    assert stab.gravity_wave_amplitude("euler", 200, theta=0.5) == pytest.approx(1.0, abs=1e-9)
    assert stab.gravity_wave_amplitude("ssp_rk2", 200, theta=0.5) < 0.91


def test_second_order_advection_is_unstable_with_forward_euler():
    # V09 の等流（Fr 0.58）でも、クーラン数 0.3 で不安定
    assert stab.spectral_radius("second_order", "euler", 0.58, 0.3, n=24) > 1.001
    # 面の値を Δt/2 上流で評価しても、移流だけを SSP-RK2 にしても、限界は 0.3 程度
    for method in ("lax_wendroff", "advection_rk2"):
        assert stab.spectral_radius("second_order", method, 0.58, 0.4, n=24) > 1.001


def test_first_order_limits_match_adr_0016():
    # 1次の移流＋1段：深い水（Fr 0.2）では C 0.5 でわずかに不安定、等流（Fr 0.58）では安定
    assert stab.spectral_radius("upwind", "euler", 0.2, 0.5, n=24) > 1.001
    assert stab.spectral_radius("upwind", "euler", 0.58, 0.5, n=24) <= 1.0 + 1e-9
    assert stab.critical_courant("upwind", "euler", 0.58) == pytest.approx(0.56, abs=0.01)


def test_unknown_names_are_rejected():
    with pytest.raises(ValueError):
        stab.spectral_radius("central", "euler", 0.5, 0.5, n=8)
    with pytest.raises(ValueError):
        stab.spectral_radius("upwind", "rk4", 0.5, 0.5, n=8)
    with pytest.raises(ValueError):
        stab.spectral_radius("upwind", "lax_wendroff", 0.5, 0.5, n=8)
    with pytest.raises(ValueError):
        stab.gravity_wave_amplitude("ab3", 100)
