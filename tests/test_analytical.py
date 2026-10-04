"""Hand-calculation values for the reference problem (spec Section 2)."""

import numpy as np
import pytest

from beamfea import analytical as an
from beamfea.config import default_problem

PROBLEM = default_problem()
P, L, EI = PROBLEM.tip_load, PROBLEM.length, PROBLEM.flexural_rigidity
SEC = PROBLEM.section


def test_second_moment_of_area() -> None:
    assert SEC.second_moment * 1e12 == pytest.approx(1666.67, abs=0.01)  # mm^4


def test_tip_deflection_12_5_mm() -> None:
    assert an.cantilever_tip_deflection(P, L, EI) * 1e3 == pytest.approx(12.5, rel=1e-12)


def test_tip_rotation() -> None:
    assert an.cantilever_tip_rotation(P, L, EI) == pytest.approx(0.0375, rel=1e-12)


def test_root_stress_150_mpa() -> None:
    m_root = float(an.cantilever_moment(0.0, P, L))
    assert m_root == pytest.approx(50.0)
    assert an.bending_stress_magnitude(m_root, SEC.extreme_fibre, SEC.second_moment) / 1e6 == pytest.approx(150.0)


def test_midspan_stress_75_mpa() -> None:
    m_mid = float(an.cantilever_moment(L / 2, P, L))
    assert an.bending_stress_magnitude(m_mid, SEC.extreme_fibre, SEC.second_moment) / 1e6 == pytest.approx(75.0)


def test_shear_stress_negligible() -> None:
    assert an.rectangular_max_shear_stress(P, SEC.area) / 1e6 == pytest.approx(0.75)


def test_beta_l_roots_match_table() -> None:
    np.testing.assert_allclose(an.cantilever_beta_l(5), an.CANTILEVER_BETA_L_TABLE, atol=5e-9)


def test_first_frequency_hand_calculation() -> None:
    # Independent hand calculation: f1 = 1.8751^2 / (2 pi) * sqrt(EI / (rho A L^4)).
    # EI = 333.33 N*m^2, rho*A = 1.57 kg/m, L^4 = 0.0625 m^4 -> f1 = 32.6 Hz.
    # NOTE: the project spec quoted ~103 Hz; that value is off by sqrt(10)
    # (see docs/ASSUMPTIONS.md). The physics gives 32.6 Hz.
    hand = 1.87510407**2 / (2 * np.pi) * np.sqrt(333.3333333 / (7850 * 2e-4 * 0.5**4))
    f = an.cantilever_natural_frequencies(EI, PROBLEM.mass_per_length, L, 3)
    assert f[0] == pytest.approx(hand, abs=0.01)
    assert f[0] == pytest.approx(32.6, abs=0.5)


def test_frequency_ratios_match_beam_theory() -> None:
    f = an.cantilever_natural_frequencies(EI, PROBLEM.mass_per_length, L, 3)
    assert f[1] / f[0] == pytest.approx(6.267, abs=1e-3)
    assert f[2] / f[0] == pytest.approx(17.548, abs=1e-3)
