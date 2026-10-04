"""Modal solver checks: accuracy, upper-bound property and convergence."""

import numpy as np
import pytest

from beamfea import analytical as an
from beamfea.assembly import assemble_beam
from beamfea.config import default_problem
from beamfea.modal import solve_modal
from beamfea.solver import cantilever_bcs

PROBLEM = default_problem()


def fe_frequencies(n: int, n_modes: int) -> np.ndarray:
    _, k, m = assemble_beam(PROBLEM.length, n, PROBLEM.flexural_rigidity, PROBLEM.mass_per_length)
    return solve_modal(k, m, cantilever_bcs(), n_modes).frequencies_hz


def exact_frequencies(n_modes: int) -> np.ndarray:
    return an.cantilever_natural_frequencies(
        PROBLEM.flexural_rigidity, PROBLEM.mass_per_length, PROBLEM.length, n_modes)


def test_first_mode_within_half_percent() -> None:
    assert fe_frequencies(16, 1)[0] == pytest.approx(exact_frequencies(1)[0], rel=5e-3)


def test_first_five_modes_are_upper_bounds() -> None:
    assert np.all(fe_frequencies(16, 5) > exact_frequencies(5))


def test_frequencies_strictly_increasing() -> None:
    assert np.all(np.diff(fe_frequencies(16, 10)) > 0.0)


def test_mode_one_error_decreases_with_refinement() -> None:
    exact = exact_frequencies(1)[0]
    errors = [fe_frequencies(n, 1)[0] - exact for n in (1, 2, 4, 8, 16, 32)]
    assert all(e > 0 for e in errors)
    assert all(b < a for a, b in zip(errors, errors[1:]))


def test_mode_shapes_mass_normalised() -> None:
    _, k, m = assemble_beam(PROBLEM.length, 10, PROBLEM.flexural_rigidity, PROBLEM.mass_per_length)
    phi = solve_modal(k, m, cantilever_bcs(), 4).mode_shapes
    np.testing.assert_allclose(phi.T @ m @ phi, np.eye(4), atol=1e-9)
    np.testing.assert_allclose(phi[:2, :], 0.0)  # clamped DOFs stay zero
