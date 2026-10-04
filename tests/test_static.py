"""Static solver checks against closed-form Euler-Bernoulli results."""

import numpy as np
import pytest

from beamfea import analytical as an
from beamfea.assembly import W, assemble_load, assemble_stiffness, generate_mesh
from beamfea.config import default_problem
from beamfea.postprocess import evaluate_fields, max_stress
from beamfea.solver import cantilever_bcs, simply_supported_bcs, solve_static

PROBLEM = default_problem()
P, L, EI = PROBLEM.tip_load, PROBLEM.length, PROBLEM.flexural_rigidity


def solve_tip(n: int):
    mesh = generate_mesh(L, n)
    k = assemble_stiffness(mesh, EI)
    f = assemble_load(mesh, point_loads=[(n, W, -P)])
    return mesh, solve_static(k, f, cantilever_bcs())


@pytest.mark.parametrize("n", [1, 2, 5, 20])
def test_tip_deflection_and_rotation_exact_for_any_mesh(n: int) -> None:
    # Exact solution is cubic -> Hermite cubics reproduce it, even with 1 element.
    _, sol = solve_tip(n)
    assert -sol.deflections[-1] == pytest.approx(an.cantilever_tip_deflection(P, L, EI), rel=1e-9)
    assert -sol.rotations[-1] == pytest.approx(an.cantilever_tip_rotation(P, L, EI), rel=1e-9)


@pytest.mark.parametrize("n", [1, 10])
def test_nodal_deflections_match_curve(n: int) -> None:
    mesh, sol = solve_tip(n)
    exact = -an.cantilever_deflection_curve(mesh.nodes, P, L, EI)
    np.testing.assert_allclose(sol.deflections, exact, rtol=1e-9, atol=1e-15)


def test_reactions_balance_load() -> None:
    _, sol = solve_tip(10)
    assert sol.reactions[(0, 0)] == pytest.approx(P, rel=1e-9)  # upward force
    assert sol.reactions[(0, 1)] == pytest.approx(P * L, rel=1e-9)  # CCW moment


def test_simply_supported_central_load() -> None:
    n = 10
    mesh = generate_mesh(L, n)
    k = assemble_stiffness(mesh, EI)
    f = assemble_load(mesh, point_loads=[(n // 2, W, -P)])
    sol = solve_static(k, f, simply_supported_bcs(mesh.n_nodes))
    assert -sol.deflections[n // 2] == pytest.approx(
        an.simply_supported_midspan_deflection(P, L, EI), rel=1e-9)
    assert sol.reactions[(0, 0)] == pytest.approx(P / 2, rel=1e-9)


@pytest.mark.parametrize("n", [2, 4, 10])
def test_udl_cantilever_tip(n: int) -> None:
    q = 200.0
    mesh = generate_mesh(L, n)
    k = assemble_stiffness(mesh, EI)
    sol = solve_static(k, assemble_load(mesh, udl_intensity=-q), cantilever_bcs())
    assert -sol.deflections[-1] == pytest.approx(an.cantilever_udl_tip_deflection(q, L, EI), rel=1e-9)


def test_root_stress_is_150_mpa() -> None:
    mesh, sol = solve_tip(PROBLEM.n_elements)
    fields = evaluate_fields(mesh, sol.displacements, PROBLEM.material.youngs_modulus,
                             PROBLEM.section.second_moment, PROBLEM.section.depth)
    peak, x_peak = max_stress(fields)
    assert peak == pytest.approx(150e6, rel=1e-9)
    assert x_peak == pytest.approx(0.0)
    assert fields.stress_top[0] > 0.0  # top fibre in tension for a downward tip load


def test_global_stiffness_singular_before_bcs() -> None:
    mesh = generate_mesh(L, 8)
    k = assemble_stiffness(mesh, EI)
    np.testing.assert_allclose(k, k.T)
    assert np.linalg.matrix_rank(k) == mesh.n_dofs - 2
