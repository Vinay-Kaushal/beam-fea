"""Element-level checks: symmetry, rigid-body modes, mass, shape functions."""

import numpy as np
import pytest

from beamfea.elements import (
    element_consistent_mass,
    element_stiffness,
    element_udl_load,
    shape_function_second_derivatives,
    shape_functions,
)

EI = 333.333  # N*m^2
RHO_A = 1.57  # kg/m
LE = 0.05  # m


def test_stiffness_is_symmetric() -> None:
    ke = element_stiffness(EI, LE)
    np.testing.assert_allclose(ke, ke.T)


def test_stiffness_has_two_rigid_body_modes() -> None:
    # A free element can translate and rotate without strain energy,
    # so exactly 2 eigenvalues must be (numerically) zero.
    eig = np.linalg.eigvalsh(element_stiffness(EI, LE))
    tol = 1e-10 * np.max(np.abs(eig))
    assert np.sum(np.abs(eig) < tol) == 2
    assert np.all(eig > -tol)


def test_rigid_translation_and_rotation_produce_no_force() -> None:
    ke = element_stiffness(EI, LE)
    translation = np.array([1.0, 0.0, 1.0, 0.0])
    rotation = np.array([0.0, 1.0, LE, 1.0])  # w = x, theta = 1
    np.testing.assert_allclose(ke @ translation, 0.0, atol=1e-9)
    np.testing.assert_allclose(ke @ rotation, 0.0, atol=1e-9)


def test_mass_is_symmetric_positive_definite() -> None:
    me = element_consistent_mass(RHO_A, LE)
    np.testing.assert_allclose(me, me.T)
    assert np.all(np.linalg.eigvalsh(me) > 0.0)


def test_mass_rigid_translation_gives_total_mass() -> None:
    v = np.array([1.0, 0.0, 1.0, 0.0])
    assert v @ element_consistent_mass(RHO_A, LE) @ v == pytest.approx(RHO_A * LE, rel=1e-12)


@pytest.mark.parametrize("xi", np.linspace(0.0, 1.0, 7))
def test_partition_of_unity_for_translation(xi: float) -> None:
    n = shape_functions(xi, LE)
    assert n[0] + n[2] == pytest.approx(1.0, abs=1e-14)


def test_shape_functions_interpolate_nodal_values() -> None:
    np.testing.assert_allclose(shape_functions(0.0, LE), [1.0, 0.0, 0.0, 0.0], atol=1e-15)
    np.testing.assert_allclose(shape_functions(1.0, LE), [0.0, 0.0, 1.0, 0.0], atol=1e-15)


def test_curvature_of_pure_bending_is_constant() -> None:
    # w = x^2/2 has curvature 1 everywhere: d = [0, 0, LE^2/2, LE].
    d = np.array([0.0, 0.0, LE**2 / 2.0, LE])
    xi = np.linspace(0.0, 1.0, 5)
    np.testing.assert_allclose(d @ shape_function_second_derivatives(xi, LE), 1.0)


def test_udl_load_sums_to_total_load() -> None:
    fe = element_udl_load(-200.0, LE)
    assert fe[0] + fe[2] == pytest.approx(-200.0 * LE)
    assert fe[1] == pytest.approx(-fe[3])
