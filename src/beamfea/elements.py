"""Element-level quantities for the 2-node Euler-Bernoulli beam element.

Role in the pipeline
--------------------
Everything that happens *inside one element* lives here: Hermite cubic
shape functions and their derivatives, the 4x4 element stiffness matrix,
the 4x4 consistent mass matrix and the consistent load vector for a
uniformly distributed load. ``assembly.py`` then scatters these small
matrices into the global system.

Element DOF ordering: d = [w1, theta1, w2, theta2]
    w     transverse deflection [m]   (positive +y)
    theta rotation dw/dx [rad]        (positive counter-clockwise)
Local coordinate: xi = x_local / Le, with xi in [0, 1].
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def shape_functions(xi: float | FloatArray, length: float) -> FloatArray:
    """Hermite cubic shape functions N1..N4 for a beam element.

    Why Hermite cubics: the beam strain energy contains w'' (curvature), so
    the interpolation must keep *both* w and dw/dx continuous between
    elements (C1 continuity). Two nodes x two DOFs = four unknowns, which
    fixes a cubic polynomial uniquely.

    Args:
        xi: Local coordinate(s) xi = x/Le in [0, 1] [-].
        length: Element length Le [m].

    Returns:
        Array of shape (4,) or (4, n_points). N1, N3 are dimensionless;
        N2, N4 have units of metres (they multiply rotations in rad).
    """
    xi = np.asarray(xi, dtype=float)
    return np.array(
        [
            1.0 - 3.0 * xi**2 + 2.0 * xi**3,
            length * (xi - 2.0 * xi**2 + xi**3),
            3.0 * xi**2 - 2.0 * xi**3,
            length * (-(xi**2) + xi**3),
        ]
    )


def shape_function_first_derivatives(
    xi: float | FloatArray, length: float
) -> FloatArray:
    """First derivatives dN/dx of the Hermite shape functions.

    Used to recover the rotation field theta(x) = dw/dx inside an element.

    Args:
        xi: Local coordinate(s) xi = x/Le in [0, 1] [-].
        length: Element length Le [m].

    Returns:
        Array of shape (4,) or (4, n_points). Entries 1 and 3 have units
        [1/m]; entries 2 and 4 are dimensionless.
    """
    xi = np.asarray(xi, dtype=float)
    return np.array(
        [
            (-6.0 * xi + 6.0 * xi**2) / length,
            1.0 - 4.0 * xi + 3.0 * xi**2,
            (6.0 * xi - 6.0 * xi**2) / length,
            -2.0 * xi + 3.0 * xi**2,
        ]
    )


def shape_function_second_derivatives(
    xi: float | FloatArray, length: float
) -> FloatArray:
    """Second derivatives d2N/dx2 of the Hermite shape functions.

    Curvature kappa = w'' = sum_j N_j'' d_j, so these are the "strain
    displacement" functions (the B matrix) of the beam element.

    Args:
        xi: Local coordinate(s) xi = x/Le in [0, 1] [-].
        length: Element length Le [m].

    Returns:
        Array of shape (4,) or (4, n_points). Entries 1 and 3 have units
        [1/m^2]; entries 2 and 4 have units [1/m].
    """
    xi = np.asarray(xi, dtype=float)
    return np.array(
        [
            (-6.0 + 12.0 * xi) / length**2,
            (-4.0 + 6.0 * xi) / length,
            (6.0 - 12.0 * xi) / length**2,
            (-2.0 + 6.0 * xi) / length,
        ]
    )


def element_stiffness(flexural_rigidity: float, length: float) -> FloatArray:
    """Element stiffness matrix Ke of the Euler-Bernoulli beam.

    Obtained from the strain energy U = 1/2 * integral(EI * (w'')^2 dx):
    substituting w'' = B d gives Ke = integral(B^T EI B dx), which evaluates
    exactly to the closed form below.

    Args:
        flexural_rigidity: Bending stiffness EI [N*m^2].
        length: Element length Le [m].

    Returns:
        4x4 symmetric matrix. Units: [N/m] for w-w terms, [N] for w-theta
        terms, [N*m] for theta-theta terms.
    """
    le = length
    return (flexural_rigidity / le**3) * np.array(
        [
            [12.0, 6.0 * le, -12.0, 6.0 * le],
            [6.0 * le, 4.0 * le**2, -6.0 * le, 2.0 * le**2],
            [-12.0, -6.0 * le, 12.0, -6.0 * le],
            [6.0 * le, 2.0 * le**2, -6.0 * le, 4.0 * le**2],
        ]
    )


def element_consistent_mass(mass_per_length: float, length: float) -> FloatArray:
    """Consistent element mass matrix Me.

    "Consistent" means it uses the *same* Hermite shape functions as the
    stiffness: from the kinetic energy T = 1/2 * integral(rho*A * w_dot^2 dx)
    we get Me = integral(rho*A * N^T N dx).

    Args:
        mass_per_length: rho*A [kg/m].
        length: Element length Le [m].

    Returns:
        4x4 symmetric positive-definite matrix. Units: [kg] for w-w terms,
        [kg*m] for w-theta terms, [kg*m^2] for theta-theta terms.
    """
    le = length
    return (mass_per_length * le / 420.0) * np.array(
        [
            [156.0, 22.0 * le, 54.0, -13.0 * le],
            [22.0 * le, 4.0 * le**2, 13.0 * le, -3.0 * le**2],
            [54.0, 13.0 * le, 156.0, -22.0 * le],
            [-13.0 * le, -3.0 * le**2, -22.0 * le, 4.0 * le**2],
        ]
    )


def element_udl_load(intensity: float, length: float) -> FloatArray:
    """Consistent nodal load vector for a uniformly distributed load.

    fe = integral(N^T q dx). Note it contains end *moments* (+-q*Le^2/12),
    not just half the load at each node; these moments are what make the
    FE nodal values exact for a UDL.

    Args:
        intensity: Distributed load q, positive in +y [N/m].
        length: Element length Le [m].

    Returns:
        Array of shape (4,): [F1 [N], M1 [N*m], F2 [N], M2 [N*m]].
    """
    return intensity * length * np.array([0.5, length / 12.0, 0.5, -length / 12.0])
