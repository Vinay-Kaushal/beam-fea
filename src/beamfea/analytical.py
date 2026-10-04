"""Closed-form Euler-Bernoulli beam results used to verify the FE solver.

Role in the pipeline
--------------------
These are the "hand calculations" at the bottom of the verification
ladder. Formulas return *magnitudes* for downward loads (the textbook
form); scripts convert to the FE sign convention (w positive +y) where
needed. All inputs and outputs are SI.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .elements import FloatArray

# Tabulated roots of the cantilever frequency equation 1 + cos(bL) cosh(bL) = 0.
CANTILEVER_BETA_L_TABLE: tuple[float, ...] = (
    1.87510407,
    4.69409113,
    7.85475744,
    10.99554073,
    14.13716839,
)


def cantilever_tip_deflection(load: float, length: float, ei: float) -> float:
    """Tip deflection of a cantilever with tip load, delta = P L^3 / (3 EI) [m].

    Args:
        load: Tip load magnitude P [N].
        length: Beam length L [m].
        ei: Bending stiffness EI [N*m^2].
    """
    return load * length**3 / (3.0 * ei)


def cantilever_tip_rotation(load: float, length: float, ei: float) -> float:
    """Tip rotation magnitude of a tip-loaded cantilever, P L^2 / (2 EI) [rad]."""
    return load * length**2 / (2.0 * ei)


def cantilever_deflection_curve(
    x: FloatArray, load: float, length: float, ei: float
) -> FloatArray:
    """Deflection magnitude w(x) = P x^2 (3L - x) / (6 EI) [m] for a tip load.

    Args:
        x: Positions along the beam [m].
        load: Tip load magnitude P [N].
        length: Beam length L [m].
        ei: Bending stiffness EI [N*m^2].
    """
    x = np.asarray(x, dtype=float)
    return load * x**2 * (3.0 * length - x) / (6.0 * ei)


def cantilever_moment(x: FloatArray, load: float, length: float) -> FloatArray:
    """Hogging bending-moment magnitude M(x) = P (L - x) [N*m] for a tip load."""
    return load * (length - np.asarray(x, dtype=float))


def bending_stress_magnitude(moment: float, c: float, second_moment: float) -> float:
    """Extreme-fibre bending stress |sigma| = M c / I [Pa]."""
    return abs(moment) * c / second_moment


def rectangular_max_shear_stress(shear_force: float, area: float) -> float:
    """Peak transverse shear stress in a rectangle, tau = 1.5 V / A [Pa]."""
    return 1.5 * abs(shear_force) / area


def cantilever_udl_tip_deflection(
    intensity: float, length: float, ei: float
) -> float:
    """Tip deflection of a cantilever under UDL, q L^4 / (8 EI) [m].

    Args:
        intensity: Distributed load magnitude q [N/m].
        length: Beam length L [m].
        ei: Bending stiffness EI [N*m^2].
    """
    return intensity * length**4 / (8.0 * ei)


def simply_supported_midspan_deflection(
    load: float, length: float, ei: float
) -> float:
    """Mid-span deflection, simply supported with central load, P L^3 / (48 EI) [m]."""
    return load * length**3 / (48.0 * ei)


def cantilever_beta_l(n_modes: int) -> FloatArray:
    """Roots beta_n L of 1 + cos(bL) cosh(bL) = 0 to full double precision.

    The tabulated values (8 decimals) are used as starting brackets and then
    refined with Brent's method. The refinement matters only for the
    convergence study, where FE errors become smaller than the rounding
    of the tabulated constants.

    Args:
        n_modes: Number of roots required [-].

    Returns:
        Array of beta_n L values [-].
    """
    roots = []
    for n in range(n_modes):
        # Beyond the table, beta_n L -> (2n - 1) pi / 2 (n starting at 1).
        guess = (
            CANTILEVER_BETA_L_TABLE[n]
            if n < len(CANTILEVER_BETA_L_TABLE)
            else (2 * (n + 1) - 1) * np.pi / 2.0
        )
        roots.append(
            brentq(
                lambda b: 1.0 + np.cos(b) * np.cosh(b),
                guess - 0.1,
                guess + 0.1,
                xtol=1e-15,
            )
        )
    return np.array(roots)


def cantilever_natural_frequencies(
    ei: float, mass_per_length: float, length: float, n_modes: int = 5
) -> FloatArray:
    """Bending natural frequencies of a cantilever [Hz].

    omega_n = (beta_n L)^2 * sqrt(EI / (rho A L^4)),  f_n = omega_n / (2 pi).

    Args:
        ei: Bending stiffness EI [N*m^2].
        mass_per_length: rho*A [kg/m].
        length: Beam length L [m].
        n_modes: Number of modes [-].
    """
    beta_l = cantilever_beta_l(n_modes)
    omega = beta_l**2 * np.sqrt(ei / (mass_per_length * length**4))
    return omega / (2.0 * np.pi)
