"""Modal analysis: the generalized eigenproblem K phi = omega^2 M phi.

Role in the pipeline
--------------------
Free vibration of the undamped beam obeys M u'' + K u = 0. Assuming
harmonic motion u = phi * sin(omega t) gives (K - omega^2 M) phi = 0, a
*generalized* symmetric eigenproblem. After eliminating the constrained
DOFs (same as in the static solve) we hand K_ff and M_ff to
``scipy.linalg.eigh``, which returns eigenvalues lambda = omega^2 in
ascending order and mass-normalised eigenvectors (phi^T M phi = 1).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh

from .elements import FloatArray
from .solver import BoundaryConditions, free_dofs


@dataclass(frozen=True)
class ModalSolution:
    """Natural frequencies and mode shapes.

    Attributes:
        eigenvalues: lambda = omega^2, ascending [rad^2/s^2].
        mode_shapes: Full-DOF mass-normalised mode shapes, shape
            (n_dofs, n_modes); constrained DOFs are zero.
    """

    eigenvalues: FloatArray
    mode_shapes: FloatArray

    @property
    def angular_frequencies(self) -> FloatArray:
        """omega = sqrt(lambda) [rad/s]."""
        return np.sqrt(self.eigenvalues)

    @property
    def frequencies_hz(self) -> FloatArray:
        """f = omega / (2 pi) [Hz]."""
        return self.angular_frequencies / (2.0 * np.pi)


def solve_modal(
    k_full: FloatArray,
    m_full: FloatArray,
    bcs: BoundaryConditions,
    n_modes: int | None = None,
) -> ModalSolution:
    """Solve the constrained generalized eigenproblem.

    Args:
        k_full: Global stiffness before BCs.
        m_full: Global mass before BCs.
        bcs: Constrained (node, local_dof) pairs.
        n_modes: Number of lowest modes to keep (None = all available).

    Returns:
        ``ModalSolution`` sorted by ascending frequency.
    """
    free = free_dofs(k_full.shape[0], bcs)
    k_ff = k_full[np.ix_(free, free)]
    m_ff = m_full[np.ix_(free, free)]
    eigenvalues, vectors = eigh(k_ff, m_ff)  # ascending, M-orthonormal
    order = np.argsort(eigenvalues)
    eigenvalues, vectors = eigenvalues[order], vectors[:, order]
    if n_modes is not None:
        eigenvalues, vectors = eigenvalues[:n_modes], vectors[:, :n_modes]

    shapes = np.zeros((k_full.shape[0], eigenvalues.size))
    shapes[free, :] = vectors
    # Eigenvectors are only defined up to sign; fix it (largest |w| > 0)
    # so plots and CSVs are deterministic between runs and machines.
    for j in range(shapes.shape[1]):
        w = shapes[0::2, j]
        if w[np.argmax(np.abs(w))] < 0.0:
            shapes[:, j] *= -1.0
    return ModalSolution(eigenvalues=eigenvalues, mode_shapes=shapes)
