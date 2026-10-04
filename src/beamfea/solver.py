"""Boundary conditions, static solution and reaction recovery.

Role in the pipeline
--------------------
The assembled K is singular: an unsupported beam can translate and rotate
freely without straining (2 rigid-body modes). Boundary conditions remove
those modes. Here they are applied by **elimination**: constrained DOFs
(prescribed as zero) are simply deleted from the system, leaving the
smaller, non-singular free-free block K_ff u_f = F_f. After solving, the
full displacement vector is rebuilt and reactions are recovered from
R = K u - F at the constrained DOFs (an equilibrium check).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .assembly import THETA, W, global_dof
from .elements import FloatArray

BoundaryConditions = list[tuple[int, int]]
"""List of (node, local_dof) pairs constrained to zero."""


def cantilever_bcs() -> BoundaryConditions:
    """Fixed support at node 0: w = 0 and theta = 0."""
    return [(0, W), (0, THETA)]


def simply_supported_bcs(n_nodes: int) -> BoundaryConditions:
    """Pin + roller: w = 0 at the first and last node, rotations free."""
    return [(0, W), (n_nodes - 1, W)]


def constrained_dofs(bcs: BoundaryConditions) -> np.ndarray:
    """Sorted, unique global DOF indices that are constrained to zero."""
    return np.array(sorted({global_dof(n, d) for n, d in bcs}), dtype=int)


def free_dofs(n_dofs: int, bcs: BoundaryConditions) -> np.ndarray:
    """Global DOF indices that remain unknown after applying ``bcs``."""
    return np.setdiff1d(np.arange(n_dofs), constrained_dofs(bcs))


@dataclass(frozen=True)
class StaticSolution:
    """Result of a static solve.

    Attributes:
        displacements: Full DOF vector u = [w0, th0, w1, th1, ...]
            ([m] for w, [rad] for theta).
        reactions: Mapping (node, local_dof) -> reaction ([N] for W,
            [N*m] for THETA), positive in +y / counter-clockwise.
    """

    displacements: FloatArray
    reactions: dict[tuple[int, int], float]

    @property
    def deflections(self) -> FloatArray:
        """Nodal deflections w [m]."""
        return self.displacements[0::2]

    @property
    def rotations(self) -> FloatArray:
        """Nodal rotations theta [rad]."""
        return self.displacements[1::2]


def compute_reactions(
    k_full: FloatArray,
    u_full: FloatArray,
    f_full: FloatArray,
    bcs: BoundaryConditions,
) -> dict[tuple[int, int], float]:
    """Recover support reactions R = K u - F at the constrained DOFs.

    Why this works: at every DOF, internal (elastic) force K u must balance
    external load F plus any support reaction R. At free DOFs R = 0; at
    constrained DOFs the residual K u - F is exactly what the support supplies.

    Args:
        k_full: Global stiffness before BCs.
        u_full: Full displacement vector [m, rad].
        f_full: Global applied load vector [N, N*m].
        bcs: Boundary conditions used in the solve.

    Returns:
        Mapping (node, local_dof) -> reaction [N or N*m].
    """
    residual = k_full @ u_full - f_full
    return {(n, d): float(residual[global_dof(n, d)]) for n, d in bcs}


def solve_static(
    k_full: FloatArray, f_full: FloatArray, bcs: BoundaryConditions
) -> StaticSolution:
    """Solve K u = F with zero-displacement BCs applied by elimination.

    Args:
        k_full: Global stiffness matrix before BCs.
        f_full: Global load vector [N, N*m].
        bcs: Constrained (node, local_dof) pairs.

    Returns:
        ``StaticSolution`` with the full displacement vector and reactions.
    """
    free = free_dofs(k_full.shape[0], bcs)
    # Partition: only the free-free block is needed because the prescribed
    # displacements are zero (so the K_fc * u_c term vanishes).
    k_ff = k_full[np.ix_(free, free)]
    u_full = np.zeros(k_full.shape[0])
    u_full[free] = np.linalg.solve(k_ff, f_full[free])
    reactions = compute_reactions(k_full, u_full, f_full, bcs)
    return StaticSolution(displacements=u_full, reactions=reactions)
