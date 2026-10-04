"""Mesh generation and assembly of the global K, M and F.

Role in the pipeline
--------------------
Takes the small element matrices from ``elements.py`` and "scatters" them
into global matrices. Each node has two DOFs, so node i owns global
indices 2*i (w) and 2*i + 1 (theta). Neighbouring elements share a node,
which is where their contributions add up - that addition *is* the
enforcement of displacement continuity and force equilibrium at the node.

Dense NumPy arrays are used because the meshes here are tiny (< 200 DOFs).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .elements import (
    FloatArray,
    element_consistent_mass,
    element_stiffness,
    element_udl_load,
)

DOFS_PER_NODE = 2
W = 0  # local DOF index for deflection w
THETA = 1  # local DOF index for rotation theta


def global_dof(node: int, local_dof: int) -> int:
    """Return the global DOF index of ``local_dof`` (W or THETA) at ``node``."""
    return DOFS_PER_NODE * node + local_dof


@dataclass(frozen=True)
class Mesh:
    """Uniform 1D beam mesh.

    Attributes:
        nodes: Node x-coordinates, shape (n_nodes,) [m].
        connectivity: Node indices of each element, shape (n_elements, 2).
    """

    nodes: FloatArray
    connectivity: np.ndarray

    @property
    def n_nodes(self) -> int:
        """Number of nodes [-]."""
        return int(self.nodes.size)

    @property
    def n_elements(self) -> int:
        """Number of elements [-]."""
        return int(self.connectivity.shape[0])

    @property
    def n_dofs(self) -> int:
        """Total number of DOFs before boundary conditions [-]."""
        return DOFS_PER_NODE * self.n_nodes

    def element_length(self, element: int) -> float:
        """Length of element ``element`` [m]."""
        i, j = self.connectivity[element]
        return float(self.nodes[j] - self.nodes[i])

    def element_dofs(self, element: int) -> list[int]:
        """Global DOF indices [w1, th1, w2, th2] of ``element``."""
        i, j = self.connectivity[element]
        return [
            global_dof(i, W),
            global_dof(i, THETA),
            global_dof(j, W),
            global_dof(j, THETA),
        ]


def generate_mesh(length: float, n_elements: int) -> Mesh:
    """Create a uniform mesh of ``n_elements`` beam elements on [0, length].

    Args:
        length: Beam length L [m].
        n_elements: Number of elements (>= 1) [-].

    Returns:
        A ``Mesh`` with n_elements + 1 equally spaced nodes.
    """
    if n_elements < 1:
        raise ValueError("n_elements must be >= 1")
    nodes = np.linspace(0.0, length, n_elements + 1)
    connectivity = np.column_stack(
        [np.arange(n_elements), np.arange(1, n_elements + 1)]
    )
    return Mesh(nodes=nodes, connectivity=connectivity)


def assemble_stiffness(mesh: Mesh, flexural_rigidity: float) -> FloatArray:
    """Assemble the global stiffness matrix K (before boundary conditions).

    Args:
        mesh: Beam mesh.
        flexural_rigidity: Bending stiffness EI [N*m^2].

    Returns:
        (n_dofs, n_dofs) symmetric, singular matrix (2 rigid-body modes).
    """
    k_global = np.zeros((mesh.n_dofs, mesh.n_dofs))
    for e in range(mesh.n_elements):
        dofs = mesh.element_dofs(e)
        # np.ix_ adds the 4x4 block into the rows/columns of this element.
        k_global[np.ix_(dofs, dofs)] += element_stiffness(
            flexural_rigidity, mesh.element_length(e)
        )
    return k_global


def assemble_mass(mesh: Mesh, mass_per_length: float) -> FloatArray:
    """Assemble the global consistent mass matrix M.

    Args:
        mesh: Beam mesh.
        mass_per_length: rho*A [kg/m].

    Returns:
        (n_dofs, n_dofs) symmetric positive-definite matrix.
    """
    m_global = np.zeros((mesh.n_dofs, mesh.n_dofs))
    for e in range(mesh.n_elements):
        dofs = mesh.element_dofs(e)
        m_global[np.ix_(dofs, dofs)] += element_consistent_mass(
            mass_per_length, mesh.element_length(e)
        )
    return m_global


def assemble_load(
    mesh: Mesh,
    point_loads: list[tuple[int, int, float]] | None = None,
    udl_intensity: float = 0.0,
) -> FloatArray:
    """Assemble the global load vector F.

    Args:
        mesh: Beam mesh.
        point_loads: List of (node, local_dof, value). For local_dof = W the
            value is a force [N] (positive +y); for THETA it is a moment
            [N*m] (positive counter-clockwise). Point loads go straight into
            F because they act exactly at a node.
        udl_intensity: Uniform distributed load q over the whole beam,
            positive +y [N/m]. Converted to consistent nodal loads per element.

    Returns:
        (n_dofs,) load vector [N or N*m].
    """
    f_global = np.zeros(mesh.n_dofs)
    for node, local_dof, value in point_loads or []:
        f_global[global_dof(node, local_dof)] += value
    if udl_intensity != 0.0:
        for e in range(mesh.n_elements):
            f_global[mesh.element_dofs(e)] += element_udl_load(
                udl_intensity, mesh.element_length(e)
            )
    return f_global


def assemble_beam(
    length: float, n_elements: int, flexural_rigidity: float, mass_per_length: float
) -> tuple[Mesh, FloatArray, FloatArray]:
    """Convenience wrapper: uniform mesh plus global K and M in one call.

    Args:
        length: Beam length L [m].
        n_elements: Number of elements [-].
        flexural_rigidity: EI [N*m^2].
        mass_per_length: rho*A [kg/m].

    Returns:
        (mesh, K, M) with K, M assembled before boundary conditions.
    """
    mesh = generate_mesh(length, n_elements)
    return (
        mesh,
        assemble_stiffness(mesh, flexural_rigidity),
        assemble_mass(mesh, mass_per_length),
    )
