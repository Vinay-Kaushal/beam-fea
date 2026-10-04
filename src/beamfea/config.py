"""Problem definition: material, cross-section and beam configuration.

Role in the pipeline
--------------------
This is the *input layer*. Every other module works purely in SI units
(m, N, Pa, kg, s). Values that engineers naturally quote in millimetres or
megapascals are converted to SI here, exactly once, via the helper functions
below. All dataclasses are frozen so a problem definition can never be
mutated part-way through an analysis.
"""

from __future__ import annotations

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Unit conversion helpers (input layer only)
# ---------------------------------------------------------------------------


def mm_to_m(value_mm: float) -> float:
    """Convert a length from millimetres [mm] to metres [m]."""
    return value_mm * 1.0e-3


def m_to_mm(value_m: float) -> float:
    """Convert a length from metres [m] to millimetres [mm] (for reporting)."""
    return value_m * 1.0e3


def gpa_to_pa(value_gpa: float) -> float:
    """Convert a modulus from gigapascals [GPa] to pascals [Pa]."""
    return value_gpa * 1.0e9


def mpa_to_pa(value_mpa: float) -> float:
    """Convert a stress from megapascals [MPa] to pascals [Pa]."""
    return value_mpa * 1.0e6


def pa_to_mpa(value_pa: float) -> float:
    """Convert a stress from pascals [Pa] to megapascals [MPa] (for reporting)."""
    return value_pa * 1.0e-6


# ---------------------------------------------------------------------------
# Problem dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Material:
    """Linear-elastic isotropic material.

    Attributes:
        name: Human-readable material name.
        youngs_modulus: Young's modulus E [Pa].
        poissons_ratio: Poisson's ratio nu [-]. Not used by Euler-Bernoulli
            theory; kept so the 1D model and the 3D ANSYS model share one
            material definition.
        density: Mass density rho [kg/m^3].
        yield_strength: Yield strength [Pa], used only for reporting a
            safety factor.
    """

    name: str
    youngs_modulus: float
    poissons_ratio: float
    density: float
    yield_strength: float


@dataclass(frozen=True)
class RectangularSection:
    """Solid rectangular cross-section.

    Bending happens in the x-y plane, so ``depth`` (h) is measured along y
    and ``width`` (b) along z (out of plane).

    Attributes:
        width: Width b, out of the bending plane [m].
        depth: Depth h, in the bending direction (y) [m].
    """

    width: float
    depth: float

    @property
    def area(self) -> float:
        """Cross-sectional area A = b*h [m^2]."""
        return self.width * self.depth

    @property
    def second_moment(self) -> float:
        """Second moment of area about z, I = b*h^3/12 [m^4]."""
        return self.width * self.depth**3 / 12.0

    @property
    def extreme_fibre(self) -> float:
        """Distance from neutral axis to extreme fibre, c = h/2 [m]."""
        return self.depth / 2.0


@dataclass(frozen=True)
class BeamProblem:
    """Complete definition of the reference cantilever problem.

    Sign convention: deflection w is positive in +y, rotation theta = dw/dx
    is positive counter-clockwise, and the tip load acts in -y (downward).
    ``tip_load`` stores the load *magnitude* P [N]; the solver applies it
    as -P on the tip w-DOF.

    Attributes:
        length: Beam length L [m].
        material: Material properties.
        section: Cross-section properties.
        tip_load: Magnitude P of the downward tip point load [N].
        n_elements: Default number of finite elements for static runs [-].
    """

    length: float
    material: Material
    section: RectangularSection
    tip_load: float
    n_elements: int = 10

    @property
    def flexural_rigidity(self) -> float:
        """Bending stiffness EI [N*m^2]."""
        return self.material.youngs_modulus * self.section.second_moment

    @property
    def mass_per_length(self) -> float:
        """Mass per unit length rho*A [kg/m]."""
        return self.material.density * self.section.area


def structural_steel() -> Material:
    """Return the structural steel used in the reference problem (SI units)."""
    return Material(
        name="Structural steel",
        youngs_modulus=gpa_to_pa(200.0),
        poissons_ratio=0.3,
        density=7850.0,
        yield_strength=mpa_to_pa(250.0),
    )


def default_problem() -> BeamProblem:
    """Return the reference cantilever of Section 2 of the spec (SI units).

    L = 500 mm, b = 20 mm, h = 10 mm, steel, P = 100 N downward at the tip.
    """
    return BeamProblem(
        length=mm_to_m(500.0),
        material=structural_steel(),
        section=RectangularSection(width=mm_to_m(20.0), depth=mm_to_m(10.0)),
        tip_load=100.0,
        n_elements=10,
    )


@dataclass(frozen=True)
class StudySettings:
    """Mesh settings used by the scripts (kept here so nothing is hard-coded twice).

    Attributes:
        static_mesh_sweep: Element counts used to show tip-load exactness [-].
        modal_elements: Element count for the main modal run [-].
        modal_n_modes: Number of modes compared with theory [-].
        convergence_meshes: Element counts for the modal h-refinement study [-].
        convergence_n_modes: Number of modes tracked in the convergence study [-].
        udl_intensity: Distributed load magnitude q for the UDL case [N/m]
            (applied downward).
        points_per_element: Sampling points per element for field plots [-].
    """

    static_mesh_sweep: tuple[int, ...] = (1, 2, 5, 10, 20, 50)
    modal_elements: int = 20
    modal_n_modes: int = 5
    convergence_meshes: tuple[int, ...] = (1, 2, 4, 8, 16, 32, 64)
    convergence_n_modes: int = 3
    udl_intensity: float = 200.0
    points_per_element: int = 11
