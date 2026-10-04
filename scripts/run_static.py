"""Static verification: FE solver vs closed-form beam theory.

Pipeline step 1 of run_all.py. Solves the reference cantilever (tip load),
shows the tip deflection is exact for any mesh, checks reactions and
stresses, and runs two extra validation cases (UDL cantilever, simply
supported beam with a central load).

Outputs: results/static_results.csv, results/static_deflection.png,
         results/static_moment_stress.png
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
RESULTS = ROOT / "results"

from beamfea import analytical as an  # noqa: E402
from beamfea.assembly import W, assemble_load, generate_mesh, assemble_stiffness  # noqa: E402
from beamfea.config import BeamProblem, StudySettings, default_problem, m_to_mm, pa_to_mpa  # noqa: E402
from beamfea.plotting import plot_moment_stress, plot_static_deflection  # noqa: E402
from beamfea.postprocess import (  # noqa: E402
    BeamFields,
    evaluate_fields,
    format_table,
    max_stress,
    relative_error_pct,
    value_at,
    write_csv,
)
from beamfea.solver import (  # noqa: E402
    StaticSolution,
    cantilever_bcs,
    simply_supported_bcs,
    solve_static,
)

CSV_HEADERS = [
    "case", "quantity", "unit", "n_elements",
    "fe_value", "analytical_value", "abs_error", "rel_error_pct",
]


@dataclass(frozen=True)
class ComparisonRow:
    """One FE-vs-analytical comparison (values in reporting units)."""

    case: str
    quantity: str
    unit: str
    n_elements: int
    fe_value: float
    analytical_value: float

    @property
    def abs_error(self) -> float:
        """|FE - analytical| [same unit as the values]."""
        return abs(self.fe_value - self.analytical_value)

    @property
    def rel_error_pct(self) -> float:
        """Signed relative error [%]."""
        return relative_error_pct(self.fe_value, self.analytical_value)


def solve_cantilever_tip(problem: BeamProblem, n_elements: int) -> tuple[StaticSolution, BeamFields]:
    """Solve the tip-loaded cantilever and recover its fields.

    Args:
        problem: Beam definition (SI).
        n_elements: Number of elements [-].

    Returns:
        (static solution, sampled fields).
    """
    mesh = generate_mesh(problem.length, n_elements)
    k = assemble_stiffness(mesh, problem.flexural_rigidity)
    # Tip load acts in -y, so the nodal force is -P on the last w-DOF.
    f = assemble_load(mesh, point_loads=[(n_elements, W, -problem.tip_load)])
    sol = solve_static(k, f, cantilever_bcs())
    fields = evaluate_fields(
        mesh, sol.displacements, problem.material.youngs_modulus,
        problem.section.second_moment, problem.section.depth,
        StudySettings().points_per_element,
    )
    return sol, fields


def tip_load_rows(problem: BeamProblem, settings: StudySettings) -> list[ComparisonRow]:
    """Tip deflection for every mesh in the sweep (shows exactness)."""
    exact_mm = m_to_mm(an.cantilever_tip_deflection(problem.tip_load, problem.length, problem.flexural_rigidity))
    rows = []
    for n in settings.static_mesh_sweep:
        sol, _ = solve_cantilever_tip(problem, n)
        # Report the downward deflection as a positive magnitude, like the formula.
        rows.append(ComparisonRow("cantilever_tip_load", "tip_deflection", "mm", n,
                                  -m_to_mm(sol.deflections[-1]), exact_mm))
    return rows


def detail_rows(problem: BeamProblem, sol: StaticSolution, fields: BeamFields, n: int) -> list[ComparisonRow]:
    """Rotation, reactions and stresses for the default mesh."""
    p, length, ei = problem.tip_load, problem.length, problem.flexural_rigidity
    c, i_z = problem.section.extreme_fibre, problem.section.second_moment
    root_moment = an.cantilever_moment(0.0, p, length)
    mid_moment = an.cantilever_moment(length / 2.0, p, length)
    case = "cantilever_tip_load"
    return [
        ComparisonRow(case, "tip_rotation", "rad", n, -sol.rotations[-1],
                      an.cantilever_tip_rotation(p, length, ei)),
        ComparisonRow(case, "reaction_force", "N", n, sol.reactions[(0, 0)], p),
        ComparisonRow(case, "reaction_moment", "N*m", n, sol.reactions[(0, 1)], float(root_moment)),
        ComparisonRow(case, "root_top_stress", "MPa", n, pa_to_mpa(fields.stress_top[0]),
                      pa_to_mpa(an.bending_stress_magnitude(float(root_moment), c, i_z))),
        ComparisonRow(case, "midspan_top_stress", "MPa", n,
                      pa_to_mpa(value_at(fields.x, fields.stress_top, length / 2.0)),
                      pa_to_mpa(an.bending_stress_magnitude(float(mid_moment), c, i_z))),
    ]


def extra_case_rows(problem: BeamProblem, settings: StudySettings, n: int) -> list[ComparisonRow]:
    """UDL cantilever and simply supported central-load validation cases."""
    ei, length = problem.flexural_rigidity, problem.length
    mesh = generate_mesh(length, n)
    k = assemble_stiffness(mesh, ei)

    q = settings.udl_intensity
    sol_udl = solve_static(k, assemble_load(mesh, udl_intensity=-q), cantilever_bcs())

    mid = n // 2  # n is even, so the mid-span point is a node
    f_ss = assemble_load(mesh, point_loads=[(mid, W, -problem.tip_load)])
    sol_ss = solve_static(k, f_ss, simply_supported_bcs(mesh.n_nodes))
    return [
        ComparisonRow("cantilever_udl", "tip_deflection", "mm", n,
                      -m_to_mm(sol_udl.deflections[-1]),
                      m_to_mm(an.cantilever_udl_tip_deflection(q, length, ei))),
        ComparisonRow("simply_supported_central_load", "midspan_deflection", "mm", n,
                      -m_to_mm(sol_ss.deflections[mid]),
                      m_to_mm(an.simply_supported_midspan_deflection(problem.tip_load, length, ei))),
    ]


def make_plots(problem: BeamProblem, sol: StaticSolution, fields: BeamFields, n: int) -> None:
    """Deflection plot and moment/stress plot for the default mesh."""
    p, length, ei = problem.tip_load, problem.length, problem.flexural_rigidity
    x = np.linspace(0.0, length, 201)
    nodes = np.linspace(0.0, length, n + 1)
    plot_static_deflection(
        m_to_mm(nodes), m_to_mm(sol.deflections),
        m_to_mm(x), -m_to_mm(an.cantilever_deflection_curve(x, p, length, ei)),
        RESULTS / "static_deflection.png",
        f"Cantilever, tip load P = {p:.0f} N: FE ({n} elements) vs analytical",
    )
    # Analytical moment is hogging (negative in the sagging-positive convention).
    m_an = -an.cantilever_moment(x, p, length)
    sigma_an = -m_an * problem.section.extreme_fibre / problem.section.second_moment
    plot_moment_stress(
        m_to_mm(fields.x), fields.moment, pa_to_mpa(fields.stress_top),
        m_to_mm(x), m_an, pa_to_mpa(sigma_an),
        RESULTS / "static_moment_stress.png",
    )


def main() -> dict[str, object]:
    """Run the static study, write CSV/plots, print a table, return key results."""
    problem, settings = default_problem(), StudySettings()
    n = problem.n_elements
    sol, fields = solve_cantilever_tip(problem, n)

    rows = tip_load_rows(problem, settings) + detail_rows(problem, sol, fields, n)
    rows += extra_case_rows(problem, settings, n)
    write_csv(RESULTS / "static_results.csv", CSV_HEADERS, [
        [r.case, r.quantity, r.unit, r.n_elements, f"{r.fe_value:.12g}",
         f"{r.analytical_value:.12g}", f"{r.abs_error:.3e}", f"{r.rel_error_pct:.3e}"]
        for r in rows
    ])
    make_plots(problem, sol, fields, n)

    peak, x_peak = max_stress(fields)
    print("\n=== STATIC VERIFICATION (FE vs Euler-Bernoulli theory) ===")
    print(format_table(
        ["case", "quantity", "n", "FE", "analytical", "unit", "rel. error [%]"],
        [[r.case, r.quantity, str(r.n_elements), f"{r.fe_value:.9g}",
          f"{r.analytical_value:.9g}", r.unit, f"{r.rel_error_pct:+.2e}"] for r in rows],
    ))
    yield_mpa = pa_to_mpa(problem.material.yield_strength)
    print(f"\nMax |sigma| (= von Mises, uniaxial) = {pa_to_mpa(peak):.3f} MPa at x = "
          f"{m_to_mm(x_peak):.1f} mm; safety factor vs yield ({yield_mpa:.0f} MPa) = "
          f"{problem.material.yield_strength / peak:.2f}")
    print("Wrote results/static_results.csv, static_deflection.png, static_moment_stress.png")
    return {"rows": rows, "max_stress_mpa": pa_to_mpa(peak), "max_stress_x_mm": m_to_mm(x_peak)}


if __name__ == "__main__":
    main()
