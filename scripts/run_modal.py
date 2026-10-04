"""Modal verification: FE natural frequencies vs cantilever beam theory.

Pipeline step 2 of run_all.py. Solves the generalized eigenproblem for the
reference cantilever with 20 elements and compares the first 5 bending
frequencies with the closed-form values.

Outputs: results/modal_results.csv, results/mode_shapes.png
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
RESULTS = ROOT / "results"

from beamfea import analytical as an  # noqa: E402
from beamfea.assembly import assemble_beam  # noqa: E402
from beamfea.config import StudySettings, default_problem, m_to_mm  # noqa: E402
from beamfea.modal import solve_modal  # noqa: E402
from beamfea.plotting import plot_mode_shapes  # noqa: E402
from beamfea.postprocess import evaluate_fields, format_table, relative_error_pct, write_csv  # noqa: E402
from beamfea.solver import cantilever_bcs  # noqa: E402

N_PLOTTED_MODES = 4


def main() -> dict[str, object]:
    """Run the modal study, write CSV/plot, print a table, return key results."""
    problem, settings = default_problem(), StudySettings()
    n, n_modes = settings.modal_elements, settings.modal_n_modes
    mesh, k, m = assemble_beam(problem.length, n, problem.flexural_rigidity, problem.mass_per_length)
    modal = solve_modal(k, m, cantilever_bcs(), n_modes)

    f_fe = modal.frequencies_hz
    f_an = an.cantilever_natural_frequencies(
        problem.flexural_rigidity, problem.mass_per_length, problem.length, n_modes)
    errors = [relative_error_pct(a, b) for a, b in zip(f_fe, f_an)]

    write_csv(RESULTS / "modal_results.csv",
              ["mode", "n_elements", "f_fe_hz", "f_analytical_hz", "error_pct"],
              [[k_ + 1, n, f"{a:.6f}", f"{b:.6f}", f"{e:.6e}"]
               for k_, (a, b, e) in enumerate(zip(f_fe, f_an, errors))])

    # Interpolate mode shapes inside elements with the Hermite functions
    # (only the deflection field is needed, so E, I, h are dummies here).
    shapes, x = [], None
    for j in range(N_PLOTTED_MODES):
        fields = evaluate_fields(mesh, modal.mode_shapes[:, j], 1.0, 1.0, 1.0,
                                 settings.points_per_element)
        shapes.append(fields.deflection)
        x = fields.x
    plot_mode_shapes(m_to_mm(x), shapes, f_fe[:N_PLOTTED_MODES], RESULTS / "mode_shapes.png")

    print(f"\n=== MODAL VERIFICATION ({n} elements, consistent mass) ===")
    print(format_table(
        ["mode", "f_FE [Hz]", "f_analytical [Hz]", "error [%]"],
        [[str(i + 1), f"{a:.4f}", f"{b:.4f}", f"{e:+.2e}"]
         for i, (a, b, e) in enumerate(zip(f_fe, f_an, errors))],
    ))
    print("All FE frequencies above theory (upper bound):", all(e > 0 for e in errors))
    print("Wrote results/modal_results.csv, mode_shapes.png")
    return {"n_elements": n, "f_fe": list(f_fe), "f_an": list(f_an), "errors_pct": errors}


if __name__ == "__main__":
    main()
