"""Modal h-refinement study: frequency error vs mesh size.

Pipeline step 3 of run_all.py. Static results are exact for the tip load
(see docs/CONCEPTS.md), so the meaningful convergence study in Python is the
modal one. For each mesh the first 3 frequencies are compared with theory
and the observed convergence rate of mode 1 is measured by a least-squares
fit of log(error) against log(h), with h = L / n.

Outputs: results/modal_convergence.csv, results/modal_convergence.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
RESULTS = ROOT / "results"

from beamfea import analytical as an  # noqa: E402
from beamfea.assembly import assemble_beam  # noqa: E402
from beamfea.config import StudySettings, default_problem, m_to_mm  # noqa: E402
from beamfea.modal import solve_modal  # noqa: E402
from beamfea.plotting import plot_modal_convergence  # noqa: E402
from beamfea.postprocess import format_table, write_csv  # noqa: E402
from beamfea.solver import cantilever_bcs  # noqa: E402

# Fit only the asymptotic range: the 1-element mesh is pre-asymptotic.
MIN_ELEMENTS_FOR_FIT = 2


def asymptotic_mask(n: np.ndarray, err: np.ndarray) -> np.ndarray:
    """Select meshes for the slope fit [bool array].

    Starts at MIN_ELEMENTS_FOR_FIT and stops at the first mesh whose error
    does *not* decrease: from there on the error is dominated by
    floating-point round-off in the eigen-solver (relative error ~1e-8),
    not by discretisation, and including it would corrupt the slope.
    """
    use = np.zeros(n.size, dtype=bool)
    previous = np.inf
    for i in range(n.size):
        if n[i] < MIN_ELEMENTS_FOR_FIT or not np.isfinite(err[i]):
            continue
        if not (0.0 < err[i] < previous):
            break
        use[i], previous = True, err[i]
    return use


def fit_slope(h: np.ndarray, err: np.ndarray) -> float:
    """Least-squares slope p of log(err) = p log(h) + c [-]."""
    p, _ = np.polyfit(np.log(h), np.log(err), 1)
    return float(p)


def main() -> dict[str, object]:
    """Run the convergence study, write CSV/plot, print table and slope."""
    problem, settings = default_problem(), StudySettings()
    meshes, n_modes = settings.convergence_meshes, settings.convergence_n_modes
    f_exact = an.cantilever_natural_frequencies(
        problem.flexural_rigidity, problem.mass_per_length, problem.length, n_modes)

    freqs = np.full((len(meshes), n_modes), np.nan)
    for i, n in enumerate(meshes):
        _, k, m = assemble_beam(problem.length, n, problem.flexural_rigidity, problem.mass_per_length)
        f = solve_modal(k, m, cantilever_bcs()).frequencies_hz[:n_modes]
        freqs[i, : f.size] = f  # 1 element has only 2 free DOFs -> 2 modes
    rel_err = (freqs - f_exact) / f_exact

    n_arr = np.array(meshes, dtype=float)
    h = problem.length / n_arr
    use = asymptotic_mask(n_arr, rel_err[:, 0])
    slope_h = fit_slope(h[use], rel_err[use, 0])

    rows = []
    for i, n in enumerate(meshes):
        row: list[object] = [n, f"{m_to_mm(h[i]):.4f}"]
        for j in range(n_modes):
            row += (["n/a", "n/a"] if np.isnan(freqs[i, j])
                    else [f"{freqs[i, j]:.8f}", f"{100 * rel_err[i, j]:.6e}"])
        rows.append(row)
    headers = ["n_elements", "element_size_mm"]
    for j in range(n_modes):
        headers += [f"f{j + 1}_fe_hz", f"f{j + 1}_error_pct"]
    write_csv(RESULTS / "modal_convergence.csv", headers, rows)

    print("\n=== MODAL h-REFINEMENT STUDY ===")
    print("Exact: " + ", ".join(f"f{j + 1} = {f:.4f} Hz" for j, f in enumerate(f_exact)))
    print(format_table(
        ["n", "h [mm]"] + [f"mode {j + 1} error [%]" for j in range(n_modes)],
        [[str(n), f"{m_to_mm(h[i]):.3f}"] +
         ["n/a" if np.isnan(rel_err[i, j]) else f"{100 * rel_err[i, j]:+.3e}" for j in range(n_modes)]
         for i, n in enumerate(meshes)],
    ))
    used = [int(v) for v in n_arr[use]]
    excluded = [int(v) for v in n_arr[~use & (n_arr >= MIN_ELEMENTS_FOR_FIT)]]
    plot_modal_convergence(meshes, rel_err, -slope_h, RESULTS / "modal_convergence.png", excluded)
    print(f"Fitted convergence rate, mode 1: error ~ h^{slope_h:.3f} "
          f"(fit over n = {used}; theory for cubic Hermite + consistent mass: ~h^4)")
    if excluded:
        print(f"Excluded from fit (round-off floor reached, error stopped decreasing): n = {excluded}")
    all_above = bool(np.all(rel_err[np.isfinite(rel_err)] > 0))
    print("FE frequencies converge from above:", all_above)
    print("Wrote results/modal_convergence.csv, modal_convergence.png")
    return {"slope_h": slope_h, "fit_meshes": used, "excluded_meshes": excluded, "converges_from_above": all_above,
            "meshes": list(meshes), "rel_err": rel_err}


if __name__ == "__main__":
    main()
