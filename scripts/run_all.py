"""Run the complete pipeline and write results/summary.md.

Order: hand calculations -> run_static -> run_modal -> run_convergence ->
compare_ansys. Every number in summary.md is produced by this run; ANSYS
numbers appear only if ansys/ansys_results.csv has been filled in by hand.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(SCRIPTS))
RESULTS = ROOT / "results"

import compare_ansys  # noqa: E402
import run_convergence  # noqa: E402
import run_modal  # noqa: E402
import run_static  # noqa: E402
from beamfea import analytical as an  # noqa: E402
from beamfea.config import BeamProblem, default_problem, m_to_mm, pa_to_mpa  # noqa: E402
from beamfea.postprocess import format_table, markdown_table  # noqa: E402


def hand_calculations(problem: BeamProblem) -> list[list[str]]:
    """Closed-form sanity values for the reference cantilever (reporting units)."""
    p, length, ei = problem.tip_load, problem.length, problem.flexural_rigidity
    sec = problem.section
    root_m = float(an.cantilever_moment(0.0, p, length))
    mid_m = float(an.cantilever_moment(length / 2.0, p, length))
    freqs = an.cantilever_natural_frequencies(ei, problem.mass_per_length, length, 3)
    rows = [
        ["Second moment of area I = b h³/12", f"{sec.second_moment * 1e12:.2f}", "mm⁴"],
        ["Tip deflection P L³/(3EI)", f"{m_to_mm(an.cantilever_tip_deflection(p, length, ei)):.4f}", "mm"],
        ["Tip rotation P L²/(2EI)", f"{an.cantilever_tip_rotation(p, length, ei):.5f}", "rad"],
        ["Root bending moment P L", f"{root_m:.2f}", "N·m"],
        ["Root max bending stress M c / I", f"{pa_to_mpa(an.bending_stress_magnitude(root_m, sec.extreme_fibre, sec.second_moment)):.2f}", "MPa"],
        ["Mid-span top-fibre stress", f"{pa_to_mpa(an.bending_stress_magnitude(mid_m, sec.extreme_fibre, sec.second_moment)):.2f}", "MPa"],
        ["Max transverse shear stress 1.5 V / A", f"{pa_to_mpa(an.rectangular_max_shear_stress(p, sec.area)):.3f}", "MPa"],
    ]
    rows += [[f"Bending natural frequency f{j + 1}", f"{f:.3f}", "Hz"] for j, f in enumerate(freqs)]
    return rows


def ansys_section(report: dict[str, object]) -> list[str]:
    """Markdown lines for the ANSYS part of the summary."""
    lines = ["## 5. ANSYS Mechanical comparison", ""]
    if report["status"] != "available":
        lines += ["**ANSYS data pending.** Follow `ansys/ANSYS_GUIDE.md`, fill in "
                  "`ansys/ansys_results.csv`, then re-run `python scripts/run_all.py`.", ""]
        return lines
    static, modal = report["static"], report["modal"]
    if static:
        hdr = ["size [mm]", "elements", "nodes", "tip ANSYS [mm]", "tip theory [mm]", "diff [%]",
               "mid σ ANSYS [MPa]", "mid σ theory [MPa]", "diff [%]", "root σeq ANSYS [MPa]",
               "root σ beam [MPa]", "diff [%]"]
        lines += ["### Static mesh convergence (values entered by the user from ANSYS)", "",
                  markdown_table(hdr, [[str(c) or "n/r" for c in r] for r in static["rows"]]), ""]
        for key, label in (("tip", "Tip deflection"), ("mid", "Mid-span top stress")):
            size, diff = static["finest"][key]
            if np.isfinite(size):
                lines += [f"- {label}, finest mesh ({size:g} mm): ANSYS differs from beam "
                          f"theory by {diff:+.3f} %."]
        lines += ["", "The root peak stress is a stress singularity and is not a converged "
                  "quantity (n/r = not recorded).", "",
                  "![ANSYS convergence](ansys_convergence.png)", ""]
    else:
        lines += ["ANSYS static data pending.", ""]
    if modal:
        lines += ["### Bending-Y modes", "",
                  markdown_table(modal["headers"], [[str(c) for c in r] for r in modal["rows"]]), ""]
    else:
        lines += ["ANSYS modal data pending.", ""]
    return lines


def write_summary(problem: BeamProblem, hand: list[list[str]], static: dict, modal: dict,
                  conv: dict, ansys: dict) -> Path:
    """Assemble results/summary.md from the script results."""
    rows = static["rows"]
    static_md = markdown_table(
        ["case", "quantity", "n", "FE", "analytical", "unit", "rel. error [%]"],
        [[r.case, r.quantity, str(r.n_elements), f"{r.fe_value:.9g}", f"{r.analytical_value:.9g}",
          r.unit, f"{r.rel_error_pct:+.2e}"] for r in rows])
    modal_md = markdown_table(
        ["mode", "f FE [Hz]", "f analytical [Hz]", "error [%]"],
        [[str(j + 1), f"{a:.4f}", f"{b:.4f}", f"{e:+.2e}"]
         for j, (a, b, e) in enumerate(zip(modal["f_fe"], modal["f_an"], modal["errors_pct"]))])
    rel = conv["rel_err"]
    conv_md = markdown_table(
        ["n elements"] + [f"mode {j + 1} error [%]" for j in range(rel.shape[1])],
        [[str(n)] + ["n/a" if np.isnan(v) else f"{100 * v:+.3e}" for v in rel[i]]
         for i, n in enumerate(conv["meshes"])])
    excluded = conv["excluded_meshes"]
    lines = [
        "# beam-fea results summary", "",
        "_Generated by `python scripts/run_all.py`. Do not edit by hand._", "",
        f"Reference problem: steel cantilever, L = {m_to_mm(problem.length):.0f} mm, "
        f"b = {m_to_mm(problem.section.width):.0f} mm, h = {m_to_mm(problem.section.depth):.0f} mm, "
        f"tip load P = {problem.tip_load:.0f} N downward.", "",
        "## 1. Hand calculations (closed-form Euler-Bernoulli)", "",
        markdown_table(["quantity", "value", "unit"], hand), "",
        "## 2. Static verification: FE vs analytical", "",
        static_md, "",
        f"Maximum |σ| (von Mises, uniaxial) = {static['max_stress_mpa']:.3f} MPa at "
        f"x = {static['max_stress_x_mm']:.1f} mm.", "",
        "The tip-load deflection is exact for every mesh (including 1 element) because the "
        "exact solution is a cubic, which Hermite elements reproduce; the tiny non-zero errors "
        "are floating-point round-off.", "",
        "![Static deflection](static_deflection.png)", "",
        "![Moment and stress](static_moment_stress.png)", "",
        f"## 3. Modal verification ({modal['n_elements']} elements)", "",
        modal_md, "",
        "![Mode shapes](mode_shapes.png)", "",
        "## 4. Modal h-refinement", "",
        conv_md, "",
        f"Fitted convergence rate for mode 1: **error ∝ h^{conv['slope_h']:.3f}** "
        f"(fit over n = {conv['fit_meshes']}).",
        (f" Excluded n = {excluded}: error stopped decreasing (floating-point round-off floor)."
         if excluded else ""), "",
        f"All FE frequencies above the exact values (converge from above): "
        f"{conv['converges_from_above']}.", "",
        "![Modal convergence](modal_convergence.png)", "",
    ]
    lines += ansys_section(ansys)
    path = RESULTS / "summary.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main() -> None:
    """Run every study in order and write the summary."""
    problem = default_problem()
    hand = hand_calculations(problem)
    print("=== HAND CALCULATIONS (closed form) ===")
    print(format_table(["quantity", "value", "unit"], hand))
    static = run_static.main()
    modal = run_modal.main()
    conv = run_convergence.main()
    ansys = compare_ansys.main()
    path = write_summary(problem, hand, static, modal, conv, ansys)
    print(f"\nWrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
