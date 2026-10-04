"""Compare hand-entered ANSYS Mechanical results with beam theory and the FE solver.

Pipeline step 4 of run_all.py (top rung of the verification ladder).
ANSYS numbers come ONLY from ansys/ansys_results.csv, which the user fills
in by hand following ansys/ANSYS_GUIDE.md. Nothing is estimated here: if
the file is missing, empty, or has no measured values, the script prints
"ANSYS data pending" and writes no output files.

Outputs (only when data exist):
    results/ansys_convergence.png      static mesh convergence (3 panels)
    results/ansys_vs_theory.csv        static ANSYS vs beam theory per mesh
    results/ansys_modal_comparison.csv bending-Y modes: ANSYS vs FE vs theory

Why ANSYS (3D solid) differs from the 1D beam model: shear deformation,
Poisson restraint at the fixed face, solid-element discretisation error,
and the stress singularity at the fixed support (see docs/CONCEPTS.md).
"""

from __future__ import annotations

import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
RESULTS = ROOT / "results"
ANSYS_CSV = ROOT / "ansys" / "ansys_results.csv"

from beamfea import analytical as an  # noqa: E402
from beamfea.assembly import assemble_beam  # noqa: E402
from beamfea.config import BeamProblem, StudySettings, default_problem, m_to_mm, pa_to_mpa  # noqa: E402
from beamfea.modal import solve_modal  # noqa: E402
from beamfea.plotting import plot_ansys_convergence  # noqa: E402
from beamfea.postprocess import format_table, relative_error_pct, write_csv  # noqa: E402
from beamfea.solver import cantilever_bcs  # noqa: E402

EXPECTED_HEADER = [
    "study", "element_size_mm", "n_elements", "n_nodes", "tip_deflection_mm",
    "midspan_top_stress_mpa", "max_equivalent_stress_root_mpa", "mode",
    "frequency_hz", "mode_type",
]
PENDING_MESSAGE = "ANSYS data pending: see ansys/ANSYS_GUIDE.md"


class AnsysDataError(ValueError):
    """Raised when ansys_results.csv exists but cannot be interpreted."""


@dataclass(frozen=True)
class StaticRecord:
    """One ANSYS static mesh (NaN = not recorded). Units: mm, MPa."""

    element_size_mm: float
    n_elements: float
    n_nodes: float
    tip_deflection_mm: float
    midspan_stress_mpa: float
    max_root_stress_mpa: float


@dataclass(frozen=True)
class ModalRecord:
    """One ANSYS mode identified as bending in Y. Units: Hz."""

    ansys_mode: int
    frequency_hz: float


@dataclass(frozen=True)
class TheoryValues:
    """Beam-theory reference values in reporting units (mm, MPa)."""

    tip_mm: float
    midspan_mpa: float
    root_mpa: float


def _number(cell: str | None) -> float:
    """Parse a CSV cell to float; blank -> NaN. Raises AnsysDataError on junk."""
    text = (cell or "").strip()
    if not text:
        return math.nan
    try:
        return float(text)
    except ValueError as exc:
        raise AnsysDataError(f"not a number: {text!r}") from exc


def load_ansys_csv(path: Path) -> tuple[list[StaticRecord], list[ModalRecord]]:
    """Read the user's ANSYS CSV, keeping only rows that contain measurements.

    Static rows need an element size and at least one measured quantity.
    Modal rows need ``mode_type == bending_y`` and a frequency; other mode
    types (bending_z, torsion, ...) are deliberately ignored.

    Args:
        path: CSV in the format of ansys/ansys_results_template.csv.

    Returns:
        (static records sorted coarse -> fine, bending-Y modal records sorted
        by frequency).
    """
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        header = [h.strip() for h in (reader.fieldnames or [])]
        if not header:
            return [], []
        if header != EXPECTED_HEADER:
            raise AnsysDataError(f"unexpected header {header}; expected {EXPECTED_HEADER}")
        static, modal = [], []
        for line_no, raw in enumerate(reader, start=2):
            row = {k.strip(): (v or "") for k, v in raw.items() if k}
            try:
                study = row["study"].strip().lower()
                if study == "static":
                    rec = StaticRecord(*(_number(row[c]) for c in EXPECTED_HEADER[1:7]))
                    measured = (rec.tip_deflection_mm, rec.midspan_stress_mpa, rec.max_root_stress_mpa)
                    if math.isfinite(rec.element_size_mm) and any(math.isfinite(v) for v in measured):
                        static.append(rec)
                elif study == "modal":
                    freq = _number(row["frequency_hz"])
                    if row["mode_type"].strip().lower() == "bending_y" and math.isfinite(freq):
                        mode = _number(row["mode"])
                        modal.append(ModalRecord(int(mode) if math.isfinite(mode) else 0, freq))
            except AnsysDataError as exc:
                raise AnsysDataError(f"line {line_no}: {exc}") from exc
    static.sort(key=lambda r: -r.element_size_mm)
    modal.sort(key=lambda r: r.frequency_hz)
    return static, modal


def theory_values(problem: BeamProblem) -> TheoryValues:
    """Beam-theory tip deflection, mid-span and root top-fibre stress."""
    p, length = problem.tip_load, problem.length
    c, i_z = problem.section.extreme_fibre, problem.section.second_moment
    return TheoryValues(
        tip_mm=m_to_mm(an.cantilever_tip_deflection(p, length, problem.flexural_rigidity)),
        midspan_mpa=pa_to_mpa(an.bending_stress_magnitude(float(an.cantilever_moment(length / 2, p, length)), c, i_z)),
        root_mpa=pa_to_mpa(an.bending_stress_magnitude(float(an.cantilever_moment(0.0, p, length)), c, i_z)),
    )


def _pct(value: float, reference: float) -> float:
    """Relative difference [%] or NaN if the value was not recorded."""
    return relative_error_pct(value, reference) if math.isfinite(value) else math.nan


def _fmt(value: float, spec: str = ".4f") -> str:
    """Format a number, or an empty cell for NaN (not recorded)."""
    return format(value, spec) if math.isfinite(value) else ""


def _finest(static: list[StaticRecord], diffs: list[float], attr: str) -> tuple[float, float]:
    """(element size [mm], % diff) of the finest mesh where ``attr`` was recorded."""
    for rec, diff in zip(reversed(static), reversed(diffs)):
        if math.isfinite(getattr(rec, attr)):
            return rec.element_size_mm, diff
    return math.nan, math.nan


def compare_static(static: list[StaticRecord], theory: TheoryValues) -> dict[str, object]:
    """Write ansys_vs_theory.csv and ansys_convergence.png; return a summary."""
    rows, table, tip_diffs, mid_diffs = [], [], [], []
    for r in static:
        tip, mid, root = abs(r.tip_deflection_mm), abs(r.midspan_stress_mpa), abs(r.max_root_stress_mpa)
        diffs = (_pct(tip, theory.tip_mm), _pct(mid, theory.midspan_mpa), _pct(root, theory.root_mpa))
        tip_diffs.append(diffs[0])
        mid_diffs.append(diffs[1])
        rows.append([_fmt(r.element_size_mm, "g"), _fmt(r.n_elements, ".0f"), _fmt(r.n_nodes, ".0f"),
                     _fmt(tip), f"{theory.tip_mm:.4f}", _fmt(diffs[0], ".3f"),
                     _fmt(mid), f"{theory.midspan_mpa:.4f}", _fmt(diffs[1], ".3f"),
                     _fmt(root), f"{theory.root_mpa:.4f}", _fmt(diffs[2], ".3f")])
        table.append([c or "n/r" for c in (rows[-1][0], rows[-1][3], rows[-1][5], rows[-1][6],
                                             rows[-1][8], rows[-1][9], rows[-1][11])])
    write_csv(RESULTS / "ansys_vs_theory.csv", [
        "element_size_mm", "n_elements", "n_nodes",
        "tip_deflection_ansys_mm", "tip_deflection_theory_mm", "tip_deflection_diff_pct",
        "midspan_stress_ansys_mpa", "midspan_stress_theory_mpa", "midspan_stress_diff_pct",
        "max_root_stress_ansys_mpa", "root_stress_beam_theory_mpa", "max_root_stress_diff_pct",
    ], rows)
    plot_ansys_convergence(
        np.array([r.element_size_mm for r in static]),
        np.abs([r.tip_deflection_mm for r in static]),
        np.abs([r.midspan_stress_mpa for r in static]),
        np.abs([r.max_root_stress_mpa for r in static]),
        theory.tip_mm, theory.midspan_mpa, theory.root_mpa,
        RESULTS / "ansys_convergence.png",
    )
    print("\n=== ANSYS STATIC vs BEAM THEORY (values as entered in ansys_results.csv) ===")
    print(format_table(["size [mm]", "tip [mm]", "diff [%]", "mid σ [MPa]", "diff [%]",
                        "root σeq [MPa]", "vs 150 [%]"], table))
    finest = {"tip": _finest(static, tip_diffs, "tip_deflection_mm"),
              "mid": _finest(static, mid_diffs, "midspan_stress_mpa")}
    for key, label in (("tip", "tip deflection"), ("mid", "mid-span top stress")):
        size, diff = finest[key]
        if math.isfinite(size):
            print(f"Finest mesh with {label} ({size:g} mm): ANSYS vs beam theory = {diff:+.3f} %")
    print("(n/r = not recorded.) The root peak stress is a singularity and is not expected to converge.")
    return {"rows": rows, "finest": finest}


def compare_modal(modal: list[ModalRecord], problem: BeamProblem) -> dict[str, object]:
    """Write ansys_modal_comparison.csv for bending-Y modes; return a summary."""
    settings = StudySettings()
    n_modes = len(modal)
    _, k, m = assemble_beam(problem.length, settings.modal_elements,
                            problem.flexural_rigidity, problem.mass_per_length)
    f_fe = solve_modal(k, m, cantilever_bcs(), n_modes).frequencies_hz
    f_an = an.cantilever_natural_frequencies(problem.flexural_rigidity, problem.mass_per_length,
                                             problem.length, n_modes)
    rows = []
    for j, rec in enumerate(modal):
        rows.append([j + 1, rec.ansys_mode, f"{rec.frequency_hz:.4f}", f"{f_fe[j]:.4f}", f"{f_an[j]:.4f}",
                     f"{relative_error_pct(rec.frequency_hz, f_an[j]):.3f}",
                     f"{relative_error_pct(rec.frequency_hz, f_fe[j]):.3f}"])
    headers = ["bending_y_mode", "ansys_mode_number", "f_ansys_hz", "f_python_fe_hz",
               "f_analytical_hz", "ansys_vs_analytical_pct", "ansys_vs_python_fe_pct"]
    write_csv(RESULTS / "ansys_modal_comparison.csv", headers, rows)
    print("\n=== ANSYS MODAL (bending-Y) vs PYTHON FE vs THEORY ===")
    print(format_table(["#", "ANSYS mode", "ANSYS [Hz]", "FE [Hz]", "theory [Hz]",
                        "vs theory [%]", "vs FE [%]"], [[str(c) for c in r] for r in rows]))
    return {"headers": headers, "rows": rows}


def main() -> dict[str, object]:
    """Compare ANSYS data if present; otherwise report 'pending'. Never fails."""
    pending: dict[str, object] = {"status": "pending", "static": None, "modal": None}
    if not ANSYS_CSV.exists():
        print(f"\n{PENDING_MESSAGE} (ansys/ansys_results.csv not found)")
        return pending
    try:
        static, modal = load_ansys_csv(ANSYS_CSV)
    except AnsysDataError as exc:
        print(f"\nCould not read ansys/ansys_results.csv: {exc}")
        print(PENDING_MESSAGE)
        return {**pending, "error": str(exc)}
    if not static and not modal:
        print(f"\n{PENDING_MESSAGE} (ansys_results.csv contains no measured values)")
        return pending

    problem = default_problem()
    report: dict[str, object] = {"status": "available", "static": None, "modal": None}
    if static:
        report["static"] = compare_static(static, theory_values(problem))
    else:
        print("\nANSYS static data pending (no static rows with values).")
    if modal:
        report["modal"] = compare_modal(modal, problem)
    else:
        print("\nANSYS modal data pending (no rows with mode_type = bending_y and a frequency).")
    return report


if __name__ == "__main__":
    main()
