"""Post-processing: fields inside elements, bending moment and stress.

Role in the pipeline
--------------------
The solver only returns nodal values. Here the Hermite shape functions are
used again to interpolate *inside* each element, giving w(x), theta(x) and
curvature kappa(x) = w''(x). Then the classic stress-recovery chain is
applied:

    curvature  ->  bending moment M = E I kappa  ->  stress sigma = -M y / I

Sign convention: M is positive when sagging (w'' > 0). With y = +h/2 the
top fibre is in compression for a sagging moment, hence the minus sign.

Von Mises for the beam model equals |sigma| because the stress state is
uniaxial (sigma_xx only); transverse shear stress is neglected - for the
reference problem it is ~0.75 MPa versus 150 MPa bending stress.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .assembly import Mesh
from .elements import (
    FloatArray,
    shape_function_first_derivatives,
    shape_function_second_derivatives,
    shape_functions,
)


@dataclass(frozen=True)
class BeamFields:
    """Fields sampled along the beam (points repeat at element boundaries).

    Attributes:
        x: Sample positions [m].
        deflection: w(x) [m].
        rotation: theta(x) [rad].
        curvature: kappa(x) = w''(x) [1/m].
        moment: Bending moment M(x), sagging positive [N*m].
        stress_top: Bending stress at y = +h/2 [Pa].
        stress_bottom: Bending stress at y = -h/2 [Pa].
        von_mises: Uniaxial von Mises stress max(|top|, |bottom|) [Pa].
    """

    x: FloatArray
    deflection: FloatArray
    rotation: FloatArray
    curvature: FloatArray
    moment: FloatArray
    stress_top: FloatArray
    stress_bottom: FloatArray
    von_mises: FloatArray


def bending_stress(moment: FloatArray, y: float, second_moment: float) -> FloatArray:
    """Flexure formula sigma = -M y / I.

    Args:
        moment: Bending moment M, sagging positive [N*m].
        y: Distance from neutral axis, positive +y [m].
        second_moment: I [m^4].

    Returns:
        Normal stress sigma_xx [Pa] (tension positive).
    """
    return -np.asarray(moment) * y / second_moment


def evaluate_fields(
    mesh: Mesh,
    displacements: FloatArray,
    youngs_modulus: float,
    second_moment: float,
    depth: float,
    points_per_element: int = 11,
) -> BeamFields:
    """Interpolate deflection, rotation, curvature, moment and stress.

    Args:
        mesh: Beam mesh.
        displacements: Full DOF vector [m, rad].
        youngs_modulus: E [Pa].
        second_moment: I [m^4].
        depth: Section depth h [m] (extreme fibres at +-h/2).
        points_per_element: Samples per element including both ends [-].

    Returns:
        ``BeamFields`` sampled at ``points_per_element`` points per element.
    """
    xi = np.linspace(0.0, 1.0, points_per_element)
    xs, ws, ths, ks = [], [], [], []
    for e in range(mesh.n_elements):
        le = mesh.element_length(e)
        d = displacements[mesh.element_dofs(e)]
        x0 = mesh.nodes[mesh.connectivity[e, 0]]
        xs.append(x0 + xi * le)
        ws.append(d @ shape_functions(xi, le))
        ths.append(d @ shape_function_first_derivatives(xi, le))
        ks.append(d @ shape_function_second_derivatives(xi, le))
    curvature = np.concatenate(ks)
    moment = youngs_modulus * second_moment * curvature
    top = bending_stress(moment, depth / 2.0, second_moment)
    bottom = bending_stress(moment, -depth / 2.0, second_moment)
    return BeamFields(
        x=np.concatenate(xs),
        deflection=np.concatenate(ws),
        rotation=np.concatenate(ths),
        curvature=curvature,
        moment=moment,
        stress_top=top,
        stress_bottom=bottom,
        von_mises=np.maximum(np.abs(top), np.abs(bottom)),
    )


def max_stress(fields: BeamFields) -> tuple[float, float]:
    """Maximum von Mises (= |sigma|) stress and where it occurs.

    Returns:
        (max stress [Pa], x location [m]). Ties resolve to the first point.
    """
    idx = int(np.argmax(fields.von_mises))
    return float(fields.von_mises[idx]), float(fields.x[idx])


def value_at(fields_x: FloatArray, values: FloatArray, x: float) -> float:
    """Return the sampled value at position ``x`` (nearest sample point).

    Useful for reading off e.g. the mid-span stress; with an even number of
    elements x = L/2 is a node, so the value is exact.

    Args:
        fields_x: Sample positions [m].
        values: Field values at those positions [any unit].
        x: Query position [m].
    """
    return float(values[int(np.argmin(np.abs(fields_x - x)))])


# ---------------------------------------------------------------------------
# Tabular reporting (console tables and CSV files)
# ---------------------------------------------------------------------------


def relative_error_pct(fe_value: float, reference: float) -> float:
    """Signed relative error 100 * (fe - ref) / |ref| [%]."""
    return 100.0 * (fe_value - reference) / abs(reference)


def format_table(headers: list[str], rows: list[list[str]]) -> str:
    """Render rows of strings as an aligned plain-text table.

    Args:
        headers: Column titles.
        rows: Table body; every row has len(headers) string cells.

    Returns:
        Multi-line string ready for printing.
    """
    widths = [
        max(len(h), *(len(r[i]) for r in rows)) if rows else len(h)
        for i, h in enumerate(headers)
    ]
    line = "  ".join("-" * w for w in widths)
    out = ["  ".join(h.ljust(w) for h, w in zip(headers, widths)), line]
    out += ["  ".join(_align(c, w) for c, w in zip(r, widths)) for r in rows]
    return "\n".join(out)


def _align(cell: str, width: int) -> str:
    """Right-align numeric cells, left-align text cells."""
    try:
        float(cell)
    except ValueError:
        return cell.ljust(width)
    return cell.rjust(width)


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    """Render rows of strings as a GitHub-flavoured Markdown table."""
    out = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def write_csv(path: Path, headers: list[str], rows: list[list[object]]) -> Path:
    """Write a CSV file (creating parent folders) and return its path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        writer.writerows(rows)
    return path
