"""All plotting functions (matplotlib, headless, PNG at 150 dpi).

Role in the pipeline
--------------------
Turns arrays produced by the solver/post-processor into figures in
``results/``. The plotting functions take plain arrays in *reporting*
units (mm, MPa, Hz) so that nothing here knows about FE internals.
Conventions: FE results are drawn as markers, analytical references as
solid lines, so "dots on the line" means agreement.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")  # headless backend: works in CI and over SSH

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401
import numpy as np  # noqa: E402

from .elements import FloatArray  # noqa: E402

DPI = 150
# Fixed categorical order (never cycled): blue, orange, aqua, yellow.
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
REFERENCE = "#52514e"  # neutral ink for analytical / reference lines

plt.rcParams.update(
    {
        "axes.grid": True,
        "grid.color": "#dddddd",
        "grid.linewidth": 0.6,
        "axes.edgecolor": "#888888",
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 9,
        "lines.linewidth": 2.0,
        "lines.markersize": 6,
        "font.family": "DejaVu Sans",
    }
)


def _save(fig: plt.Figure, path: Path) -> Path:
    """Save ``fig`` to ``path`` as PNG at the project DPI and close it."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def plot_static_deflection(
    x_nodes_mm: FloatArray,
    w_nodes_mm: FloatArray,
    x_curve_mm: FloatArray,
    w_curve_mm: FloatArray,
    path: Path,
    title: str,
) -> Path:
    """FE nodal deflections (markers) against the analytical curve (line).

    Args:
        x_nodes_mm: Node positions [mm].
        w_nodes_mm: FE nodal deflections, +y positive [mm].
        x_curve_mm: Positions for the analytical curve [mm].
        w_curve_mm: Analytical deflection, +y positive [mm].
        path: Output PNG path.
        title: Figure title.
    """
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.plot(x_curve_mm, w_curve_mm, color=REFERENCE, label="Analytical  w(x)")
    ax.plot(
        x_nodes_mm, w_nodes_mm, "o", color=SERIES[0], mfc="white", mew=2,
        label="FE nodes",
    )
    ax.set_title(title)
    ax.set_xlabel("Position along beam x [mm]")
    ax.set_ylabel("Deflection w [mm]  (+y up)")
    ax.legend(loc="lower left")
    return _save(fig, path)


def plot_moment_stress(
    x_fe_mm: FloatArray,
    moment_fe: FloatArray,
    stress_fe_mpa: FloatArray,
    x_an_mm: FloatArray,
    moment_an: FloatArray,
    stress_an_mpa: FloatArray,
    path: Path,
) -> Path:
    """Bending moment and top-fibre stress, FE vs analytical (two panels).

    Args:
        x_fe_mm: FE sample positions [mm].
        moment_fe: FE bending moment, sagging positive [N*m].
        stress_fe_mpa: FE top-fibre bending stress [MPa].
        x_an_mm: Analytical sample positions [mm].
        moment_an: Analytical moment, same sign convention [N*m].
        stress_an_mpa: Analytical top-fibre stress [MPa].
        path: Output PNG path.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 6.4), sharex=True)
    step = max(1, x_fe_mm.size // 40)  # thin the markers so the line stays visible
    for ax, fe, an, ylabel, title in (
        (ax1, moment_fe, moment_an, "Bending moment M [N·m]",
         "Bending moment (sagging positive)"),
        (ax2, stress_fe_mpa, stress_an_mpa, "Top-fibre stress σ [MPa]",
         "Bending stress at top fibre (y = +h/2, tension positive)"),
    ):
        ax.plot(x_an_mm, an, color=REFERENCE, label="Analytical")
        ax.plot(x_fe_mm[::step], fe[::step], "o", color=SERIES[1], mfc="white",
                mew=1.6, label="FE (interpolated in elements)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(loc="best")
    ax2.set_xlabel("Position along beam x [mm]")
    return _save(fig, path)


def plot_mode_shapes(
    x_mm: FloatArray,
    shapes: Sequence[FloatArray],
    frequencies_hz: Sequence[float],
    path: Path,
) -> Path:
    """First few mode shapes, each normalised to unit peak |w|.

    Args:
        x_mm: Sample positions [mm].
        shapes: One deflection array per mode [arbitrary units].
        frequencies_hz: Natural frequency of each mode [Hz].
        path: Output PNG path.
    """
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    for k, (shape, f) in enumerate(zip(shapes, frequencies_hz)):
        normalised = shape / np.max(np.abs(shape))
        ax.plot(x_mm, normalised, color=SERIES[k % len(SERIES)],
                label=f"Mode {k + 1}: {f:.1f} Hz")
    ax.axhline(0.0, color="#888888", linewidth=0.8)
    ax.set_title("Cantilever bending mode shapes (FE, 20 elements)")
    ax.set_xlabel("Position along beam x [mm]")
    ax.set_ylabel("Normalised deflection [-]")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    return _save(fig, path)


def plot_modal_convergence(
    n_elements: Sequence[int],
    relative_errors: FloatArray,
    fitted_slope_vs_n: float,
    path: Path,
    excluded_n: Sequence[int] = (),
) -> Path:
    """Relative frequency error vs number of elements on log-log axes.

    Args:
        n_elements: Element counts [-].
        relative_errors: Array (n_meshes, n_modes) of (f_FE - f_exact)/f_exact
            [-]; NaN where the mesh has too few DOFs for that mode.
        fitted_slope_vs_n: Fitted slope of log(error) vs log(n) for mode 1 [-].
        path: Output PNG path.
        excluded_n: Meshes left out of the mode-1 fit (round-off floor) [-].
    """
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    n = np.asarray(n_elements, dtype=float)
    for k in range(relative_errors.shape[1]):
        ok = np.isfinite(relative_errors[:, k]) & (relative_errors[:, k] > 0)
        ax.loglog(n[ok], 100.0 * relative_errors[ok, k], "o-",
                  color=SERIES[k % len(SERIES)], label=f"Mode {k + 1}")
    for n_ex in excluded_n:
        i = list(n_elements).index(n_ex)
        ax.annotate("round-off floor\n(excluded from fit)",
                    xy=(n_ex, 100.0 * relative_errors[i, 0]), xytext=(-150, -6),
                    textcoords="offset points", fontsize=8, color=REFERENCE,
                    arrowprops={"arrowstyle": "->", "color": REFERENCE})
    ax.set_title(
        "Modal h-refinement (consistent mass)\n"
        f"mode 1 fitted slope = {fitted_slope_vs_n:.2f} vs n, i.e. error ∝ h^{-fitted_slope_vs_n:.2f}"
    )
    ax.set_xlabel("Number of elements n [-]")
    ax.set_ylabel("Relative frequency error (FE − exact)/exact [%]")
    ax.legend(loc="upper right")
    return _save(fig, path)


def plot_ansys_convergence(
    element_size_mm: FloatArray,
    tip_deflection_mm: FloatArray,
    midspan_stress_mpa: FloatArray,
    max_root_stress_mpa: FloatArray,
    analytical_tip_mm: float,
    analytical_midspan_mpa: float,
    beam_root_stress_mpa: float,
    path: Path,
) -> Path:
    """ANSYS mesh-convergence plots with beam-theory reference lines.

    Three panels: tip deflection and mid-span top stress (which should
    converge) and peak equivalent stress at the fixed support (which should
    keep rising - a stress singularity). NaN entries are skipped.

    Args:
        element_size_mm: ANSYS element sizes [mm].
        tip_deflection_mm: |tip deflection| from ANSYS [mm].
        midspan_stress_mpa: Mid-span top-surface normal stress [MPa].
        max_root_stress_mpa: Max equivalent stress (at the support) [MPa].
        analytical_tip_mm: Beam-theory tip deflection [mm].
        analytical_midspan_mpa: Beam-theory mid-span top stress [MPa].
        beam_root_stress_mpa: Beam-theory root stress [MPa].
        path: Output PNG path.
    """
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    panels = (
        (tip_deflection_mm, analytical_tip_mm, "Tip deflection |w| [mm]",
         "Tip deflection (should converge)"),
        (midspan_stress_mpa, analytical_midspan_mpa, "Normal stress σx [MPa]",
         "Mid-span top stress (should converge)"),
        (max_root_stress_mpa, beam_root_stress_mpa, "Equivalent stress [MPa]",
         "Peak stress at fixed support\n(singularity: does NOT converge)"),
    )
    for k, (ax, (values, ref, ylabel, title)) in enumerate(zip(axes, panels)):
        ok = np.isfinite(values)
        ax.semilogx(element_size_mm[ok], values[ok], "o-", color=SERIES[k],
                    label="ANSYS")
        ax.axhline(ref, color=REFERENCE, linestyle="--", linewidth=1.5,
                   label="Beam theory")
        # Label the actual mesh sizes instead of default log-scale ticks.
        ax.set_xticks(element_size_mm)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlim(element_size_mm.max() * 1.25, element_size_mm.min() / 1.25)  # finer → right
        ax.set_xlabel("Element size [mm] (finer →)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(loc="best")
    return _save(fig, path)
