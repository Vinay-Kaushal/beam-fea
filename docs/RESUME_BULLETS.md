# Resume bullets: placeholders and where their numbers come from

Never type a number on your resume from memory. Each placeholder below maps to a generated
file and column. Regenerate with `python scripts/run_all.py` before copying.

## Placeholder map

| Placeholder | Meaning | Source file | Column / how to read | Status |
|---|---|---|---|---|
| `{TIP_ERR_PY}` | Python FE tip-deflection error vs theory | `results/static_results.csv` | `rel_error_pct`, rows `case = cantilever_tip_load`, `quantity = tip_deflection` (take the largest magnitude across n) | Available now |
| `{N_MESH_PY}` | Number of meshes in the Python modal convergence study | `results/modal_convergence.csv` | number of data rows | Available now |
| `{SLOPE}` | Fitted modal convergence rate (error ∝ h^SLOPE) | `results/summary.md` §4 (also printed by `run_convergence.py`) | "Fitted convergence rate for mode 1" | Available now |
| `{N_MODES_PY}` | Modes compared, Python vs theory | `results/modal_results.csv` | number of data rows | Available now |
| `{MAX_FREQ_ERR_PY}` | Max frequency error, Python FE (20 elements) vs theory | `results/modal_results.csv` | max of `error_pct` | Available now |
| `{TIP_ERR_ANSYS}` | ANSYS tip-deflection error vs theory at the finest mesh | `results/ansys_vs_theory.csv` | `tip_deflection_diff_pct`, last row with a value | **Pending ANSYS CSV** |
| `{MID_ERR_ANSYS}` | ANSYS mid-span stress error vs theory at the finest mesh | `results/ansys_vs_theory.csv` | `midspan_stress_diff_pct`, last row with a value | **Pending ANSYS CSV** |
| `{N_MESH_ANSYS}` | Number of ANSYS mesh sizes actually solved | `results/ansys_vs_theory.csv` | rows with a tip deflection value | **Pending ANSYS CSV** |
| `{CONV_TOL}` | Mesh-convergence tolerance achieved in ANSYS: % change in tip deflection between the two finest meshes | `results/ansys_vs_theory.csv` | compute from `tip_deflection_ansys_mm` of the last two rows: 100·\|a−b\|/b | **Pending ANSYS CSV** |
| `{N_MODES_ANSYS}` | Bending-Y modes compared with ANSYS | `results/ansys_modal_comparison.csv` | number of data rows | **Pending ANSYS CSV** |
| `{MAX_FREQ_ERR_ANSYS}` | Max ANSYS vs analytical frequency difference | `results/ansys_modal_comparison.csv` | max \|`ansys_vs_analytical_pct`\| | **Pending ANSYS CSV** |

## Fill-in template (three bullets)

- Developed a **from-scratch 1D Euler-Bernoulli beam finite element solver** in Python
  (NumPy/SciPy): Hermite cubic elements, consistent mass, elimination BCs, reaction
  recovery and stress recovery. It matched closed-form beam theory to within
  `{TIP_ERR_PY}` % in tip deflection, with 40+ automated pytest checks in GitHub Actions CI.
- Performed **modal analysis and h-refinement studies** over `{N_MESH_PY}` meshes. The
  first `{N_MODES_PY}` natural frequencies were within `{MAX_FREQ_ERR_PY}` % of theory, and
  I measured the convergence rate as h^`{SLOPE}`, confirming the theoretical upper-bound,
  4th-order behaviour.
- **Validated against ANSYS Mechanical** 3D solid models through a `{N_MESH_ANSYS}`-level
  mesh-convergence study (converged to `{CONV_TOL}` %). Tip deflection and mid-span stress
  agreed with theory within `{TIP_ERR_ANSYS}` % / `{MID_ERR_ANSYS}` %, and `{N_MODES_ANSYS}`
  bending modes within `{MAX_FREQ_ERR_ANSYS}` %. I identified the fixed-support stress
  singularity and explained the 1D-vs-3D differences (shear, Poisson restraint,
  discretisation).

> The third bullet depends entirely on `ansys/ansys_results.csv`. Do not use it until you
> have run the ANSYS study yourself and `compare_ansys.py` has produced the files above.
> The "40+" test count comes from the `pytest -q` summary line; re-check it if you add or
> remove tests.
