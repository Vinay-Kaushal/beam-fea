# Workflow and architecture

## The problem this project solves

Commercial FE tools (ANSYS, HyperWorks/OptiStruct) produce results whether or not the model
is right. An engineer has to know *why* a result is trustworthy. This project builds that
trust from first principles. It writes the finite element method for a beam from scratch,
proves it against exact theory and convergence behaviour, and then uses it to validate (and
be validated by) a commercial 3D solid model.

## Architecture

```
                 ┌────────────┐
                 │  config.py │  frozen dataclasses, mm→SI conversion
                 └─────┬──────┘
                       │ BeamProblem (E, I, ρA, L, P)
                 ┌─────▼──────┐
                 │ elements.py│  Hermite N, N', N''; Ke; Me; UDL fe   (4×4 / 4×1)
                 └─────┬──────┘
                       │ element matrices
                 ┌─────▼──────┐
                 │ assembly.py│  mesh, scatter → global K, M, F        (ndof × ndof)
                 └─────┬──────┘
                       │ K, M, F, BCs
                 ┌─────▼──────┐
                 │  solver.py │  elimination, K_ff u_f = F_f, R = K u − F
                 └──┬──────┬──┘
                    │      │
      ┌─────────────▼┐   ┌─▼───────────┐
      │postprocess.py│   │  modal.py   │  eigh(K_ff, M_ff) → f, φ
      │ κ → M → σ    │   └─┬───────────┘
      └──────┬───────┘     │
             │   ┌─────────▼──┐   ┌──────────────┐
             └──►│ plotting.py│   │ analytical.py│  closed-form references
                 └─────┬──────┘   └──────┬───────┘
                       │                 │
                 ┌─────▼─────────────────▼───────────────────────────┐
                 │ scripts/: run_static, run_modal, run_convergence, │
                 │           compare_ansys (reads ansys/*.csv),      │
                 │           run_all → summary.md                    │
                 └─────────────────────┬─────────────────────────────┘
                                       ▼
                                  results/  (CSV, PNG, summary.md)
```

## Data flow: one static run (`run_static.py`, default case)

1. `default_problem()` builds the problem in SI (L = 0.5 m, I = 1.6667e−9 m⁴, P = 100 N).
2. `generate_mesh(L, 10)` → 11 nodes, 10 elements.
3. `assemble_stiffness` → 22×22 K (singular).
4. `assemble_load` puts −100 N on DOF 20 (w at node 10).
5. `solve_static(K, F, cantilever_bcs())` eliminates DOFs 0 and 1, solves the 20×20 system,
   rebuilds u (22), and computes reactions R = K u − F at DOFs 0, 1.
6. `evaluate_fields` interpolates w, θ, κ inside each element (11 points per element) →
   M = EIκ → σ_top, σ_bottom, σ_vm.
7. Compare with `analytical.py`, write `static_results.csv`, plot two PNGs.

## Data flow: one modal run (`run_modal.py`)

1. `assemble_beam(L, 20, EI, ρA)` → mesh, 42×42 K and M.
2. `solve_modal` extracts K_ff, M_ff (40×40), calls `scipy.linalg.eigh`, keeps 5 modes,
   expands them to full DOF vectors and fixes their signs.
3. f = √λ / 2π compared with `cantilever_natural_frequencies` → `modal_results.csv`.
4. Mode shapes interpolated with Hermite functions → `mode_shapes.png`.

## The verification ladder

```
  Rung 3  ANSYS Mechanical 3D solid      (validation: different, richer physics)
            ▲   compare_ansys.py
  Rung 2  Own Python FE solver           (verification: exact match, convergence rate)
            ▲   run_static / run_modal / run_convergence + pytest
  Rung 1  Hand calculation / closed form (analytical.py, summary.md §1)
```

Each rung is trusted only after it agrees with the one below. Where it disagrees, the
difference must be explained (see `docs/CONCEPTS.md` §17).

## Module responsibilities

| Module | Responsibility |
|---|---|
| `config.py` | Frozen problem dataclasses, study settings, mm/MPa ↔ SI helpers |
| `elements.py` | Hermite shape functions and derivatives, Ke, consistent Me, UDL load vector |
| `assembly.py` | Uniform mesh, DOF numbering, global K, M, F assembly |
| `solver.py` | BC definitions, elimination, static solve, reaction recovery |
| `postprocess.py` | In-element fields, moment, bending/von Mises stress, tables and CSV output |
| `modal.py` | Constrained generalized eigenproblem, frequencies, mass-normalised modes |
| `analytical.py` | Closed-form deflections, moments, stresses, cantilever frequencies |
| `plotting.py` | All figures (headless Agg backend, 150 dpi) |
| `scripts/run_static.py` | Static verification table and plots |
| `scripts/run_modal.py` | Modal verification table and mode-shape plot |
| `scripts/run_convergence.py` | Modal h-refinement and fitted convergence slope |
| `scripts/compare_ansys.py` | Reads the hand-filled ANSYS CSV and compares (or reports "pending") |
| `scripts/run_all.py` | Runs everything, writes `results/summary.md` |

## Reproduce everything

```bash
git clone https://github.com/Vinay-Kaushal/beam-fea.git && cd beam-fea
python -m pip install -r requirements.txt     # or: make install
python -m pytest -q                           # or: make test
python scripts/run_all.py                     # or: make run
# after ANSYS (see ansys/ANSYS_GUIDE.md):
cp ansys/ansys_results_template.csv ansys/ansys_results.csv   # then fill it in
python scripts/compare_ansys.py
python scripts/run_all.py                     # refreshes summary.md with ANSYS data
```

`make clean` removes generated CSV/PNG files and caches.

## Technologies and why

| Technology | Why |
|---|---|
| **Python 3.10+** | Readable, standard in engineering automation (ANSYS and HyperWorks both script in Python) |
| **NumPy** | Dense linear algebra and vectorised shape-function evaluation |
| **SciPy** | `linalg.eigh` for the symmetric generalized eigenproblem; `brentq` to refine βₙL roots |
| **Matplotlib** | Publication-quality static plots, headless in CI |
| **pytest** | Executable verification: every claim (exactness, rank, upper bounds) is a test |
| **GitHub Actions** | Proves on every push that tests pass and the pipeline runs headless on Python 3.10 and 3.12 |
| **ANSYS Mechanical** | Industry-standard commercial FE tool; an independent 3D solid model to validate against |
