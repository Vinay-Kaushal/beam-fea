# Assumptions and deviations

## Deviation: the spec's frequencies are off by √10
The spec quoted f₁ ≈ 103 Hz (and 646 / 1810 Hz). For this beam:
EI = 333.33 N·m², ρA = 1.57 kg/m, so √(EI/(ρA·L⁴)) = 58.28 s⁻¹ and
f₁ = 1.8751² × 58.28 / 2π = **32.6 Hz**. The ratios between modes are right, but every
value is √10 too high. All static sanity values match exactly. The tests check f₁ against
this hand calculation, not against 103 Hz. In the ANSYS guide, bending Z ≈ 2× bending Y.

## Decisions
- `tests/conftest.py` (not in the spec's tree) adds `src/` to the import path; the scripts
  do the same. No install step is needed.
- Sign convention: w and loads positive in +y, θ = dw/dx, M = EIw'' (sagging positive),
  σ = −My/I. Tables report downward deflection as a positive magnitude.
- Table/CSV helpers live in `postprocess.py`; the `assemble_beam()` helper lives in `assembly.py`.
- UDL case: q = 200 N/m. Simply supported case: 100 N at mid-span. Both use 10 elements.
- βₙL roots are refined with `brentq`, so the reference values aren't the limiting error.
- `eigh(K_ff, M_ff)` is kept as specified. Its round-off floor is about 1e−8 relative, so the
  slope fit runs from n = 2 to the last mesh before the error stops decreasing. Excluded
  meshes are reported on the console, the plot and in `summary.md`.
- ANSYS: the script uses magnitudes of deflection and stress, compares only `bending_y`
  modes, writes no files while data are pending, and treats a malformed CSV as pending
  (exit 0).
- `results/` is tracked in git so the README images render on GitHub.
- Docs quote only inputs and closed-form hand calculations. FE and ANSYS numbers are
  referenced by file.

## Modelling assumptions
Euler-Bernoulli (no shear deformation, no Poisson effect), linear elastic, small
deflection, undamped. Von Mises = |σ| (uniaxial); peak shear stress (0.75 MPa) is neglected.
