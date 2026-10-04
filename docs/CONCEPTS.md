# FEA concepts, explained through this project

Every section ends with a **worked example** from the reference problem: a steel cantilever
with L = 500 mm, b = 20 mm (width), h = 10 mm (depth), E = 200 GPa, ρ = 7850 kg/m³, and a
100 N downward tip load. The worked examples use closed-form hand calculations, the same
ones the code prints in `results/summary.md` §1. For FE outputs (errors, slopes) this file
points to the generated files instead of quoting numbers.

---

## 1. What FEA is and why it exists

Most real structures are governed by differential equations (equilibrium, compatibility,
material law) that have no closed-form solution once the geometry, loads or supports get
complicated. The **finite element method** replaces the unknown continuous field (for
example the deflection w(x)) with a finite number of unknowns at points called **nodes**.
Between nodes the field is interpolated by simple functions. The differential equation then
becomes a set of linear algebraic equations, **K u = F**, which a computer can solve.

*Worked example:* our cantilever has a known exact answer (δ = PL³/3EI = 12.5 mm), which
is exactly why we use it: we can check that the machinery works before trusting it on
problems with no exact answer.

## 2. Discretisation: nodes, elements, DOFs

- **Node:** a point where unknowns live.
- **Element:** the piece of structure between nodes (here, a beam segment of length Le).
- **Degree of freedom (DOF):** an independent unknown at a node. A planar beam node has two:
  deflection **w** (m) and rotation **θ = dw/dx** (rad).

With n elements there are n + 1 nodes and 2(n + 1) DOFs. Node i owns global DOFs
2i (w) and 2i + 1 (θ).

*Worked example:* 10 elements → 11 nodes → 22 DOFs; the clamp removes 2, so 20 unknowns.

## 3. Euler-Bernoulli assumptions, and when they fail

1. Plane sections remain plane **and perpendicular to the neutral axis**, so there is no
   shear deformation.
2. Small deflections and rotations (linear geometry).
3. Linear-elastic, isotropic material.
4. Stress is uniaxial (σ_xx only); the Poisson effect across the section is ignored.

They fail for **short/deep beams** (L/h below about 10, where shear deformation matters:
use Timoshenko theory), large deflections (geometric non-linearity), yielding, and near
supports or loads where the 3D stress state matters (Saint-Venant zones).

*Worked example:* L/h = 500/10 = 50, so the beam is slender and Euler-Bernoulli is an
excellent model. Peak transverse shear stress is τ = 1.5 V/A = 1.5 × 100 / 200 mm² =
0.75 MPa, against 150 MPa bending stress. Shear is negligible.

## 4. Shape functions, and why Hermite cubics for beams

Inside an element, w(x) = Σ N_j(x)·d_j, where d = [w1, θ1, w2, θ2]. The beam's strain energy
contains the **second derivative** w''. For the energy to be finite across element
boundaries, both w **and** w' must be continuous (**C¹ continuity**). Four nodal values
fix a cubic polynomial uniquely, giving the **Hermite cubics** (ξ = x/Le):

    N1 = 1 − 3ξ² + 2ξ³          N2 = Le(ξ − 2ξ² + ξ³)
    N3 = 3ξ² − 2ξ³              N4 = Le(−ξ² + ξ³)

Check: N1 + N3 = 1, so a rigid translation is reproduced exactly (tested in
`tests/test_elements.py`).

*Worked example:* the exact cantilever deflection under a tip load,
w(x) = P x²(3L − x)/(6EI), **is a cubic**. A cubic element can represent it exactly, so
even **one element** gives the exact tip deflection (12.5 mm) and rotation (0.0375 rad).
This is expected behaviour, not a bug. See `results/static_results.csv`, rows n = 1…50.

## 5. Deriving the element matrices

**Stiffness.** Strain energy U = ½ ∫ EI (w'')² dx. Substituting w'' = B·d with
B = [N1'', N2'', N3'', N4''] gives U = ½ dᵀ (∫ Bᵀ EI B dx) d, so

    Ke = ∫₀^Le Bᵀ EI B dx = (EI/Le³) [[12, 6Le, −12, 6Le], [6Le, 4Le², −6Le, 2Le²], …]

**Consistent mass.** Kinetic energy T = ½ ∫ ρA ẇ² dx with ẇ = N·ḋ gives

    Me = ∫₀^Le ρA Nᵀ N dx = (ρA Le/420) [[156, 22Le, 54, −13Le], …]

It is called "consistent" because it uses the same shape functions as K.

**Consistent load (UDL).** fe = ∫ Nᵀ q dx = q Le [½, Le/12, ½, −Le/12]. Note the end
**moments**. Just lumping half the load on each node would lose accuracy.

*Worked example:* with Le = 50 mm, the rigid translation vector v = [1, 0, 1, 0] gives
vᵀ Me v = ρA·Le = 1.57 kg/m × 0.05 m = 0.0785 kg, the element's real mass (tested).

## 6. Assembly

Each element's 4×4 matrix is added into the global matrix at the rows/columns of its
DOFs. Where two elements share a node, their contributions **add**. That addition is what
enforces continuity of w and θ and equilibrium of forces/moments at the node.

*Worked example:* 2 elements → global K is 6×6; the middle node's 2×2 block is the sum of
element 1's lower-right block and element 2's upper-left block.

## 7. Boundary conditions, and why K is singular without them

A free beam can translate and rotate as a rigid body with **zero strain energy**. Those
motions are null vectors of K, so K is singular (not invertible). In planar bending there
are **2 rigid-body modes** (vertical translation and rotation), so rank(K) = ndof − 2
(tested). Supports remove them. This project applies BCs by **elimination**: it deletes the
constrained rows/columns and solves the reduced system. That is exact, unlike the penalty
method, which adds a huge stiffness and introduces a small error.

*Worked example:* cantilever: w0 = θ0 = 0 (2 DOFs removed). Simply supported: w at both
ends = 0, rotations free (also 2 DOFs, enough to stop both rigid modes).

## 8. Solving K u = F

After elimination, K_ff is symmetric positive definite, so a direct solve
(`numpy.linalg.solve`) is robust. For large models, commercial codes use sparse direct
solvers or iterative solvers (PCG), but the idea is the same.

*Worked example:* the tip w-DOF gets F = −100 N; the solution gives w_tip = −12.5 mm.

## 9. Reactions as an equilibrium check

At every DOF, internal force K u = external load F + reaction R. So **R = K u − F**,
evaluated at the constrained DOFs. If the reactions don't balance the applied loads,
something is wrong (a load in the wrong DOF, a sign error, a missing element).

*Worked example:* R_w = +100 N (upward) and R_θ = +P·L = +50 N·m (counter-clockwise), which
balance the −100 N tip load and its moment about the root (tested).

## 10. Stress recovery: curvature → moment → stress

1. Curvature κ(x) = w''(x) = Σ N_j''(x) d_j (from the shape-function second derivatives)
2. Bending moment M = EI·κ (sagging positive)
3. Bending stress σ = −M·y/I, with y = ±h/2 at the extreme fibres

*Worked example:* at the root, M = −P·L = −50 N·m (hogging). Top fibre (y = +5 mm):
σ = −(−50)(0.005)/(1.6667e−9) = **+150 MPa (tension)**. At mid-span, M = −25 N·m, so
σ_top = **75 MPa**. Note that stresses come from *derivatives* of the displacement, so in
general they converge more slowly than displacements.

## 11. Von Mises stress

Von Mises combines a 3D stress state into one scalar to compare with uniaxial yield:

    σ_vm = √(½[(σ1−σ2)² + (σ2−σ3)² + (σ3−σ1)²])

For a uniaxial state (σ_xx only) this reduces to **|σ_xx|**. In this beam model, shear is
neglected, so σ_vm = |σ|.

*Worked example:* max σ_vm = 150 MPa at the root, so the safety factor against yield is
250/150 = 1.67.

## 12. Mesh convergence and h-refinement

**h-refinement:** make elements smaller (more of them). **p-refinement:** raise the
polynomial order of each element. A result is "converged" when refining further changes it
less than a tolerance you set in advance. On a log-log plot of error vs element size h,
the slope is the **convergence rate**.

*Worked example:* static results here are exact for any mesh (§4), so the meaningful
convergence study is the **modal** one. `scripts/run_convergence.py` measures the rate for
mode 1. Theory predicts eigenvalue error ∝ h^(2(p+1−m)) = h^4 for cubics (p = 3) and a
2nd-order problem (m = 2). The fitted value is in `results/summary.md` §4 and
`results/modal_convergence.png`. At the finest meshes the error reaches the floating-point
**round-off floor** (about 1e−8 relative), where it stops falling. Those points are
excluded from the fit and labelled on the plot.

## 13. Stress singularities

At a re-entrant corner, a point load or a fully fixed face, linear elasticity predicts
**infinite** stress. In FE, the peak stress there grows every time you refine the mesh. It
never converges. Such peaks are not "results". Engineers instead (a) read stresses away
from the singularity (Saint-Venant), (b) model the real geometry (fillets) or the real
support stiffness, or (c) use fracture/fatigue methods built for singular fields.

*Worked example:* in the ANSYS solid model, the clamp at x = 0 stops lateral (Poisson)
contraction, so the corners of the fixed face are singular. Their von Mises peak rises
with refinement. The mid-span stress (75 MPa by theory) is the converged quantity to
compare. The 1D beam has no singularity: it reports the nominal 150 MPa.

## 14. Modal analysis and the generalized eigenproblem

Free undamped vibration: M ü + K u = 0. Try u = φ sin(ωt):

    (K − ω² M) φ = 0   →   K φ = λ M φ,  λ = ω²

After eliminating the constrained DOFs, `scipy.linalg.eigh(K_ff, M_ff)` returns λ ascending
and **mass-normalised** modes (φᵀ M φ = 1). Then f = √λ / 2π.

*Worked example (hand calc):* ω₁ = (β₁L)² √(EI/(ρA L⁴)) with β₁L = 1.8751. EI = 333.33 N·m²,
ρA = 1.57 kg/m, L⁴ = 0.0625 m⁴ → √(…) = 58.28 s⁻¹ → ω₁ = 204.9 rad/s → **f₁ ≈ 32.6 Hz**.
Higher modes follow the ratios (β_nL/β₁L)²: 6.27× and 17.55× (≈ 204 Hz and ≈ 572 Hz).

## 15. Why FE overestimates frequencies

The FE solution is a **Rayleigh-Ritz** approximation. It only allows deflected shapes that
can be built from the shape functions. Restricting the possible shapes acts like extra
constraints, making the model **stiffer**. With a *consistent* mass matrix, the Rayleigh
quotient ωₙ² = φᵀKφ / φᵀMφ is therefore an **upper bound**, and refining the mesh brings the
frequencies down toward the exact values. (A *lumped* mass matrix breaks this guarantee;
its frequencies can be below the exact values.)

*Worked example:* every FE frequency in `results/modal_results.csv` and
`results/modal_convergence.csv` has a positive error (FE > exact).

## 16. Verification vs validation

- **Verification:** "are we solving the equations right?" Compare the code with exact
  solutions, check convergence rates, run unit tests. (Hand calc → own solver.)
- **Validation:** "are we solving the right equations?" Compare the model with physical
  reality, or with a higher-fidelity model. (Own 1D solver → ANSYS 3D solid → ideally a
  physical test.)

This project's **verification ladder**: hand calculation → own Python FE → ANSYS
Mechanical.

## 17. Sources of error between the 1D beam and the 3D solid

| Source | Direction / size |
|---|---|
| Shear deformation (absent in Euler-Bernoulli) | Solid slightly more flexible; frequencies slightly lower; grows with mode number and with h/L |
| Poisson restraint at the clamp | Local stiffening and a stress singularity at the root |
| 3D end effects (load spread on the end face) | Local only (Saint-Venant) |
| Solid discretisation (element type, tet vs hex, size) | Coarse/linear/tet meshes too stiff; shrinks with refinement |
| Extra modes in 3D (bending Z, torsion, axial) | Not errors, but they must be filtered out before comparing |

## 18. Limits of this project

- Linear, static or free-vibration only. No damping, forced response, buckling, contact or
  plasticity.
- 1D, single bending plane. No axial or torsional DOFs, no out-of-plane bending.
- Uniform prismatic section, uniform mesh, a single material.
- Dense matrices, fine for < 1000 DOFs but not for industrial models.
- ANSYS validation is against another model, not a physical test.
