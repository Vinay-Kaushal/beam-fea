# Interview Q&A: beam-fea

Short, defensible answers tied to this project. Numbers quoted are hand calculations
(reproduced in `results/summary.md` §1). For FE errors and slopes, quote the generated files.

**1. Why are the element matrices 4×4?**
The element has 2 nodes × 2 DOFs (deflection w, rotation θ) = 4 DOFs. Ke and Me relate
the 4 nodal displacements to 4 nodal forces/moments. Rotation is needed because beam
bending energy depends on w'', which requires C¹ (slope) continuity between elements.

**2. Why is the global K singular before boundary conditions? What is the rank deficiency?**
An unsupported beam can move as a rigid body with zero strain energy, and those motions are
null vectors of K. In planar bending there are two: vertical translation and rigid rotation.
So rank(K) = ndof − 2. The test `test_global_stiffness_singular_before_bcs` checks this.

**3. Why is the tip deflection exact with just one element?**
The exact deflection under a tip load is w = Px²(3L − x)/6EI, a cubic. Hermite elements
interpolate with complete cubics, so the exact solution lies inside the FE solution space,
and the Galerkin method recovers it exactly (12.5 mm to machine precision for n = 1, 2, 5,
…). More generally, 1D FE nodal values are exact for nodal loads. For a UDL the nodal values
are still exact with consistent loads, but the field between nodes is not.

**4. Then how did you show convergence?**
With the modal problem. The exact mode shapes are trigonometric/hyperbolic, not
polynomial, so the FE frequencies are approximate. I refined from 1 to 64 elements and
fitted the slope of log(error) vs log(h) for mode 1 (value in `results/summary.md` §4;
theory predicts about 4). At the finest meshes the error hits the round-off floor, so I
excluded those points from the fit and said so.

**5. Why do FE frequencies converge from above?**
FE is a Rayleigh-Ritz method. It restricts the deformation to a finite set of shapes,
which acts like extra constraints, so the model is too stiff. With a consistent mass
matrix the Rayleigh quotient gives an upper bound for each mode. Refining enlarges the
space and lowers the frequencies toward the exact ones.

**6. What is a stress singularity? Why is the fixed-support peak stress meaningless?**
A point where linear elasticity predicts infinite stress: sharp re-entrant corners, point
loads, a fully fixed face that blocks Poisson contraction. In FE the peak there keeps
rising with every refinement, so it never converges and has no physical value. I compare
the mid-span stress (75 MPa by theory) instead. The real fix is to model the actual
geometry (fillet) or support stiffness, or to read stresses at a distance (Saint-Venant).

**7. h-refinement vs p-refinement?**
h: smaller elements of the same order. p: same mesh, higher polynomial order. p converges
faster for smooth solutions; h is more robust near singularities and is what ANSYS/HyperMesh
users do day to day (ANSYS SOLID186 is quadratic, i.e. p = 2, and is then h-refined).

**8. Why von Mises?**
Ductile metals yield when the distortion energy reaches its uniaxial-test value. Von Mises
turns a 3D stress state into one equivalent stress to compare with the yield strength. In
the 1D beam model the stress is uniaxial, so σ_vm = |σ_xx|: 150 MPa at the root, a safety
factor of 1.67 against 250 MPa yield.

**9. Consistent vs lumped mass?**
Consistent mass uses the same shape functions as K. It is a full matrix with rotational
inertia coupling, and it gives upper-bound frequencies. Lumped mass puts the mass on the
diagonal (translational DOFs). It is cheaper and ideal for explicit dynamics, but
frequencies may lie below the exact ones, and accuracy per element is usually lower for
beams.

**10. Why does ANSYS differ from beam theory?**
ANSYS solves 3D elasticity. The differences come from shear deformation (absent in
Euler-Bernoulli), Poisson restraint at the clamp, solid-element discretisation error, how
the load is applied on the end face, and the singularity at the support. For a slender
beam (L/h = 50) the converged deflection and mid-span stress should be close; the
measured differences will go in `results/ansys_vs_theory.csv` once the ANSYS study is done.

**11. How do you check a mesh is converged?**
Choose the quantity of interest and a tolerance first (e.g. a 1 % change). Refine
systematically (halve the element size) and plot the quantity vs element size. Converged =
the change between successive meshes is below the tolerance. Use quantities that *can*
converge: displacements, stresses away from singularities, frequencies. Never use a
singular peak.

**12. Verification vs validation?**
Verification: solving the equations right (code vs exact solutions, convergence rates,
unit tests). Validation: solving the right equations (model vs reality or a richer model).
My ladder: hand calc → own solver (verification) → ANSYS 3D (validation of the 1D model).

**13. How did you apply boundary conditions, and why elimination over penalty?**
I deleted the constrained rows/columns and solved the free-DOF block K_ff u_f = F_f. That
is exact for zero prescribed displacements. The penalty method adds a large number to the
diagonal, which is simple but leaves a small error and can worsen conditioning.

**14. How do you recover reactions, and why bother?**
R = K u − F at the constrained DOFs. It is a free equilibrium check: for the cantilever,
R_w = +100 N and R_θ = +50 N·m must balance the load. A mismatch reveals sign, DOF-mapping
or load errors.

**15. How do you get stress from displacements?**
Curvature κ = Σ N_j'' d_j → M = EIκ → σ = −My/I at y = ±h/2. Stress depends on second
derivatives, so it is less accurate than displacement in general. Here it is exact
because the exact moment is linear and the Hermite cubic's κ is linear.

**16. Why did your ANSYS modal analysis show modes the Python model doesn't have?**
The solid can bend in Z (about 2× the Y frequency, since I_y/I_z = (b/h)² = 4 and f ∝ √I),
twist, and stretch axially. My 1D planar model has only bending in Y. I identified modes by
their shapes and compared only `bending_y` modes.

**17. What would make this model fail?**
Short/deep beams (shear matters: use Timoshenko), large deflections (geometric
non-linearity), yielding (material non-linearity), local 3D effects at supports/holes, or
dynamic loads near resonance without damping.

**18. Why did you choose this geometry orientation in ANSYS?**
X along the length, the 10 mm depth along Y, the 20 mm width along Z, with the load in −Y.
Then the beam bends about Z with I = b h³/12 = 1666.67 mm⁴, matching the Python model.
Swapping Y/Z gives I 4× larger and a deflection 4× too small.

**19. What was the trickiest numerical issue?**
Floating-point round-off. With very fine meshes the eigenvalue error stops decreasing
(it reaches about 1e−8 relative) because the K and M entries span many orders of magnitude.
I detect that point automatically and exclude it from the slope fit instead of
reporting a corrupted rate. I also refined the tabulated βₙL constants to full precision
with a root finder, so the reference itself was not the limiting error.

**20. What would you do next?**
2D plane-stress elements (CST/Q4) to model the root singularity myself; Timoshenko beam
elements to quantify shear deformation; a lumped-mass option for comparison; frame
elements (axial + bending) for 2D trusses/frames; and repeating the study in a
**HyperMesh → OptiStruct** workflow (CBEAM/CBAR and solid CHEXA, a SOL 101 static and an
EIGRL modal card) to compare against ANSYS.

---

## Things I must be ready to explain from my resume

- [ ] Write Ke from memory and explain where it comes from (∫ Bᵀ EI B dx)
- [ ] Why 2 DOFs per node; Hermite cubics; C¹ continuity
- [ ] Why the tip deflection is exact even with 1 element, and why that isn't a bug
- [ ] Rank of K = ndof − 2 before BCs; elimination vs penalty
- [ ] Reaction recovery R = K u − F and its signs (+100 N, +50 N·m)
- [ ] Hand-calc: δ = 12.5 mm, σ_root = 150 MPa, σ_mid = 75 MPa, f₁ ≈ 32.6 Hz
- [ ] Why frequencies converge from above (Rayleigh-Ritz, consistent mass) and the h⁴ rate
- [ ] The fitted slope from `results/summary.md` and why the finest mesh was excluded
- [ ] Stress singularity at the clamp, and why mid-span stress is the converged comparison
- [ ] Every ANSYS number on my resume, traced to `ansys/ansys_results.csv` and the screenshots
- [ ] Mesh method used in ANSYS (hex sweep vs tet), element counts, student node limit
- [ ] Bending-Y vs bending-Z vs torsion mode identification
- [ ] Verification vs validation, and the verification ladder
