# ANSYS Mechanical guide: validating the beam solver with a 3D solid model

This guide walks you through reproducing the reference cantilever in **ANSYS Mechanical
(Workbench, Student version)**, recording the results in a CSV file, and comparing them
with the Python solver and beam theory.

> **Honesty rule.** The comparison script only reads numbers that *you* type into
> `ansys/ansys_results.csv`. Nothing is estimated. If a cell is blank, it is reported as
> "not recorded"; if the file is missing, everything is reported as "pending".

The beam-theory values quoted below as "expected" are the hand calculations printed in
`results/summary.md` §1 (regenerate with `python scripts/run_all.py`). They are there to
help you spot set-up mistakes. They are **not** ANSYS results.

---

## 0. Before you start

- Install **ANSYS Student** (free) from the ANSYS website. It runs on Windows only.
- Check your version's **node/element limit** (it is printed on the download page and in
  the Mechanical status bar; recent student releases allow roughly 128k nodes/elements for
  structural physics). The finest mesh in this study may go over the limit. That's
  expected; Step 4 tells you what to do.
- Create a folder `results/ansys_screenshots/` in this repo for your screenshots.
- Units: in Mechanical choose **Units → Metric (mm, kg, N, s, mV, mA)**, so that you
  can type every value below directly in mm, N and MPa.

---

## 1. Geometry: a 500 × 10 × 20 mm box

1. In Workbench, drag **Static Structural** onto the Project Schematic.
2. Open **Geometry** in SpaceClaim or DesignModeler and make a box:
   - **X (length):** 0 → 500 mm
   - **Y (depth h, bending direction):** −5 → +5 mm (10 mm)
   - **Z (width b, out of plane):** −10 → +10 mm (20 mm)

   Centre the box on the X axis (as above). The neutral axis is then y = 0 and the top
   surface is y = +5 mm, which makes the stress probe in Step 3 easy.

**Why this orientation?** The load acts in −Y, so the beam bends about the **Z axis**.
The relevant second moment of area is about Z:

    I_z = b·h³/12 = 20 × 10³ / 12 = 1666.67 mm⁴

This is the same I used by the Python model. If you swapped Y and Z, the beam would bend
about its *stiff* axis (I = h·b³/12 = 6666.67 mm⁴, 4× larger), and the deflection would
come out 4× too small. That is the most common set-up mistake in this exercise.

3. **Optional but recommended:** split the body (or just the top face) with a plane at
   x = 250 mm. This puts an edge exactly at mid-span on the top surface, so you can scope
   a stress result to it. A Construction-Geometry **Path** (Step 3) works too if you
   prefer not to split.

---

## 2. Material: Structural Steel

1. In **Engineering Data**, keep the default **Structural Steel** and check it:
   - Young's modulus **E = 200 GPa** (2e5 MPa)
   - Poisson's ratio **ν = 0.3**
   - Density **ρ = 7850 kg/m³** (7.85e-6 kg/mm³)
2. If any value differs (versions change defaults), edit it to the values above. The
   comparison is only fair if both models use the same material.
3. In Mechanical, check that the body's **Material Assignment** is Structural Steel.

---

## 3. Static Structural set-up

1. **Fixed Support:** select the face at **x = 0**.
2. **Force:** select the face at **x = 500 mm**. Define by **Components**:
   X = 0, **Y = −100 N**, Z = 0.
3. **Solution results to insert:**

   | Result | Scope | What to record |
   |---|---|---|
   | **Directional Deformation, Y axis** | whole body (or a vertex/edge at the tip) | Minimum value = tip deflection (negative = downward). Expected magnitude ≈ 12.5 mm |
   | **Equivalent (von Mises) Stress** | whole body | Maximum value. It sits at the fixed face (see Step 4) |
   | **Normal Stress, X axis** | the mid-span top edge from the split (x = 250, y = +5), or a **Path** from (250, 5, −10) to (250, 5, 10) | Average value along the edge/path. Expected ≈ +75 MPa (tension) |

   To create the path: **Model → Construction Geometry → Path**, Path Type = Two Points,
   type the coordinates above, then scope the Normal Stress result to that path. Use the
   **average** across the width; it is the quantity beam theory predicts.

4. **Solve.**

**Why the mid-span stress, not the maximum?** The fixed face is a stress singularity (Step
4). The mid-span value is far from the support and the load (Saint-Venant's principle), so
it converges with mesh refinement and can be compared fairly with σ = M·c/I.

---

## 4. Mesh convergence study

Repeat the static solve for five mesh sizes. For each one:

1. **Mesh → Insert → Sizing** (Body Sizing) on the body, with element size
   **20, 10, 5, 2.5, 1.25 mm**.
2. **Element shape:** use hexahedral elements if possible. Insert **Method → Sweep**
   (or **MultiZone**) on the body; the box sweeps cleanly. If that fails, use the default
   (tetrahedral) and **record which method you used**, because it affects the results.
   (A note in the CSV is not supported, so write it in your README.)
3. **Generate Mesh**, then read **Mesh → Statistics**: number of **Nodes** and **Elements**.
4. Solve and record:
   - `n_elements`, `n_nodes`
   - `tip_deflection_mm`: Directional Deformation Y at the tip (you may type it negative;
     the script uses the magnitude)
   - `midspan_top_stress_mpa`: Normal Stress X at mid-span top (average)
   - `max_equivalent_stress_root_mpa`: maximum Equivalent Stress anywhere (it will be at
     the fixed face)
5. **Node limit:** if a mesh exceeds the student limit, Mechanical will refuse to solve.
   Leave that row's result cells **blank** and note "exceeded student limit at X mm" in your
   README. That is a legitimate result; don't fake the row.

**What you should see:**
- Tip deflection and mid-span stress **settle** as the mesh is refined. That is
  convergence.
- The maximum equivalent stress at the fixed face **keeps rising** with refinement. That is
  a **stress singularity**: the fixed support stops the face from contracting sideways
  (the Poisson effect), so the face corners carry a theoretically infinite stress in linear
  elasticity. That peak is a modelling artefact, not a design value. Beam theory's 150 MPa
  at the root is the meaningful nominal value.

Screenshot the mesh, the deformation and the stress contour for at least one mesh, and
the final convergence table (Step 6).

---

## 5. Modal analysis

1. In Workbench, drag **Modal** onto the project and drop it on the **Model** cell of the
   Static Structural system, so it shares geometry, material and mesh. Or create it
   separately with the same set-up.
2. Mesh: Body Sizing **2.5 mm** (same method as above).
3. **Fixed Support** on the x = 0 face. **No load** (a modal analysis has no loads).
4. **Analysis Settings → Max Modes to Find = 8** (or more).
5. Solve and open each **Total Deformation** mode, using the animation to see the shape.

**Which modes are which.** A 3D solid can vibrate in ways a 2D beam model cannot.
Expect, in roughly this order:
- **Bending in Y** (the plane of the Python model). The first one is near the hand-calc
  value in `results/summary.md` (≈ 32.6 Hz).
- **Bending in Z** (sideways). This is about **2× the Y frequency** for each order, because
  I_y = h·b³/12 is 4× larger than I_z and frequency ∝ √I.
- Higher bending modes in Y and Z, then **torsion** (twisting about X) and eventually an
  axial mode.

> **Note on the spec's numbers.** The original project brief quoted ~103 Hz / ~206 Hz for
> the first Y/Z bending modes. Those values are off by a factor of √10 (see
> `docs/ASSUMPTIONS.md`). Don't use them to identify modes; use the mode *shapes*.

6. For each ANSYS mode, record a `modal` row in the CSV: `mode` = ANSYS mode number,
   `frequency_hz`, and `mode_type` = **`bending_y`**, `bending_z`, `torsion` or `other`.
   **Only `bending_y` rows are compared.** The others are ignored by the script, but
   recording them shows you understood the full spectrum.

---

## 6. Screenshots to save in `results/ansys_screenshots/`

| File name (suggested) | Content |
|---|---|
| `mesh_2p5mm.png` | Mesh view with statistics visible |
| `deformation_y.png` | Directional Deformation Y contour (true scale or noted scale) |
| `stress_equivalent.png` | Equivalent stress contour, zoomed on the fixed face |
| `stress_normal_x_midspan.png` | Normal Stress X with the mid-span path/probe visible |
| `convergence_table.png` | Your filled table (e.g. screenshot of the CSV or of Mechanical's results) |
| `mode1_bending_y.png`, `mode2_bending_y.png`, `mode3_bending_y.png` | First three bending-Y mode shapes |

---

## 7. Filling the CSV and running the comparison

1. Copy the template:

   ```bash
   cp ansys/ansys_results_template.csv ansys/ansys_results.csv
   ```

2. Open `ansys/ansys_results.csv` in a text editor or Excel (save as **CSV**, keep the
   header **exactly**). Fill in only what you measured:
   - `static` rows: `n_elements`, `n_nodes`, `tip_deflection_mm`,
     `midspan_top_stress_mpa`, `max_equivalent_stress_root_mpa`. Leave `mode`,
     `frequency_hz`, `mode_type` blank.
   - `modal` rows: the template has 8 rows numbered 1–8 (ANSYS mode numbers). Fill in
     `frequency_hz` and `mode_type`. Leave the static columns blank. Add rows if you
     extracted more modes.
3. Run:

   ```bash
   python scripts/compare_ansys.py     # comparison only
   python scripts/run_all.py           # everything, including results/summary.md
   ```

4. Outputs: `results/ansys_convergence.png`, `results/ansys_vs_theory.csv`,
   `results/ansys_modal_comparison.csv`, and the ANSYS section of `results/summary.md`.

---

## 8. Expected differences vs beam theory, and what to write afterwards

The Python model is a 1D Euler-Bernoulli beam; ANSYS solves 3D elasticity. Small
differences are expected and are *explainable*:

| Source | Effect |
|---|---|
| **Shear deformation** (ignored by Euler-Bernoulli) | Makes the solid slightly more flexible, so deflection is a little larger and frequencies a little lower. Small here because L/h = 50 (a slender beam); grows for higher modes, whose half-wavelengths are shorter |
| **Poisson restraint at the fixed face** | The clamp prevents lateral contraction, which stiffens the root region locally |
| **Solid discretisation error** | Coarse meshes (especially linear or tetrahedral elements) are too stiff; this shrinks with refinement |
| **Load application** | A force spread over the end face, rather than a point load, changes the local stress near the tip only (Saint-Venant) |
| **Stress singularity at the support** | Peak von Mises at the fixed face grows without bound with refinement. Compare mid-span stress, not the peak |

**After the study, add to the README** (replace the brackets with numbers from
`results/summary.md` / `results/ansys_vs_theory.csv`, never from memory):

> ANSYS Mechanical (3D solid, [hex/tet] mesh, [element type]) was used to validate the
> solver. At the finest converged mesh ([size] mm) the tip deflection differed from beam
> theory by [x] % and the mid-span bending stress by [y] %. The first [k] bending-Y
> frequencies agreed within [z] %. The peak stress at the fixed support rose from [a] to
> [b] MPa with refinement, confirming a stress singularity. [Mesh N mm exceeded the
> student node limit.]
