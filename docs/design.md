# flatmat — design

Status: **design phase, phase-0 spike done.** The spike
([`spike/`](../spike/README.md)) runs §1 end to end on OCP 7.8 and 8.0. The
results are folded into this document (marked *spike*). Anything still
marked *verify* has not been checked yet.

Research behind every decision below, with sources, is in
[`research.md`](research.md).

## 0. Goal and decisions

flatmat reads a sheet metal part from STEP, recognises it as sheet metal
(thickness, flanges, bends) and produces a flat pattern with bend lines, a
bend table and a DXF that a laser cutter and a press brake can use.

| Decision | Choice | Why (short) |
|---|---|---|
| Licence | Apache-2.0 | Same as build123d and the OCP bindings; patent grant; usable from proprietary software. FreeCAD (LGPL) code is read for ideas, never copied. |
| Language | Python | Run time is dominated by OCCT's STEP reader and healing, which is C++ in every option. Unfolding is O(faces): milliseconds. |
| Kernel | `cadquery-ocp-novtk` | Headless OCCT bindings, wheels for CPython 3.11–3.14 on Linux/macOS/Windows. Target 8.0.x, CI also on 7.8.x. |
| Runtime dependencies | OCP, ezdxf (+ pyclipper or shapely for 2D checks) | No build123d (unreleased sheet metal, pins OCP 8.0, pulls scipy/sympy/ipython), no networkx (a plain BFS suffices). |
| Concurrency | one process per part | OCCT's STEP reader was not reliably thread-safe up to 7.9; OCP has no free-threaded wheels. |
| C++ | only if the phase-0 benchmark shows a need | A C++ extension can only exchange `TopoDS_Shape` with OCP when built against the identical OCCT build. |

Public API rule: **no OCP types in the public API** (escape hatch:
`SheetPart.occ_shape`). Everything a user reads is frozen dataclasses of
floats, tuples and enums.

---

## 1. Algorithm

### 1.0 Input and preparation (`flatmat/io.py`)

1. Read with `STEPControl_Reader`, unit forced to mm. Healing is
   configurable (`heal=True|False|"light"`): on perforated parts healing is
   about 90 % of import time (Quaoar).
2. Collect solids (`TopExp_Explorer(SOLID)`). Exactly one expected, else
   `MULTIPLE_SOLIDS` — or `solid=<index>`; assemblies see §6.4.
3. `ShapeUpgrade_UnifySameDomain` (UnifyFaces, UnifyEdges, no BSpline
   concatenation) merges coplanar split faces and seam-split cylinders. It
   only merges faces on the *same* underlying surface, so plane/cylinder
   boundaries survive (*verify*). Keep its history so diagnostics can name
   original faces.
4. v0.2: canonical recognition of B-spline faces as plane/cylinder/cone
   (§5.5).

### 1.1 Face table and attributed adjacency graph (`flatmat/recognize.py`)

Per face: id, surface type (`BRepAdaptor_Surface.GetType`),
orientation-corrected normal (flip if `REVERSED`); for cylinders the axis
(`gp_Ax1`), radius R and convexity (normal pointing away from or towards the
axis); `BRepTools.UVBounds`; area (`BRepGProp.SurfaceProperties`).

Per edge (`TopExp.MapShapesAndAncestors(EDGE, FACE)`):

| ancestors | meaning |
|---|---|
| 1, and `BRep_Tool.IsClosed(edge, face)` | seam → ignore |
| 1, otherwise | `OPEN_SHELL` |
| > 2 | `NON_MANIFOLD` |
| 2 | dihedral at the edge midpoint from both pcurves (`BRepLProp_SLProps` normals). Angle < `tol.angular` → **SMOOTH**; else sign of `(n1 × n2) · t` (t oriented along face 1) → **CONVEX**/**CONCAVE**. |

### 1.2 Thickness

- **Cylinder pairs:** axes parallel and within `tol.linear` of each other,
  overlapping axial spans (rules out counterbores), ΔR > tol → candidate ΔR,
  weight = area. Candidates are **bucketed by axis line** (canonical
  direction + foot point, rounded). *Spike*: comparing all pairs was O(n²)
  on the hole walls of perforated parts.
- **Plane pairs, analytically:** for the largest skin planes, take a point
  inside the face, its distance to each antiparallel plane of the other
  skin, and test whether its foot lies inside that face. *Spike*: a ray cast
  against the whole solid (`IntCurvesFace_ShapeIntersector`) was O(n²) on
  perforated parts. No O(n²) booleans either (build123d's shadow solids).
- **Point-in-face tests** use one `IntTools_FClass2d` built per face and
  reused. *Spike*: `BRepClass_FaceClassifier` rebuilds per call (~9 ms on a
  face with 1000 holes against 0.8 ms); `BRepTopAdaptor_FClass2d` no longer
  exists in OCP 8.0.
- `t` = weighted mode of the clustered candidates; `confidence` = share of
  skin area that agrees; < 0.95 → `VARIABLE_THICKNESS`. A user-given `t` is
  validated, not trusted.

### 1.3 Skins, thickness faces, bends

- Connected components over SMOOTH edges. The two largest by area are the
  **skins**; every other component is a **thickness face** (side wall or hole
  wall).
- Checks:
  - every thickness component touches both skins via sharp edges — otherwise
    `EDGE_FEATURE` (chamfer, countersink; see §5.3);
  - area identity `A_outer − A_inner = Σ_bends θ · t · L` → else
    `NOT_TWO_SKINS`.
- **Reference skin A**: default the skin containing the largest plane
  (option `reference="auto"|"inner"|"outer"`). In the flat, A faces the
  viewer (+Z); `FlatPattern.mirrored()` gives the other side and swaps
  up/down.
- **Bends** = cylinder faces of A:
  - `r_in` = R if A is concave on that face, else R − t — decided **per
    cylinder face**, never per skin (joggles, §5.1);
  - `angle` (deflection) = u-span; `opening_angle = 180° − angle`;
  - `direction` = UP if A is the inside of the bend (flange turns towards
    the viewer);
  - `is_hem` if angle ≥ 179°; angle ≥ 359.9° → `BEND_TOO_LARGE`.
- **Grouping**: non-adjacent cylinder faces with the same axis line, r and
  angle are one press-brake bend with several `segments` (a bend
  interrupted by a cut-out). Bend lines per segment, one table row per bend.
- **Flanges**: components of A once bend faces are removed.
- Cones, tori, B-splines (before canonical recognition) →
  `UNSUPPORTED_SURFACE` with face ids.

### 1.4 Developing each face on its own (`flatmat/unfold.py`)

Every face of skin A is developed independently into a local 2D frame,
following build123d's design and Quaoar part 5 (not FreeCAD's chained bend
matrices):

- **Plane**: pcurves on a `gp_Pln` are already metric 2D curves in the
  plane's frame. A `REVERSED` face is mirrored (u → −u) so every developed
  face has normal +Z. Circles stay exact `Geom2d_Circle`.
- **Cylinder**: map `(u, v) → (R_n · u, v)` with `R_n = r_in + K·t` (sign
  handling §1.7). Lines stay lines. B-splines: scale the poles, keep degree,
  knots, multiplicities and weights — exact, because an affine map preserves
  NURBS (FreeCAD's unfolder drops degree/knots/mults here). Other 2D curve
  types are converted first with `Geom2dConvert.CurveToBSplineCurve`.
  Periodic pcurves are shifted by k·2π into `[umin, umax]`.
- Output per face: 2D wires (outer + inner); every 2D segment keeps its
  **`source_edge`**.

### 1.5 Placement by vertex identity

- Root: the largest plane of A (option `root=`); identity transform.
- BFS over SMOOTH edges of A. For tree edge parent P → child C over shared
  edge e:
  - take the 2D points of the **same `TopoDS_Vertex`** V1, V2 in both
    developments (`BRep_Tool.Parameter(vertex, edge, face)` evaluated on the
    pcurve) — the correspondence is unambiguous, no UV-corner case analysis;
  - rigid 2D transform (rotation + translation, **never** reflection, since
    all faces are +Z) from the two point pairs;
  - check both images of e have equal length (`DEVELOP_FAILED`); tangent
    plane/cylinder edges are rulings, straight and unaffected by the u
    scaling;
  - check C lies on the opposite side of e from P (`DEVELOP_FOLDBACK`).
- **Non-tree edges (cycles)**: both images must coincide within
  `tol.flat`, else `CYCLE_INCONSISTENT` naming the edge as a suggested cut
  (tubes, §5.6).

### 1.6 Flat pattern (`flatmat/flat2d.py`)

- **Outline without booleans**: boundary edges of A (edges with only one
  face in A, i.e. against a thickness face) are mapped to 2D and chained
  into loops **by 3D vertex identity** — exact, no fuzzy endpoint matching.
  Outer loop = largest |signed area|; the rest are holes (checked by
  containment).
- **Bend lines**: iso-u at `u_mid` of each segment (x = R_n·u_mid in the
  developed frame), clipped to the developed face → inner intervals. The
  tangent edges (strip borders) go to an optional `BEND_EXTENT` layer.
- **Checks**:
  - self-overlap: tessellate loops at 0.01 mm, test with pyclipper/shapely →
    `FLAT_OVERLAP` naming the flanges;
  - area identity `A_flat = Σ A_plane + Σ A_cyl · R_n / R_A` within 1e-6
    relative → `AREA_MISMATCH`. A strong self-check that costs nothing.
- Output geometry stays exact (`Line`, `Arc`, `Spline`). Fitting to
  lines/arcs/biarcs is an export option with an explicit tolerance.

### 1.7 Bend model (`flatmat/bendmodel.py`)

- `BA = θ · (r_in + K·t)`
- `OSSB = tan(θ/2) · (r_in + t)` — undefined for hems (θ ≥ 179°), returned
  as `None`
- `BD = 2·OSSB − BA`
- K models:
  - scalar;
  - `KTable` r/t → K, piecewise linear with bisection, clamped, warning
    `K_EXTRAPOLATED`;
  - `KTable.din6935()`: `K = k_DIN / 2` with `k_DIN = 0.65 + 0.5·lg(r/s)`.
    **For r/s > 5 DIN sets k = 1** (K = 0.5); the switch is continuous
    (0.9995 at 5).
    - DIN 6935 does **not cover r/s ≤ 0.65**: its Table 5 starts at
      "over 0.65". Below that flatmat clamps at k ≥ 0 and warns
      `K_OUTSIDE_STANDARD`. The formula reaches 0 at r/s ≈ 0.05; the clamp
      is a practitioner convention.
    - The DIN compensation value `v` (legs dimensioned to the outer vertex,
      β = opening angle = 180° − bend angle):
      - β 0–90°: `v = π(180°−β)/180°·(r + s/2·k) − 2(r + s)`
      - β 90–165°: `v = π(180°−β)/180°·(r + s/2·k) − 2(r + s)·tan((180°−β)/2)`
      - β 165–180°: v = 0
    - These were rebuilt from a near-verbatim English rendering of the
      2011-10 text and match all 7 of its worked examples to 0.01 mm. Some
      German sources flip the sign (`L = a + b − v`).
  - callable `k_for(bend, t) -> float` for material-, angle- or
    hem-dependent values.
- Neutral fibre sits K·t from the **inside** skin. Measured from reference
  skin A: K·t if A is the inside, else (1 − K)·t (build123d's
  `surface_offset` logic).

---

## 2. API and data model

```python
import flatmat as fm

pat = fm.unfold_step("bracket.step", k=0.44)                      # one call

model = fm.load_step("bracket.step", heal=True)                   # -> StepModel
sheet = fm.recognize(model, thickness=None, tol=fm.Tolerances())  # -> SheetPart
pat   = fm.unfold(sheet, k=fm.KTable.din6935(), root="largest", reference="auto")
report = fm.check(pat, rules=fm.DfmRules.default(), tooling=None) # §7, optional

fm.write_dxf(pat, "bracket.dxf", layers=fm.LayerMap.default(), labels="delem", version="R2010")
fm.write_bend_table(pat, "bracket.json")                          # or .csv
fm.write_svg(pat, "bracket.svg")
```

CLI:

```
flatmat unfold part.step -k 0.44 -o part.dxf --table part.json [--svg part.svg]
               [--thickness 2] [--reference inner] [--check]
```

Non-zero exit code on refusal; the diagnostic goes to stderr as JSON.

### Types (frozen dataclasses)

- `Tolerances(linear=1e-3, angular=1e-4, thickness_rel=1e-3, flat=1e-2)` —
  mm and rad; see §6.3 for how they scale.
- `SheetPart`: `thickness`, `thickness_confidence`, `faces: dict[FaceId,
  FaceInfo]`, `skins: tuple[SkinId, SkinId]`, `reference: SkinId`,
  `flanges`, `bends`, `thickness_faces`, `features` (§5.3), `diagnostics`,
  `source: Provenance` (file, solid index, units, healing applied).
- `Bend`: `id`, `segments: list[BendSegment(face_id, length, axis_p,
  axis_d)]`, `r_in`, `angle_deg`, `opening_angle_deg`, `direction:
  Literal["up", "down"]`, `kind: Literal["bend", "hem", "sharp", "joggle"]`,
  `flanges: tuple[FlangeId, FlangeId]`; after unfolding `k`, `ba`, `bd`,
  `ossb`.
- `Flange`: `id`, `face_ids`, `area`, `normal`.
- `FlatPattern`: `outer: Loop2D`, `holes: list[Loop2D]`, `bend_lines:
  list[BendLine]`, `bend_zones`, `thickness`, `k_model`, `bbox`, `area`,
  `face_transforms`, `bends`, `features`, `diagnostics`; `mirrored()`,
  `to_polylines(tol)`.
- `Loop2D` = list of `Seg2D` = `Line(p0, p1)` | `Arc(c, r, a0, a1, ccw)` |
  `Spline(deg, poles, knots, mults, weights)`, each with `source_edge`.
- `BendLine`: `bend_id`, `segment_idx`, `p0`, `p1`, `direction`,
  `angle_deg`, `r_in`.
- `Diagnostic(code, severity, message, face_ids, edge_ids, data)`. Severity
  `error` raises `FlatmatError(diagnostic)`; `warning` is attached to the
  result.

### Diagnostic codes

| group | codes |
|---|---|
| input | `NO_SOLID`, `MULTIPLE_SOLIDS`, `OPEN_SHELL`, `NON_MANIFOLD`, `HEALING_CHANGED_TOPOLOGY`, `DIRTY_TOLERANCE`, `MISALIGNED_BODIES` |
| thickness/skins | `NO_THICKNESS`, `VARIABLE_THICKNESS`, `NOT_TWO_SKINS` |
| geometry | `EDGE_FEATURE`, `UNSUPPORTED_SURFACE`, `BEND_TOO_LARGE`, `BEND_UNPAIRED`, `FORMED_FEATURE`, `HEM_MERGED` |
| unfolding | `CYCLE_INCONSISTENT`, `DEVELOP_FAILED`, `DEVELOP_FOLDBACK` |
| flat checks | `FLAT_OVERLAP`, `AREA_MISMATCH` |
| warnings | `K_EXTRAPOLATED`, `COUNTERSINK_PROJECTED`, `SHARP_BEND_ASSUMED_RADIUS`, `SURFACE_CONVERTED`, `NEAR_TANGENT_KINK`, `HEM_CONTACT`, `HEM_DEFAULT_ALLOWANCE`, `INSIDE_OUT_FIXED`, `SLIVER_FACE` |
| DFM (§7) | `DFM_*` — see §7 |

### Bend table (JSON, versioned schema)

```json
{"schema": "flatmat.bendtable/1", "units": "mm", "thickness": 2.0, "k_model": "din6935",
 "bends": [{"id": 1, "direction": "up", "angle_deg": 90.0, "opening_angle_deg": 90.0,
            "r_in": 2.0, "k": 0.42, "ba": 4.46, "bd": 3.54, "ossb": 4.0, "kind": "bend",
            "segments": [{"line": [[0, 25.1], [120, 25.1]], "length": 120.0}]}]}
```

CSV: the same content, flattened to one row per segment.

### DXF

- R2010 by default (Amada Dr.Abe reportedly needs ≤ 2010), `$INSUNITS = 4`
  (mm), everything in the XY plane at 1:1.
- Default layers `OUTER`, `INNER`, `BEND_UP`, `BEND_DOWN`, `BEND_EXTENT`
  (off), `ETCH`, `TEXT` — all configurable through `LayerMap`.
- Bend lines in linetype `CENTER`, never hidden (SendCutSend does not detect
  hidden bend lines); continuous across cut-outs where the shop wants it
  (option), otherwise per segment.
- Arcs native; splines as `SPLINE` or fitted.
- `labels="delem"`: `+90` / `-90` (up/down), `H` for hems, `R2` for radius
  next to each bend line, as Delem controls read them. Also `"plain"` or
  `None`.

---

## 3. Testing and benchmark

### Fixture generator (`tests/fixtures/gen.py`, OCP only, no licence issues)

A 2D profile (flat lengths, angles, r, t) with inner fillet r and outer
fillet r + t, extruded by a width. Shapes: L, U, Z, hat, hem, joggle, box
with corner reliefs. Variants: holes in a flange, **a hole through a bend**,
a cut-out splitting a bend (grouping), seam-split cylinders
(`ShapeUpgrade_FaceDivide`), NURBS-converted (`BRepBuilderAPI_NurbsConvert`,
v0.2), each round-tripped through STEP.

### Tests

- **Analytic**: flat length = Σ inner flat lengths + Σ BA within 1e-6;
  recovered r, angle, direction and t match the generator input; bend model
  formulas against textbook values.
- **Property-based (hypothesis)**: random profiles with 1–5 bends, angles
  10–170°, r 0.5–10, t 0.5–6, K 0.3–0.5 → flat length, area identity, bend
  count, angle, r, direction, no overlap.
- **Refusals**: two solids, tapered plate, cone bend, chamfered edge, open
  shell, tube with an inconsistent cycle.
- **Broken STEP corpus** (§6): synthetic defects (enlarged tolerances,
  stripped pcurves, split faces, reversed faces) each expecting either a
  correct result after healing or a specific code.

### BenDFM (`bench/bendfm.py`)

20 000 synthetic bent parts with folded STEP, unfolded STEP, bend-sequence
JSON and labels JSON. GPL-3.0 data: downloaded and cached at run time,
**never committed**. Compare our unfold of the folded STEP with their
unfolded STEP on invariants (area, perimeter, hole count, sorted hole
areas) and bend count against the labels. Report success rate and a
histogram of refusal codes. *Verify* first which K and units BenDFM uses.

### FreeCAD oracle

Optional CI job in a FreeCAD + SheetMetal container comparing flat areas on
the fixture set.

### Benchmark (`bench/run.py`)

Times per stage (read, heal, unify, recognize, unfold, flat2d, export).
Scaling series: perforated plate with N = 10 … 5000 holes, to show linear
behaviour. 500 BenDFM parts for p50/p95. Target: recognize + unfold < 50 ms
for ≤ 500 faces. Results stored as a CI JSON artifact; no hard threshold at
first.

### CI

ubuntu / macOS / windows × Python 3.11 and 3.14 × OCP 7.8 and 8.0 (reduced
matrix); ruff, mypy, pytest.

---

## 4. Order of work

1. ~~**Phase-0 spike**~~ **done** ([`spike/README.md`](../spike/README.md)).
   - 11 generated parts are exact to ≤ 3.5e-14.
   - BenDFM samples 0–5 match their own unfolds to ≤ 2.9e-5 (limited by
     sampling).
   - OCP 7.8 and 8.0 give identical results.
   - Typical parts take 10–90 ms including STEP read. 1000 holes take 1.5 s,
     linear after six hotspot fixes.
2. Repo skeleton: `pyproject.toml` (Apache-2.0), CI, ruff/mypy, fixture
   generator.
3. `recognize` → `bendmodel` → `unfold` → `flat2d` → `export` → `cli`, each
   with tests.
4. Benchmark and BenDFM evaluation; v0.1 on PyPI (check the name is free).
5. v0.2: canonical recognition, sharp bends, countersinks, DFM checks.
6. Consumers: Snest calls flatmat from its STEP worker where it currently
   refuses bent parts.

### Scope by version

- **v0.1**: planes + cylindrical bends; one solid of constant thickness;
  seam-split bends; holes and cut-outs crossing bends; hems with r > 0; DXF
  + bend table.
- **v0.2**: canonical recognition of B-spline faces; sharp bends; countersinks
  and counterbores; warn-and-skip instead of abort; DFM checks (§7).
- **later**: cones; formed features beyond detection; tubes with a
  user-chosen seam; assemblies; bend sequence feasibility.

## 5. Hard geometry

Research: [`research.md` §8](research.md#8-hard-geometry). No open tool
handles these cases, so this is flatmat's own work.

### 5.1 Neutral radius per face, not per skin

Joggles show that one skin can be the concave side of one bend and the
convex side of the next. The neutral radius is therefore always derived
from **each cylinder face's own convexity**:

- concave on A: `R_n = R + K·t`
- convex on A: `R_n = R − t + K·t`

This replaces any global "inner skin" assumption in §1.3 and §1.7.

### 5.2 Hems (v0.1 if the inner skin exists, v0.2 otherwise)

- **Detection**: a cylinder with Δu ≥ π − ε.
  - Δu ≈ π: flat or open hem.
  - Δu > π: teardrop or rope.
  - Δu > π with no flat leg after it: rolled.

  The legs are antiparallel planes ≈ 2·r_in apart. `Bend.kind = "hem"`
  with a sub-kind.
- **Thickness pairing uses a material test**: the midpoint between the
  candidate faces must classify IN (`BRepClass3d_SolidClassifier`).
  Otherwise facing inner legs of a closed hem (gap ≈ 0.002 mm) are paired
  as "thickness".
- **Contact**: skin faces that are not adjacent but lie within 0.1·t
  (`BRepExtrema_DistShapeShape`), non-manifold edges, or two lumps give
  `HEM_CONTACT` (a warning if the topology is still manifold, otherwise an
  error).
- **Allowance**:
  - Never use k(r/t) as r → 0, because DIN's k drops to ≤ 0 below
    r/t ≈ 0.05.
  - Use a dedicated `HemTable(t, gap) -> allowance`, Creo style: the flat
    length is the sum of the lip lengths plus the allowance.
  - Fallback: `k_hem` (default 0.43, the shop rule) with a
    `HEM_DEFAULT_ALLOWANCE` warning.
- **Merged-skin hems** (fused r = 0: a 2t planar pair plus a half-cylinder
  R ≈ t plus a t-high step face) give `HEM_MERGED` (error). An opt-in
  reconstruction splits the 2t block at its mid-plane (v0.2+).

### 5.3 Sharp bends (v0.2)

Classify each straight concave edge E between two skin faces A and B that
have thickness partners. Cut-out corner edges of length ≈ t between
thickness faces are excluded. Offset A and B by t; their intersection gives
the outer virtual-sharp line L.

| Case | Outer side | Handling |
|---|---|---|
| (a) sharp/sharp | a convex edge on L | synthetic bend: develop the inner skin and insert a strip along E of width `gap = BA(r_eff) − 2·r_eff·tan(φ/2)`, equivalently `2t·tan(φ/2) − BD(φ)`; with r_eff = 0 that is φ·K·t |
| (b) sharp inside, round outside | a cylinder R ≈ t, axis on E | develop the **outer** skin; the normal cylinder path with r_in = R − t = 0 works unchanged |
| (c) other | — | `VARIABLE_THICKNESS` at that bend |

- `r_eff` comes from `Tooling` (air bending: p·V, see §7). Without tooling
  it is 0, with `SHARP_BEND_ASSUMED_RADIUS` (warning).
- A bend-deduction table, if given, plugs straight into case (a).
- Collinear split edges are grouped into one bend. Sharp hems are handled
  as hems.

### 5.4 Holes, countersinks, formed features (detection v0.2, flat output per type)

For each inner wire of a planar skin face, classify the faces across its
edges:

| Pattern | Feature | Flat output |
|---|---|---|
| one through-cylinder | hole | the circle (exact) |
| coaxial cone (half-angle 41–60°) + cylinder | countersink | through-diameter as the cut, annotation (diameter, angle, side), `COUNTERSINK_PROJECTED` |
| cone through the full thickness | knife-edge countersink | the smaller circle, warning |
| cylinder + annular plane + cylinder | counterbore | through-diameter, annotation |
| tangent torus → coaxial cylinder pair with ΔR ≈ t ending in a free edge | extruded hole | pre-hole diameter from a material/tool table (published formulas disagree), else `FORMED_FEATURE` |
| tangent faces returning to a parallel plane at offset h, no thickness faces | emboss, dimple, bead | marking-layer outline plus tool id, no cut |
| thickness faces inside a formed region | louvre, lance, bridge | slit lines plus marker plus tool id |

All of these become `SheetPart.features` (type, faces, parameters, side).
Unrecognised formed regions give `FORMED_FEATURE` (an error by default, a
warning with `partial=True`).

### 5.5 B-spline faces (canonical recognition, v0.2)

- Recognition is **per face** with `GeomConvert_SurfToAnaSurf(surface)`,
  calling `ConvertToAnalytical(tol, *BRepTools.UVBounds_s(face))` and then
  `Gap()`.
  - Do not use the shell mode of `ShapeAnalysis_CanonicalRecognition`: in
    7.8.0 and 8.0.0 it compares the first face with itself, so coaxiality
    is never checked. We compare axes ourselves.
  - Note that `SurfToAnaSurf` clips the tolerance to bbox diagonal / 1000.
- The converted surface has a new parameterisation, so the original
  pcurves are useless. For cylinders and cones we compute (u, v) ourselves
  from the analytic axes (exact, unwrapping u) instead of rebuilding faces.
- Faces that are not developable (Gaussian curvature ≠ 0 on a grid via
  `GeomLProp_SLProps`) give `UNSUPPORTED_SURFACE` (non-developable).
- A converted face carries `SURFACE_CONVERTED` (the gap). A converted
  cylinder whose sagitta is below the unify tolerance is protected with
  `KeepShapes` (§6.1).

### 5.6 Cones (later)

- **Development**:
  - OCCT `P(u,v) = O + (R + v·sinα)(cos u·X + sin u·Y) + v·cosα·Z`;
  - flat `ρ = v + R/sinα`, `ψ = u·sinα`, point `(ρ cosψ, ρ sinψ)`;
  - sector angle Δu·|sinα|.
- **Neutral layer**: shift ρ by `K·t·cotα` from the concave skin; ψ is
  unchanged.
- Not an affine u-scale: u = const maps to rays and v = const to arcs.
  Other pcurves are sampled and refit (`Geom2dAPI_PointsToBSpline`) with a
  stated tolerance.
- K is constant over the cone (as SolidWorks does).
- **Refuse**: inner and outer cones sharing an apex (thickness not
  constant), and inexact or oblique cones.

### 5.7 Joggles, curved flanges, rolled parts, tubes

- **Joggle**:
  - two parallel-axis cylinders of opposite convexity, adjacent or
    separated by a short flat;
  - works with §5.1 automatically;
  - `Bend.kind = "joggle"` groups the pair;
  - a `JoggleTable(offset)` (Creo-style `offset_allow`) can override the
    allowance.
- **Curved flange** (a flange along a curved edge, i.e. a torus bend): not
  developable, so `UNSUPPORTED_SURFACE`. A flange on a straight generatrix
  of a rolled base is fine (all cylinders).
- **Rolled part vs tube** (Quaoar part 8): Σ Δu ≈ 2π across the cylinder
  chain.
  - Genus 0 (a rolled part with its slit modelled): fine. Seam edges are
    never glued (they are already excluded in §1.1).
  - Genus 1 (a closed tube): `CYCLE_INCONSISTENT` unless `cut_edge=` is
    given. Then the output is a rectangle 2π(R_in + K·t) wide with the seam
    marked.

### 5.8 Reliefs

- Relief side faces are thickness planes perpendicular to the bend axis.
  They develop to straight lines perpendicular to the bend line.
- **Bend line extent = the trimmed v-range of the bend face.** Notched bend
  faces from partial reliefs are developed exactly from their pcurves.
- Pitfalls handled explicitly:
  - sliver faces where a relief barely clips a bend (area below a tolerance
    gives `SLIVER_FACE`, a warning);
  - closed corners where flanges touch (the contact check from §5.2).
- A missing relief (the bend runs out with material on both sides) is a DFM
  finding, `DFM_RELIEF_MISSING` (§7), not a recognition error.

## 6. Robustness on real-world STEP files

Research: [`research.md` §9](research.md#9-robustness-on-real-world-step-files).
Principle (Quaoar): **detect and report, don't guess.** Healing that changes
topology in ways we cannot trace is worse than a clear refusal.

### 6.1 Import pipeline

```
read ─▶ (heal: full | fast | off) ─▶ validate ─▶ canonicalize ─▶ unify ─▶ validate ─▶ recognize
```

1. **Read**:
   - Plain parts use `STEPControl_Reader`; assemblies use
     `STEPCAFControl_Reader` (§6.4).
   - Set the unit to mm before the transfer (`xstep.cascade.unit`;
     `XCAFDoc_DocumentTool::SetLengthUnit` for XCAF).
   - Set `read.surfacecurve.mode=3`, so missing pcurves are rebuilt from 3D.
   - Set `read.maxprecision.val = 0.1`; the default is 1 mm, which can
     swallow thin features.
   - Set `read.step.tessellated = 0`. Mesh-only faces mean
     `UNSUPPORTED_SURFACE`.
2. **Healing modes** (*spike*, measured on an L with 1000 holes):
   - **The reader's built-in healing runs `ShapeFix_Shell`, which is
     quadratic in the face count:** 4.1 s of a 4.4 s transfer. Face- and
     wire-level fixes on the same part take only ~130 ms.
   - `fast` (default): the reader's `FixShape` with `FixShellMode = 0` and
     `FixSolidMode = 0`; orientation is checked by us (volume > 0). Transfer
     drops to **0.18 s (8.0) / 0.25 s (7.8)** and the solid stays valid.
     - **OCP 8.0**: `reader.SetShapeFixParameters(map)` **after**
       `ReadFile`, because the actor that stores the map only exists then.
       Pass a complete key→value map: seed it with
       `SetShapeFixParameters(DE_ShapeFixParameters())`, read it back,
       `Bind` the overrides, set it again. The overload with "additional
       parameters" does not override existing keys.
       `SetShapeProcessFlags` cannot be called from Python, because
       `std::bitset` is not bound.
     - **OCP 7.8**: a resource file named `STEP` in a directory pointed to
       by `$CSF_STEPDefaults`, set before the first read in the process,
       with `FromSTEP.exec.op : FixShape` and the two modes at 0.
   - `full`: the reader's default healing. On 7.8 that needs a process
     without the resource override, so the worker pool keeps one 7.8
     process per mode.
   - **`FixShellMode = 0` also skips all face- and wire-level fixes**
     (*spike*, on a real SolidWorks file). So `fast` must add missing
     pcurves itself: `ShapeFix_Wire(w, f, tol).FixEdgeCurves()` on every
     wire with an edge whose `BRep_Tool.CurveOnSurface` is null. It works
     in place and took a few ms. SolidWorks 2023 writes **no pcurves for
     cylindrical faces at all**.
   - Before touching any pcurve, check it exists: `BRepAdaptor_Curve2d` on
     a missing pcurve raises on 7.8 but **segfaults on 8.0**, which would
     take down the worker process.
   - Fallback: if any later stage fails on a `fast` import (outline not
     closed, area identity broken), re-import with `full` once.
3. **Validate**:
   - **No `BRepCheck_Analyzer` by default** (*spike*): on a face with 1000
     holes it costs 1.3 s on 8.0 and **13 s** on 7.8, even topology-only.
     It runs only as a diagnostic when a later stage fails, or with
     `strict=True`. The self-checks (§1.5, §1.6) are cheaper and more to
     the point.
   - `ShapeAnalysis_ShapeTolerance`: tolerance above 5 % of t means
     `DIRTY_TOLERANCE` (a warning, or an error if `strict`).
   - Volume above 0, otherwise `BRepLib::OrientClosedSolid`, then
     `INSIDE_OUT_FIXED` (a warning).
   - A null surface means `UNSUPPORTED_SURFACE`.
4. **Canonicalize** (v0.2):
   - `ShapeAnalysis_CanonicalRecognition` (present in 7.8): `IsPlane`,
     `IsCylinder`, `IsCone`.
   - `GeomLib_IsPlanarSurface`.
   - An extrusion of a circle is a cylinder.
   - Each conversion gets a `SURFACE_CONVERTED` warning with the deviation.
5. **Unify**:
   - `ShapeUpgrade_UnifySameDomain` with UnifyFaces, UnifyEdges,
     `AllowInternalEdges=false`, and a small `SetLinearTolerance`.
   - **Every flange–bend edge goes into `KeepShapes`.** Otherwise a flat
     B-spline bend whose sagitta is below the tolerance would be merged
     into the flange.
   - The `History()` is kept, so diagnostics can name the original face
     ids.

### 6.2 Tolerances (`Tolerances.adaptive(shape)`)

- `linear = clamp(max(1e-4, 1e-6 · bbox_diag, 2 · max entity tolerance),
  upper = 0.05 · t)`. Hitting the clamp means `DIRTY_TOLERANCE`, never a
  wider tolerance.
- Thickness bins of `max(linear, 0.005 · t)`, using the area-weighted mode.
  It is accepted when ≥ 95 % of the skin area lies within `±max(0.01 mm,
  0.01 · t)`.
- Edge classes:
  - **SMOOTH** below 0.1°; **SHARP** above 1°.
  - **NEAR_TANGENT** in between, resolved by an exact surface-pair test.
    For a plane and a cylinder: |dist(axis, plane) − R| < linear and axis ⟂
    normal. For B-splines: normals at 3–5 points along the edge.
  - **SAME_DOMAIN**: same underlying surface (arbitrary splits); merged or
    grouped.
- A genuine small kink (e.g. 0.5°) is **not** snapped to tangent, because
  that would mean 8.7 mm lateral error per metre of flange. It is developed
  as a zero-radius bend with a `NEAR_TANGENT_KINK` warning.

### 6.3 Defect detection

| Defect | Detection | Result |
|---|---|---|
| misaligned halves (Quaoar 0.24 mm case) | thickness outliers plus short step faces next to a bend | `MISALIGNED_BODIES`, with the faces |
| bend without a coaxial partner, inner = outer radius | cylinder pairing in §1.2 fails for that bend | `BEND_UNPAIRED` |
| welded or fused flanges ("fake origami") | a thickness component touching more than two skin regions, or a cycle in the bend graph | `CYCLE_INCONSISTENT` / `NOT_TWO_SKINS` |
| out-of-period UV on cylinders | normalise every UV into `[umin, umax]` before developing | — |

A partial result is never returned silently. With `partial=True` the
caller gets a `FlatPattern` with `complete=False` and the unidentified faces
listed.

### 6.4 Assemblies, multiple bodies, quantities

- `fm.load_step_assembly(path)` returns a list of `PartInstance(name, path_id,
  location, solid)`.
  - It walks from `XCAFDoc_ShapeTool.GetFreeShapes`, following references,
    components and locations, with names from `TDataStd_Name`.
  - Reference implementation: build123d's `importers.py`, or
    `XCAFPrs_DocumentExplorer` (`OnlyLeafNodes`).
- `fm.group_identical(instances)` returns `[(prototype, quantity,
  instances)]`.
  - It groups first by referred label, then by the geometric hash
    `V:E:F:Sh:So:area:volume:OBB` with a handedness suffix (`h+`, `h−`,
    `h0`).
  - **Mirrored parts stay separate**, because their flats are mirrored.
- A hidden flat-pattern body is detected heuristically: planar faces only,
  two large antiparallel faces at distance t, the same t, and a large-face
  area roughly equal to the developed area, and/or a name matching
  `flat|abwicklung|unfold`. It is reported as a `reference_flat` and used
  as a free cross-check.

### 6.5 Determinism

- Stable ids come from `TopExp::MapShapes` index order. Never iterate
  hash-keyed maps or sets for anything that affects output.
- Root choice: the largest area, ties broken by the rounded centroid.
  Orientation: minimum-area rectangle with a fixed tie-break, bbox minimum
  at the origin.
- DXF:
  - coordinates rounded to 1e-6 mm, −0.0 → 0, arcs normalised to CCW;
  - entities sorted by (layer, type, coordinates);
  - `ezdxf.options.write_fixed_meta_data_for_testing` in tests.
- Output metadata records the flatmat, OCP and OCCT versions, the healing
  mode and the tolerances. Golden-file tests run per pinned OCP version.

## 7. Manufacturability checks (`flatmat/dfm.py`, v0.2)

Research: [`research.md` §10](research.md#10-manufacturability-rules-for-laser--press-brake).
Sources disagree, so **every threshold is a parameter**. The defaults below
are a reasonable middle and are documented with their source.

### 7.1 Inputs

- `FlatPattern` and `SheetPart` (bends, flanges, holes, bend zones).
- `Material(name, rm_mpa, family, min_r_over_t, air_radius_pct_of_v,
  hem_allowed)`. A small built-in table: mild steel, stainless 304,
  aluminium 5052-H32 and 6061-T6, galvanised.
- `Tooling` (optional), a parametric schema:
  - `Punch(height, tip_angle, tip_radius, gooseneck_throat_depth,
    gooseneck_throat_height)`
  - `Die(v, angle, shoulder_radius, width)`
  - Default: BenDFM's punch (10 mm thick, 90° tip, 300 mm long) and a V
    chosen by the Amada rule below.
- `rolling_direction` (optional vector in the flat).

### 7.2 Derived quantities

- **V when no die is given** (Amada): V = 6t (t ≤ 2.5 mm), 8t (3–8 mm), 10t
  (9–10 mm), 12t (≥ 12 mm).
- **Expected air-bend radius:** r_air ≈ p·V, with p = 0.16 (mild steel),
  0.21 (stainless), 0.12 (aluminium). If the model's r_in differs by more
  than a tolerance, warn `DFM_RADIUS_NOT_AIR_BENDABLE`: the shop's radius,
  and therefore the K-factor and flat size, will differ from the model.
- **Tonnage:** `F[kN] = 1.42 · Rm[MPa] · t² · L / (1000 · V)` per bend (mm),
  × 1 for air bending. Reported, and checked against `Machine.max_force`
  when one is given.

### 7.3 Rules

| Code | Check | Default threshold |
|---|---|---|
| `DFM_FLANGE_TOO_SHORT` | each flange length measured from the bend (outer mould line), over ≥ 50 % of the bend length | `max(0.7·V, 2.5·t + r)` |
| `DFM_RADIUS_TOO_SMALL` | r_in vs the material's minimum radius; the parallel-to-grain value when the bend line is within 45° of the rolling direction or the direction is unknown | material table |
| `DFM_HOLE_NEAR_BEND` | distance from a hole or slot edge to the bend tangent line and to the centreline | `max(2.5·t + r from tangent, V/2 from centreline)`; slots `4·t + r` |
| `DFM_CUT_IN_BEND_ZONE` | a cut-out inside the keep-out band that does not cross the whole bend | band = BA strip ± `max(2.5·t + r, V/2 − BA/2)` |
| `DFM_HOLE_TOO_SMALL` | hole diameter | ≥ t (cut by laser: ≥ 0.5·t, configurable) |
| `DFM_HOLE_NEAR_EDGE` | hole to the outer contour | ≥ 2·t |
| `DFM_HOLES_TOO_CLOSE` | web between holes | ≥ 2·t |
| `DFM_FEATURE_TOO_SMALL` | narrowest gap or slot in the outline | ≥ 1.5·t (laser) |
| `DFM_RELIEF_MISSING` | a bend end where material continues past the end without a cut-out ≥ t wide and ≥ r + t deep past the tangent | relief width t, depth r + t |
| `DFM_CORNER_GAP` | two bend zones meeting at a corner | gap ≥ 0.4 mm or a corner relief |
| `DFM_HEM_TOO_SHORT` | hem return length | ≥ 4·t |
| `DFM_HEM_MATERIAL` | hem on a material or thickness where it is not allowed | material table, Protocase maximum thicknesses |
| `DFM_U_CHANNEL` | U with base < 2 × flange height, or flange taller than the punch clearance or gooseneck throat | — |
| `DFM_TOOL_COLLISION` | 2D section collision (§7.4) | — |
| `DFM_BEND_ANGLE` | angle beyond what the shop forms in one bend | 130° (SendCutSend), configurable |

Every finding is a `Diagnostic` with severity `warning` (or `error` if the
caller sets `strict=True`), pointing to the bends, faces and 2D locations
involved, so a UI can draw it.

### 7.4 Collision and sequence (later)

- **Section test:** for each bend, cut a section perpendicular to the bend
  line. Represent the punch (straight or gooseneck) and the die as
  polygons. Sweep the bend angle in 5° steps with the already-formed
  flanges attached, testing at the left, centre and right of the bend (as
  BenDFM does).
- **Sequence hint** (Duflou 1999):
  - search backwards from the folded part
  - collisions become precedence constraints
  - heuristics, by weight:
    - shape-defining bends last (4)
    - leaves of the connectivity graph first (2)
    - shorter bends first (1)
  - collinear bends are grouped
- **Output:** a feasible order, or the bend that cannot be formed, with
  the reason.

## OCP binding pitfalls (*spike*)

- **Iterating an OCP shape list from Python costs ~0.3 ms per list** (7.8
  and 8.0), about 1 s for the edge→face table of a 3000-edge part. Use
  `First()`/`Last()`/`Size()`; an edge has at most two faces.
- C++ `double&` out-parameters do not come back (e.g. `First`/`Last` of
  `BRep_Tool.CurveOnSurface`). Use `BRepAdaptor_Curve2d`.
- Cache one `BRepAdaptor_Surface` per face. Rebuilding it per edge made
  edge classification take 9 s on a 1000-hole part.
- 7.8 vs 8.0 names: maps and lists moved to `OCP.collections` with long
  template names; `TopoDS.Face_s` became `TopoDS.Face`;
  `BRepTopAdaptor_FClass2d` is gone. A ~20-line `_occ.py` shim covers it
  (see `spike/occ_compat.py`).

## Open points

- ~~PyPI name~~: `flatmat` is free (404 on 2026-10-11).
- ~~DIN 6935 range~~: resolved (§1.7).
- ~~BenDFM K~~: K = 0.45 (ANSI) for every bend checked, in mm, true
  cylinders, made with PythonOCC. Confirmed by the spike reproducing their
  flat areas.
- ~~OCP 7.8 vs 8.0~~: resolved (see the pitfalls above and `spike/README.md`).
- ~~Fast healing on 7.8~~: resolved via `$CSF_STEPDefaults` (§6.1).
- Is OCCT 8.0's STEP reader thread-safe? (Irrelevant while we use one
  process per part.)
- Hem allowance and extruded-hole pre-hole tables: no open data found.
- Real STEP exports: one SolidWorks 2023 part tried (missing cylinder
  pcurves, now handled). Inventor, Fusion, Onshape, Creo and Solid Edge are
  still to collect.
- `ShapeUpgrade_UnifySameDomain` with `KeepShapes` and canonical recognition
  are untested.
