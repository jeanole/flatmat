# flatmat — design

Status: **v0.1 implemented, v0.2 in progress** (see the
[README](../README.md)). This document began as the design. The phase-0
spike ([`spike/`](../spike/README.md)) ran §1 end to end on OCP 7.8 and 8.0;
its results are folded in (marked *spike*). Where the code differs from the
plan, the section says so (*as built*); the modules are listed in
[§4](#code-layout-as-built). Anything still marked *verify* has not been
checked yet.

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

### 1.6 Flat pattern (`flatmat/unfold.py`)

*As built*: there is no separate `flat2d.py`; the outline, the bend lines and
both checks are part of `unfold.py`. The overlap check samples the loops at
0.005 mm and tests them with shapely; `AREA_MISMATCH` uses
`Tolerances.area_rel` (1e-4 by default).

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
| warnings | `K_EXTRAPOLATED`, `COUNTERSINK_PROJECTED`, `ENGRAVED_RECESS`, `SHARP_BEND_ASSUMED_RADIUS`, `SURFACE_CONVERTED`, `NEAR_TANGENT_KINK`, `HEM_CONTACT`, `HEM_DEFAULT_ALLOWANCE`, `INSIDE_OUT_FIXED`, `SLIVER_FACE` |
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

### Fixture generator (`tests/parts.py`, OCP only, no licence issues)

A 2D profile (flat lengths, angles, r, t) with inner fillet r and outer
fillet r + t, extruded by a width. Shapes: L, U, Z, hat, hem, joggle, box
with corner reliefs. Variants: holes in a flange, **a hole through a bend**,
a cut-out splitting a bend (grouping), seam-split cylinders
(`ShapeUpgrade_FaceDivide`), NURBS-converted (`BRepBuilderAPI_NurbsConvert`,
v0.2), each round-tripped through STEP.

*As built*: `tests/parts.py` builds L, U, Z, hat, obtuse, acute, open hem
and joggle profiles, holes in a flange and through a bend, sharp bends,
countersinks and counterbores, embosses and extruded holes, louvres, bridges
and lances, one part per DFM rule (cut-outs as boxes), and boxes with and
without corner reliefs (`box_solid`). The NURBS variants are made in
`tests/test_canonical.py`, the engraved recesses in `tests/test_engrave.py`.
Not built yet: a cut-out splitting a bend, and seam-split cylinders.

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

### BenDFM (planned: `bench/bendfm.py`)

20 000 synthetic bent parts with folded STEP, unfolded STEP, bend-sequence
JSON and labels JSON. GPL-3.0 data: downloaded and cached at run time,
**never committed**. Compare our unfold of the folded STEP with their
unfolded STEP on invariants (area, perimeter, hole count, sorted hole
areas) and bend count against the labels. Report success rate and a
histogram of refusal codes. *Verify* first which K and units BenDFM uses.

*As built*: only the six example parts in the BenDFM repository are used, by
`tests/test_external.py` when `FLATMAT_BENDFM_DIR` points at them (CI
downloads them). It checks thickness, bend count and flat area against their
unfolded STEP, with K = 0.45. The evaluation on the full 20 000 parts is not
built yet.

### FreeCAD oracle

Optional CI job in a FreeCAD + SheetMetal container comparing flat areas on
the fixture set.

### Benchmark (planned: `bench/run.py`)

Times per stage (read, heal, unify, recognize, unfold, flat2d, export).
Scaling series: perforated plate with N = 10 … 5000 holes, to show linear
behaviour. 500 BenDFM parts for p50/p95. Target: recognize + unfold < 50 ms
for ≤ 500 faces. Results stored as a CI JSON artifact; no hard threshold at
first.

*As built*: not built yet. The phase-0 spike measured the scaling series
with `spike/run.py`.

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
2. ~~Repo skeleton~~ **done**: `pyproject.toml` (Apache-2.0), CI matrix
   (Linux/macOS/Windows × OCP 7.8/8.0), ruff, mypy, the fixture generator.
3. ~~`recognize` → `bendmodel` → `unfold` → `export` → `cli`~~ **done for
   v0.1**, in `src/flatmat/`. Changes against the plan above:
   - the 2D outline lives in `unfold.py` (no separate `flat2d.py`);
   - loop areas are exact (Green's theorem per line and arc; splines by
     Gauss quadrature per knot span);
   - collinear lines and co-circular arcs are joined (SolidWorks splits
     circles into halves);
   - periodic B-spline pcurves are clamped with `SetNotPeriodic()` before
     export (DXF needs a clamped knot vector);
   - bend-line end points are found by bisection on the face classifier
     (to ~1e-9 of the bend length).

   Tests: 77 in total, including hypothesis profiles, the BenDFM samples
   and private parts via environment variables.
4. Benchmark and BenDFM evaluation; v0.1 on PyPI (check the name is free).
5. v0.2: canonical recognition, sharp bends, countersinks, DFM checks;
   embosses and extruded holes, the corner rule, tooling, collisions and the
   bend order (done).
6. Consumers: Snest calls flatmat from its STEP worker where it currently
   refuses bent parts.

### Code layout (as built)

All in `src/flatmat/`:

| module | contents |
|---|---|
| `io.py` | STEP import and healing (§1.0, §6.1) → `StepModel` |
| `canonical.py` | B-spline faces and edges that are really planes, cylinders, cones, lines or circles → analytic (§5.5) |
| `recognize.py` | face graph, skins, thickness, sharp bends, bends, flanges, hole and formed features, engravings (§1.1–1.3, §5.3–5.4) → `SheetPart` |
| `_topology.py` | internal: the OCP-side context `recognize()` hands to `unfold()` (faces, bend faces, sharp edges, pcurves, feature data) |
| `unfold.py` | per-face development, placement, outline, bend lines, self-checks (§1.4–1.6) → `FlatPattern` |
| `bendmodel.py` | BA, BD, OSSB, K models, DIN 6935, extruded-hole pre-holes (§1.7) |
| `dfm.py` | manufacturability checks, one function per rule group (§7) → `DfmReport` |
| `tooling.py` | punches, dies and the die-opening rule (§7.5) |
| `fold.py`, `sequence.py` | the folded part, tool collisions and a bend order (§7.4) → `BendPlan` |
| `export.py` | DXF, bend table (JSON/CSV), SVG (§2) |
| `model.py`, `diagnostics.py` | the public dataclasses, `FEATURE_KINDS` (how each stage treats each kind of feature) and the diagnostic codes |
| `_occ.py` | the differences between OCP 7.8 and 8.0 |
| `cli.py` | `flatmat unfold` |

### Scope by version

- **v0.1**: planes + cylindrical bends; one solid of constant thickness;
  seam-split bends; holes and cut-outs crossing bends; hems with r > 0; DXF
  + bend table.
- **v0.2**: canonical recognition of B-spline faces; sharp bends; countersinks
  and counterbores; embosses, extruded holes, louvres, bridges and lances
  (§5.4); warn-and-skip instead
  of abort; DFM checks (§7) with the corner rule, tooling (§7.5), tool
  collisions and a bend order (§7.4).
- **later**: cones; tubes with a user-chosen seam; assemblies.

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

### 5.3 Sharp bends (implemented)

**As built:**

- The thickness is estimated **before** the skins are known
  (`_plane_partners`): planes are bucketed by normal, and each plane gets as
  partner the nearest antiparallel plane *behind* it whose face contains the
  foot of an inner point. A sharp bend splits a skin into separate smooth
  components, so the skins can only be found after the sharp bends.
- A sharp-bend edge is a straight, non-tangent edge between two partnered
  planes whose partners either meet at a sharp edge (sharp on both sides) or
  touch a common cylinder of radius ≈ t (sharp inside, round outside). Skins
  are the components over smooth edges *and* sharp-bend edges.
- Placement across a sharp edge moves the child face perpendicular to the
  edge by `Δ = BA(r) − 2·s·tan(φ/2)`. Here `s = r` if the reference skin is
  the inside of the corner and `s = r + t` if it is the outside, with
  `r = sharp_radius` (default 0, which raises the warning
  `SHARP_BEND_ASSUMED_RADIUS`).
- The outline is closed at the two vertices of the edge (`_close`):
  - for `Δ ≥ 0`, a bridging line;
  - for `Δ < 0`, both outline lines are trimmed to their intersection (or
    joined if collinear);
  - anything other than straight lines there gives `SHARP_BEND_UNSUPPORTED`.

  Any other gap in the outline is now an error (`OUTLINE_OPEN`).
- **Sharp inside, round outside, with the reference skin outside:** the
  cylinder carries the corner. Its neutral radius is
  `(BA(r) − 2·r·tan(φ/2)) / φ`, so a `sharp_radius` applies there too.
- Bend lines lie in the middle of the inserted strip. The bend table reports
  the radius actually used.
- Tests:
  - 7 generated parts (L, U, Z, obtuse, mixed; sharp inside or on both
    sides) exact to 1e-9;
  - random profiles that mix round, sharp-inside and sharp-both bends;
  - the bend-line position and the `sharp_radius` formula.
- **A bug found on the way:** canonicalising a direction's sign *before*
  rounding let −1e-17 noise split parallel planes into different buckets.
  Now it rounds first (`_canonical_dir`).

**Original plan:**

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

### 5.4 Holes, countersinks, formed features

**Countersinks and counterbores: implemented.**

- **Detection:** non-skin faces are grouped over *all* edges. A group of
  coaxial cylinders, cones and annular planes that touches both skins is a
  hole feature:
  - with a cone → countersink (included angle from the cone);
  - with two or more cylinder radii → counterbore (depth from the large
    circles).
- **Flat output:** the through (smallest) circle is cut. Where the large rim
  is on the reference skin, its loop is replaced by the through circle
  (`COUNTERSINK_PROJECTED` for countersinks).
- **Reporting:** the feature appears in `FlatPattern.features` with its flat
  centre, in the bend table (`features`) and as DXF text (`CSK 90° Ø10 TOP`).
- **Unchanged:** groups that are not hole features still go through the
  per-component rule, so a chamfered outer edge is still `EDGE_FEATURE`.

**Embosses and extruded holes: implemented** (`recognize.py`, `unfold.py`).

- **Detection, per skin:**
  - *Seeds* are the skin's tori, spheres, cones and freeform faces. They grow
    over cylinders tangent to them (an emboss's wall, a collar's bore).
  - Removing the seeds leaves the main sheet (the largest part) and what the
    seeds enclose (an emboss's raised top).
  - A region is accepted only if it meets the main sheet along inner wires of
    one plane face (the host). Otherwise `FORMED_FEATURE` names its faces.
    Cones alone are not seeds, so conical bends still give
    `UNSUPPORTED_SURFACE`.
- **Classification:**
  - no edges to thickness faces → **emboss** (dimple, bead, raised pad);
  - edges to thickness faces, and every face a torus or cylinder coaxial on
    an axis normal to the host → **extruded hole** (Durchzug);
  - anything else with cuts (louvre, lance, bridge) → `FORMED_FEATURE`.
- **Measures:**
  - the height comes from the region's vertices and inner sample points,
    along the host's normal; the side is the side it stands out on;
  - for a collar: bore, outside height (from the far side of the sheet to
    the rim) and the fillet's inner radius.
- **Pre-hole** (`bendmodel.extruded_hole_prehole`):
  - laid flat, the neutral fibre from the rim of the formed region runs
    through the fillet and the wall: `d = 2·(ρ − θ·(r + K·t) − h)`;
  - with K = 0.5 and a 90° fillet this is exactly the published
    `d = D − 2·(H − 0.43·r − 0.72·t)` (tested);
  - circumferential stretching of the rim is ignored. `unfold(pre_hole=…)`
    takes a diameter or `(Feature) -> diameter` for the tool maker's table;
  - no material left → `PRE_HOLE` warning, pierced (no cut).
- **Flat output:**
  - the formed faces leave the reference skin. Their rim becomes a loop of
    the host face;
  - an emboss's rim is not cut. The area inside it is material, and it is
    returned as `Feature.outline`;
  - an extruded hole cuts its pre-hole; its rim is the outline;
  - the area self-check adds what lies inside the rims;
  - DXF: rims on layer `FORM`, text `EMBOSS H5 TOP` or
    `EXTR Ø6 H4.5 BOTTOM (PRE Ø2.708)`. The bend table has `height`,
    `radius` and `rim`.
- **Sharp-edged collars** (no fillet) are not supported. A counterbore
  deeper than t is not taken for one; they end as `EDGE_FEATURE`.
- **Tests:** generated by revolving the cross-section: raised pads up and
  down, collars up, down and two on one part, a collar too tall for a
  pre-hole. All are exact in area and round-trip through STEP.

**Louvres, bridges and lances: implemented**, bent or formed, with straight or
curved slits.

- **Detection** (`_formed_regions`, `_classify`): the opening's inner wire on
  the host face is walked in order (`BRepTools_WireExplorer`). Each edge is
  *attached* (to the region) or *cut* (to a thickness face); runs of each are
  counted.
  - **Formed** (seeded, as for embosses: tori, spheres, freeform faces):
    - held on one side and cut along one run: a **louvre** (Kieme) if the
      cut is a straight slit (within t of its chord) or takes at most half the
      rim, else a **lance** (Lasche: a tongue cut free around most of its
      rim);
    - held on two sides and cut along two runs: a **bridge** (Brücke).
  - **Bent** (planes and cylinders only, hanging in an inner wire of a plane
    face and touching nothing else):
    - held on one side: a **lance**, developed like any flange;
    - held on two sides: a **bridge** with bent ramps (its development
      would not close).
  - Anything else (held on more sides, cuts inside the region) →
    `FORMED_FEATURE`.
- **Slit width** (`_slit_width`):
  - the region's own cut edges (its edges along thickness faces) are
    projected onto the host plane;
  - the median distance to the sheet's cut edges is the width. The ends,
    where both meet, do not count;
  - in 3D the two sides only lie in one plane where the slit leaves the
    sheet, so a plain 3D distance would give 0.
- **Flat output:**
  - a formed louvre, bridge or lance is not developed, as for an emboss.
    Its rim is a marking (`Feature.outline`, layer `FORM`) and the inside is
    material.
  - **Slits with width:** each slit is the sheet's cut edge (exact: lines
    and arcs of the rim) closed by the region's cut edge projected onto the
    sheet.
    - Projected edges that are lines or circles (a cut face normal to the
      sheet) become exact `Line`/`Arc` segments; others stay a polyline
      within 1e-4 mm.
    - If no projected edge joins a run's ends, the rim and the projected
      edges are polygonized, and the faces along the run are the slits.
    - The projection is where the region's edge lies once formed; for a cut
      face normal to the sheet, that is also where the blank is cut.
  - **Slits without width:** cut lines in `FlatPattern.slits` (layer
    `INNER`); the fold model and the checks treat them as cuts.
  - **Bent lance:** the tab is developed like any flange
    (`Feature.bend_id`), and the opening around it is cut as it is.
  - `Feature.holes` lists the holes each feature cuts (also pre-holes and
    countersinks). `orientation_deg` is the direction of a bent lance's bend
    line, otherwise of the chord of the longest slit. `size` is the extent
    along and across it. DXF text: `LOUVRE 20x10.5 A0 H5 TOP`.
- **Checks:**
  - their holes are judged with the feature, not as holes or narrow slots:
    rim (a bent lance: its opening) to the outline ≥ 3t, to other bends
    ≥ 4t + r;
  - a bent lance's flange is measured to its opening;
  - lances are formed with the blank in a punching tool, not on the press
    brake. A bent lance's bend gets no press-brake rules (flange to V, air
    radius, tooling, tonnage, U channel), and the bend order treats it like
    a hem (made first, `BendPlan.unchecked`). Air-bending a tab inside a
    frame would push the frame into the die, and the simulation shows
    exactly that.
- **Tests** (`tests/test_cut_forms.py`), all slits 0.5 mm wide and
  round-tripped through STEP:
  - a raised half pad open along its diameter (straight louvre, up and
    down);
  - the hood of a raised round pad outside a circle, slit along it (curved
    louvre; the slot is four exact arcs);
  - a tongue and a strip cut from a raised round pad (formed lance, formed
    bridge);
  - a strip raised by two joggles at each end (bent bridge);
  - tabs bent 30° and 90° out of openings as long as their flat (bent
    lance);
  - flat areas exact to 1e-9 for straight cuts and to 1e-5 for curved ones
    against a polygonised reference;
  - widths, sizes, directions and heights as built;
  - a synthetic slit without width through DXF, the fold model and the
    checks.
- **Not supported:** cuts inside a formed region; regions held on more than
  two sides of their opening. The material a formed feature stretches is
  not accounted for: inside its rim the blank is flat.

**Engravings: implemented, behind `recognize(engrave=True)`** (`recognize.py`,
`unfold.py`). A recess that does not go through the sheet — a pocket, a
groove, engraved text, a blind hole — is machined on a mill; on a laser the
shop marks its rim instead. Without the switch such a part is refused as
before (`EDGE_FEATURE`).

- **Detection** (`_recess`): a group of non-skin faces (grouped over all
  edges, as for hole features) that touches only one skin is a recess when
  - it opens into a single plane face of that skin (the host), through inner
    wires of that face only: it reaches no edge of the part and no bend;
  - every edge of those wires is its own (none is shared with a through
    cut);
  - it stays less than t below the host (vertices and inner sample points,
    along the host's normal); that is its `depth`.
  Faces of the group flush with the host are islands its floor leaves
  standing (the middle of an O); their wires are rims too.
- **Features:** one `engraving` per rim, so a ring groove is two. `side` is
  the skin it opens into, `outer_diameter` the rim's diameter if it is one
  circle, `depth` the recess's. One `ENGRAVED_RECESS` warning per part
  counts them and gives the deepest: the depth is not made.
- **Flat output:** the rim is `Feature.outline`, never a hole: the flat is
  that of the part without the recess (area and size tested exact).
  - a rim in a face of the reference skin is a loop of that face and is
    chained with the outline, exactly (lines, arcs, splines); the area
    self-check adds what lies inside it, as for an emboss;
  - a rim on the far skin, or round an island, is projected along the
    normal onto the plane of the reference skin behind it
    (`_projected_rim`): circles stay circles, other edges are fitted as a
    line or an arc, else polylined at 1e-4 mm;
  - DXF: rims on layer `ENGRAVE`, no text per rim (text is many rims). The
    bend table has `depth` and `rim`.
- **Not supported:** a recess reaching the outline (a rebate) or running
  into a bend: `EDGE_FEATURE`, as without the switch. A recess so large that
  its floor outweighs the skin it is cut into is read as a skin.

### 5.5 B-spline faces (canonical recognition, implemented)

**As built** (`src/flatmat/canonical.py`, called from `recognize()`; switch off
with `canonical=False`):

- **Edges:** freeform edges are tested with
  `ShapeAnalysis_CanonicalRecognition(edge).IsLine` / `.IsCircle`; these
  fill a mutable `gp_Lin` / `gp_Circ`. They are replaced by
  `Geom_Line` / `Geom_Circle` with the original parameter direction (the
  circle's axis is flipped if needed), `Range` is set from the vertices,
  and the edge's pcurves are dropped. `GeomConvert_CurveToAnaCurve` was
  not usable from Python: it returns the new curve and range through
  `Handle&` / `double&` out-parameters, and OCP does not hand those back.
- **Faces:** `GeomConvert_SurfToAnaSurf` per face, with the face's UV bounds
  first. *Spike finding:* with bounds the conversion sometimes fails on
  NURBS cylinders, so it retries without them. The target is a plane,
  cylinder or cone. The new surface is `UReverse()`d if its normal
  disagrees with the old one, so the wire orientation stays valid. The
  face is updated in place with `BRep_Builder.UpdateFace`, and the pcurves
  of its edges are removed (`UpdateEdge(e, None, f, tol)` works from OCP).
- **pcurves** are then rebuilt by `ShapeFix_Wire.FixEdgeCurves`.
- It always works on a `BRepBuilderAPI_Copy`, never on the input. A single
  `SURFACE_CONVERTED` warning reports the counts and the largest deviation.
- **Tests:** all 20 test parts (round bends, sharp bends, countersinks and
  counterbores), written as pure NURBS with `BRepBuilderAPI_NurbsConvert`
  and round-tripped through STEP, unfold to the same exact areas (≤ 1e-14).
  Their outlines come back as `LINE` / `ARC`. The user's SolidWorks parts
  are unaffected: their B-spline relief walls are not analytic, stay as
  they are, and are thickness faces.

**Original plan:**

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

## 7. Manufacturability checks (`flatmat/dfm.py`, implemented)

**As built** (`fm.check(pattern, material=, die_v=, rolling_direction_deg=,
machine_max_force_kn=, rules=DfmRules(...), strict=)`, the CLI `--check`):

- Everything is checked on the flat pattern (2D, shapely) and the bend data.
- **Every finding** is a `Diagnostic` with its location in the flat
  (`data["at"]`). It can be drawn as a marker on layer `DFM` in the DXF and
  goes into the bend table JSON (`dfm`).
- **Severity:** findings are warnings, or errors with `strict`. The
  air-bend radius mismatch is only `info`; `report.ok` ignores infos.
- **Changes against the plan below**, all from running the rules on real
  parts (the user's SolidWorks parts and the BenDFM samples):
  - *Unknown rolling direction* uses the normal minimum radius. The
    along-the-grain minimum applies only when a direction is given; anything
    else flagged every SolidWorks default part.
  - *Minimum radii follow practice:* steel 0.5t, stainless 1t. Xometry's
    0.8t and 2t flagged every CAD default bend.
  - *Bend ends at another bend* (box corners) are not judged as end
    reliefs. They need a corner rule (3D, still to do), and every corner of
    the SolidWorks box sheets was flagged otherwise. Narrow slots near bend
    ends are reliefs, not `FEATURE_TOO_SMALL`.
  - *The relief check fills a band per side:* width `relief_width_t·t`
    (default 0.5t, as in SendCutSend and SolidWorks), reaching `r + t` past
    the tangent line. It reports only when material fills more than half of
    it, so an outline that merely flares out at an angle from the end is not
    flagged.
  - *Narrow slots* need a morphological closing; the opening only finds
    narrow webs. `DFM_HOLE_TOO_SMALL` is for true circles only, with the
    exact diameter.
  - Duplicates are suppressed: a narrow web next to a hole already reported
    near an edge or another hole is not reported again, and there is one
    `CUT_IN_BEND_ZONE` per hole and bend.
- **Tests:** 10 parts with exactly one deliberate problem each (or none);
  the expected values (flange OML, required minimum, force, air radius);
  `strict`, configurable and disabled rules; grain; tonnage; hem material;
  exports.
- **On real parts:** the 4 SolidWorks parts give no problems except the 8
  holes of the bracket, which sit 3.3–3.6 mm from the tangent line where
  2.5t + r = 4.75 mm is required. BenDFM gives 1 U channel and 2 reliefs.

**Original plan:**

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
| `DFM_CORNER_RELIEF` / `DFM_CORNER_GAP` | two bend zones meeting at a corner (§7.6) | relief ⌀ ≥ t clear; folded gap ≥ 0.4 mm |
| `DFM_HEM_TOO_SHORT` | hem return length | ≥ 4·t |
| `DFM_HEM_MATERIAL` | hem on a material or thickness where it is not allowed | material table, Protocase maximum thicknesses |
| `DFM_U_CHANNEL` | U with base < 2 × flange height, or flange taller than the punch clearance or gooseneck throat | — |
| `DFM_TOOL_COLLISION` | collision with punch, die, clamp or beam in every order tried (§7.4) | — |
| `DFM_TOOL_ANGLE` / `DFM_TOOL_LENGTH` | the bend's opening angle below the punch or die angle; bend longer than a tool (§7.5) | 2° overbend allowed |
| `DFM_FORM_NEAR_BEND` / `DFM_FORM_NEAR_EDGE` | rim of an emboss, extruded hole, louvre or bridge (a lance: its opening) to the bend tangent line / to the outline | 4·t + r (CADEX: extruded holes 4t) / 3·t |
| `DFM_BEND_ANGLE` | angle beyond what the shop forms in one bend | 130° (SendCutSend), configurable |

Every finding is a `Diagnostic` with severity `warning` (or `error` if the
caller sets `strict=True`), pointing to the bends, faces and 2D locations
involved, so a UI can draw it.

### 7.4 Collision and sequence (implemented: `flatmat/fold.py`, `flatmat/sequence.py`)

**As built** (`fm.plan_bends(pattern, tooling, order=None)`,
`fm.check(..., sequence=True | [ids])`, CLI `--sequence [ORDER]`):

- **Folding the flat back up** (`FoldModel`):
  - the flat is cut along the BA strips (extended by t past the bend line's
    ends, so the corners separate too); what is left are the flanges;
  - a spanning tree over the bend lines from the largest region fixes which
    side moves;
  - a bend at angle φ turns the far side rigidly. Each region is the sheet's
    mid-surface, so the strip's ends map to the ends of an arc of radius
    r + t/2 (length kept while the bend is part way). This is exact for any
    K and either skin inside: folded STEP parts land t/2 inside the original
    solid to 1e-10 mm (generated parts, a box, the six BenDFM samples);
  - bends whose two sides are one region (no relief) or that have no
    flange between them (a joggle with a zero flat) cannot be placed. They
    are `BendPlan.unchecked`.
- **Press brake:**
  - the bend's inside faces the punch; both legs rise by φ/2;
  - the die is placed so the outside of the legs rests on its shoulders
    (tangent to the shoulder radius);
  - the punch's tip radius is limited to the bend's radius.
- **Tools are prisms along the bend line:**
  - punch and die: segmented to the bend's length unless a length is
    given;
  - the clamp above the punch and the lower beam below the die run the
    machine's full length.
- **Test per state:**
  - every flange's part within a tool's length, seen along the bend line, is
    an affine image of its outline;
  - buffered by t/2 minus a clearance (2 % of t, at least 0.02 mm), it must
    not meet the tool's section;
  - the bend is swept in 5° steps, the full angle first;
  - the part may lie either way round (rotated 180°), which matters for a
    gooseneck: `BendStep.turned`.
- **Order:**
  - the backward search follows Duflou, with weights 4 / 2 / 1 for
    shape-defining (the smaller side's area share), not at a leaf, and
    length;
  - results are memoised per (bend, bends done), and sets of remaining
    bends that dead-end are remembered (no subset is searched twice);
  - a node limit of 2000;
  - turning the part over is `BendStep.flip`;
  - hems (pre-bent, closed with a hemming tool) and lances (formed in a
    punching tool) are made first and not simulated.
- **Output:**
  - a feasible order, or the longest feasible tail plus one `Collision`
    per bend that could not be placed (tool, angle reached, flange
    location, bends done);
  - with `check(sequence=…)`: `DFM_TOOL_COLLISION` findings, the plan in
    the bend table JSON (`dfm.sequence`), each bend's `step` (JSON and
    CSV) and `(n)` before the DXF bend labels.
- **Not modelled:**
  - the bend zones themselves, formed features standing out of the sheet;
  - springback, the back gauge, the operator.
- **Validation:**
  - textbook cases: a deep U (20 wide, 60 legs) is infeasible with a
    straight punch and feasible with a gooseneck, one leg turned into the
    throat;
  - a Z whose long first flange hangs into the lower beam if bent first;
  - a U with a lip where some orders fail;
  - a box with 70 mm walls hits the punch.
  - BenDFM samples with BenDFM's tools: orders found for 0, 1, 3 and 5 in
    < 1 s together. BenDFM's labels are not comparable: they judge three
    sections per bend (a collision only if all three collide), with a die
    block of unpublished size, for one random order.

**Original plan:**

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

### 7.5 Tooling (`flatmat/tooling.py`, implemented)

No open format exists (research §10.8), so tools are a few catalogue numbers,
loaded and saved as JSON (`Tooling.load/save`, CLI `--tooling FILE`):

- **`Punch`**: height (tip to clamp), tip angle, tip radius, body width,
  length (None: segmented to each bend), clamp width.
  - A **gooseneck** keeps its nose up to `throat_bottom`, then leans back at
    45° by `throat_depth`, so a 90° bend's leg clears its underside.
  - The recess on its +y side, up to `throat_top`, is the throat.
- **`Die`**: V (between the sharp edges, before the shoulder radius), angle,
  shoulder radius, block width and height, length (None: segmented), and the
  width of the lower beam below it.
- **Catalogue** (`PUNCHES`, `DIES`; CLI `--punch`, `--die`, `--die V<v>`):
  - punches: straight 88°, gooseneck 88°, acute 30°;
  - dies: V6–V50 at 88°, and a V16 at 30°;
  - BenDFM's tools: a 90° punch 10 × 300 mm, a V40 die.

  These are typical values entered by hand, not any maker's data.
  `Tooling.default(t)` is the straight punch and a V from the Amada rule.
- **In the checks:**
  - V comes from the die unless `die_v` is given;
  - the expected air-bend radius is at least the punch's tip radius;
  - `DFM_TOOL_ANGLE` when the bend's opening angle is below the die's or
    the punch's angle (2° overbend allowed);
  - `DFM_TOOL_LENGTH` when a tool with a length is shorter than the bend;
  - the tools are used by the collision check (§7.4).

### 7.6 Corner rule (implemented)

Where two bends meet (a bend line ends within BA/2 + t of another bend's
line), the end is judged by the corner rule instead of the end-relief rule:

- **`DFM_CORNER_RELIEF`** (2D):
  - the point that tears is the end of the bend's tangent line on the side
    towards the other bend. For a square box corner, that is where both
    tangent lines cross;
  - material may fill at most `corner_fill` (60 %) of a disc of radius
    `corner_relief_t·t` (0.5 t, i.e. a relief ⌀ ≥ t) around it;
  - a bare square notch fills 75 % and is reported. A relief cut past the
    tangent lines fills nothing.
- **`DFM_CORNER_GAP`** / **`DFM_CORNER_OVERLAP`** (3D, on the `FoldModel`):
  - every two folded flanges are plates t thick around their mid-surface;
    their smallest distance comes from outline samples (three layers, at
    most 0.25 mm apart);
  - closer than `corner_gap` (0.4 mm, SendCutSend's 0.015") is
    `DFM_CORNER_GAP`, intersecting is `DFM_CORNER_OVERLAP`;
  - flanges joined by one bend, or by two parallel bends (a joggle, a U),
    are skipped.
- **Tests:** a box built from two U profiles:
  - with square reliefs in the base's corners: clean;
  - short bends ending 0.2 mm before the long ones: 4 × `DFM_CORNER_RELIEF`;
  - ending 2 mm before, without a relief: the end-relief rule, as before;
  - closed corners (tabs on all walls) 0.1 mm apart: 4 × `DFM_CORNER_GAP`
    at exactly 0.1 mm. At 1 mm: none.
- **Caveat:** these corner ends were not judged before. On the BenDFM
  samples the corner rule now reports their bare square corner notches
  (fill 75 %); the private SolidWorks boxes are still to be re-checked.

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
- B-spline recognition is tested only on exact (rational) NURBS from
  `BRepBuilderAPI_NurbsConvert`. Exporters that *approximate* cylinders with
  polynomial splines (deviations ~1e-6 to 1e-3) are untested; the
  tolerance is `Tolerances.linear` (1e-4 mm).
- ~~Fast healing on 7.8~~: resolved via `$CSF_STEPDefaults` (§6.1).
- Is OCCT 8.0's STEP reader thread-safe? (Irrelevant while we use one
  process per part.)
- Hem allowance and extruded-hole pre-hole tables: no open data found. The
  pre-hole uses the neutral fibre (§5.4), which matches the one published
  formula; `unfold(pre_hole=…)` takes a shop's table.
- Real punch and die dimensions: the built-in catalogue is typical values,
  not checked against a maker's catalogue.
- Real STEP exports: one SolidWorks 2023 part tried (missing cylinder
  pcurves, now handled). Inventor, Fusion, Onshape, Creo and Solid Edge are
  still to collect.
- `ShapeUpgrade_UnifySameDomain` with `KeepShapes` and canonical recognition
  are untested.
