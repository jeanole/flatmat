# Phase-0 spike

A throwaway prototype that checks the design
([`docs/design.md`](../docs/design.md) §1) end to end on real OCP before the
library is written. It covers planes and cylindrical bends only, in one file
with no error taxonomy. **It is not the library.**

| File | Purpose |
|---|---|
| `api_probe.py` | checks that every OCP name the design uses exists in the installed OCP |
| `occ_compat.py` | the differences between OCP 7.8 and 8.0 that the probe and the runs found |
| `gen.py` | generates test parts (2D profile + bends, extruded, optional holes) with exactly known flat lengths, written to STEP |
| `unfold.py` | the pipeline: STEP → adjacency graph → skins → thickness → bends → per-face development → placement by shared vertices → outline by vertex chaining → bend lines → self-checks |
| `run.py` | generated parts, BenDFM samples, scaling series, DXF output, timings |
| `results/` | output of the final run on OCP 8.0.1.1 (Python 3.13) and 7.8.1.1 (Python 3.12) |
| `*.png` | DXF previews (BenDFM part 0, the hole through a bend, the hat profile) |

## Running it

```sh
uv venv -p 3.13 .venv && uv pip install -p .venv "cadquery-ocp-novtk==8.0.1.1.0" ezdxf
.venv/bin/python spike/api_probe.py
# BenDFM sample files (MIT repo; the files are not committed here):
mkdir -p /tmp/bendfm && cd /tmp/bendfm && for i in 0 1 2 3 4 5; do for f in $i.stp ${i}_unfolded.stp \
  ${i}_sequence.json ${i}_labels.json; do curl -sSfO https://raw.githubusercontent.com/UGent-CVAMO/bendfm/main/example_files/$f; done; done
.venv/bin/python spike/run.py --bendfm /tmp/bendfm --scaling
```

## Results

### Correctness: everything passes on both OCP versions

**Generated parts** (K = 0.44, each one round-tripped through STEP):

- Covered: L (up and down), U, Z, hat (4 bends), obtuse 135°, acute 45°, open
  hem 180°, joggle (two tangent cylinders of opposite convexity), holes in
  flanges, and a hole drilled radially through a bend.
- Flat area against the analytic `width × (Σ flats + Σ θ·(r + K·t))`: relative
  error **≤ 3.5e-14**.
- The two parts with holes show ~2e-5, but only because the spike measures
  area on circles sampled with 64 points; the library writes exact arcs.
- Thickness, inner radius, angle and the up/down pattern are recovered for
  all 11 parts.
- The B-spline pcurve of the hole through the bend was developed by mapping
  its poles (affine map, weights kept). It matches pointwise mapping to
  **2.7e-15**, so developing B-splines exactly works as designed.

**BenDFM samples 0–5** (K = 0.45, their convention, see `docs/research.md`
§5):

- All 6 recognised.
- Thickness equals `part_thickness` from their JSON.
- Bend count equals the number of entries in their sequence JSON (2–8 bends).
- Our flat area against the top-skin area of *their* unfolded STEP: relative
  error **1.5e-9 to 2.9e-5**, again limited by the sampling of their rounded
  flanges.

**OCP 7.8.1.1 and 8.0.1.1 give identical numbers.**

### Speed (OCP 8.0, Python 3.13, this container)

| Part | Faces | Total incl. STEP read |
|---|---|---|
| generated L / U / Z / hat | 10–22 | 9–25 ms |
| BenDFM samples | 30–86 | 22–87 ms |
| L with 100 holes | 110 | 133 ms |
| L with 1000 holes | 1010 | 1.47 s, of which STEP read 0.90 s |

7.8 is about 1.0–1.5× slower. The scaling series is now linear; the first
version needed **19 s** for 1000 holes. The fixes are findings for the design:

| Hotspot (1000 holes) | Before | After | Cause → fix |
|---|---|---|---|
| STEP import | 4.1 s (8.0), 4.9 s (7.8) | 0.18 s / 0.25 s transfer | The reader's built-in healing runs `ShapeFix_Shell`, which is **quadratic in the face count**. Switch off `FixShellMode` and `FixSolidMode` and check the orientation ourselves (volume > 0). 8.0: `SetShapeFixParameters` **after** `ReadFile`, as a full key→value map. 7.8: a resource file named `STEP` in `$CSF_STEPDefaults`. |
| `BRepCheck_Analyzer` | 1.3 s (8.0), **13 s** (7.8) | not run | Quadratic on a face with many wires, even topology-only. Our own self-checks (area identity, closed outline, cycle consistency) cover what matters. BRepCheck is for diagnostics only. |
| edge → face table | 0.9 s | 4 ms | Iterating an OCP `TopTools_ListOfShape` from Python costs ~0.3 ms per list. Use `First()`/`Last()`; an edge has at most 2 faces. |
| edge classification | 9 s | 55 ms | `BRepAdaptor_Surface` was rebuilt per edge. Cache one per face. |
| point-in-face | 9 ms per call | 0.8 ms | `BRepClass_FaceClassifier` rebuilds per call. Build one `IntTools_FClass2d` per face (21 ms) and reuse it. |
| thickness | 1.5 s | 70 ms | (a) A ray cast against the whole solid → analytic antiparallel plane pairs. (b) Comparing all cylinder pairs → bucket by axis line. |

What is left is linear. STEP parsing (`ReadFile`, ~0.4 s for 1000 holes) is
the largest single item, and it is OCCT C++ in any language.

### OCP 7.8 vs 8.0 API differences found

| | 7.8 | 8.0 |
|---|---|---|
| shape maps and lists | `TopTools_IndexedMapOfShape`, `TopTools_IndexedDataMapOfShapeListOfShape`, `TopTools_ListOfShape` | `OCP.collections.IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher`, `IndexedDataMap_TopoDS_Shape_List_TopoDS_Shape_TopTools_ShapeMapHasher`, `List_TopoDS_Shape` |
| downcasts | `TopoDS.Face_s(...)` | `TopoDS.Face(...)` |
| face classifier | `BRepTopAdaptor_FClass2d` exists | gone; use `IntTools_FClass2d` (both) |
| healing control | resource file via `$CSF_STEPDefaults` | `SetShapeFixParameters(map)` after `ReadFile` |
| `SetShapeProcessFlags` | missing | present, but **not callable from Python** (`std::bitset` is not bound) |

Present in both: `ShapeAnalysis_CanonicalRecognition`,
`GeomConvert_SurfToAnaSurf`, `GeomLib_IsPlanarSurface`,
`IntCurvesFace_ShapeIntersector`, `ShapeUpgrade_UnifySameDomain` (with
`KeepShapes` and `History`), `BRepClass3d_SolidClassifier`,
`BRepExtrema_DistShapeShape`.

General OCP pitfall: C++ `double&` out-parameters do not come back to Python
(for example the `First`/`Last` of `BRep_Tool.CurveOnSurface`). Use
`BRepAdaptor_Curve2d`, which exposes the parameter range.

### First real CAD export: SolidWorks 2023 (AP214), user-provided

The file is not committed. It is a U-shaped holder: one solid, 56 faces (18
planes, 38 cylinders), 2 mm thick, 15 holes.

- **SolidWorks writes no pcurves for cylindrical faces**: all 152 edge uses on
  cylinders lacked them. Planes are fine, because OCCT computes their pcurves
  on the fly.
- The reader's full healing adds the missing pcurves, but the fast mode
  above skipped them: **`FixShellMode = 0` also skips every face- and
  wire-level fix below the shell.** So the fast path now calls
  `ShapeFix_Wire.FixEdgeCurves()` itself on every wire that has an edge
  without a pcurve. This runs in place and took 4 ms (8.0) / 6 ms (7.8).
- Without that fix, `BRepAdaptor_Curve2d` on a missing pcurve **raises
  `Standard_NullObject` on 7.8 and segfaults the interpreter on 8.0**. The
  spike now checks `BRep_Tool.CurveOnSurface` first and refuses with
  `MISSING_PCURVE` instead of crashing.
- After the fix, OCP 7.8 and 8.0 give identical results:
  - thickness 2.000 (confidence 1.0)
  - two 90° bends, r_in 1.0, both bending the same way
  - flat 250.000 × 241.906 mm at K = 0.44
  - area identity 1.1e-5 (sampled circles)
  - read 27 ms + unfold 79 ms on 8.0
- The regression set (11 generated parts, 6 BenDFM samples) still passes on
  both versions.

## Not covered by the spike

- Exact `Line`/`Arc`/`Spline` output (the spike writes sampled polylines into
  the DXF).
- `UnifySameDomain`, canonical recognition, sharp bends, cones and the
  diagnostics taxonomy.
- Self-overlap check, assemblies and DFM.
- The full BenDFM set (only samples 0–5) and real-world STEP files from
  commercial CAD systems.
