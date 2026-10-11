# flatmat — research findings

Collected 2026-10-11 while scoping flatmat. Desk research only: nothing
below was run or benchmarked unless it says so. Each section ends with what
could not be verified. The design that follows from it is in
[`design.md`](design.md).

## Contents

1. [Landscape: what already exists](#1-landscape-what-already-exists)
2. [How the open unfolders work](#2-how-the-open-unfolders-work)
3. [Language and performance: Python, C++ or Rust](#3-language-and-performance-python-c-or-rust)
4. [Kernel, bindings, licensing](#4-kernel-bindings-licensing)
5. [Test data](#5-test-data)
6. [Bend allowance and K-factor](#6-bend-allowance-and-k-factor)
7. [Flat-pattern DXF and press-brake conventions](#7-flat-pattern-dxf-and-press-brake-conventions)
8. [Hard geometry](#8-hard-geometry)
9. [Robustness on real-world STEP files](#9-robustness-on-real-world-step-files)
10. [Manufacturability rules for laser + press brake](#10-manufacturability-rules-for-laser--press-brake)

---

## 1. Landscape: what already exists

**Finding: there is no standalone, permissively licensed, pip-installable
library that takes an arbitrary STEP solid and returns thickness, a bend
table (radius, angle, up/down, axis, K) and a flat-pattern DXF with bend
lines.** That gap is what flatmat is for.

### Open source

| Project | Licence | State | Notes |
|---|---|---|---|
| [build123d](https://github.com/gumyr/build123d) sheet metal ([PR #1434](https://github.com/gumyr/build123d/pull/1434)) | Apache-2.0 | merged into `dev` 2026-09-29, **not** in release 0.13.0 (2026-09-21) | `sheet_shells(solid)` works on imported STEP solids; `unfold(shell, params)`. Planes and cylinders only; no bend lines or bend table; one solid; requires `cadquery-ocp-novtk >= 8.0, < 8.1` and pulls scipy, sympy, ipython, requests, lib3mf. |
| [FreeCAD SheetMetal](https://github.com/shaise/FreeCAD_SheetMetal) (new unfolder, ex [YouCanNotUnfold](https://github.com/alexneufeld/YouCanNotUnfold)) | LGPL-2.1 (README still says GPLv3; relicensed in v0.4.00) | V0.8.24, 2026-09-30 | The most capable open unfolder: bend info, DXF/SVG with bend lines and angle labels, K-factor spreadsheets. Needs FreeCAD (`Part`, `TechDraw`) and networkx; no CLI; headless use via `freecadcmd` looks possible but is untested. |
| [tmkaio/progressive-die-designer](https://github.com/tmkaio/progressive-die-designer) | MIT | 3 commits, 0 stars | Proof of concept: STEP via CadQuery/OCP, bend cylinders, θ·(r + K·t), DXF via ezdxf. |
| [ulrich1a/sheet_ufo](https://github.com/ulrich1a/sheet_ufo) | — | historic | The original FreeCAD unfolder macro. |
| `cadunpack` (PyPI) | unknown | **withdrawn** | Described as "STEP assembly extraction + sheet-metal unfolding with raw DXF output"; the PyPI JSON and GitHub repo now return 404. |
| [quiddity](https://github.com/pzfreo/quiddity) | Apache-2.0 | — | STEP feature recognition for machined parts; no unfolding. |

Not usable: KoksaSa/sheet-metal-unfolder (CC BY-NC), browser and TypeScript
calculators, quirxi/SheetMetalCone (Inkscape cone extension). Unfolders
inside applications (SmartTwinModeller, on-demand-craft-greece, forge,
emithran) were found by code search but not read.

**Rust**: no sheet metal unfolder found. [Fornjot](https://github.com/hannobraun/fornjot)
was archived on 2026-06-19.

### Commercial

| Product | What it offers | Price |
|---|---|---|
| [CADEX MTK](https://docs.cadexsoft.com/mtk/mtk_sheet_metal_unfolding) | unfolds shells/solids with simple, curved and hem bends; `SetKFactor` (default 0.5); C++/C#/Python; 2025.8 added complex formed features | separate licence, not published |
| [OCCT Sheet Metal Operations](https://www.opencascade.com/components/sheet-metal-operations/) / [Unfolding](https://www.opencascade.com/components/unfolding-component/) | recognises bends, hems, jogs and bridges; flat pattern with K-factor; bend data for costing; C++ | quote only |
| HOOPS Exchange | reads native flat-pattern data from CAD files; no unfolder found | — |
| Lantek Flex3d | CAD add-ins on native formats; STEP not confirmed | — |
| Bendline (hosted, built on the OCC/FreeCAD stack) | STEP → DXF with cut/bend layers, PDF bend table, press-brake CSV; API | ~£0.90 per download |
| [SendCutSend](https://sendcutsend.com/faq/how-to-configure-bends-for-my-step-file) | unfolds uploads internally; needs one solid built to their radius and K | — |

Not verified: what cadunpack did and why it vanished; unfolder.nl's offer;
Spatial CGM; Xometry; CADEX/OCCT prices; FreeCAD's unfolder under
`freecadcmd`; build123d on real STEP parts with hems or reliefs.

---

## 2. How the open unfolders work

Read from source (`raw.githubusercontent.com`), not executed.

### build123d (dev branch)

| File | LOC | Content |
|---|---|---|
| `src/build123d/sheet_utils.py` | 568 | `SheetMetalParameters`, `neutral_radius`, `bend_allowance`, `is_positive_bend`, unfold core `_unfold_shell` (~370 LOC) |
| `src/build123d/operations_sheet.py` | 3135 | `unfold()` (l.1619), `sheet_shells` (l.1674), `_bend_pairs` (l.1775), `_flat_pairs` (l.1802), plus flange/hem/jog/miter/reliefs |
| `src/build123d/build_sheet.py` | 359 | `BuildSheet` |
| `src/build123d/topology/two_d.py` | — | `Shell.unfold` (l.3221) |

- **Thickness / skins** (`sheet_shells`):
  - Bend pairs are coaxial cylinders with overlapping axial spans and
    different radii.
  - Flat pairs come from extruding a "shadow" solid along −normal and taking
    the first antiparallel face it hits. That is one boolean per candidate,
    O(n²).
  - Thickness is the median of the largest distance cluster: weighted by
    count when there are bends, otherwise by area.
  - Faces are sewn into exactly two shells; it does not say which one is
    inside.
- **Unfold** (`_unfold_shell`):
  - Only planes and cylinders, at least one plane, no non-manifold edges.
  - The face adjacency is a dict built from
    `TopExp.MapShapesAndAncestors(EDGE→FACE)`.
  - Each face is developed on its own from its pcurves
    (`BRep_Tool.CurveOnSurface` → `BRepBuilderAPI_MakeEdge` on the XY
    plane), closed with `ShapeFix_Wire`/`ShapeFix_Face`.
  - Cylinder u is scaled by the neutral radius with `BRepBuilderAPI_GTransform`
    (diag(R, 1, 1)).
  - BFS from the largest plane; each child is rigidly snapped to its placed
    parent along the shared edge, with the side taken from
    `_topods_material_side`.
- **Cycles**: a non-tree link must reproduce the same developed endpoints,
  otherwise it raises *"The sheet contains a cycle that requires a cut"*.
- **K**: one scalar for the part (default 0.5); `BA = rad(θ)·(r_in + k·t)`;
  neutral offset depends on the reference surface (inside, outside, mid,
  neutral).
- **Holes across bends**: exact via pcurves, provided they were cut by a
  solid cutter.
- **Bend lines**: not emitted; which faces came from bends is lost.
- **Limits**: planes and cylinders only, min bend radius 1e-3, one K, no
  overlap check.

### FreeCAD SheetMetal new unfolder (`SheetMetalNewUnfolder.py`, 1610 LOC, © 2025 Alex Neufeld)

- **Thickness**, trying in order:
  1. the mode of straight edge lengths normal to the root face;
  2. the parallel face closest in area;
  3. the mode of coaxial cylinder radius differences.
- **Graph**:
  - networkx; a link is added only between two faces that share an edge
    and are tangent.
  - The component containing the user-picked root face is the skin.
  - Seam edges (one ancestor) drop out; seam-split cylinders are joined by
    the "coaxial, same radius" rule.
- **Traversal and placement**:
  - `nx.minimum_spanning_tree` with an effectively arbitrary weight; cycles
    are dropped silently.
  - Each cylinder gets a 4×4 unbend matrix, `A · translate(allowance) ·
    rotate · A⁻¹`. The local frame comes from four UV-corner cases
    ("figured out by brute force").
  - Transforms are chained along the shortest path from the root.
- **Cylinder unroll**: pcurve `(u, v) → (v − vmin, (u − umin) · BA/θ)`.
  Lines and B-splines are supported; other curve types turn into a
  warning and the face is skipped.
- **Suspected bugs** (from reading the code, not run):
  - The B-spline rebuild calls `buildFromPolesMultsKnots(poles, weights)`
    without degree, mults or knots, which distorts non-default curves.
  - In the K-table lookup the `while r_over_t <= values[i]` loop never runs,
    so it always interpolates between rows 0 and 1.
- **Cones**: they are recognised as tangent, so they enter the graph, but
  they are carried along rigidly instead of being unrolled, which gives
  wrong output.
- **Output**: bend lines at the midline of each strip, with angle and
  radius labels. The 2D cleanup converts B-splines to lines, arcs and
  biarcs at 0.1 mm. DXF and SVG go through FreeCAD Sketcher objects.
- **Porting to raw OCCT**:
  - The algorithm maps 1:1 onto OCP plus a graph:
    - `TopExp`, `BRepAdaptor_Surface`, `BRepTools::UVBounds`,
      `BRep_Tool::CurveOnSurface`, `BRepLProp_SLProps`
    - `gp_Trsf`/`gp_Ax3`
    - `ShapeAnalysis_FreeBounds` + `ShapeFix_Wire`, `BRepAlgoAPI_Common`
  - Not in OCCT: `FaceMakerBullseye` (needs our own wire nesting) and
    `toBiArcs`.

### Quaoar / Analysis Situs blog series

[Part 1](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-1),
[2](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-2),
[3](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-3),
[4](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-4-1),
[5](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-5),
[6](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-6),
[12](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-12-searching-for-feasible-bends),
plus [the CADA paper post](https://quaoar.su/blog/page/our-paper-in-computer-aided-design-and-applications)
and [UV mapping for sheet metals](https://quaoar.su/blog/page/uv-mapping-for-sheet-metals).
Analysis Situs is BSD-3 C++; the sheet metal module (SMRU) source was not
found publicly.

- **Recognition: attributed adjacency graph (AAG).**
  - Faces are nodes; each adjacency is labelled convex, concave or smooth.
  - Seed from the largest face; thickness by ray cast to the opposite face.
  - Grow across smooth edges until sharp borders are reached.
  - Removing the thickness faces (a ring of degree-2 nodes) leaves exactly
    two skin components.
- **Unfolding (part 5, the revised approach):**
  - Each face is unrolled individually.
  - Neighbours attach at "anchor" shared edges via trihedrons.
  - No special status for bends, the result does not depend on traversal
    order, and gaps are never filled. This is the design build123d also
    uses.
- **Cylinder unroll (part 2):** u is scaled by R; pcurves are
  re-approximated with `AdvApprox_ApproxAFunction` +
  `GeomLib_MakeCurvefromApprox`. Circles then come out as B-splines and
  need recognising again.
- **K (part 3):** `R_N = r + k·T`, with k depending on r/t, angle and
  forming method → a table.
- **Pitfalls they list:**
  - rolled parts without flanges defeat rotation-based unfolding;
  - bend compensation computed twice makes results fragile;
  - angle round-off leaves the flat slightly non-planar;
  - bend faces glued to other bend faces are ambiguous;
  - tubes fail the path-independence check;
  - self-overlapping flats have to be detected;
  - mismatched radii and variable thickness appear in real parts (prefer
    warn-and-skip over abort);
  - two halves misaligned by 0.24 mm break recognition;
  - periodic UV goes out of range on cylinders;
  - flipped orientation flips up/down.

### Comparison

| Aspect | build123d | FreeCAD new unfolder | Quaoar |
|---|---|---|---|
| Thickness | coaxial ΔR, else area-weighted shadow pairs (O(n²) booleans) | normal edges → parallel face → coaxial ΔR | ray cast |
| Skin split | sew faces at distance t → 2 shells | tangent component of the root | AAG minus thickness ring |
| Placement | per-face develop + rigid snap on shared edge | chained 4×4 unbend matrices | per-face unroll + anchors |
| Cycles | checked, raises | dropped silently | path-independence check |
| K | scalar | scalar or r/t table (buggy interpolation) | table |
| Cones | rejected | wrong output | not covered |
| Bend lines | none | midline + labels | stored on links |
| Core LOC | ~550 | ~1100 | closed |

### Lessons adopted in the design

1. Split skins by tangency, then confirm the two components are offsets at
   distance t.
2. Develop every face independently from its pcurves and snap neighbours
   on shared edges — order-independent, handles seam-split and multi-face
   bends without case analysis.
3. Check path independence on non-tree links.
4. Keep provenance: every flat edge maps back to its source edge and bend.
5. Compute K and the neutral radius from the inside radius with explicit
   sign handling; use an r/t table with correct interpolation; record K per
   bend.
6. Keep bend development exact: preserve degree, knots and mults (or use
   GTransform). Fit arcs only at export, with an explicit tolerance.
7. Refuse or warn explicitly on unsupported geometry.
8. Post-check the flat: self-overlap, `BRepCheck`, area identity.
9. Thickness by ray casts (`IntCurvesFace_ShapeIntersector`), not O(n²)
   booleans.
10. Depend only on OCP; write DXF with ezdxf.

---

## 3. Language and performance: Python, C++ or Rust

**Conclusion: Python on OCP.** STEP import and healing dominate the run
time and are C++ in every option. Unfolding is O(faces) and costs
milliseconds even with Python call overhead.

### Where the time goes

- On perforated sheet metal, **shape healing takes ~90 % of STEP import
  time**. Disabling the self-intersection checks made some of ~2,600 test
  cases several times faster with no effect on recognition
  ([Quaoar](https://quaoar.su/blog/page/perforated-sheet-metals)). The same
  post found the STEP writer slowed by `BRepAdaptor_Surface` computing UV
  bounds per edge.
- **Large-file import reports:**
  - 34 s TransferRoot on a 4.5 MB file (2011)
  - a 400 MB file at ~57 s read + 44 s transfer in a release build
  - 1 GB at ~6 min (GCC) vs ~19 min (MSVC)
    ([tracker 32014](https://tracker.dev.opencascade.org/view.php?id=32014))
  - an 8 MB file taking > 10 min and > 10 GB on 7.9.0, said to be fixed
    after 8.0.0rc3
    ([forum](https://occt3d.com/dev/content/slow-import-step-files-and-high-memory-usage-occt-790/index.html))
- **OCCT releases:**
  - 7.8 claims up to 40 % faster large reads
  - 8.0.0 has a faster BSpline cache and TopExp, and bit-packed TopoDS
    ([release](https://dev.opencascade.org/content/open-cascade-technology-800-release))
  - multithreaded STEP transfer was postponed in Aug 2025 as "not stable"
- **No published timings** were found for topology exploration,
  BRepAdaptor queries, booleans or BRepMesh on typical small parts.
  Estimate: for 50–500 faces, read + heal dominate and unfold is
  milliseconds.

### Binding overhead

- pybind11 ≈ 200 ns per call vs ≈ 30 ns for a plain C call
  ([measurement](https://gitlab.mpcdf.mpg.de/mtr/pypocketfft/-/merge_requests/32));
  older numbers put it ~4× slower than Cython
  ([#1227](https://github.com/pybind/pybind11/issues/1227)).
- nanobind claims ~3× lower overhead than pybind11, ~10× when classes are
  passed ([benchmark](https://nanobind.readthedocs.io/en/latest/benchmark.html)).
- No OCP-specific benchmark found. No evidence that OCP releases the GIL; in
  [OCP #220](https://github.com/CadQuery/OCP/issues/220) the free-threading
  contributor recommends single-thread access, though 4 threads gave a
  348 % speedup on independent work.

### Existing unfolders' speed

- FreeCAD's new unfolder is pure Python and v0.8.24 lists an "unfolder
  speedup" without numbers. An Ondsel post calls unfolding "rather slow",
  worth a C++ rewrite ([post](https://ondsel.com/blog/sheetmetal-strategy)).
- build123d reports only that its new thickening is > 10× faster.
- Quaoar part 12 says "tens of seconds" for its bend search is not
  acceptable.
- CADEX: no timings published.

### Rust

| Option | State |
|---|---|
| [opencascade-rs](https://github.com/bschwind/opencascade-rs) | "a major work in progress", one spare-time maintainer, 52 open issues. `Shape` has `read_step`; low-level coverage only where someone hand-bound it. BRepAdaptor, pcurves etc. would need our own C++ glue. LGPL-2.1; default build links OCCT statically. |
| [truck](https://github.com/ricosjp/truck) / truck-stepio 0.3.0 (2024) | Open issues: an edge silently dropped (#133), parabola/hyperbola handling (#132), panics on NaN or degenerate placements (#128, #123). Fork monstertruck 0.4.1. Apache-2.0. |
| [foxtrot](https://github.com/formlabs/foxtrot) | 2021 hackathon project, "nowhere near industrial strength". |
| [ruststep](https://github.com/ricosjp/ruststep) | "experimental stage, should not be used for products". |
| step-io 0.2.5 | experimental. |
| Fornjot | archived 2026-06-19. |

No evidence that any pure-Rust stack reads real AP203/AP214 with B-spline
faces robustly enough for production.

### C++ route

- nanobind: 2.7–4.4× faster compiles, 3–5× smaller binaries, stable ABI
  (one wheel across Python versions). pybind11 3.0 has free-threading
  support.
- **The catch**: a C++ extension can only exchange `TopoDS_Shape` with OCP
  when compiled against the *identical* OCCT build. OCP wheels ship the
  shared libraries but not the headers. In practice the C++ route means
  building OCCT yourself on every platform (OCP's Windows build alone is
  ~5.5 h on GitHub runners) and giving up interop with build123d and
  CadQuery.

### Concurrency

- Separate objects on separate threads is the usual advice; the same object
  across threads is not safe.
- STEP reader thread safety:
  - Crashes in the `STEPControl_Reader` constructor with 16 threads on
    7.8.1 ([forum](https://occt3d.com/dev/content/stepcontrolreader-access-violation-reading-location/index.html)).
  - A mutex added to `XSControl_WorkSession` in Jan 2025, yet random crashes
    in `STEPCAFControl_Reader::ReadStream` were still reported on 7.9.0
    ([forum](https://occt3d.com/dev/content/read-access-violation-stepcafcontrolreaderreadstream-concurrent-step-file-reading-occt-790/index.html)).
  - 8.0 replaced `Standard_Mutex` with `std::mutex`, made error handlers
    thread-local, and 8.0.0p1 fixed `Standard_Strtod` sharing a pool across
    threads
    ([notes](https://newreleases.io/project/github/Open-Cascade-SAS/OCCT/release/V8_0_0_p1)).
- **Practice: one process per worker.**

Not verified: per-call costs of OCCT operations on small parts; any OCP
overhead or GIL measurement; seconds per part for any unfolder; whether
OCCT 8.0's STEP reader is fully thread-safe.

---

## 4. Kernel, bindings, licensing

| Option | Facts | Licence |
|---|---|---|
| [cadquery-ocp](https://pypi.org/project/cadquery-ocp/) | 8.0.1.1.0 (2026-10-05), OCCT 8.0.1, CPython 3.11–3.14, Windows x86-64, manylinux_2_28 x86-64/aarch64, macOS x86-64/arm64; ~48–68 MB wheels; built with VTK 9.6.2 | bindings Apache-2.0; bundled OCCT LGPL-2.1 + exception |
| [cadquery-ocp-novtk](https://pypi.org/project/cadquery-ocp-novtk/) | same version and platforms, no VTK → best fit for a headless library | same |
| [pythonocc-core](https://github.com/tpaviot/pythonocc-core) | 8.0.1.1 (2026-10-05), SWIG, conda-only in practice (PyPI has a 2014 stub) | LGPL-3.0 |
| [opencascade-rs](https://github.com/bschwind/opencascade-rs) | see §3 | LGPL-2.1 |
| truck | see §3 | Apache-2.0 |

Wheels are built by [ocp-build-system](https://github.com/CadQuery/ocp-build-system)
(VTK, novtk and proxy variants). There are no free-threaded (3.13t/3.14t)
OCP wheels; [#220](https://github.com/CadQuery/OCP/issues/220) is open.

**OCCT's LGPL-2.1 with
[exception](https://github.com/Open-Cascade-SAS/OCCT/blob/master/OCCT_LGPL_EXCEPTION.txt):**

- A Python package that imports OCP is a "work that uses the library", so
  flatmat keeps its own licence (Apache-2.0). Add an OCCT notice and point
  to the OCCT source.
- Proprietary users (such as Snest) are fine with dynamic linking: they
  must include the LGPL text and notice, provide or offer the OCCT source,
  allow the libraries to be replaced, and not forbid reverse engineering
  for debugging those changes.
- Static linking (for example opencascade-rs `builtin`) would also require
  shipping relinkable object files.

This is a reading of the licence text, not legal advice.

Not verified: the pybind11/nanobind licences (BSD-3 from memory); the
pythonocc-core LICENSE file itself; truck-stepio's robustness.

---

## 5. Test data

| Dataset | Sheet metal? | Format | Licence |
|---|---|---|---|
| **BenDFM** ([paper](https://arxiv.org/html/2603.13102v1), [repo](https://github.com/UGent-CVAMO/bendfm), [Zenodo](https://zenodo.org/records/18622958)) | **yes** — 20 000 synthetic bent parts (BenDFM 14k + BenDFM-U 6k), train/val/test | per part: folded STEP, unfolded STEP, bend-sequence JSON, labels JSON; one 3.1 GB archive | data **GPL-3.0-or-later**, code MIT |
| [ABC](https://deep-geometry.github.io/abc-dataset/) | no sheet metal subset found | STEP, Parasolid | model copyright stays with the Onshape creators |
| [Fusion 360 Gallery](https://github.com/AutodeskAILab/Fusion360GalleryDataset/blob/master/LICENSE.md) | no | STEP | non-commercial research only, no redistribution |
| [MFCAD](https://github.com/hducg/MFCAD) | no (machining features) | STEP + labels | MIT |
| [MFCAD++](https://pure.qub.ac.uk/en/datasets/mfcad-dataset-dataset-for-paper-hierarchical-cadnet-learning-from/) | no | 59 665 STEP | CC BY |
| MFInstSeg ([AAGNet](https://github.com/whjdark/AAGNet)) | no | ~62k STEP + AAG | code MIT, dataset licence not stated |
| SMCAD (Ma & Yang 2024) | yes (hinges, stiffeners) | unknown | unknown, availability not found |
| FreeCAD SheetMetal `SMTests` | scripts only seen (`testFolder.py`, `testKfactor.py`) | — | LGPL-2.1 |

**Consequence:** GPL-3.0 test data can be downloaded and used for testing
but must not be committed to this Apache-2.0 repository. Fixtures in the
repo are generated by our own code (see the design, §3).

Not verified: the MFInstSeg licence; whether SMCAD is available; whether
`SMTests` contains models (the listing was truncated); which K-factor and
units BenDFM uses.

---

## 6. Bend allowance and K-factor

- **Common ("ANSI") convention:** K = distance of the neutral fibre from
  the inside surface divided by t, 0 ≤ K ≤ 0.5. This is a CAD-vendor
  convention; no actual ANSI standard number was found behind it.
  - `BA = π·(R + K·T)·A / 180`
    ([SolidWorks](https://help.solidworks.com/2022/English/SolidWorks/sldworks/c_K_Factor.htm))
  - `OSSB = tan(A/2)·(R + T)`, `BD = 2·OSSB − BA` (textbook, from memory)
- **DIN 6935:**
  - correction factor `k = 0.65 + 0.5·lg(r/s)`, range 0–2
    ([PTC Creo](https://support.ptc.com/help/creo/ced_modeling/r20.7.0.0/en/ced_modeling/OSDM_Modules/SheetMetal_adm_allow.html))
  - DIN measures in units of s/2, so **K_ANSI = k_DIN / 2**; CATIA's DIN
    default `K = (0.65 + log(R/T)/2)/2` confirms it
  - From memory, standard paywalled
    ([Beiblatt 1](https://din.de/en/wdc-beuth:din21:144420042)):
    - k = 1 for r/s ≥ 5.
    - Compensation `v = π(180−β)/180·(r + k·s/2) − 2(r + s)` for β 0–90°.
    - For β 90–165° the last term is multiplied by `tan((180−β)/2)`.
    - v ≈ 0 for β 165–180°.
- **Typical K (ANSI convention, air bending):**

  | Source | Values |
  |---|---|
  | [Durmapress](https://www.durmapress.com/k-factor-in-press-brake-bending-formula-chart-practical-guide/) soft / medium / hard | r ≤ t: 0.33 / 0.38 / 0.40 · 1–3t: 0.40 / 0.43 / 0.45 · > 3t: 0.50 · bottoming 0.42–0.46 · coining 0.38–0.44 |
  | [Xometry](https://xometry.pro/en/?p=1993) | air, R ≤ T: 0.36–0.45 (stainless high) · bottoming 0.44–0.48 · R > 3T: 0.5 |
  | FreeCAD README example | r/t 1 → 0.38, 3 → 0.43, 99 → 0.5 |

- **Air bending:** the inside radius comes from the die opening, not the
  punch ([The Fabricator](https://www.thefabricator.com/article/bending/how-to-predict-an-air-formed-inside-bend-radius-with-precision)):
  - mild cold-rolled steel 15–17 % of V
  - 304 stainless 20–22 %
  - hot-rolled pickled and oiled 14–16 %
  - 5xxx aluminium 10–12 %

  So the radius should come from the model (or the die), not be assumed.
- **Open tables:** none formally licensed. FreeCAD's spreadsheet format
  (Radius/Thickness | K-factor ANSI/DIN) is the closest reusable schema.

---

## 7. Flat-pattern DXF and press-brake conventions

- **No standard layer naming exists.** CAD exporters let the user map
  bend-up and bend-down lines to layers; SolidWorks 2022 dropped its default
  `BEND` layer
  ([SolidWorks](https://help.solidworks.com/2012/english/solidworks/sldworks/c_mapping_bend_directions_when_exporting_to_dxf_dwg_files.htm)).
- **SendCutSend**
  ([blog](https://sendcutsend.com/blog/how-to-map-bend-lines-into-layers-and-linetypes-in-solidworks/)):
  - one `BEND` layer for up and down
  - a non-hidden linetype (hidden bend lines are not detected)
  - 1:1 in mm or inches
  - bend lines continuous across cut-outs (SolidWorks splits them)
- **Amada Dr.Abe**: users report it needs DXF 2010 or older
  ([forum](https://forums.autodesk.com/t5/inventor-forum/remove-bend-lines-in-flat-pattern-for-dxf/td-p/11636991)).
- **Delem DA-66T**
  ([manual](https://machinemfg.com/wp-content/uploads/2017/09/Press-Brake-Delem-DA-66T-Operation-Manual.pdf)):
  - Layers are assigned to "bendline" or "contour", or bend lines are
    detected automatically.
  - **The bend angle is a text label next to the line: positive = up,
    negative = down, `H±angle` for a hem, `R value` for the radius.**
  - Thickness and material can be supplied as text.
  - Handles both cutting and projection dimensions; accepts a per-bend
    allowance.
  - Entities: LINE, ARC, CIRCLE, LWPOLYLINE, POLYLINE, SPLINE, TEXT, MTEXT,
    INSERT, in the XY plane.
- **Trumpf TruTops** uses its proprietary GEO format (contours plus bend
  data, material, thickness); DXF carries geometry only
  ([GEO Translator](https://www.cadstudio.cz/en/apps/geo-translator/)).

**Adopted:**
- configurable layers with defaults `OUTER` / `INNER` / `BEND_UP` /
  `BEND_DOWN` / `ETCH` / `TEXT`
- optional Delem-style labels
- bend table as separate JSON/CSV
- R2010 by default

Not verified: documented DXF conventions for Trumpf TruTops Boost and Amada
VPSS 3i / Dr.Abe; whether Amada has a native bend-data format.

---

## 8. Hard geometry

**Summary:**
- None of the open tools develops cones, B-spline bends, sharp bends,
  formed features or curved flanges. FreeCAD SheetMetal (old and new
  unfolder), build123d `dev` and Quaoar's published methods all stop short
  of these, so flatmat has to build most of the hard cases itself.
- **Sharp bends:** if one skin has a cylinder of radius ≈ t, develop that
  skin. The existing cylinder code then works unchanged, with r_in = R − t
  = 0. Only corners that are sharp on both skins need a synthetic bend
  strip.
- **B-spline bends:** `ShapeAnalysis_CanonicalRecognition` and
  `GeomConvert_SurfToAnaSurf` are exposed in OCP 7.8 and 8.0 (checked in
  the `cadquery-ocp-stubs` 7.8.1.1 and 8.0.1.1.0 wheels, not imported at
  runtime). There is a shell-level bug, and the converted surface has a new
  parameterisation, so the old pcurves are useless.

### 8.1 Hems

**What existing tools do:**

- **build123d `dev`:**
  - `HemType.FLAT/OPEN/TEARDROP/ROLLED`.
  - A flat hem is modelled with `MIN_BEND_RADIUS = 1e-3`, because r = 0
    would leave "legs lying on one another, which neither a surface nor a
    solid can represent".
  - Flat and open hems: a 180° bend with r = max(0.5·opening, 1e-3).
  - Rolled hem: up to 270° + asin(r/(r+t)).
  - Teardrop: solved with `brentq`; with zero gap the angle is
    180° + 2·atan(r/leg).
  - No hem-specific allowance: `BA = |φ|·(r_in + k·t)`.
- **FreeCAD new unfolder:** no hem logic. A hem is a cylinder with Δu ≈ π,
  and only Δu > 359.9° is rejected.
- **[CADEX MTK](https://docs.cadexsoft.com/mtk/classcadex_1_1_sheet_metal___hem_bend):**
  `SheetMetal_HemBend(radius, angle, width)` with types flattened, open,
  teardrop, rope and rolled.
- **[Creo](https://support.ptc.com/help/creo/ced_modeling/r20.7.0.0/en/ced_modeling/OSDM_Modules/SheetMetal_adm_allow.html):**
  an empirical `hem_allow` table keyed on t and the lip gap. Example for
  1.25 mm: −0.70 mm at gap 0.001, −0.20 mm at gap 1.25.
- **[Shop rule](https://sheetmetal.me/tooling-terminology/sheet-metal-hems/):**
  allowance ≈ 0.43·t, "not 100 % accurate". With the base measured to the
  outer fold, the theory gives flat = B + L − 2t + π·k·t ≈ B + L − 0.6t at
  k ≈ 0.45, which is consistent with both the 0.43t rule and Creo (our
  inference).

**Detection:**

- A cylinder with Δu ≥ π − ε is a hem candidate:
  - Δu ≈ π → flat or open hem;
  - Δu > π → teardrop or rope;
  - Δu > π with no flat leg after it → rolled.
- The legs are antiparallel planes at distance ≈ 2·r_in.
- **Touching skins:**
  - skin faces that are not adjacent but closer than ~0.1·t
    (`BRepExtrema_DistShapeShape`) with facing normals;
  - non-manifold edges (more than 2 ancestor faces);
  - `BRepCheck_Analyzer`;
  - the solid arriving as two lumps.
- **Thickness pairing needs a material-side test:** classify the midpoint
  between the two faces with `BRepClass3d_SolidClassifier` and require IN.
  Otherwise facing inner legs with a 0.002 mm gap pair up wrongly.
- **Merged skins** (a true r = 0 boolean fuse): the inner faces vanish,
  leaving a planar pair 2t apart, a half-cylinder of R ≈ t and a t-high
  step face.

**Unfolding:**

- If the inner skin exists, no special case: neutral radius r + k_hem·t,
  with a separate k_hem or a hem table.
- **Do not** use k(r/t) as r → 0: DIN's curve goes ≤ 0 below r/t ≈ 0.05.

**Refuse:** merged-skin hems (optional opt-in reconstruction splitting the
2t block at mid-plane), non-manifold contact, flat self-overlap.

### 8.2 Sharp bends (r = 0)

- **[SolidWorks](https://www.javelin-tech.com/blog/2018/02/solidworks-zero-radius-bend-values/):**
  the minimum radius is 0.001 mm, but flattening often fails at that value;
  Javelin suggests ≈ 0.025 mm.
- **Inventor:** stores r = 0 as 0.001 mm (forum, unofficial).
- **Onshape:** nothing official found.
- **Shop reality:** the radius comes from the tooling (air bending
  ≈ 15–17 % of V for cold-rolled steel, 20–22 % for 304). A bend counts as
  "sharp" below ≈ 0.63·t (rule of thumb). For the flat, replace r = 0 by
  r_eff from tooling.

**Pairing:**

- Candidate: a straight concave edge E between skin faces A and B that
  both have thickness partners. Exclude cut-out corner edges of length ≈ t
  that lie between thickness faces.
- Offset A and B by t into the material; their intersection is the outer
  virtual-sharp line L. Then:
  - (a) a convex edge within tol of L → **sharp on both sides**;
  - (b) a cylinder R ≈ t whose axis lies on E, tangent to A′ and B′ →
    **inner sharp, outer round**;
  - (c) anything else → inconsistent thickness.
- Collinear split edges are grouped into one bend.

**Unfolding:**

- With φ = π − dihedral angle, develop the inner skin and insert a strip
  along E of width
  `gap = BA(r_eff) − 2·r_eff·tan(φ/2)`, with `BA = φ·(r_eff + k·t)`.
- This equals `2t·tan(φ/2) − BD(φ)`, so bend-deduction tables plug in
  directly. With r_eff = 0 the gap is φ·k·t.
- Case (b): develop the outer skin; the existing code works with r_in = 0.
- Sharp hems (tan(φ/2) → ∞): handled as hems.

**DIN 6935 bend deduction** (de.wikipedia "Biegeverkürzung"; β = opening
angle):

- 90° < β ≤ 165°: `v = 2(r+t)·tan((180°−β)/2) − π·(180°−β)/180°·(r + t/2·k)`
- 0 < β ≤ 90°: `v = 2(r+t) − π·(180°−β)/180°·(r + t/2·k)`
- β > 165°: v ≈ 0
- `k_DIN = 0.65 + 0.5·lg(r/t)`; K = k_DIN/2 (FreeCAD's
  `_convert_to_ansi_kfactor` halves it too).
- The validity range (k = 1 for r/t > 5) is still unverified.

### 8.3 Formed features and holes

- **[CADEX MTK](https://docs.cadexsoft.com/mtk/mtk_sheet_metal_features)**
  classes:
  - `SheetMetal_Hole`
  - `SheetMetal_ComplexHole` (pierced or extruded, with radius, depth and
    axis)
  - `SheetMetal_Bead`, `SheetMetal_Louver` (depth), `SheetMetal_Bridge`
  - `SheetMetal_Cutout`, `SheetMetal_Notch` / `StraightNotch` / `VNotch`
  - `SheetMetal_Tab`, `SheetMetal_FormingFeature`

  The unfolder "ignores" beads, louvres and bridges. The relief DFM check
  is W_r ≥ T·Rt_w and L_r ≥ R_b·Rt_l (defaults 0.1 / 0.1).
- **Quaoar
  [part 7](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-7-drilled-holes) /
  [part 13](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-13-hole-mapping):**
  machined holes are recognised and stored as AAG attributes. In the flat,
  the main diameter is substituted and the countersink or counterbore is
  dropped; it is reported in JSON (`diameter`, `countersinkDiameter/Angle`,
  `...Opposite`, `faceIds`).
- **Old FreeCAD unfolder:** `handle_hole()` / `compute_replacement_circle()`
  turn cone + cylinder and counterbore stacks into a circle of the smaller
  radius.
- **SolidWorks forming tools:** the flat shows the punch, the profile, or a
  centre mark.

**Cheap detection** (per inner wire of a planar skin face, look at the
faces across its edges):

- Thickness faces that reach the matching loop on the other skin → **cut**:
  - one cylinder → plain hole;
  - cone + cylinder, coaxial, cone half-angle 41–60° → countersink;
  - cone only, through the full thickness → knife-edge countersink (flag
    it);
  - cylinder + annular plane + cylinder → counterbore.
- Tangent curved faces (e.g. a torus) → **formed feature**:
  - coaxial cylinder pair with ΔR ≈ t ending in a free edge → extruded
    hole;
  - returns to a plane parallel to the skin at offset h, with no thickness
    faces → emboss, dimple or bead (bead if elongated);
  - thickness faces inside the region → louvre, lance or bridge.

**Output:**

- Holes: the through (minimum) diameter as the cut, plus annotations.
- Extruded holes: a pre-hole from a material or tool table. The two
  published formulas, `d = D − (2H − 0.86R_d − 1.43T)` (Fabricator snippet,
  low-carbon steel) and `pre-hole = D + 1.8T − 2H` (conic), use different
  definitions of D, so neither is used directly.
- Louvre, lance, bridge: slit lines, a marker and a tool id.
- Emboss, dimple, bead: a mark on a marking layer and a tool id, no cut.

### 8.4 Cones

- **OCCT parameterisation** (`Geom_ConicalSurface.hxx`):
  `P(u,v) = O + (R + v·sinα)(cos u·X + sin u·Y) + v·cosα·Z`. v is the slant
  length; α can be negative (apex on the +Z side).
- **Development:** `ρ = v + R/sinα` (slant distance to the apex),
  `ψ = u·sinα`, flat point `(ρ cosψ, ρ sinψ)`. The sector angle is
  Δu·|sinα|; for a full cone that is 2π·r/s.
- **Neutral layer** (our derivation):
  - an offset d along the normal keeps α and moves the apex:
    `ρ' = ρ + d·cotα`, circle radius r + d·cosα, ψ unchanged;
  - so shift ρ by k·t·cotα from the concave skin.
  - Not a u-scale: u = const maps to rays, v = const to arcs about the apex
    image, and other pcurves must be sampled and refit
    (`Geom2dAPI_PointsToBSpline`).
  - r varies along a conical bend, so k(r/t) varies too; tools use a
    constant k.
- **Tool support:**
  - [SolidWorks](https://help.solidworks.com/2019/English/SolidWorks/sldworks/t_creating_sheet_metal_parts_with_conical_faces.htm)
    requires K-factor rather than bend tables for cones;
  - [Onshape](https://www.onshape.com/en/resource-center/tech-tips/cylindrical-conical-sheet-metal-parts)
    supports cylinders and cones;
  - no open tool unrolls cones (FreeCAD warns, build123d rejects).
- **Refuse:** inner and outer cones sharing an apex (thickness not
  constant); oblique or inexact cones.

### 8.5 B-spline bend surfaces

Present in OCP 7.8.1.1 and 8.0.1.1.0 (stubs inspected):

- **`ShapeAnalysis_CanonicalRecognition(shape)`**: `IsPlane` / `IsCylinder`
  / `IsCone` / `IsSphere(tol, gp_X)`, `IsLine` / `IsCircle` / `IsEllipse`,
  `GetGap()`, `GetStatus()`. Since OCCT 7.7.0; no `IsTorus`.
- **`GeomConvert_SurfToAnaSurf(S)`**: `SetConvType(Target|Simplest|MinGap)`,
  `SetTarget(GeomAbs_*)`, `ConvertToAnalytical(tol[, u0, u1, v0, v1])`,
  `Gap()`, `IsCanonical_s`, `IsSame_s`. Handles tori.
- Also: `GeomConvert_CurveToAnaCurve`, `GeomLib_IsPlanarSurface(S, 1e-7)`,
  `ShapeCustom_Surface.ConvertToAnalytical`,
  `ShapeCustom.SweptToElementary_s`, `BRepLib_FindSurface`.

**Pitfalls**
([source](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V7_8_0/src/ShapeAnalysis/ShapeAnalysis_CanonicalRecognition.cxx)):

- `SurfToAnaSurf` clips the tolerance to bbox diagonal / 1000 and measures
  the gap on a 21×21 grid.
- For faces and shells, CanonicalRecognition uses the full surface, with no
  least-squares fallback (that exists only for edges and wires). Better to
  call `SurfToAnaSurf` directly with `BRepTools.UVBounds_s(face)`.
- **Shell bug in V7_8_0 and V8_0_0:** `SetSurfParams(theTarget,
  anElemSurf, …)` instead of `anElemSurf1`, so it compares the first face
  with itself and never checks coaxiality across the shell. → Recognise per
  face and compare axes ourselves.
- `ShapeCustom_Surface` handles only B-spline and Bezier input, using
  heuristics.
- The converted surface has a new parameterisation, so the old pcurves are
  useless. Either compute (u, v) ourselves from the `gp_Cylinder` axes
  (exact, unwrapping u), or rebuild with `BRepBuilderAPI_MakeFace(surf,
  wire)` + `ShapeFix_Face`.
- Refuse non-developable faces: Gaussian curvature ≠ 0 on a grid
  (`GeomLProp_SLProps`).

### 8.6 Joggles, curved flanges, rolled tubes

- **Joggles:**
  - Detection: two cylinders with parallel axes and opposite convexity,
    either adjacent or separated by a short flat.
  - **Neutral radius per face from its own convexity:** convex R − t + k·t,
    concave R + k·t. One skin carries the concave face of bend 1 and the
    convex face of bend 2, so a global "inner skin" assumption breaks.
  - build123d's minimum jog offset is (r1 + r2)(1 − cosθ); Creo uses an
    `offset_allow` table up to an offset of 2t.
- **Curved flanges:**
  - A flange on a straight generatrix of a rolled base is fine (all
    cylinders, parallel axes).
  - A flange along a curved edge has a **torus bend**, which is not
    developable (MTK: "curved bend"). Refuse it, or approximate and flag
    it.
- **Rolled vs tube**
  ([Quaoar part 8](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-8-rolled-parts-vs-tubes)):
  - a rolled part has genus 0, a tube genus 1;
  - check Σ Δu ≈ 2π, then check topologically that thickness faces
    separate the sheet.
  - A full OCCT cylinder has a seam edge (`BRep_Tool.IsClosed_s(edge,
    face)`, two pcurves). Developing from pcurves opens it naturally, but
    the placement must not glue the seam to itself.
  - Split halves need a chosen cut edge (user choice, or one that avoids
    holes).
  - Output: a rectangle 2π(R_in + k·t) wide with the seam marked.
  - SolidWorks requires a modelled slit.

### 8.7 Bend and corner reliefs

- Relief side faces are thickness planes perpendicular to the bend axis.
  They meet the bend cylinder in v = const arcs, which develop as straight
  lines perpendicular to the bend line.
- A round relief end adds a thickness cylinder in the adjacent flat.
- **The bend line's extent is the trimmed v-range of the bend face.** A
  partial-depth relief leaves a notched cylinder face.
- Cut-outs can split a bend into several faces, which are grouped
  (coaxial, same R and φ) into one bend with segments (Quaoar part 11).
- Corner reliefs (round, square, obround, constant width) fall out of the
  face-by-face development.
- **Pitfalls:**
  - sliver faces where a relief barely clips a bend;
  - closed corners where flanges touch (the same contact checks as for
    hems).
- No relief: the bend runs to the edge with material on both sides.
  build123d raises "needs the sheet ripped rather than notched"; flatmat
  reports it as a DFM issue.

Not verified:

- Inventor and Onshape behaviour at r = 0 (forum sources only)
- the SolidWorks sharp-bend help page
- MTK hem enum identifiers and what MTK draws in the flat for holes and
  formed features
- the full Fabricator articles (403)
- the measuring convention behind the 0.43t hem rule
- the validity range of DIN's k
- Quaoar parts 2, 4, 6, 9 and 12 in detail

[Analysis Situs SMRU K-factor](https://analysissitus.org/featuresext/features_smru-k-factor.html):
default k 0.4468, with linear presets k = a·r/t + b.

## 9. Robustness on real-world STEP files

Checked against OCCT 7.8.0/7.9/8.0.1 source, the OCCT 8.0.1 STEP guide,
FreeCAD SheetMetal and build123d source, and the Quaoar posts. Anything
not confirmed is marked *[unverified]*.

### 9.1 How real exports break sheet metal recognition

| Defect | Evidence | Remedy |
|---|---|---|
| **Cylinders split into two halves** (holes, tubes, full bends) | SolidWorks has a "Split periodic faces" export option (OCCT tracker #33368); same for Creo ([CS274017](https://www.ptc.com/en/support/article/CS274017)); [NX KB](https://support.sw.siemens.com/en-US/okba/PL8010682/Imported-STEP-splits-cylindrical-faces-in-two/index.html) | `ShapeUpgrade_UnifySameDomain`, or a separate "same-domain" edge class in the graph |
| **Planes or cylinders written as B-splines** | fabrication importers reject all-B-spline solids ([McNeel / 247TailorSteel](https://discourse.mcneel.com/t/247tailorsteel-stepfiles/105631)); some writers have a "BSpline" mode ([Altair](https://2017.help.altair.com/2017/engsol/step_export_options.htm)) | canonical recognition first: `ShapeAnalysis_CanonicalRecognition` (**present in 7.8**: `IsPlane`/`IsCylinder`/`IsCone`/`IsLine`/`IsCircle` with tolerance), `GeomLib_IsPlanarSurface`; treat `Geom_SurfaceOfLinearExtrusion` of a circle as a cylinder |
| **Missing pcurves** | optional in STEP; FreeCAD writes them by default ([#2744](https://tracker.freecad.org/view.php?id=2744)); OCCT once proposed not writing them ([#25654](https://tracker.dev.opencascade.org/view.php?id=25654)) | `read.surfacecurve.mode=3`; `ShapeFix_Edge::FixAddPCurve`, `BRepLib::BuildPCurveForEdgeOnPlane` |
| **Large tolerances and gaps** | ACIS/ShapeManager tolerant modelling (Inventor, Fusion) exports bodies that "will continue to have gaps" ([Kubotek](https://help.kubotekkosmos.com/keycreator/2024SP1/A10196)). OCCT reader healing may raise tolerances up to `read.maxprecision.val` = **1 mm** by default. | measure with `ShapeAnalysis_ShapeTolerance`; flag if max vertex tolerance > ~5 % of t; consider `read.maxprecision.val` ≈ 0.1 |
| **Misaligned halves joined by a boolean** | Quaoar: a 0.24 mm shift killed recognition (§9.3) | detect thickness outliers and short step faces next to bends; report, don't guess |
| **Several solids; surfaces instead of solids** | `shell_based_surface_model` → shells in a compound ([STEP guide](https://occt3d.com/dev/doc/overview/html/occt_user_guides__step.html)); Fusion users get "surfaces" ([forum](https://forums.autodesk.com/t5/fusion-design-validate-document/export-sheet-metal-file-to-solidworks/m-p/10252190)) | iterate `TopAbs_SOLID`; with no solid, `BRepBuilderAPI_Sewing` + `ShapeFix_Solid` |
| **Inside-out solid** | wrong `same_sense` or unoriented solid | check `GProp` volume > 0, or `BRepLib::OrientClosedSolid`; otherwise convex and concave swap |
| **Faces split arbitrarily** (coplanar flange pieces, a bend split along its arc) | feature history, exporter behaviour | fourth edge class "same-domain"; merge or group |
| **Units** | the reader converts to `xstep.cascade.unit` (default mm); `STEPCAFControl_Reader::prepareUnits` sets the document unit if missing (7.8 source); `read.precision.val` is always in mm | `XCAFDoc_DocumentTool::SetLengthUnit(doc, 1, MM)` before transfer; log `FileUnits`; plausibility-check t against gauges |
| **AP242 tessellated geometry** | `read.step.tessellated` defaults to On; 2 = only where no B-rep | reject or flag faces with a null `BRep_Tool::Surface` (mesh-only) |
| **Hidden flat-pattern body** | no evidence that any major CAD system writes it by default (Inventor exports flat patterns only to SAT/DXF/DWG, [help](https://help.autodesk.com/cloudhelp/2020/ENU/Inventor-Help/files/GUID-5EE0993A-643E-4575-9D20-3803804753E1.htm)); possible through user multibody exports *[unverified]* | heuristic in §9.5 |

### 9.2 OCCT healing tools

**What the STEP reader already runs (7.8).** The `FromSTEP` sequence
(`XSTEPResource/STEP`) is a single `FixShape` with nearly every mode at −1
("fix if needed"). Self-intersection, lacking and add-pcurve fixes are −1;
`ModifyTopologyMode` is 0.

**Switching import healing off:**

- **7.8:** no clean switch. Without a configured sequence,
  `XSAlgo_AlgoContainer::ProcessShape` falls back to a default
  `ShapeFix_Shape` anyway (source).
- **7.9 / 8.0:** `XSControl_Reader::SetShapeProcessFlags(ShapeProcess::OperationsFlags())`.
  Empty flags count as "set", so nothing runs (source).
- Quaoar's `read.step.nohealing` exists only in their fork
  ([post](https://quaoar.su/blog/page/perforated-sheet-metals)): healing was
  ~90 % of import time on perforated parts, the hotspot was the wire
  self-intersection checks, and switching it off was "several times faster"
  on ~2,600 cases with no recognition loss.

| Tool | Use |
|---|---|
| `ShapeFix_Shape` | umbrella Solid → Shell → Face → Wire → Edge; for unknown input; slow on perforated parts, don't run it twice |
| `ShapeFix_Wire` | `FixReorder`, `FixConnected`, `FixLacking`, `FixDegenerated`, `FixSelfIntersection`; the last is the expensive one, switch it off for large hole counts |
| `ShapeFix_Face` | `FixMissingSeam`, `FixOrientation`, `FixSmallAreaWire`, `FixIntersectingWires`, `FixSplitFace` |
| `ShapeFix_Wireframe` | `FixWireGaps`, `FixSmallEdges` (`ModeDropSmallEdges`, `LimitAngle`); dropping small edges can delete relief edges, so leave drop off |
| `BRepLib::BuildCurves3d`, `BuildPCurveForEdgeOnPlane` | rebuild missing curves (pcurves on planes are cheap and exact); follow with `SameParameter` and `UpdateTolerances` |
| `ShapeAnalysis_ShapeTolerance` | `Tolerance(shape, mode)`, `OverTolerance`, `InTolerance`; drives the adaptive tolerances |
| `BRepCheck_Analyzer` | always; topology-only is cheap, geometry checks cost more |
| `BOPAlgo_ArgumentAnalyzer` | the most expensive; diagnostics or opt-in only |
| `BRepLib::ContinuityOfFaces`, `EncodeRegularity` | cross-check for smooth edges; the default angle 1e-10 is far too strict |

**`ShapeUpgrade_UnifySameDomain` (7.8 source):**

- Merges only neighbouring faces whose surfaces coincide (`IsSameDomain`):
  - two planar faces: `GeomLib_IsPlanarSurface(S, LinTol)` — this includes
    B-spline and Bezier surfaces — with normals parallel within AngTol and
    plane distance < LinTol;
  - two elementary surfaces that `IntPatch_ImpImpIntersection` reports as
    `TangentFaces()`;
  - a cylinder and a swept surface with equal radius and axis.
- **A plane and a tangent cylinder are never "same domain"**, so bend
  boundaries survive.
- **Exception:** a large-radius or small-angle bend stored as a B-spline
  whose sagitta R(1 − cos(θ/2)) is below LinTol passes the planarity test
  and gets merged into the flange. Pass all flange–bend edges to
  `KeepShapes` and keep `SetLinearTolerance` small (the defaults are
  1e-7 / 1e-12).
- Options: `UnifyEdges`, `UnifyFaces`, `ConcatBSplines`,
  `AllowInternalEdges` (default false), `SetSafeInputMode`. Faces in
  different shells are never merged.
- `History()` returns a `BRepTools_History` (`Modified`, `IsRemoved`),
  collected by default.
- Known bugs: [#30158](https://tracker.dev.opencascade.org/view.php?id=30158)
  (seams), [#30099](https://tracker.dev.opencascade.org/view.php?id=30099)
  (cylinders). Always run `BRepCheck` afterwards.

### 9.3 Quaoar on broken parts

- **[Another case with a broken sheet metal part](https://quaoar.su/blog/page/another-case-with-a-broken-sheet-metal-part):**
  two halves joined by a boolean, one shifted 0.24 mm, leaving a tiny
  "elevation" face near a bend, so the part was not recognised. The repair
  was manual (partition, realign, fuse). → detect and report, don't
  auto-heal.
- **[Inspection of a FreeCAD sheet metal model (RU)](http://quaoar.su/blog/page/inspekcija-odnoj-listovoj-modeli):**
  bends without a coaxial partner, inner and outer radius equal, a flange
  welded onto another lobe ("fake origami"). Result: a partial unfold with
  the unidentified faces marked.
- **[Thickness heuristic](https://quaoar.su/blog/page/on-thickness-heuristic-for-sheet-metal-parts):**
  - false positives such as a U-profile with nearly parallel planes;
  - slanted cuts make side faces "thicker" than t, so correct for the
    dihedral angle;
  - ray probes on faces and edges;
  - all heuristics in one function, no magic numbers.
- **[Planarity check](https://quaoar.su/blog/page/on-planarity-check-for-a-sheet-panel):**
  "fail gracefully when garbage comes in".
- **[Part 5](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-5):**
  - bends glued to bends add spurious graph links;
  - rotating flanges one after another accumulates angle round-off;
  - remedy: unroll each face and attach it via anchor edges (close to our
    shared-vertex placement); compute compensation once.
- **[Part 11, bend grouping](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-11-bend-grouping):**
  group by angle, axis, inner radius **and** up/down direction; geometry
  alone fails.
- **[Part 13, hole mapping](https://quaoar.su/blog/page/on-sheet-metal-unfolding-part-13-hole-mapping):**
  history keeps "dead" shapes; check every active leaf still exists.
- **[UV mapping](https://quaoar.su/blog/page/uv-mapping-for-sheet-metals):**
  `ShapeAnalysis_Surface::ValueOfUV` can return out-of-period values on
  cylinders, so normalise them. Reusing a merged face's UV space for the
  flat layout is wrong.

### 9.4 Tolerances

**Reference values** (verified in source):

- OCCT: Confusion 1e-7, Approximation 1e-6, Angular 1e-12.
- FreeCAD new unfolder: `eps` 1e-6, `eps_angular` 1e-12 rad, `fuzz` 1e-3 mm,
  `spline2arc_tol` 0.1 mm; thickness by count-mode.
- `smIsParallel`: |1 − |cos θ|| < 1e-6 (θ < ≈ 0.08°).
- Legacy unfolder: `cFaceTol` 0.002; first thickness estimate = volume /
  (area / 2).
- build123d: `TOLERANCE = 1e-6`.

**Recommendations** (engineering suggestions, not from a source):

- `tol_lin = max(1e-4 mm, 1e-6 · bbox_diag, 2 · max BRep_Tool::Tolerance of
  the entities compared)`, capped at 0.05·t. Hitting the cap marks the part
  as dirty rather than widening further.
- Thickness:
  - candidates: antiparallel plane pairs with overlapping projections, and
    coaxial ΔR;
  - bin width max(tol_lin, 0.5 % of t);
  - area-weighted mode;
  - accept if ≥ 95 % of sheet area lies within ±max(0.01 mm, 1 % of t).
- Smooth vs sharp:
  - analytic pairs are tested on the surfaces themselves, e.g. a plane and
    a cylinder: |dist(axis, plane) − R| < tol_lin and axis ⟂ normal;
  - B-splines: compare normals at 3–5 points along the edge;
  - smooth < 0.1°, sharp > 1°, anything between is near-tangent and
    resolved by the surface test.
- **A 0.5° kink must not be snapped to tangent.** That would carry a
  lateral error of L·sin 0.5° ≈ 8.7 mm per metre of flange. Develop it as a
  zero-radius 0.5° bend, with a warning.

### 9.5 Assemblies, bodies, quantities

- **Reading:**
  - create the XCAF document, set its length unit, then
    `STEPCAFControl_Reader` with `SetNameMode(True)` and `Transfer(doc)`;
  - walk from `XCAFDoc_ShapeTool.GetFreeShapes` (`IsReference` →
    `GetReferredShape`, `IsAssembly` → `GetComponents`), names from
    `TDataStd_Name`, placement from `GetLocation`, multiplied down the
    tree;
  - Python reference: build123d
    [`importers.py`](https://github.com/gumyr/build123d/blob/dev/src/build123d/importers.py);
    alternatively `XCAFPrs_DocumentExplorer` with `OnlyLeafNodes` (stable
    path ids).
- **Quantities:**
  - count instances per referred label, multiplied through nested
    assemblies;
  - for duplicated products, a geometric hash
    `V:E:F:Sh:So:area:volume:xOBB:yOBB:zOBB`
    ([post](https://quaoar.su/blog/page/hash-codes-for-identical-cad-geometries))
    plus a handedness suffix `h+` / `h−` / `h0`
    ([post](https://quaoar.su/blog/page/unique-hashes-for-mirrored-parts)).
  - Mirrored parts stay separate items: their flat is mirrored, which
    matters for countersinks and forming direction.
- **Flat-pattern body heuristic:** a planar-only solid with two large
  antiparallel faces t apart and no bend cylinders, the same t, a face area
  ≈ the developed area of the folded body, and/or a name matching
  `flat|abwicklung|unfold`. Report it as a reference flat; it is also a
  free cross-check.

### 9.6 Determinism

- **Stable:** the same file with the same OCCT build gives the same import
  order and `TopExp_Explorer` traversal. `TopExp::MapShapes` index maps
  follow insertion order, so they are the way to get stable ids
  ([forum](https://dev.opencascade.org/content/hashcode-topodsshape)).
- **Unstable:**
  - iteration over `TopTools_MapOfShape` / `DataMap` and Python sets/dicts
    keyed by shape hash (the hash is a pointer address);
  - parallel booleans and meshing;
  - healing results across OCCT versions (7.9 moved to
    `XSAlgo_ShapeProcessor`, 8.0 renamed types).
  - → pin OCP and write the OCCT version into the output metadata.
- **Stable DXF:**
  - deterministic root (largest area, then rounded centroid);
  - minimum-area-rectangle orientation with a fixed tie-break, bbox
    minimum at the origin;
  - round to 1e-6 mm, −0.0 → 0, arcs normalised to CCW;
  - entities sorted by (layer, type, coordinates);
  - `ezdxf.options.write_fixed_meta_data_for_testing = True` fixes
    `$TDUPDATE` and the GUIDs (ezdxf source).

Not verified:

- which CAD systems put a flat-pattern body into STEP
- Inventor- and Fusion-specific B-spline or tolerance behaviour
- the default of SolidWorks' "Split periodic faces"
- switching off healing in 7.8 via a custom resource file
- import healing speed across OCCT versions
- detailed AP203/214/242 differences
- Quaoar unfolding parts 7–10

## 10. Manufacturability rules for laser + press brake

Sources disagree on exact numbers, so every threshold in flatmat must be a
configurable parameter. Tags:

- **[V]** verified on the cited page
- **[S]** search snippet only; this includes all thefabricator.com figures,
  because that site returned 403
- **[GK]** general knowledge or our own derivation

Symbols: t = thickness, r = inner radius, V = die opening, Rm = tensile
strength.

### 10.1 V-die, minimum flange, radius, tonnage

**V-die selection**

- [Amada catalog](https://online.flippingbook.com/view/795100862/114/) [V]:
  V = 6t (t 0.5–2.5 mm), 8t (3–8), 10t (9–10), 12t (≥ 12), for mild steel.
- [Wila](https://wilatooling.com/knowledge-innovation/knowledge-articles/v-opening)
  [V]: 6–8t for steel up to ~3 mm; 12–15t+ for high-strength steel; below
  6t is strongly discouraged.
- [SSAB Strenx](https://www.ssab.com/en/brands-and-products/strenx/guarantee/thickness-flatness-bending)
  [V]: minimum W/t 10–16.
- [Durma](https://www.durmapress.com/press-brake-v-die-opening-chart-die-size-calculation-guide/)
  [V]: 8t mild steel, 10–12t stainless and heavy plate.

**Minimum flange** (no Trumpf, Bystronic or Salvagnini formula found):

| Source | Rule |
|---|---|
| Wila [V] | slightly more than V/2 |
| Durma [V] | ≈ 0.77·V (shop rule 0.7·V) |
| Fabricator [S] | 0.7·V standard die, ~1.1·V acute die; geometric minimum 2t + r |
| [Protolabs toolkit](https://protolabs.com/resources/design-for-sheet-metal-fabrication-toolkit) [V] | 2.5–4.0·t + r by material and thickness (table below) |
| [Protolabs one-pager](https://protolabs.com/media/roognect/sheet-metal-tolerances-one-pager-2.pdf) [V] | 4t |
| [Creo](https://support.ptc.com/help/creo/creo_plus/usascii/sheetmetal/sheetmetaldesign/About_Design_Rules.html) [V] | MIN_WALL_HEIGHT = 1.5t + r |
| [CADEX MTK defaults](https://docs.cadexsoft.com/mtk/classcadex_1_1_d_f_m_sheet_metal___analyzer_parameters) [V] | 3t (bend + flange) |
| [SendCutSend](https://sendcutsend.com/guidelines/bending) [V] | ~50–60 % of the bend length must meet the minimum; U-channel base ≥ 2 × flange; max bend angle 130° |

Protolabs, minimum flange as a multiple of t plus r:

| t (in) | Al/Cu | CRS/galv | Stainless |
|---|---|---|---|
| .024–.050 | 3.0t+r | 3.5t+r | 3.5t+r |
| .060–.090 | 2.5t+r | 2.75t+r | 2.75t+r |
| .100–.134 | 2.5t+r | 2.5t+r | 2.5t+r |
| ≥ .160 | 3.0t+r | 3.5t+r | 4.0t+r |

Steel, 90° air bend (Durma [V]; matches the
[Xometry chart](https://xometry.pro/en/?p=22145) [V]):

| t (mm) | 1 | 1.5 | 2 | 2.5 | 3 | 4 | 5 | 6 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| V | 8 | 12 | 16 | 20 | 25 | 32 | 40 | 50 | 63 | 80 |
| b_min | 5.5 | 8.5 | 11 | 14 | 17.5 | 22 | 28 | 35 | 45 | 55 |
| r_i | 1.3 | 2.0 | 2.6 | 3.3 | 4.0 | 5.0 | 6.5 | 8 | 10 | 13 |

In this table b_min ≈ 0.7·V and r_i ≈ 0.16·V. **Proposed rule:**
`b_min = max(0.7·V, k·t + r)`.

**Air-bend inside radius as a share of V (the "20 % rule")** [S]:

- 304 stainless 20–22 %
- cold-rolled steel 15–17 %
- hot-rolled pickled and oiled 12–14 %
- 5052-H32 9–14 % (sources conflict)
- other materials: scale by Rm / 60 ksi × 16 %

Durma [V] gives 16 % mild steel, 20–22 % stainless, 13–15 % aluminium. No
figure was found for 6061-T6.

**Minimum bend radius** (Xometry [V], by thickness band 1–6 / 6–12 / 12–25
mm):

- aluminium: 1t / 1.5t / 2–3t
- steel: 0.8t / 1.2t / 1.5–2.5t
- stainless: 2t / 2.5t / 3–4t

5052-H32 ≈ 1t and 6061-T6 ≈ 3–6t (crack-prone) come from forum and vendor
sources only [S]. CADEX default r ≥ 1t [V].

**Tonnage**

- [Bystronic](https://www.bystronic.com/zaf/en/news/how-calculate-press-brake-tonnage)
  [V]: tons = 575·t²/V/12 × length (inches, US tons), × (Rm_psi / 60 000) ×
  method factor (air 1, bottoming ≥ 5, coining ≥ 10).
- Metric (Durma [V]): `F = 1.42·Rm·t²·L / (1000·V)`, which comes out in kN
  with mm and MPa [GK check]. For t = 2, V = 16, Rm = 450 that is ≈ 160 kN/m,
  consistent with the charts.

### 10.2 Holes, slots, edges

**Hole to bend**

| Source | Rule |
|---|---|
| Creo MIN_CUT_TO_BEND [V] | 2.5t + r |
| [Fabcon checklist](https://fabcon.com/articles/?p=534) [V] | 2.5t + r from the bend tangent line |
| Xometry guide [S] | 2.5t + r; slots 4t + r |
| Yijin [S] | holes < 25 mm: 2t + r; larger: 2.5t + r |
| Protolabs one-pager [V] | 4t to the outside of the bend |
| CADEX defaults [V] | 3t; extruded holes 4t; notches 4t |
| [SendCutSend](https://sendcutsend.com/guidelines/bend-deformation/) [V] | ≥ V/2 from the bend centreline |

**Proposed check:** the larger of 2.5t + r from the tangent line and V/2
from the centreline.

**Other hole and edge rules**

- Hole diameter ≥ t (Protolabs, CADEX [V]). Laser-cut holes ≥ 0.5t
  (SendCutSend [V]).
- Hole to edge:
  - 2t (Creo, CADEX [V]; Xometry [S])
  - Protolabs [V]: 0.062" for t ≤ 0.036", else 0.125"
- Hole to hole: 2t (CADEX [V]); 5t (Creo MIN_DIST_BTWN_CUTS [V]).
- Smallest laser feature (arc diameter or gap between parallel edges): 1.5t
  (Creo MIN_LASER_DIM [V]).

### 10.3 Bend and corner reliefs

**When a relief is needed**
([SendCutSend](https://sendcutsend.com/blog/guide-to-designing-bend-reliefs) [V]):

- a bend that stops partway across the part needs a relief at that end
- a bend in the middle of an edge, or a window bend, needs both ends relieved
- a full-width bend needs none, unless the material is thick or stiff or the
  bend is out of square with the edge
- strongly recommended at ≥ 0.187"
- flanges meeting at a corner need a gap of ≥ 0.015" (0.4 mm)

**Relief size**

| Source | Width | Depth |
|---|---|---|
| Fabcon [V] | ≥ t | ≥ r + t |
| SendCutSend [V] | ≥ 0.5t | ≥ t + r + 0.5 mm |
| Protolabs [V] | ≥ 0.030" | total ≤ r + t, extending 0.015" past the radius |
| CADEX defaults [V] | ≥ 0.1t | ≥ 0.1r (permissive) |

No verified rule was found for the diameter of round corner reliefs.

**Detecting a missing relief from the flat** [GK]:

1. Each bend's zone strip has width BA and is bounded by the tangent lines.
2. At each end of the clipped bend line: if the strip runs out to a free
   edge, that end is fine. If flat material continues past the end, a
   cut-out ≥ t wide and reaching ≥ r + t past the tangent line is required,
   otherwise `DFM_RELIEF_MISSING`.
3. Corners: if two strips' end regions come within ~t of each other,
   require a corner relief or a gap of ≥ 0.4 mm.

### 10.4 Deformation zone, cut-outs across bends

- SendCutSend [V]:
  - the tooling contact zone is ±V/2 from the centreline, and cuts inside it
    distort
  - a bend must cross the whole bent area
  - intersecting bends without separate flanges cannot be made
  - window cut-outs ≥ 0.45" wide for V ≤ 0.984"
- [OSHCut](https://www.oshcut.com/materialdetails/aluminum-6061-t6) [V]:
  avoid holes inside the bend allowance unless they are centred on the bend.
- **Proposed keep-out band** [GK]: the BA strip plus max(2.5t + r, V/2 −
  BA/2) on each side.

### 10.5 Hems

- [Xometry](https://www.xometry.com/blog/understanding-sheet-metal-tolerances)
  [V]:
  - open hem: inside diameter ≥ t, return ≥ 4t
  - teardrop: inside diameter ≥ t, opening ≥ 0.25t, return ≥ 4t
  - avoid flat (closed) hems
  - curl: outside radius ≥ 2t; a bend must be ≥ 6t + curl radius away
- CADEX defaults [V]: open and flattened hem flange ≥ 4t; open hem radius
  0.5t.
- [Protocase](https://www.protocase.com/resources/bends/hem-bends) [V]: hem
  ≥ 0.5" (12.7 mm). Maximum thickness for hems: aluminium 0.051", CRS 0.075",
  galvanneal 0.079", stainless 0.063".
- "No hems on 6061-T6" was found only in forum and vendor reports
  (cracking; 5052-H32 preferred); not verified.

### 10.6 Grain direction

- Bending across the grain (bend line perpendicular to the rolling
  direction) is the safer case [S]. Bending along the grain cracks more,
  especially in aluminium and stainless.
- SSAB Strenx [V], minimum Ri/t transverse vs along the rolling direction:
  - class A, t < 8 mm: 1.5 vs 1.8
  - class B: 2.3 vs 2.5
  - class C: 2.8 vs 3.3
- **Proposed check** [GK]: if the rolling direction is known, bends whose
  line is within ~45° of it need r ≥ r_min_parallel. If it is unknown,
  assume the parallel case [S].

### 10.7 Bend sequence and collisions

**[BenDFM](https://arxiv.org/html/2603.13102v1)** (Ballegeer & Benoit;
[repo](https://github.com/UGent-CVAMO/bendfm), MIT) [V]:

- Taxonomy: configuration-dependent vs -independent × feasibility vs
  complexity.
- Labels:
  - punch–part and die–part collision
  - unfolding overlap
  - reorientation distance and angle (a "flip" when > 90°)
  - unfolded bounding-box area; bend count, volume, mass
- Collision simulation:
  - fixed punch: 10 mm thick, 90° tip, 300 mm long
  - fixed die: 40 mm opening
  - tested at the left, centre and right of each bend (a collision only if
    all three collide)
  - 5° steps, with downstream flanges still unbent
- Not covered: tool library, goosenecks, sequence search (sequences are
  sampled randomly), hole/flange/relief rules.
- Parameters: t 2–6 mm, r = 1–1.5t, K 0.3–0.5, flanges ≥ 75 mm.
- It covers collision and overlap only; its overlap check is directly
  reusable.

**Duflou et al. 1999** (IJPR 37(14),
[doi:10.1080/002075499190239](https://doi.org/10.1080/002075499190239)) [V]:

- branch-and-bound TSP with a dynamically updated penalty matrix
- searches **backwards** from the folded part, because first bends never
  collide and only a few bends can be the last
- collisions become hard precedence constraints
- weighted heuristics:
  - shape-defining bends after non-shape-defining ones (weight 4)
  - work from the leaves of the connectivity graph towards the root
    (weight 2)
  - shorter bends before longer ones (weight 1)
- collinear bends are grouped into one operation

**Further work** [S]: Gupta, Bourne, Kim & Krishnan, J. Manuf. Systems 17(5),
1998 (state-space search; A* is confirmed only in Amada patent US6507767);
Zhang et al. 2011 (A*); Rico et al. 2003 (channels and spirals); Márkus et
al. 2002; De Vin et al. 1994.

**Cheap heuristics** [GK]:

- A 2D cross-section collision test perpendicular to each bend: punch as a
  polygon (straight or gooseneck), die as a V block, angle swept in 5°
  steps with already-formed flanges attached (as BenDFM does).
- U-channel with base < 2 × flange, or flange taller than the punch
  clearance or gooseneck throat → "needs a gooseneck, or infeasible".
- Box: the last flange is feasible only if box depth < usable punch height
  and inner length ≥ the smallest segmented-punch width.
- Sequence hint: Duflou's backward search with the three heuristics.

### 10.8 Tooling data

No open standard for punch and die geometry was found.

- Wila provides tool data files through its app; the format and openness
  are unknown.
- Mate publishes DXF files of its die inserts [S].
- Langmuir Systems has community die profiles under CC-BY (hobby machines)
  [S].
- Radan, SigmaBEND and Almacam hold full tool libraries but are closed.

→ flatmat should define a small parametric punch/die schema, seeded with
BenDFM's defaults and a few catalog values entered by hand.

Not verified:

- any Trumpf, Bystronic or Salvagnini minimum-flange formula
- thefabricator.com page bodies
- exact air-bend radius percentages for aluminium, including 6061
- the 6061-T6 hem ban
- an authoritative 6061 vs 5052 radius table
- a rule for the diameter of round corner reliefs
- whether Gupta 1998 itself uses A*
- BenDFM label field names beyond `y_tool_collision`,
  `y_unfolding_collision` and `bbox_area_unfolded`
