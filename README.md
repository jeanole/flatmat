# flatmat

Sheet metal recognition and unfolding for STEP files, in Python on
OpenCascade.

flatmat is meant to take a bent sheet metal part from a STEP file and:

1. check that it really is sheet metal and measure its thickness;
2. recognise flanges and bends (inner radius, angle, up/down, axis), including
   bends interrupted by cut-outs;
3. unfold it with a K-factor model (scalar, r/t table or DIN 6935) into an
   exact flat pattern: lines, arcs and splines, with holes that cross bends;
4. write a DXF for laser cutting and press-brake software, plus a bend table
   (JSON/CSV);
5. optionally check manufacturability: minimum flange, hole-to-bend
   distance, reliefs and similar rules.

Every result carries provenance: each flat edge maps back to the STEP edge it
came from. Anything flatmat cannot handle is refused with a structured
diagnostic that names the faces involved, never silently unfolded wrong.

## Status

**Design phase.** There is no library code yet, only a phase-0 spike that
runs the core algorithm end to end and checks it:

- 11 generated parts come out exact to ≤ 3.5e-14.
- The BenDFM samples match their own unfolds.
- OCP 7.8 and 8.0 give identical results.
- Typical parts take 10–90 ms; a part with 1000 holes takes 1.5 s, and the
  time grows linearly.

See [`spike/README.md`](spike/README.md).

- [`docs/design.md`](docs/design.md): algorithm, API and data model, test and
  benchmark plan, scope per version.
- [`docs/research.md`](docs/research.md): the research behind it, with sources:
  existing tools, how the build123d, FreeCAD and Analysis Situs unfolders
  work, Python vs C++ vs Rust, licensing, test data, K-factor standards, DXF
  conventions, hard geometry, robustness and manufacturing rules.

## Why another unfolder

No standalone, permissively licensed library does STEP → bend table → DXF
today:

- build123d's new sheet metal module (Apache-2.0) has not been released yet.
  It supports only planes and cylinders, and outputs no bend lines or bend
  table.
- FreeCAD SheetMetal (LGPL) needs the FreeCAD runtime.
- The commercial SDKs (CADEX MTK, the OCCT Sheet Metal component) are priced
  per quote.

## Planned scope

- **v0.1:** flat faces and cylindrical bends; one solid of constant
  thickness; bends whose face is split in two; holes and cut-outs across
  bends; hems; DXF and bend table.
- **v0.2:** recognising B-spline faces as cylinders or planes; sharp bends;
  countersinks and counterbores; manufacturability checks.
- **Later:** cones, other formed features, tubes, assemblies, bend-sequence
  feasibility.

## Stack

Python ≥ 3.11 on [`cadquery-ocp-novtk`](https://pypi.org/project/cadquery-ocp-novtk/)
(OCCT 8.0, also tested on 7.8) and [ezdxf](https://ezdxf.mozman.at/). No
FreeCAD and no build123d at runtime.

## License

[Apache-2.0](LICENSE).

flatmat uses Open CASCADE Technology (OCCT) through the OCP bindings. OCCT
is licensed under the LGPL-2.1 with the Open CASCADE exception; its source is
available at <https://github.com/Open-Cascade-SAS/OCCT>.
