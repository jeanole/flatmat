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

**v0.1 in development.** The package works end to end for planes and
cylindrical bends. Tested on:

- generated parts and random profiles (hypothesis): exact to 1e-9;
- 6 BenDFM samples: they match BenDFM's own unfolds;
- 4 real SolidWorks 2023/2026 exports.

Results are identical on OCP 7.8 and 8.0.

| | |
|---|---|
| [`docs/design.md`](docs/design.md) | algorithm, API, diagnostics, tests, hard cases, robustness, manufacturing checks, roadmap |
| [`docs/research.md`](docs/research.md) | the research behind it, with sources |
| [`spike/`](spike/README.md) | the phase-0 prototype and its measurements |

### Usage

```sh
pip install -e .            # needs Python >= 3.11; pulls cadquery-ocp-novtk
flatmat unfold part.step -k 0.44 -o part.dxf --table part.json --svg part.svg
flatmat unfold part.step --din  # K per bend from DIN 6935
flatmat unfold part.step --sharp-radius 1.6   # radius the die gives on sharp (r = 0) bends
```

```python
import flatmat as fm

sheet = fm.recognize(fm.load_step("part.step"))      # thickness, bends, flanges
pat = fm.unfold(sheet, k=fm.KTable.from_rows([(1, 0.38), (3, 0.43), (99, 0.5)]))
print(pat.size, [(b.direction, b.angle_deg, b.r_in, b.ba) for b in pat.bends])
fm.write_dxf(pat, "part.dxf")          # exact LINE / ARC / CIRCLE / SPLINE, bend lines, Delem-style labels
fm.write_bend_table(pat, "part.json")  # or .csv
```

A part flatmat cannot handle raises `fm.FlatmatError`. Its `.diagnostic` has a
code (`VARIABLE_THICKNESS`, `EDGE_FEATURE`, `UNSUPPORTED_SURFACE`,
`CYCLE_INCONSISTENT`, …) and names the faces involved.

### Tests

```sh
pip install -e ".[test]" && pytest
FLATMAT_BENDFM_DIR=/path/to/bendfm/example_files pytest tests/test_external.py
FLATMAT_PARTS_DIR=/path/to/private/step/files pytest tests/test_external.py   # never committed
```

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

- **v0.1 (implemented, not released yet):** flat faces and cylindrical bends;
  one solid of constant thickness; holes and cut-outs across bends; hems; K as
  a scalar, an r/t table or DIN 6935; exact DXF and bend table.
- **v0.2 (in progress):** sharp bends (r = 0, inside or both sides, with an
  optional tooling radius) — done; countersinks and counterbores — done;
  B-spline faces and edges that are really planes, cylinders, cones, lines
  or circles — done; manufacturability checks.
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
