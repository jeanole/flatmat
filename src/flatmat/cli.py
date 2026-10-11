"""Command line: ``flatmat unfold part.step -k 0.44 -o part.dxf --table part.json``."""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .bendmodel import KDin6935
from .diagnostics import FlatmatError
from .export import write_bend_table, write_dxf, write_svg
from .io import load_step
from .recognize import recognize
from .unfold import unfold


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="flatmat", description="Sheet metal unfolding for STEP files")
    p.add_argument("--version", action="version", version=f"flatmat {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)
    u = sub.add_parser("unfold", help="unfold a sheet metal part")
    u.add_argument("step", help="input STEP file")
    g = u.add_mutually_exclusive_group()
    g.add_argument("-k", "--k-factor", type=float, default=0.44, help="K-factor (default 0.44)")
    g.add_argument("--din", action="store_true", help="K from DIN 6935 per bend")
    u.add_argument("-o", "--dxf", help="output DXF")
    u.add_argument("--table", help="bend table, .json or .csv")
    u.add_argument("--svg", help="SVG preview")
    u.add_argument("--thickness", type=float, help="expected thickness (checked, not trusted)")
    u.add_argument("--sharp-radius", type=float, default=0.0, help="inner radius to use for sharp bends (default 0)")
    u.add_argument("--solid", type=int, help="index of the solid if the file has several")
    u.add_argument("--heal", choices=["fast", "full"], default="fast")
    u.add_argument("--labels", choices=["delem", "plain", "none"], default="delem")
    u.add_argument("--check", action="store_true", help="manufacturability checks (laser + press brake)")
    u.add_argument("--material", default="steel", help="material for --check (steel, galvanised, stainless, ...)")
    u.add_argument("--die-v", type=float, help="die opening V in mm for --check (default: Amada rule)")
    u.add_argument("--rolling-direction", type=float, help="rolling direction in the flat, degrees from +X")
    u.add_argument("--max-force", type=float, help="press brake capacity in kN for --check")
    return p


def main(argv: list[str] | None = None) -> int:
    a = _parser().parse_args(argv)
    try:
        sheet = recognize(load_step(a.step, heal=a.heal), solid=a.solid, thickness=a.thickness)
        pat = unfold(sheet, k=KDin6935() if a.din else a.k_factor, sharp_radius=a.sharp_radius)
    except FlatmatError as e:
        print(json.dumps(e.diagnostic.to_dict()), file=sys.stderr)
        return 2
    w, h = pat.size
    print(f"{a.step}: t={pat.thickness:.3f} mm, {len(pat.bends)} bends, {len(pat.holes)} holes, "
          f"flat {w:.3f} x {h:.3f} mm ({pat.k_model})")  # fmt: skip
    for b in pat.bends:
        print(f"  bend {b.id}: {b.direction:4} {b.angle_deg:7.2f}°  r_in {b.r_in:g}  K {b.k:.4f}  BA {b.ba:.4f}"
              f"  segments {len(b.segments)}")  # fmt: skip
    for f in pat.features:
        print(f"  {f.kind} {f.side}: through Ø{f.through_diameter:g}, outer Ø{f.outer_diameter:g}")
    for d in pat.diagnostics:
        print(f"  {d.severity.value}: {d.code.value}: {d.message}", file=sys.stderr)
    report = None
    if a.check:
        from .dfm import check

        report = check(pat, material=a.material, die_v=a.die_v, rolling_direction_deg=a.rolling_direction,
                       machine_max_force_kn=a.max_force)  # fmt: skip
        print(f"  check ({report.material.name}, V={report.die_v:g}): "
              f"{len(report.problems)} problems, {len(report.findings) - len(report.problems)} notes")  # fmt: skip
        for d in report.findings:
            print(f"    {d.severity.value:7} {d.code.value}: {d.message}")
    if a.dxf:
        write_dxf(pat, a.dxf, labels=None if a.labels == "none" else a.labels, dfm=report)
    if a.table:
        write_bend_table(pat, a.table, dfm=report)
    if a.svg:
        write_svg(pat, a.svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
