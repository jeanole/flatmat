"""Command line: ``flatmat unfold part.step -k 0.44 -o part.dxf --table part.json``."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Callable

from . import __version__
from .bendmodel import KDin6935
from .dfm import MATERIALS, check
from .diagnostics import FlatmatError
from .export import write_bend_table, write_dxf, write_svg
from .io import load_step
from .model import Feature
from .recognize import recognize
from .tooling import DIES, PUNCHES, Tooling
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
    u.add_argument("--material", default="steel", choices=sorted(MATERIALS), help="material for --check")
    u.add_argument("--die-v", type=float, help="die opening V in mm for --check (default: Amada rule)")
    u.add_argument("--rolling-direction", type=float, help="rolling direction in the flat, degrees from +X")
    u.add_argument("--max-force", type=float, help="press brake capacity in kN for --check")
    u.add_argument("--tooling", help="punch and die as JSON (see flatmat.tooling) for --check")
    u.add_argument("--punch", help="built-in punch for --check: " + ", ".join(PUNCHES))
    u.add_argument("--die", help="built-in die for --check: " + ", ".join(DIES) + " or V<opening>")
    u.add_argument(
        "--sequence",
        nargs="?",
        const="search",
        metavar="ORDER",
        help="with --check: find a collision-free bend order, or check ORDER (bend ids, e.g. 3,1,2)",
    )
    u.add_argument("--pre-hole", type=float, help="pre-hole diameter for extruded holes (default: neutral fibre)")
    u.add_argument(
        "--engrave",
        action="store_true",
        help="take a recess that does not go through (pocket, groove, text) as an engraving of its rim",
    )
    return p


def _tooling(a, t: float):
    if a.tooling:
        return Tooling.load(a.tooling)
    if a.punch or a.die:
        default = Tooling.default(t)
        punch = a.punch or default.punch.name
        die = a.die or default.die.name
        return Tooling.named(punch, die)
    return None


def _hole_details(f: Feature) -> str:
    return f": through Ø{f.through_diameter:g}, outer Ø{f.outer_diameter:g}"


def _oriented_details(f: Feature) -> str:
    if not f.size:
        return ""
    return f": {f.size[0]:g} x {f.size[1]:g} at {f.orientation_deg:g}°, height {f.height:g}"


# what the summary prints after "<kind> <side>" (one entry per kind in model.FEATURE_KINDS)
_DETAILS: dict[str, Callable[[Feature], str]] = {
    "countersink": _hole_details,
    "counterbore": _hole_details,
    "emboss": lambda f: f": height {f.height:g}",
    "extruded_hole": lambda f: (
        f": bore Ø{f.outer_diameter:g}, height {f.height:g}, pre-hole Ø{f.through_diameter:.3f}"
    ),
    "louvre": _oriented_details,
    "bridge": _oriented_details,
    "lance": _oriented_details,
    "engraving": lambda f: f": depth {f.depth:.3g}" + (f", Ø{f.outer_diameter:g}" if f.outer_diameter else ""),
}


def main(argv: list[str] | None = None) -> int:
    a = _parser().parse_args(argv)
    try:
        sheet = recognize(load_step(a.step, heal=a.heal), solid=a.solid, thickness=a.thickness, engrave=a.engrave)
        pat = unfold(sheet, k=KDin6935() if a.din else a.k_factor, sharp_radius=a.sharp_radius, pre_hole=a.pre_hole)
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
        print(f"  {f.kind_info.name} {f.side}{_DETAILS[f.kind](f)}")
    for d in pat.diagnostics:
        print(f"  {d.severity.value}: {d.code.value}: {d.message}", file=sys.stderr)
    report = None
    if a.check:
        try:
            tooling = _tooling(a, pat.thickness)
        except (KeyError, ValueError, OSError) as e:
            print(f"flatmat: bad tooling: {e}", file=sys.stderr)
            return 2
        seq: bool | list[int] = False
        if a.sequence == "search":
            seq = True
        elif a.sequence:
            seq = [int(x) for x in a.sequence.split(",")]
        report = check(pat, material=a.material, die_v=a.die_v, tooling=tooling, sequence=seq,
                       rolling_direction_deg=a.rolling_direction, machine_max_force_kn=a.max_force)  # fmt: skip
        tools = f", {tooling.punch.name} / {tooling.die.name}" if tooling else ""
        print(f"  check ({report.material.name}, V={report.die_v:g}{tools}): "
              f"{len(report.problems)} problems, {len(report.findings) - len(report.problems)} notes")  # fmt: skip
        for d in report.findings:
            print(f"    {d.severity.value:7} {d.code.value}: {d.message}")
        plan = report.sequence
        if plan is not None:
            if plan.feasible:
                steps = ", ".join(f"{s.bend_id}" + (" (flip)" if s.flip else "") + (" (turned)" if s.turned else "")
                                  for s in plan.steps)  # fmt: skip
                print(f"  bend order: {steps}")
            else:
                print("  bend order: none found" + (f"; longest tail {plan.order}" if plan.order else ""))
            if plan.unchecked:
                ids = ", ".join(map(str, plan.unchecked))
                print(f"  not simulated: bends {ids} (hems and lances, made with their own tools; or not separable)")
    if a.dxf:
        write_dxf(pat, a.dxf, labels=None if a.labels == "none" else a.labels, dfm=report)
    if a.table:
        write_bend_table(pat, a.table, dfm=report)
    if a.svg:
        write_svg(pat, a.svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
