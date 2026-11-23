"""Sheet metal test parts with holes across bends and box corners, with flats worked out by hand.

Every part is built in OCP from primitives, and every expected flat is a closed form:

- **windows, notches and slots across a bend** (boxes and cylinders cut through the
  folded part): their walls across the bend are planes normal to the bend axis, so the
  cut is l1 + BA + l2 long in the flat at any K;
- **round holes cut in the flat and then bent** (laser cut, then folded at K_f): the hole
  wall in the bend is radial, a ruled surface between the hole's images on the two skins.
  Unfolded at K_f the hole is the circle again; at another K its bend part is stretched
  along the flat by (r + K t) / (r + K_f t);
- **round holes drilled radially through a bend** (a cylinder on the bend's bisector): the
  skin the unfolder develops sees the hole as r_s^2 sin^2(phi) + z^2 = a^2, which develops
  (s = rn * phi) to the area 4 rn r_s (E(m) - (1 - m) K(m)), m = (a / r_s)^2, with E and K
  the complete elliptic integrals. It depends on the skin: ``flat_area`` gives both;
- **open boxes** (``parts.box_solid``): base + 4 bends + 4 walls, less the corner reliefs,
  plus the corner tabs of a closed corner.

    python tests/handworked_parts.py OUTDIR     # writes <name>.step and <name>.json per part
"""

from __future__ import annotations

import json
import math
import os
import sys
from collections.abc import Callable
from dataclasses import dataclass

from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakeSolid,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_Sewing,
)
from OCP.BRepFill import BRepFill
from OCP.BRepLib import BRepLib
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCP.Geom import Geom_CylindricalSurface
from OCP.Geom2d import Geom2d_Ellipse
from OCP.gp import gp_Ax2, gp_Ax2d, gp_Ax3, gp_Dir, gp_Dir2d, gp_Elips2d, gp_Pnt, gp_Pnt2d
from OCP.ShapeFix import ShapeFix_Solid
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from OCP.TopoDS import TopoDS

sys.path.insert(0, os.path.dirname(__file__))
from parts import Part, box_solid, make_solid, write_step  # noqa: E402

K = 0.44


def ba(angle_deg: float, r: float, t: float, k: float) -> float:
    return math.radians(angle_deg) * (r + k * t)


def ellipk_ellipe(m: float) -> tuple[float, float]:
    """Complete elliptic integrals K(m), E(m) (parameter m = k^2) by the AGM."""
    a, b, s, p = 1.0, math.sqrt(1.0 - m), 0.5 * m, 0.5
    for _ in range(40):  # quadratic convergence: a handful of steps reach machine precision
        if abs(a - b) <= 1e-15 * a:
            break
        a, b, c = (a + b) / 2, math.sqrt(a * b), (a - b) / 2
        p *= 2
        s += p * c * c
    kk = math.pi / (2 * a)
    return kk, kk * (1 - s)


def radial_hole_area(rn: float, rs: float, a: float) -> float:
    """Area of a radius-a radial hole seen on a skin of radius rs, developed at radius rn."""
    m = (a / rs) ** 2
    kk, ee = ellipk_ellipe(m)
    return 4 * rn * rs * (ee - (1 - m) * kk)


@dataclass
class HandPart:
    name: str
    t: float
    build: Callable[[], object]
    bends: list[tuple[float, float]]  # (angle in degrees, inner radius)
    flat_area: Callable[[float], float | dict[str, float]]  # a dict: per developed skin ("inner", "outer")
    flat_size: Callable[[float], tuple[float, float]]
    holes: int = 0  # inner loops of the flat pattern
    k: float = K
    volume_rel: float | None = 1e-9  # flat area (K 0.5) x t = volume within this; None: not an identity
    circle: tuple[float, float] | None = None  # (radius, at K): the flat hole is this circle
    skin_cylinder: Callable[[str], float] | None = None  # radial holes: the skin's bend face area
    note: str = ""


# ---- an L bend: flange A along x 0..la (centre line y 0), bend 90 deg left about (la, rc),
# flange B at x la + rc, y rc..rc + lb; extruded along z 0..w


def _l_part(name, t, r, w, la, lb, **kw):
    return Part(name, t=t, width=w, flats=[la, lb], bends=[(90, r)], **kw)


def _window_box(p: Part, l1, l2, z0, z1):
    """A box through the bend of an L: l1 into flange A, l2 into flange B, z0..z1 across."""
    la, t, rc = p.flats[0], p.t, p.bends[0][1] + p.t / 2
    return BRepPrimAPI_MakeBox(gp_Pnt(la - l1, -t, z0), gp_Pnt(la + rc + t, rc + l2, z1)).Shape()


def _cut(solid, *tools):
    for tool in tools:
        solid = BRepAlgoAPI_Cut(solid, tool).Shape()
    u = ShapeUpgrade_UnifySameDomain(solid, True, True, False)
    u.Build()
    return u.Shape()


def _refold_tool(p: Part, u_c: float, z_c: float, a: float, k_fold: float):
    """The radial wall of a radius-a hole cut in the flat and folded into the bend of an L at
    K k_fold, centred u_c radians after the bend start: a ruled solid between its images on
    two cylinders just inside and outside the skins."""
    la, t, r = p.flats[0], p.t, p.bends[0][1]
    rc, rn = r + t / 2, r + k_fold * t
    ax = gp_Ax3(gp_Pnt(la, rc, 0), gp_Dir(0, 0, 1), gp_Dir(0, -1, 0))  # u = 0 at the bend start

    def edge(radius):
        surf = Geom_CylindricalSurface(ax, radius)
        el = gp_Elips2d(gp_Ax2d(gp_Pnt2d(u_c, z_c), gp_Dir2d(0, 1)), a, a / rn)
        e = BRepBuilderAPI_MakeEdge(Geom2d_Ellipse(el), surf).Edge()
        BRepLib.BuildCurves3d_s(e, 1e-9)
        return e, surf

    (e1, s1), (e2, s2) = edge(r - 0.5), edge(r + t + 0.5)
    sew = BRepBuilderAPI_Sewing(1e-6)
    sew.Add(BRepFill.Face_s(e1, e2))
    sew.Add(BRepBuilderAPI_MakeFace(s1, BRepBuilderAPI_MakeWire(e1).Wire()).Face())
    sew.Add(BRepBuilderAPI_MakeFace(s2, BRepBuilderAPI_MakeWire(e2).Wire()).Face())
    sew.Perform()
    fix = ShapeFix_Solid(BRepBuilderAPI_MakeSolid(TopoDS.Shell(sew.SewedShape())).Solid())
    fix.Perform()
    return fix.Solid()


def _refold_at_tangent(p: Part, z_c: float, a: float, k_fold: float):
    """A radius-a hole cut in the flat on the line where flange A meets the bend, then folded:
    half a straight hole in the flange, half a radial wall in the bend."""
    la, t, r = p.flats[0], p.t, p.bends[0][1]
    rc = r + t / 2
    big = r + 2 * t + 2
    bend_half = BRepAlgoAPI_Common(
        _refold_tool(p, 0.0, z_c, a, k_fold),
        BRepPrimAPI_MakeBox(gp_Pnt(la, -big, z_c - 2 * a), gp_Pnt(la + big, rc, z_c + 2 * a)).Shape(),
    ).Shape()
    hole = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(la, -t, z_c), gp_Dir(0, 1, 0)), a, 2 * t).Shape()
    flange_half = BRepAlgoAPI_Common(
        hole, BRepPrimAPI_MakeBox(gp_Pnt(la - 2 * a, -t, z_c - 2 * a), gp_Pnt(la, t, z_c + 2 * a)).Shape()
    ).Shape()
    return BRepAlgoAPI_Fuse(bend_half, flange_half).Shape()


def _l_cases() -> list[HandPart]:
    out = []

    # a 20 x (8 + BA + 10) window across the bend, closed on all sides
    p = _l_part("hb_window", 2.0, 3.0, 80, 40, 35)
    out.append(HandPart(
        "hb_window", p.t, lambda p=p: _cut(make_solid(p), _window_box(p, 8, 10, 30, 50)), [(90, 3.0)],
        flat_area=lambda k, p=p: p.width * p.flat_length(k) - 20 * (8 + ba(90, 3.0, 2.0, k) + 10),
        flat_size=lambda k, p=p: tuple(sorted((p.width, p.flat_length(k)))), holes=1,
        note="a 20 wide window 8 into one flange and 10 into the other",
    ))  # fmt: skip

    # a notch from the part's edge across the bend: 12 deep, 6 into each flange
    p = _l_part("hb_notch", 2.0, 3.0, 80, 40, 35)
    out.append(HandPart(
        "hb_notch", p.t, lambda p=p: _cut(make_solid(p), _window_box(p, 6, 6, -1, 12)), [(90, 3.0)],
        flat_area=lambda k, p=p: p.width * p.flat_length(k) - 12 * (6 + ba(90, 3.0, 2.0, k) + 6),
        flat_size=lambda k, p=p: tuple(sorted((p.width, p.flat_length(k)))),
        note="a bend relief cut in from the edge: the outline steps in, no hole",
    ))  # fmt: skip

    # a slot 8 wide, its round ends (r 4) 6 into each flange
    p = _l_part("hb_slot", 2.0, 2.0, 60, 30, 30)

    def slot(p=p):
        la, rc = p.flats[0], p.bends[0][1] + p.t / 2
        end_a = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(la - 6, -p.t, 30), gp_Dir(0, 1, 0)), 4, 2 * p.t).Shape()
        end_b = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(la + rc - p.t, rc + 6, 30), gp_Dir(1, 0, 0)), 4, 2 * p.t)
        return _cut(make_solid(p), _window_box(p, 6, 6, 26, 34), end_a, end_b.Shape())

    out.append(HandPart(
        "hb_slot", p.t, slot, [(90, 2.0)],
        flat_area=lambda k, p=p: p.width * p.flat_length(k) - 8 * (12 + ba(90, 2.0, 2.0, k)) - math.pi * 16,
        flat_size=lambda k, p=p: tuple(sorted((p.width, p.flat_length(k)))), holes=1,
        note="a slot across the bend: two straight sides, its round ends in the flanges",
    ))  # fmt: skip

    # a d10 hole drilled radially through the middle of an R8 bend, t 3
    p = _l_part("hb_radial", 3.0, 8.0, 60, 40, 30, holes=[("bend", 0, 0.5, 30.0, 10.0)])
    skins = {"inner": 8.0, "outer": 11.0}
    out.append(HandPart(
        "hb_radial", p.t, lambda p=p: make_solid(p), [(90, 8.0)],
        flat_area=lambda k, p=p: {
            s: p.width * p.flat_length(k) - radial_hole_area(8 + k * 3, rs, 5.0) for s, rs in skins.items()
        },
        flat_size=lambda k, p=p: tuple(sorted((p.width, p.flat_length(k)))), holes=1, volume_rel=None,
        skin_cylinder=lambda s: skins[s] * math.pi / 2 * 60 - radial_hole_area(skins[s], skins[s], 5.0),
        note="the flat hole is not a circle and its area depends on the skin developed",
    ))  # fmt: skip

    # a d8 hole laser cut in the flat at the middle of an R10 bend, then folded at K 0.44
    p = _l_part("hb_refold_mid", 2.0, 10.0, 60, 40, 30)
    th = math.pi / 2
    out.append(HandPart(
        "hb_refold_mid", p.t, lambda p=p: _cut(make_solid(p), _refold_tool(p, th / 2, 30, 4.0, K)), [(90, 10.0)],
        flat_area=lambda k, p=p: p.width * p.flat_length(k) - math.pi * 16 * (10 + k * 2) / (10 + K * 2),
        flat_size=lambda k, p=p: tuple(sorted((p.width, p.flat_length(k)))), holes=1, volume_rel=1e-6,
        circle=(4.0, K), note="cut flat, folded at K 0.44: unfolded at 0.44 it is the d8 circle again",
    ))  # fmt: skip

    # a d6 hole laser cut in the flat on the tangent line of an R6 bend, then folded at K 0.44
    p = _l_part("hb_refold_tangent", 2.0, 6.0, 60, 40, 30)
    out.append(HandPart(
        "hb_refold_tangent", p.t, lambda p=p: _cut(make_solid(p), _refold_at_tangent(p, 30, 3.0, K)), [(90, 6.0)],
        flat_area=lambda k, p=p: p.width * p.flat_length(k) - math.pi * 9 / 2 * (1 + (6 + k * 2) / (6 + K * 2)),
        flat_size=lambda k, p=p: tuple(sorted((p.width, p.flat_length(k)))), holes=1, volume_rel=1e-6,
        circle=(3.0, K), note="half the hole in the flange, half in the bend",
    ))  # fmt: skip
    return out


# ---- open boxes (parts.box_solid): base length x width, walls of height h, bends of radius r


def _box_area(k, length, width, height, t, r, gap, corner_gap=None, relief=None):
    b = ba(90, r, t, k)
    rc = r + t / 2
    a = length * width + 2 * (b + height) * width + 2 * (b + height) * (length - 2 * gap)
    if relief:
        a -= 4 * relief * (gap + relief)
    if corner_gap is not None:
        a += 4 * height * (rc + t / 2) + 4 * height * (rc + gap - t / 2 - corner_gap)
    return a


def _box_size(k, length, width, height, t, r, **_):
    b = ba(90, r, t, k)
    return tuple(sorted((length + 2 * (b + height), width + 2 * (b + height))))


def _box(name, note, window=None, **kw):
    def build():
        solid = box_solid(**kw)
        if window:  # through the bend at y = 0: l1 into the wall, l2 into the base, z0..z1
            l1, l2, z0, z1 = window
            h, t, rc = kw["height"], kw["t"], kw["r"] + kw["t"] / 2
            solid = _cut(solid, BRepPrimAPI_MakeBox(gp_Pnt(h - l1, -t, z0), gp_Pnt(h + rc + t, rc + l2, z1)).Shape())
        return solid

    def area(k):
        a = _box_area(k, **kw)
        if window:
            l1, l2, z0, z1 = window
            a -= (z1 - z0) * (l1 + ba(90, kw["r"], kw["t"], k) + l2)
        return a

    return HandPart(
        name, kw["t"], build, [(90, kw["r"])] * 4, flat_area=area, flat_size=lambda k: _box_size(k, **kw),
        holes=1 if window else 0, note=note,
    )  # fmt: skip


def _box_cases() -> list[HandPart]:
    base = dict(length=80, width=60, height=30, t=2.0, r=2.0, gap=2.0)
    return [
        _box("box_open", "corners open: the short walls stop 2 short of the long walls' bends", **base),
        _box("box_relief", "5 x 5 square reliefs cut into the base at the corners", relief=5.0, **base),
        _box(
            "box_closed", "closed corners: tabs on the walls, 1 mm apart; reliefs 5", corner_gap=1.0, relief=5.0, **base
        ),  # fmt: skip
        _box(
            "box_thin",
            "t 1.5, R 1, 120 x 70 x 25, gap 1, reliefs 3",
            length=120,
            width=70,
            height=25,
            t=1.5,
            r=1.0,
            gap=1.0,
            relief=3.0,
        ),  # fmt: skip
        _box(
            "box_deep_closed",
            "t 1, R 1.5, 60 x 40 x 40, closed corners 0.5 apart, reliefs 3",
            length=60,
            width=40,
            height=40,
            t=1.0,
            r=1.5,
            gap=1.5,
            corner_gap=0.5,
            relief=3.0,
        ),  # fmt: skip
        _box(
            "box_window",
            "box_relief with a 20 wide window across a long bend, 10 into the wall and 15 into the base",
            window=(10, 15, 20, 40),
            relief=5.0,
            **base,
        ),  # fmt: skip
    ]


def cases() -> list[HandPart]:
    return _l_cases() + _box_cases()


def expected(c: HandPart, k: float | None = None) -> dict:
    k = c.k if k is None else k
    area = c.flat_area(k)
    flat = {"area_mm2": area, "size": list(c.flat_size(k)), "holes": c.holes}
    if isinstance(area, dict):
        flat["area_mm2_by_developed_skin"] = flat.pop("area_mm2")
    if c.circle:
        flat["hole_circle_radius"] = c.circle[0]
    return {
        "source": "flatmat tests/handworked_parts.py: flat worked out by hand",
        "thickness": c.t,
        "k": k,
        "bends": [{"angle_deg": a, "r_in": r, "ba": ba(a, r, c.t, k)} for a, r in c.bends],
        "flat": flat,
        "note": c.note,
    }


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(out, exist_ok=True)
    for c in cases():
        write_step(c.build(), os.path.join(out, f"{c.name}.step"))
        with open(os.path.join(out, f"{c.name}.json"), "w") as fh:
            json.dump(expected(c), fh, indent=2)
        print(c.name)
