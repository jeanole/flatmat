"""SolidCraft's sheet metal reference parts, rebuilt in OCP, with hand-worked flat patterns.

SolidCraft (github.com/tangency-labs/solidcraft, MIT OR Apache-2.0) keeps analytic
reference cases in ``tests/oracle/<id>/``: a recipe and the closed-form volume and area
of the part it makes, "worked out by hand", never measured in a commercial CAD product.
Its sheet metal cases are 67 to 76 (69 is not a sheet metal part there). Each one here
is the same part, built from the cross-section its ``expected.json`` describes: its
volume and area must match SolidCraft's numbers, which proves it is the same part
before flatmat unfolds it.

SolidCraft checks its folded parts only (volume, area, centre of mass, bounding box);
its flat pattern tests check bend allowances (pi/2)(R + K t) and that the flat area times
the thickness is the folded volume. ``flat_area(k)``, ``flat_size(k)`` and the bends are
worked out the same way, by hand, for flatmat's unfold.

    python tests/oracle_parts.py OUTDIR     # writes <name>.step and <name>.json per part
"""

from __future__ import annotations

import json
import math
import os
import sys
from collections.abc import Callable
from dataclasses import dataclass

from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakePrism
from OCP.GC import GC_MakeArcOfCircle
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt, gp_Vec
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain

sys.path.insert(0, os.path.dirname(__file__))
from parts import Part, make_solid, write_step  # noqa: E402


def ba(angle_deg: float, r: float, t: float, k: float) -> float:
    return math.radians(angle_deg) * (r + k * t)


@dataclass
class OraclePart:
    name: str
    oracle: str  # SolidCraft's tests/oracle directory
    t: float
    k: float  # the K-factor of SolidCraft's rule for the part
    volume: float  # SolidCraft's expected.json
    area: float
    build: Callable[[], object]
    bends: list[tuple[float, float]]  # (angle in degrees, inner radius)
    flat_area: Callable[[float], float]
    flat_size: Callable[[float], tuple[float, float]]
    holes: int = 0  # inner loops of the flat pattern
    note: str = ""


def _profile(name, t, width, flats, bends, cuts=()):
    return Part(name, t=t, width=width, flats=flats, bends=bends, cuts=list(cuts))


def _rect(p: Part) -> dict:
    """The fields of a straight profile part, a rectangle width x flat length when flat."""
    return dict(
        build=lambda: make_solid(p),
        bends=[(abs(a), r) for a, r in p.bends],
        flat_area=lambda k: p.width * p.flat_length(k),
        flat_size=lambda k: tuple(sorted((p.width, p.flat_length(k)))),
    )


def _plate():
    return BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 0), gp_Pnt(100, 60, 2.5)).Shape()


def _quarter_annulus_x(y, z, r0, r1, x0, x1):
    """A 90 degree bend about an axis along x at (y, z), the quadrant y < centre, z < centre."""
    return _annulus(gp_Pnt(x0, y, z), gp_Dir(1, 0, 0), (0, -1, 0), (0, 0, -1), r0, r1, x1 - x0)


def _quarter_annulus_y(x, z, r0, r1, y0, y1):
    """A 90 degree bend about an axis along y at (x, z), the quadrant x > centre, z < centre."""
    return _annulus(gp_Pnt(x, y0, z), gp_Dir(0, 1, 0), (1, 0, 0), (0, 0, -1), r0, r1, y1 - y0)


def _annulus(c: gp_Pnt, axis: gp_Dir, u, v, r0, r1, length):
    """The quarter annulus between radii r0 and r1 from direction u to direction v about c."""

    def at(r, a):
        return gp_Pnt(*(c.Coord(i + 1) + r * (math.cos(a) * u[i] + math.sin(a) * v[i]) for i in range(3)))

    w = BRepBuilderAPI_MakeWire()
    m = math.pi / 4
    w.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(at(r1, 0), at(r1, m), at(r1, 2 * m)).Value()).Edge())
    w.Add(BRepBuilderAPI_MakeEdge(at(r1, 2 * m), at(r0, 2 * m)).Edge())
    w.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(at(r0, 2 * m), at(r0, m), at(r0, 0)).Value()).Edge())
    w.Add(BRepBuilderAPI_MakeEdge(at(r0, 0), at(r1, 0)).Edge())
    face = BRepBuilderAPI_MakeFace(w.Wire()).Face()
    d = gp_Vec(axis.X(), axis.Y(), axis.Z()) * length
    return BRepPrimAPI_MakePrism(face, d).Shape()


def _box(x0, y0, z0, x1, y1, z1):
    return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape()


# 75: each mitred panel runs on past its bend by 0.732233047, SolidCraft's number (its
# formula text, (R + t) - (t + gap/sqrt2) with gap 0.01, does not give it: 5 - 2.5 - 2.5/sqrt2 does)
EXT75 = 0.732233047


def _bracket():
    """Oracle 75: base 75 x 35 (y 5..40), a 25 mm flange on y = 0 and a 10 mm one on x = 80,
    mitred at the corner, two d6 holes in the base."""
    pieces = [
        _box(0, 5, 0, 75, 40, 2.5),
        _quarter_annulus_x(5, 5, 2.5, 5, 0, 75),
        _box(0, 0, 5, 75 + EXT75, 2.5, 30),
        _quarter_annulus_y(75, 5, 2.5, 5, 5, 40),
        _box(77.5, 5 - EXT75, 5, 80, 40, 15),
    ]
    solid = pieces[0]
    for p in pieces[1:]:
        solid = BRepAlgoAPI_Fuse(solid, p).Shape()
    for x in (25, 55):
        hole = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x, 20, -1), gp_Dir(0, 0, 1)), 3, 5).Shape()
        solid = BRepAlgoAPI_Cut(solid, hole).Shape()
    u = ShapeUpgrade_UnifySameDomain(solid, True, True, False)
    u.Build()
    return u.Shape()


def cases() -> list[OraclePart]:
    out = []
    # 67: a plate 100 x 60 x 2.5, no bend
    out.append(OraclePart(
        "sc67_base_flange", "67-sm-base-flange", t=2.5, k=0.44, volume=15000.0, area=12800.0,
        build=_plate, bends=[], flat_area=lambda k: 6000.0, flat_size=lambda k: (60.0, 100.0),
    ))  # fmt: skip
    # 68: base y 5..60, 90 deg R 2.5, flange z 5..20; along x 100
    out.append(OraclePart(
        "sc68_edge_flange", "68-sm-edge-flange", t=2.5, k=0.44, volume=18972.621556370217, area=16057.549676223574,
        **_rect(_profile("sc68", 2.5, 100, [55, 15], [(90, 2.5)])),
    ))  # fmt: skip
    # 70: 45 deg flange: base 60 - 2.071068, flat leg 15 - 2.071068
    sb = 5 * math.tan(math.radians(22.5))
    out.append(OraclePart(
        "sc70_flange_45", "70-sm-flange-45", t=2.5, k=0.44, volume=18450.776872252372, area=15629.637035246944,
        **_rect(_profile("sc70", 2.5, 100, [60 - sb, 15 - sb], [(45, 2.5)])),
    ))  # fmt: skip
    # 71: U channel t 1.5, R 3, K 0.33: base 41, flanges 15.5; along x 80
    out.append(OraclePart(
        "sc71_custom_rule", "71-sm-custom-rule", t=1.5, k=0.33, volume=10053.716694115408, area=13896.29850950676,
        **_rect(_profile("sc71", 1.5, 80, [15.5, 41, 15.5], [(90, 3.0), (90, 3.0)])),
    ))  # fmt: skip
    # 72: Z contour flange: plate 35, wall 22.5, plate 27.5, both R 2.5, opposite ways; width 50
    out.append(OraclePart(
        "sc72_contour_flange", "72-sm-contour-flange", t=2.5, k=0.44, volume=12097.621556370217,
        area=10412.00210735098, **_rect(_profile("sc72", 2.5, 50, [35, 22.5, 27.5], [(90, 2.5), (-90, 2.5)])),
    ))  # fmt: skip
    # 73: 68 plus a flat hem: 180 deg, inner radius 0.005 (half the 0.01 hem gap), leg 5.495
    out.append(OraclePart(
        "sc73_hem", "73-sm-hem", t=2.5, k=0.44, volume=21332.046251434014, area=17992.27792617589,
        **_rect(_profile("sc73", 2.5, 100, [55, 15, 5.495], [(-90, 2.5), (-180, 0.005)])),
        note="the hem leg lies 0.01 mm from the flange it folds over",
    ))  # fmt: skip
    # 74: U channel 80 x 50 (base 40, flanges 15, R 2.5, K 0.44), a 10 x 22 slot cut across
    # the first bend while flat and refolded: it takes 5 of the base, the bend, and 22 - 5 - BA
    # of the flange. Profile coords: flange 1 along x 0..15, bend about (15, 3.75), base x 18.75.
    cut74 = 22 - 5 - ba(90, 2.5, 2.5, 0.44)
    p74 = _profile("sc74", 2.5, 80, [15, 40, 15], [(90, 2.5), (90, 2.5)],
                   cuts=[((15 - cut74, -2, 35), (25, 3.75 + 5, 45))])  # fmt: skip
    out.append(OraclePart(
        "sc74_unfold_refold", "74-sm-unfold-refold", t=2.5, k=0.44, volume=15800.304003966867,
        area=13610.326162673397, build=lambda: make_solid(p74), bends=[(90, 2.5), (90, 2.5)],
        # the slot is 22 long in the flat at K 0.44 only; at another K its bend part is BA(K)
        flat_area=lambda k: 80 * p74.flat_length(k) - 10 * (5 + ba(90, 2.5, 2.5, k) + cut74),
        flat_size=lambda k: tuple(sorted((80, p74.flat_length(k)))), holes=1,
    ))  # fmt: skip
    # 75: bracket, see _bracket
    out.append(OraclePart(
        "sc75_bracket_holes", "75-sm-bracket-holes", t=2.5, k=0.44, volume=13667.58243421114,
        area=11819.540919701749, build=_bracket, bends=[(90, 2.5), (90, 2.5)],
        flat_area=lambda k: (75 * 35 - 2 * math.pi * 9 + ba(90, 2.5, 2.5, k) * (75 + 35)
                             + 25 * (75 + EXT75) + 10 * (35 + EXT75)),
        flat_size=lambda k: tuple(sorted((85 + ba(90, 2.5, 2.5, k), 60 + ba(90, 2.5, 2.5, k)))), holes=2,
    ))  # fmt: skip
    # 76: t 3, R 3: base 54, flange 14; along x 100
    out.append(OraclePart(
        "sc76_thickness_variant", "76-sm-thickness-variant", t=3.0, k=0.44, volume=22520.575041173113,
        area=16064.12819493887, **_rect(_profile("sc76", 3.0, 100, [54, 14], [(90, 3.0)])),
    ))  # fmt: skip
    return out


def expected(c: OraclePart, k: float | None = None) -> dict:
    k = c.k if k is None else k
    return {
        "source": f"SolidCraft tests/oracle/{c.oracle} (folded volume and area); flat by hand",
        "thickness": c.t,
        "k": k,
        "folded": {"volume_mm3": c.volume, "area_mm2": c.area},
        "bends": [{"angle_deg": a, "r_in": r, "ba": ba(a, r, c.t, k)} for a, r in c.bends],
        "flat": {"area_mm2": c.flat_area(k), "size": list(c.flat_size(k)), "holes": c.holes},
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
