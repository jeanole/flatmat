"""Bent sheet metal test parts with known flat lengths (shared by the tests).

A part is a 2D profile, flats and bends along a centre line, extruded by a
width. A bend turns the centre line by a signed angle (positive = left) with a
given inner radius, so both skins are concentric arcs r and r + t. Straight
segments have the same length on both skins, which makes the expected flat
length exact:

    flat_length = sum(flats) + sum(|theta| * (r + K * t))
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakePrism
from OCP.GC import GC_MakeArcOfCircle, GC_MakeSegment
from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt, gp_Vec
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer

from flatmat._occ import ShapeList as TopTools_ListOfShape


@dataclass
class Part:
    name: str
    t: float
    width: float
    flats: list[float]
    bends: list[tuple[float, float]]  # (signed angle in degrees, inner radius)
    holes: list[tuple[str, int, float, float, float]] = field(default_factory=list)
    # ("flat", index, fraction along, z, diameter) or ("bend", index, 0.5, z, diameter)

    def flat_length(self, k: float) -> float:
        return sum(self.flats) + sum(math.radians(abs(a)) * (r + k * self.t) for a, r in self.bends)


def _profile(part: Part):
    """Centre-line turtle -> left and right boundary pieces, plus hole placements."""
    t = part.t
    p = (0.0, 0.0)
    h = 0.0
    left, right, anchors = [], [], {}

    def nrm(h):  # left normal
        return (-math.sin(h), math.cos(h))

    def off(pt, h, d):
        n = nrm(h)
        return (pt[0] + n[0] * d, pt[1] + n[1] * d)

    for i, length in enumerate(part.flats):
        q = (p[0] + length * math.cos(h), p[1] + length * math.sin(h))
        if length > 0:  # zero-length flat = two tangent bends (joggle)
            left.append(("line", off(p, h, t / 2), off(q, h, t / 2)))
            right.append(("line", off(p, h, -t / 2), off(q, h, -t / 2)))
        anchors[("flat", i)] = (p, h, length)
        p = q
        if i < len(part.bends):
            ang, r = part.bends[i]
            th = math.radians(ang)
            s = 1.0 if th > 0 else -1.0
            rc = r + t / 2
            n = nrm(h)
            c = (p[0] + s * rc * n[0], p[1] + s * rc * n[1])

            def on_arc(phi, radius, c=c, s=s, h=h):
                # point at turned heading h + phi, at `radius` from the centre
                hh = h + phi
                nn = nrm(hh)
                return (c[0] - s * radius * nn[0], c[1] - s * radius * nn[1])

            r_left = rc - s * t / 2
            r_right = rc + s * t / 2
            for side, rad in ((left, r_left), (right, r_right)):
                side.append(("arc", on_arc(0, rad), on_arc(th / 2, rad), on_arc(th, rad)))
            anchors[("bend", i)] = (c, h, th, rc, s)
            p = on_arc(th, rc)
            h += th
    return left, right, anchors


def _p3(xy, z=0.0):
    return gp_Pnt(xy[0], xy[1], z)


def make_solid(part: Part):
    left, right, anchors = _profile(part)
    edges = []
    for kind, *pts in left:
        edges.append(_edge(kind, pts))
    edges.append(_edge("line", [left[-1][-1], right[-1][-1]]))
    for kind, *pts in reversed(right):
        edges.append(_edge(kind, list(reversed(pts))))
    edges.append(_edge("line", [right[0][1], left[0][1]]))
    wire = BRepBuilderAPI_MakeWire()
    for e in edges:
        wire.Add(e)
    face = BRepBuilderAPI_MakeFace(wire.Wire()).Face()
    solid = BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, part.width)).Shape()

    if part.holes:
        tools = TopTools_ListOfShape()
        for kind, idx, frac, z, d in part.holes:
            if kind == "flat":
                p0, h, length = anchors[("flat", idx)]
                c = (p0[0] + frac * length * math.cos(h), p0[1] + frac * length * math.sin(h))
                axis = (-math.sin(h), math.cos(h))
            else:
                cen, h, th, rc, s = anchors[("bend", idx)]
                hh = h + th / 2
                axis = (-math.sin(hh), math.cos(hh))
                c = (cen[0] - s * rc * axis[0], cen[1] - s * rc * axis[1])
            base = gp_Pnt(c[0] - axis[0] * 5 * part.t, c[1] - axis[1] * 5 * part.t, z)
            ax = gp_Ax2(base, gp_Dir(axis[0], axis[1], 0))
            tools.Append(BRepPrimAPI_MakeCylinder(ax, d / 2, 10 * part.t).Shape())
        args = TopTools_ListOfShape()
        args.Append(solid)
        cut = BRepAlgoAPI_Cut()
        cut.SetArguments(args)
        cut.SetTools(tools)
        cut.Build()
        solid = cut.Shape()
    return solid


def _edge(kind, pts):
    if kind == "line":
        return BRepBuilderAPI_MakeEdge(GC_MakeSegment(_p3(pts[0]), _p3(pts[1])).Value()).Edge()
    return BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(_p3(pts[0]), _p3(pts[1]), _p3(pts[2])).Value()).Edge()


def write_step(shape, path: str) -> None:
    w = STEPControl_Writer()
    w.Transfer(shape, STEPControl_AsIs)
    if w.Write(path) != IFSelect_RetDone:
        raise RuntimeError(f"STEP write failed: {path}")


def cases() -> list[Part]:
    return [
        Part("L", t=2.0, width=100, flats=[50, 40], bends=[(90, 2.0)]),
        Part("L_down", t=2.0, width=100, flats=[50, 40], bends=[(-90, 2.0)]),
        Part("U", t=1.5, width=80, flats=[30, 60, 30], bends=[(90, 1.5), (90, 1.5)]),
        Part("Z", t=3.0, width=120, flats=[40, 30, 40], bends=[(90, 4.0), (-90, 4.0)]),
        Part("hat", t=1.0, width=60, flats=[20, 25, 40, 25, 20], bends=[(-90, 1.0), (90, 1.0), (90, 1.0), (-90, 1.0)]),
        Part("obtuse", t=2.0, width=50, flats=[40, 40], bends=[(135, 3.0)]),
        Part("acute", t=2.0, width=50, flats=[40, 40], bends=[(45, 6.0)]),
        Part("hem_open", t=1.0, width=50, flats=[40, 15], bends=[(180, 1.0)]),
        Part("joggle", t=2.0, width=80, flats=[40, 0.0, 40], bends=[(30, 3.0), (-30, 3.0)]),
        Part(
            "L_holes",
            t=2.0,
            width=100,
            flats=[50, 40],
            bends=[(90, 2.0)],
            holes=[("flat", 0, 0.4, 30.0, 8.0), ("flat", 1, 0.6, 70.0, 10.0)],
        ),
        Part(
            "L_hole_in_bend", t=2.0, width=100, flats=[50, 40], bends=[(90, 5.0)], holes=[("bend", 0, 0.5, 50.0, 6.0)]
        ),
    ]
