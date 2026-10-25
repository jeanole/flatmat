"""Bent sheet metal test parts with known flat lengths (shared by the tests).

A part is a 2D profile, flats and bends along a centre line, extruded by a
width. A bend turns the centre line by a signed angle (positive = left) with a
given inner radius, so both skins are concentric arcs r and r + t. Straight
segments have the same length on both skins, which makes the expected flat
length exact:

    flat_length = sum(flats) + sum(|theta| * (r + K * t))

A bend with r = 0 has a sharp inner corner and an outer radius t; listing its
index in ``sharp_outer`` makes the outer corner sharp too. The flat length is the
same formula with r = 0 (inner legs plus the neutral fibre through the corner).

``features`` adds countersinks and counterbores to flats:
("csk", flat, frac, z, d_through, d_top, side) with a 90° cone, or
("cbore", flat, frac, z, d_through, d_bore, side, depth); side +1 = left skin.

``formed`` adds formed features to flats, made by revolving their cross-section:
("emboss", flat, frac, z, rim_radius, height, fillet_in, fillet_top, side): a round
pad raised by ``height`` with a vertical wall; ("collar", flat, frac, z, rim_radius,
height, fillet, side): an extruded hole whose outside height (from the far skin) is
``height``. ``side`` +1 = it stands out on the left skin's side; ``fillet`` radii are
inner radii.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_Transform,
)
from OCP.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeCone,
    BRepPrimAPI_MakeCylinder,
    BRepPrimAPI_MakePrism,
    BRepPrimAPI_MakeRevol,
)
from OCP.GC import GC_MakeArcOfCircle, GC_MakeSegment
from OCP.gp import gp_Ax1, gp_Ax2, gp_Ax3, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec
from OCP.IFSelect import IFSelect_RetDone
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
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
    sharp_outer: tuple[int, ...] = ()
    features: list[tuple] = field(default_factory=list)
    cuts: list[tuple] = field(default_factory=list)  # boxes ((x0, y0, z0), (x1, y1, z1)) in profile coords
    formed: list[tuple] = field(default_factory=list)
    adds: list[tuple] = field(default_factory=list)  # boxes fused on after the cuts

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
            for side_list, rad, outer in ((left, r_left, s < 0), (right, r_right, s > 0)):
                if rad < 1e-12:
                    continue  # sharp inner corner: the neighbouring lines meet at the centre
                if outer and i in part.sharp_outer:
                    side_list.append(("corner",))
                else:
                    side_list.append(("arc", on_arc(0, rad), on_arc(th / 2, rad), on_arc(th, rad)))
            anchors[("bend", i)] = (c, h, th, rc, s)
            p = on_arc(th, rc)
            h += th
    return _corners(left), _corners(right), anchors


def _corners(side):
    """Replace ("corner",) markers by extending the neighbouring lines to their intersection."""
    out = list(side)
    while any(e[0] == "corner" for e in out):
        k = next(i for i, e in enumerate(out) if e[0] == "corner")
        (_, a0, a1), (_, b0, b1) = out[k - 1], out[k + 1]
        d1 = (a1[0] - a0[0], a1[1] - a0[1])
        d2 = (b1[0] - b0[0], b1[1] - b0[1])
        den = d1[0] * d2[1] - d1[1] * d2[0]
        u = ((b0[0] - a0[0]) * d2[1] - (b0[1] - a0[1]) * d2[0]) / den
        x = (a0[0] + u * d1[0], a0[1] + u * d1[1])
        out[k - 1 : k + 2] = [("line", a0, x), ("line", x, b1)]
    return out


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
    # a sharp inner corner leaves two lines meeting without an arc: nothing to add
    edges.append(_edge("line", [right[0][1], left[0][1]]))
    wire = BRepBuilderAPI_MakeWire()
    for e in edges:
        wire.Add(e)
    face = BRepBuilderAPI_MakeFace(wire.Wire()).Face()
    solid = BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, part.width)).Shape()

    for lo, hi in part.cuts:
        tool = BRepPrimAPI_MakeBox(gp_Pnt(*lo), gp_Pnt(*hi)).Shape()
        solid = BRepAlgoAPI_Cut(solid, tool).Shape()
    if part.holes or part.features:
        tools = TopTools_ListOfShape()
        for kind, idx, frac, z, d, big, side, *rest in part.features:
            p0, h, length = anchors[("flat", idx)]
            c = (p0[0] + frac * length * math.cos(h), p0[1] + frac * length * math.sin(h))
            n = (-math.sin(h), math.cos(h))
            eps = 0.05
            top = gp_Pnt(c[0] + side * n[0] * (part.t / 2 + eps), c[1] + side * n[1] * (part.t / 2 + eps), z)
            down = gp_Dir(-side * n[0], -side * n[1], 0)
            if kind == "csk":  # 90° included angle: the radius shrinks by 1 per unit depth
                tools.Append(
                    BRepPrimAPI_MakeCone(gp_Ax2(top, down), big / 2 + eps, d / 2, big / 2 + eps - d / 2).Shape()
                )
            else:
                tools.Append(BRepPrimAPI_MakeCylinder(gp_Ax2(top, down), big / 2, rest[0] + eps).Shape())
            base = gp_Pnt(c[0] - n[0] * 5 * part.t, c[1] - n[1] * 5 * part.t, z)
            tools.Append(BRepPrimAPI_MakeCylinder(gp_Ax2(base, gp_Dir(n[0], n[1], 0)), d / 2, 10 * part.t).Shape())
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
    for lo, hi in part.adds:
        solid = BRepAlgoAPI_Fuse(solid, BRepPrimAPI_MakeBox(gp_Pnt(*lo), gp_Pnt(*hi)).Shape()).Shape()
    for kind, idx, frac, z, *args in part.formed:
        side = args[-1]
        p0, h, length = anchors[("flat", idx)]
        c = (p0[0] + frac * length * math.cos(h), p0[1] + frac * length * math.sin(h))
        n = (-math.sin(h), math.cos(h))
        hole, add = _formed_tools(kind, part.t, *args[:-1])
        # local z = 0 is the far skin, local +z the side the feature stands out on, local x along the flat
        org = gp_Pnt(c[0] - side * n[0] * part.t / 2, c[1] - side * n[1] * part.t / 2, z)
        tr = gp_Trsf()
        tr.SetTransformation(gp_Ax3(org, gp_Dir(side * n[0], side * n[1], 0), gp_Dir(math.cos(h), math.sin(h), 0)),
                             gp_Ax3())  # fmt: skip
        hole = BRepBuilderAPI_Transform(hole, tr, True).Shape()
        add = BRepBuilderAPI_Transform(add, tr, True).Shape()
        solid = BRepAlgoAPI_Fuse(BRepAlgoAPI_Cut(solid, hole).Shape(), add).Shape()
    if part.formed or part.adds:
        u = ShapeUpgrade_UnifySameDomain(solid, True, True, False)
        u.Build()
        solid = u.Shape()
    return solid


def _formed_tools(kind: str, t: float, *a):
    """(material to cut away, material to add) for a formed feature, in its local frame.

    Local frame: the far skin at z = 0, the near skin at z = t, the feature stands
    out towards +z; x runs along the flat, y across it. The origin is the centre.
    """
    deep = (-1.0, 3 * t + 100.0)  # z range of the cutting tools

    def prism(x0, x1, y0, y1):
        return BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, deep[0]), gp_Pnt(x1, y1, deep[1])).Shape()

    def revolve(prof, angle):
        face = BRepBuilderAPI_MakeFace(_profile_wire(prof)).Face()
        return BRepPrimAPI_MakeRevol(face, gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), angle).Shape()

    def cylinder(r):
        return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0, 0, deep[0]), gp_Dir(0, 0, 1)), r, deep[1] - deep[0]).Shape()

    if kind in ("emboss", "collar"):
        rim, height, *f = a
        prof = (
            _emboss_profile(t, rim, height, f[0], f[1]) if kind == "emboss" else _collar_profile(t, rim, height, f[0])
        )
        return cylinder(rim), revolve(prof, 2 * math.pi)
    if kind == "louvre":
        # a raised half pad (y >= 0), open along its straight edge: a slit `gap` wide at y in [-gap, 0]
        rim, height, a1, a2, gap = a
        half = BRepAlgoAPI_Common(cylinder(rim), prism(-rim - 1, rim + 1, 0.0, rim + 1)).Shape()
        hole = BRepAlgoAPI_Fuse(half, prism(-rim, rim, -gap, 0.0)).Shape()
        return hole, revolve(_emboss_profile(t, rim, height, a1, a2), math.pi)
    if kind == "louvre_arc":
        # the part of a raised round pad outside a circle (radius rc, centre at x = -cx) is the hood;
        # the sheet comes back inside that circle, `gap` short of it: a slit along an arc
        rim, height, a1, a2, cx, rc, gap = a
        pad = revolve(_emboss_profile(t, rim, height, a1, a2), 2 * math.pi)
        hood = BRepAlgoAPI_Cut(pad, _cyl_at(-cx, rc, deep)).Shape()
        hole = BRepAlgoAPI_Cut(cylinder(rim), _cyl_at(-cx, rc - gap, deep)).Shape()
        return hole, hood
    if kind in ("lance_formed", "bridge_formed"):
        # a tongue cut from a raised round pad: |y| < width/2, x < tip; held only at the pad's rim at -x;
        # slits `gap` wide on the other three sides. Without a tip (tip >= rim): a strip held at both ends
        rim, height, a1, a2, width, tip, gap = a
        pad = revolve(_emboss_profile(t, rim, height, a1, a2), 2 * math.pi)
        tongue = BRepAlgoAPI_Common(pad, prism(-rim - 1, tip, -width / 2, width / 2)).Shape()
        hole = BRepAlgoAPI_Common(cylinder(rim), prism(-rim - 1, tip + gap, -width / 2 - gap, width / 2 + gap)).Shape()
        return hole, tongue
    if kind == "bridge":
        # a strip `width` wide raised by two 30° joggles at each end; slits `gap` wide along both sides
        span, width, height, r, gap = a
        ang = 30.0
        rise = height - 2 * (1 - math.cos(math.radians(ang))) * (r + t / 2)  # height from the sloped flats
        fr = rise / math.sin(math.radians(ang))
        e = 5.0
        prof = Part(
            "strip", t=t, width=width, flats=[e, fr, 1.0, fr, e], bends=[(ang, r), (-ang, r), (-ang, r), (ang, r)]
        )
        _, _, anchors = _profile(prof)
        x_end = anchors[("flat", 4)][0][0]
        # make the top flat as long as needed for the requested span
        top = span - (x_end - e - 1.0)
        prof.flats[2] = top
        _, _, anchors = _profile(prof)
        x_end = anchors[("flat", 4)][0][0]
        dx = -(e + x_end) / 2
        strip = _place_profile(make_solid(prof), dx, width, t)
        return prism(e + dx, x_end + dx, -width / 2 - gap, width / 2 + gap), strip
    if kind == "lance":
        # a tab `width` wide bent up by `angle` (radius r) along y at x = -length/2, cut free on the
        # other three sides with slits `gap` wide; the opening is as long as the flat tab (K = 0.44)
        length, width, angle, r, gap = a
        e = 5.0
        ba = math.radians(angle) * (r + 0.44 * t)
        tab = Part("tab", t=t, width=width, flats=[e, length - ba], bends=[(angle, r)])
        dx = -length / 2 - e
        return (
            prism(-length / 2, length / 2 + gap, -width / 2 - gap, width / 2 + gap),
            _place_profile(make_solid(tab), dx, width, t),
        )
    raise ValueError(kind)


def _cyl_at(x: float, r: float, deep):
    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x, 0, deep[0]), gp_Dir(0, 0, 1)), r, deep[1] - deep[0]).Shape()


def _place_profile(solid, dx: float, width: float, t: float):
    """A profile solid (profile in xy, extruded along z by width) into the local frame of _formed_tools."""
    tr = gp_Trsf()
    tr.SetValues(1, 0, 0, dx, 0, 0, -1, width / 2, 0, 1, 0, t / 2)  # x -> x, y -> z, z -> -y
    return BRepBuilderAPI_Transform(solid, tr, True).Shape()


def _arc3(c, r, a0, a1):
    def at(a):
        return (c[0] + r * math.cos(a), c[1] + r * math.sin(a))

    return (at(a0), at((a0 + a1) / 2), at(a1))


def _emboss_profile(t, rim, height, a1, a2):
    """(radius, z) cross-section of a raised round pad; far skin at z = 0, near skin at z = t."""
    pi = math.pi
    c1, ra = (rim, t + a1), rim - a1
    zw = height - a2
    c2 = (ra - t - a2, zw)
    return [
        _arc3(c1, a1, -pi / 2, -pi),
        ((ra, t + a1), (ra, zw)),
        _arc3(c2, a2 + t, 0, pi / 2),
        ((c2[0], t + height), (0, t + height)),
        ((0, t + height), (0, height)),
        ((0, height), (c2[0], height)),
        _arc3(c2, a2, pi / 2, 0),
        ((ra - t, zw), (ra - t, t + a1)),
        _arc3(c1, a1 + t, -pi, -pi / 2),
        ((rim, 0), (rim, t)),
    ]


def _collar_profile(t, rim, height, a1):
    pi = math.pi
    c1 = (rim, t + a1)
    ra, rb = rim - a1, rim - a1 - t
    return [
        _arc3(c1, a1, -pi / 2, -pi),
        ((ra, t + a1), (ra, height)),
        ((ra, height), (rb, height)),
        ((rb, height), (rb, t + a1)),
        _arc3(c1, a1 + t, -pi, -pi / 2),
        ((rim, 0), (rim, t)),
    ]


def _profile_wire(segs):
    w = BRepBuilderAPI_MakeWire()
    for sg in segs:
        pts = [gp_Pnt(r, 0, zz) for r, zz in sg]
        if len(pts) == 2:
            w.Add(BRepBuilderAPI_MakeEdge(GC_MakeSegment(*pts).Value()).Edge())
        else:
            w.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(*pts).Value()).Edge())
    return w.Wire()


def _edge(kind, pts):
    if kind == "line":
        return BRepBuilderAPI_MakeEdge(GC_MakeSegment(_p3(pts[0]), _p3(pts[1])).Value()).Edge()
    return BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(_p3(pts[0]), _p3(pts[1]), _p3(pts[2])).Value()).Edge()


def write_step(shape, path: str) -> None:
    w = STEPControl_Writer()
    w.Transfer(shape, STEPControl_AsIs)
    if w.Write(path) != IFSelect_RetDone:
        raise RuntimeError(f"STEP write failed: {path}")


def sharp_cases() -> list[Part]:
    return [
        Part("sharp_in_L", t=2.0, width=80, flats=[50, 40], bends=[(90, 0.0)]),
        Part("sharp_both_L", t=2.0, width=80, flats=[50, 40], bends=[(90, 0.0)], sharp_outer=(0,)),
        Part("sharp_both_U", t=1.5, width=60, flats=[30, 50, 30], bends=[(90, 0.0), (90, 0.0)], sharp_outer=(0, 1)),
        Part("sharp_both_Z", t=3.0, width=70, flats=[40, 30, 40], bends=[(90, 0.0), (-90, 0.0)], sharp_outer=(0, 1)),
        Part("sharp_in_Z", t=3.0, width=70, flats=[40, 30, 40], bends=[(90, 0.0), (-90, 0.0)]),
        Part("sharp_both_obtuse", t=2.0, width=50, flats=[40, 40], bends=[(120, 0.0)], sharp_outer=(0,)),
        Part("sharp_mixed", t=2.0, width=60, flats=[40, 30, 40], bends=[(90, 0.0), (90, 2.0)], sharp_outer=(0,)),
    ]


def feature_cases() -> list[Part]:
    return [
        Part(
            "csk_L",
            t=3.0,
            width=80,
            flats=[60, 40],
            bends=[(90, 2.0)],
            features=[("csk", 0, 0.4, 20.0, 6.0, 10.0, +1), ("csk", 0, 0.4, 60.0, 6.0, 10.0, -1)],
        ),
        Part(
            "cbore_L",
            t=4.0,
            width=80,
            flats=[60, 40],
            bends=[(90, 2.0)],
            features=[("cbore", 0, 0.4, 20.0, 5.0, 9.0, +1, 2.0), ("cbore", 1, 0.5, 60.0, 5.0, 9.0, -1, 2.0)],
        ),
    ]


def formed_cases() -> list[Part]:
    """t = 2, an L with r = 2; the base is the left skin's flat 0 (50 long), 100 wide."""
    L = dict(t=2.0, width=100, flats=[50, 40], bends=[(90, 2.0)])
    return [
        Part("emboss_up", formed=[("emboss", 0, 0.5, 50.0, 10.0, 5.0, 1.0, 1.0, +1)], **L),
        Part("emboss_down", formed=[("emboss", 0, 0.5, 50.0, 10.0, 5.0, 1.0, 1.0, -1)], **L),
        Part("collar_up", formed=[("collar", 0, 0.5, 50.0, 6.0, 4.5, 1.0, +1)], **L),
        Part(
            "collar_down_two",
            formed=[("collar", 0, 0.4, 30.0, 8.0, 8.0, 1.0, -1), ("collar", 1, 0.5, 70.0, 7.0, 6.0, 1.5, +1)],
            **L,
        ),
    ]


def cut_formed_cases() -> list[Part]:
    """Louvres, bridges and lances on the 60 mm base of an L (t = 2, r = 2, 100 wide); slits 0.5 wide.

    ("louvre", flat, frac, z, rim_radius, height, fillet_in, fillet_top, gap, side): a raised
    half pad open along its straight edge; ("bridge", flat, frac, z, span, width, height, r,
    gap, side): a strip raised by two joggles, slit along both sides; ("lance", flat, frac,
    z, length, width, angle, r, gap, side): a tab bent out of an opening as long as its flat.
    Cut from a raised round pad (rim, height, fillet_in, fillet_top): ("louvre_arc", ...,
    cx, rc, gap, side): the hood outside a circle of radius rc centred cx before the pad's
    centre, slit along that circle; ("lance_formed" / "bridge_formed", ..., width, tip, gap,
    side): the pad's strip |y| < width/2 up to x = tip, held at the rim at -x (a tongue), or
    at both ends if tip >= rim.
    """
    L = dict(t=2.0, width=100, flats=[60, 40], bends=[(90, 2.0)])
    return [
        Part("louvre_up", formed=[("louvre", 0, 0.5, 50.0, 10.0, 5.0, 1.0, 1.0, 0.5, +1)], **L),
        Part("louvre_down", formed=[("louvre", 0, 0.5, 30.0, 8.0, 4.5, 1.0, 1.0, 0.5, -1)], **L),
        Part("bridge", formed=[("bridge", 0, 0.5, 50.0, 30.0, 10.0, 4.0, 2.0, 0.5, +1)], **L),
        Part("lance", formed=[("lance", 0, 0.5, 50.0, 25.0, 12.0, 30.0, 2.0, 0.5, +1)], **L),
        Part("lance_90", formed=[("lance", 0, 0.4, 30.0, 20.0, 10.0, 90.0, 2.0, 0.5, -1)], **L),
        # cut from a raised round pad (rim 12, 5 high, fillets 1): stretched, not developable
        Part("louvre_arc", formed=[("louvre_arc", 0, 0.5, 50.0, 12.0, 5.0, 1.0, 1.0, 20.0, 25.0, 0.5, +1)], **L),
        Part("lance_formed", formed=[("lance_formed", 0, 0.5, 50.0, 12.0, 5.0, 1.0, 1.0, 10.0, 6.0, 0.5, -1)], **L),
        Part("bridge_formed", formed=[("bridge_formed", 0, 0.5, 50.0, 12.0, 5.0, 1.0, 1.0, 8.0, 13.0, 0.5, +1)], **L),
    ]


def dfm_cases() -> list[Part]:
    """Parts with one deliberate manufacturability problem each (t = 2, r = 2, V = 12).

    Profile coordinates: the base runs along +x at y in [-1, 1]; the first bend (left,
    90°, r = 2) has its tangent line at x = 50 and the flange rises at x in [52, 54].
    """
    L = dict(t=2.0, width=100, bends=[(90, 2.0)])
    return [
        Part("dfm_ok", flats=[50, 40], **L),
        Part("dfm_short_flange", flats=[50, 3], **L),
        Part("dfm_hole_near_bend", flats=[50, 40], holes=[("flat", 0, 0.92, 50.0, 4.0)], **L),
        Part("dfm_small_hole", flats=[50, 40], holes=[("flat", 0, 0.4, 50.0, 1.2)], **L),
        Part("dfm_hole_near_edge", flats=[50, 40], holes=[("flat", 0, 0.4, 3.5, 4.0)], **L),
        Part(
            "dfm_holes_close",
            flats=[50, 40],
            holes=[("flat", 0, 0.4, 40.0, 6.0), ("flat", 0, 0.4, 47.0, 6.0)],
            **L,
        ),
        # flange and bend removed for z > 60: the bend ends inside the base, no relief
        Part("dfm_no_relief", flats=[50, 40], cuts=[((50.0, -2.0, 60.0), (60.0, 50.0, 101.0))], **L),
        # the same with a relief slot 2.5 wide, 4.5 deep into the base
        Part(
            "dfm_relief",
            flats=[50, 40],
            cuts=[((50.0, -2.0, 60.0), (60.0, 50.0, 101.0)), ((45.5, -2.0, 57.5), (50.5, 2.0, 62.5))],
            **L,
        ),
        # a 1 mm slot through the base, 20 mm long
        Part("dfm_narrow_slot", flats=[50, 40], cuts=[((10.0, -2.0, 40.0), (30.0, 2.0, 41.0))], **L),
        Part("dfm_short_hem", t=1.0, width=60, flats=[40, 1.5], bends=[(180, 1.0)]),
    ]


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


def box_solid(length: float, width: float, height: float, t: float, r: float, gap: float, corner_gap=None,
              relief=None):  # fmt: skip
    """An open box: a base length x width, four walls of height ``height`` bent up with radius r.

    Two U profiles share the base. The long walls run the full width; the short walls
    stop ``gap`` short of the long walls' bend tangent lines, which leaves a corner
    relief. With ``corner_gap`` the short walls get tabs that reach to within that
    distance of the long walls' inner faces (a closed corner).
    The base's mid-plane is x = height + r + t/2; the walls stand towards -x.
    """
    rc = r + t / 2
    u1 = make_solid(Part("box_u1", t=t, width=width, flats=[height, length, height], bends=[(90, r), (90, r)]))
    u2 = make_solid(
        Part("box_u2", t=t, width=length - 2 * gap, flats=[height, width, height], bends=[(90, r), (90, r)])
    )
    tr = gp_Trsf()
    # local x -> x, local y (along the base) -> z, local z (extrusion) -> -y; a rotation
    tr.SetValues(1, 0, 0, 0, 0, 0, -1, rc + length - gap, 0, 1, 0, -rc)
    u2 = BRepBuilderAPI_Transform(u2, tr, True).Shape()
    solid = BRepAlgoAPI_Fuse(u1, u2).Shape()
    if corner_gap is not None:
        tabs = []
        for yc in (0.0, length + 2 * rc):  # the long walls cover the short walls' ends
            tabs.append(((0, yc - t / 2, -rc - t / 2), (height, yc + t / 2, 0)))
            tabs.append(((0, yc - t / 2, width), (height, yc + t / 2, width + rc + t / 2)))
        for zc in (-rc, width + rc):  # the short walls reach to within corner_gap of them
            for y0, y1 in ((t / 2 + corner_gap, rc + gap), (rc + length - gap, length + 2 * rc - t / 2 - corner_gap)):
                tabs.append(((0, y0, zc - t / 2), (height, y1, zc + t / 2)))
        for lo, hi in tabs:
            solid = BRepAlgoAPI_Fuse(solid, BRepPrimAPI_MakeBox(gp_Pnt(*lo), gp_Pnt(*hi)).Shape()).Shape()
    if relief:
        x = height + rc
        for y0, y1 in ((rc, rc + gap + relief), (rc + length - gap - relief, rc + length)):
            for z0, z1 in ((0, relief), (width - relief, width)):
                solid = BRepAlgoAPI_Cut(
                    solid, BRepPrimAPI_MakeBox(gp_Pnt(x - t, y0, z0), gp_Pnt(x + t, y1, z1)).Shape()
                ).Shape()
    u = ShapeUpgrade_UnifySameDomain(solid, True, True, False)
    u.Build()
    return u.Shape()
