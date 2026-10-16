"""Sheet metal recognition (docs/design.md §1.1–1.3).

Faces and their adjacency form an attributed adjacency graph. Edges along which two
faces meet tangentially are SMOOTH; the connected components over SMOOTH edges are the
two skins (the two largest) and the thickness faces. Bends are the cylinders of the
reference skin.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field, replace
from typing import Any, Literal, cast

from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.BRepLProp import BRepLProp_SLProps
from OCP.BRepTools import BRepTools
from OCP.GeomAbs import (
    GeomAbs_Circle,
    GeomAbs_Cone,
    GeomAbs_Cylinder,
    GeomAbs_Line,
    GeomAbs_Plane,
    GeomAbs_Sphere,
    GeomAbs_Torus,
)
from OCP.gp import gp_Pnt, gp_Pnt2d, gp_Vec
from OCP.GProp import GProp_GProps
from OCP.IntTools import IntTools_FClass2d
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_IN, TopAbs_REVERSED, TopAbs_VERTEX
from OCP.TopExp import TopExp

from ._occ import EdgeFaceMap, ShapeMap, as_edge, as_face
from .bendmodel import extruded_hole_prehole
from .canonical import canonicalize
from .diagnostics import Code, Diagnostic, refuse, warn
from .io import StepModel
from .model import Bend, BendSegment, FaceInfo, Feature, Flange, SheetPart, Tolerances

# ------------------------------------------------------------------------------ internals


@dataclass
class _Face:
    idx: int
    shape: Any
    kind: str
    area: float
    normal: gp_Vec  # outward (face orientation applied), at the UV centre
    mirror: float  # +1 if the surface normal equals the outward normal, else -1
    uv: tuple[float, float, float, float]
    ad: Any  # BRepAdaptor_Surface, one per face
    radius: float = 0.0
    axis_p: Any = None  # gp_Pnt for cylinders
    axis_d: Any = None  # gp_Vec for cylinders and cones
    semi_angle: float = 0.0  # cones
    minor: float = 0.0  # tori
    _cls: Any = None

    def inside(self, u: float, v: float) -> bool:
        # one classifier per face: BRepClass_FaceClassifier rebuilds per call (~9 ms on
        # a face with 1000 holes against 0.8 ms here)
        if self._cls is None:
            self._cls = IntTools_FClass2d(self.shape, 1e-7)
        return self._cls.Perform(gp_Pnt2d(u, v)) == TopAbs_IN


@dataclass
class _BendFace:
    face: int
    bend_id: int
    r_in: float
    concave: bool  # reference skin is the inside of the bend at this face
    angle: float  # rad


@dataclass
class _SharpEdge:
    """A sharp bend edge between two faces of the reference skin."""

    edge: int
    faces: tuple[int, int]
    angle: float  # deflection, rad
    concave: bool  # reference skin is the inside of the bend here
    length: float
    vertices: tuple[int, int]  # vertex map indices
    bend_id: int = 0


@dataclass
class Context:
    """Everything unfold() needs; kept out of the public model."""

    solid: Any
    fmap: Any
    vmap: Any
    faces: list[_Face]
    edges: list[tuple[Any, tuple[int, ...]]]
    adj: dict[int, list[tuple[int, int]]]  # face -> [(neighbour, edge index)] over SMOOTH edges
    skin_a: set[int]
    root: int
    thickness: float
    bend_faces: dict[int, _BendFace] = field(default_factory=dict)
    sharp_edges: dict[int, _SharpEdge] = field(default_factory=dict)
    features: dict[int, tuple[int, ...]] = field(default_factory=dict)  # feature id -> faces
    feature_through_r: dict[int, float] = field(default_factory=dict)
    feature_on_a: dict[int, bool] = field(default_factory=dict)
    feature_kind: dict[int, str] = field(default_factory=dict)
    feature_edges: dict[int, set[int]] = field(default_factory=dict)  # edges of a formed feature's opening
    feature_cuts: dict[int, list[list[int]]] = field(default_factory=dict)  # its slits, as runs of edges
    feature_gap: dict[int, float] = field(default_factory=dict)  # their width
    developed: set[int] = field(default_factory=set)  # lances laid flat by bending
    feature_free: dict[int, tuple[int, list[int]]] = field(default_factory=dict)  # host face, free edges
    # rims of engravings not in a face of the reference skin (the far skin, islands):
    # feature id -> (the reference-skin plane to project onto, edges of the rim in wire order)
    projected_rims: dict[int, tuple[int, list[int]]] = field(default_factory=dict)
    tol: Tolerances = field(default_factory=Tolerances)


def _props(face, ad, u, v):
    pr = BRepLProp_SLProps(ad, u, v, 1, 1e-7)
    n = gp_Vec(pr.Normal())
    sn = pr.D1U().Crossed(pr.D1V())
    if face.Orientation() == TopAbs_REVERSED:
        n.Reverse()
    return n, (1.0 if sn.Dot(n) > 0 else -1.0)


def pcurve(e, f) -> BRepAdaptor_Curve2d:
    """BRepAdaptor_Curve2d, refusing instead of crashing when the pcurve is missing.

    On a missing pcurve the adaptor raises on OCP 7.8 and segfaults on 8.0.
    """
    if BRep_Tool.CurveOnSurface_s(e, f, 0.0, 0.0) is None:
        raise refuse(Code.MISSING_PCURVE, "an edge has no pcurve on its face")
    return BRepAdaptor_Curve2d(e, f)


def _face_table(solid) -> tuple[Any, list[_Face]]:
    fmap = ShapeMap()
    TopExp.MapShapes_s(solid, TopAbs_FACE, fmap)
    faces = []
    for i in range(1, fmap.Extent() + 1):
        f = as_face(fmap.FindKey(i))
        ad = BRepAdaptor_Surface(f, True)
        u0, u1, v0, v1 = BRepTools.UVBounds_s(f)
        n, mirror = _props(f, ad, (u0 + u1) / 2, (v0 + v1) / 2)
        g = GProp_GProps()
        BRepGProp.SurfaceProperties_s(f, g)
        kind = {
            GeomAbs_Plane: "plane",
            GeomAbs_Cylinder: "cylinder",
            GeomAbs_Cone: "cone",
            GeomAbs_Torus: "torus",
            GeomAbs_Sphere: "sphere",
        }.get(ad.GetType(), "other")
        rec = _Face(i, f, kind, g.Mass(), n, mirror, (u0, u1, v0, v1), ad)
        if kind == "cylinder":
            c = ad.Cylinder()
            rec.radius = c.Radius()
            rec.axis_p = c.Location()
            rec.axis_d = gp_Vec(c.Axis().Direction())
        elif kind == "cone":
            k = ad.Cone()
            rec.radius = k.RefRadius()
            rec.axis_p = k.Location()
            rec.axis_d = gp_Vec(k.Axis().Direction())
            rec.semi_angle = k.SemiAngle()
        elif kind == "torus":
            tor = ad.Torus()
            rec.radius = tor.MajorRadius()
            rec.minor = tor.MinorRadius()
            rec.axis_p = tor.Location()
            rec.axis_d = gp_Vec(tor.Axis().Direction())
        elif kind == "sphere":
            sph = ad.Sphere()
            rec.radius = sph.Radius()
            rec.axis_p = sph.Location()
        faces.append(rec)
    return fmap, faces


def _edge_table(solid, fmap) -> list[tuple[Any, tuple[int, ...]]]:
    efm = EdgeFaceMap()
    TopExp.MapShapesAndAncestors_s(solid, TopAbs_EDGE, TopAbs_FACE, efm)
    edges = []
    for i in range(1, efm.Extent() + 1):
        e = as_edge(efm.FindKey(i))
        lst = efm.FindFromIndex(i)
        # Iterating an OCP shape list from Python costs ~0.3 ms per list; an edge has at
        # most two faces, so First()/Last() cover the manifold case.
        if lst.Size() <= 2:
            idx = sorted({fmap.FindIndex(lst.First()), fmap.FindIndex(lst.Last())})
        else:
            idx = sorted({fmap.FindIndex(f) for f in lst})
        edges.append((e, tuple(idx)))
    return edges


def _is_smooth(e, f1: _Face, f2: _Face, tol: Tolerances) -> bool:
    c = BRepAdaptor_Curve(e)
    a, b = c.FirstParameter(), c.LastParameter()
    worst = 0.0
    for s in (0.25, 0.5, 0.75):
        tm = a + s * (b - a)
        n = []
        for f in (f1, f2):
            uv = pcurve(e, f.shape).Value(tm)
            n.append(_props(f.shape, f.ad, uv.X(), uv.Y())[0])
        worst = max(worst, n[0].Angle(n[1]))
    return worst < tol.smooth_angle


def _components(n: int, links) -> list[list[int]]:
    parent = list(range(n + 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in links:
        parent[find(a)] = find(b)
    groups = defaultdict(list)
    for i in range(1, n + 1):
        groups[find(i)].append(i)
    return list(groups.values())


def _canonical_dir(v: gp_Vec, digits: int = 4) -> tuple[tuple[float, float, float], gp_Vec]:
    """Unit direction with a canonical sign: rounded first, so -1e-17 noise cannot flip it."""
    d = v.Normalized()
    r = [round(c, digits) + 0.0 for c in (d.X(), d.Y(), d.Z())]
    first = next(c for c in r if c != 0.0)
    if first < 0:
        d.Reverse()
        r = [-c + 0.0 for c in r]
    return (r[0], r[1], r[2]), d


def _line_key(p: gp_Pnt, v: gp_Vec) -> tuple[float, ...]:
    """A hashable key for the infinite line through p along v."""
    key, d = _canonical_dir(v)
    o = gp_Vec(p.XYZ())
    foot = o - d.Multiplied(o.Dot(d))
    return (*key, *(round(c, 3) + 0.0 for c in (foot.X(), foot.Y(), foot.Z())))


def _axis_key(f: _Face) -> tuple[float, ...]:
    return _line_key(f.axis_p, f.axis_d)


def _plane_uv(f: _Face, p: gp_Pnt) -> tuple[float, float]:
    pos = f.ad.Plane().Position()
    d = gp_Vec(pos.Location(), p)
    return d.Dot(gp_Vec(pos.XDirection())), d.Dot(gp_Vec(pos.YDirection()))


def _inner_point(f: _Face):
    u0, u1, v0, v1 = f.uv
    for su in (0.5, 0.25, 0.75, 0.1, 0.9):
        for sv in (0.5, 0.25, 0.75, 0.1, 0.9):
            u, v = u0 + su * (u1 - u0), v0 + sv * (v1 - v0)
            if f.inside(u, v):
                return f.ad.Value(u, v)
    return None


def _thickness(faces: list[_Face], skins: list[set[int]], tol: Tolerances, max_planes: int = 12):
    cands: list[tuple[float, float]] = []
    # coaxial cylinder pairs, bucketed by axis line (all pairs was O(n^2) on hole walls)
    buckets = defaultdict(list)
    for f in faces:
        if f.kind == "cylinder" and (f.idx in skins[0] or f.idx in skins[1]):
            buckets[_axis_key(f)].append(f)
    for group in buckets.values():
        for i, a in enumerate(group):
            for b in group[i + 1 :]:
                if (a.idx in skins[0]) != (b.idx in skins[0]):
                    dr = abs(a.radius - b.radius)
                    if dr > tol.linear:
                        cands.append((dr, min(a.area, b.area)))
    # plane pairs, analytically: distance to antiparallel planes of the other skin,
    # with the foot inside that face (a ray cast against the solid was O(n^2))
    planes = [sorted((faces[i - 1] for i in s if faces[i - 1].kind == "plane"), key=lambda f: -f.area) for s in skins]
    for side in (0, 1):
        for f in planes[side][:max_planes]:
            p = _inner_point(f)
            if p is None:
                continue
            best = None
            for g in planes[1 - side][:max_planes]:
                if g.normal.Angle(f.normal) < math.pi - 1e-3:
                    continue
                d = abs(gp_Vec(g.ad.Plane().Location(), p).Dot(g.normal) / g.normal.Magnitude())
                if d <= tol.linear or (best is not None and d >= best):
                    continue
                foot = gp_Pnt(p.XYZ() - f.normal.Normalized().XYZ() * d)
                if g.inside(*_plane_uv(g, foot)):
                    best = d
            if best is not None:
                cands.append((best, f.area))
    if not cands:
        raise refuse(Code.NO_THICKNESS, "no pair of opposite faces found")
    bins: dict[float, float] = defaultdict(float)
    for v, w in cands:
        bins[round(v, 3)] += w
    t0 = max(bins, key=lambda k: bins[k])
    band = max(0.01, tol.thickness_rel * t0)
    near = [(v, w) for v, w in cands if abs(v - t0) <= band]
    t = sum(v * w for v, w in near) / sum(w for _, w in near)
    conf = sum(w for _, w in near) / sum(w for _, w in cands)
    return t, conf


# ------------------------------------------------------------------------------ partners


def _plane_partners(faces: list[_Face], tol: Tolerances) -> tuple[float | None, dict[int, int]]:
    """Global thickness estimate and, for every plane, the opposite plane at that distance.

    Planes are bucketed by their (unsigned) normal direction; the partner of f is the
    nearest antiparallel plane *behind* f (into the material) whose face contains the
    foot of an inner point of f. Skins are not known yet: sharp bends split them.
    """
    planes = [f for f in faces if f.kind == "plane"]
    buckets: dict[tuple, list[_Face]] = defaultdict(list)
    for f in planes:
        buckets[_canonical_dir(f.normal)[0]].append(f)
    found: dict[int, tuple[float, int]] = {}
    for group in buckets.values():
        if len(group) < 2:
            continue
        for f in group:
            p = _inner_point(f)
            if p is None:
                continue
            nf = f.normal.Normalized()
            best = None
            for g in group:
                if g is f or g.normal.Angle(f.normal) < math.pi - 1e-3:
                    continue
                d = gp_Vec(g.ad.Plane().Location(), p).Dot(nf)  # > 0: g lies behind f
                if d <= tol.linear or (best is not None and d >= best[0]):
                    continue
                foot = gp_Pnt(p.XYZ() - nf.XYZ() * d)
                if g.inside(*_plane_uv(g, foot)):
                    best = (d, g.idx)
            if best is not None:
                found[f.idx] = best
    if not found:
        return None, {}
    bins: dict[float, float] = defaultdict(float)
    for i, (d, _) in found.items():
        bins[round(d, 3)] += faces[i - 1].area
    t0 = max(bins, key=lambda k: bins[k])
    band = max(0.01, tol.thickness_rel * t0)
    return t0, {i: g for i, (d, g) in found.items() if abs(d - t0) <= band}


def _edge_line(e):
    c = BRepAdaptor_Curve(e)
    if c.GetType() != GeomAbs_Line:
        return None
    a, b = c.FirstParameter(), c.LastParameter()
    return c.Value(a), c.Value(b)


def _concave(e, f: _Face, h: _Face, scale: float) -> bool:
    """Is the edge between f and h concave (h rises on the outward side of f)?"""
    ends = _edge_line(e)
    p0, p1 = ends
    m = gp_Pnt((p0.X() + p1.X()) / 2, (p0.Y() + p1.Y()) / 2, (p0.Z() + p1.Z()) / 2)
    te = gp_Vec(p0, p1).Normalized()
    d = h.normal.Normalized().Crossed(te)
    q = gp_Pnt(m.XYZ() + d.XYZ() * scale)
    if not h.inside(*_plane_uv(h, q)):
        d.Reverse()
    return d.Dot(f.normal) > 0


def _sharp_bend_edges(faces, edges, partners, nb_all, smooth_nb, t, tol) -> dict[int, tuple]:
    """Sharp edges between skin faces that are bends (docs/design.md §5.3).

    Edge E between planes f, h is a sharp bend if both have partners at distance t and
    the partners meet either at a sharp edge (sharp on both sides) or through a
    cylinder of radius ~t (sharp inside, round outside).
    """
    out = {}
    for ei, (e, fi) in enumerate(edges):
        if len(fi) != 2:
            continue
        a, b = fi
        fa, fb = faces[a - 1], faces[b - 1]
        if fa.kind != "plane" or fb.kind != "plane" or b in smooth_nb[a]:
            continue
        if a not in partners or b not in partners or _edge_line(e) is None:
            continue
        pa, pb = partners[a], partners[b]
        if pa in (a, b) or pb in (a, b) or pa == pb:
            continue
        ok = pb in nb_all[pa]  # sharp on both sides
        if not ok:  # sharp inside, round (R ~ t) outside
            for z in smooth_nb[pa] & smooth_nb[pb]:
                fz = faces[z - 1]
                if fz.kind == "cylinder" and abs(fz.radius - t) <= max(0.01, tol.thickness_rel * t):
                    ok = True
                    break
        if ok:
            out[ei] = (a, b)
    return out


# ------------------------------------------------------------------------------ hole features


def _circle_radii(face: _Face) -> list[tuple[float, Any]]:
    from OCP.TopAbs import TopAbs_EDGE as _E
    from OCP.TopExp import TopExp_Explorer

    out = []
    ex = TopExp_Explorer(face.shape, _E)
    while ex.More():
        c = BRepAdaptor_Curve(as_edge(ex.Current()))
        if c.GetType() == GeomAbs_Circle:
            circ = c.Circle()
            out.append((circ.Radius(), circ.Location()))
        ex.Next()
    return out


def _hole_feature(group: list[int], faces: list[_Face], tol: Tolerances, t: float = math.inf):
    """A coaxial stack of cylinders, cones and annular planes -> (kind, through_r, outer_r, angle, depth)."""
    fs = [faces[i - 1] for i in group]
    axial = [f for f in fs if f.kind in ("cylinder", "cone")]
    if not axial or any(f.kind not in ("cylinder", "cone", "plane") for f in fs):
        return None
    key = _axis_key(axial[0])
    if any(
        _axis_key(f)[:3] != key[:3] or max(abs(x - y) for x, y in zip(_axis_key(f)[3:], key[3:], strict=True)) > 1e-2
        for f in axial
    ):
        return None
    ax = axial[0].axis_d.Normalized()
    for f in fs:
        if f.kind == "plane" and abs(abs(f.normal.Normalized().Dot(ax)) - 1) > 1e-6:
            return None
    cones = [f for f in axial if f.kind == "cone"]
    cyl_r = sorted({round(f.radius, 6) for f in axial if f.kind == "cylinder"})
    circles = [c for f in fs for c in _circle_radii(f)]
    if not circles:
        return None
    rmin = min(r for r, _ in circles)
    rmax = max(r for r, _ in circles)
    if cones:
        return "countersink", rmin, rmax, math.degrees(2 * abs(cones[0].semi_angle)), None
    if len(cyl_r) >= 2:
        big = [p for r, p in circles if abs(r - rmax) < 1e-6]
        depth = max((big[i].Distance(big[j]) for i in range(len(big)) for j in range(i + 1, len(big))), default=None)
        if depth is not None and depth > t * (1 + tol.thickness_rel):
            return None  # deeper than the sheet: a sharp-edged collar, not a counterbore
        return "counterbore", rmin, rmax, None, depth
    return None


# ------------------------------------------------------------------------------ formed features


@dataclass
class _Formed:
    """A region of one skin that is not laid flat by bending (docs/design.md §5.4).

    Formed regions (emboss, extruded hole, louvre, bridge) are left out of the
    development; a lance is a tab bent out of an opening, developed like any flange.
    """

    skin: int
    faces: set[int]
    host: int | None = None  # the plane face whose inner wire bounds the region
    kind: str | None = None  # "emboss" | "extruded_hole" | "louvre" | "bridge" | "lance" | None (not supported)
    reason: str = ""
    wire: list[int] = field(default_factory=list)  # edge indices of that inner wire
    cuts: list[list[int]] = field(default_factory=list)  # runs of the wire along cuts (slits), in wire order
    developed: bool = False  # a lance bent from planes and cylinders: laid flat like any flange
    gap: float = 0.0  # width of the slits, seen on the sheet
    free: list[int] = field(default_factory=list)  # the region's edges along thickness faces (its cut edge)


_FORM_KINDS = ("torus", "sphere", "cone", "other")


def _wires(face: _Face):
    """(outer wire, [inner wires]) of a face."""
    from OCP.TopAbs import TopAbs_WIRE
    from OCP.TopExp import TopExp_Explorer

    outer = BRepTools.OuterWire_s(face.shape)
    inner = []
    ex = TopExp_Explorer(face.shape, TopAbs_WIRE)
    while ex.More():
        if not ex.Current().IsSame(outer):
            inner.append(ex.Current())
        ex.Next()
    return outer, inner


def _wire_edges(wire, face: _Face, host: int, edges) -> list[int]:
    """Edge indices of a wire, in order along it."""
    from OCP.BRepTools import BRepTools_WireExplorer

    from ._occ import as_wire

    mine = [(ei, e) for ei, (e, fi) in enumerate(edges) if host in fi]
    out = []
    ex = BRepTools_WireExplorer(as_wire(wire), face.shape)
    while ex.More():
        cur = ex.Current()
        out += [ei for ei, e in mine if e.IsSame(cur)][:1]
        ex.Next()
    return out


def _runs(flags: list[bool]) -> list[list[int]]:
    """Maximal runs of True in a cyclic list, as lists of positions."""
    n = len(flags)
    if all(flags):
        return [list(range(n))]
    start = next(i for i in range(n) if not flags[i])
    runs: list[list[int]] = []
    cur: list[int] = []
    for k in range(1, n + 1):
        i = (start + k) % n
        if flags[i]:
            cur.append(i)
        elif cur:
            runs.append(cur)
            cur = []
    if cur:
        runs.append(cur)
    return runs


def _thin(chain: list[int], edges, t: float) -> bool:
    """Do a cut's edges stay within t of the chord between its ends (a straight slit or slot)?"""
    pts = []
    for ei in chain:
        c = BRepAdaptor_Curve(edges[ei][0])
        a, b = c.FirstParameter(), c.LastParameter()
        pts += [c.Value(a + (b - a) * k / 8) for k in range(9)]
    p, q = max(((p, q) for p in pts for q in pts), key=lambda pq: pq[0].Distance(pq[1]))
    if p.Distance(q) < 1e-9:
        return True
    d = gp_Vec(p, q).Normalized()
    return max(gp_Vec(p, x).Crossed(d).Magnitude() for x in pts) <= t * 1.01 + 1e-6


def _classify(reg: _Formed, faces: list[_Face], edges, nb_all, in_skin, t: float, wire, seeded: bool) -> None:
    assert reg.host is not None
    host = faces[reg.host - 1]
    reg.wire = _wire_edges(wire, host, reg.host, edges)
    other = [next(iter(set(edges[ei][1]) - {reg.host}), 0) for ei in reg.wire]
    attached = [o in reg.faces for o in other]
    cut = [o not in in_skin for o in other]
    if not all(a or c for a, c in zip(attached, cut, strict=True)):
        reg.reason = "a formed region meets other skin faces in its opening"
        return
    n_att = len(_runs(attached))
    reg.cuts = [[reg.wire[i] for i in run] for run in _runs(cut)] if any(cut) else []
    thin = [_thin(c, edges, t) for c in reg.cuts]
    reg.gap = _slit_width(reg, faces, edges, in_skin)
    if not seeded:  # planes and cylinders only
        if n_att == 1:
            reg.kind = "lance"  # a tab bent out of an opening: developed like any flange
            reg.developed = True
        elif n_att == 2 and len(reg.cuts) == 2:
            reg.kind = "bridge"
        else:
            reg.reason = "a region bent out of an opening and held on more than two sides is not supported"
        return
    if not reg.cuts:
        cut_inside = any(n not in in_skin for i in reg.faces for n in nb_all[i])
        if not cut_inside:
            reg.kind = "emboss"
        elif _collar(reg, faces) is not None:
            reg.kind = "extruded_hole"
        else:
            reg.reason = "a formed feature with cuts inside it is not supported"
    elif n_att == 1 and len(reg.cuts) == 1:
        # held on one side: a louvre is cut along a slit (straight or curved), a formed lance
        # (a tongue) is cut free around most of its rim
        cut_len = sum(_edge_length(edges[ei][0]) for ei in reg.cuts[0])
        rim_len = sum(_edge_length(edges[ei][0]) for ei in reg.wire)
        reg.kind = "louvre" if thin[0] or cut_len <= 0.5 * rim_len else "lance"
    elif n_att == 2 and len(reg.cuts) == 2:
        reg.kind = "bridge"
    else:
        reg.reason = "a formed feature held on more than two sides of its opening is not supported"


def _edge_length(e) -> float:
    g = GProp_GProps()
    BRepGProp.LinearProperties_s(e, g)
    return g.Mass()


def _slit_width(reg: _Formed, faces: list[_Face], edges, in_skin) -> float:
    """Width of the slits, seen on the sheet: the region's free edges projected onto the host
    plane, and their distance from the sheet's cut edges (the median, so the ends where the two
    meet do not count).
    """
    from OCP.BRep import BRep_Builder
    from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
    from OCP.BRepExtrema import BRepExtrema_DistShapeShape
    from OCP.TopoDS import TopoDS_Compound

    cut = [ei for c in reg.cuts for ei in c]
    free = [  # the region's edges along thickness faces
        ei
        for ei, (_, fi) in enumerate(edges)
        if len(fi) == 2 and len(set(fi) & reg.faces) == 1 and next(iter(set(fi) - reg.faces)) not in in_skin
    ]
    reg.free = free
    if not cut or not free or reg.host is None:
        return 0.0
    bb = BRep_Builder()
    comp = TopoDS_Compound()
    bb.MakeCompound(comp)
    for ei in cut:
        bb.Add(comp, edges[ei][0])
    h = faces[reg.host - 1]
    n = h.normal.Normalized()
    o = h.ad.Plane().Location()
    ds = []
    for ei in free:
        c = BRepAdaptor_Curve(edges[ei][0])
        a, b = c.FirstParameter(), c.LastParameter()
        for k in range(1, 8):
            p = c.Value(a + (b - a) * k / 8)
            q = gp_Pnt(p.XYZ() - n.XYZ() * gp_Vec(o, p).Dot(n))
            d = BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(q).Vertex(), comp)
            if d.IsDone():
                ds.append(d.Value())
    if not ds:
        return 0.0
    ds.sort()
    return ds[len(ds) // 2]


def _formed_regions(
    faces: list[_Face], skin: set[int], s: int, smooth_nb, nb_all, edges, in_skin, t: float
) -> list[_Formed]:
    """Formed regions and lances of one skin.

    Formed regions: seeds are the skin's tori, spheres, cones and freeform faces,
    grown over cylinders tangent to them. Removing the seeds leaves the main sheet (the
    largest part) and the faces the seeds enclose (an emboss's raised top). A region is
    valid only if it meets the main sheet along an inner wire of a single plane face.

    Without seeds, a region of planes and cylinders that hangs in an inner wire of a
    plane face (it touches nothing else) is a tab bent out of an opening (a lance) if it
    is held on one side, a bridge if on two.
    """
    out: list[_Formed] = []
    seeds = {i for i in skin if faces[i - 1].kind in _FORM_KINDS}
    seeded_faces: set[int] = set()
    if any(faces[i - 1].kind in ("torus", "sphere", "other") for i in seeds):  # cones alone: conical bends
        grow = True
        while grow:
            grow = False
            for i in skin:
                if i not in seeds and faces[i - 1].kind == "cylinder" and smooth_nb[i] & seeds:
                    seeds.add(i)
                    grow = True
        rest = skin - seeds
        rest_comps = [
            c
            for c in _components(len(faces), [(a, b) for a in rest for b in smooth_nb[a] if b in rest and a < b])
            if c[0] in rest
        ]
        main = set(max(rest_comps, key=lambda c: sum(faces[i - 1].area for i in c))) if rest_comps else set()
        formed = skin - main if main else set()  # nothing flat left: unsupported surfaces, not formed features
        comps = _components(len(faces), [(a, b) for a in formed for b in smooth_nb[a] if b in formed and a < b])
        for c in comps:
            if c[0] not in formed:
                continue
            reg = _Formed(s, set(c))
            out.append(reg)
            seeded_faces |= reg.faces
            hosts = {b for a in c for b in smooth_nb[a] if b in main}
            if len(hosts) != 1 or faces[next(iter(hosts)) - 1].kind != "plane":
                reg.reason = "a formed region must lie inside one flat face"
                continue
            reg.host = next(iter(hosts))
            outer, inner = _wires(faces[reg.host - 1])
            bound = [e for e, fi in edges if len(fi) == 2 and reg.host in fi and (set(fi) - {reg.host}) <= reg.faces]
            if any(_has_edge(outer, e) for e in bound):
                reg.reason = "a formed region reaches the outline of its face"
                continue
            wire = next((w for w in inner if any(_has_edge(w, e) for e in bound)), None)
            _classify(reg, faces, edges, nb_all, in_skin, t, wire, seeded=True)

    # regions of planes and cylinders hanging in an inner wire: lances and bridges
    largest = max((i for i in skin if faces[i - 1].kind == "plane"), key=lambda i: faces[i - 1].area, default=None)
    for h in sorted(skin - seeded_faces):
        fh = faces[h - 1]
        if fh.kind != "plane":
            continue
        _, inner = _wires(fh)
        for w in inner:
            start = {
                o
                for ei in _wire_edges(w, fh, h, edges)
                for o in edges[ei][1]
                if o != h and o in skin and o in smooth_nb[h] and o not in seeded_faces
            }
            if not start:
                continue
            region, todo = set(start), list(start)
            while todo:
                a = todo.pop()
                for b in smooth_nb[a]:
                    if b != h and b in skin and b not in region:
                        region.add(b)
                        todo.append(b)
            touches = {b for a in region for b in smooth_nb[a] if b not in region}
            if touches != {h} or largest in region or region & seeded_faces:
                continue
            reg = _Formed(s, region, host=h)
            _classify(reg, faces, edges, nb_all, in_skin, t, w, seeded=False)
            out.append(reg)
    return out


def _has_edge(wire, e) -> bool:
    from OCP.TopExp import TopExp_Explorer

    ex = TopExp_Explorer(wire, TopAbs_EDGE)
    while ex.More():
        if ex.Current().IsSame(e):
            return True
        ex.Next()
    return False


def _collar(reg: _Formed, faces: list[_Face]):
    """(cylinder, tori) if the region is a coaxial collar on an axis normal to its host."""
    fs = [faces[i - 1] for i in reg.faces]
    if any(f.kind not in ("cylinder", "torus") for f in fs) or reg.host is None:
        return None
    cyls = [f for f in fs if f.kind == "cylinder"]
    tori = [f for f in fs if f.kind == "torus"]
    if not cyls or len({round(f.radius, 6) for f in cyls}) != 1:
        return None
    n = faces[reg.host - 1].normal.Normalized()
    key = _axis_key(cyls[0])
    for f in fs:
        k = _axis_key(f)
        if k[:3] != key[:3] or max(abs(x - y) for x, y in zip(k[3:], key[3:], strict=True)) > 1e-2:
            return None
    if abs(abs(cyls[0].axis_d.Normalized().Dot(n)) - 1) > 1e-6:
        return None
    return cyls, tori


def _region_extent(reg: _Formed, faces: list[_Face]) -> tuple[float, float]:
    """Smallest and largest signed distance of the region from its host plane (along its normal)."""
    assert reg.host is not None
    h = faces[reg.host - 1]
    n = h.normal.Normalized()
    o = h.ad.Plane().Location()
    from OCP.TopExp import TopExp_Explorer

    from ._occ import as_vertex

    pts = []
    for i in reg.faces:
        f = faces[i - 1]
        ex = TopExp_Explorer(f.shape, TopAbs_VERTEX)
        while ex.More():
            pts.append(BRep_Tool.Pnt_s(as_vertex(ex.Current())))
            ex.Next()
        u0, u1, v0, v1 = f.uv
        for a in range(1, 8):
            for b in range(1, 8):
                u, v = u0 + (u1 - u0) * a / 8, v0 + (v1 - v0) * b / 8
                if f.inside(u, v):
                    pts.append(f.ad.Value(u, v))
    ds = [gp_Vec(o, p).Dot(n) for p in pts] or [0.0]
    return min(0.0, *ds), max(0.0, *ds)


def _host_circle(reg: _Formed, faces: list[_Face], edges) -> tuple[float, Any] | None:
    """Radius and centre of the region's boundary on its host, if it is a circle."""
    out = None
    for e, fi in edges:
        if len(fi) == 2 and reg.host in fi and (set(fi) - {reg.host}) <= reg.faces:
            c = BRepAdaptor_Curve(e)
            if c.GetType() != GeomAbs_Circle:
                return None
            circ = c.Circle()
            if out is not None and abs(out[0] - circ.Radius()) > 1e-6:
                return None
            out = (circ.Radius(), circ.Location())
    return out


def _formed_feature(fid: int, reg: _Formed, faces: list[_Face], edges, t: float, diags: list) -> Feature:
    """A Feature for a formed region of the reference skin."""
    lo, hi = _region_extent(reg, faces)
    up = hi >= -lo  # displaced towards the reference skin's outside: towards the viewer of the flat
    h = hi if up else -lo
    side: Literal["top", "bottom"] = "top" if up else "bottom"
    if reg.kind in ("louvre", "bridge", "lance"):
        return Feature(fid, reg.kind, 0.0, 0.0, side, face_ids=tuple(sorted(reg.faces)), height=h)  # type: ignore[arg-type]
    circ = _host_circle(reg, faces, edges)
    if reg.kind == "emboss":
        size = 2 * circ[0] if circ else 0.0
        return Feature(fid, "emboss", 0.0, size, side, face_ids=tuple(sorted(reg.faces)), height=h)
    cyls, tori = _collar(reg, faces)  # type: ignore[misc]
    if circ is None:
        raise refuse(Code.FORMED_FEATURE, "an extruded hole whose rim is not a circle", faces=sorted(reg.faces))
    rho = circ[0]
    r_cyl = cyls[0].radius
    straight = max(f.uv[3] - f.uv[2] for f in cyls)  # v runs along the axis
    if tori:
        tor = tori[0]
        sweep = max(f.uv[3] - f.uv[2] for f in tori)
        r_in = tor.minor if up else tor.minor - t  # the reference skin is the inside of the fillet if up
    else:
        sweep, r_in = 0.0, 0.0
    bore = r_cyl - t if up else r_cyl
    height = h + t if up else h  # outside height, from the far side of the sheet to the rim
    d_pre = extruded_hole_prehole(rho, r_in, sweep, straight, t)
    if d_pre <= 0:
        diags.append(warn(Code.PRE_HOLE, f"extruded hole ⌀{2 * bore:g}: no material for a pre-hole; pierced",
                          faces=sorted(reg.faces)))  # fmt: skip
    return Feature(
        fid,
        "extruded_hole",
        max(0.0, d_pre),
        2 * bore,
        side,
        face_ids=tuple(sorted(reg.faces)),
        height=height,
        radius=r_in,
    )


# ------------------------------------------------------------------------------ engravings


@dataclass
class _Recess:
    """A recess that does not go through the sheet: a pocket, a groove, engraved text, a blind hole."""

    faces: list[int]
    host: int  # the skin plane it opens into
    rims: list[tuple[int, list[int]]]  # (the face it bounds, its edge indices in wire order) per rim
    depth: float


def _recess(group: list[int], skin: set[int], faces: list[_Face], edges, t: float, tol: Tolerances):
    """The group as a recess into one plane of ``skin``, or None.

    It must open into a single plane face of the skin only through inner wires of that face
    (so it reaches no edge of the part), and stay less than the thickness below the plane.
    Faces of the group flush with that plane are islands its floor leaves standing (the
    middle of an O): their rims are rims of the recess too. Every edge of a rim must be the
    recess's own (none is shared with a through cut).
    """
    gs = set(group)
    hosts = {j for _, fi in edges if len(fi) == 2 and (fi[0] in gs) != (fi[1] in gs) for j in fi if j not in gs}
    if len(hosts) != 1:
        return None
    h = next(iter(hosts))
    host = faces[h - 1]
    if h not in skin or host.kind != "plane":
        return None
    n = host.normal.Normalized()
    o = host.ad.Plane().Location()
    islands = {
        i
        for i in gs
        if faces[i - 1].kind == "plane"
        and faces[i - 1].normal.Angle(host.normal) < tol.smooth_angle
        and abs(gp_Vec(o, faces[i - 1].ad.Plane().Location()).Dot(n)) <= tol.linear
    }
    walls = gs - islands
    if not walls:
        return None
    mouth = {ei for ei, (_, fi) in enumerate(edges) if len(fi) == 2 and (fi[0] in walls) != (fi[1] in walls)}
    outer, inner = _wires(host)
    if any(_has_edge(outer, edges[ei][0]) for ei in mouth if h in edges[ei][1]):
        return None
    rims: list[tuple[int, list[int]]] = []
    surfaces = [(h, inner)]
    for i in sorted(islands):
        o_i, inner_i = _wires(faces[i - 1])
        surfaces.append((i, [o_i, *inner_i]))
    for i, wires in surfaces:
        for w in wires:
            if any(_has_edge(w, edges[ei][0]) for ei in mouth if i in edges[ei][1]):
                rims.append((i, _wire_edges(w, faces[i - 1], i, edges)))
    if not rims or {ei for _, rim in rims for ei in rim} != mouth:
        return None
    lo, hi = _region_extent(_Formed(-1, walls, host=h), faces)
    depth = -lo
    if hi > tol.linear or depth <= tol.linear or depth >= t * (1 - tol.thickness_rel):
        return None
    return _Recess(sorted(group), h, rims, depth)


def _rim_diameter(rim: list[int], edges) -> float:
    """The diameter of a rim that is one full circle, else 0."""
    if len(rim) != 1:
        return 0.0
    c = BRepAdaptor_Curve(edges[rim[0]][0])
    if c.GetType() != GeomAbs_Circle or abs(c.LastParameter() - c.FirstParameter() - 2 * math.pi) > 1e-9:
        return 0.0
    return 2 * c.Circle().Radius()


# ------------------------------------------------------------------------------ public


def recognize(
    model: StepModel | Any,
    *,
    thickness: float | None = None,
    solid: int | None = None,
    tol: Tolerances | None = None,
    canonical: bool = True,
    engrave: bool = False,
) -> SheetPart:
    """Recognise a solid as sheet metal: thickness, skins, bends, flanges, hole features.

    ``model`` is a :class:`~flatmat.io.StepModel` or an OCP solid. With ``canonical``
    (default), B-spline faces and edges that are really planes, cylinders, cones, lines
    or circles are converted first (on a copy).

    With ``engrave``, a recess that does not go through the sheet (a pocket, a groove,
    engraved text, a blind hole) is an ``engraving`` feature: its rim is marked on the flat
    pattern and its depth is not made. Without it such a part is refused (``EDGE_FEATURE``).
    """
    tol = tol or Tolerances()
    diags: list[Diagnostic] = []
    source = ""
    if isinstance(model, StepModel):
        source = model.path
        diags.extend(model.diagnostics)
        shape = model.solids[solid] if solid is not None else model.solid
    else:
        shape = model

    if canonical:
        shape, converted = canonicalize(shape, tol.linear)
        if converted is not None:
            diags.append(converted)

    fmap, faces = _face_table(shape)
    edges = _edge_table(shape, fmap)
    vmap = ShapeMap()
    TopExp.MapShapes_s(shape, TopAbs_VERTEX, vmap)

    smooth = []
    smooth_nb: dict[int, set[int]] = defaultdict(set)
    nb_all: dict[int, set[int]] = defaultdict(set)
    sharp_nb: dict[int, set[int]] = defaultdict(set)
    for ei, (e, fi) in enumerate(edges):
        if BRep_Tool.Degenerated_s(e):
            continue
        if len(fi) == 1:
            if not BRep_Tool.IsClosed_s(e, faces[fi[0] - 1].shape):
                raise refuse(Code.OPEN_SHELL, "an edge borders only one face", faces=fi, edges=(ei,))
            continue
        if len(fi) > 2:
            raise refuse(Code.NON_MANIFOLD, "an edge borders more than two faces", faces=fi, edges=(ei,))
        a, b = fi
        nb_all[a].add(b)
        nb_all[b].add(a)
        if _is_smooth(e, faces[a - 1], faces[b - 1], tol):
            smooth.append((a, b, ei))
            smooth_nb[a].add(b)
            smooth_nb[b].add(a)
        else:
            sharp_nb[a].add(b)
            sharp_nb[b].add(a)

    # sharp bends join faces of one skin across a non-tangent edge
    t_est, partners = _plane_partners(faces, tol)
    sharp = _sharp_bend_edges(faces, edges, partners, nb_all, smooth_nb, t_est, tol) if t_est else {}
    links = [(a, b) for a, b, _ in smooth] + list(sharp.values())

    comps = sorted(_components(len(faces), links), key=lambda c: -sum(faces[i - 1].area for i in c))
    if len(comps) < 3:
        raise refuse(Code.NOT_TWO_SKINS, f"{len(comps)} face groups; a sheet has two skins and edge faces")
    skins = [set(comps[0]), set(comps[1])]
    in_skin = skins[0] | skins[1]
    regions: list[_Formed] = []
    t_guess = thickness or t_est or 1.0  # to tell a slit from an opening; refined below
    for s, sk in enumerate(skins):
        regions += _formed_regions(faces, sk, s, smooth_nb, nb_all, edges, in_skin, t_guess)
    for reg in regions:
        if reg.kind is None:
            raise refuse(Code.FORMED_FEATURE, reg.reason, faces=sorted(reg.faces))
    formed = {i for reg in regions if not reg.developed for i in reg.faces}
    for sk in skins:
        bad = [i for i in sk if faces[i - 1].kind in _FORM_KINDS and i not in formed]
        if bad:
            raise refuse(Code.UNSUPPORTED_SURFACE, "skin faces must be planes or cylinders", faces=bad)

    t, conf = _thickness(faces, skins, tol)
    if thickness is not None:
        if abs(thickness - t) > max(0.01, tol.thickness_rel * thickness):
            raise refuse(Code.VARIABLE_THICKNESS, f"given thickness {thickness} but measured {t:.4f}", measured=t)
        t = thickness
    if conf < tol.thickness_confidence:
        raise refuse(
            Code.VARIABLE_THICKNESS, f"only {conf:.0%} of the sheet agrees on t = {t:.4f}", thickness=t, confidence=conf
        )

    # faces between the skins: hole features (coaxial stacks) or edge faces touching both skins
    between = [i for i in range(1, len(faces) + 1) if i not in in_skin]
    groups = _components(len(faces), [(a, b) for a in between for b in nb_all[a] if b not in in_skin and a < b])
    groups = [g for g in groups if g[0] not in in_skin]
    feature_groups: list[tuple[list[int], tuple]] = []
    recesses: list[_Recess] = []
    featured: set[int] = set()
    for g in groups:
        touches = {s for s in (0, 1) if any(n in skins[s] for i in g for n in nb_all[i])}
        hf = _hole_feature(g, faces, tol, t) if touches == {0, 1} else None
        if hf is not None:
            feature_groups.append((g, hf))
            featured |= set(g)
        elif engrave and len(touches) == 1:
            rec = _recess(g, skins[next(iter(touches))], faces, edges, t, tol)
            if rec is not None:
                recesses.append(rec)
                featured |= set(g)
    for comp in comps[2:]:
        if set(comp) & featured:
            continue
        touches = {s for s in (0, 1) if any(n in skins[s] for i in comp for n in sharp_nb[i])}
        if touches != {0, 1}:
            raise refuse(
                Code.EDGE_FEATURE,
                "a face group between the skins touches only one skin (chamfer, formed feature?)",
                faces=comp,
            )

    planes = [f for f in faces if f.kind == "plane" and f.idx in in_skin and f.idx not in formed]
    if not planes:
        raise refuse(Code.UNSUPPORTED_SURFACE, "no planar skin face to start from")
    root = max(planes, key=lambda f: (f.area, -f.idx))
    a_side = 0 if root.idx in skins[0] else 1
    skin_a = skins[a_side] - formed

    # bends: cylinders of the reference skin, grouped by axis line, radius and angle
    groups_b: dict[tuple, list] = defaultdict(list)
    bend_faces: dict[int, _BendFace] = {}
    for i in sorted(skin_a):
        f = faces[i - 1]
        if f.kind != "cylinder":
            continue
        u0, u1, v0, v1 = f.uv
        p = f.ad.Value((u0 + u1) / 2, (v0 + v1) / 2)
        ap = gp_Vec(f.axis_p, p)
        radial = ap - f.axis_d.Multiplied(ap.Dot(f.axis_d) / f.axis_d.SquareMagnitude())
        concave = f.normal.Dot(radial) < 0  # normal points at the axis: A is the inside
        r_in = f.radius if concave else f.radius - t
        if abs(r_in) < tol.linear * 10:
            r_in = 0.0
        ang = u1 - u0
        if math.degrees(ang) >= 359.9:
            raise refuse(Code.BEND_TOO_LARGE, f"bend of {math.degrees(ang):.1f}°", faces=(i,))
        bf = _BendFace(i, 0, r_in, concave, ang)
        key = ("cyl", *_axis_key(f), round(r_in, 3), round(math.degrees(ang), 2), concave)
        groups_b[key].append(bf)
        bend_faces[i] = bf

    # sharp bends on the reference skin, grouped by edge line and angle
    sharp_edges: dict[int, _SharpEdge] = {}
    for ei, (a, b) in sharp.items():
        if a not in skin_a:
            continue
        e = edges[ei][0]
        p0, p1 = _edge_line(e)
        fa, fb = faces[a - 1], faces[b - 1]
        ang = fa.normal.Angle(fb.normal)
        concave = _concave(e, fa, fb, max(1e-3, 1e-3 * t))
        vs = _edge_vertex_ids(e, vmap)
        se = _SharpEdge(ei, (a, b), ang, concave, p0.Distance(p1), vs)
        key = ("sharp", *_line_key(p0, gp_Vec(p0, p1)), round(math.degrees(ang), 2), concave)
        groups_b[key].append(se)
        sharp_edges[ei] = se

    # flanges: components of skin A once the bend faces are removed (sharp edges separate too)
    flat_ids = {i for i in skin_a if i not in bend_faces}
    flange_links = [(a, b) for a, b, _ in smooth if a in flat_ids and b in flat_ids]
    flange_comps = sorted(sorted(c) for c in _components(len(faces), flange_links) if c[0] in flat_ids)
    flange_of = {i: n + 1 for n, c in enumerate(flange_comps) for i in c}
    flanges = tuple(Flange(n + 1, tuple(c), sum(faces[i - 1].area for i in c)) for n, c in enumerate(flange_comps))

    def first_id(members) -> int:
        return min(m.face if isinstance(m, _BendFace) else 100000 + m.edge for m in members)

    bends = []
    for n, (key, members) in enumerate(sorted(groups_b.items(), key=lambda kv: first_id(kv[1]))):
        bid = n + 1
        segs = []
        adj_flanges: set[int] = set()
        if key[0] == "cyl":
            for m in sorted(members, key=lambda m: m.face):
                m.bend_id = bid
                f = faces[m.face - 1]
                u0, u1, v0, v1 = f.uv
                segs.append(
                    BendSegment(
                        m.face,
                        v1 - v0,
                        (f.axis_p.X(), f.axis_p.Y(), f.axis_p.Z()),
                        (f.axis_d.X(), f.axis_d.Y(), f.axis_d.Z()),
                    )
                )
                adj_flanges |= {flange_of[nb] for nb in smooth_nb[m.face] if nb in flange_of}
            m0 = members[0]
            deg = math.degrees(m0.angle)
            kind: Literal["bend", "hem", "sharp"] = "sharp" if m0.r_in == 0.0 else ("hem" if deg >= 179.0 else "bend")
            bends.append(Bend(bid, tuple(segs), m0.r_in, deg, "up" if m0.concave else "down", kind,
                              tuple(sorted(adj_flanges))))  # fmt: skip
        else:
            for m in sorted(members, key=lambda m: m.edge):
                m.bend_id = bid
                p0, p1 = _edge_line(edges[m.edge][0])
                segs.append(
                    BendSegment(
                        0, m.length, (p0.X(), p0.Y(), p0.Z()), tuple(gp_Vec(p0, p1).Normalized().Coord()), m.edge
                    )  # fmt: skip
                )
                adj_flanges |= {flange_of[x] for x in m.faces if x in flange_of}
            m0 = members[0]
            bends.append(Bend(bid, tuple(segs), 0.0, math.degrees(m0.angle), "up" if m0.concave else "down",
                              "sharp", tuple(sorted(adj_flanges))))  # fmt: skip

    features = []
    feat_faces: dict[int, tuple[int, ...]] = {}
    feat_r: dict[int, float] = {}
    feat_a: dict[int, bool] = {}
    for n, (g, (fkind, rmin, rmax, angle, depth)) in enumerate(feature_groups):
        fid = n + 1
        # the side carrying the large rim: the skin touched by the face holding the largest circle
        big_face = max(g, key=lambda i: max((r for r, _ in _circle_radii(faces[i - 1])), default=0))
        on_a = any(nb in skin_a for nb in nb_all[big_face])
        features.append(Feature(fid, fkind, 2 * rmin, 2 * rmax, "top" if on_a else "bottom", angle, depth, tuple(g)))
        feat_faces[fid] = tuple(g)
        feat_r[fid] = rmin
        feat_a[fid] = on_a
        if fkind == "countersink" and on_a:
            diags.append(warn(Code.COUNTERSINK_PROJECTED, f"countersink ⌀{2 * rmax:g} cut as ⌀{2 * rmin:g}", faces=g))

    feat_edges: dict[int, set[int]] = {}
    feat_cuts: dict[int, list[list[int]]] = {}
    feat_gap: dict[int, float] = {}
    feat_free: dict[int, tuple[int, list[int]]] = {}
    developed: set[int] = set()
    for reg in (r for r in regions if r.skin == a_side):
        fid = len(features) + 1
        feat = _formed_feature(fid, reg, faces, edges, t, diags)
        feat_gap[fid] = reg.gap
        feat_free[fid] = (reg.host, reg.free)  # type: ignore[assignment]
        if reg.developed:
            hinge = next((b for b in bends for sg in b.segments if sg.face_id in reg.faces
                          and reg.host in smooth_nb[sg.face_id]), None)  # fmt: skip
            feat = replace(feat, bend_id=hinge.id if hinge else None)
            developed.add(fid)
        else:
            feat_faces[fid] = tuple(sorted(reg.faces))
            feat_r[fid] = feat.through_diameter / 2
            feat_a[fid] = True
        features.append(feat)
        feat_edges[fid] = set(reg.wire)
        feat_cuts[fid] = reg.cuts

    # engravings: one per rim (a ring-shaped groove has two). A rim in a face of the reference
    # skin is chained with the outline; any other is projected onto the plane of that skin behind it.
    projected_rims: dict[int, tuple[int, list[int]]] = {}
    for rec in recesses:
        on_a = rec.host in skin_a
        behind = rec.host if on_a else partners.get(rec.host)
        if behind is None or behind not in skin_a:
            raise refuse(Code.EDGE_FEATURE, "a recess in a face with no face of the other skin behind it",
                         faces=rec.faces)  # fmt: skip
        for face, rim in rec.rims:
            fid = len(features) + 1
            features.append(Feature(fid, "engraving", 0.0, _rim_diameter(rim, edges), "top" if on_a else "bottom",
                                    depth=rec.depth, face_ids=tuple(rec.faces)))  # fmt: skip
            if face in skin_a:
                feat_edges[fid] = set(rim)
            else:
                projected_rims[fid] = (behind, rim)
    if recesses:
        deepest = max(r.depth for r in recesses)
        diags.append(warn(Code.ENGRAVED_RECESS,
                          f"{len(recesses)} recess(es) up to {deepest:.3g} deep engraved along their rims; "
                          "their depth is not made",
                          faces=sorted(i for r in recesses for i in r.faces),
                          recesses=len(recesses), depth=deepest))  # fmt: skip

    adj: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for a, b, ei in smooth:
        adj[a].append((b, ei))
        adj[b].append((a, ei))
    for ei, (a, b) in sharp.items():
        if a in skin_a:
            adj[a].append((b, ei))
            adj[b].append((a, ei))

    role = {i: "skin_a" for i in skin_a} | {i: "skin_b" for i in skins[1 - a_side]} | {i: "formed" for i in formed}
    face_infos = tuple(
        FaceInfo(
            f.idx,
            cast(Literal["plane", "cylinder", "other"], f.kind if f.kind in ("plane", "cylinder") else "other"),
            f.area,
            cast(Literal["skin_a", "skin_b", "thickness", "formed"], role.get(f.idx, "thickness")),
        )
        for f in faces
    )
    ctx = Context(shape, fmap, vmap, faces, edges, adj, skin_a, root.idx, t, bend_faces, tol=tol)
    ctx.sharp_edges = sharp_edges
    ctx.features, ctx.feature_through_r, ctx.feature_on_a = feat_faces, feat_r, feat_a
    ctx.feature_kind = {f.id: f.kind for f in features}
    ctx.feature_edges, ctx.feature_cuts = feat_edges, feat_cuts
    ctx.feature_gap, ctx.developed, ctx.feature_free = feat_gap, developed, feat_free
    ctx.projected_rims = projected_rims
    return SheetPart(
        thickness=t,
        thickness_confidence=conf,
        faces=face_infos,
        bends=tuple(bends),
        flanges=flanges,
        reference_root=root.idx,
        diagnostics=tuple(diags),
        source=source,
        features=tuple(features),
        _ctx=ctx,
    )


def _edge_vertex_ids(e, vmap) -> tuple[int, int]:
    from OCP.TopExp import TopExp_Explorer

    from ._occ import as_vertex

    vs = []
    ex = TopExp_Explorer(e, TopAbs_VERTEX)
    while ex.More():
        vs.append(vmap.FindIndex(as_vertex(ex.Current())))
        ex.Next()
    return vs[0], vs[-1]
