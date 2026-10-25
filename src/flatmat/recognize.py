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

from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepGProp import BRepGProp
from OCP.BRepTools import BRepTools, BRepTools_WireExplorer
from OCP.GeomAbs import (
    GeomAbs_Circle,
    GeomAbs_Cone,
    GeomAbs_Cylinder,
    GeomAbs_Line,
    GeomAbs_Plane,
    GeomAbs_Sphere,
    GeomAbs_Torus,
)
from OCP.gp import gp_Pnt, gp_Vec
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX, TopAbs_WIRE
from OCP.TopExp import TopExp, TopExp_Explorer
from OCP.TopoDS import TopoDS_Compound

from ._occ import EdgeFaceMap, ShapeMap, as_edge, as_face, as_vertex, as_wire
from ._topology import BendFace, Context, Face, SharpEdge, edge_vertices, inner_uv, pcurve, surface_props
from .bendmodel import extruded_hole_prehole
from .canonical import canonicalize
from .diagnostics import Code, Diagnostic, refuse, warn
from .io import StepModel
from .model import Bend, BendSegment, FaceInfo, Feature, Flange, SheetPart, Tolerances

# ------------------------------------------------------------------------------ tables


def _face_table(solid) -> tuple[Any, list[Face]]:
    fmap = ShapeMap()
    TopExp.MapShapes_s(solid, TopAbs_FACE, fmap)
    faces = []
    for i in range(1, fmap.Extent() + 1):
        f = as_face(fmap.FindKey(i))
        ad = BRepAdaptor_Surface(f, True)
        u0, u1, v0, v1 = BRepTools.UVBounds_s(f)
        n, mirror = surface_props(f, ad, (u0 + u1) / 2, (v0 + v1) / 2)
        g = GProp_GProps()
        BRepGProp.SurfaceProperties_s(f, g)
        kind = {
            GeomAbs_Plane: "plane",
            GeomAbs_Cylinder: "cylinder",
            GeomAbs_Cone: "cone",
            GeomAbs_Torus: "torus",
            GeomAbs_Sphere: "sphere",
        }.get(ad.GetType(), "other")
        rec = Face(i, f, kind, g.Mass(), n, mirror, (u0, u1, v0, v1), ad)
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


# ------------------------------------------------------------------------------ graph


@dataclass
class _Graph:
    """Face adjacency over the manifold edges, split by how the two faces meet."""

    smooth: list[tuple[int, int, int]]  # (face, face, edge) where the faces meet tangentially
    smooth_nb: dict[int, set[int]]
    sharp_nb: dict[int, set[int]]
    nb_all: dict[int, set[int]]


def _is_smooth(e, f1: Face, f2: Face, tol: Tolerances) -> bool:
    c = BRepAdaptor_Curve(e)
    a, b = c.FirstParameter(), c.LastParameter()
    worst = 0.0
    for s in (0.25, 0.5, 0.75):
        tm = a + s * (b - a)
        n = []
        for f in (f1, f2):
            uv = pcurve(e, f.shape).Value(tm)
            n.append(surface_props(f.shape, f.ad, uv.X(), uv.Y())[0])
        worst = max(worst, n[0].Angle(n[1]))
    return worst < tol.smooth_angle


def _graph(faces: list[Face], edges, tol: Tolerances) -> _Graph:
    """The adjacency graph; refuses open shells and non-manifold edges."""
    g = _Graph([], defaultdict(set), defaultdict(set), defaultdict(set))
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
        g.nb_all[a].add(b)
        g.nb_all[b].add(a)
        if _is_smooth(e, faces[a - 1], faces[b - 1], tol):
            g.smooth.append((a, b, ei))
            g.smooth_nb[a].add(b)
            g.smooth_nb[b].add(a)
        else:
            g.sharp_nb[a].add(b)
            g.sharp_nb[b].add(a)
    return g


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


def _skins(faces: list[Face], links) -> tuple[list[list[int]], list[set[int]]]:
    """Face groups by area (largest first) and the two skins."""
    comps = sorted(_components(len(faces), links), key=lambda c: -sum(faces[i - 1].area for i in c))
    if len(comps) < 3:
        raise refuse(Code.NOT_TWO_SKINS, f"{len(comps)} face groups; a sheet has two skins and edge faces")
    return comps, [set(comps[0]), set(comps[1])]


def _check_skin_surfaces(faces: list[Face], skins: list[set[int]], formed: set[int]) -> None:
    """Outside formed features, the skins must be planes or cylinders."""
    for sk in skins:
        bad = [i for i in sk if faces[i - 1].kind in _FORM_KINDS and i not in formed]
        if bad:
            raise refuse(Code.UNSUPPORTED_SURFACE, "skin faces must be planes or cylinders", faces=bad)


# ------------------------------------------------------------------------------ geometry keys


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


def _axis_key(f: Face) -> tuple[float, ...]:
    return _line_key(f.axis_p, f.axis_d)


def _plane_uv(f: Face, p: gp_Pnt) -> tuple[float, float]:
    pos = f.ad.Plane().Position()
    d = gp_Vec(pos.Location(), p)
    return d.Dot(gp_Vec(pos.XDirection())), d.Dot(gp_Vec(pos.YDirection()))


def _inner_point(f: Face):
    uv = inner_uv(f)
    return None if uv is None else f.ad.Value(*uv)


# ------------------------------------------------------------------------------ thickness


def _mode(cands: list[tuple[float, float]], tol: Tolerances) -> tuple[float, float]:
    """The (value, weight) bin with the most weight, and the band around it that counts as equal."""
    bins: dict[float, float] = defaultdict(float)
    for v, w in cands:
        bins[round(v, 3)] += w
    t0 = max(bins, key=lambda k: bins[k])
    return t0, max(0.01, tol.thickness_rel * t0)


def _plane_partners(faces: list[Face], tol: Tolerances) -> dict[int, tuple[float, int]]:
    """For every plane, the opposite plane: plane -> (distance, partner).

    Planes are bucketed by their (unsigned) normal direction; the partner of f is the
    nearest antiparallel plane *behind* f (into the material) whose face contains the
    foot of an inner point of f. Skins are not known yet: sharp bends split them.
    """
    buckets: dict[tuple, list[Face]] = defaultdict(list)
    for f in faces:
        if f.kind == "plane":
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
    return found


def _thickness(
    faces: list[Face], skins: list[set[int]], partners: dict[int, tuple[float, int]], tol: Tolerances
) -> tuple[float, float]:
    """Thickness and the share of the sheet (by area) that agrees with it.

    Candidates: coaxial cylinders on opposite skins, and planes whose partner lies on the
    other skin.
    """
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
    for i, (d, g) in partners.items():
        if (i in skins[0] and g in skins[1]) or (i in skins[1] and g in skins[0]):
            cands.append((d, faces[i - 1].area))
    if not cands:
        raise refuse(Code.NO_THICKNESS, "no pair of opposite faces found")
    t0, band = _mode(cands, tol)
    near = [(v, w) for v, w in cands if abs(v - t0) <= band]
    t = sum(v * w for v, w in near) / sum(w for _, w in near)
    conf = sum(w for _, w in near) / sum(w for _, w in cands)
    return t, conf


def _checked_thickness(t: float, conf: float, given: float | None, tol: Tolerances) -> float:
    if given is not None:
        if abs(given - t) > max(0.01, tol.thickness_rel * given):
            raise refuse(Code.VARIABLE_THICKNESS, f"given thickness {given} but measured {t:.4f}", measured=t)
        t = given
    if conf < tol.thickness_confidence:
        raise refuse(
            Code.VARIABLE_THICKNESS, f"only {conf:.0%} of the sheet agrees on t = {t:.4f}", thickness=t, confidence=conf
        )
    return t


# ------------------------------------------------------------------------------ sharp bends


def _edge_line(e):
    c = BRepAdaptor_Curve(e)
    if c.GetType() != GeomAbs_Line:
        return None
    a, b = c.FirstParameter(), c.LastParameter()
    return c.Value(a), c.Value(b)


def _concave(e, f: Face, h: Face, scale: float) -> bool:
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


def _sharp_bend_edges(faces, edges, partners, g: _Graph, t, tol) -> dict[int, tuple[int, int]]:
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
        if fa.kind != "plane" or fb.kind != "plane" or b in g.smooth_nb[a]:
            continue
        if a not in partners or b not in partners or _edge_line(e) is None:
            continue
        pa, pb = partners[a], partners[b]
        if pa in (a, b) or pb in (a, b) or pa == pb:
            continue
        ok = pb in g.nb_all[pa]  # sharp on both sides
        if not ok:  # sharp inside, round (R ~ t) outside
            for z in g.smooth_nb[pa] & g.smooth_nb[pb]:
                fz = faces[z - 1]
                if fz.kind == "cylinder" and abs(fz.radius - t) <= max(0.01, tol.thickness_rel * t):
                    ok = True
                    break
        if ok:
            out[ei] = (a, b)
    return out


# ------------------------------------------------------------------------------ hole features


def _circle_radii(face: Face) -> list[tuple[float, Any]]:
    out = []
    ex = TopExp_Explorer(face.shape, TopAbs_EDGE)
    while ex.More():
        c = BRepAdaptor_Curve(as_edge(ex.Current()))
        if c.GetType() == GeomAbs_Circle:
            circ = c.Circle()
            out.append((circ.Radius(), circ.Location()))
        ex.Next()
    return out


def _hole_feature(group: list[int], faces: list[Face], tol: Tolerances, t: float = math.inf):
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


def _between_skins(
    faces: list[Face], edges, comps, skins: list[set[int]], g: _Graph, t: float, engrave: bool, tol: Tolerances
):
    """Hole features among the faces between the skins, [(faces, feature tuple)], and with
    ``engrave`` the recesses.

    The other face groups between the skins must touch both (edge faces); one that
    touches only one is refused.
    """
    in_skin = skins[0] | skins[1]
    between = [i for i in range(1, len(faces) + 1) if i not in in_skin]
    groups = _components(len(faces), [(a, b) for a in between for b in g.nb_all[a] if b not in in_skin and a < b])
    groups = [grp for grp in groups if grp[0] not in in_skin]
    found: list[tuple[list[int], tuple]] = []
    recesses: list[_Recess] = []
    featured: set[int] = set()
    for grp in groups:
        touches = {s for s in (0, 1) if any(n in skins[s] for i in grp for n in g.nb_all[i])}
        hf = _hole_feature(grp, faces, tol, t) if touches == {0, 1} else None
        if hf is not None:
            found.append((grp, hf))
            featured |= set(grp)
        elif engrave and len(touches) == 1:
            rec = _recess(grp, skins[next(iter(touches))], faces, edges, t, tol)
            if rec is not None:
                recesses.append(rec)
                featured |= set(grp)
    for comp in comps[2:]:
        if set(comp) & featured:
            continue
        touches = {s for s in (0, 1) if any(n in skins[s] for i in comp for n in g.sharp_nb[i])}
        if touches != {0, 1}:
            raise refuse(
                Code.EDGE_FEATURE,
                "a face group between the skins touches only one skin (chamfer, formed feature?)",
                faces=comp,
            )
    return found, recesses


@dataclass
class _Features:
    """Feature records as they are found, and the per-feature data unfold() needs."""

    items: list[Feature] = field(default_factory=list)
    faces: dict[int, tuple[int, ...]] = field(default_factory=dict)  # faces the outline cuts as a hole
    through_r: dict[int, float] = field(default_factory=dict)
    on_a: dict[int, bool] = field(default_factory=dict)
    edges: dict[int, set[int]] = field(default_factory=dict)
    cuts: dict[int, list[list[int]]] = field(default_factory=dict)
    gap: dict[int, float] = field(default_factory=dict)
    free: dict[int, tuple[int, list[int]]] = field(default_factory=dict)
    developed: set[int] = field(default_factory=set)
    projected_rims: dict[int, tuple[int, list[int]]] = field(default_factory=dict)
    diags: list[Diagnostic] = field(default_factory=list)

    @property
    def next_id(self) -> int:
        return len(self.items) + 1


def _hole_features(table: _Features, feature_groups, faces: list[Face], g: _Graph, skin_a: set[int]) -> None:
    """Countersinks and counterbores, with a warning per projected countersink."""
    for grp, (fkind, rmin, rmax, angle, depth) in feature_groups:
        fid = table.next_id
        # the side carrying the large rim: the skin touched by the face holding the largest circle
        big_face = max(grp, key=lambda i: max((r for r, _ in _circle_radii(faces[i - 1])), default=0))
        on_a = any(nb in skin_a for nb in g.nb_all[big_face])
        table.items.append(
            Feature(fid, fkind, 2 * rmin, 2 * rmax, "top" if on_a else "bottom", angle, depth, tuple(grp))
        )
        table.faces[fid] = tuple(grp)
        table.through_r[fid] = rmin
        table.on_a[fid] = on_a
        if fkind == "countersink" and on_a:
            table.diags.append(
                warn(Code.COUNTERSINK_PROJECTED, f"countersink ⌀{2 * rmax:g} cut as ⌀{2 * rmin:g}", faces=grp)
            )


def _formed_features(
    table: _Features, regions: list[_Formed], a_side: int, faces: list[Face], edges, t: float, bends, g: _Graph
) -> None:
    """The formed regions of the reference skin: left out of the development, or (lances) developed."""
    for reg in (r for r in regions if r.skin == a_side):
        fid = table.next_id
        feat = _formed_feature(fid, reg, faces, edges, t, table.diags)
        table.gap[fid] = reg.gap
        table.free[fid] = (reg.host, reg.free)  # type: ignore[assignment]
        if reg.developed:
            hinge = next((b for b in bends for sg in b.segments if sg.face_id in reg.faces
                          and reg.host in g.smooth_nb[sg.face_id]), None)  # fmt: skip
            feat = replace(feat, bend_id=hinge.id if hinge else None)
            table.developed.add(fid)
        else:
            table.faces[fid] = tuple(sorted(reg.faces))
            table.through_r[fid] = feat.through_diameter / 2
            table.on_a[fid] = True
        table.items.append(feat)
        table.edges[fid] = set(reg.wire)
        table.cuts[fid] = reg.cuts


def _engravings(table: _Features, recesses: list[_Recess], skin_a: set[int], partners: dict[int, int], edges) -> None:
    """One engraving per rim (a ring-shaped groove has two).

    A rim in a face of the reference skin is chained with the outline; any other is
    projected onto the plane of that skin behind it.
    """
    for rec in recesses:
        on_a = rec.host in skin_a
        behind = rec.host if on_a else partners.get(rec.host)
        if behind is None or behind not in skin_a:
            raise refuse(Code.EDGE_FEATURE, "a recess in a face with no face of the other skin behind it",
                         faces=rec.faces)  # fmt: skip
        for face, rim in rec.rims:
            fid = table.next_id
            table.items.append(Feature(fid, "engraving", 0.0, _rim_diameter(rim, edges), "top" if on_a else "bottom",
                                       depth=rec.depth, face_ids=tuple(rec.faces)))  # fmt: skip
            if face in skin_a:
                table.edges[fid] = set(rim)
            else:
                table.projected_rims[fid] = (behind, rim)
    if recesses:
        deepest = max(r.depth for r in recesses)
        table.diags.append(warn(Code.ENGRAVED_RECESS,
                                f"{len(recesses)} recess(es) up to {deepest:.3g} deep engraved along their rims; "
                                "their depth is not made",
                                faces=sorted(i for r in recesses for i in r.faces),
                                recesses=len(recesses), depth=deepest))  # fmt: skip


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


def _wires(face: Face):
    """(outer wire, [inner wires]) of a face."""
    outer = BRepTools.OuterWire_s(face.shape)
    inner = []
    ex = TopExp_Explorer(face.shape, TopAbs_WIRE)
    while ex.More():
        if not ex.Current().IsSame(outer):
            inner.append(ex.Current())
        ex.Next()
    return outer, inner


def _wire_edges(wire, face: Face, host: int, edges) -> list[int]:
    """Edge indices of a wire, in order along it."""
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


def _classify(reg: _Formed, faces: list[Face], edges, nb_all, in_skin, t: float, wire, seeded: bool) -> None:
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


def _slit_width(reg: _Formed, faces: list[Face], edges, in_skin) -> float:
    """Width of the slits, seen on the sheet: the region's free edges projected onto the host
    plane, and their distance from the sheet's cut edges (the median, so the ends where the two
    meet do not count).
    """
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
    faces: list[Face], skin: set[int], s: int, smooth_nb, nb_all, edges, in_skin, t: float
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
    ex = TopExp_Explorer(wire, TopAbs_EDGE)
    while ex.More():
        if ex.Current().IsSame(e):
            return True
        ex.Next()
    return False


def _collar(reg: _Formed, faces: list[Face]):
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


def _region_extent(reg: _Formed, faces: list[Face]) -> tuple[float, float]:
    """Smallest and largest signed distance of the region from its host plane (along its normal)."""
    assert reg.host is not None
    h = faces[reg.host - 1]
    n = h.normal.Normalized()
    o = h.ad.Plane().Location()
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


def _host_circle(reg: _Formed, faces: list[Face], edges) -> tuple[float, Any] | None:
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


def _formed_feature(fid: int, reg: _Formed, faces: list[Face], edges, t: float, diags: list) -> Feature:
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


def _recess(group: list[int], skin: set[int], faces: list[Face], edges, t: float, tol: Tolerances):
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


# ------------------------------------------------------------------------------ bends and flanges


def _bend_groups(faces: list[Face], edges, skin_a: set[int], sharp, vmap, t: float, tol: Tolerances):
    """The bend faces and sharp bend edges of the reference skin, grouped into bends.

    Cylinders group by axis line, radius and angle; sharp edges by edge line and angle.
    Returns [(kind, members)] in bend id order, the members carrying their bend id.
    """
    groups: dict[tuple, list[BendFace | SharpEdge]] = defaultdict(list)
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
        key = ("cyl", *_axis_key(f), round(r_in, 3), round(math.degrees(ang), 2), concave)
        groups[key].append(BendFace(i, 0, r_in, concave, ang))

    for ei, (a, b) in sharp.items():
        if a not in skin_a:
            continue
        e = edges[ei][0]
        p0, p1 = _edge_line(e)
        fa, fb = faces[a - 1], faces[b - 1]
        ang = fa.normal.Angle(fb.normal)
        concave = _concave(e, fa, fb, max(1e-3, 1e-3 * t))
        vs = edge_vertices(e)
        se = SharpEdge(ei, (a, b), ang, concave, p0.Distance(p1), (vmap.FindIndex(vs[0]), vmap.FindIndex(vs[-1])))
        key = ("sharp", *_line_key(p0, gp_Vec(p0, p1)), round(math.degrees(ang), 2), concave)
        groups[key].append(se)

    def first(members) -> tuple[int, int]:  # bends are numbered by their first face, sharp ones after
        return min((0, m.face) if isinstance(m, BendFace) else (1, m.edge) for m in members)

    ordered = sorted(groups.items(), key=lambda kv: first(kv[1]))
    return [(key[0], [replace(m, bend_id=n + 1) for m in members]) for n, (key, members) in enumerate(ordered)]


def _flanges(faces: list[Face], skin_a: set[int], bend_faces, smooth) -> tuple[tuple[Flange, ...], dict[int, int]]:
    """Components of skin A once the bend faces are removed (sharp edges separate too)."""
    flat_ids = {i for i in skin_a if i not in bend_faces}
    links = [(a, b) for a, b, _ in smooth if a in flat_ids and b in flat_ids]
    comps = sorted(sorted(c) for c in _components(len(faces), links) if c[0] in flat_ids)
    flange_of = {i: n + 1 for n, c in enumerate(comps) for i in c}
    flanges = tuple(Flange(n + 1, tuple(c), sum(faces[i - 1].area for i in c)) for n, c in enumerate(comps))
    return flanges, flange_of


def _bend(kind: str, members, faces: list[Face], edges, flange_of: dict[int, int], g: _Graph) -> Bend:
    segs = []
    adj_flanges: set[int] = set()
    m0 = members[0]
    if kind == "cyl":
        for m in sorted(members, key=lambda m: m.face):
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
            adj_flanges |= {flange_of[nb] for nb in g.smooth_nb[m.face] if nb in flange_of}
        deg = math.degrees(m0.angle)
        bkind: Literal["bend", "hem", "sharp"] = "sharp" if m0.r_in == 0.0 else ("hem" if deg >= 179.0 else "bend")
        r_in = m0.r_in
    else:
        for m in sorted(members, key=lambda m: m.edge):
            p0, p1 = _edge_line(edges[m.edge][0])
            segs.append(
                BendSegment(0, m.length, (p0.X(), p0.Y(), p0.Z()), tuple(gp_Vec(p0, p1).Normalized().Coord()), m.edge)
            )
            adj_flanges |= {flange_of[x] for x in m.faces if x in flange_of}
        deg, bkind, r_in = math.degrees(m0.angle), "sharp", 0.0
    return Bend(m0.bend_id, tuple(segs), r_in, deg, "up" if m0.concave else "down", bkind, tuple(sorted(adj_flanges)))


def _skin_adjacency(smooth, sharp, skin_a: set[int]) -> dict[int, list[tuple[int, int]]]:
    """Face -> [(neighbour, edge)] over smooth edges and the sharp bends of the reference skin."""
    adj: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for a, b, ei in smooth:
        adj[a].append((b, ei))
        adj[b].append((a, ei))
    for ei, (a, b) in sharp.items():
        if a in skin_a:
            adj[a].append((b, ei))
            adj[b].append((a, ei))
    return adj


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
    g = _graph(faces, edges, tol)

    # sharp bends join faces of one skin across a non-tangent edge
    plane_partners = _plane_partners(faces, tol)
    t_est: float | None = None
    partners: dict[int, int] = {}
    sharp: dict[int, tuple[int, int]] = {}
    if plane_partners:
        t_est, band = _mode([(d, faces[i - 1].area) for i, (d, _) in plane_partners.items()], tol)
        partners = {i: p for i, (d, p) in plane_partners.items() if abs(d - t_est) <= band}
        sharp = _sharp_bend_edges(faces, edges, partners, g, t_est, tol)

    comps, skins = _skins(faces, [(a, b) for a, b, _ in g.smooth] + list(sharp.values()))
    in_skin = skins[0] | skins[1]
    regions: list[_Formed] = []
    t_guess = thickness or t_est or 1.0  # to tell a slit from an opening; refined below
    for s, sk in enumerate(skins):
        regions += _formed_regions(faces, sk, s, g.smooth_nb, g.nb_all, edges, in_skin, t_guess)
    for reg in regions:
        if reg.kind is None:
            raise refuse(Code.FORMED_FEATURE, reg.reason, faces=sorted(reg.faces))
    formed = {i for reg in regions if not reg.developed for i in reg.faces}
    _check_skin_surfaces(faces, skins, formed)

    t, conf = _thickness(faces, skins, plane_partners, tol)
    t = _checked_thickness(t, conf, thickness, tol)
    feature_groups, recesses = _between_skins(faces, edges, comps, skins, g, t, engrave, tol)

    planes = [f for f in faces if f.kind == "plane" and f.idx in in_skin and f.idx not in formed]
    if not planes:
        raise refuse(Code.UNSUPPORTED_SURFACE, "no planar skin face to start from")
    root = max(planes, key=lambda f: (f.area, -f.idx))
    a_side = 0 if root.idx in skins[0] else 1
    skin_a = skins[a_side] - formed

    bend_groups = _bend_groups(faces, edges, skin_a, sharp, vmap, t, tol)
    bend_faces = {m.face: m for kind, ms in bend_groups if kind == "cyl" for m in ms if isinstance(m, BendFace)}
    sharp_edges = {m.edge: m for kind, ms in bend_groups if kind == "sharp" for m in ms if isinstance(m, SharpEdge)}
    flanges, flange_of = _flanges(faces, skin_a, bend_faces, g.smooth)
    bends = tuple(_bend(kind, ms, faces, edges, flange_of, g) for kind, ms in bend_groups)
    table = _Features()
    _hole_features(table, feature_groups, faces, g, skin_a)
    _formed_features(table, regions, a_side, faces, edges, t, bends, g)
    _engravings(table, recesses, skin_a, partners, edges)
    diags.extend(table.diags)
    features = table.items

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
    ctx = Context(
        shape,
        fmap,
        vmap,
        faces,
        edges,
        _skin_adjacency(g.smooth, sharp, skin_a),
        skin_a,
        root.idx,
        t,
        bend_faces,
        sharp_edges,
        features=table.faces,
        feature_through_r=table.through_r,
        feature_on_a=table.on_a,
        feature_kind={f.id: f.kind for f in features},
        feature_edges=table.edges,
        feature_cuts=table.cuts,
        feature_gap=table.gap,
        developed=table.developed,
        feature_free=table.free,
        projected_rims=table.projected_rims,
        tol=tol,
    )
    return SheetPart(
        thickness=t,
        thickness_confidence=conf,
        faces=face_infos,
        bends=bends,
        flanges=flanges,
        reference_root=root.idx,
        diagnostics=tuple(diags),
        source=source,
        features=tuple(features),
        _ctx=ctx,
    )
