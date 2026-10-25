"""Sheet metal recognition (docs/design.md §1.1–1.3).

Faces and their adjacency form an attributed adjacency graph. Edges along which two
faces meet tangentially are SMOOTH; the connected components over SMOOTH edges are the
two skins (the two largest) and the thickness faces. Bends are the cylinders of the
reference skin.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, replace
from typing import Any, Literal, cast

from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.BRepTools import BRepTools
from OCP.GeomAbs import GeomAbs_Circle, GeomAbs_Cone, GeomAbs_Cylinder, GeomAbs_Line, GeomAbs_Plane
from OCP.gp import gp_Pnt, gp_Vec
from OCP.GProp import GProp_GProps
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX
from OCP.TopExp import TopExp, TopExp_Explorer

from ._occ import EdgeFaceMap, ShapeMap, as_edge, as_face
from ._topology import BendFace, Context, Face, SharpEdge, edge_vertices, inner_uv, pcurve, surface_props
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
        kind = {GeomAbs_Plane: "plane", GeomAbs_Cylinder: "cylinder", GeomAbs_Cone: "cone"}.get(ad.GetType(), "other")
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
    """Face groups by area (largest first) and the two skins, which must be planes or cylinders."""
    comps = sorted(_components(len(faces), links), key=lambda c: -sum(faces[i - 1].area for i in c))
    if len(comps) < 3:
        raise refuse(Code.NOT_TWO_SKINS, f"{len(comps)} face groups; a sheet has two skins and edge faces")
    skins = [set(comps[0]), set(comps[1])]
    for sk in skins:
        bad = [i for i in sk if faces[i - 1].kind in ("other", "cone")]
        if bad:
            raise refuse(Code.UNSUPPORTED_SURFACE, "skin faces must be planes or cylinders", faces=bad)
    return comps, skins


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


def _hole_feature(group: list[int], faces: list[Face], tol: Tolerances):
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
        return "counterbore", rmin, rmax, None, depth
    return None


def _between_skins(faces: list[Face], comps, skins: list[set[int]], g: _Graph, tol: Tolerances):
    """Hole features among the faces between the skins: [(faces, feature tuple)].

    The other face groups between the skins must touch both (edge faces); one that
    touches only one is refused.
    """
    in_skin = skins[0] | skins[1]
    between = [i for i in range(1, len(faces) + 1) if i not in in_skin]
    groups = _components(len(faces), [(a, b) for a in between for b in g.nb_all[a] if b not in in_skin and a < b])
    groups = [grp for grp in groups if grp[0] not in in_skin]
    found: list[tuple[list[int], tuple]] = []
    featured: set[int] = set()
    for grp in groups:
        touches = {s for s in (0, 1) if any(n in skins[s] for i in grp for n in g.nb_all[i])}
        hf = _hole_feature(grp, faces, tol) if touches == {0, 1} else None
        if hf is not None:
            found.append((grp, hf))
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
    return found


def _features(feature_groups, faces: list[Face], g: _Graph, skin_a: set[int]):
    """Feature records, the per-feature data unfold() needs, and a warning per projected countersink."""
    features = []
    feat_faces: dict[int, tuple[int, ...]] = {}
    feat_r: dict[int, float] = {}
    feat_a: dict[int, bool] = {}
    diags = []
    for n, (grp, (fkind, rmin, rmax, angle, depth)) in enumerate(feature_groups):
        fid = n + 1
        # the side carrying the large rim: the skin touched by the face holding the largest circle
        big_face = max(grp, key=lambda i: max((r for r, _ in _circle_radii(faces[i - 1])), default=0))
        on_a = any(nb in skin_a for nb in g.nb_all[big_face])
        features.append(Feature(fid, fkind, 2 * rmin, 2 * rmax, "top" if on_a else "bottom", angle, depth, tuple(grp)))
        feat_faces[fid] = tuple(grp)
        feat_r[fid] = rmin
        feat_a[fid] = on_a
        if fkind == "countersink" and on_a:
            diags.append(warn(Code.COUNTERSINK_PROJECTED, f"countersink ⌀{2 * rmax:g} cut as ⌀{2 * rmin:g}", faces=grp))
    return features, feat_faces, feat_r, feat_a, diags


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
) -> SheetPart:
    """Recognise a solid as sheet metal: thickness, skins, bends, flanges, hole features.

    ``model`` is a :class:`~flatmat.io.StepModel` or an OCP solid. With ``canonical``
    (default), B-spline faces and edges that are really planes, cylinders, cones, lines
    or circles are converted first (on a copy).
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
    sharp: dict[int, tuple[int, int]] = {}
    if plane_partners:
        t_est, band = _mode([(d, faces[i - 1].area) for i, (d, _) in plane_partners.items()], tol)
        partners = {i: p for i, (d, p) in plane_partners.items() if abs(d - t_est) <= band}
        sharp = _sharp_bend_edges(faces, edges, partners, g, t_est, tol)

    comps, skins = _skins(faces, [(a, b) for a, b, _ in g.smooth] + list(sharp.values()))
    t, conf = _thickness(faces, skins, plane_partners, tol)
    t = _checked_thickness(t, conf, thickness, tol)
    feature_groups = _between_skins(faces, comps, skins, g, tol)

    planes = [f for f in faces if f.kind == "plane" and f.idx in skins[0] | skins[1]]
    if not planes:
        raise refuse(Code.UNSUPPORTED_SURFACE, "no planar skin face to start from")
    root = max(planes, key=lambda f: (f.area, -f.idx))
    a_side = 0 if root.idx in skins[0] else 1
    skin_a = skins[a_side]

    bend_groups = _bend_groups(faces, edges, skin_a, sharp, vmap, t, tol)
    bend_faces = {m.face: m for kind, ms in bend_groups if kind == "cyl" for m in ms if isinstance(m, BendFace)}
    sharp_edges = {m.edge: m for kind, ms in bend_groups if kind == "sharp" for m in ms if isinstance(m, SharpEdge)}
    flanges, flange_of = _flanges(faces, skin_a, bend_faces, g.smooth)
    bends = tuple(_bend(kind, ms, faces, edges, flange_of, g) for kind, ms in bend_groups)
    features, feat_faces, feat_r, feat_a, feat_diags = _features(feature_groups, faces, g, skin_a)
    diags.extend(feat_diags)

    role = {i: "skin_a" for i in skin_a} | {i: "skin_b" for i in skins[1 - a_side]}
    face_infos = tuple(
        FaceInfo(
            f.idx,
            cast(Literal["plane", "cylinder", "other"], f.kind if f.kind != "cone" else "other"),
            f.area,
            cast(Literal["skin_a", "skin_b", "thickness"], role.get(f.idx, "thickness")),
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
        feat_faces,
        feat_r,
        feat_a,
        tol,
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
