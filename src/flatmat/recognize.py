"""Sheet metal recognition (docs/design.md §1.1–1.3).

Faces and their adjacency form an attributed adjacency graph. Edges along which two
faces meet tangentially are SMOOTH; the connected components over SMOOTH edges are the
two skins (the two largest) and the thickness faces. Bends are the cylinders of the
reference skin.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Literal, cast

from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.BRepLProp import BRepLProp_SLProps
from OCP.BRepTools import BRepTools
from OCP.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane
from OCP.gp import gp_Pnt, gp_Pnt2d, gp_Vec
from OCP.GProp import GProp_GProps
from OCP.IntTools import IntTools_FClass2d
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_IN, TopAbs_REVERSED, TopAbs_VERTEX
from OCP.TopExp import TopExp

from ._occ import EdgeFaceMap, ShapeMap, as_edge, as_face
from .diagnostics import Code, Diagnostic, refuse
from .io import StepModel
from .model import Bend, BendSegment, FaceInfo, Flange, SheetPart, Tolerances

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
    axis_d: Any = None  # gp_Vec for cylinders
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
        kind = {GeomAbs_Plane: "plane", GeomAbs_Cylinder: "cylinder"}.get(ad.GetType(), "other")
        rec = _Face(i, f, kind, g.Mass(), n, mirror, (u0, u1, v0, v1), ad)
        if kind == "cylinder":
            c = ad.Cylinder()
            rec.radius = c.Radius()
            rec.axis_p = c.Location()
            rec.axis_d = gp_Vec(c.Axis().Direction())
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


def _axis_key(f: _Face) -> tuple[float, ...]:
    d = f.axis_d.Normalized()
    if (d.X(), d.Y(), d.Z()) < (0, 0, 0):
        d.Reverse()
    o = gp_Vec(f.axis_p.XYZ())
    foot = o - d.Multiplied(o.Dot(d))
    return tuple(round(x, 3) for x in (d.X(), d.Y(), d.Z(), foot.X(), foot.Y(), foot.Z()))


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


# ------------------------------------------------------------------------------ public


def recognize(
    model: StepModel | Any,
    *,
    thickness: float | None = None,
    solid: int | None = None,
    tol: Tolerances | None = None,
) -> SheetPart:
    """Recognise a solid as sheet metal: thickness, skins, bends, flanges.

    ``model`` is a :class:`~flatmat.io.StepModel` or an OCP solid.
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

    fmap, faces = _face_table(shape)
    edges = _edge_table(shape, fmap)
    vmap = ShapeMap()
    TopExp.MapShapes_s(shape, TopAbs_VERTEX, vmap)

    smooth, adj = [], defaultdict(list)
    sharp_nb = defaultdict(set)
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
        if _is_smooth(e, faces[a - 1], faces[b - 1], tol):
            smooth.append((a, b))
            adj[a].append((b, ei))
            adj[b].append((a, ei))
        else:
            sharp_nb[a].add(b)
            sharp_nb[b].add(a)

    comps = sorted(_components(len(faces), smooth), key=lambda c: -sum(faces[i - 1].area for i in c))
    if len(comps) < 3:
        raise refuse(Code.NOT_TWO_SKINS, f"{len(comps)} smooth face groups; a sheet has two skins and edge faces")
    skins = [set(comps[0]), set(comps[1])]
    for s in skins:
        bad = [i for i in s if faces[i - 1].kind == "other"]
        if bad:
            raise refuse(Code.UNSUPPORTED_SURFACE, "skin faces must be planes or cylinders in v0.1", faces=bad)

    t, conf = _thickness(faces, skins, tol)
    if thickness is not None:
        if abs(thickness - t) > max(0.01, tol.thickness_rel * thickness):
            raise refuse(Code.VARIABLE_THICKNESS, f"given thickness {thickness} but measured {t:.4f}", measured=t)
        t = thickness
    if conf < tol.thickness_confidence:
        raise refuse(
            Code.VARIABLE_THICKNESS, f"only {conf:.0%} of the sheet agrees on t = {t:.4f}", thickness=t, confidence=conf
        )

    for comp in comps[2:]:
        touches = {s for s in (0, 1) if any(n in skins[s] for i in comp for n in sharp_nb[i])}
        if touches != {0, 1}:
            raise refuse(
                Code.EDGE_FEATURE,
                "a face group between the skins touches only one skin (chamfer, countersink, formed feature?)",
                faces=comp,
            )

    planes = [f for f in faces if f.kind == "plane" and (f.idx in skins[0] or f.idx in skins[1])]
    if not planes:
        raise refuse(Code.UNSUPPORTED_SURFACE, "no planar skin face to start from")
    root = max(planes, key=lambda f: (f.area, -f.idx))
    a_side = 0 if root.idx in skins[0] else 1
    skin_a = skins[a_side]

    # bends: cylinders of the reference skin, grouped by axis line, radius and angle
    groups: dict[tuple, list[_BendFace]] = defaultdict(list)
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
        ang = u1 - u0
        if math.degrees(ang) >= 359.9:
            raise refuse(Code.BEND_TOO_LARGE, f"bend of {math.degrees(ang):.1f}°", faces=(i,))
        bf = _BendFace(i, 0, r_in, concave, ang)
        key = (*_axis_key(f), round(r_in, 3), round(math.degrees(ang), 2), concave)
        groups[key].append(bf)
        bend_faces[i] = bf

    # flanges: components of skin A once the bend faces are removed
    flat_ids = {i for i in skin_a if i not in bend_faces}
    flange_links = [(a, b) for a, b in smooth if a in flat_ids and b in flat_ids]
    flange_comps = sorted(sorted(c) for c in _components(len(faces), flange_links) if c[0] in flat_ids)
    flange_of = {i: n + 1 for n, c in enumerate(flange_comps) for i in c}
    flanges = tuple(Flange(n + 1, tuple(c), sum(faces[i - 1].area for i in c)) for n, c in enumerate(flange_comps))

    bends = []
    for n, (_key, members) in enumerate(sorted(groups.items(), key=lambda kv: min(m.face for m in kv[1]))):
        bid = n + 1
        segs = []
        adj_flanges = set()
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
            adj_flanges |= {flange_of[nb] for nb, _ in adj[m.face] if nb in flange_of}
        m0 = members[0]
        deg = math.degrees(m0.angle)
        bends.append(
            Bend(
                id=bid,
                segments=tuple(segs),
                r_in=m0.r_in,
                angle_deg=deg,
                direction="up" if m0.concave else "down",
                kind="hem" if deg >= 179.0 else "bend",
                flanges=tuple(sorted(adj_flanges)),
            )
        )

    role = {i: "skin_a" for i in skin_a} | {i: "skin_b" for i in skins[1 - a_side]}
    face_infos = tuple(
        FaceInfo(
            f.idx,
            cast(Literal["plane", "cylinder", "other"], f.kind),
            f.area,
            cast(Literal["skin_a", "skin_b", "thickness"], role.get(f.idx, "thickness")),
        )
        for f in faces
    )
    ctx = Context(shape, fmap, vmap, faces, edges, adj, skin_a, root.idx, t, bend_faces, tol)
    return SheetPart(
        thickness=t,
        thickness_confidence=conf,
        faces=face_infos,
        bends=tuple(bends),
        flanges=flanges,
        reference_root=root.idx,
        diagnostics=tuple(diags),
        source=source,
        _ctx=ctx,
    )
