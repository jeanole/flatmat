"""Phase-0 spike of the flatmat pipeline (docs/design.md §1).

Not the library: one file, no error taxonomy, planes and cylinders only. It
exists to check the design end to end on real OCP:

    read STEP -> adjacency graph -> skins -> thickness -> bends
    -> develop each face on its own -> place by shared vertices
    -> outline by vertex chaining -> bend lines -> self-checks -> DXF

Every stage is timed.
"""

from __future__ import annotations

import math
import os
import tempfile
import time
from collections import defaultdict, deque
from dataclasses import dataclass, field

from occ_compat import EdgeFaceMap, ShapeMap, as_edge, as_face, as_solid, as_vertex
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.BRepLProp import BRepLProp_SLProps
from OCP.BRepTools import BRepTools
from OCP.GeomAbs import GeomAbs_BSplineCurve, GeomAbs_Cylinder, GeomAbs_Line, GeomAbs_Plane
from OCP.gp import gp_Pnt, gp_Pnt2d, gp_Vec
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.IntTools import IntTools_FClass2d
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_IN, TopAbs_REVERSED, TopAbs_SOLID, TopAbs_VERTEX
from OCP.TopExp import TopExp, TopExp_Explorer

SMOOTH_DEG = 0.1  # design §6.2: smooth < 0.1°, sharp > 1°


class Refused(Exception):
    pass


@dataclass
class Face:
    idx: int
    shape: object
    kind: str  # "plane" | "cylinder" | "other"
    area: float
    normal: gp_Vec  # outward, orientation applied, at the UV centre
    mirror: float  # +1 if surface normal == outward normal, else -1
    uv: tuple[float, float, float, float]
    ad: object = None
    _cls: object = None
    radius: float = 0.0
    axis_p: gp_Pnt | None = None
    axis_d: gp_Vec | None = None


@dataclass
class Bend:
    face: int
    r_in: float
    angle_deg: float
    direction: str
    k: float
    rn: float
    ba: float
    line: tuple = ()
    segments: list = field(default_factory=list)


@dataclass
class Result:
    thickness: float
    thickness_confidence: float
    bends: list[Bend]
    outer: list[tuple[float, float]]
    holes: list[list[tuple[float, float]]]
    flat_area: float
    expected_area: float
    timings: dict[str, float]
    n_faces: int
    checks: dict[str, float] = field(default_factory=dict)
    bspline_max_err: float = 0.0


# ----------------------------------------------------------------------------- io


# STEP import healing (docs/design.md §6.1). The reader's ShapeFix pass runs
# ShapeFix_Shell, which is quadratic in the face count (~4 s for a plate with 1000
# holes). "fast" switches the shell and solid fixes off and validates afterwards;
# an invalid result is re-read with full healing.
#   OCP 8.0: XSControl_Reader.SetShapeFixParameters, set after ReadFile (the actor
#            that stores them only exists then), as a complete key -> value map.
#   OCP 7.8: no such method; the FromSTEP sequence comes from a resource file named
#            "STEP" in $CSF_STEPDefaults, which must be set before the first read.
_FAST_OVERRIDES = {"FixShape.FixShellMode": "0", "FixShape.FixSolidMode": "0"}
_RESOURCE_78 = """FromSTEP.exec.op : FixShape
FromSTEP.FixShape.Tolerance3d : &Runtime.Tolerance
FromSTEP.FixShape.MaxTolerance3d : &Runtime.MaxTolerance
FromSTEP.FixShape.MinTolerance3d : 1.e-7
FromSTEP.FixShape.FixShellMode : 0
FromSTEP.FixShape.FixSolidMode : 0
ToSTEP.exec.op : SplitCommonVertex,DirectFaces
"""


def _has_fix_parameters():
    return hasattr(STEPControl_Reader, "SetShapeFixParameters")


if not _has_fix_parameters() and "CSF_STEPDefaults" not in os.environ:
    _res_dir = tempfile.mkdtemp(prefix="flatmat_occt78_")
    with open(os.path.join(_res_dir, "STEP"), "w") as fh:
        fh.write(_RESOURCE_78)
    os.environ["CSF_STEPDefaults"] = _res_dir
    _FAST_78 = True
else:
    _FAST_78 = False


def _apply_fast(reader):
    from OCP.DE import DE_ShapeFixParameters
    from OCP.TCollection import TCollection_AsciiString as A

    reader.SetShapeFixParameters(DE_ShapeFixParameters())  # defaults into the map
    params = reader.GetShapeFixParameters()
    for k, v in _FAST_OVERRIDES.items():
        params.Bind(A(k), A(v))
    reader.SetShapeFixParameters(params)


def read_step(path: str, heal: str = "fast"):
    r = STEPControl_Reader()
    if r.ReadFile(path) != IFSelect_RetDone:
        raise Refused(f"cannot read {path}")
    if heal == "fast" and _has_fix_parameters():
        _apply_fast(r)
    elif heal == "full" and _FAST_78:
        raise Refused("full healing on OCP 7.8 needs a fresh process without CSF_STEPDefaults")
    r.TransferRoots()
    shape = r.OneShape()
    solids = []
    ex = TopExp_Explorer(shape, TopAbs_SOLID)
    while ex.More():
        solids.append(as_solid(ex.Current()))
        ex.Next()
    if len(solids) != 1:
        raise Refused(f"MULTIPLE_SOLIDS: {len(solids)}")
    solid = solids[0]
    if heal == "fast":
        # No BRepCheck_Analyzer here: on a face with 1000 holes it costs 1.3 s (OCP 8.0)
        # and 13 s (7.8) even topology-only. Our own self-checks (area identity, closed
        # outline, cycle consistency) catch what matters; BRepCheck is for diagnostics.
        g = GProp_GProps()
        BRepGProp.VolumeProperties_s(solid, g)
        if g.Mass() < 0:  # the solid fix would have reoriented it
            solid = as_solid(solid.Reversed())
    return solid


# ------------------------------------------------------------------- face table


def _props(face, ad, u, v):
    pr = BRepLProp_SLProps(ad, u, v, 1, 1e-7)
    n = gp_Vec(pr.Normal())
    sn = pr.D1U().Crossed(pr.D1V())
    if face.Orientation() == TopAbs_REVERSED:
        n.Reverse()
    return n, (1.0 if sn.Dot(n) > 0 else -1.0)


def face_table(solid):
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
        st = ad.GetType()
        kind = {GeomAbs_Plane: "plane", GeomAbs_Cylinder: "cylinder"}.get(st, "other")
        rec = Face(i, f, kind, g.Mass(), n, mirror, (u0, u1, v0, v1), ad)
        if kind == "cylinder":
            c = ad.Cylinder()
            rec.radius = c.Radius()
            rec.axis_p = c.Location()
            rec.axis_d = gp_Vec(c.Axis().Direction())
        faces.append(rec)
    return fmap, faces


def edge_table(solid, fmap):
    efm = EdgeFaceMap()
    TopExp.MapShapesAndAncestors_s(solid, TopAbs_EDGE, TopAbs_FACE, efm)
    edges = []  # (edge, [face idx ...]) with 1-based face idx
    for i in range(1, efm.Extent() + 1):
        e = as_edge(efm.FindKey(i))
        lst = efm.FindFromIndex(i)
        # Iterating an OCP shape list from Python costs ~0.3 ms per list (both 7.8 and
        # 8.0): a second on 3000 edges. First()/Last() cover the manifold case.
        if lst.Size() <= 2:
            idx = sorted({fmap.FindIndex(lst.First()), fmap.FindIndex(lst.Last())})
        else:
            idx = sorted({fmap.FindIndex(f) for f in lst})
        edges.append((e, idx))
    return edges


def classify_edge(e, f1: Face, f2: Face) -> bool:
    """True if the two faces meet tangentially (SMOOTH) along e."""
    c = BRepAdaptor_Curve(e)
    tm = 0.5 * (c.FirstParameter() + c.LastParameter())
    normals = []
    for f in (f1, f2):
        c2 = BRepAdaptor_Curve2d(e, f.shape)
        uv = c2.Value(tm)
        n, _ = _props(f.shape, f.ad, uv.X(), uv.Y())
        normals.append(n)
    return normals[0].Angle(normals[1]) < math.radians(SMOOTH_DEG)


# -------------------------------------------------------------- recognition


def components(n, links):
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


def thickness(solid, fmap, faces, skin_ids, tol=1e-4, max_planes=12):
    cands = []  # (value, weight)
    # Coaxial cylinder pairs. Bucket by the axis line (canonical direction + foot point
    # nearest the origin, rounded) so only candidates on the same line are compared:
    # comparing all pairs was O(n^2) for the hole walls of a perforated part.
    buckets = defaultdict(list)
    for f in faces:
        if f.kind != "cylinder":
            continue
        d = f.axis_d.Normalized()
        if (d.X(), d.Y(), d.Z()) < (0, 0, 0):
            d.Reverse()
        o = gp_Vec(f.axis_p.XYZ())
        foot = o - d.Multiplied(o.Dot(d))
        key = tuple(round(x, 3) for x in (d.X(), d.Y(), d.Z(), foot.X(), foot.Y(), foot.Z()))
        buckets[key].append(f)
    for group in buckets.values():
        for i, a in enumerate(group):
            for b in group[i + 1 :]:
                dr = abs(a.radius - b.radius)
                if dr > tol:
                    cands.append((dr, min(a.area, b.area)))
    # Plane pairs, analytically: a point inside a skin plane, its distance to each
    # antiparallel plane of the other skin, and whether its foot lies inside that face.
    # (A ray cast against the whole solid was O(n^2) on perforated parts.)
    planes = sorted((f for f in faces if f.kind == "plane" and f.idx in skin_ids), key=lambda f: -f.area)
    planes = planes[:max_planes]
    for f in planes:
        p = _inner_point(f)
        if p is None:
            continue
        best = None
        for g in planes:
            if g is f or g.normal.Angle(f.normal) < math.pi - 1e-3:
                continue
            d = abs(gp_Vec(g.ad.Plane().Location(), p).Dot(g.normal) / g.normal.Magnitude())
            if d <= tol or (best is not None and d >= best):
                continue
            foot = gp_Pnt(p.XYZ() - f.normal.Normalized().XYZ() * d)
            if _contains(g, foot):
                best = d
        if best is not None:
            cands.append((best, f.area))
    if not cands:
        raise Refused("NO_THICKNESS")
    bins = defaultdict(float)
    for v, w in cands:
        bins[round(v, 3)] += w
    t = max(bins, key=bins.get)
    total = sum(bins.values())
    agree = sum(w for v, w in cands if abs(v - t) <= max(0.01, 0.01 * t))
    return t, agree / total


def _plane_uv(f: Face, p: gp_Pnt):
    pos = f.ad.Plane().Position()
    d = gp_Vec(pos.Location(), p)
    return d.Dot(gp_Vec(pos.XDirection())), d.Dot(gp_Vec(pos.YDirection()))


def inside(f: Face, u: float, v: float) -> bool:
    """Point-in-face test with a classifier built once per face (IntTools_FClass2d).

    BRepClass_FaceClassifier rebuilds everything per call: ~9 ms on a face with
    1000 holes, against 0.8 ms per Perform after a 21 ms build here.
    """
    if f._cls is None:
        f._cls = IntTools_FClass2d(f.shape, 1e-7)
    return f._cls.Perform(gp_Pnt2d(u, v)) == TopAbs_IN


def _contains(f: Face, p: gp_Pnt) -> bool:
    return inside(f, *_plane_uv(f, p))


def _inner_point(f: Face):
    u0, u1, v0, v1 = f.uv
    for su in (0.5, 0.25, 0.75):
        for sv in (0.5, 0.25, 0.75):
            u, v = u0 + su * (u1 - u0), v0 + sv * (v1 - v0)
            if inside(f, u, v):
                return f.ad.Value(u, v)
    return None


# -------------------------------------------------------------- development


class Dev:
    """Map (u, v) on one face of skin A to that face's own developed 2D frame."""

    def __init__(self, f: Face, rn: float | None):
        self.f = f
        self.s = rn if f.kind == "cylinder" else 1.0

    def __call__(self, u, v):
        return (self.f.mirror * self.s * u, v)


def rigid(a1, a2, b1, b2):
    """Rotation + translation taking b1->a1 and b2->a2 (no reflection)."""
    ang = math.atan2(a2[1] - a1[1], a2[0] - a1[0]) - math.atan2(b2[1] - b1[1], b2[0] - b1[0])
    c, s = math.cos(ang), math.sin(ang)
    tx = a1[0] - (c * b1[0] - s * b1[1])
    ty = a1[1] - (s * b1[0] + c * b1[1])
    return (c, s, tx, ty)


def apply(T, p):
    c, s, tx, ty = T
    return (c * p[0] - s * p[1] + tx, s * p[0] + c * p[1] + ty)


def edge_vertices(e):
    vs = []
    ex = TopExp_Explorer(e, TopAbs_VERTEX)
    while ex.More():
        vs.append(as_vertex(ex.Current()))
        ex.Next()
    return vs


def vertex_uv(v, e, f):
    t = BRep_Tool.Parameter_s(v, e, f)
    p = BRepAdaptor_Curve2d(e, f).Value(t)
    return p.X(), p.Y()


def sample_edge(e, f, n):
    c2 = BRepAdaptor_Curve2d(e, f)
    a, b = c2.FirstParameter(), c2.LastParameter()
    return [(lambda p: (p.X(), p.Y()))(c2.Value(a + (b - a) * i / n)) for i in range(n + 1)], c2


def samples_for(c2) -> int:
    return 1 if c2.GetType() == GeomAbs_Line else 64


def bspline_check(c2, dev: Dev) -> float:
    """Exactness of the pole-scaling development of a B-spline pcurve.

    Builds the developed curve by mapping the poles (affine map, weights kept)
    and compares it with mapping evaluated points. Returns the max deviation.
    """
    if c2.GetType() != GeomAbs_BSplineCurve:
        return 0.0
    bs = c2.BSpline()
    cp = bs.Copy()
    from OCP.Geom2d import Geom2d_BSplineCurve

    cp = Geom2d_BSplineCurve.DownCast(cp) if hasattr(Geom2d_BSplineCurve, "DownCast") else cp
    for i in range(1, bs.NbPoles() + 1):
        p = bs.Pole(i)
        x, y = dev(p.X(), p.Y())
        cp.SetPole(i, gp_Pnt2d(x, y))
    err = 0.0
    a, b = c2.FirstParameter(), c2.LastParameter()
    for i in range(51):
        t = a + (b - a) * i / 50
        p = bs.Value(t)
        x, y = dev(p.X(), p.Y())
        q = cp.Value(t)
        err = max(err, math.hypot(q.X() - x, q.Y() - y))
    return err


# ----------------------------------------------------------------- pipeline


def unfold(solid, k: float = 0.44) -> Result:
    T = {}
    t0 = time.perf_counter()

    def lap(name):
        nonlocal t0
        now = time.perf_counter()
        T[name] = (now - t0) * 1000
        t0 = now

    fmap, faces = face_table(solid)
    edges = edge_table(solid, fmap)
    vmap = ShapeMap()
    TopExp.MapShapes_s(solid, TopAbs_VERTEX, vmap)
    lap("tables")

    smooth, adj = [], defaultdict(list)
    for ei, (e, fi) in enumerate(edges):
        if BRep_Tool.Degenerated_s(e):
            continue
        if len(fi) == 1:
            if not BRep_Tool.IsClosed_s(e, faces[fi[0] - 1].shape):
                raise Refused("OPEN_SHELL")
            continue
        if len(fi) > 2:
            raise Refused("NON_MANIFOLD")
        a, b = fi
        if classify_edge(e, faces[a - 1], faces[b - 1]):
            smooth.append((a, b))
            adj[a].append((b, ei))
            adj[b].append((a, ei))
    lap("classify")

    comps = sorted(components(len(faces), smooth), key=lambda c: -sum(faces[i - 1].area for i in c))
    if len(comps) < 2:
        raise Refused("NOT_TWO_SKINS")
    skins = comps[:2]
    skin_of = {i: s for s, comp in enumerate(skins) for i in comp}
    lap("skins")

    t, conf = thickness(solid, fmap, faces, set(skins[0]) | set(skins[1]))
    if conf < 0.95:
        raise Refused(f"VARIABLE_THICKNESS (confidence {conf:.2f})")
    lap("thickness")

    planes = [f for f in faces if f.kind == "plane" and f.idx in skin_of]
    root = max(planes, key=lambda f: f.area)
    A = set(skins[skin_of[root.idx]])
    for i in A:
        if faces[i - 1].kind == "other":
            raise Refused(f"UNSUPPORTED_SURFACE face {i}")

    bends = {}
    for i in A:
        f = faces[i - 1]
        if f.kind != "cylinder":
            continue
        u0, u1, v0, v1 = f.uv
        p = f.ad.Value((u0 + u1) / 2, (v0 + v1) / 2)
        ap = gp_Vec(f.axis_p, p)
        radial = ap - f.axis_d.Multiplied(ap.Dot(f.axis_d) / f.axis_d.SquareMagnitude())
        concave = f.normal.Dot(radial) < 0  # normal points to the axis: A is the inside
        r_in = f.radius if concave else f.radius - t
        rn = r_in + k * t
        ang = u1 - u0
        bends[i] = Bend(i, r_in, math.degrees(ang), "up" if concave else "down", k, rn, ang * rn)
    devs = {i: Dev(faces[i - 1], bends[i].rn if i in bends else None) for i in A}
    lap("bends")

    # place: BFS from root over smooth edges inside A
    place = {root.idx: (1.0, 0.0, 0.0, 0.0)}
    q = deque([root.idx])
    nontree = []
    seen_edges = set()
    while q:
        p = q.popleft()
        for c, ei in adj[p]:
            if c not in A or ei in seen_edges:
                continue
            seen_edges.add(ei)
            e = edges[ei][0]
            vs = edge_vertices(e)
            if len(vs) < 2:
                continue
            fp, fc = faces[p - 1].shape, faces[c - 1].shape
            a1, a2 = (apply(place[p], devs[p](*vertex_uv(v, e, fp))) for v in vs[:2])
            b1, b2 = (devs[c](*vertex_uv(v, e, fc)) for v in vs[:2])
            if abs(math.dist(a1, a2) - math.dist(b1, b2)) > 1e-6:
                raise Refused(f"DEVELOP_FAILED edge {ei}: {math.dist(a1, a2)} vs {math.dist(b1, b2)}")
            Tc = rigid(a1, a2, b1, b2)
            if c in place:
                nontree.append((ei, max(math.dist(apply(place[c], b), apply(Tc, b)) for b in (b1, b2))))
                continue
            place[c] = Tc
            q.append(c)
    if len(place) != len(A):
        raise Refused(f"DEVELOP_FAILED: placed {len(place)} of {len(A)} faces")
    worst_cycle = max((d for _, d in nontree), default=0.0)
    if worst_cycle > 1e-3:
        raise Refused(f"CYCLE_INCONSISTENT {worst_cycle:.4f}")
    lap("place")

    # outline: edges of A with exactly one face in A, chained by vertex identity
    bspline_err = 0.0
    pieces = []  # (v_start, v_end, points)
    for e, fi in edges:
        inside = [i for i in fi if i in A]
        if len(inside) != 1 or BRep_Tool.Degenerated_s(e):
            continue
        i = inside[0]
        f = faces[i - 1].shape
        if BRep_Tool.IsClosed_s(e, f):
            continue
        pts, c2 = sample_edge(e, f, 1)
        n = samples_for(c2)
        if n > 1:
            pts, _ = sample_edge(e, f, n)
        bspline_err = max(bspline_err, bspline_check(c2, devs[i]))
        pts = [apply(place[i], devs[i](*uv)) for uv in pts]
        vs = edge_vertices(e)
        a = vmap.FindIndex(vs[0])
        b = vmap.FindIndex(vs[-1])
        # orient the samples so they run from vertex a to vertex b
        pa = apply(place[i], devs[i](*vertex_uv(vs[0], e, f)))
        if math.dist(pts[0], pa) > math.dist(pts[-1], pa):
            pts.reverse()
        pieces.append((a, b, pts))
    loops = chain(pieces)
    lap("outline")

    loops.sort(key=lambda L: -abs(shoelace(L)))
    outer, holes = loops[0], loops[1:]
    flat_area = abs(shoelace(outer)) - sum(abs(shoelace(h)) for h in holes)
    expected = sum(faces[i - 1].area * (bends[i].rn / faces[i - 1].radius if i in bends else 1.0) for i in A)

    for i, b in bends.items():
        b.segments = bend_line(faces[i - 1], devs[i], place[i])
    lap("bendlines")

    return Result(
        t,
        conf,
        list(bends.values()),
        outer,
        holes,
        flat_area,
        expected,
        T,
        len(faces),
        checks={"cycle": worst_cycle},
        bspline_max_err=bspline_err,
    )


def chain(pieces):
    by_v = defaultdict(list)
    for n, (a, b, _) in enumerate(pieces):
        by_v[a].append(n)
        by_v[b].append(n)
    used, loops = set(), []
    for start in range(len(pieces)):
        if start in used:
            continue
        a, b, pts = pieces[start]
        used.add(start)
        loop = list(pts)
        cur = b
        while cur != a:
            nxt = next((n for n in by_v[cur] if n not in used), None)
            if nxt is None:
                raise Refused("OUTLINE_OPEN")
            used.add(nxt)
            na, nb, npts = pieces[nxt]
            if na == cur:
                loop += npts[1:]
                cur = nb
            else:
                loop += list(reversed(npts))[1:]
                cur = na
        loops.append(loop)
    return loops


def shoelace(pts):
    return 0.5 * sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1]))


def bend_line(f: Face, dev: Dev, T, steps=400):
    """Iso-u at u_mid, clipped to the face (interval search along v)."""
    u0, u1, v0, v1 = f.uv
    um = 0.5 * (u0 + u1)
    flags = []
    for i in range(steps + 1):
        v = v0 + (v1 - v0) * i / steps
        flags.append((v, inside(f, um, v) or i in (0, steps)))
    segs, start, prev = [], None, None
    for v, ok in flags:
        if ok and start is None:
            start = v
        if not ok and start is not None:
            segs.append((start, prev))
            start = None
        prev = v
    if start is not None:
        segs.append((start, prev))
    return [(apply(T, dev(um, a)), apply(T, dev(um, b))) for a, b in segs]
