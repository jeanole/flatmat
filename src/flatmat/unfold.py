"""Unfolding (docs/design.md §1.4–1.6).

Every face of the reference skin is developed on its own into a local 2D frame
(planes: their own (u, v); cylinders: (R_n·u, v)), then placed by a rigid 2D
transform that maps the shared vertices of an edge onto the already placed
neighbour. The outline is chained from the skin's boundary edges by vertex
identity, so no booleans and no fuzzy matching are involved.
"""

from __future__ import annotations

import math
from collections import defaultdict, deque

from OCP.BRep import BRep_Tool
from OCP.Geom2d import Geom2d_TrimmedCurve
from OCP.Geom2dConvert import Geom2dConvert
from OCP.GeomAbs import GeomAbs_BSplineCurve, GeomAbs_Circle, GeomAbs_Line
from OCP.TopAbs import TopAbs_VERTEX
from OCP.TopExp import TopExp_Explorer
from shapely.geometry import Polygon
from shapely.validation import explain_validity

from ._occ import OCP_VERSION, as_vertex
from .bendmodel import as_k_model, bend_allowance, bend_deduction, outside_setback
from .diagnostics import Code, Diagnostic, refuse, warn
from .model import Arc, BendLine, FlatPattern, Line, Loop2D, Point, SheetPart, Spline
from .recognize import Context, _Face, pcurve

Rigid = tuple[float, float, float, float]  # cos, sin, tx, ty
_IDENTITY: Rigid = (1.0, 0.0, 0.0, 0.0)


def _apply(T: Rigid, p: Point) -> Point:
    c, s, tx, ty = T
    return (c * p[0] - s * p[1] + tx, s * p[0] + c * p[1] + ty)


def _rigid(a1: Point, a2: Point, b1: Point, b2: Point) -> Rigid:
    """Rotation + translation taking b1 -> a1 and b2 -> a2 (never a reflection)."""
    ang = math.atan2(a2[1] - a1[1], a2[0] - a1[0]) - math.atan2(b2[1] - b1[1], b2[0] - b1[0])
    c, s = math.cos(ang), math.sin(ang)
    return (c, s, a1[0] - (c * b1[0] - s * b1[1]), a1[1] - (s * b1[0] + c * b1[1]))


class _Dev:
    """(u, v) on one face of the reference skin -> that face's developed 2D frame."""

    def __init__(self, f: _Face, rn: float | None):
        self.m = f.mirror
        self.s = rn if (f.kind == "cylinder" and rn is not None) else 1.0

    def __call__(self, u: float, v: float) -> Point:
        return (self.m * self.s * u, v)

    @property
    def flips(self) -> bool:  # does the map reverse orientation?
        return self.m < 0


def _vertices(e) -> list:
    out = []
    ex = TopExp_Explorer(e, TopAbs_VERTEX)
    while ex.More():
        out.append(as_vertex(ex.Current()))
        ex.Next()
    return out


def _vertex_uv(v, e, f) -> tuple[float, float]:
    p = pcurve(e, f).Value(BRep_Tool.Parameter_s(v, e, f))
    return p.X(), p.Y()


def unfold(
    sheet: SheetPart, k=0.44, *, sharp_radius: float = 0.0, check_overlap: bool = True, pre_hole=None
) -> FlatPattern:
    """Develop the recognised sheet into a flat pattern with bend lines.

    ``k`` is a K-factor (float), a :class:`~flatmat.bendmodel.KTable`,
    ``KTable.din6935()`` or a callable ``(r_in, t, angle_rad) -> K``.

    ``sharp_radius`` is the inner radius the shop actually gets where the model has a
    sharp bend (r = 0), e.g. from the die: air bending gives ~16 % of V on mild steel.

    ``pre_hole`` sets the pre-hole of extruded holes: a diameter, or a callable
    ``(Feature) -> diameter`` (e.g. a lookup in the tool maker's table). By default it
    comes from the neutral fibre (:func:`~flatmat.bendmodel.extruded_hole_prehole`).
    """
    ctx: Context = sheet._ctx  # type: ignore[assignment]
    if ctx is None:
        raise ValueError("SheetPart has no recognition context; call flatmat.recognize() first")
    kmodel = as_k_model(k)
    t = ctx.thickness
    faces = ctx.faces
    diags: list[Diagnostic] = list(sheet.diagnostics)

    # K per bend; a sharp bend (r = 0 in the model) is made with sharp_radius
    k_of_bend: dict[int, float] = {}
    r_of_bend: dict[int, float] = {}
    for b in sheet.bends:
        r = sharp_radius if b.kind == "sharp" else b.r_in
        kk, d = kmodel(r, t, math.radians(b.angle_deg))
        k_of_bend[b.id], r_of_bend[b.id] = kk, r
        if d is not None:
            diags.append(d)
    if sharp_radius == 0.0 and any(b.kind == "sharp" for b in sheet.bends):
        diags.append(warn(Code.SHARP_BEND_ASSUMED_RADIUS, "sharp bends developed with inner radius 0"))

    def sharp_delta(bend_id: int, angle: float, concave: bool) -> float:
        """Strip width between the two faces of a sharp corner (negative: overlap).

        The faces of the reference skin run to the virtual sharp. A real bend of radius r
        has BA between tangent points that sit r·tan(φ/2) (reference skin inside) or
        (r + t)·tan(φ/2) (outside) back from that sharp.
        """
        r = r_of_bend[bend_id]
        setback = (r if concave else r + t) * math.tan(angle / 2)
        return bend_allowance(angle, r, t, k_of_bend[bend_id]) - 2 * setback

    # neutral radius per bend face, from its own convexity (joggles: design §5.1)
    devs: dict[int, _Dev] = {}
    for i in ctx.skin_a:
        bf = ctx.bend_faces.get(i)
        rn = None
        if bf is not None:
            if bf.r_in == 0.0:
                # sharp inside, round (R = t) outside, reference skin outside: the flanges end
                # at the cylinder's tangent lines, which sit on the inner sharp, so the
                # cylinder carries the inside-sharp strip BA(r) - 2·r·tan(φ/2)
                rn = sharp_delta(bf.bend_id, bf.angle, True) / bf.angle
            else:
                rn = bf.r_in + k_of_bend[bf.bend_id] * t
        devs[i] = _Dev(faces[i - 1], rn)
    shifts = {ei: sharp_delta(se.bend_id, se.angle, se.concave) for ei, se in ctx.sharp_edges.items()}

    place, sharp_geo = _place(ctx, devs, shifts)
    sharp_v = {v: shifts[ei] for ei, se in ctx.sharp_edges.items() for v in se.vertices}
    sheet_features = list(sheet.features)
    if pre_hole is not None:
        diags = [d for d in diags if d.code != Code.PRE_HOLE]  # the recognised pre-hole is replaced
        sheet_features = [
            _with(f, through_diameter=float(pre_hole(f) if callable(pre_hole) else pre_hole))
            if f.kind == "extruded_hole"
            else f
            for f in sheet_features
        ]
    pre_r = {f.id: f.through_diameter / 2 for f in sheet_features if f.kind == "extruded_hole"}
    outer, holes, feature_centers, outlines, hole_fid, slits, orient = _outline(ctx, devs, place, sharp_v, diags, pre_r)

    # self-checks
    area = outer.area() - sum(-h.area() for h in holes)
    expected = 0.0
    for i in ctx.skin_a:
        f = faces[i - 1]
        expected += f.area * (devs[i].s / f.radius if i in ctx.bend_faces else 1.0)
    expected += sum(shifts[ei] * se.length for ei, se in ctx.sharp_edges.items())
    expected += sum(extra for _, extra in feature_centers.values())  # rims replaced by through holes
    rel = abs(area - expected) / expected
    if rel > ctx.tol.area_rel:
        raise refuse(Code.AREA_MISMATCH, f"flat area {area:.4f} vs developed {expected:.4f}", rel_err=rel)
    if check_overlap:
        poly = Polygon(outer.polyline(0.005), [h.polyline(0.005) for h in holes])
        if not poly.is_valid:
            raise refuse(Code.FLAT_OVERLAP, f"the flat pattern overlaps itself: {explain_validity(poly)}")

    # bend lines and bend values
    lines: list[BendLine] = []
    bends = []
    for b in sheet.bends:
        kk, r = k_of_bend[b.id], r_of_bend[b.id]
        rad = math.radians(b.angle_deg)
        bends.append(
            _with(
                b,
                r_in=r,  # for a sharp bend: the radius it is made with
                k=kk,
                ba=bend_allowance(rad, r, t, kk),
                bd=bend_deduction(rad, r, t, kk),
                ossb=outside_setback(rad, r, t),
            )
        )
        for seg in b.segments:
            if seg.edge_id is not None:  # sharp bend: the middle of the inserted strip
                a1, a2, n, d = sharp_geo[seg.edge_id]
                h = d / 2
                p0, p1 = (a1[0] + h * n[0], a1[1] + h * n[1]), (a2[0] + h * n[0], a2[1] + h * n[1])
                lines.append(BendLine(b.id, 0, p0, p1, b.direction, b.angle_deg, r, b.kind))
                continue
            f = faces[seg.face_id - 1]
            for p0, p1 in _bend_line(f, devs[seg.face_id], place[seg.face_id]):
                lines.append(BendLine(b.id, seg.face_id, p0, p1, b.direction, b.angle_deg, r, b.kind))

    # normalise: bbox minimum at the origin
    pts = outer.polyline(0.01)
    dx, dy = -min(p[0] for p in pts), -min(p[1] for p in pts)
    shift = (1.0, 0.0, dx, dy)
    outer = _move_loop(outer, shift)
    holes = [_move_loop(h, shift) for h in holes]
    lines = [_with(bl, p0=_apply(shift, bl.p0), p1=_apply(shift, bl.p1)) for bl in lines]
    slit_loops = [(_move_loop(L, shift), fid) for L, fid in slits]
    features = tuple(
        _cut_feature(
            _with(
                f,
                center=_apply(shift, feature_centers[f.id][0]) if f.id in feature_centers else None,
                outline=_move_loop(outlines[f.id], shift) if f.id in outlines else None,
                holes=tuple(i for i, fid in enumerate(hole_fid) if fid == f.id),
            ),
            holes,
            lines,
            orient.get(f.id),
        )
        for f in sheet_features
    )

    return FlatPattern(
        outer=outer,
        holes=tuple(holes),
        bend_lines=tuple(lines),
        bends=tuple(bends),
        thickness=t,
        k_model=kmodel.describe(),
        area=area,
        diagnostics=tuple(diags),
        meta={
            "source": sheet.source,
            "ocp": OCP_VERSION,
            "area_rel_err": rel,
            "sharp_radius": sharp_radius,
            "root_face": ctx.root,
            "root_offset": (dx, dy),  # flat = root face (u, v) developed + this
        },
        features=features,
        slits=tuple(L for L, _ in slit_loops),
    )


def _cut_feature(f, holes: list[Loop2D], lines: list[BendLine], angle: float | None):
    """Louvres, bridges, lances: where they point and how big they are in the flat.

    The direction is that of the slits (louvre, bridge) or of the bend line (lance).
    """
    if f.kind not in ("louvre", "bridge", "lance"):
        return f
    if f.kind == "lance":
        pts = [p for i in f.holes for p in holes[i].polyline(0.01)]
        bl = next((b for b in lines if b.bend_id == f.bend_id), None)
        if bl is not None:
            angle = math.atan2(bl.p1[1] - bl.p0[1], bl.p1[0] - bl.p0[0])
    else:
        pts = f.outline.polyline(0.01) if f.outline is not None else []
    if not pts or angle is None:
        return f
    deg = math.degrees(angle) % 180.0
    if deg > 180.0 - 1e-9:
        deg = 0.0
    a = math.radians(deg)
    u, v = (math.cos(a), math.sin(a)), (-math.sin(a), math.cos(a))
    pa = [x * u[0] + y * u[1] for x, y in pts]
    pb = [x * v[0] + y * v[1] for x, y in pts]
    if f.center is None:  # a lance: the middle of its opening
        ca, cb = (max(pa) + min(pa)) / 2, (max(pb) + min(pb)) / 2
        f = _with(f, center=(ca * u[0] + cb * v[0], ca * u[1] + cb * v[1]))
    return _with(f, orientation_deg=deg, size=(max(pa) - min(pa), max(pb) - min(pb)))


def _with(obj, **kw):
    from dataclasses import replace

    return replace(obj, **kw)


# ------------------------------------------------------------------------------ placement


def _inner_uv(f: _Face) -> tuple[float, float]:
    u0, u1, v0, v1 = f.uv
    for su in (0.5, 0.25, 0.75, 0.1, 0.9):
        for sv in (0.5, 0.25, 0.75, 0.1, 0.9):
            u, v = u0 + su * (u1 - u0), v0 + sv * (v1 - v0)
            if f.inside(u, v):
                return u, v
    return 0.5 * (u0 + u1), 0.5 * (v0 + v1)


def _place(ctx: Context, devs: dict[int, _Dev], shifts: dict[int, float]):
    """Rigid placement of every reference-skin face; a sharp edge inserts a strip of shifts[ei]."""
    faces, edges = ctx.faces, ctx.edges
    place = {ctx.root: _IDENTITY}
    sharp_geo: dict[int, tuple[Point, Point, Point, float]] = {}
    q = deque([ctx.root])
    seen = set()
    worst = 0.0
    while q:
        p = q.popleft()
        for c, ei in ctx.adj[p]:
            if c not in ctx.skin_a or ei in seen:
                continue
            seen.add(ei)
            e = edges[ei][0]
            vs = _vertices(e)
            if len(vs) < 2 or vs[0].IsSame(vs[-1]):
                continue  # closed edge: carries no placement information
            fp, fc = faces[p - 1].shape, faces[c - 1].shape
            a1, a2 = (_apply(place[p], devs[p](*_vertex_uv(v, e, fp))) for v in (vs[0], vs[-1]))
            b1, b2 = (devs[c](*_vertex_uv(v, e, fc)) for v in (vs[0], vs[-1]))
            if abs(math.dist(a1, a2) - math.dist(b1, b2)) > ctx.tol.flat:
                raise refuse(
                    Code.DEVELOP_FAILED,
                    "a shared edge has different lengths in its two developments",
                    faces=(p, c),
                    edges=(ei,),
                )
            Tc = _rigid(a1, a2, b1, b2)
            if ei in shifts:
                # move C away from P, perpendicular to the edge, by the strip width
                L = math.dist(a1, a2)
                n = (-(a2[1] - a1[1]) / L, (a2[0] - a1[0]) / L)
                qc = _apply(Tc, devs[c](*_inner_uv(faces[c - 1])))
                if (qc[0] - a1[0]) * n[0] + (qc[1] - a1[1]) * n[1] < 0:
                    n = (-n[0], -n[1])
                d = shifts[ei]
                Tc = (Tc[0], Tc[1], Tc[2] + d * n[0], Tc[3] + d * n[1])
                sharp_geo[ei] = (a1, a2, n, d)
            if c in place:
                worst = max(worst, *(math.dist(_apply(place[c], b), _apply(Tc, b)) for b in (b1, b2)))
                if worst > ctx.tol.flat:
                    raise refuse(
                        Code.CYCLE_INCONSISTENT,
                        "faces around a closed loop do not develop consistently (tube? needs a cut)",
                        faces=(p, c),
                        edges=(ei,),
                        gap=worst,
                    )
                continue
            place[c] = Tc
            q.append(c)
    if len(place) != len(ctx.skin_a):
        missing = sorted(ctx.skin_a - set(place))
        raise refuse(Code.DEVELOP_FAILED, "part of the skin is not connected to the root", faces=missing)
    return place, sharp_geo


# ------------------------------------------------------------------------------ outline


def _edge_segment(e, ei: int, f: _Face, dev: _Dev, T: Rigid):
    """The developed, placed, exact 2D geometry of edge e on face f."""
    c2 = pcurve(e, f.shape)
    a, b = c2.FirstParameter(), c2.LastParameter()

    def m(p) -> Point:
        return _apply(T, dev(p.X(), p.Y()))

    kind = c2.GetType()
    if kind == GeomAbs_Line:
        return Line(m(c2.Value(a)), m(c2.Value(b)), ei)
    if kind == GeomAbs_Circle and f.kind == "plane":
        circ = c2.Circle()
        cen = m(circ.Location())
        p0, pm = m(c2.Value(a)), m(c2.Value((a + b) / 2))
        a0 = math.atan2(p0[1] - cen[1], p0[0] - cen[0])
        am = math.atan2(pm[1] - cen[1], pm[0] - cen[0])
        sweep = abs(b - a)
        # direction: does the arc reach the mid point turning counter-clockwise?
        ccw = ((am - a0) % (2 * math.pi)) < math.pi if sweep < 2 * math.pi - 1e-9 else _circle_ccw(m, c2, a, b, cen)
        return Arc(cen, circ.Radius(), a0, a0 + (sweep if ccw else -sweep), ccw, ei)
    # everything else: an exact B-spline, developed by mapping its poles
    if kind == GeomAbs_BSplineCurve:
        bs = c2.BSpline()
        bs = bs.Copy() if hasattr(bs, "Copy") else bs
        bs = _downcast_bspline(bs)
        if bs.IsPeriodic():
            bs.SetNotPeriodic()  # DXF needs a clamped knot vector; this is exact (adds poles)
        if a > bs.FirstParameter() + 1e-12 or b < bs.LastParameter() - 1e-12:
            bs.Segment(a, b)
    else:
        curve = BRep_Tool.CurveOnSurface_s(e, f.shape, 0.0, 0.0)
        bs = Geom2dConvert.CurveToBSplineCurve_s(Geom2d_TrimmedCurve(curve, a, b))
    if bs.IsPeriodic():
        bs.SetNotPeriodic()
    n = bs.NbPoles()
    poles = tuple(m(bs.Pole(i)) for i in range(1, n + 1))
    weights = tuple(bs.Weight(i) for i in range(1, n + 1)) if bs.IsRational() else None
    knots = tuple(bs.Knot(i) for i in range(1, bs.NbKnots() + 1))
    mults = tuple(bs.Multiplicity(i) for i in range(1, bs.NbKnots() + 1))
    return Spline(bs.Degree(), poles, knots, mults, weights, ei)


def _downcast_bspline(h):
    from OCP.Geom2d import Geom2d_BSplineCurve

    if isinstance(h, Geom2d_BSplineCurve):
        return h
    dc = getattr(Geom2d_BSplineCurve, "DownCast", None)
    return dc(h) if dc else h


def _circle_ccw(m, c2, a, b, cen) -> bool:
    p0, p1 = m(c2.Value(a)), m(c2.Value(a + (b - a) * 0.25))
    cross = (p0[0] - cen[0]) * (p1[1] - cen[1]) - (p0[1] - cen[1]) * (p1[0] - cen[0])
    return cross > 0


def _outline(ctx: Context, devs, place, sharp_v: dict[int, float], diags: list, pre_r: dict[int, float]):
    faces = ctx.faces
    feature_of_face = {i: fid for fid, fs in ctx.features.items() for i in fs}
    feature_of_edge: dict[int, int] = {}
    pieces = []  # (v_start, v_end, segment)
    for ei, (e, fi) in enumerate(ctx.edges):
        inside = [i for i in fi if i in ctx.skin_a]
        if len(inside) != 1 or BRep_Tool.Degenerated_s(e):
            continue
        i = inside[0]
        f = faces[i - 1]
        if BRep_Tool.IsClosed_s(e, f.shape):
            continue  # seam
        other = [j for j in fi if j != i]
        if other and other[0] in feature_of_face:
            feature_of_edge[ei] = feature_of_face[other[0]]
        seg = _edge_segment(e, ei, f, devs[i], place[i])
        vs = _vertices(e)
        va, vb = ctx.vmap.FindIndex(vs[0]), ctx.vmap.FindIndex(vs[-1])
        pa = _apply(place[i], devs[i](*_vertex_uv(vs[0], e, f.shape)))
        if va != vb and math.dist(seg.start, pa) > math.dist(seg.end, pa):
            seg = seg.reversed()
        pieces.append((va, vb, seg))
    for opening, es in ctx.feature_edges.items():  # openings of louvres, bridges, lances
        for ei in es:
            feature_of_edge.setdefault(ei, opening)
    cut_edges = {ei for chains in ctx.feature_cuts.values() for c in chains for ei in c}
    loops = [
        _merge_collinear(_close(segs, jv, sharp_v, ctx.tol.flat), keep_apart=cut_edges) for segs, jv in _chain(pieces)
    ]

    # hole features: report the centre; where the large rim is on the reference skin,
    # cut the through hole instead (the countersink / counterbore is made afterwards).
    # Formed features: an emboss's rim is a marking, not a cut (the area inside it is
    # material); an extruded hole cuts its pre-hole; a louvre or bridge cuts its slits,
    # the rest of its rim is a marking; a lance's opening is cut as it is.
    centers: dict[int, tuple[Point, float]] = {}
    outlines: dict[int, Loop2D] = {}
    slits: list[tuple[Loop2D, int]] = []
    orient: dict[int, float] = {}  # louvres, bridges: direction of the slits
    out: list[tuple[Loop2D, int | None]] = []
    for L in loops:
        fids = {feature_of_edge[s.source_edge] for s in L.segments if s.source_edge in feature_of_edge}
        fid: int | None = next(iter(fids)) if len(fids) == 1 else None
        if fid is None:
            out.append((L, None))
            continue
        kind = ctx.feature_kind.get(fid)
        if kind == "lance":
            out.append((L, fid))
            continue
        if kind in ("louvre", "bridge"):
            rim = L if L.area() > 0 else _reverse_loop(L)
            outlines[fid] = rim
            cut_area = 0.0
            runs = _cut_runs(rim, cut_edges)
            if runs:
                c = max(runs, key=lambda ch: math.dist(ch[0].start, ch[-1].end))
                orient[fid] = math.atan2(c[-1].end[1] - c[0].start[1], c[-1].end[0] - c[0].start[0])
            for chain in runs:
                closing = Line(chain[-1].end, chain[0].start, chain[0].source_edge)
                slot = Loop2D((*chain, closing))
                if abs(slot.area()) <= 1e-9 * max(1.0, closing.length()) ** 2:
                    slits.append((Loop2D(tuple(chain)), fid))  # no width: a cut line
                else:
                    out.append((slot, fid))
                    cut_area += abs(slot.area())
            centers[fid] = (_centroid(rim), rim.area() - cut_area)
            continue
        if kind == "emboss":
            rim = L if L.area() > 0 else _reverse_loop(L)
            outlines[fid] = rim
            centers[fid] = (_centroid(rim), rim.area())
            continue
        if kind == "extruded_hole":
            seg = L.segments[0]
            if len(L.segments) != 1 or not isinstance(seg, Arc) or not seg.is_circle:
                raise refuse(Code.FORMED_FEATURE, "the rim of an extruded hole is not a circle in the flat")
            outlines[fid] = Loop2D((Arc(seg.center, seg.radius, 0.0, 2 * math.pi, True, seg.source_edge),))
            r = pre_r.get(fid, 0.0)
            centers[fid] = (seg.center, math.pi * (seg.radius**2 - r**2))
            if r > 0:
                out.append((Loop2D((Arc(seg.center, r, 0.0, -2 * math.pi, False, seg.source_edge),)), fid))
            continue
        if len(L.segments) != 1 or not isinstance(L.segments[0], Arc) or not L.segments[0].is_circle:
            diags.append(warn(Code.FEATURE_NOT_CIRCULAR, "a hole feature is not a circle in the flat; kept as is"))
            out.append((L, fid))
            continue
        arc = L.segments[0]
        if ctx.feature_on_a.get(fid):
            r = ctx.feature_through_r[fid]
            centers[fid] = (arc.center, math.pi * (arc.radius**2 - r**2))
            out.append((Loop2D((Arc(arc.center, r, arc.a0, arc.a1, arc.ccw, arc.source_edge),)), fid))
        else:
            centers[fid] = (arc.center, 0.0)
            out.append((L, fid))
    out.sort(key=lambda lf: -abs(lf[0].area()))
    outer = out[0][0]
    if outer.area() < 0:
        outer = _reverse_loop(outer)
    holes = [_reverse_loop(h) if h.area() > 0 else h for h, _ in out[1:]]
    hole_fid = [f for _, f in out[1:]]
    return outer, holes, centers, outlines, hole_fid, slits, orient


def _cut_runs(L: Loop2D, cut_edges: set[int]) -> list[list]:
    """Runs of consecutive segments of a closed loop that come from cut edges."""
    segs = list(L.segments)
    flags = [s.source_edge in cut_edges for s in segs]
    if all(flags) or not any(flags):
        return [segs] if all(flags) else []
    k = flags.index(False)
    segs, flags = segs[k:] + segs[:k], flags[k:] + flags[:k]
    runs: list[list] = []
    cur: list = []
    for sg, f in zip(segs, flags, strict=True):
        if f:
            cur.append(sg)
        elif cur:
            runs.append(cur)
            cur = []
    if cur:
        runs.append(cur)
    return runs


def _centroid(L: Loop2D) -> Point:
    pts = L.polyline(0.001)
    a = cx = cy = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1], strict=True):
        c = x0 * y1 - x1 * y0
        a += c
        cx += (x0 + x1) * c
        cy += (y0 + y1) * c
    return (cx / (3 * a), cy / (3 * a))


def _chain(pieces) -> list[tuple[list, list[int]]]:
    """Chain pieces into loops by vertex identity: (segments, vertex at the start of each)."""
    by_v = defaultdict(list)
    for n, (a, b, _) in enumerate(pieces):
        by_v[a].append(n)
        by_v[b].append(n)
    used: set[int] = set()
    loops = []
    for start in range(len(pieces)):
        if start in used:
            continue
        a, b, seg = pieces[start]
        used.add(start)
        segs, jv = [seg], [a]
        cur = b
        while cur != a:
            nxt = next((n for n in by_v[cur] if n not in used), None)
            if nxt is None:
                raise refuse(Code.OUTLINE_OPEN, "the flat outline does not close", edges=(seg.source_edge,))
            used.add(nxt)
            na, nb, s = pieces[nxt]
            jv.append(cur)
            if na == cur:
                segs.append(s)
                cur = nb
            else:
                segs.append(s.reversed())
                cur = na
        loops.append((segs, jv))
    return loops


def _close(segs: list, jv: list[int], sharp_v: dict[int, float], tol: float) -> Loop2D:
    """Make consecutive segments meet.

    At a sharp-bend vertex the two sides were moved apart (strip, delta > 0: bridge the
    gap with a line) or into each other (delta < 0: trim both lines to their
    intersection, or join them if they are collinear). Any other gap is an error.
    """
    segs = list(segs)
    if len(segs) == 1:
        return Loop2D(tuple(segs))
    out: list = []
    n = len(segs)
    for i in range(n):
        prev, cur, v = segs[i - 1], segs[i], jv[i]
        if math.dist(prev.end, cur.start) <= tol:
            continue
        if v not in sharp_v:
            raise refuse(Code.OUTLINE_OPEN, "gap in the flat outline", edges=(cur.source_edge,))
        if sharp_v[v] >= 0:
            continue  # bridged below, in order
        if not (isinstance(prev, Line) and isinstance(cur, Line)):
            raise refuse(
                Code.SHARP_BEND_UNSUPPORTED,
                "the outline at a sharp bend is not straight; cannot trim the overlap",
                edges=(prev.source_edge, cur.source_edge),
            )
        x = _intersect(prev, cur)
        if x is None:  # collinear: drop the overlap by joining the two
            x = cur.start if _along(prev, cur.start) < _along(prev, prev.end) else prev.end
            segs[i - 1] = Line(prev.p0, x, prev.source_edge)
            segs[i] = Line(x, cur.p1, cur.source_edge)
        else:
            segs[i - 1] = Line(prev.p0, x, prev.source_edge)
            segs[i] = Line(x, cur.p1, cur.source_edge)
    for i in range(n):
        prev, cur = segs[i - 1], segs[i]
        if math.dist(prev.end, cur.start) > tol:  # a strip (delta > 0): bridge it
            out.append(Line(prev.end, cur.start, cur.source_edge))
        out.append(cur)
    return Loop2D(tuple(s for s in out if not (isinstance(s, Line) and s.length() <= tol)))


def _along(line: Line, p: Point) -> float:
    d = (line.p1[0] - line.p0[0], line.p1[1] - line.p0[1])
    return (p[0] - line.p0[0]) * d[0] + (p[1] - line.p0[1]) * d[1]


def _intersect(a: Line, b: Line) -> Point | None:
    d1 = (a.p1[0] - a.p0[0], a.p1[1] - a.p0[1])
    d2 = (b.p1[0] - b.p0[0], b.p1[1] - b.p0[1])
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) <= 1e-9 * (math.hypot(*d1) * math.hypot(*d2) or 1.0):
        return None
    u = ((b.p0[0] - a.p0[0]) * d2[1] - (b.p0[1] - a.p0[1]) * d2[0]) / den
    return (a.p0[0] + u * d1[0], a.p0[1] + u * d1[1])


def _merge_collinear(L: Loop2D, tol: float = 1e-9, keep_apart: set[int] | None = None) -> Loop2D:
    """Join consecutive collinear lines and co-circular arcs (split at bend tangents, circle halves).

    Segments from ``keep_apart`` edges (slits) are not joined with others.
    """
    segs = list(L.segments)
    if len(segs) < 2:
        return L
    apart = keep_apart or set()

    def joinable(a, b) -> bool:
        if (a.source_edge in apart) != (b.source_edge in apart):
            return False
        if isinstance(a, Line) and isinstance(b, Line):
            (x0, y0), (x1, y1), (x2, y2) = a.p0, a.p1, b.p1
            cross = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
            dot = (x1 - x0) * (x2 - x1) + (y1 - y0) * (y2 - y1)
            return abs(cross) <= tol * max(1.0, math.dist(a.p0, b.p1)) ** 2 and dot > 0
        if isinstance(a, Arc) and isinstance(b, Arc):
            # CAD systems split circles into arcs (SolidWorks: two halves)
            return (
                a.ccw == b.ccw
                and math.dist(a.center, b.center) <= 1e-7 * max(1.0, a.radius)
                and abs(a.radius - b.radius) <= 1e-7 * max(1.0, a.radius)
                and a.sweep + b.sweep <= 2 * math.pi + 1e-9
            )
        return False

    def join(a, b):
        if isinstance(a, Line):
            return Line(a.p0, b.p1, a.source_edge)
        d = b.sweep if a.ccw else -b.sweep
        return Arc(a.center, a.radius, a.a0, a.a1 + d, a.ccw, a.source_edge)

    out = [segs[0]]
    for s in segs[1:]:
        if joinable(out[-1], s):
            out[-1] = join(out[-1], s)
        else:
            out.append(s)
    while len(out) > 1 and joinable(out[-1], out[0]):
        out[0] = join(out[-1], out[0])
        out.pop()
    return Loop2D(tuple(out))


def _reverse_loop(L: Loop2D) -> Loop2D:
    return Loop2D(tuple(s.reversed() for s in reversed(L.segments)))


def _move_loop(L: Loop2D, T: Rigid) -> Loop2D:
    out = []
    for s in L.segments:
        if isinstance(s, Line):
            out.append(_with(s, p0=_apply(T, s.p0), p1=_apply(T, s.p1)))
        elif isinstance(s, Arc):
            out.append(_with(s, center=_apply(T, s.center)))
        else:
            out.append(_with(s, poles=tuple(_apply(T, p) for p in s.poles)))
    return Loop2D(tuple(out))


# ------------------------------------------------------------------------------ bend lines


def _bend_line(f: _Face, dev: _Dev, T: Rigid, steps: int = 200) -> list[tuple[Point, Point]]:
    """The iso-u line at mid angle, clipped to the face.

    Coarse scan along v for inside/outside, then bisection on each transition so the
    end points are exact to ~1e-9 of the bend length.
    """
    u0, u1, v0, v1 = f.uv
    um = 0.5 * (u0 + u1)
    eps = (v1 - v0) * 1e-9

    def ins(v: float) -> bool:
        return f.inside(um, v)

    def edge(lo: float, hi: float) -> float:  # ins(lo) != ins(hi)
        want = ins(lo)
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if ins(mid) == want:
                lo = mid
            else:
                hi = mid
            if hi - lo < eps:
                break
        return 0.5 * (lo + hi)

    vs = [v0 + (v1 - v0) * i / steps for i in range(steps + 1)]
    flags = [ins(v) for v in vs]
    # the boundary itself classifies as ON; treat the ends as inside if their neighbour is
    flags[0] = flags[0] or flags[1]
    flags[-1] = flags[-1] or flags[-2]
    segs = []
    start = vs[0] if flags[0] else None
    for i in range(1, len(vs)):
        if flags[i] and not flags[i - 1]:
            start = edge(vs[i - 1], vs[i])
        elif not flags[i] and flags[i - 1] and start is not None:
            segs.append((start, edge(vs[i - 1], vs[i])))
            start = None
    if start is not None:
        segs.append((start, vs[-1]))
    return [(_apply(T, dev(um, a)), _apply(T, dev(um, b))) for a, b in segs if b - a > eps]
