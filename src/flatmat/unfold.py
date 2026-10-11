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
from .diagnostics import Code, Diagnostic, refuse
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


def unfold(sheet: SheetPart, k=0.44, *, check_overlap: bool = True) -> FlatPattern:
    """Develop the recognised sheet into a flat pattern with bend lines.

    ``k`` is a K-factor (float), a :class:`~flatmat.bendmodel.KTable`,
    ``KTable.din6935()`` or a callable ``(r_in, t, angle_rad) -> K``.
    """
    ctx: Context = sheet._ctx  # type: ignore[assignment]
    if ctx is None:
        raise ValueError("SheetPart has no recognition context; call flatmat.recognize() first")
    kmodel = as_k_model(k)
    t = ctx.thickness
    faces = ctx.faces
    diags: list[Diagnostic] = list(sheet.diagnostics)

    # neutral radius per bend face, from its own convexity (joggles: design §5.1)
    k_of_bend: dict[int, float] = {}
    for b in sheet.bends:
        kk, d = kmodel(b.r_in, t, math.radians(b.angle_deg))
        k_of_bend[b.id] = kk
        if d is not None:
            diags.append(d)
    devs: dict[int, _Dev] = {}
    for i in ctx.skin_a:
        bf = ctx.bend_faces.get(i)
        rn = None
        if bf is not None:
            rn = bf.r_in + k_of_bend[bf.bend_id] * t
        devs[i] = _Dev(faces[i - 1], rn)

    place = _place(ctx, devs)
    outer, holes = _outline(ctx, devs, place)

    # self-checks
    area = outer.area() - sum(-h.area() for h in holes)
    expected = 0.0
    for i in ctx.skin_a:
        f = faces[i - 1]
        bf = ctx.bend_faces.get(i)
        expected += f.area * ((bf.r_in + k_of_bend[bf.bend_id] * t) / f.radius if bf else 1.0)
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
        kk = k_of_bend[b.id]
        rad = math.radians(b.angle_deg)
        bends.append(
            _with(
                b,
                k=kk,
                ba=bend_allowance(rad, b.r_in, t, kk),
                bd=bend_deduction(rad, b.r_in, t, kk),
                ossb=outside_setback(rad, b.r_in, t),
            )
        )
        for seg in b.segments:
            f = faces[seg.face_id - 1]
            for p0, p1 in _bend_line(f, devs[seg.face_id], place[seg.face_id]):
                lines.append(BendLine(b.id, seg.face_id, p0, p1, b.direction, b.angle_deg, b.r_in, b.kind))

    # normalise: bbox minimum at the origin
    pts = outer.polyline(0.01)
    dx, dy = -min(p[0] for p in pts), -min(p[1] for p in pts)
    shift = (1.0, 0.0, dx, dy)
    outer = _move_loop(outer, shift)
    holes = [_move_loop(h, shift) for h in holes]
    lines = [_with(bl, p0=_apply(shift, bl.p0), p1=_apply(shift, bl.p1)) for bl in lines]

    return FlatPattern(
        outer=outer,
        holes=tuple(holes),
        bend_lines=tuple(lines),
        bends=tuple(bends),
        thickness=t,
        k_model=kmodel.describe(),
        area=area,
        diagnostics=tuple(diags),
        meta={"source": sheet.source, "ocp": OCP_VERSION, "area_rel_err": rel},
    )


def _with(obj, **kw):
    from dataclasses import replace

    return replace(obj, **kw)


# ------------------------------------------------------------------------------ placement


def _place(ctx: Context, devs: dict[int, _Dev]) -> dict[int, Rigid]:
    faces, edges = ctx.faces, ctx.edges
    place = {ctx.root: _IDENTITY}
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
    return place


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


def _outline(ctx: Context, devs: dict[int, _Dev], place: dict[int, Rigid]) -> tuple[Loop2D, list[Loop2D]]:
    faces = ctx.faces
    pieces = []  # (v_start, v_end, segment)
    for ei, (e, fi) in enumerate(ctx.edges):
        inside = [i for i in fi if i in ctx.skin_a]
        if len(inside) != 1 or BRep_Tool.Degenerated_s(e):
            continue
        i = inside[0]
        f = faces[i - 1]
        if BRep_Tool.IsClosed_s(e, f.shape):
            continue  # seam
        seg = _edge_segment(e, ei, f, devs[i], place[i])
        vs = _vertices(e)
        va, vb = ctx.vmap.FindIndex(vs[0]), ctx.vmap.FindIndex(vs[-1])
        pa = _apply(place[i], devs[i](*_vertex_uv(vs[0], e, f.shape)))
        if va != vb and math.dist(seg.start, pa) > math.dist(seg.end, pa):
            seg = seg.reversed()
        pieces.append((va, vb, seg))
    loops = [_merge_collinear(L) for L in _chain(pieces)]
    loops.sort(key=lambda L: -abs(L.area()))
    outer, holes = loops[0], loops[1:]
    if outer.area() < 0:
        outer = _reverse_loop(outer)
    holes = [_reverse_loop(h) if h.area() > 0 else h for h in holes]
    return outer, holes


def _chain(pieces) -> list[Loop2D]:
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
        segs = [seg]
        cur = b
        while cur != a:
            nxt = next((n for n in by_v[cur] if n not in used), None)
            if nxt is None:
                raise refuse(Code.OUTLINE_OPEN, "the flat outline does not close", edges=(seg.source_edge,))
            used.add(nxt)
            na, nb, s = pieces[nxt]
            if na == cur:
                segs.append(s)
                cur = nb
            else:
                segs.append(s.reversed())
                cur = na
        loops.append(Loop2D(tuple(segs)))
    return loops


def _merge_collinear(L: Loop2D, tol: float = 1e-9) -> Loop2D:
    """Join consecutive collinear lines and co-circular arcs (split at bend tangents, circle halves)."""
    segs = list(L.segments)
    if len(segs) < 2:
        return L

    def joinable(a, b) -> bool:
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
