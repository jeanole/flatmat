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
from shapely.geometry import Polygon
from shapely.validation import explain_validity

from ._occ import OCP_VERSION
from ._topology import Context, Face, edge_vertices, inner_uv, pcurve
from .bendmodel import as_k_model, bend_allowance, bend_deduction, outside_setback
from .diagnostics import Code, Diagnostic, refuse, warn
from .geom2d import IDENTITY, Rigid, cross, dot, heading, offset, perp, rigid_fit, sub, transform, unit
from .model import FEATURE_KINDS, Arc, BendLine, FlatPattern, Line, Loop2D, Point, SheetPart, Spline


class _Dev:
    """(u, v) on one face of the reference skin -> that face's developed 2D frame."""

    def __init__(self, f: Face, rn: float | None):
        self.m = f.mirror
        self.s = rn if (f.kind == "cylinder" and rn is not None) else 1.0

    def __call__(self, u: float, v: float) -> Point:
        return (self.m * self.s * u, v)

    @property
    def flips(self) -> bool:  # does the map reverse orientation?
        return self.m < 0


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
            if f.kind_info.flat == "pre_hole"
            else f
            for f in sheet_features
        ]
    pre_r = {f.id: f.through_diameter / 2 for f in sheet_features if f.kind_info.flat == "pre_hole"}
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
                p0, p1 = offset(a1, n, d / 2), offset(a2, n, d / 2)
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
    lines = [_with(bl, p0=transform(shift, bl.p0), p1=transform(shift, bl.p1)) for bl in lines]
    slit_loops = [(_move_loop(L, shift), fid) for L, fid in slits]
    features = tuple(
        _cut_feature(
            _with(
                f,
                center=transform(shift, feature_centers[f.id][0]) if f.id in feature_centers else None,
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

    The direction is that of the bend line (a bent lance) or of the chord of the longest
    slit (the others; for a formed lance, the side it is held on).
    """
    if not f.kind_info.oriented:
        return f
    if f.outline is not None:
        pts = f.outline.polyline(0.01)
    else:  # a bent lance: its opening, along its bend line
        pts = [p for i in f.holes for p in holes[i].polyline(0.01)]
        bl = next((b for b in lines if b.bend_id == f.bend_id), None)
        if bl is not None:
            angle = heading(bl.p0, bl.p1)
    if not pts or angle is None:
        return f
    deg = math.degrees(angle) % 180.0
    if deg > 180.0 - 1e-9:
        deg = 0.0
    a = math.radians(deg)
    u, v = (math.cos(a), math.sin(a)), (-math.sin(a), math.cos(a))
    pa = [dot(p, u) for p in pts]
    pb = [dot(p, v) for p in pts]
    if f.center is None:  # a lance: the middle of its opening
        ca, cb = (max(pa) + min(pa)) / 2, (max(pb) + min(pb)) / 2
        f = _with(f, center=(ca * u[0] + cb * v[0], ca * u[1] + cb * v[1]))
    return _with(f, orientation_deg=deg, size=(max(pa) - min(pa), max(pb) - min(pb)))


def _with(obj, **kw):
    from dataclasses import replace

    return replace(obj, **kw)


# ------------------------------------------------------------------------------ placement


def _place(ctx: Context, devs: dict[int, _Dev], shifts: dict[int, float]):
    """Rigid placement of every reference-skin face; a sharp edge inserts a strip of shifts[ei]."""
    faces, edges = ctx.faces, ctx.edges
    place = {ctx.root: IDENTITY}
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
            vs = edge_vertices(e)
            if len(vs) < 2 or vs[0].IsSame(vs[-1]):
                continue  # closed edge: carries no placement information
            fp, fc = faces[p - 1].shape, faces[c - 1].shape
            a1, a2 = (transform(place[p], devs[p](*_vertex_uv(v, e, fp))) for v in (vs[0], vs[-1]))
            b1, b2 = (devs[c](*_vertex_uv(v, e, fc)) for v in (vs[0], vs[-1]))
            if abs(math.dist(a1, a2) - math.dist(b1, b2)) > ctx.tol.flat:
                raise refuse(
                    Code.DEVELOP_FAILED,
                    "a shared edge has different lengths in its two developments",
                    faces=(p, c),
                    edges=(ei,),
                )
            Tc = rigid_fit(a1, a2, b1, b2)
            if ei in shifts:
                # move C away from P, perpendicular to the edge, by the strip width
                u, _ = unit(a1, a2)
                n = perp(u)
                fc_rec = faces[c - 1]
                qc = transform(Tc, devs[c](*(inner_uv(fc_rec) or fc_rec.uv_centre)))
                if dot(sub(qc, a1), n) < 0:
                    n = (-n[0], -n[1])
                d = shifts[ei]
                Tc = (Tc[0], Tc[1], Tc[2] + d * n[0], Tc[3] + d * n[1])
                sharp_geo[ei] = (a1, a2, n, d)
            if c in place:
                worst = max(worst, *(math.dist(transform(place[c], b), transform(Tc, b)) for b in (b1, b2)))
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


def _edge_segment(e, ei: int, f: Face, dev: _Dev, T: Rigid):
    """The developed, placed, exact 2D geometry of edge e on face f."""
    c2 = pcurve(e, f.shape)
    a, b = c2.FirstParameter(), c2.LastParameter()

    def m(p) -> Point:
        return transform(T, dev(p.X(), p.Y()))

    kind = c2.GetType()
    if kind == GeomAbs_Line:
        return Line(m(c2.Value(a)), m(c2.Value(b)), ei)
    if kind == GeomAbs_Circle and f.kind == "plane":
        circ = c2.Circle()
        cen = m(circ.Location())
        p0, pm = m(c2.Value(a)), m(c2.Value((a + b) / 2))
        a0, am = heading(cen, p0), heading(cen, pm)
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
    return cross(sub(p0, cen), sub(p1, cen)) > 0


def _outline(ctx: Context, devs, place, sharp_v: dict[int, float], diags: list, pre_r: dict[int, float]):
    faces = ctx.faces
    feature_of_face = {i: fid for fid, fd in ctx.features.items() for i in fd.outline_faces}
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
        vs = edge_vertices(e)
        va, vb = ctx.vmap.FindIndex(vs[0]), ctx.vmap.FindIndex(vs[-1])
        pa = transform(place[i], devs[i](*_vertex_uv(vs[0], e, f.shape)))
        if va != vb and math.dist(seg.start, pa) > math.dist(seg.end, pa):
            seg = seg.reversed()
        pieces.append((va, vb, seg))
    for opening, fd in ctx.features.items():  # openings of louvres, bridges, lances; rims of engravings
        for ei in fd.edges:
            feature_of_edge.setdefault(ei, opening)
    cut_edges = {ei for fd in ctx.features.values() for c in fd.cuts for ei in c}
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
        fd = ctx.features[fid]
        flat = FEATURE_KINDS[fd.kind].flat
        if fd.developed:  # a lance laid flat: its opening is cut as it is
            out.append((L, fid))
            continue
        if flat == "slits":
            rim = L if L.area() > 0 else _reverse_loop(L)
            outlines[fid] = rim
            cut_area = 0.0
            runs = _cut_runs(rim, cut_edges)
            if runs:
                c = max(runs, key=lambda ch: math.dist(ch[0].start, ch[-1].end))
                orient[fid] = heading(c[0].start, c[-1].end)
            host, free = fd.free
            if fd.gap <= 1e-6 or not free:
                slits.extend((Loop2D(tuple(chain)), fid) for chain in runs)  # no width: cut lines
            else:
                edge_pts = [_on_host(ctx, devs, place, host, ctx.edges[ei][0]) for ei in free]
                for slot in _slots(rim, runs, edge_pts):
                    out.append((slot, fid))
                    cut_area += abs(slot.area())
            centers[fid] = (_centroid(rim), rim.area() - cut_area)
            continue
        if flat == "rim":
            rim = L if L.area() > 0 else _reverse_loop(L)
            outlines[fid] = rim
            centers[fid] = (_centroid(rim), rim.area())
            continue
        if flat == "pre_hole":
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
        if fd.on_a:
            r = fd.through_r
            centers[fid] = (arc.center, math.pi * (arc.radius**2 - r**2))
            out.append((Loop2D((Arc(arc.center, r, arc.a0, arc.a1, arc.ccw, arc.source_edge),)), fid))
        else:
            centers[fid] = (arc.center, 0.0)
            out.append((L, fid))
    # rims of engravings on the far skin or round an island: projected onto the plane behind them
    for fid, fd in ctx.features.items():
        if fd.projected_rim is None:
            continue
        behind, rim_edges = fd.projected_rim
        rim = _projected_rim(ctx, devs, place, behind, rim_edges)
        rim = rim if rim.area() > 0 else _reverse_loop(rim)
        outlines[fid] = rim
        centers[fid] = (_centroid(rim), 0.0)
    out.sort(key=lambda lf: -abs(lf[0].area()))
    outer = out[0][0]
    if outer.area() < 0:
        outer = _reverse_loop(outer)
    holes = [_reverse_loop(h) if h.area() > 0 else h for h, _ in out[1:]]
    hole_fid = [f for _, f in out[1:]]
    return outer, holes, centers, outlines, hole_fid, slits, orient


def _projected_rim(ctx: Context, devs, place, host: int, rim: tuple[int, ...]) -> Loop2D:
    """A closed wire of edges, in wire order, projected onto the plane ``host`` in the flat."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve

    pieces: list[list] = []
    for ei in rim:
        e = ctx.edges[ei][0]
        c = BRepAdaptor_Curve(e)
        if c.GetType() == GeomAbs_Circle and abs(c.LastParameter() - c.FirstParameter() - 2 * math.pi) < 1e-9:
            circ = c.Circle()
            centre = _on_host_point(ctx, devs, place, host, circ.Location())
            pieces.append([Arc(centre, circ.Radius(), 0.0, 2 * math.pi, True, ei)])
            continue
        segs = [_with(s, source_edge=ei) for s in _on_host(ctx, devs, place, host, e)]
        if segs:
            pieces.append(segs)
    if not pieces:
        raise refuse(Code.OUTLINE_OPEN, "the rim of an engraving is empty in the flat", edges=tuple(rim))
    # the edges are in wire order, each running either way: turn each to start where the last ended
    if len(pieces) > 1:
        first, second = pieces[0], pieces[1]
        ends = (second[0].start, second[-1].end)
        if min(math.dist(first[0].start, p) for p in ends) < min(math.dist(first[-1].end, p) for p in ends):
            pieces[0] = [s.reversed() for s in reversed(first)]
        for n in range(1, len(pieces)):
            prev, cur = pieces[n - 1][-1].end, pieces[n]
            if math.dist(cur[-1].end, prev) < math.dist(cur[0].start, prev):
                pieces[n] = [s.reversed() for s in reversed(cur)]
    return Loop2D(tuple(s for piece in pieces for s in piece))


def _on_host_point(ctx: Context, devs, place, host: int, p) -> Point:
    """A 3D point projected onto the plane ``host``, in the flat."""
    pos = ctx.faces[host - 1].ad.Plane().Position()
    o, xd, yd = pos.Location(), pos.XDirection(), pos.YDirection()
    d = (p.X() - o.X(), p.Y() - o.Y(), p.Z() - o.Z())
    uv = (d[0] * xd.X() + d[1] * xd.Y() + d[2] * xd.Z(), d[0] * yd.X() + d[1] * yd.Y() + d[2] * yd.Z())
    return transform(place[host], devs[host](*uv))


def _on_host(ctx: Context, devs, place, host: int, e) -> list:
    """An edge of a formed region projected onto its host plane, in the flat."""
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.GCPnts import GCPnts_QuasiUniformDeflection

    c = BRepAdaptor_Curve(e)
    sampler = GCPnts_QuasiUniformDeflection(c, 1e-4)
    params = [sampler.Parameter(i) for i in range(1, sampler.NbPoints() + 1)] if sampler.IsDone() else []
    if len(params) < 2:
        params = [c.FirstParameter(), c.LastParameter()]
    return _fit([_on_host_point(ctx, devs, place, host, c.Value(u)) for u in params])


def _fit(pts: list[Point], tol: float = 1e-7) -> list:
    """A projected edge as one Line or Arc if it is one (a cut face normal to the sheet), else a polyline."""
    if len(pts) < 2 or math.dist(pts[0], pts[-1]) < 1e-12 and len(pts) < 3:
        return []
    a, m, b = pts[0], pts[len(pts) // 2], pts[-1]
    ab = math.dist(a, b)
    if ab > 1e-12 and all(abs(cross(sub(b, a), sub(p, a))) / ab <= tol for p in pts):
        return [Line(a, b)]
    den = 2 * (a[0] * (m[1] - b[1]) + m[0] * (b[1] - a[1]) + b[0] * (a[1] - m[1]))
    if abs(den) > 1e-12:
        ux = ((a[0] ** 2 + a[1] ** 2) * (m[1] - b[1]) + (m[0] ** 2 + m[1] ** 2) * (b[1] - a[1])
              + (b[0] ** 2 + b[1] ** 2) * (a[1] - m[1])) / den  # fmt: skip
        uy = ((a[0] ** 2 + a[1] ** 2) * (b[0] - m[0]) + (m[0] ** 2 + m[1] ** 2) * (a[0] - b[0])
              + (b[0] ** 2 + b[1] ** 2) * (m[0] - a[0])) / den  # fmt: skip
        c = (ux, uy)
        r = math.dist(c, a)
        if all(abs(math.dist(c, p) - r) <= tol for p in pts):
            a0 = heading(c, a)
            angles = [heading(c, p) for p in pts]
            sweep = 0.0
            for u, v in zip(angles, angles[1:], strict=False):
                sweep += (v - u + math.pi) % (2 * math.pi) - math.pi
            return [Arc(c, r, a0, a0 + sweep, sweep > 0)]
    return [Line(p, q) for p, q in zip(pts, pts[1:], strict=False) if math.dist(p, q) > 1e-12]


def _slots(rim: Loop2D, runs: list[list], free: list[list]) -> list[Loop2D]:
    """The slits of a formed region, with their width: between the sheet's cut edge (a run of
    the rim) and the region's own cut edge projected onto the sheet (``free``: per edge, its
    segments).

    Each run is closed by the region's edge that joins its two ends. Where none does, the
    rim and the projected edges are polygonized and the faces along a run are the slits.
    """
    chains = _join([p for p in free if p])
    out = []
    for run in runs:
        a, b = run[0].start, run[-1].end
        best = None
        for ch in chains:
            for segs in (ch, [s.reversed() for s in reversed(ch)]):
                err = max(math.dist(segs[0].start, b), math.dist(segs[-1].end, a))
                if best is None or err < best[0]:
                    best = (err, segs)
        if best is not None and best[0] <= 2e-3:
            segs = list(run) + best[1]
            loop = []
            for prev, cur in zip(segs[-1:] + segs[:-1], segs, strict=True):
                if math.dist(prev.end, cur.start) > 1e-9:  # the region's edge meets the rim within 2e-3
                    loop.append(Line(prev.end, cur.start, cur.source_edge))
                loop.append(cur)
            out.append(_merge_collinear(Loop2D(tuple(loop))))
        else:
            pts = [[sg.start for sg in p] + [p[-1].end] for p in free if p]
            out.extend(_slots_by_faces(rim, [run], pts))
    return out


def _join(pieces: list[list], tol: float = 1e-6) -> list[list]:
    """Segment chains joined end to end where their ends meet."""

    def rev(segs):
        return [s.reversed() for s in reversed(segs)]

    todo = [list(p) for p in pieces]
    chains: list[list] = []
    while todo:
        cur = todo.pop()
        grown = True
        while grown:
            grown = False
            for i, p in enumerate(todo):
                if math.dist(cur[-1].end, p[0].start) <= tol:
                    cur = cur + p
                elif math.dist(cur[-1].end, p[-1].end) <= tol:
                    cur = cur + rev(p)
                elif math.dist(cur[0].start, p[-1].end) <= tol:
                    cur = p + cur
                elif math.dist(cur[0].start, p[0].start) <= tol:
                    cur = rev(p) + cur
                else:
                    continue
                todo.pop(i)
                grown = True
                break
        chains.append(cur)
    return chains


def _slots_by_faces(rim: Loop2D, runs: list[list], free: list[list[Point]]) -> list[Loop2D]:
    """Fallback: the faces cut from the rim by the projected region edges that lie along a run."""
    from shapely import set_precision, snap
    from shapely.geometry import LineString
    from shapely.ops import polygonize, unary_union

    rim_pts = rim.polyline(1e-4)
    rim_line = LineString([*rim_pts, rim_pts[0]])
    # the region's edges end on the rim (a polyline within 1e-4 of it) and on each other
    # (projected from different edges): snap to the rim and to a 1e-6 grid so they meet
    lines = [rim_line] + [snap(LineString(p), rim_line, 2e-4) for p in free if len(p) > 1]
    lines = [set_precision(g, 1e-6) for g in lines]
    cut_lines = []
    for chain in runs:
        pts: list[Point] = []
        for sg in chain:
            sp = sg.sample(1e-4)
            pts.extend(sp if not pts else sp[1:])
        cut_lines.append(LineString(pts))
    out = []
    for face in polygonize(unary_union(lines)):
        if any(face.exterior.intersection(cl.buffer(1e-5)).length > 1e-3 for cl in cut_lines):
            ring = list(face.exterior.coords)[:-1]
            loop = Loop2D(tuple(Line(p, q) for p, q in zip(ring, ring[1:] + ring[:1], strict=True)))
            out.append(_merge_collinear(loop, 1e-9))
    return out


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
    return dot(sub(p, line.p0), sub(line.p1, line.p0))


def _intersect(a: Line, b: Line) -> Point | None:
    d1, d2 = sub(a.p1, a.p0), sub(b.p1, b.p0)
    den = cross(d1, d2)
    if abs(den) <= 1e-9 * (math.hypot(*d1) * math.hypot(*d2) or 1.0):
        return None
    return offset(a.p0, d1, cross(sub(b.p0, a.p0), d2) / den)


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
            out.append(_with(s, p0=transform(T, s.p0), p1=transform(T, s.p1)))
        elif isinstance(s, Arc):
            out.append(_with(s, center=transform(T, s.center)))
        else:
            out.append(_with(s, poles=tuple(transform(T, p) for p in s.poles)))
    return Loop2D(tuple(out))


# ------------------------------------------------------------------------------ bend lines


def _bend_line(f: Face, dev: _Dev, T: Rigid, steps: int = 200) -> list[tuple[Point, Point]]:
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
    return [(transform(T, dev(um, a)), transform(T, dev(um, b))) for a, b in segs if b - a > eps]
