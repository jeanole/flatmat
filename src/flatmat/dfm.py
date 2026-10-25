"""Manufacturability checks for laser cutting + press-brake air bending (docs/design.md §7).

Every rule works on the flat pattern (2D, shapely) plus the bend data, and every
threshold is a parameter: published shop rules disagree (see docs/research.md §10),
so the defaults below are a documented middle, not the truth.

    report = fm.check(pattern, material="steel", die_v=None)
    for d in report.findings: print(d.code, d.message, d.data["at"])
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import cached_property

from shapely import STRtree
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import unary_union

from .diagnostics import Code, Diagnostic, Severity
from .model import Arc, Bend, BendLine, FlatPattern

# ------------------------------------------------------------------------------- inputs


@dataclass(frozen=True)
class Material:
    """Material properties the rules need. Values are typical, not guaranteed."""

    name: str
    rm_mpa: float  # tensile strength, for tonnage
    air_radius_pct: float  # air-bent inside radius as a fraction of V ("20 % rule")
    min_r_over_t: float  # minimum inner radius / t, bend line across the rolling direction
    min_r_over_t_parallel: float  # ... bend line along the rolling direction
    hem_allowed: bool = True
    hem_max_t: float | None = None  # mm


# Sources: The Fabricator / Durma (air radius % of V), Protocase (hem thickness limits),
# forum/vendor reports for 6061-T6. See research.md §10. Minimum radii: Xometry quotes
# 0.8t (steel) and 2t (stainless) for 1-6 mm; in practice cold-rolled mild steel is
# bent at 0.5t and 304 at ~1t (CAD defaults, e.g. SolidWorks r = 1 mm at t = 2 mm), so
# the defaults below follow practice; the "along the grain" values only apply when the
# rolling direction is given.
MATERIALS: dict[str, Material] = {
    "steel": Material("steel", 400, 0.16, 0.5, 1.0, True, 1.9),
    "galvanised": Material("galvanised", 400, 0.16, 0.5, 1.0, True, 2.0),
    "stainless": Material("stainless", 600, 0.21, 1.0, 1.5, True, 1.6),
    "aluminium_5052": Material("aluminium_5052", 230, 0.12, 1.0, 1.5, True, 1.3),
    "aluminium_6061_t6": Material("aluminium_6061_t6", 310, 0.12, 3.0, 4.0, False, None),
}


@dataclass(frozen=True)
class DfmRules:
    min_flange_v: float = 0.7  # min flange >= this * V ...
    min_flange_t: float = 2.5  # ... and >= this * t + r (outer mould line)
    flange_fraction: float = 0.5  # share of the bend length that must meet the minimum
    hole_to_bend_t: float = 2.5  # hole edge >= this * t + r from the bend tangent line ...
    hole_to_bend_v: float = 0.5  # ... and >= this * V from the bend centre line
    min_hole_t: float = 1.0  # hole diameter >= this * t (laser: 0.5 is common)
    hole_to_edge_t: float = 2.0
    hole_to_hole_t: float = 2.0
    min_feature_t: float = 1.5  # narrowest web / slot of the outline
    relief_width_t: float = 0.5  # relief width >= this * t (SendCutSend; Fabcon says 1.0); depth >= r + t
    relief_fill: float = 0.5  # report when material fills more than this share of the band beyond a bend end
    hem_return_t: float = 4.0
    max_bend_angle: float | None = None  # e.g. 130 for shops that form in one hit
    u_channel_ratio: float = 2.0  # base >= ratio * flange
    radius_tolerance: float = 0.25  # model r vs air-bend r (relative), for the info warning
    disabled: frozenset[str] = frozenset()  # rule codes to skip, e.g. {"DFM_U_CHANNEL"}


def die_opening(t: float) -> float:
    """V from the Amada catalogue rule: 6t (t 0.5-2.5), 8t (3-8), 10t (9-10), 12t (>= 12).

    The catalogue's bands have gaps; the switch points lie in the middle of them, so a
    measured 2.9999 mm sheet gets 8t like a nominal 3 mm one.
    """
    if t < 2.75:
        return 6.0 * t
    if t < 8.5:
        return 8.0 * t
    if t < 11.0:
        return 10.0 * t
    return 12.0 * t


def bending_force_kn(rm_mpa: float, t: float, length: float, v: float) -> float:
    """Air bending force F = 1.42·Rm·t²·L / (1000·V) in kN (mm, MPa)."""
    return 1.42 * rm_mpa * t * t * length / (1000.0 * v)


# ------------------------------------------------------------------------------- report


@dataclass(frozen=True)
class BendCheck:
    bend_id: int
    die_v: float
    air_radius: float
    force_kn: float
    min_flange: float  # required
    flanges: tuple[float, float]  # measured, outer mould line, each side (median along the bend)


@dataclass(frozen=True)
class DfmReport:
    material: Material
    die_v: float
    findings: tuple[Diagnostic, ...]
    bends: tuple[BendCheck, ...] = ()
    rules: DfmRules = field(default_factory=DfmRules)

    @property
    def ok(self) -> bool:
        """No warnings or errors (informational findings do not count)."""
        return not self.problems

    @property
    def problems(self) -> tuple[Diagnostic, ...]:
        return tuple(d for d in self.findings if d.severity != Severity.INFO)

    def codes(self) -> set[str]:
        return {d.code.value for d in self.findings}

    def to_dict(self) -> dict:
        return {
            "material": self.material.name,
            "die_v": self.die_v,
            "findings": [d.to_dict() for d in self.findings],
            "bends": [
                {
                    "id": b.bend_id,
                    "die_v": b.die_v,
                    "air_radius": b.air_radius,
                    "force_kn": b.force_kn,
                    "min_flange": b.min_flange,
                    "flanges": list(b.flanges),
                }
                for b in self.bends
            ],
        }


# ------------------------------------------------------------------------------- helpers


def _finding(code: Code, msg: str, sev: Severity, at=None, bend: int | None = None, **data) -> Diagnostic:
    d = dict(data)
    if at is not None:
        d["at"] = [round(at[0], 4), round(at[1], 4)]
    if bend is not None:
        d["bend"] = bend
    return Diagnostic(code, sev, msg, (), (), d)


def _setback(b: Bend, t: float) -> float:
    return b.ossb if b.ossb is not None else b.r_in + t


def _unit(p0, p1):
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    n = math.hypot(dx, dy)
    return (dx / n, dy / n), n


# ------------------------------------------------------------------------------- check


@dataclass
class _Sheet:
    """The flat pattern as the rules see it: shapely geometry, bends and their lines."""

    pattern: FlatPattern
    t: float
    v: float  # die opening
    outer: Polygon
    holes: list[Polygon]
    part: Polygon  # outer minus holes
    bends: dict[int, Bend]
    lines_by_bend: dict[int, list[BendLine]]
    all_lines: list[tuple[BendLine, LineString]]

    @cached_property
    def line_tree(self) -> STRtree:
        return STRtree([ls for _, ls in self.all_lines])

    @cached_property
    def hole_tree(self) -> STRtree:
        return STRtree(self.holes)

    @cached_property
    def ends(self) -> list[tuple[tuple[float, float], BendLine, bool]]:
        """Bend-line ends, and whether each sits on another bend (a corner where two bends meet)."""
        reach = max(((b.ba or 0.0) / 2 + self.t for b in self.bends.values()), default=0.0)
        out = []
        for bl, _ls in self.all_lines:
            for end in (bl.p0, bl.p1):
                corner = any(
                    other is not bl and ls2.distance(Point(end)) <= (self.bends[other.bend_id].ba or 0.0) / 2 + self.t
                    for other, ls2 in (self.all_lines[j] for j in _near(self.line_tree, (*end, *end), reach))
                )
                out.append((end, bl, corner))
        return out

    @cached_property
    def end_tree(self) -> STRtree:
        return STRtree([Point(e) for e, _, _ in self.ends])


def _near(tree: STRtree, bounds, margin: float) -> list[int]:
    """Indices (ascending) of the tree's geometries that can come within margin of the box ``bounds``.

    A box query: the exact test is the caller's. Distances never undercut the distance
    between bounding boxes, so nothing within the margin is missed.
    """
    x0, y0, x1, y1 = bounds
    m = margin + 1e-6
    return sorted(int(j) for j in tree.query(box(x0 - m, y0 - m, x1 + m, y1 + m)))


class _Findings:
    """Collects findings; a rule that is switched off adds none, and can ask to skip its work."""

    def __init__(self, rules: DfmRules, sev: Severity):
        self.rules, self.sev = rules, sev
        self.items: list[Diagnostic] = []

    def on(self, *codes: Code) -> bool:
        return any(c.value not in self.rules.disabled for c in codes)

    def add(self, code: Code, msg: str, at=None, bend=None, info: bool = False, **data) -> None:
        if self.on(code):
            self.items.append(_finding(code, msg, Severity.INFO if info else self.sev, at, bend, **data))


def check(
    pattern: FlatPattern,
    *,
    material: str | Material = "steel",
    die_v: float | None = None,
    rolling_direction_deg: float | None = None,
    machine_max_force_kn: float | None = None,
    rules: DfmRules | None = None,
    strict: bool = False,
) -> DfmReport:
    """Check a flat pattern against laser and press-brake rules.

    ``rolling_direction_deg`` is the rolling direction in the flat pattern's frame
    (0 = +X). Bends within 45° of it must meet the larger "along the grain" minimum
    radius; without it, every bend is held to the normal minimum.
    With ``strict`` the findings are errors instead of warnings.
    """
    rules = rules or DfmRules()
    if isinstance(material, str):
        if material not in MATERIALS:
            raise ValueError(f"unknown material {material!r}; known: {', '.join(sorted(MATERIALS))}")
        mat = MATERIALS[material]
    else:
        mat = material
    t = pattern.thickness
    v = die_v or die_opening(t)
    out = _Findings(rules, Severity.ERROR if strict else Severity.WARNING)

    outer = Polygon(pattern.outer.polyline(0.01))
    holes = [Polygon(h.polyline(0.01)) for h in pattern.holes]
    lines_by_bend: dict[int, list] = {}
    for bl in pattern.bend_lines:
        lines_by_bend.setdefault(bl.bend_id, []).append(bl)
    sheet = _Sheet(
        pattern,
        t,
        v,
        outer,
        holes,
        Polygon(outer.exterior.coords, [h.exterior.coords for h in holes]),
        {b.id: b for b in pattern.bends},
        lines_by_bend,
        [(bl, LineString([bl.p0, bl.p1])) for bl in pattern.bend_lines],
    )

    checks = _bend_rules(sheet, mat, rules, out, rolling_direction_deg, machine_max_force_kn)
    if out.on(Code.DFM_U_CHANNEL):
        _u_channels(pattern, sheet.bends, checks, rules, out.add)
    flagged = _hole_rules(sheet, rules, out)
    if out.on(Code.DFM_FEATURE_TOO_SMALL):
        _web_and_slot_rules(sheet, rules, out, flagged)
    if out.on(Code.DFM_RELIEF_MISSING):
        _relief_rules(sheet, rules, out)
    return DfmReport(mat, v, tuple(out.items), tuple(checks), rules)


def _bend_rules(
    sheet: _Sheet,
    mat: Material,
    rules: DfmRules,
    out: _Findings,
    rolling_direction_deg: float | None,
    machine_max_force_kn: float | None,
) -> list[BendCheck]:
    """Per bend: flanges, hem, radius, angle and tonnage; always measured for the report."""
    add = out.add
    t, v, bends = sheet.t, sheet.v, sheet.bends
    checks: list[BendCheck] = []
    for bid, b in bends.items():
        bls = sheet.lines_by_bend.get(bid, [])
        if not bls:
            continue
        ba = b.ba or 0.0
        need = max(rules.min_flange_v * v, rules.min_flange_t * t + b.r_in)
        legs: list[list[float]] = [[], []]
        worst: list[tuple[float, tuple]] = [(math.inf, (0, 0)), (math.inf, (0, 0))]
        for bl in bls:
            u, length = _unit(bl.p0, bl.p1)
            n = (-u[1], u[0])
            for k in range(11):
                s = 0.05 + 0.9 * k / 10
                p = (bl.p0[0] + s * length * u[0], bl.p0[1] + s * length * u[1])
                for side, sg in ((0, 1.0), (1, -1.0)):
                    leg = _leg(p, (sg * n[0], sg * n[1]), sheet, bl, ba, b)
                    legs[side].append(leg)
                    if leg < worst[side][0]:
                        worst[side] = (leg, p)
        med = tuple(sorted(x)[len(x) // 2] for x in legs)
        for side in (0, 1):
            if b.kind == "hem" and med[side] <= med[1 - side]:
                continue  # the return of a hem has its own rule
            ok = sum(1 for x in legs[side] if x >= need - 1e-9) / len(legs[side])
            if ok < rules.flange_fraction:
                add(
                    Code.DFM_FLANGE_TOO_SHORT,
                    f"bend {bid}: flange {med[side]:.2f} mm < {need:.2f} mm (V={v:g})",
                    worst[side][1],
                    bid,
                    flange=med[side],
                    required=need,
                )
        if b.kind == "hem":
            ret = min(med)
            if ret < rules.hem_return_t * t:
                add(Code.DFM_HEM_TOO_SHORT, f"bend {bid}: hem return {ret:.2f} mm < {rules.hem_return_t * t:.2f} mm",
                    bls[0].p0, bid, hem_return=ret)  # fmt: skip
            if not mat.hem_allowed or (mat.hem_max_t is not None and t > mat.hem_max_t):
                add(Code.DFM_HEM_MATERIAL, f"bend {bid}: no hems in {mat.name} at t={t:g} mm", bls[0].p0, bid)
        # radius vs material (grain-aware) and vs what air bending gives
        u, _ = _unit(bls[0].p0, bls[0].p1)
        along = False  # unknown grain: the normal minimum; the along-the-grain one needs the direction
        if rolling_direction_deg is not None:
            rd = (math.cos(math.radians(rolling_direction_deg)), math.sin(math.radians(rolling_direction_deg)))
            along = abs(u[0] * rd[0] + u[1] * rd[1]) > math.cos(math.radians(45))
        r_min = (mat.min_r_over_t_parallel if along else mat.min_r_over_t) * t
        if b.r_in < r_min - 1e-9 and b.kind != "hem":
            add(Code.DFM_RADIUS_TOO_SMALL, f"bend {bid}: r={b.r_in:g} < {r_min:g} mm for {mat.name}"
                + (" (along the grain)" if along else ""), bls[0].p0, bid, r_min=r_min)  # fmt: skip
        r_air = mat.air_radius_pct * v
        if b.kind != "hem" and abs(b.r_in - r_air) > rules.radius_tolerance * max(r_air, t):
            add(Code.DFM_RADIUS_NOT_AIR_BENDABLE, f"bend {bid}: model r={b.r_in:g}, air bending in V{v:g} gives "
                f"~{r_air:.2f} mm; the flat size will differ", bls[0].p0, bid, info=True, air_radius=r_air)  # fmt: skip
        if rules.max_bend_angle is not None and b.kind != "hem" and b.angle_deg > rules.max_bend_angle:
            add(Code.DFM_BEND_ANGLE, f"bend {bid}: {b.angle_deg:.1f}° > {rules.max_bend_angle:g}°", bls[0].p0, bid)
        total_len = sum(math.dist(bl.p0, bl.p1) for bl in bls)
        force = bending_force_kn(mat.rm_mpa, t, total_len, v)
        if machine_max_force_kn is not None and force > machine_max_force_kn:
            add(Code.DFM_TONNAGE, f"bend {bid}: {force:.0f} kN > machine {machine_max_force_kn:g} kN",
                bls[0].p0, bid, force_kn=force)  # fmt: skip
        checks.append(BendCheck(bid, v, r_air, force, need, (med[0], med[1])))
    return checks


def _hole_rules(sheet: _Sheet, rules: DfmRules, out: _Findings) -> set[int]:
    """Holes: size, distance to the edge, to each other and to the bends.

    Returns the holes that sit in a narrow web (near the edge or another hole), which the
    web rule does not report a second time; so those distances are measured whenever
    the web rule is on, even when their own rules are off.
    """
    add = out.add
    t, v, holes, bends = sheet.t, sheet.v, sheet.holes, sheet.bends
    webs = out.on(Code.DFM_FEATURE_TOO_SMALL)
    to_edge = webs or out.on(Code.DFM_HOLE_NEAR_EDGE)
    to_holes = webs or out.on(Code.DFM_HOLES_TOO_CLOSE)
    to_bends = out.on(Code.DFM_CUT_IN_BEND_ZONE, Code.DFM_HOLE_NEAR_BEND)
    tree = sheet.hole_tree if holes and to_holes else None
    # how far from a bend line a hole can still be reported, for the box query below
    reach = max(
        (
            max(rules.hole_to_bend_t * t + b.r_in + (b.ba or 0.0) / 2, rules.hole_to_bend_v * v, (b.ba or 0.0) / 2)
            for b in bends.values()
        ),
        default=0.0,
    )
    flagged: set[int] = set()
    for i, h in enumerate(holes):
        c = h.centroid
        loop = sheet.pattern.holes[i]
        seg0 = loop.segments[0]
        if len(loop.segments) == 1 and isinstance(seg0, Arc) and seg0.is_circle:
            dia = 2 * seg0.radius
            if dia < rules.min_hole_t * t - 1e-9:
                add(Code.DFM_HOLE_TOO_SMALL, f"hole ⌀{dia:g} < {rules.min_hole_t * t:g} mm", (c.x, c.y),
                    diameter=dia)  # fmt: skip
        if to_edge:
            d_edge = h.distance(sheet.outer.exterior)
            if d_edge < rules.hole_to_edge_t * t - 1e-9:
                flagged.add(i)
                add(Code.DFM_HOLE_NEAR_EDGE, f"hole {d_edge:.2f} mm from the edge < {rules.hole_to_edge_t * t:g}",
                    (c.x, c.y), distance=d_edge)  # fmt: skip
        if tree is not None:
            for j in tree.query(h.buffer(rules.hole_to_hole_t * t)):
                j = int(j)
                if j <= i:
                    continue
                d = h.distance(holes[j])
                if d < rules.hole_to_hole_t * t - 1e-9:
                    flagged.update((i, j))
                    add(Code.DFM_HOLES_TOO_CLOSE, f"web between holes {d:.2f} mm < {rules.hole_to_hole_t * t:g}",
                        (c.x, c.y), distance=d)  # fmt: skip
        if not to_bends:
            continue
        in_zone: set[int] = set()
        for j in _near(sheet.line_tree, h.bounds, reach):
            bl, ls = sheet.all_lines[j]
            b = bends[bl.bend_id]
            ba = b.ba or 0.0
            d = h.distance(ls)
            if d <= ba / 2 + 1e-9:
                if b.id not in in_zone:
                    in_zone.add(b.id)
                    add(Code.DFM_CUT_IN_BEND_ZONE, f"hole inside the bend zone of bend {b.id}", (c.x, c.y), b.id)
                continue
            if b.id in in_zone:
                continue
            need = max(rules.hole_to_bend_t * t + b.r_in + ba / 2, rules.hole_to_bend_v * v)
            if d < need - 1e-9 and _beside(h, bl):
                add(Code.DFM_HOLE_NEAR_BEND, f"hole {d - ba / 2:.2f} mm from the tangent line of bend {b.id}; "
                    f"needs {need - ba / 2:.2f}", (c.x, c.y), b.id, distance=d, required=need)  # fmt: skip
    return flagged


def _web_and_slot_rules(sheet: _Sheet, rules: DfmRules, out: _Findings, flagged: set[int]) -> None:
    """Narrow webs (an opening removes them) and narrow slots (a closing fills them)."""
    add = out.add
    t, part, bends = sheet.t, sheet.part, sheet.bends
    rad = rules.min_feature_t * t / 2
    near_flagged = unary_union([sheet.holes[i].buffer(rules.hole_to_edge_t * t) for i in flagged]) if flagged else None
    relief_margin = max((b.ba or 0) / 2 + b.r_in + 2 * t for b in bends.values()) if bends else 0.0
    opened = part.buffer(-rad, join_style="round").buffer(rad, join_style="round")
    closed = part.buffer(rad, join_style="round").buffer(-rad, join_style="round")
    for what, diff in (("web", part.difference(opened)), ("slot", closed.difference(part))):
        for piece in getattr(diff, "geoms", [diff]):
            # a convex (opening) or concave (closing) corner alone changes < ~0.22 r²
            if piece.is_empty or piece.area <= 0.6 * rad * rad:
                continue
            if near_flagged is not None and piece.within(near_flagged):
                continue  # the same narrow web as a hole-to-edge / hole-to-hole finding
            if what == "slot" and any(
                Point(sheet.ends[j][0]).distance(piece) <= relief_margin
                for j in _near(sheet.end_tree, piece.bounds, relief_margin)
            ):
                continue  # a relief: judged by the relief rule, not as a slot
            c = piece.representative_point()
            add(Code.DFM_FEATURE_TOO_SMALL, f"{what} narrower than {2 * rad:.2f} mm"
                + (" (or a very acute corner)" if what == "web" else ""), (c.x, c.y), area=piece.area)  # fmt: skip


def _relief_rules(sheet: _Sheet, rules: DfmRules, out: _Findings) -> None:
    """Reliefs at bend ends (not at corners where two bends meet: those need a corner relief)."""
    add = out.add
    t, part = sheet.t, sheet.part
    corner_ends = {(id(bl), e) for e, bl, corner in sheet.ends if corner}
    for bl, _ls in sheet.all_lines:
        b = sheet.bends[bl.bend_id]
        ba = b.ba or 0.0
        u, _ = _unit(bl.p0, bl.p1)
        n = (-u[1], u[0])
        depth = ba / 2 + b.r_in + t
        width = rules.relief_width_t * t
        for end, sign in ((bl.p1, 1.0), (bl.p0, -1.0)):
            if (id(bl), end) in corner_ends:
                continue
            if any(
                sheet.holes[j].exterior.distance(Point(end)) < 1e-3 for j in _near(sheet.hole_tree, (*end, *end), 1e-3)
            ):
                continue  # the bend line stops at a hole in the bend zone: DFM_CUT_IN_BEND_ZONE
            q0 = (end[0] + sign * u[0] * 1e-3, end[1] + sign * u[1] * 1e-3)
            q1 = (end[0] + sign * u[0] * width, end[1] + sign * u[1] * width)
            # one band per side beyond the end: from the bend line to r + t past the tangent line.
            # A real missing relief fills a band; an outline that merely flares out at an
            # angle from the end touches a thin wedge of it.
            fill = 0.0
            for sg in (1.0, -1.0):
                band = Polygon(
                    [
                        q0,
                        q1,
                        (q1[0] + sg * n[0] * depth, q1[1] + sg * n[1] * depth),
                        (q0[0] + sg * n[0] * depth, q0[1] + sg * n[1] * depth),
                    ]
                )
                fill = max(fill, part.intersection(band).area / band.area)
            if fill > rules.relief_fill:
                add(Code.DFM_RELIEF_MISSING, f"bend {b.id} ends in material: needs a relief ≥ {width:g} mm wide, "
                    f"≥ {b.r_in + t:g} mm deep past the tangent line", end, b.id)  # fmt: skip


def _beside(h: Polygon, bl) -> bool:
    """Is the hole beside the bend line (projects onto it), not off its end?"""
    u, length = _unit(bl.p0, bl.p1)
    minx, miny, maxx, maxy = h.bounds
    corners = [(minx, miny), (minx, maxy), (maxx, miny), (maxx, maxy)]
    s = [(p[0] - bl.p0[0]) * u[0] + (p[1] - bl.p0[1]) * u[1] for p in corners]
    return max(s) > 0 and min(s) < length


def _leg(p, d, sheet: _Sheet, this_line: BendLine, ba: float, b: Bend) -> float:
    """Outer-mould-line length of the flange from bend line point p in direction d.

    Stops at the outline or at the next bend line, whichever comes first.
    """
    far = 1e6
    ray = LineString([p, (p[0] + d[0] * far, p[1] + d[1] * far)])
    best = math.inf
    stop_bend: Bend | None = None
    t, bends = sheet.t, sheet.bends
    hit = ray.intersection(sheet.outer.exterior)
    for g in getattr(hit, "geoms", [hit]):
        if g.is_empty:
            continue
        for q in getattr(g, "coords", []):
            dd = math.dist(p, q)
            if 1e-6 < dd < best:
                best = dd
    for j in sorted(int(j) for j in sheet.line_tree.query(ray, predicate="intersects")):
        bl, ls = sheet.all_lines[j]
        if bl is this_line:
            continue
        x = ray.intersection(ls)
        if x.is_empty:
            continue
        pts = [x] if isinstance(x, Point) else list(getattr(x, "geoms", []))
        for q in pts:
            if not isinstance(q, Point):
                continue
            dd = math.dist(p, (q.x, q.y))
            if 1e-6 < dd < best:
                best, stop_bend = dd, bends[bl.bend_id]
    if best is math.inf:
        return 0.0
    leg = best - ba / 2 + _setback(b, t)
    if stop_bend is not None:
        leg += -(stop_bend.ba or 0.0) / 2 + _setback(stop_bend, t)
    return leg


def _u_channels(pattern: FlatPattern, bends: dict[int, Bend], checks: list[BendCheck], rules: DfmRules, add) -> None:
    by_id = {c.bend_id: c for c in checks}
    lines: dict[int, BendLine] = {}
    for bl in pattern.bend_lines:
        lines.setdefault(bl.bend_id, bl)
    ids = sorted(lines)
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            A, B = bends[a], bends[b]
            if A.direction != B.direction or abs(A.angle_deg - 90) > 15 or abs(B.angle_deg - 90) > 15:
                continue
            la, lb = lines[a], lines[b]
            ua, _ = _unit(la.p0, la.p1)
            ub, _ = _unit(lb.p0, lb.p1)
            if abs(abs(ua[0] * ub[0] + ua[1] * ub[1]) - 1) > 1e-6:
                continue
            n = (-ua[1], ua[0])
            base = abs((lb.p0[0] - la.p0[0]) * n[0] + (lb.p0[1] - la.p0[1]) * n[1])
            # base outer mould line between the two bends
            base_oml = (
                base
                - (A.ba or 0) / 2
                - (B.ba or 0) / 2
                + _setback(A, pattern.thickness)
                + _setback(B, pattern.thickness)
            )
            if a not in by_id or b not in by_id:
                continue

            def leg_away(line, other, check) -> float:
                # flanges[0] is on the +normal side of the line; the leg is the side away from the other bend
                u, _ = _unit(line.p0, line.p1)
                nn = (-u[1], u[0])
                towards = (other.p0[0] - line.p0[0]) * nn[0] + (other.p0[1] - line.p0[1]) * nn[1] > 0
                return check.flanges[1] if towards else check.flanges[0]

            flange = max(leg_away(la, lb, by_id[a]), leg_away(lb, la, by_id[b]))
            if base_oml < rules.u_channel_ratio * flange - 1e-9:
                add(Code.DFM_U_CHANNEL, f"bends {a} and {b} form a U: base {base_oml:.1f} mm < "
                    f"{rules.u_channel_ratio:g} × flange {flange:.1f} mm (gooseneck punch?)",
                    la.p0, a, base=base_oml, flange=flange)  # fmt: skip


__all__ = [
    "MATERIALS",
    "BendCheck",
    "DfmReport",
    "DfmRules",
    "Material",
    "bending_force_kn",
    "check",
    "die_opening",
]
