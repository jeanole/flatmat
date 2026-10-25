"""Tool collisions and a bend sequence (docs/design.md §7.4).

Each bend is simulated in the press brake: the part is folded to the state "the bends
before it at their full angle, the rest flat", the bend itself is swept from 0 to its
angle in steps, and at each step every flange is tested against the punch (with its
clamp) and the die (with the lower beam). The tools are prisms along the bend line:
a flange collides if its points within the tool's length, seen along the bend line,
fall inside the tool's section. The part may lie either way round (rotated 180° about
the vertical), which matters for a gooseneck.

The order is searched backwards from the finished part (Duflou et al. 1999): the last
bend must be makeable with all others done; then the one before, and so on. At each
step the candidates are tried in order of three weighted heuristics: shape-defining
bends late (4), bends at the leaves of the flange tree early (2), short bends early
(1). A dead end backtracks, up to a node limit.

What is not modelled: the bend zones themselves (only flanges), formed features
standing out of the sheet, springback, the back gauge, and the operator's hands.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import shapely
from shapely.geometry import LineString, Polygon

from .fold import FoldModel, Hinge, apply, hinge_frame
from .model import FlatPattern
from .tooling import Tooling

_BIG = 1e6


@dataclass(frozen=True)
class Collision:
    """Why bend ``bend_id`` cannot be made in a given state."""

    bend_id: int
    tool: Literal["punch", "clamp", "die", "beam"]  # clamp: upper beam; beam: lower beam
    angle_deg: float  # how far the bend had got
    at: tuple[float, float]  # in the flat: the flange that hits the tool
    done: tuple[int, ...]  # the bends made before it

    def describe(self) -> str:
        before = f" after bends {', '.join(map(str, self.done))}" if self.done else ""
        what = {"clamp": "punch clamp", "beam": "lower beam"}.get(self.tool, self.tool)
        return f"bend {self.bend_id}{before}: the part hits the {what} at {self.angle_deg:.0f}°"


@dataclass(frozen=True)
class BendStep:
    bend_id: int
    flip: bool  # the part is turned over (the other skin up) before this bend
    turned: bool  # the part lies the other way round (rotated 180° about the vertical)


@dataclass(frozen=True)
class BendPlan:
    """A bend order, or the reason there is none.

    ``feasible`` is True when every bend in ``steps`` can be made in that order without
    a collision. When no order is found, ``steps`` is the longest tail found and
    ``collisions`` say what blocks the bends that could not be placed.
    """

    steps: tuple[BendStep, ...]
    feasible: bool
    collisions: tuple[Collision, ...] = ()
    unchecked: tuple[int, ...] = ()  # hems and lances (made first, with their own tools) and bends not placed
    tooling: Tooling | None = None
    nodes: int = 0

    @property
    def order(self) -> tuple[int, ...]:
        return tuple(s.bend_id for s in self.steps)

    @property
    def flips(self) -> int:
        return sum(s.flip for s in self.steps)

    def to_dict(self) -> dict:
        return {
            "feasible": self.feasible,
            "order": list(self.order),
            "steps": [{"bend": s.bend_id, "flip": s.flip, "turned": s.turned} for s in self.steps],
            "collisions": [
                {"bend": c.bend_id, "tool": c.tool, "angle_deg": c.angle_deg, "at": list(c.at), "done": list(c.done)}
                for c in self.collisions
            ],
            "unchecked": list(self.unchecked),
            "tooling": self.tooling.to_dict() if self.tooling else None,
        }


@dataclass
class _Sim:
    model: FoldModel
    tooling: Tooling
    step: float
    clearance: float
    exteriors: dict[int, Polygon] = field(default_factory=dict)
    points: dict[int, np.ndarray] = field(default_factory=dict)
    memo: dict = field(default_factory=dict)
    punch_cache: dict = field(default_factory=dict)

    def punch(self, r_in: float) -> tuple[Polygon, Polygon]:
        rp = min(self.tooling.punch.tip_radius, r_in)
        key = round(rp, 6)
        if key not in self.punch_cache:
            self.punch_cache[key] = self.tooling.punch.sections(rp)
        return self.punch_cache[key]

    def attempt(self, bend: int, done: frozenset[int]) -> tuple[Collision | None, bool]:
        """Can ``bend`` be made with ``done`` already bent? -> (first collision or None, turned)."""
        key = (bend, done)
        if key not in self.memo:
            first = None
            for turned in (False, True):
                c = self._sweep(bend, done, turned)
                if c is None:
                    self.memo[key] = (None, turned)
                    break
                first = first or c
            else:
                self.memo[key] = (first, False)
        return self.memo[key]

    def _sweep(self, bend: int, done: frozenset[int], turned: bool) -> Collision | None:
        m = self.model
        hs = m.hinges_of(bend)
        h = hs[0]
        full = {x.bend_id: x.angle for x in m.hinges if x.bend_id in done}
        n = max(1, math.ceil(math.degrees(h.angle) / self.step))
        # the full angle first: most collisions show there, and a hit ends the sweep
        for i in [n, *range(1, n)]:
            phi = h.angle * i / n
            angles = dict(full)
            angles[bend] = phi
            c = self._state(h, hs, angles, phi, turned)
            if c is not None:
                return Collision(bend, c[0], math.degrees(phi), c[1], tuple(sorted(done)))
        return None

    def _state(self, h: Hinge, hs: list[Hinge], angles: dict[int, float], phi: float, turned: bool):
        m = self.model
        T = m.transforms(angles)
        F, r_in = hinge_frame(h, phi, m.t)
        M = F @ np.linalg.inv(T[h.parent])  # base frame -> machine
        if turned:
            M = np.diag([-1.0, -1.0, 1.0, 1.0]) @ M
        # the bend's extent along the machine x axis
        xs: list[float] = []
        for x in hs:
            Tx = T[x.parent]
            xs.extend(apply(M @ Tx, np.array([x.p0, x.p1]))[:, 0])
        x0, x1 = min(xs), max(xs)
        tools = self.tooling
        xm = (x0 + x1) / 2
        pl = tools.punch.length
        px = (x0, x1) if pl is None else (xm - pl / 2, xm + pl / 2)
        dl = tools.die.length
        dx = (x0, x1) if dl is None else (xm - dl / 2, xm + dl / 2)
        everywhere = (-_BIG, _BIG)
        punch, clamp = self.punch(h.r_in)
        die, beam = self._die(phi, r_in, m.t)
        rad = max(m.t / 2 - self.clearance, 1e-3)
        checks = (("punch", punch, px), ("die", die, dx), ("clamp", clamp, everywhere), ("beam", beam, everywhere))
        for r, Tr in T.items():
            A = M @ Tr
            pts = self.points[r] @ A[:3, :2].T + A[:3, 3]
            lo3, hi3 = pts.min(axis=0) - rad, pts.max(axis=0) + rad
            for tool, poly, (lo, hi) in checks:
                b = poly.bounds  # (ymin, zmin, ymax, zmax)
                if hi3[0] < lo or lo3[0] > hi or hi3[1] < b[0] or lo3[1] > b[2] or hi3[2] < b[1] or lo3[2] > b[3]:
                    continue
                img = _section(self.exteriors[r], A, lo, hi)
                if img is None:
                    continue
                geom, flat = img
                g = geom.buffer(rad, quad_segs=4)
                if g.intersects(poly):
                    c = flat.representative_point()
                    return tool, (round(c.x, 3), round(c.y, 3))
        return None

    def _die(self, phi: float, r_in: float, t: float) -> tuple[Polygon, Polygon]:
        """The die placed so the outside of both legs rests on its shoulders."""
        die = self.tooling.die
        sy, sz = die.shoulder_center()
        rs = die.shoulder_radius
        half = phi / 2
        ro = r_in + t  # outside radius, centre at (0, r_in)
        # the leg's outside line is tangent to the shoulder circle (sy, top + sz)
        top = r_in + (sy * math.sin(half) - ro - rs) / math.cos(half) - sz if half < math.pi / 2 - 1e-6 else -_BIG
        block, beam = die.sections()
        return shapely.affinity.translate(block, 0.0, top), shapely.affinity.translate(beam, 0.0, top)


def _section(poly: Polygon, A: np.ndarray, lo: float, hi: float):
    """The part of a flat region within machine x in [lo, hi], seen along x -> (yz geometry, flat piece)."""
    g = A[0, :2]
    c0 = A[0, 3]
    ng = math.hypot(g[0], g[1])
    if ng < 1e-12:
        if not lo - 1e-9 <= c0 <= hi + 1e-9:
            return None
        piece = poly
    else:
        d = g / ng
        q = np.array([-d[1], d[0]])
        s0, s1 = (max(lo, -_BIG) - c0) / ng, (min(hi, _BIG) - c0) / ng
        L = _BIG
        strip = Polygon([s0 * d - L * q, s1 * d - L * q, s1 * d + L * q, s0 * d + L * q])
        piece = poly.intersection(strip)
        if piece.is_empty:
            return None
    geoms = [p for p in getattr(piece, "geoms", [piece]) if isinstance(p, Polygon) and not p.is_empty]
    if not geoms:
        return None
    parts = []
    for p in geoms:
        pts = np.asarray(p.exterior.coords)[:, :2]
        yz = pts @ A[1:3, :2].T + A[1:3, 3]
        img = Polygon(yz)
        parts.append(img if img.is_valid and img.area > 1e-9 else LineString(yz))
    geom = parts[0] if len(parts) == 1 else shapely.union_all(parts)
    flat = geoms[0] if len(geoms) == 1 else shapely.union_all(geoms)
    return geom, flat


def plan_bends(
    pattern: FlatPattern,
    tooling: Tooling | None = None,
    *,
    order: list[int] | tuple[int, ...] | None = None,
    step_deg: float = 5.0,
    clearance: float | None = None,
    max_nodes: int = 2000,
    model: FoldModel | None = None,
) -> BendPlan:
    """Find a collision-free bend order, or check a given ``order``.

    ``tooling`` defaults to :meth:`Tooling.default` for the sheet's thickness. Each bend
    is swept in ``step_deg`` steps; ``clearance`` (default 2 % of t, at least 0.02 mm)
    is how close a flange may come to a tool.
    """
    t = pattern.thickness
    tooling = tooling or Tooling.default(t)
    model = model or FoldModel(pattern)
    sim = _Sim(model, tooling, step_deg, clearance if clearance is not None else max(0.02, 0.02 * t))
    sim.exteriors = {i: Polygon(r.exterior.coords) for i, r in enumerate(model.regions)}
    sim.points = {i: np.asarray(p.exterior.coords)[:, :2] for i, p in sim.exteriors.items()}
    bends = {b.id: b for b in pattern.bends}
    # hems (pre-bent, closed with a flattening tool) and lances (formed with the blank, in a
    # punching tool) are not made on the press brake: not simulated, made first
    lances = {f.bend_id for f in pattern.features if f.kind == "lance"}
    hems = frozenset(b for b in model.bend_ids() if bends[b].kind == "hem" or b in lances)
    ids = [b for b in model.bend_ids() if b not in hems]
    unchecked = tuple(sorted(set(bends) - set(ids)))
    direction = {b: bends[b].direction for b in ids}

    def steps_for(seq: list[int], turned: dict[int, bool]) -> tuple[BendStep, ...]:
        out = []
        for i, b in enumerate(seq):
            flip = i > 0 and direction[b] != direction[seq[i - 1]]
            out.append(BendStep(b, flip, turned.get(b, False)))
        return tuple(out)

    if order is not None:
        seq = [b for b in order if b in ids]
        cols, turned = [], {}
        for i, b in enumerate(seq):
            c, tn = sim.attempt(b, frozenset(seq[:i]) | hems)
            turned[b] = tn
            if c is not None:
                cols.append(c)
        return BendPlan(steps_for(seq, turned), not cols, tuple(cols), unchecked, tooling, len(sim.memo))

    score = _scores(model, pattern, ids)
    best: list[int] = []  # longest tail found (in reverse)
    blocked: dict[int, Collision] = {}
    dead: set[frozenset[int]] = set()  # sets of bends that cannot all be made first, whatever the tail
    nodes = 0

    def search(remaining: frozenset[int], tail: list[int], turned: dict[int, bool]) -> list[int] | None:
        nonlocal nodes, best
        if not remaining:
            return tail
        if remaining in dead:
            return None
        if len(tail) > len(best):
            best = list(tail)
        for b in sorted(remaining, key=lambda b: (-score[b], b)):
            if nodes >= max_nodes:
                return None
            nodes += 1
            c, tn = sim.attempt(b, (remaining - {b}) | hems)
            if c is not None:
                if b not in blocked or len(c.done) < len(blocked[b].done):
                    blocked[b] = c
                continue
            turned[b] = tn
            got = search(remaining - {b}, tail + [b], turned)
            if got is not None:
                return got
        if nodes < max_nodes:
            dead.add(remaining)
        return None

    turned_found: dict[int, bool] = {}
    found = search(frozenset(ids), [], turned_found)
    if found is not None:
        seq = list(reversed(found))
        return BendPlan(steps_for(seq, turned_found), True, (), unchecked, tooling, nodes)
    tail = list(reversed(best))
    stuck = [blocked[b] for b in ids if b not in tail and b in blocked]
    return BendPlan(steps_for(tail, turned_found), False, tuple(stuck), unchecked, tooling, nodes)


def _scores(model: FoldModel, pattern: FlatPattern, ids: list[int]) -> dict[int, float]:
    """Lateness of each bend (higher: later), Duflou's heuristics with weights 4, 2, 1."""
    area = {i: r.area for i, r in enumerate(model.regions)}
    total = sum(area.values()) or 1.0
    hinge_count: dict[int, int] = {}
    for h in model.hinges:
        hinge_count[h.parent] = hinge_count.get(h.parent, 0) + 1
        hinge_count[h.child] = hinge_count.get(h.child, 0) + 1
    length = {b: 0.0 for b in ids}
    for bl in pattern.bend_lines:
        if bl.bend_id in length:
            length[bl.bend_id] += math.dist(bl.p0, bl.p1)
    longest = max(length.values(), default=1.0) or 1.0
    out = {}
    for b in ids:
        hs = model.hinges_of(b)
        moved = sum(area[r] for h in hs for r in model.subtree(h.child))
        shape = 2 * min(moved, total - moved) / total  # 1: splits the part in half
        leaf = 0.0 if any(hinge_count.get(h.child, 0) == 1 for h in hs) else 1.0
        out[b] = 4 * shape + 2 * leaf + length[b] / longest
    return out


__all__ = ["BendPlan", "BendStep", "Collision", "plan_bends"]
