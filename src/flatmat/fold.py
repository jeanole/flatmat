"""Folding a flat pattern back up (docs/design.md §7.4).

The flat is cut into regions along the bend zones (the BA strips around each bend
line): what is left are the flanges. A bend line joins the region on each side of its
strip; from the largest region (the base) a spanning tree over these hinges fixes
which side of each bend moves. Folding a bend by an angle φ turns everything beyond
its strip rigidly, so a state "these bends at these angles" is one rigid 3D transform
per region.

The regions stand for the sheet's mid-surface (t/2 from either skin), which is exact
whatever the K-factor and whichever skin is inside each bend: the strip ends at the
two ends of a mid-surface arc of radius r + t/2, whose length θ·(r + t/2) is kept
while the bend is part way (radius θ·(r + t/2)/φ). The bend zones themselves are not
drawn. That is
what the corner and collision checks need: flanges, where they are, and how far
apart. Coordinates: the flat's x, y and z = 0 for the base; +z is the viewer's side,
towards which an "up" bend turns.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass

import numpy as np
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

from .geom2d import cut_slits, midpoint, offset, perp, unit
from .model import FlatPattern


@dataclass(frozen=True)
class Hinge:
    """One bend line segment joining region ``parent`` (nearer the base) to ``child``."""

    bend_id: int
    parent: int
    child: int
    a: tuple[float, float]  # on the parent's tangent line, at the middle of the segment
    u: tuple[float, float]  # along the bend line
    n: tuple[float, float]  # from parent to child, in the flat
    ba: float  # strip width in the flat
    angle: float  # full bend angle, rad
    sigma: float  # +1: the child turns towards +z ("up")
    r_in: float
    t: float
    p0: tuple[float, float]  # the bend line segment
    p1: tuple[float, float]


class FoldModel:
    """Regions of a flat pattern and the hinges between them."""

    def __init__(self, pattern: FlatPattern, tol: float = 0.01):
        t = pattern.thickness
        self.t = t
        self.pattern = pattern
        bends = {b.id: b for b in pattern.bends}
        part = Polygon(pattern.outer.polyline(tol), [h.polyline(tol) for h in pattern.holes]).buffer(0)
        part = cut_slits(part, pattern.slits, tol)
        strips = []
        for bl in pattern.bend_lines:
            b = bends[bl.bend_id]
            u, length = unit(bl.p0, bl.p1)
            if length <= 1e-9:
                continue
            ext = t  # past the line's ends, so the strip separates its two sides in the corners too
            half = max((b.ba or 0.0) / 2, 1e-3)
            n = perp(u)
            q0, q1 = offset(bl.p0, u, -ext), offset(bl.p1, u, ext)
            strips.append(Polygon([offset(q, n, s * half) for q, s in ((q0, -1), (q1, -1), (q1, 1), (q0, 1))]))
        rest = part.difference(unary_union(strips)) if strips else part
        pieces = [g for g in getattr(rest, "geoms", [rest]) if isinstance(g, Polygon) and g.area > 1e-6 * part.area]
        self.regions: list[Polygon] = pieces
        self.unseparated: set[int] = set()  # bends whose two sides are one region (no relief?)

        def region_at(p) -> int | None:
            pt = Point(p)
            best = min(range(len(pieces)), key=lambda i: pieces[i].distance(pt), default=None)
            if best is None or pieces[best].distance(pt) > 0.05 + 0.05 * t:
                return None
            return best

        raw = []
        for bl in pattern.bend_lines:
            b = bends[bl.bend_id]
            u, length = unit(bl.p0, bl.p1)
            if length <= 1e-9:
                continue
            n = perp(u)
            m = midpoint(bl.p0, bl.p1)
            d = (b.ba or 0.0) / 2 + 0.02 + 0.01 * t
            ra = region_at(offset(m, n, -d))
            rb = region_at(offset(m, n, d))
            if ra is None or rb is None:
                continue
            if ra == rb:
                self.unseparated.add(b.id)
                continue
            raw.append((b, bl, ra, rb, m, u, n))
        self.base = max(range(len(pieces)), key=lambda i: pieces[i].area) if pieces else 0
        # spanning tree from the base
        by_region: dict[int, list] = {}
        for h in raw:
            by_region.setdefault(h[2], []).append(h)
            by_region.setdefault(h[3], []).append(h)
        self.hinges: list[Hinge] = []
        self.parent: dict[int, int] = {}  # region -> index of the hinge to its parent
        self.depth: dict[int, int] = {self.base: 0}
        q = deque([self.base])
        while q:
            r = q.popleft()
            for b, bl, ra, rb, m, u, n in by_region.get(r, []):
                c = rb if ra == r else ra
                if c in self.depth:
                    continue
                nn = n if ra == r else (-n[0], -n[1])
                ba = b.ba or 0.0
                a = offset(m, nn, -ba / 2)
                self.hinges.append(
                    Hinge(
                        b.id,
                        r,
                        c,
                        a,
                        u,
                        nn,
                        ba,
                        math.radians(b.angle_deg),
                        1.0 if b.direction == "up" else -1.0,
                        b.r_in,
                        t,
                        bl.p0,
                        bl.p1,
                    )  # fmt: skip
                )
                self.parent[c] = len(self.hinges) - 1
                self.depth[c] = self.depth[r] + 1
                q.append(c)
        self.reached = set(self.depth)

    # ------------------------------------------------------------------ structure

    def bend_ids(self) -> list[int]:
        return sorted({h.bend_id for h in self.hinges})

    def hinges_of(self, bend_id: int) -> list[Hinge]:
        return [h for h in self.hinges if h.bend_id == bend_id]

    def subtree(self, region: int) -> set[int]:
        """The region and everything that hangs from it (away from the base)."""
        out = {region}
        grow = True
        while grow:
            grow = False
            for h in self.hinges:
                if h.parent in out and h.child not in out:
                    out.add(h.child)
                    grow = True
        return out

    # ------------------------------------------------------------------ folding

    def transforms(self, angles: dict[int, float]) -> dict[int, np.ndarray]:
        """Region -> 4x4 transform from the flat (z = 0) to the folded base frame.

        ``angles`` maps bend ids to their current angle in rad; missing bends are flat.
        """
        out = {self.base: np.eye(4)}
        order = sorted(self.parent, key=lambda r: self.depth[r])
        for r in order:
            h = self.hinges[self.parent[r]]
            out[r] = out[h.parent] @ fold_matrix(h, angles.get(h.bend_id, 0.0))
        return out

    def folded(self, angles: dict[int, float] | None = None) -> dict[int, np.ndarray]:
        """All bends at their full angle unless ``angles`` says otherwise."""
        full = {h.bend_id: h.angle for h in self.hinges}
        full.update(angles or {})
        return self.transforms(full)


def fold_matrix(h: Hinge, phi: float) -> np.ndarray:
    """Fold the child side of hinge ``h`` by ``phi`` (in the parent's flat frame)."""
    U = np.array([h.u[0], h.u[1], 0.0])
    N = np.array([h.n[0], h.n[1], 0.0])
    Z = np.array([0.0, 0.0, 1.0])
    a = np.array([h.a[0], h.a[1], 0.0])
    if abs(phi) < 1e-12:  # not bent yet: the blank as cut
        return np.eye(4)
    c, s = math.cos(phi), math.sin(phi)
    N2 = c * N + h.sigma * s * Z
    Z2 = -h.sigma * s * N + c * Z
    R = np.outer(U, U) + np.outer(N2, N) + np.outer(Z2, Z)
    rn = h.angle * mid_radius(h) / phi  # the mid-surface arc keeps its length while it bends
    start = a + h.ba * N  # the child's tangent line in the flat
    end = a + rn * s * N + h.sigma * rn * (1 - c) * Z
    M = np.eye(4)
    M[:3, :3] = R
    M[:3, 3] = end - R @ start
    return M


def mid_radius(h: Hinge) -> float:
    return h.r_in + h.t / 2


def hinge_frame(h: Hinge, phi: float, t: float) -> tuple[np.ndarray, float]:
    """Press-brake frame of hinge ``h`` bent to ``phi``, in the parent's flat frame.

    Returns (4x4 from the parent's flat frame to machine coordinates, inner radius):
    x along the bend line, z up through the arc's middle (the inside of the bend faces
    up, towards the punch), origin at the inside surface's lowest point.
    """
    U = np.array([h.u[0], h.u[1], 0.0])
    N = np.array([h.n[0], h.n[1], 0.0])
    Z = np.array([0.0, 0.0, 1.0])
    a = np.array([h.a[0], h.a[1], 0.0])
    phi = max(phi, 1e-6)
    rn = h.angle * mid_radius(h) / phi
    r_in = max(rn - t / 2, 0.0)
    zm = -math.sin(phi / 2) * N + h.sigma * math.cos(phi / 2) * Z
    xm = U
    ym = np.cross(zm, xm)
    center = a + h.sigma * rn * Z
    origin = center - r_in * zm
    M = np.eye(4)
    M[0, :3], M[1, :3], M[2, :3] = xm, ym, zm
    M[:3, 3] = -M[:3, :3] @ origin
    return M, r_in


def apply(M: np.ndarray, pts2: np.ndarray) -> np.ndarray:
    """4x4 transform of flat points (N x 2, z = 0) -> N x 3."""
    return pts2 @ M[:3, :2].T + M[:3, 3]


def boundary_points(poly: Polygon, step: float) -> np.ndarray:
    """Points along a polygon's outline, at most ``step`` apart."""
    ring = LineString(poly.exterior.coords).segmentize(step)
    return np.asarray(ring.coords)[:-1, :2]


__all__ = ["FoldModel", "Hinge", "apply", "boundary_points", "fold_matrix", "hinge_frame"]
