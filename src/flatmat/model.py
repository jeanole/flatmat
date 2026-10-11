"""Public data model. Plain frozen dataclasses: no OCP types cross this boundary."""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Any, Literal

from .diagnostics import Diagnostic

Point = tuple[float, float]


# ---------------------------------------------------------------------------- 2D geometry


@dataclass(frozen=True)
class Line:
    p0: Point
    p1: Point
    source_edge: int = 0

    @property
    def start(self) -> Point:
        return self.p0

    @property
    def end(self) -> Point:
        return self.p1

    def reversed(self) -> Line:
        return replace(self, p0=self.p1, p1=self.p0)

    def sample(self, tol: float = 0.01) -> list[Point]:
        return [self.p0, self.p1]

    def length(self) -> float:
        return math.dist(self.p0, self.p1)


@dataclass(frozen=True)
class Arc:
    """Circular arc from ``a0`` to ``a1`` (radians); counter-clockwise if ``ccw``.

    A full circle has ``|a1 - a0| == 2*pi``.
    """

    center: Point
    radius: float
    a0: float
    a1: float
    ccw: bool = True
    source_edge: int = 0

    def _at(self, a: float) -> Point:
        return (self.center[0] + self.radius * math.cos(a), self.center[1] + self.radius * math.sin(a))

    @property
    def start(self) -> Point:
        return self._at(self.a0)

    @property
    def end(self) -> Point:
        return self._at(self.a1)

    @property
    def sweep(self) -> float:
        return abs(self.a1 - self.a0)

    @property
    def is_circle(self) -> bool:
        return abs(self.sweep - 2 * math.pi) < 1e-9

    def reversed(self) -> Arc:
        return replace(self, a0=self.a1, a1=self.a0, ccw=not self.ccw)

    def sample(self, tol: float = 0.01) -> list[Point]:
        # chord error r(1 - cos(d/2)) <= tol
        step = 2 * math.acos(max(-1.0, 1 - tol / self.radius)) if self.radius > tol else math.pi / 4
        n = max(4, math.ceil(self.sweep / max(step, 1e-6)))
        return [self._at(self.a0 + (self.a1 - self.a0) * i / n) for i in range(n + 1)]

    def length(self) -> float:
        return self.radius * self.sweep


@dataclass(frozen=True)
class Spline:
    """Clamped (possibly rational) B-spline, exactly as developed from the pcurve."""

    degree: int
    poles: tuple[Point, ...]
    knots: tuple[float, ...]
    mults: tuple[int, ...]
    weights: tuple[float, ...] | None = None
    source_edge: int = 0

    @property
    def flat_knots(self) -> list[float]:
        return [k for k, m in zip(self.knots, self.mults, strict=True) for _ in range(m)]

    def point(self, u: float) -> Point:
        """De Boor evaluation (rational if weighted)."""
        p, kv = self.degree, self.flat_knots
        n = len(self.poles)
        u = min(max(u, kv[p]), kv[n])
        k = p
        while k < n - 1 and u >= kv[k + 1]:
            k += 1
        w = self.weights or (1.0,) * n
        d = [(self.poles[j][0] * w[j], self.poles[j][1] * w[j], w[j]) for j in range(k - p, k + 1)]
        for r in range(1, p + 1):
            for j in range(p, r - 1, -1):
                i = j + k - p
                den = kv[i + p - r + 1] - kv[i]
                a = 0.0 if den == 0 else (u - kv[i]) / den
                d[j] = (
                    (1 - a) * d[j - 1][0] + a * d[j][0],
                    (1 - a) * d[j - 1][1] + a * d[j][1],
                    (1 - a) * d[j - 1][2] + a * d[j][2],
                )
        x, y, ww = d[p]
        return (x / ww, y / ww)

    @property
    def start(self) -> Point:
        return self.point(self.knots[0])

    @property
    def end(self) -> Point:
        return self.point(self.knots[-1])

    def reversed(self) -> Spline:
        lo, hi = self.knots[0], self.knots[-1]
        return replace(
            self,
            poles=tuple(reversed(self.poles)),
            knots=tuple(lo + hi - k for k in reversed(self.knots)),
            mults=tuple(reversed(self.mults)),
            weights=tuple(reversed(self.weights)) if self.weights else None,
        )

    def sample(self, tol: float = 0.01) -> list[Point]:
        # control polygon length bounds the curve length; dense enough for checks
        ctrl = sum(math.dist(a, b) for a, b in zip(self.poles, self.poles[1:], strict=False))
        n = max(8 * len(self.poles), min(4000, math.ceil(ctrl / max(tol, 1e-6) ** 0.5 * 2)))
        lo, hi = self.knots[0], self.knots[-1]
        return [self.point(lo + (hi - lo) * i / n) for i in range(n + 1)]


Seg2D = Line | Arc | Spline


@dataclass(frozen=True)
class Loop2D:
    segments: tuple[Seg2D, ...]

    def polyline(self, tol: float = 0.01) -> list[Point]:
        pts: list[Point] = []
        for s in self.segments:
            sp = s.sample(tol)
            pts.extend(sp if not pts else sp[1:])
        return pts[:-1] if len(pts) > 1 and math.dist(pts[0], pts[-1]) < 1e-9 else pts

    def area(self) -> float:
        """Signed area (counter-clockwise positive), exact for lines and arcs.

        Green's theorem per segment: 0.5 * integral of (x dy - y dx). Splines are
        integrated on a dense sampling of the exact curve.
        """
        total = 0.0
        for s in self.segments:
            if isinstance(s, Line):
                (x0, y0), (x1, y1) = s.p0, s.p1
                total += x0 * y1 - x1 * y0
            elif isinstance(s, Arc):
                (cx, cy), (xs, ys), (xe, ye) = s.center, s.start, s.end
                total += s.radius**2 * (s.a1 - s.a0) + cx * (ye - ys) - cy * (xe - xs)
            else:
                p = s.sample(1e-5)
                total += sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(p, p[1:], strict=False))
        return 0.5 * total


# ---------------------------------------------------------------------------- recognition


@dataclass(frozen=True)
class Tolerances:
    """Linear in mm, angular in rad."""

    linear: float = 1e-4
    smooth_angle: float = math.radians(0.1)
    thickness_rel: float = 0.01
    thickness_confidence: float = 0.95
    flat: float = 1e-3
    area_rel: float = 1e-4


@dataclass(frozen=True)
class FaceInfo:
    id: int
    kind: Literal["plane", "cylinder", "other"]
    area: float
    role: Literal["skin_a", "skin_b", "thickness"]


@dataclass(frozen=True)
class BendSegment:
    face_id: int
    length: float
    axis_point: tuple[float, float, float]
    axis_dir: tuple[float, float, float]


@dataclass(frozen=True)
class Bend:
    id: int
    segments: tuple[BendSegment, ...]
    r_in: float
    angle_deg: float
    direction: Literal["up", "down"]
    kind: Literal["bend", "hem"]
    flanges: tuple[int, ...] = ()
    # filled by unfold()
    k: float | None = None
    ba: float | None = None
    bd: float | None = None
    ossb: float | None = None

    @property
    def opening_angle_deg(self) -> float:
        return 180.0 - self.angle_deg


@dataclass(frozen=True)
class Flange:
    id: int
    face_ids: tuple[int, ...]
    area: float


@dataclass(frozen=True)
class SheetPart:
    thickness: float
    thickness_confidence: float
    faces: tuple[FaceInfo, ...]
    bends: tuple[Bend, ...]
    flanges: tuple[Flange, ...]
    reference_root: int
    diagnostics: tuple[Diagnostic, ...] = ()
    source: str = ""
    # internal: the recognition context unfold() works on (OCP objects); not public API
    _ctx: Any = field(default=None, repr=False, compare=False)

    @property
    def occ_shape(self):
        """Escape hatch: the OCP solid this was recognised from."""
        return getattr(self._ctx, "solid", None)


# ---------------------------------------------------------------------------- flat pattern


@dataclass(frozen=True)
class BendLine:
    bend_id: int
    segment_face: int
    p0: Point
    p1: Point
    direction: Literal["up", "down"]
    angle_deg: float
    r_in: float
    kind: Literal["bend", "hem"] = "bend"


@dataclass(frozen=True)
class FlatPattern:
    outer: Loop2D
    holes: tuple[Loop2D, ...]
    bend_lines: tuple[BendLine, ...]
    bends: tuple[Bend, ...]
    thickness: float
    k_model: str
    area: float
    diagnostics: tuple[Diagnostic, ...] = ()
    meta: dict = field(default_factory=dict)

    @property
    def bbox(self) -> tuple[float, float, float, float]:
        pts = self.outer.polyline(0.01)
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        return (min(xs), min(ys), max(xs), max(ys))

    @property
    def size(self) -> tuple[float, float]:
        x0, y0, x1, y1 = self.bbox
        return (x1 - x0, y1 - y0)
