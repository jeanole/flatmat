"""Small 2D helpers for the flat pattern: vectors as (x, y) tuples, rigid transforms, slits.

Each helper does its arithmetic in the order the callers used to write it out, so
results are the same to the last bit.
"""

from __future__ import annotations

import math
from collections.abc import Iterable

from shapely.geometry import LineString

from .model import Loop2D, Point

Vec = tuple[float, float]


def sub(a: Point, b: Point) -> Vec:
    """a - b."""
    return (a[0] - b[0], a[1] - b[1])


def dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1]


def cross(a: Vec, b: Vec) -> float:
    """z of a × b: > 0 if b turns counter-clockwise from a."""
    return a[0] * b[1] - a[1] * b[0]


def perp(u: Vec) -> Vec:
    """u turned 90° counter-clockwise (the left normal of a direction)."""
    return (-u[1], u[0])


def offset(p: Point, d: Vec, s: float) -> Point:
    """p moved by s along d."""
    return (p[0] + d[0] * s, p[1] + d[1] * s)


def midpoint(a: Point, b: Point) -> Point:
    return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)


def unit(p0: Point, p1: Point) -> tuple[Vec, float]:
    """Unit direction from p0 to p1, and the distance (+x for coincident points)."""
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    n = math.hypot(dx, dy)
    return ((dx / n, dy / n) if n > 0 else (1.0, 0.0)), n


def heading(a: Point, b: Point) -> float:
    """Angle of the direction from a to b, in radians from +x."""
    return math.atan2(b[1] - a[1], b[0] - a[0])


# ------------------------------------------------------------------------------ rigid transforms

Rigid = tuple[float, float, float, float]  # cos, sin, tx, ty
IDENTITY: Rigid = (1.0, 0.0, 0.0, 0.0)


def transform(T: Rigid, p: Point) -> Point:
    c, s, tx, ty = T
    return (c * p[0] - s * p[1] + tx, s * p[0] + c * p[1] + ty)


def rigid_fit(a1: Point, a2: Point, b1: Point, b2: Point) -> Rigid:
    """Rotation + translation taking b1 -> a1 and b2 -> a2 (never a reflection)."""
    ang = heading(a1, a2) - heading(b1, b2)
    c, s = math.cos(ang), math.sin(ang)
    return (c, s, a1[0] - (c * b1[0] - s * b1[1]), a1[1] - (s * b1[0] + c * b1[1]))


# ------------------------------------------------------------------------------ slits

SLIT_WIDTH = 1e-3  # how wide the 2D checks take a cut without width (FlatPattern.slits) to be


def cut_slits(part, slits: Iterable[Loop2D], tol: float):
    """The shapely ``part`` with the slits cut out, so they separate what they cut."""
    for sl in slits:
        part = part.difference(LineString(sl.polyline(tol)).buffer(SLIT_WIDTH))
    return part
