"""Bend allowance, K-factor models and DIN 6935.

Conventions: ``K`` is the "ANSI" K-factor, the neutral fibre's distance from the
inside surface divided by the thickness (0 <= K <= 0.5 in practice). ``angle``
is the bend (deflection) angle; DIN 6935's opening angle is ``180° - angle``.
See docs/research.md §6 for sources.
"""

from __future__ import annotations

import bisect
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

from .diagnostics import Code, Diagnostic, warn


def bend_allowance(angle_rad: float, r_in: float, t: float, k: float) -> float:
    """Length of the neutral fibre through the bend: ``θ·(r + K·t)``."""
    return abs(angle_rad) * (r_in + k * t)


def outside_setback(angle_rad: float, r_in: float, t: float) -> float | None:
    """``tan(θ/2)·(r + t)``; undefined for hems (θ ≥ 179°)."""
    if abs(angle_rad) >= math.radians(179.0):
        return None
    return math.tan(abs(angle_rad) / 2) * (r_in + t)


def bend_deduction(angle_rad: float, r_in: float, t: float, k: float) -> float | None:
    """``2·OSSB − BA`` (legs measured to the outer mould lines)."""
    ossb = outside_setback(angle_rad, r_in, t)
    return None if ossb is None else 2 * ossb - bend_allowance(angle_rad, r_in, t, k)


# ------------------------------------------------------------------------------- DIN 6935


def din6935_k(r_over_s: float) -> float:
    """DIN 6935 correction factor ``k = 0.65 + 0.5·lg(r/s)``, ``k = 1`` for ``r/s > 5``.

    The standard does not cover ``r/s <= 0.65`` (its Table 5 starts "over 0.65"); the
    formula is evaluated there anyway and clamped at 0. Callers that care use
    :func:`din6935_k_checked` to get a warning.
    """
    if r_over_s > 5:
        return 1.0
    if r_over_s <= 0:
        return 0.0
    return max(0.0, 0.65 + 0.5 * math.log10(r_over_s))


def din6935_v(opening_deg: float, r: float, s: float) -> float:
    """DIN 6935 compensation value ``v`` (Ausgleichswert): ``L = a + b + v``.

    Legs ``a``, ``b`` are dimensioned to the outer vertex; ``v`` is normally negative.
    ``opening_deg`` is the opening angle β between the legs (180° = flat).
    """
    beta = opening_deg
    k = din6935_k(r / s)
    arc = math.pi * (180.0 - beta) / 180.0 * (r + s / 2 * k)
    if beta <= 90.0:
        return arc - 2 * (r + s)
    if beta <= 165.0:
        return arc - 2 * (r + s) * math.tan(math.radians(180.0 - beta) / 2)
    return 0.0


# ------------------------------------------------------------------------------- K models


@dataclass(frozen=True)
class KScalar:
    k: float

    def __call__(self, r_in: float, t: float, angle_rad: float) -> tuple[float, Diagnostic | None]:
        return self.k, None

    def describe(self) -> str:
        return f"K={self.k:g}"


@dataclass(frozen=True)
class KTable:
    """K as a function of r/t, piecewise linear between rows, clamped outside them."""

    r_over_t: tuple[float, ...]
    k: tuple[float, ...]
    name: str = "table"

    def __post_init__(self) -> None:
        if len(self.r_over_t) != len(self.k) or len(self.k) < 1:
            raise ValueError("r_over_t and k must have the same, non-zero length")
        if list(self.r_over_t) != sorted(self.r_over_t):
            raise ValueError("r_over_t must be ascending")

    @classmethod
    def from_rows(cls, rows: Sequence[tuple[float, float]], name: str = "table") -> KTable:
        rows = sorted(rows)
        return cls(tuple(r for r, _ in rows), tuple(k for _, k in rows), name)

    @classmethod
    def din6935(cls) -> KDin6935:
        return KDin6935()

    def __call__(self, r_in: float, t: float, angle_rad: float) -> tuple[float, Diagnostic | None]:
        x = r_in / t
        xs, ks = self.r_over_t, self.k
        if x <= xs[0] or x >= xs[-1]:
            k = ks[0] if x <= xs[0] else ks[-1]
            outside = x < xs[0] or x > xs[-1]
            d = warn(Code.K_EXTRAPOLATED, f"r/t={x:.3g} outside the {self.name} table", r_over_t=x) if outside else None
            return k, d
        i = bisect.bisect_right(xs, x)
        x0, x1, k0, k1 = xs[i - 1], xs[i], ks[i - 1], ks[i]
        return k0 + (k1 - k0) * (x - x0) / (x1 - x0), None

    def describe(self) -> str:
        return self.name


@dataclass(frozen=True)
class KDin6935:
    """K from DIN 6935: ``K = k_DIN / 2``."""

    def __call__(self, r_in: float, t: float, angle_rad: float) -> tuple[float, Diagnostic | None]:
        x = r_in / t
        d = None
        if x <= 0.65:
            d = warn(Code.K_OUTSIDE_STANDARD, f"DIN 6935 does not cover r/s={x:.3g} <= 0.65", r_over_t=x)
        return din6935_k(x) / 2, d

    def describe(self) -> str:
        return "DIN 6935"


KCallable = Callable[[float, float, float], float]
KModel = KScalar | KTable | KDin6935


def as_k_model(k: float | KModel | KCallable) -> KModel | _KFunc:
    if isinstance(k, int | float):
        if not 0.0 <= float(k) <= 1.0:
            raise ValueError(f"K-factor {k} outside [0, 1]")
        return KScalar(float(k))
    if isinstance(k, KScalar | KTable | KDin6935):
        return k
    if callable(k):
        return _KFunc(k)
    raise TypeError(f"unsupported K model: {k!r}")


@dataclass(frozen=True)
class _KFunc:
    fn: KCallable

    def __call__(self, r_in: float, t: float, angle_rad: float) -> tuple[float, Diagnostic | None]:
        return float(self.fn(r_in, t, angle_rad)), None

    def describe(self) -> str:
        return getattr(self.fn, "__name__", "callable")
