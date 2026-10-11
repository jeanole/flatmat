"""Press-brake tooling: a small parametric punch and die schema (docs/design.md §7.5).

There is no open standard for tool geometry (docs/research.md §10.8), so a tool is a
handful of numbers that every catalogue gives: the punch's height, tip angle and
radius, body width and, for a gooseneck, its throat; the die's V opening, angle,
shoulder radius and block size. :class:`Tooling` loads and saves them as JSON::

    {"punch": {"name": "PS.88", "height": 100, "tip_angle_deg": 88, "tip_radius": 0.8},
     "die":   {"name": "V12",   "v": 12, "angle_deg": 88, "shoulder_radius": 1.5}}

The built-in :data:`PUNCHES` and :data:`DIES` are typical values for European-style
tooling, entered by hand for examples and defaults; they are not any maker's data.

Section coordinates (used by the collision check): ``y`` across the bend line, ``z``
up. The punch's lowest point (the bottom of its tip radius) is at the origin; the die's
top face is at ``z = 0`` with the V centred on ``y = 0``.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import asdict, dataclass, field, fields, replace

from shapely.geometry import Polygon, box

_FAR = 2000.0  # how far the clamp and the lower beam reach (mm)


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


@dataclass(frozen=True)
class Punch:
    """Upper tool. Lengths in mm, angles in degrees.

    A straight punch is a wedge of ``tip_angle_deg`` that widens to ``width`` and runs
    up to ``height``, where the clamp (``clamp_width`` wide) begins. A gooseneck keeps
    that nose up to ``throat_bottom``, then leans back at 45° by ``throat_depth``; the
    recess on its +y side, up to ``throat_top`` (default: the full height), is the
    throat where a returned flange fits. ``length`` None means the punch is made up of segments as long as
    each bend. The clamp runs the machine's full length.
    """

    name: str = "straight"
    height: float = 100.0
    tip_angle_deg: float = 88.0
    tip_radius: float = 0.8
    width: float = 12.0
    throat_depth: float = 0.0
    throat_bottom: float = 0.0
    throat_top: float | None = None
    length: float | None = None
    clamp_width: float = 60.0

    @property
    def is_gooseneck(self) -> bool:
        return self.throat_depth > 0

    def section(self, tip_radius: float | None = None) -> Polygon:
        """The punch and its clamp in section coordinates, the tip's lowest point at the origin."""
        body, clamp = self.sections(tip_radius)
        return body.union(clamp)

    def sections(self, tip_radius: float | None = None) -> tuple[Polygon, Polygon]:
        """(punch, clamp) in section coordinates, the tip's lowest point at the origin."""
        rp = self.tip_radius if tip_radius is None else max(0.0, tip_radius)
        half = math.radians(self.tip_angle_deg) / 2
        w = self.width
        hw = (w / 2) / math.tan(half)  # where the flanks reach the body width
        H = max(self.height, hw + 1.0)
        if not self.is_gooseneck:
            pts = [(0.0, 0.0), (w / 2, hw), (w / 2, H), (-w / 2, H), (-w / 2, hw)]
        else:
            # the neck leans back at 45° from the top of the nose, so a 90° bend's leg
            # (also at 45°) clears its underside; the recess above it is the throat
            d = self.throat_depth
            z1 = max(self.throat_bottom, hw + 1.0)
            zf = z1 + d + w / 2  # where the throat's floor meets its back wall
            z2 = min(self.throat_top if self.throat_top is not None else H, H)
            if z2 <= zf:
                raise ValueError(
                    f"punch {self.name}: throat_top must lie above {zf:g} (throat_bottom + depth + width/2)"
                )
            pts = [(0.0, 0.0), (w / 2, hw), (w / 2, z1), (-d, zf), (-d, z2)]
            pts += [(w / 2, z2), (w / 2, H)] if z2 < H else [(-d, H)]
            pts += [(-d - w, H), (-d - w, zf), (-w / 2, z1), (-w / 2, hw)]
        body = Polygon(pts)
        if rp > 0:
            # round the tip, then put its lowest point back at the origin
            body = body.buffer(-rp, join_style="round").buffer(rp, join_style="round")
            lift = body.bounds[1]
            from shapely.affinity import translate

            body = translate(body, 0.0, -lift)
            H = H - lift
        ys = [p[0] for p in pts]
        mid = (min(ys) + max(ys)) / 2
        clamp = box(mid - self.clamp_width / 2, H, mid + self.clamp_width / 2, H + _FAR)
        return body, clamp


@dataclass(frozen=True)
class Die:
    """Lower tool: a block ``width`` wide and ``height`` high with a V of ``angle_deg``.

    ``v`` is measured between the sharp edges where the V's flanks meet the top face
    (before the ``shoulder_radius`` rounds them). Below the block, a holder and the
    lower beam ``base_width`` wide, the machine's full length. ``length`` None means
    the die is made up of segments as long as each bend.
    """

    name: str = "V"
    v: float = 12.0
    angle_deg: float = 88.0
    shoulder_radius: float = 1.0
    width: float | None = None
    height: float = 60.0
    length: float | None = None
    base_width: float = 60.0

    @property
    def top_width(self) -> float:
        return self.width if self.width is not None else max(self.v + 20.0, 1.6 * self.v)

    @property
    def depth(self) -> float:
        """Depth of the V below the top face."""
        return (self.v / 2) / math.tan(math.radians(self.angle_deg) / 2)

    def shoulder_center(self) -> tuple[float, float]:
        """Centre of the right shoulder's rounding, top face at z = 0."""
        half = math.radians(self.angle_deg) / 2
        e1 = (1.0, 0.0)  # along the top face, outwards
        e2 = (-math.sin(half), -math.cos(half))  # down the flank
        b = (e1[0] + e2[0], e1[1] + e2[1])
        nb = math.hypot(*b)
        gamma = math.acos(e1[0] * e2[0] + e1[1] * e2[1])  # interior angle at the edge
        d = self.shoulder_radius / math.sin(gamma / 2)
        return (self.v / 2 + d * b[0] / nb, d * b[1] / nb)

    def section(self) -> Polygon:
        """The die block and the beam below it, top face at z = 0."""
        block, beam = self.sections()
        return block.union(beam)

    def sections(self) -> tuple[Polygon, Polygon]:
        """(die block, lower beam), top face at z = 0."""
        W = max(self.top_width, self.v + 4 * self.shoulder_radius + 2.0)
        H = max(self.height, self.depth + 1.0)
        over = self.v / 2 + math.tan(math.radians(self.angle_deg) / 2)  # the flanks continued to z = 1
        vee = Polygon([(-over, 1.0), (0.0, -self.depth), (over, 1.0)])
        block = box(-W / 2, -H, W / 2, 0.0).difference(vee)
        if self.shoulder_radius > 0:
            # round the convex edges (the shoulders) only: an opening by r
            block = block.buffer(-self.shoulder_radius, join_style="round").buffer(
                self.shoulder_radius, join_style="round"
            )
        beam = box(-self.base_width / 2, -H - _FAR, self.base_width / 2, -H)
        return block, beam


@dataclass(frozen=True)
class Tooling:
    punch: Punch = field(default_factory=Punch)
    die: Die = field(default_factory=Die)

    @classmethod
    def default(cls, t: float) -> Tooling:
        """A straight 88° punch and an 88° die with V from the Amada rule."""
        v = die_opening(t)
        return cls(PUNCHES["straight_88"], Die(f"V{v:g}", v, 88.0, _shoulder(v)))

    @classmethod
    def bendfm(cls) -> Tooling:
        """BenDFM's fixed tools: a 90° punch 10 mm thick and 300 mm high, a 40 mm die."""
        return cls(PUNCHES["bendfm"], DIES["bendfm"])

    @classmethod
    def named(cls, punch: str, die: str) -> Tooling:
        """From the built-in catalogues; ``die`` may also be ``"V<opening>"`` (88°)."""
        if punch not in PUNCHES:
            raise KeyError(f"unknown punch {punch!r}; built in: {', '.join(PUNCHES)}")
        if die not in DIES and die.upper().startswith("V") and die[1:].replace(".", "", 1).isdigit():
            v = float(die[1:])
            return cls(PUNCHES[punch], Die(f"V{v:g}", v, 88.0, _shoulder(v)))
        if die not in DIES:
            raise KeyError(f"unknown die {die!r}; built in: {', '.join(DIES)} or V<opening>")
        return cls(PUNCHES[punch], DIES[die])

    def to_dict(self) -> dict:
        return {"punch": asdict(self.punch), "die": asdict(self.die)}

    @classmethod
    def from_dict(cls, d: dict) -> Tooling:
        return cls(_build(Punch, d.get("punch", {})), _build(Die, d.get("die", {})))

    @classmethod
    def load(cls, path: str | os.PathLike) -> Tooling:
        with open(os.fspath(path)) as fh:
            return cls.from_dict(json.load(fh))

    def save(self, path: str | os.PathLike) -> None:
        with open(os.fspath(path), "w") as fh:
            json.dump(self.to_dict(), fh, indent=1)

    def with_die(self, **kw) -> Tooling:
        return replace(self, die=replace(self.die, **kw))


def _build(cls, d: dict):
    names = {f.name for f in fields(cls)}
    unknown = set(d) - names
    if unknown:
        raise ValueError(f"unknown {cls.__name__} field(s): {', '.join(sorted(unknown))}")
    return cls(**d)


def _shoulder(v: float) -> float:
    """A typical shoulder radius for a V opening (catalogues: ~0.1 V, at least 0.5 mm)."""
    return round(max(0.5, 0.1 * v), 1)


PUNCHES: dict[str, Punch] = {
    "straight_88": Punch("straight_88", height=100.0, tip_angle_deg=88.0, tip_radius=0.8, width=12.0),
    "gooseneck_88": Punch(
        "gooseneck_88",
        height=110.0,
        tip_angle_deg=88.0,
        tip_radius=0.8,
        width=10.0,
        throat_depth=40.0,
        throat_bottom=20.0,
        throat_top=100.0,
    ),  # fmt: skip
    "acute_30": Punch("acute_30", height=100.0, tip_angle_deg=30.0, tip_radius=0.6, width=10.0),
    # BenDFM tests 2D sections: its "10 mm thick, 300 mm long" punch is a 10 x 300 slab
    "bendfm": Punch("bendfm", height=300.0, tip_angle_deg=90.0, tip_radius=1.0, width=10.0, clamp_width=10.0),
}

DIES: dict[str, Die] = {
    **{f"V{v:g}": Die(f"V{v:g}", float(v), 88.0, _shoulder(v)) for v in (6, 8, 10, 12, 16, 20, 24, 32, 40, 50)},
    "V16_30": Die("V16_30", 16.0, 30.0, 1.6),
    "bendfm": Die("bendfm", 40.0, 90.0, 2.0),
}

__all__ = ["DIES", "PUNCHES", "Die", "Punch", "Tooling", "die_opening"]
