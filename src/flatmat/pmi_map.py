"""PMI on the recognised sheet: the path each tolerance takes, threads, the general tolerance.

:func:`map_pmi` turns the :class:`~flatmat.pmi.PmiItem` s of one solid into a
:class:`SheetPmi`. For each item it gives the :class:`PmiPath` through the part: which
flanges its sides sit on and which bends lie between them. flatmat_web uses the path
to decide which process holds the tolerance (laser, laser plus bends, a bend angle);
the engine does not judge tolerances itself.

A face's flange:

- a face of the reference skin (``skin_a``) is in :attr:`Flange.face_ids`;
- a plane of the other skin (``skin_b``) takes the flange of the reference-skin plane
  parallel to it at the sheet thickness (the partner recognition found);
- a thickness face (a cut edge or a hole wall) takes the flanges of the skin faces it
  shares an edge with, looking through further thickness faces when it touches none.
  A side of an item sits on the flanges all its faces share. A cut side running along
  a bent part touches several flanges, and a side on more than one is unmapped (which
  flange the dimension starts from cannot be told). A diameter's side referenced
  through the hole's rim edge (the hole wall and a skin plane) is the hole wall.
- a diameter or radius whose text names a thread (``M6x1.0 - 6H THRU``) is unmapped:
  the tap holds the tapped size, not the laser; the thread is in :attr:`SheetPmi.threads`.
  Opposite skins of one flange, as two sides (faces or edges on them) or as one side
  holding both (a size dimension), are the sheet thickness when the nominal is the
  sheet's thickness or there is none: unmapped too. Flatness of a cut face is unmapped.
- a bend face, either skin's, and a hole wall that touches one are "on a bend"; a
  face of a formed feature (emboss, extruded hole, louvre, bridge, lance), of an
  engraving, or a countersink's cone and a counterbore's bore and step, is on that
  feature. Both are unmapped. The through-bore of a countersink or counterbore is
  laser-cut: a hole of its flange, with :attr:`PmiPath.feature_id` set.

Between two flanges the path is the shortest chain of bends in the flange/bend graph
(one edge per bend between its flanges), ties broken by the lower bend ids. GD&T
follows the v1 list: position (the path from the feature to its farthest datum),
perpendicularity, parallelism and angularity between adjacent flanges (a bend angle),
and flatness of one flange. Every other type is unmapped, with the reason in plain
English. An item whose tolerance could not be read keeps its path: whether it is
checked is decided downstream.

Threads (``M6``, ``M8x1``, ``M6 THRU``, ``4x M5``: ISO 261 sizes and pitches, see
:func:`parse_thread`) and the general tolerance (``ISO 2768-mK``, ``DIN 6930-m``) are
parsed from every text the reader yields: dimension callouts, presentation names and
notes. :attr:`Thread.pmi_ids` names the item a thread came from (none for a note), so a
caller can tell a thread's own dimension from the thread. A thread's count is the larger of the callout's
``4x`` and the number of holes it points at (coaxial walls of one flange are one hole).
"""

from __future__ import annotations

import heapq
import re
from collections import defaultdict
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Literal

from .diagnostics import Code, warn
from .pmi import PmiData, PmiItem, PmiSide

if TYPE_CHECKING:
    from .model import Feature, SheetPart

PathKind = Literal["laser", "laser_hole", "bend_chain", "bend_angle", "flatness", "unmapped"]


@dataclass(frozen=True)
class PmiPath:
    kind: PathKind
    flanges: tuple[int, ...] = ()  # flange ids the sides sit on
    bends: tuple[int, ...] = ()  # bend ids between them (bend_chain) or the one bend (bend_angle)
    feature_id: int | None = None  # the hole feature, when a side is one
    reason: str = ""  # unmapped only, plain English


@dataclass(frozen=True)
class Thread:
    designation: str  # "M6", "M8x1"
    nominal_mm: float
    pitch_mm: float | None  # None = coarse
    through: bool
    face_ids: tuple[int, ...]  # the hole walls it points at
    count: int = 1  # how many threads: the callout's "4x", or the holes it points at if more
    pmi_ids: tuple[int, ...] = ()  # the PmiItem ids it came from; () for a note


@dataclass(frozen=True)
class GeneralTolerance:
    standard: Literal["ISO 2768", "DIN 6930"]
    linear_class: str  # "f" | "m" | "c" | "v"
    geometric_class: str | None = None  # "H" | "K" | "L"
    text: str = ""


@dataclass(frozen=True)
class SheetPmi:
    status: Literal["none", "semantic", "graphical_only", "failed"]
    items: tuple[tuple[PmiItem, PmiPath], ...] = ()
    threads: tuple[Thread, ...] = ()
    general: GeneralTolerance | None = None
    notes: tuple[str, ...] = ()
    dropped: int = 0
    graphical: int = 0  # graphical-only annotations nobody checks (PmiData.graphical)


# ------------------------------------------------------------------------------ texts

_NOT_A_LETTER = r"(?![^\W\d_])"  # no letter follows, umlauts included: the "f" of "für" is no class
_GENERAL = re.compile(
    r"(?<![^\W_])(?:(?:DIN\s*)?(?:EN\s*)?ISO|DIN)\s*(2768|6930)"  # the standard, "DIN ISO 2768" too
    r"(?:\s*-\s*[12](?!\d))?"  # a part number: ISO 2768-1
    rf"\s*[-–:]?\s*([fmcv])(?:[\s-]?([hkl]){_NOT_A_LETTER})?{_NOT_A_LETTER}",  # the classes: m, mK, m-K
    re.IGNORECASE,
)
_STANDARDS: dict[str, Literal["ISO 2768", "DIN 6930"]] = {"2768": "ISO 2768", "6930": "DIN 6930"}
_THREAD = re.compile(  # M6, M8x1, 4x M5, 4xM5: an M not inside a word, maybe after a count "4x"
    r"(?<![^\W_])(?<![.,])(?:(\d+)\s*[xX×]\s*)?M(\d+(?:[.,]\d+)?)(?:\s*[xX×]\s*(\d+(?:[.,]\d+)?))?(?![\d.,]*\d)"
)
_THROUGH = re.compile(r"\b(?:thru|through|durchgehend|durchg\.?)(?![A-Za-z])", re.IGNORECASE)
# a fastener named before the size: a press-in nut's or a screw's hole is no tapped hole
_FASTENER = ("mutter", "schraube", "nut", "bolt", "screw", "insert", "pem")
_WORD_BEFORE = re.compile(r"([^\W\d_]+)[\s\-–/:]*$")
# ISO 261 metric sizes to M64: (coarse pitch or None for a fine-only size, fine pitches)
_ISO_261: dict[float, tuple[float | None, tuple[float, ...]]] = {
    1.0: (0.25, (0.2,)),
    1.1: (0.25, (0.2,)),
    1.2: (0.25, (0.2,)),
    1.4: (0.3, (0.2,)),
    1.6: (0.35, (0.2,)),
    1.8: (0.35, (0.2,)),
    2.0: (0.4, (0.25,)),
    2.2: (0.45, (0.25,)),
    2.5: (0.45, (0.35,)),
    3.0: (0.5, (0.35,)),
    3.5: (0.6, (0.35,)),
    4.0: (0.7, (0.5,)),
    4.5: (0.75, (0.5,)),
    5.0: (0.8, (0.5,)),
    5.5: (None, (0.5,)),
    6.0: (1.0, (0.75,)),
    7.0: (1.0, (0.75,)),
    8.0: (1.25, (1.0, 0.75)),
    9.0: (1.25, (1.0, 0.75)),
    10.0: (1.5, (1.25, 1.0, 0.75)),
    11.0: (1.5, (1.0, 0.75)),
    12.0: (1.75, (1.5, 1.25, 1.0)),
    14.0: (2.0, (1.5, 1.25, 1.0)),
    15.0: (None, (1.5, 1.0)),
    16.0: (2.0, (1.5, 1.0)),
    17.0: (None, (1.5, 1.0)),
    18.0: (2.5, (2.0, 1.5, 1.0)),
    20.0: (2.5, (2.0, 1.5, 1.0)),
    22.0: (2.5, (2.0, 1.5, 1.0)),
    24.0: (3.0, (2.0, 1.5, 1.0)),
    25.0: (None, (2.0, 1.5, 1.0)),
    26.0: (None, (1.5,)),
    27.0: (3.0, (2.0, 1.5, 1.0)),
    28.0: (None, (2.0, 1.5, 1.0)),
    30.0: (3.5, (3.0, 2.0, 1.5, 1.0)),
    32.0: (None, (2.0, 1.5)),
    33.0: (3.5, (3.0, 2.0, 1.5)),
    35.0: (None, (1.5,)),
    36.0: (4.0, (3.0, 2.0, 1.5)),
    38.0: (None, (1.5,)),
    39.0: (4.0, (3.0, 2.0, 1.5)),
    40.0: (None, (3.0, 2.0, 1.5)),
    42.0: (4.5, (4.0, 3.0, 2.0, 1.5)),
    45.0: (4.5, (4.0, 3.0, 2.0, 1.5)),
    48.0: (5.0, (4.0, 3.0, 2.0, 1.5)),
    50.0: (None, (3.0, 2.0, 1.5)),
    52.0: (5.0, (4.0, 3.0, 2.0, 1.5)),
    55.0: (None, (4.0, 3.0, 2.0, 1.5)),
    56.0: (5.5, (4.0, 3.0, 2.0, 1.5)),
    58.0: (None, (4.0, 3.0, 2.0, 1.5)),
    60.0: (5.5, (4.0, 3.0, 2.0, 1.5)),
    62.0: (None, (4.0, 3.0, 2.0, 1.5)),
    64.0: (6.0, (4.0, 3.0, 2.0, 1.5)),
}


def parse_general_tolerance(text: str) -> GeneralTolerance | None:
    """The general tolerance a text states (``ISO 2768-mK``, ``DIN 6930-m``), or None."""
    m = _GENERAL.search(text)
    if m is None:
        return None
    geometric = m.group(3).upper() if m.group(3) else None
    return GeneralTolerance(_STANDARDS[m.group(1)], m.group(2).lower(), geometric, text.strip())


def parse_thread(text: str) -> tuple[str, float, float | None, bool, int] | None:
    """``(designation, nominal mm, pitch mm or None for coarse, through, count)`` of a metric
    thread in a text; the count is the "4x" before it, else 1.

    Only an ISO 261 size is a thread, and not one named after a fastener (``Einpressmutter
    M6``, ``Schraube M6x20``, ``PEM M4``), after a number (``25 M2``: square metres) or
    before a colon (``M1:2``: a scale). A pitch is one only if ISO 261 lists it for the
    size, coarse or fine; any other number after the size (``M6 x 12 tief``, ``M6x20``) is
    a depth or a length, and the thread is coarse. A fine-only size needs its pitch."""
    for m in _THREAD.finditer(text):
        got = _thread(text, m)
        if got is not None:
            return got
    return None


def _thread(text: str, m: re.Match[str]) -> tuple[str, float, float | None, bool, int] | None:
    nominal_text = m.group(2).replace(",", ".")
    sizes = _ISO_261.get(float(nominal_text))
    if sizes is None:
        return None
    before, after = text[: m.start()], text[m.end() :]
    word = _WORD_BEFORE.search(before)
    if word is not None and any(f in word.group(1).lower() for f in _FASTENER):
        return None
    if m.group(1) is None and re.search(r"\d\s*$", before):
        return None  # "25 M2": a unit after a number
    if re.match(r"\s*:", after):
        return None  # "M1:2": a scale
    coarse, fine = sizes
    pitch_text = m.group(3).replace(",", ".") if m.group(3) else None
    pitch = float(pitch_text) if pitch_text else None
    if pitch is not None and pitch not in (coarse, *fine):
        pitch_text = pitch = None  # a depth or a length, not a pitch
    if pitch is None and coarse is None:
        return None
    designation = f"M{nominal_text}" + (f"x{pitch_text}" if pitch_text else "")
    count = max(1, int(m.group(1))) if m.group(1) else 1
    return designation, float(nominal_text), pitch, _THROUGH.search(text) is not None, count


def _gdt_name(text: str) -> str:
    """``profile of surface 0.1 | A`` -> ``profile of surface``."""
    head = text.split(" | ")[0].split(";")[0].strip()
    m = re.match(r"(.*?)\s+Ø?[-+]?\d", head)
    return (m.group(1) if m else head).strip() or "this GD&T type"


# ------------------------------------------------------------------------------ faces


class _Unmapped(Exception):
    """An item (or one of its sides) has no path; the message is the reason."""


@dataclass(frozen=True)
class _On:
    """Where a face or a side sits: the flanges it may belong to."""

    flanges: frozenset[int]
    feature: int | None = None  # the hole feature (countersink, counterbore)
    roles: frozenset[str] = frozenset()


def _a(name: str) -> str:
    return f"{'an' if name[:1] in 'aeiou' else 'a'} {name}"


def _listed(ids) -> str:
    ids = [str(i) for i in sorted(ids)]
    return ", ".join(ids[:-1]) + f" and {ids[-1]}" if len(ids) > 1 else "".join(ids)


class _Sheet:
    """What recognize() found, indexed for mapping faces: nothing is recomputed."""

    def __init__(self, sheet: SheetPart) -> None:
        ctx = sheet._ctx
        if ctx is None:
            raise ValueError("the sheet has no recognition context")
        self.thickness = sheet.thickness
        self.role = {f.id: f.role for f in sheet.faces}
        self.kind = {f.id: f.kind for f in sheet.faces}
        self.flange_of = {i: fl.id for fl in sheet.flanges for i in fl.face_ids}
        self.bend_faces = {s.face_id for b in sheet.bends for s in b.segments if s.face_id}
        self.sharp = [frozenset(se.faces) for se in ctx.sharp_edges.values()]
        self.feature: dict[int, Feature] = {i: feat for feat in sheet.features for i in feat.face_ids}
        # the laser-cut through-bore of a countersink or counterbore: its cylinder of the through diameter
        self.through_bore = {
            i
            for i, feat in self.feature.items()
            if feat.kind_info.flat == "hole"
            and ctx.faces[i - 1].kind == "cylinder"
            and abs(2 * ctx.faces[i - 1].radius - feat.through_diameter) <= 1e-3
        }
        self.partners: dict[int, int] = ctx.partners
        # bends made in the punching tool (a lance's tab), not on the press brake
        self.tool_bent = {
            feat.bend_id: feat.kind_info.name for feat in sheet.features if feat.bend_id and feat.kind_info.tool_bent
        }
        self.nb: dict[int, set[int]] = defaultdict(set)
        for _, fs in ctx.edges:
            for a in fs:
                self.nb[a].update(b for b in fs if b != a)
        # the flange/bend graph: a bend joins its flanges, and two bends whose faces touch
        # (a joggle's zero-length flat) join each other
        self.links: dict[tuple[str, int], set[tuple[str, int]]] = defaultdict(set)
        for b in sheet.bends:
            for f in b.flanges:
                self._link(("flange", f), ("bend", b.id))
        bend_of = {s.face_id: b.id for b in sheet.bends for s in b.segments if s.face_id}
        for face, bend in bend_of.items():
            for n in self.nb[face]:
                if bend_of.get(n, bend) != bend:
                    self._link(("bend", bend), ("bend", bend_of[n]))

    def _link(self, a: tuple[str, int], b: tuple[str, int]) -> None:
        self.links[a].add(b)
        self.links[b].add(a)

    def thick(self, nominal: float | None) -> bool:
        """A nominal that may be the sheet thickness (None: no nominal to tell)."""
        return nominal is None or abs(nominal - self.thickness) <= max(0.01, 0.01 * self.thickness)

    def _on_bend(self, i: int) -> bool:
        return i in self.bend_faces or (self.role.get(i) == "skin_b" and self.kind.get(i) == "cylinder")

    def _behind(self, i: int) -> int | None:
        """The flange of a skin_b plane: that of its partner on the reference skin."""
        p = self.partners.get(i)
        return self.flange_of.get(p) if p is not None else None

    def _edge_flanges(self, i: int) -> frozenset[int]:
        """The flanges of the skin faces a thickness face shares an edge with, looking
        through further thickness faces when it touches none."""
        seen, frontier, first = {i}, [i], True
        formed = ""  # a formed feature it touches (an extruded hole's rim touches nothing else)
        while frontier:
            found: set[int] = set()
            bend = False
            nxt = []
            for f in frontier:
                for n in sorted(self.nb[f] - seen):
                    seen.add(n)
                    role = self.role.get(n)
                    feat = self.feature.get(n)
                    if self._on_bend(n):
                        bend = True
                    elif feat is not None and feat.kind_info.flat != "hole":
                        formed = formed or f"on {_a(feat.kind_info.name)}"
                    elif role == "formed":
                        formed = formed or "on a formed feature"
                    elif role == "skin_a" and n in self.flange_of:
                        found.add(self.flange_of[n])
                    elif role == "skin_b" and (fl := self._behind(n)) is not None:
                        found.add(fl)
                    elif role == "thickness":
                        nxt.append(n)
            if bend and (not found or (first and self.kind.get(i) == "cylinder")):
                raise _Unmapped("on a bend")  # a hole wall touching a bend is in the bend zone
            if found:
                return frozenset(found)
            frontier, first = nxt, False
        raise _Unmapped(formed or f"face {i} is on no flange")

    def face(self, i: int) -> _On:
        if self._on_bend(i):
            raise _Unmapped("on a bend")
        feat = self.feature.get(i)
        if feat is not None and i not in self.through_bore:  # a cone, a counterbore's bore or step, a form
            raise _Unmapped(f"on {_a(feat.kind_info.name)}")
        role = self.role.get(i)
        if role is None:
            raise _Unmapped(f"face {i} is not on the part")
        if role == "formed":
            raise _Unmapped("on a formed feature")
        if role == "skin_a":
            fl = self.flange_of.get(i)
        elif role == "skin_b":
            fl = self._behind(i)
        else:
            return _On(self._edge_flanges(i), feat.id if feat is not None else None, frozenset({role}))
        if fl is None:
            raise _Unmapped(f"face {i} is on no flange")
        return _On(frozenset({fl}), None, frozenset({role}))

    def side(self, side: PmiSide) -> _On:
        if not side.face_ids:
            raise _Unmapped("refers to the whole part, not to faces")
        ons = [self.face(i) for i in side.face_ids]
        common = frozenset.intersection(*(o.flanges for o in ons))
        if not common:
            ids = set(side.face_ids)
            if any(pair <= ids for pair in self.sharp):
                raise _Unmapped("on a bend")  # the two faces of a sharp bend edge
            raise _Unmapped(f"a side spans flanges {_listed(frozenset().union(*(o.flanges for o in ons)))}")
        feature = next((o.feature for o in ons if o.feature is not None), None)
        return _On(common, feature, frozenset().union(*(o.roles for o in ons)))

    def route(self, a_set: frozenset[int], b_set: frozenset[int]) -> tuple[int, int, tuple[int, ...]] | None:
        """``(flange a, flange b, bends from a to b)``: the fewest bends, then the lower bend ids.

        A side on several flanges (a cut side along a bent part) is unmapped: which of its
        flanges the dimension starts from cannot be told, and picking one could count too
        few bends. A path across a bend made in the punching tool (a lance) is unmapped."""
        for ids in (a_set, b_set):
            _one(ids)
        best: tuple[int, tuple[int, ...], int, int] | None = None
        for a in sorted(a_set):
            reached: dict[tuple[str, int], tuple[int, ...]] = {}  # node -> the bends on the way to it
            heap: list[tuple[int, tuple[int, ...], tuple[str, int]]] = [(0, (), ("flange", a))]
            while heap:
                _, bends, node = heapq.heappop(heap)
                if node in reached:
                    continue
                reached[node] = bends
                for nxt in self.links[node]:
                    if nxt not in reached:
                        way = (*bends, nxt[1]) if nxt[0] == "bend" else bends
                        heapq.heappush(heap, (len(way), way, nxt))
            for b in sorted(b_set):
                if ("flange", b) in reached:
                    bends = reached[("flange", b)]
                    cand = (len(bends), bends, a, b)
                    if best is None or cand < best:
                        best = cand
        if best is None:
            return None
        for bend in best[1]:
            if bend in self.tool_bent:
                raise _Unmapped(f"across {_a(self.tool_bent[bend])}")
        return best[2], best[3], best[1]


def _one(flanges: frozenset[int]) -> int:
    """The one flange a side sits on; a side on several is unmapped."""
    if len(flanges) != 1:
        raise _Unmapped(f"a side spans flanges {_listed(flanges)}")
    return next(iter(flanges))


# ------------------------------------------------------------------------------ paths

_ORIENTATION = ("perpendicularity", "parallelism", "angularity")
_SKIN_A, _SKIN_B = frozenset({"skin_a"}), frozenset({"skin_b"})
_SKINS = _SKIN_A | _SKIN_B


def _chain(on: _Sheet, a: _On, b: _On, feature: int | None, what: str = "") -> PmiPath:
    r = on.route(a.flanges, b.flanges)
    if r is None:
        raise _Unmapped(f"no bend path between the flanges{what}")
    fa, fb, bends = r
    if not bends:
        return PmiPath("laser", (fa,), feature_id=feature)
    return PmiPath("bend_chain", (fa, fb), bends, feature)


def _bend_angle(on: _Sheet, a: _On, b: _On, what: str) -> PmiPath:
    r = on.route(a.flanges, b.flanges)
    if r is None:
        raise _Unmapped(f"no bend path between the flanges of the {what}")
    fa, fb, bends = r
    if not bends:
        raise _Unmapped(f"{what} within one flange")
    if len(bends) > 1:
        raise _Unmapped(f"{what} across {len(bends)} bends")
    return PmiPath("bend_angle", (fa, fb), bends)


def _datum(on: _Sheet, name: str, datums: dict[str, PmiItem]) -> _On:
    item = datums.get(name)
    if item is None:
        raise _Unmapped(f"datum {name} not found")
    if not item.sides:
        raise _Unmapped(f"datum {name} refers to no faces")
    try:
        return on.side(item.sides[0])
    except _Unmapped as e:
        raise _Unmapped(f"datum {name} {e}") from None


def _path(on: _Sheet, item: PmiItem, datums: dict[str, PmiItem]) -> PmiPath:
    kind = item.kind
    if kind == "datum":
        raise _Unmapped("datum")  # a reference, not a tolerance: flatmat_web treats it as a marker
    if kind == "other_gdt":
        raise _Unmapped(f"{_gdt_name(item.text)} is not checked")
    if kind in ("diameter", "radius") and parse_thread(item.text) is not None:
        raise _Unmapped("a thread")  # the tap holds the tapped size, not the laser; reported in threads
    pmi_sides = [_hole_side(on, s) for s in item.sides] if kind in ("diameter", "radius") else item.sides
    sides = [on.side(s) for s in pmi_sides]
    if not sides:
        raise _Unmapped("refers to no faces")
    feature = next((s.feature for s in sides if s.feature is not None), None)
    if kind in ("diameter", "radius"):
        if any(s.roles - {"thickness"} for s in sides):
            raise _Unmapped(f"a {kind} on a skin face")
        if any(on.kind.get(i) == "plane" for s in pmi_sides for i in s.face_ids):
            raise _Unmapped(f"a {kind} on a flat face")
        common = frozenset.intersection(*(s.flanges for s in sides))
        if not common:
            raise _Unmapped(f"a {kind} across flanges")
        return PmiPath("laser_hole", (_one(common),), feature_id=feature)
    if kind == "flatness":
        if any(s.roles - _SKINS for s in sides):
            raise _Unmapped("flatness of a cut face")  # a flange's flatness is that of its skins
        common = frozenset.intersection(*(s.flanges for s in sides))
        if len(common) != 1:
            raise _Unmapped("flatness of a face on several flanges")
        return PmiPath("flatness", tuple(common))
    if kind == "angle":
        if len(sides) != 2:
            raise _Unmapped("an angle needs two sides")
        return _bend_angle(on, sides[0], sides[1], "angle")
    if kind in _ORIENTATION:
        if not item.datums:
            raise _Unmapped(f"{kind} without a datum")
        return _bend_angle(on, sides[0], _datum(on, item.datums[0], datums), kind)
    if kind == "position":
        if not item.datums:
            return _chain(on, sides[0], sides[0], feature)
        paths = [_chain(on, sides[0], _datum(on, d, datums), feature, f" and datum {d}") for d in item.datums]
        return max(paths, key=lambda p: len(p.bends))  # the farthest datum; the first of equals
    # linear. Opposite skins of one flange at the sheet thickness (or with no nominal to
    # tell): the sheet thickness. As two sides, by their faces or by edges on them, or as
    # one side holding both (a dimensional_size).
    if len(sides) == 1:
        if _SKINS <= sides[0].roles and on.thick(item.nominal):
            raise _Unmapped("the sheet thickness")
        return PmiPath("laser", (_one(sides[0].flanges),), feature_id=feature)
    a, b = sides[0], sides[1]
    if {a.roles & _SKINS, b.roles & _SKINS} == {_SKIN_A, _SKIN_B} and a.flanges & b.flanges and on.thick(item.nominal):
        raise _Unmapped("the sheet thickness")
    return _chain(on, a, b, feature)


def _hole_side(on: _Sheet, side: PmiSide) -> PmiSide:
    """A diameter's side referenced through an edge: the reader gives the edge's two faces,
    the hole wall and a skin plane. Reduced to the one hole wall."""
    walls = [i for i in side.face_ids if on.role.get(i) == "thickness" and on.kind.get(i) == "cylinder"]
    rest = [i for i in side.face_ids if i not in walls]
    if len(walls) == 1 and rest and all(on.role.get(i) in _SKINS and on.kind.get(i) == "plane" for i in rest):
        return PmiSide((walls[0],), side.points)
    return side


def _hole_walls(item: PmiItem, sheet: SheetPart) -> tuple[int, ...]:
    walls = {f.id for f in sheet.faces if f.role == "thickness" and f.kind == "cylinder"}
    return tuple(sorted({i for s in item.sides for i in s.face_ids if i in walls}))


def _holes(walls: tuple[int, ...], sheet: SheetPart, on: _Sheet | None) -> int:
    """How many holes the walls are: coaxial walls of one flange are one hole (a wall split
    into two faces, a counterbore's bore over its through-bore). Holes on one axis in two
    flanges (both legs of a U) are two."""
    ctx = sheet._ctx
    if ctx is None or on is None:
        return len(walls)
    axes = set()
    for i in walls:
        try:
            flanges: frozenset[int] | None = on._edge_flanges(i)
        except _Unmapped:
            flanges = None  # on a bend or a form: its axis alone
        f = ctx.faces[i - 1]
        d = f.axis_d.Normalized()
        dx, dy, dz = (round(c, 4) + 0.0 for c in (d.X(), d.Y(), d.Z()))
        if next(c for c in (dx, dy, dz) if c != 0.0) < 0:
            dx, dy, dz = -dx, -dy, -dz
        px, py, pz = f.axis_p.X(), f.axis_p.Y(), f.axis_p.Z()
        t = px * dx + py * dy + pz * dz  # the foot of the axis nearest the origin names the line
        axes.add((dx, dy, dz, *(round(c, 3) + 0.0 for c in (px - t * dx, py - t * dy, pz - t * dz)), flanges))
    return len(axes)


def map_pmi(sheet: SheetPart, data: PmiData) -> SheetPmi:
    """The PMI of ``data`` on ``sheet`` (read from the same solid): see the module docstring.

    An item that fails to map (an error in flatmat, not an unmapped item) is unmapped with
    the reason ``"not mapped: an error in flatmat"``; :func:`with_pmi` also warns."""
    return _map(sheet, data)[0]


def _map(sheet: SheetPart, data: PmiData) -> tuple[SheetPmi, list[Exception]]:
    """:func:`map_pmi`, and the errors of the items that failed to map."""
    items = [i for i in data.items if i.kind != "note"]
    pairs: list[tuple[PmiItem, PmiPath]] = []
    failed: list[Exception] = []
    on = _Sheet(sheet) if items else None
    if on is not None:
        datums: dict[str, PmiItem] = {}
        for i in items:
            if i.kind == "datum" and i.datums:
                datums.setdefault(i.datums[0], i)
        for i in items:
            try:
                path = _path(on, i, datums)
            except _Unmapped as e:
                path = PmiPath("unmapped", reason=str(e))
            except Exception as e:  # one item loses its path, not the others theirs
                failed.append(e)
                path = PmiPath("unmapped", reason="not mapped: an error in flatmat")
            pairs.append((i, path))
    threads = []
    general = None
    for i in data.items:
        try:
            got = parse_thread(i.text)
            if got is not None:
                designation, nominal, pitch, through, count = got
                walls = _hole_walls(i, sheet)
                if on is None and sheet._ctx is not None:
                    on = _Sheet(sheet)
                holes = max(count, _holes(walls, sheet, on), 1)
                ids = () if i.kind == "note" else (i.id,)
                threads.append(Thread(designation, nominal, pitch, through, walls, holes, ids))
        except Exception as e:
            failed.append(e)
        if general is None:
            general = parse_general_tolerance(i.text)
    pmi = SheetPmi(
        data.status,
        tuple(pairs),
        tuple(threads),
        general,
        tuple(i.text for i in data.items if i.kind == "note"),
        data.dropped,
        data.graphical,
    )
    return pmi, failed


def with_pmi(sheet: SheetPart, data: PmiData) -> SheetPart:
    """``sheet`` with its PMI mapped and the read's diagnostics added; never raises."""
    if data.status == "none" and not data.diagnostics:
        return sheet  # a file without PMI: nothing changes
    try:
        pmi, failed = _map(sheet, data)
        diags = data.diagnostics
        if failed:
            msg = f"{len(failed)} PMI item(s) not mapped onto the sheet: {failed[0]!r}"
            diags = (*diags, warn(Code.PMI_READ_FAILED, msg, count=len(failed)))
    except Exception as e:  # PMI never refuses a part
        pmi = SheetPmi("failed", dropped=data.dropped, graphical=data.graphical)
        diags = (*data.diagnostics, warn(Code.PMI_READ_FAILED, f"PMI not mapped onto the sheet: {e!r}"))
    return replace(sheet, pmi=pmi, diagnostics=(*sheet.diagnostics, *diags))
