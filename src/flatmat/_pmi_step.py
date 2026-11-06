"""The plus-minus tolerances of a STEP file as the file states them, for :mod:`flatmat.pmi`.

OCCT's XCAF dimension objects lose two things that a tolerance check needs:

- the sign of the lower bound (the reader stores ``Abs(lower_bound)``), and
- the unit of angular tolerances (the reader converts an angle's nominal to degrees but
  leaves its tolerance in the file's angle unit).

So :func:`step_dimensions` reads every ``dimensional_location`` / ``dimensional_size``
of the transferred StepData model, in file order, with its ``plus_minus_tolerance``
bounds as written and the size of each bound's own unit (mm or degrees per unit).
:class:`Matcher` then pairs the document's semantic dimensions with these, in order:
OCCT creates one dimension label per such entity, in file order, unless it drops an
entity it cannot read (one whose shapes it cannot find, for one). A pair is used only
when what XCAF kept agrees with the file (the nominal, the upper bound and the lower
bound's magnitude). When OCCT dropped an entity, or a label disagrees with the entity
in its place, the order no longer holds, and a pair is used only if every remaining
entity that agrees states the same signed bounds; otherwise the dimension is unmatched
and its tolerance unreadable. Entities are related by their entity number in the file.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from OCP.StepBasic import (
    StepBasic_ConversionBasedUnit,
    StepBasic_SiPrefix,
    StepBasic_SiUnit,
    StepBasic_SiUnitName,
)
from OCP.StepShape import (
    StepShape_DimensionalCharacteristicRepresentation,
    StepShape_DimensionalLocation,
    StepShape_DimensionalSize,
    StepShape_PlusMinusTolerance,
)

_PREFIXES = {
    "Exa": 1e18,
    "Peta": 1e15,
    "Tera": 1e12,
    "Giga": 1e9,
    "Mega": 1e6,
    "Kilo": 1e3,
    "Hecto": 1e2,
    "Deca": 1e1,
    "Deci": 1e-1,
    "Centi": 1e-2,
    "Milli": 1e-3,
    "Micro": 1e-6,
    "Nano": 1e-9,
    "Pico": 1e-12,
    "Femto": 1e-15,
    "Atto": 1e-18,
}
_METRE = StepBasic_SiUnitName.StepBasic_sunMetre
_RADIAN = StepBasic_SiUnitName.StepBasic_sunRadian


@dataclass(frozen=True)
class Bound:
    """A value as the file writes it, and the size of its unit (mm or degrees), None when unknown."""

    raw: float
    size: float | None

    @property
    def value(self) -> float | None:
        return None if self.size is None else self.raw * self.size


@dataclass(frozen=True)
class StepDimension:
    nominal: Bound | None  # None when the file gives none flatmat can find
    upper: Bound | None  # the plus_minus_tolerance bounds, signed as written; None without one
    lower: Bound | None


def _unit_size(unit, depth: int = 0) -> float | None:
    """mm per unit for a length unit, degrees per unit for a plane angle unit, else None."""
    named = unit.NamedUnit() if hasattr(unit, "NamedUnit") else unit
    if named is None or depth > 4:
        return None
    if isinstance(named, StepBasic_ConversionBasedUnit):
        factor = _measure(named.ConversionFactor())
        if factor is None:
            return None
        inner = _unit_size(factor[1], depth + 1)
        return None if inner is None else factor[0] * inner
    if isinstance(named, StepBasic_SiUnit):
        scale: float | None = 1.0
        if named.HasPrefix():
            scale = _PREFIXES.get(_enum_tail(named.Prefix(), "StepBasic_sp"))
        if scale is None:
            return None
        if named.Name() == _METRE:
            return 1000.0 * scale
        if named.Name() == _RADIAN:
            return math.degrees(1.0) * scale
    return None


def _enum_tail(value: StepBasic_SiPrefix, prefix: str) -> str:
    return getattr(value, "name", str(value)).split(".")[-1].removeprefix(prefix)


def _measure(t) -> tuple[float, Any] | None:
    """``(value, unit)`` of a measure_with_unit, plain or inside a representation item."""
    if t is None:
        return None
    if hasattr(t, "GetMeasureWithUnit"):
        t = t.GetMeasureWithUnit()
    elif hasattr(t, "Measure") and not hasattr(t, "ValueComponent"):
        t = t.Measure()
    if t is None or not hasattr(t, "ValueComponent"):
        return None
    return t.ValueComponent(), t.UnitComponent()


def _bound(t) -> Bound | None:
    m = _measure(t)
    return None if m is None else Bound(m[0], _unit_size(m[1]))


_HARMLESS = (1, 4)  # a value_qualifier's precision_qualifier and value_format_type_qualifier: display only


def _plain(item) -> bool:
    """A measure item without qualifiers, or with ones that only say how it is displayed.

    NX, Creo and SolidWorks qualify a nominal with its decimal places
    (``VALUE_FORMAT_TYPE_QUALIFIER('NR2 2.2')``); that does not change the value."""
    qri = item.GetQualifiedRepresentationItem() if hasattr(item, "GetQualifiedRepresentationItem") else None
    if qri is None:
        return not hasattr(item, "Qualifiers")
    for k in range(1, qri.NbQualifiers() + 1):
        if qri.QualifiersValue(k).CaseNumber() not in _HARMLESS:
            return False
    return True


def _nominal(rep) -> Bound | None:
    """The first plain measure item of a shape_dimension_representation (what OCCT takes as the nominal)."""
    for i in range(1, rep.NbItems() + 1):
        item = rep.ItemsValue(i)
        name = item.Name().ToCString().lower() if item.Name() is not None else ""
        if any(w in name for w in ("upper", "lower", "minimum", "maximum")):
            continue
        if not _plain(item):
            continue
        b = _bound(item)
        if b is not None:
            return b
    return None


def _numbering(model, entities: Sequence[Any]):
    """The entity number of an entity the model holds, 0 when not found.

    ``model.Number`` when the binding finds the entity. OCP 8.0.1.1 finds none (``Contains``
    is False for every entity of a model read from a file), so the number comes from the
    enumeration: ``entities`` holds every wrapper of the model, and pybind hands out the one
    wrapper it holds for an entity."""
    by_wrapper = {id(e): i for i, e in enumerate(entities, 1)}

    def number(e) -> int:
        n = model.Number(e)
        return n if n else by_wrapper.get(id(e), 0)

    return number


def step_dimensions(model) -> tuple[StepDimension, ...]:
    """Every dimensional_location / dimensional_size of ``model`` (a StepData_StepModel), in file order."""
    entities = [model.Value(i) for i in range(1, model.NbEntities() + 1)]
    number = _numbering(model, entities)
    dims: list[int] = []  # entity numbers
    tolerances: dict[int, Any] = {}  # dimension entity number -> its tolerance_value
    nominals: dict[int, Bound | None] = {}
    for i, e in enumerate(entities, 1):
        if isinstance(e, (StepShape_DimensionalLocation, StepShape_DimensionalSize)):
            dims.append(i)
        elif isinstance(e, StepShape_PlusMinusTolerance):
            rng = e.Range()
            if rng.CaseNumber() == 1:  # a tolerance_value, not limits and fits
                tolerances[number(e.TolerancedDimension().Value())] = rng.ToleranceValue()
        elif isinstance(e, StepShape_DimensionalCharacteristicRepresentation):
            nominals[number(e.Dimension().Value())] = _nominal(e.Representation())
    out = []
    for n in dims:
        tv = tolerances.get(n)
        upper = _bound(tv.UpperBound()) if tv is not None else None
        lower = _bound(tv.LowerBound()) if tv is not None else None
        out.append(StepDimension(nominals.get(n), upper, lower))
    return tuple(out)


def _close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-9)


def _agrees(kept: float, b: Bound) -> bool:
    """What XCAF kept equals the file's value converted (lengths), or as written (angular tolerances)."""
    v = b.value
    return _close(kept, abs(b.raw)) or (v is not None and _close(kept, abs(v)))


@dataclass(frozen=True)
class Kept:
    """What OCCT's dimension object kept: the nominal, and the tolerance as XCAF stores it."""

    nominal: float
    upper: float | None  # None: no plus-minus tolerance; for a negative upper bound, None too
    lower_magnitude: float | None
    negative_upper: bool = False  # the reader's 2-value "range": nominal, |lower|


def _matches(k: Kept, s: StepDimension) -> bool:
    if s.nominal is not None and s.nominal.value is not None and not _close(k.nominal, s.nominal.value):
        return False
    if k.lower_magnitude is None:
        return s.upper is None
    if s.upper is None or s.lower is None or not _agrees(k.lower_magnitude, s.lower):
        return False
    if k.negative_upper:
        return s.upper.raw < 0
    return k.upper is not None and s.upper.raw >= 0 and _agrees(k.upper, s.upper)


class Matcher:
    """Pairs the document's semantic dimensions, in label order, with the file's (module docstring).

    ``labels`` is the number of semantic dimension labels in the document. When it equals
    the number of the file's dimension entities, OCCT made one label per entity and the
    k-th label is the k-th entity: a dimension takes its own entity when they agree. At
    the first label that disagrees with its entity the order no longer holds, and from
    there on (that label included) the search rule applies. When the counts differ, OCCT
    dropped some entity, and the search rule applies from the start: a dimension takes the
    first remaining entity that agrees, only if every remaining entity that agrees states
    the same signed bounds.
    """

    def __init__(self, dims: Sequence[StepDimension], labels: int) -> None:
        self.dims = list(dims)
        self.in_order = labels == len(self.dims)
        self.next = 0
        self.turns = 0  # match calls so far

    def match(self, k: Kept) -> StepDimension | None:
        self.turns += 1
        if self.in_order:
            j = self.next
            if j < len(self.dims) and _matches(k, self.dims[j]):
                self.next = j + 1
                return self.dims[j]
            self.in_order = False  # the order no longer holds: search, from this label on
        agree = [j for j in range(self.next, len(self.dims)) if _matches(k, self.dims[j])]
        if not agree or len({_said(self.dims[j]) for j in agree}) > 1:
            return None  # nothing agrees, or entities that agree say different things: no way to tell which
        self.next = agree[0] + 1
        return self.dims[agree[0]]

    def skip(self) -> None:
        """A label that could not be read: in file order it still uses its entity."""
        self.turns += 1
        if self.in_order:
            self.next += 1


def _said(s: StepDimension) -> tuple:
    """What an entity says beyond what XCAF keeps: two candidates that differ here are ambiguous."""
    return (s.upper, s.lower)
