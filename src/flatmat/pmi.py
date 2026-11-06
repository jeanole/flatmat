"""Semantic PMI (AP242 GD&T) from the XCAF document :func:`flatmat.load_step` keeps.

:func:`read_pmi` turns the document's dimensions, geometric tolerances, datums and
notes into :class:`PmiItem` s that name their faces by flatmat face id (the
``TopExp.MapShapes`` index on the solid, 1-based). Mapping the items onto flanges
and bends happens later (:mod:`flatmat.pmi_map`).

What OCCT's reader gives, and what flatmat makes of it (verified on OCP 8.0.1.1):

- Lengths come in mm and angles in degrees. OCCT scales values, tolerances and points
  by the file's length unit and converts angular nominals to degrees, but keeps only
  the *magnitude* of a lower bound and leaves angular tolerances in the file's angle
  unit. So a document read from STEP takes each tolerance from the file itself
  (:mod:`flatmat._pmi_step`): the signed bounds as written, converted by their own
  unit. A tolerance flatmat cannot find in the file, or whose unit it does not know,
  is unreadable: ``upper`` and ``lower`` are None, and ``text`` keeps what was read.
- ``upper`` and ``lower`` are signed deviations from the nominal: +0.1/+0.05,
  0/-0.2 and -0.05/-0.1 keep their signs. An inverted pair (upper below lower) or a NaN
  is unreadable. A dimension given only as lower and upper limits has no nominal:
  ``nominal``, ``upper`` and ``lower`` are None and ``text`` keeps both limits.
- :attr:`PmiItem.tolerance` says which of these a dimension is: ``"none"`` (no
  tolerance of its own, so the general tolerance applies), ``"plus_minus"`` (``upper``
  and ``lower`` hold it), ``"limits"`` (lower and upper limits, or a max / min
  qualifier: no nominal), ``"fit"`` (an ISO 286 class, ``Ø8 H7`` in ``text``) or
  ``"unreadable"``. OCCT's reader drops a STEP file's qualifiers, and OCP cannot return
  XCAF's class of tolerance, so both come from the file: an untoleranced dimension
  flatmat cannot find in the file is unreadable too, since a qualifier or a fit cannot
  be ruled out. GD&T items are ``"none"`` (the zone is ``upper``) or ``"unreadable"``
  (no zone).
- Graphical-only PMI (a presentation without semantic data) is read as dimensions
  of type ``DimensionPresentation`` or ``CommonLabel``. They never become items: a
  dimension written both ways is one item. A file with no semantic PMI at all has the
  status ``"graphical_only"``. Beside semantic PMI, a graphical annotation whose faces are
  not among one semantic item's faces (or that refers to no faces) is a callout nobody checks:
  :attr:`PmiData.graphical` counts them, with one ``PMI_GRAPHICAL_ONLY`` warning. One on
  another solid is not counted; one that repeats a semantic item's faces is silent.
- Notes do not survive OCCT's STEP round trip (its writer drops them and its reader
  fills no notes table); they are read when a document has them from another source.
  A dimension's ``text`` carries the callout OCCT gives it: the values, its
  descriptions and its presentation's name.

A referenced edge or vertex becomes the faces it bounds. A face is found on the solid
as it is, and else on the solid without its location (PMI defined on a part that the
file places as an assembly instance). An item with any face that belongs to neither is
counted in :attr:`PmiData.dropped` and left out. A reference to the whole part (a
solid or compound sharing faces with the solid) is kept with no faces on that side.
Item ids count the semantic items in document order, so an item keeps its id whichever
solid it is read for.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any, Literal

from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.TCollection import TCollection_AsciiString
from OCP.TDF import TDF_Label, TDF_Tool
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX
from OCP.TopExp import TopExp, TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.XCAFDimTolObjects import (
    XCAFDimTolObjects_DimensionType as DimType,
)
from OCP.XCAFDimTolObjects import (
    XCAFDimTolObjects_GeomToleranceType as TolType,
)
from OCP.XCAFDoc import (
    XCAFDoc_AssemblyItemRef,
    XCAFDoc_Datum,
    XCAFDoc_Dimension,
    XCAFDoc_DimTolTool,
    XCAFDoc_DocumentTool,
    XCAFDoc_GeomTolerance,
    XCAFDoc_NoteBinData,
    XCAFDoc_NoteComment,
    XCAFDoc_ShapeTool,
)

from ._occ import EdgeFaceMap, LabelSequence, ShapeMap, as_face
from ._pmi_step import Kept, Matcher, StepDimension
from .diagnostics import Code, Diagnostic, warn

PmiKind = Literal[
    "linear",
    "diameter",
    "radius",
    "angle",
    "position",
    "perpendicularity",
    "parallelism",
    "angularity",
    "flatness",
    "other_gdt",
    "datum",
    "note",
]
Point = tuple[float, float, float]
Tolerance = Literal["none", "plus_minus", "unreadable", "limits", "fit"]


@dataclass(frozen=True)
class PmiSide:
    face_ids: tuple[int, ...]  # MapShapes index on the solid, 1-based; edges/vertices reduced to faces
    points: tuple[Point, ...] = ()  # the dimension's own points on this side, mm


@dataclass(frozen=True)
class PmiItem:
    id: int  # 1-based, in document label order
    kind: PmiKind
    nominal: float | None  # mm or deg; None for datums, notes, GD&T without one
    upper: float | None  # signed deviation from nominal (dims) or zone width (GD&T, upper only)
    lower: float | None
    datums: tuple[str, ...] = ()  # "A", "B" ...
    sides: tuple[PmiSide, ...] = ()  # 1 or 2
    anchor: Point | None = None
    text: str = ""  # raw callout / note text, never lost
    tolerance: Tolerance = "none"  # what kind of tolerance a dimension has (module docstring)


@dataclass(frozen=True)
class PmiData:
    status: Literal["none", "semantic", "graphical_only", "failed"]
    items: tuple[PmiItem, ...]
    dropped: int = 0  # items whose faces belong to another solid
    diagnostics: tuple[Diagnostic, ...] = ()
    graphical: int = 0  # graphical-only annotations, not checked: no semantic item refers to their faces


_DIM_KINDS: dict[Any, PmiKind] = {
    DimType.XCAFDimTolObjects_DimensionType_Location_LinearDistance: "linear",
    DimType.XCAFDimTolObjects_DimensionType_Location_LinearDistance_FromCenterToOuter: "linear",
    DimType.XCAFDimTolObjects_DimensionType_Location_LinearDistance_FromCenterToInner: "linear",
    DimType.XCAFDimTolObjects_DimensionType_Location_LinearDistance_FromOuterToCenter: "linear",
    DimType.XCAFDimTolObjects_DimensionType_Location_LinearDistance_FromOuterToOuter: "linear",
    DimType.XCAFDimTolObjects_DimensionType_Location_LinearDistance_FromOuterToInner: "linear",
    DimType.XCAFDimTolObjects_DimensionType_Location_LinearDistance_FromInnerToCenter: "linear",
    DimType.XCAFDimTolObjects_DimensionType_Location_LinearDistance_FromInnerToOuter: "linear",
    DimType.XCAFDimTolObjects_DimensionType_Location_LinearDistance_FromInnerToInner: "linear",
    DimType.XCAFDimTolObjects_DimensionType_Size_Thickness: "linear",
    DimType.XCAFDimTolObjects_DimensionType_Size_Diameter: "diameter",
    DimType.XCAFDimTolObjects_DimensionType_Size_Radius: "radius",
    DimType.XCAFDimTolObjects_DimensionType_Location_Angular: "angle",
    DimType.XCAFDimTolObjects_DimensionType_Size_Angular: "angle",
}  # every other semantic dimension type (curved, with path, spherical ...) is "other_gdt"
_GRAPHICAL = (
    DimType.XCAFDimTolObjects_DimensionType_DimensionPresentation,
    DimType.XCAFDimTolObjects_DimensionType_CommonLabel,
)
_TOL_KINDS: dict[Any, PmiKind] = {
    TolType.XCAFDimTolObjects_GeomToleranceType_Position: "position",
    TolType.XCAFDimTolObjects_GeomToleranceType_Perpendicularity: "perpendicularity",
    TolType.XCAFDimTolObjects_GeomToleranceType_Parallelism: "parallelism",
    TolType.XCAFDimTolObjects_GeomToleranceType_Angularity: "angularity",
    TolType.XCAFDimTolObjects_GeomToleranceType_Flatness: "flatness",
}
_PREFIX = {"diameter": "Ø", "radius": "R"}


def _enum_name(value: Any, prefix: str) -> str:
    """``XCAFDimTolObjects_GeomToleranceType_ProfileOfSurface`` -> ``profile of surface``."""
    name = getattr(value, "name", str(value)).split(".")[-1].removeprefix(prefix)
    out = ""
    for ch in name.replace("_", " "):
        out += " " + ch.lower() if ch.isupper() and out and out[-1] != " " else ch.lower()
    return out.strip()


def _num(v: float) -> str:
    return f"{v:.6g}"


def tolerance(nominal: float | None, upper: float | None, lower: float | None) -> tuple[float | None, float | None]:
    """Signed deviations ``(upper, lower)`` as they are kept: an inverted pair or a NaN gives ``(None, None)``."""
    if upper is None or lower is None:
        return None, None
    if any(v is not None and math.isnan(v) for v in (nominal, upper, lower)) or upper < lower:
        return None, None
    return upper, lower


def _seq(labels) -> list[TDF_Label]:
    return [labels.Value(i) for i in range(1, labels.Length() + 1)]


def _text(s) -> str:
    """A TCollection string (or None) as Python text."""
    if s is None:
        return ""
    if hasattr(s, "ToExtString"):
        return s.ToExtString()
    return s.ToCString()


class _Faces:
    """flatmat face ids of the shapes PMI refers to."""

    def __init__(self, solid) -> None:
        self.solid = solid
        self._maps: list[tuple[Any, Any, Any]] = []

    @property
    def maps(self) -> list[tuple[Any, Any, Any]]:
        """(faces, edge -> faces, vertex -> faces) of the solid, then of the solid without its location."""
        if not self._maps:  # built at the first reference: a file without PMI costs nothing
            for s in (self.solid, self.solid.Located(TopLoc_Location())):
                fmap, efm, vfm = ShapeMap(), EdgeFaceMap(), EdgeFaceMap()
                TopExp.MapShapes_s(s, TopAbs_FACE, fmap)
                TopExp.MapShapesAndAncestors_s(s, TopAbs_EDGE, TopAbs_FACE, efm)
                TopExp.MapShapesAndAncestors_s(s, TopAbs_VERTEX, TopAbs_FACE, vfm)
                self._maps.append((fmap, efm, vfm))
        return self._maps

    def _face_id(self, face) -> int:
        for fmap, _, _ in self.maps:
            i = fmap.FindIndex(face)
            if i:
                return i
        return 0

    def of(self, shape) -> tuple[int, ...] | None:
        """The face ids ``shape`` reduces to: () for the whole part, None when it is not on the solid."""
        t = shape.ShapeType()
        if t == TopAbs_FACE:
            i = self._face_id(shape)
            return (i,) if i else None
        if t in (TopAbs_EDGE, TopAbs_VERTEX):
            for fmap, efm, vfm in self.maps:
                m = efm if t == TopAbs_EDGE else vfm
                j = m.FindIndex(shape)
                if j:
                    return tuple(sorted({fmap.FindIndex(f) for f in m.FindFromIndex(j)}))
            return None
        ex = TopExp_Explorer(shape, TopAbs_FACE)  # a solid, shell or compound: the whole part if it shares a face
        while ex.More():
            if self._face_id(ex.Current()):
                return ()
            ex.Next()
        return None

    def side(self, labels: Iterable[TDF_Label]) -> tuple[int, ...] | None:
        ids: set[int] = set()
        for lab in labels:
            got = self.of(XCAFDoc_ShapeTool.GetShape_s(lab))
            if got is None:
                return None
            ids.update(got)
        return tuple(sorted(ids))

    def centroid(self, face_ids: Iterable[int]) -> Point | None:
        g, fmap = GProp_GProps(), self.maps[0][0]
        for i in face_ids:
            one = GProp_GProps()  # SurfaceProperties overwrites: add face by face
            BRepGProp.SurfaceProperties_s(as_face(fmap.FindKey(i)), one)
            g.Add(one)
        if g.Mass() <= 0:
            return None
        return g.CentreOfMass().Coord()


def _anchor(text_point: Point | None, sides: tuple[PmiSide, ...], faces: _Faces) -> Point | None:
    if text_point is not None:
        return text_point
    pts = [p for s in sides for p in s.points]
    if pts:
        return (
            sum(p[0] for p in pts) / len(pts),
            sum(p[1] for p in pts) / len(pts),
            sum(p[2] for p in pts) / len(pts),
        )
    return faces.centroid({i for s in sides for i in s.face_ids})


def _refs(label) -> tuple[list[TDF_Label], list[TDF_Label]]:
    first, second = LabelSequence(), LabelSequence()
    XCAFDoc_DimTolTool.GetRefShapeLabel_s(label, first, second)
    return _seq(first), _seq(second)


def _sides(faces: _Faces, refs: list[list[TDF_Label]], points: list[Point | None]):
    """The sides of an item, or None when one of them is not on the solid."""
    sides = []
    for k, labels in enumerate(refs):
        p = points[k] if k < len(points) else None
        if not labels and p is None:
            continue
        ids = faces.side(labels)
        if ids is None:
            return None
        sides.append(PmiSide(ids, (p,) if p is not None else ()))
    return tuple(sides)


def _signed(v: float) -> str:
    return "0" if v == 0 else f"{v:+.6g}"


_QUALIFIERS = {"maximum": "max", "minimum": "min", "average": "avg"}
_FORM = re.compile(r"[A-Za-z]{1,2}")  # an ISO 286 fundamental deviation: H, h, JS, zc ...
_GRADE = re.compile(r"0|01|[1-9]|1[0-8]")  # an IT grade


def _fit(found: StepDimension | None) -> str:
    """The fit class as the file writes it (``H7``, ``h6``), or why it cannot be given."""
    if found is None or found.fit is None:
        return "fit (class not read)"
    form, grade = found.fit
    if _FORM.fullmatch(form) and _GRADE.fullmatch(grade):
        return f"{form}{grade}"
    return f"fit (class unreadable: {form} {grade})".replace("  ", " ")


def _xcaf_qualifier(obj) -> str:
    if not obj.HasQualifier():
        return ""
    name = _enum_name(obj.GetQualifier(), "XCAFDimTolObjects_DimensionQualifier_")
    return "" if name == "none" else name


def _values(obj, matcher: Matcher | None, unit: str) -> tuple[float | None, float | None, float | None, Tolerance, str]:
    """``(nominal, upper, lower, tolerance, value text)`` of a dimension object.

    With a matcher (a document read from STEP) the signed bounds, the fit class and the
    qualifier come from the file (:mod:`flatmat._pmi_step`); without one (a document built
    in memory) from XCAF.
    """
    if obj.IsDimWithPlusMinusTolerance():
        kept = Kept(obj.GetValue(), obj.GetUpperTolValue(), obj.GetLowerTolValue())
    elif obj.IsDimWithRange():  # OCCT's (nominal, |lower|) for a negative upper bound, or true limits
        kept = Kept(obj.GetLowerBound(), None, obj.GetUpperBound(), negative_upper=True)
    else:
        kept = Kept(obj.GetValue(), None, None, fit=obj.IsDimWithClassOfTolerance())
    found = matcher.match(kept) if matcher is not None else None  # every dimension takes its turn
    if matcher is None:
        qualifier = _xcaf_qualifier(obj)
    else:
        qualifier = _QUALIFIERS.get(found.qualifier, found.qualifier) if found is not None else ""
    nominal = kept.nominal
    base = f"{_num(nominal)}{unit}"
    q = f"{qualifier} " if qualifier else ""
    if kept.lower_magnitude is None:  # no plus-minus tolerance and no limits
        if matcher is not None and found is None:
            return nominal, None, None, "unreadable", f"{base} (not found in the file: {_RULED_OUT})"
        if kept.fit:
            return nominal, None, None, "fit", f"{q}{base} {_fit(found)}"
        if qualifier in ("max", "min"):  # an upper or a lower limit
            return None, None, None, "limits", f"{qualifier} {base}"
        if qualifier:
            return nominal, None, None, "unreadable", f"{qualifier} {base}"
        return nominal, None, None, "none", base
    if kept.negative_upper and (found is None or found.upper is None):  # limits, or a tolerance not found
        lo: float | None
        hi: float | None
        if found is not None and found.limits is not None:  # the file's own
            lo, hi = found.limits[0].value, found.limits[1].value
            if lo is None or hi is None:
                raw = f"{_num(found.limits[0].raw)}–{_num(found.limits[1].raw)}"
                return None, None, None, "unreadable", f"{raw} (limits, unit unknown)"
        else:
            lo, hi = nominal, kept.lower_magnitude
        limits = f"{_num(lo)}{unit}–{_num(hi)}{unit}"
        if lo <= hi:
            return None, None, None, "limits", f"{limits} (limits)"
        if matcher is None or found is not None:  # true limits, upper below lower
            return None, None, None, "unreadable", f"{limits} (limits inverted)"
        # inverted limits, or a tolerance with a negative upper bound that the file does not confirm
        return None, None, None, "unreadable", f"{limits} (tolerance unreadable: not found in the file)"
    if matcher is None:
        if kept.upper is None:  # in memory a range is limits, handled above
            return nominal, None, None, "unreadable", f"{base} (tolerance unreadable)"
        raw_up, raw_lo = kept.upper, 0.0 - kept.lower_magnitude
    elif found is None or found.upper is None or found.lower is None:
        # XCAF's magnitudes: the lower bound's sign is lost, and an angle's are in the file's unit
        up = _signed(kept.upper) if kept.upper is not None else "?"
        where = ", in the file's angle unit" if unit else ""
        return (
            nominal,
            None,
            None,
            "unreadable",
            f"{base} {up}/±{_num(kept.lower_magnitude)} (sign unknown: not found in the file{where})",
        )
    elif found.upper.value is None or found.lower.value is None:
        tol = f"{_signed(found.upper.raw)}/{_signed(found.lower.raw)}"
        return nominal, None, None, "unreadable", f"{q}{base} {tol} (unit unknown)"
    else:
        raw_up, raw_lo = found.upper.value, found.lower.value
    upper, lower = tolerance(nominal, raw_up, raw_lo)
    text = f"{q}{base} {_signed(raw_up)}/{_signed(raw_lo)}"
    if upper is None:
        return nominal, None, None, "unreadable", text + " (tolerance unreadable)"
    return nominal, upper, lower, "plus_minus", text


_RULED_OUT = "a fit or a max/min qualifier cannot be ruled out"


def _dimension(label, faces: _Faces, matcher: Matcher | None):
    """``(fields, graphical)``; fields is None when the item is on another solid."""
    obj = XCAFDoc_Dimension.Set_s(label).GetObject()
    dtype = obj.GetType()
    if dtype in _GRAPHICAL:
        return None, True
    kind = _DIM_KINDS.get(dtype, "other_gdt")
    unit = "°" if kind == "angle" else ""
    nominal, upper, lower, tol, value_text = _values(obj, matcher, unit)  # before the sides: each takes its turn
    first, second = _refs(label)
    pts = [obj.GetPoint().Coord() if obj.HasPoint() else None, obj.GetPoint2().Coord() if obj.HasPoint2() else None]
    sides = _sides(faces, [first, second], pts)
    if sides is None:
        return None, False
    name = _text(obj.GetSemanticName())
    if kind == "other_gdt":
        name = name or _enum_name(dtype, "XCAFDimTolObjects_DimensionType_")
    text = f"{_PREFIX.get(kind, '')}{value_text}"
    if kind == "other_gdt":
        text = f"{name} {text}"
    for i in range(1, obj.NbDescriptions() + 1):
        text += f"; {_text(obj.GetDescription(i))}"
    if _text(obj.GetPresentationName()):
        text += f"; {_text(obj.GetPresentationName())}"
    text_point = obj.GetPointTextAttach().Coord() if obj.HasTextPoint() else None
    value = None if nominal is None or math.isnan(nominal) else nominal
    if value is None and tol in ("none", "fit", "plus_minus"):
        tol = "unreadable"  # a NaN nominal
    fields = dict(kind=kind, nominal=value, upper=upper, lower=lower, sides=sides, text=text, tolerance=tol)
    return {**fields, "text_point": text_point}, False


def _datum_names(label) -> tuple[str, ...]:
    seq = LabelSequence()
    XCAFDoc_DimTolTool.GetDatumWithObjectOfTolerLabels_s(label, seq)
    named = []
    for k, lab in enumerate(_seq(seq)):
        obj = XCAFDoc_Datum.Set_s(lab).GetObject()
        named.append((obj.GetPosition() or k + 1, k, _text(obj.GetName())))
    return tuple(n for _, _, n in sorted(named))


def _tolerance_item(label, faces: _Faces):
    obj = XCAFDoc_GeomTolerance.Set_s(label).GetObject()
    ttype = obj.GetType()
    kind = _TOL_KINDS.get(ttype, "other_gdt")
    first, _ = _refs(label)
    sides = _sides(faces, [first], [])
    if sides is None:
        return None
    zone = obj.GetValue()
    datums = _datum_names(label)
    name = _enum_name(ttype, "XCAFDimTolObjects_GeomToleranceType_")
    dia = "Ø" if "Diameter" in str(obj.GetTypeOfValue()) else ""
    text = " | ".join([f"{name} {dia}{_num(zone)}", *datums])
    text_point = obj.GetPointTextAttach().Coord() if obj.HasPointText() else None
    upper = None if math.isnan(zone) else zone
    return dict(
        kind=kind,
        nominal=None,
        upper=upper,
        lower=None,
        datums=datums,
        sides=sides,
        text=text,
        text_point=text_point,
        tolerance="none" if upper is not None else "unreadable",
    )


def _datum_item(label, faces: _Faces):
    obj = XCAFDoc_Datum.Set_s(label).GetObject()
    first, _ = _refs(label)
    sides = _sides(faces, [first], [])
    if sides is None:
        return None
    name = _text(obj.GetName())
    text_point = obj.GetPointTextAttach().Coord() if obj.HasPointText() else None
    return dict(
        kind="datum",
        nominal=None,
        upper=None,
        lower=None,
        datums=(name,),
        sides=sides,
        text=name,
        text_point=text_point,
    )


def _notes(doc, faces: _Faces):
    """``(text, sides)`` per note, or ``(text, None)`` when it points at another solid."""
    tool = XCAFDoc_DocumentTool.NotesTool_s(doc.Main())
    notes = LabelSequence()
    tool.GetNotes(notes)
    if not notes.Length():
        return []
    targets: dict[str, list] = {}  # note entry -> annotated shapes
    annotated = LabelSequence()
    tool.GetAnnotatedItems(annotated)
    for lab in _seq(annotated):
        ref = XCAFDoc_AssemblyItemRef.Get_s(lab)
        if ref is None or ref.IsGUID() or ref.IsSubshapeIndex():
            continue  # notes on attributes, or on a sub-shape index: no face to name
        item_id = ref.GetItem()
        item = TDF_Label()
        TDF_Tool.Label_s(doc.GetData(), item_id.ToString(), item, False)
        if item.IsNull():
            continue
        found = LabelSequence()
        tool.GetNotes(item_id, found)
        for note in _seq(found):
            targets.setdefault(_entry(note), []).append(XCAFDoc_ShapeTool.GetShape_s(item))
    out = []
    for note in _seq(notes):
        comment = XCAFDoc_NoteComment.Get_s(note)
        if comment is not None:
            text = _text(comment.Comment())
        else:
            bin_data = XCAFDoc_NoteBinData.Get_s(note)
            text = _text(bin_data.Title()) if bin_data is not None else ""
        ids: set[int] = set()
        foreign = False
        for shape in targets.get(_entry(note), []):
            got = faces.of(shape)
            if got is None:
                foreign = True
                break
            ids.update(got)
        out.append((text, None if foreign else ((PmiSide(tuple(sorted(ids))),) if ids else ())))
    return out


def _entry(label) -> str:
    s = TCollection_AsciiString()
    TDF_Tool.Entry_s(label, s)
    return s.ToCString()


def read_pmi(doc, solid, step: Sequence[StepDimension] = (), *, trust_xcaf: bool = False) -> PmiData:
    """The semantic PMI of ``doc`` that refers to ``solid`` (see the module docstring).

    For a STEP file use :meth:`flatmat.StepModel.pmi_for`, which passes the file's
    dimensions. ``step`` is what :func:`flatmat._pmi_step.step_dimensions` read from the
    file the document came from, and tolerances come from there: without it every
    tolerance is unreadable. ``trust_xcaf`` takes XCAF's own values instead, which is
    right only for a document built in memory (OCCT's STEP reader loses signs and units).

    Never raises: any error gives status ``"failed"`` with a ``PMI_READ_FAILED`` warning.
    """
    try:
        return _read(doc, solid, step, trust_xcaf)
    except Exception as e:  # PMI never refuses a part
        return PmiData("failed", (), diagnostics=(warn(Code.PMI_READ_FAILED, f"PMI not read: {e!r}"),))


def _read(doc, solid, step: Sequence[StepDimension], trust_xcaf: bool) -> PmiData:
    faces = _Faces(solid)
    tool = XCAFDoc_DocumentTool.DimTolTool_s(doc.Main())
    dims, tols, datums = LabelSequence(), LabelSequence(), LabelSequence()
    tool.GetDimensionLabels(dims)
    tool.GetGeomToleranceLabels(tols)
    tool.GetDatumLabels(datums)
    semantic = sum(XCAFDoc_Dimension.Set_s(lab).GetObject().GetType() not in _GRAPHICAL for lab in _seq(dims))
    matcher = None if trust_xcaf else Matcher(step, semantic)
    labelled = sorted(
        [(lab.Tag(), "dim", lab) for lab in _seq(dims)]
        + [(lab.Tag(), "tol", lab) for lab in _seq(tols)]
        + [(lab.Tag(), "datum", lab) for lab in _seq(datums)],
        key=lambda t: t[0],
    )
    items: list[PmiItem] = []
    n_semantic = dropped = graphical = 0
    presented: list[frozenset[int] | None] = []  # the faces of each graphical annotation
    for _, what, lab in labelled:
        if what == "dim":
            fields, is_graphical = _dimension(lab, faces, matcher)
            if is_graphical:
                graphical += 1
                presented.append(_presented(lab, faces))
                continue
        else:
            fields = _tolerance_item(lab, faces) if what == "tol" else _datum_item(lab, faces)
        n_semantic += 1
        if fields is None:
            dropped += 1
            continue
        text_point = fields.pop("text_point")
        items.append(PmiItem(n_semantic, anchor=_anchor(text_point, fields["sides"], faces), **fields))
    for text, sides in _notes(doc, faces):
        n_semantic += 1
        if sides is None:
            dropped += 1
            continue
        items.append(
            PmiItem(n_semantic, "note", None, None, None, sides=sides, anchor=_anchor(None, sides, faces), text=text)
        )
    if n_semantic:
        semantic_faces = [{f for s in i.sides for f in s.face_ids} for i in items]
        orphans = sum(1 for k in presented if k is not None and not (k and any(k <= f for f in semantic_faces)))
        if not orphans:
            return PmiData("semantic", tuple(items), dropped)
        return PmiData("semantic", tuple(items), dropped, (_graphical_warning(orphans),), graphical=orphans)
    if graphical:
        return PmiData("graphical_only", (), diagnostics=(_graphical_warning(graphical),), graphical=graphical)
    return PmiData("none", ())


def _graphical_warning(n: int) -> Diagnostic:
    msg = f"{n} PMI annotation(s) are graphical only, without semantic data; they are not checked"
    return warn(Code.PMI_GRAPHICAL_ONLY, msg, count=n)


def _presented(label, faces: _Faces) -> frozenset[int] | None:
    """The face ids a graphical annotation refers to; None when they are on another solid.

    OCCT's STEP reader keeps only the first shape of a presentation, so an annotation
    repeats a semantic item when its faces are among the item's."""
    ids = faces.side([lab for labels in _refs(label) for lab in labels])
    return None if ids is None else frozenset(ids)
