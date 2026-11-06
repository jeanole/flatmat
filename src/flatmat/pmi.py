"""Semantic PMI (AP242 GD&T) from the XCAF document :func:`flatmat.load_step` keeps.

:func:`read_pmi` turns the document's dimensions, geometric tolerances, datums and
notes into :class:`PmiItem` s that name their faces by flatmat face id (the
``TopExp.MapShapes`` index on the solid, 1-based). Mapping the items onto flanges
and bends happens later (:mod:`flatmat.pmi_map`).

What OCCT's reader gives, and what flatmat makes of it (verified on OCP 8.0.1.1):

- Lengths come in mm and angles in degrees: the reader scales values, tolerances and
  points by the file's length unit and converts angular nominals to degrees. Angular
  *tolerances* are left in the file's angle unit; that is right for the usual
  files in degrees, and wrong for a file in radians.
- A plus-minus tolerance is stored as the upper deviation and the *magnitude* of the
  lower bound (``Abs``). flatmat gives both as signed deviations from the nominal:
  ``upper`` as read, ``lower = -magnitude``. A lower bound above the nominal
  therefore reads as below it; the reader has already lost the sign.
- A negative upper bound (an inverted pair, or both bounds below the nominal) makes
  the reader store two values (the nominal and the lower magnitude) as a "range" with
  the lower bound above the upper. Such a pair, and any NaN, is unreadable: ``upper``
  and ``lower`` are None and ``text`` keeps what was read.
- Graphical-only PMI (a presentation without semantic data) is read as dimensions
  of type ``DimensionPresentation`` or ``CommonLabel``. They never become items: a
  dimension written both ways is one item, and graphical PMI only sets the status
  ``"graphical_only"`` (with one ``PMI_GRAPHICAL_ONLY`` warning) when the file has no
  semantic PMI at all.
- OCCT's STEP reader does not fill the notes table; notes are read when the document
  has them (from another source).

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
from collections.abc import Iterable
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


@dataclass(frozen=True)
class PmiData:
    status: Literal["none", "semantic", "graphical_only", "failed"]
    items: tuple[PmiItem, ...]
    dropped: int = 0  # items whose faces belong to another solid
    diagnostics: tuple[Diagnostic, ...] = ()


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

    def __init__(self, solid):
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


def _dimension(label, faces: _Faces):
    """``(fields, graphical)``; fields is None when the item is on another solid."""
    obj = XCAFDoc_Dimension.Set_s(label).GetObject()
    dtype = obj.GetType()
    if dtype in _GRAPHICAL:
        return None, True
    kind = _DIM_KINDS.get(dtype, "other_gdt")
    first, second = _refs(label)
    pts = [obj.GetPoint().Coord() if obj.HasPoint() else None, obj.GetPoint2().Coord() if obj.HasPoint2() else None]
    sides = _sides(faces, [first, second], pts)
    if sides is None:
        return None, False
    unit = "°" if kind == "angle" else ""
    if obj.IsDimWithPlusMinusTolerance():
        nominal, raw_up, raw_lo = obj.GetValue(), obj.GetUpperTolValue(), 0.0 - obj.GetLowerTolValue()
        upper, lower = tolerance(nominal, raw_up, raw_lo)
        tol = f" +{_num(raw_up)}/{_num(raw_lo)}" if raw_lo else f" +{_num(raw_up)}/-0"
        tol += "" if upper is not None else " (tolerance unreadable)"
    elif obj.IsDimWithRange():
        lo_b, up_b = obj.GetLowerBound(), obj.GetUpperBound()
        if tolerance(None, up_b, lo_b) == (None, None):  # what the reader makes of an inverted pair
            nominal, upper, lower = lo_b, None, None
            tol = f" (tolerance unreadable: {_num(lo_b)}, {_num(up_b)})"
        else:
            nominal = (lo_b + up_b) / 2
            upper, lower = up_b - nominal, lo_b - nominal
            tol = f" ({_num(lo_b)}..{_num(up_b)})"
    else:
        nominal, upper, lower, tol = obj.GetValue(), None, None, ""
    name = _text(obj.GetSemanticName())
    if kind == "other_gdt":
        name = name or _enum_name(dtype, "XCAFDimTolObjects_DimensionType_")
    text = f"{_PREFIX.get(kind, '')}{_num(nominal)}{unit}{tol}"
    if kind == "other_gdt":
        text = f"{name} {text}"
    for i in range(1, obj.NbDescriptions() + 1):
        text += f"; {_text(obj.GetDescription(i))}"
    if math.isnan(nominal):
        nominal = None
    text_point = obj.GetPointTextAttach().Coord() if obj.HasTextPoint() else None
    return (
        dict(kind=kind, nominal=nominal, upper=upper, lower=lower, sides=sides, text=text, text_point=text_point),
        False,
    )


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
        kind=kind, nominal=None, upper=upper, lower=None, datums=datums, sides=sides, text=text, text_point=text_point
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


def read_pmi(doc, solid) -> PmiData:
    """The semantic PMI of ``doc`` that refers to ``solid`` (see the module docstring).

    Never raises: any error gives status ``"failed"`` with a ``PMI_READ_FAILED`` warning.
    """
    try:
        return _read(doc, solid)
    except Exception as e:  # PMI never refuses a part
        return PmiData("failed", (), diagnostics=(warn(Code.PMI_READ_FAILED, f"PMI not read: {e!r}"),))


def _read(doc, solid) -> PmiData:
    faces = _Faces(solid)
    tool = XCAFDoc_DocumentTool.DimTolTool_s(doc.Main())
    dims, tols, datums = LabelSequence(), LabelSequence(), LabelSequence()
    tool.GetDimensionLabels(dims)
    tool.GetGeomToleranceLabels(tols)
    tool.GetDatumLabels(datums)
    labelled = sorted(
        [(lab.Tag(), "dim", lab) for lab in _seq(dims)]
        + [(lab.Tag(), "tol", lab) for lab in _seq(tols)]
        + [(lab.Tag(), "datum", lab) for lab in _seq(datums)],
        key=lambda t: t[0],
    )
    items: list[PmiItem] = []
    n_semantic = dropped = graphical = 0
    for _, what, lab in labelled:
        if what == "dim":
            fields, is_graphical = _dimension(lab, faces)
            if is_graphical:
                graphical += 1
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
        return PmiData("semantic", tuple(items), dropped)
    if graphical:
        msg = f"{graphical} PMI annotation(s) are graphical only, without semantic data; they are not checked"
        return PmiData("graphical_only", (), diagnostics=(warn(Code.PMI_GRAPHICAL_ONLY, msg, count=graphical),))
    return PmiData("none", ())
