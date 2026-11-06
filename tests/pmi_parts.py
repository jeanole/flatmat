"""AP242 test files with semantic PMI, written with OCCT's own XCAF writer.

:func:`write_step_with_pmi` takes a shape and a list of items, each a dict:

``kind``
    ``"linear"``, ``"diameter"``, ``"radius"``, ``"angle"`` (a dimension); a GD&T type,
    ``"position"``, ``"flatness"``, ``"profile_of_surface"`` ... (any
    ``XCAFDimTolObjects_GeomToleranceType`` in snake case); ``"datum"``; ``"graphical"``
    (a presentation without semantic data); ``"note"`` (in-memory documents only:
    OCCT's STEP writer drops notes).
``faces`` / ``edges``
    One list per side of face (edge) ids, the ``TopExp.MapShapes`` index on ``shape``.
``nominal``, ``upper``, ``lower``
    A dimension's value and signed deviations, **in the file's unit** (what the drawing
    says: ``nominal=1.0, upper=0.005, lower=-0.005, unit="INCH"``).
``zone``, ``datums``
    A GD&T item's zone width (file unit) and datum names; a ``datum`` item has ``name``.
``points``, ``presentation``, ``text``
    The dimension's own point per side (mm, model space); ``presentation=True`` adds a
    graphical presentation (a line between the sides) with its text point; ``text`` is
    a note's text.

OCCT 8.0.1.1's writer has three quirks this module works around, so that the file says
what the item says:

- the AP242 schema must be set after the STEP controller is initialised (a first writer);
- dimension values and PMI points are written as they are, labelled with the file's
  unit (no conversion from the document's mm), so values are given in the file unit and
  points are divided by the unit's size here;
- a geometric tolerance in a file without a semantic dimension is labelled in metres;
  each file is read back and the zones are rewritten in metres when they came back
  1000 times too large.

A text point is only written with a presentation: without one, the writer crashes on
a document with GD&T in it. Presentations are for dimensions only: one on a GD&T item
in a file without a dimension makes the write fail.
"""

from __future__ import annotations

import math

from OCP.BRep import BRep_Builder
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeVertex
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepGProp import BRepGProp
from OCP.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane
from OCP.gp import gp_Pnt
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.STEPCAFControl import STEPCAFControl_Reader, STEPCAFControl_Writer
from OCP.TCollection import TCollection_ExtendedString, TCollection_HAsciiString
from OCP.TDocStd import TDocStd_Document
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE
from OCP.TopExp import TopExp
from OCP.TopoDS import TopoDS_Compound
from OCP.XCAFApp import XCAFApp_Application
from OCP.XCAFDimTolObjects import (
    XCAFDimTolObjects_DatumObject,
    XCAFDimTolObjects_DimensionObject,
    XCAFDimTolObjects_GeomToleranceObject,
)
from OCP.XCAFDimTolObjects import XCAFDimTolObjects_DimensionType as DimType
from OCP.XCAFDimTolObjects import XCAFDimTolObjects_GeomToleranceType as TolType
from OCP.XCAFDoc import (
    XCAFDoc_Datum,
    XCAFDoc_Dimension,
    XCAFDoc_DocumentTool,
    XCAFDoc_GeomTolerance,
)

from flatmat._occ import LabelSequence, ShapeMap, as_face

UNITS = {"MM": 1.0, "CM": 10.0, "M": 1000.0, "INCH": 25.4}
_DIMS = {
    "linear": DimType.XCAFDimTolObjects_DimensionType_Location_LinearDistance,
    "diameter": DimType.XCAFDimTolObjects_DimensionType_Size_Diameter,
    "radius": DimType.XCAFDimTolObjects_DimensionType_Size_Radius,
    "angle": DimType.XCAFDimTolObjects_DimensionType_Location_Angular,
    "graphical": DimType.XCAFDimTolObjects_DimensionType_DimensionPresentation,
}
SURFACES = {"plane": GeomAbs_Plane, "cylinder": GeomAbs_Cylinder}


def _tol_type(kind: str):
    camel = "".join(w.capitalize() for w in kind.split("_"))
    for name in dir(TolType):
        if name.lower() == f"xcafdimtolobjects_geomtolerancetype_{camel}".lower():
            return getattr(TolType, name)
    raise ValueError(f"unknown GD&T kind {kind!r}")


def _map(shape, kind=TopAbs_FACE):
    m = ShapeMap()
    TopExp.MapShapes_s(shape, kind, m)
    return m


def _at(shape, kind, point, tol: float = 1e-6) -> list[int]:
    m, v = _map(shape, kind), BRepBuilderAPI_MakeVertex(gp_Pnt(*point)).Vertex()
    return [i for i in range(1, m.Size() + 1) if BRepExtrema_DistShapeShape(m.FindKey(i), v).Value() < tol]


def face_at(shape, point, surface: str | None = None) -> int:
    """The id of the one face of ``shape`` through ``point`` (of the given surface type)."""
    m = _map(shape)
    found = [
        i
        for i in _at(shape, TopAbs_FACE, point)
        if surface is None or BRepAdaptor_Surface(as_face(m.FindKey(i))).GetType() == SURFACES[surface]
    ]
    if len(found) != 1:
        raise ValueError(f"{len(found)} faces through {point}")
    return found[0]


def edge_at(shape, point) -> int:
    """The id of the one edge of ``shape`` through ``point`` (``TopExp.MapShapes`` edge index)."""
    found = _at(shape, TopAbs_EDGE, point)
    if len(found) != 1:
        raise ValueError(f"{len(found)} edges through {point}")
    return found[0]


def centroid(shape, face_ids) -> tuple[float, float, float]:
    m, g = _map(shape), GProp_GProps()
    for i in face_ids:
        one = GProp_GProps()
        BRepGProp.SurfaceProperties_s(as_face(m.FindKey(i)), one)
        g.Add(one)
    return g.CentreOfMass().Coord()


def _seq(labels):
    s = LabelSequence()
    for lab in labels:
        s.Append(lab)
    return s


def make_doc_with_pmi(shape, items: list[dict], unit: str = "MM", *, zone_scale: float = 1.0):
    """An XCAF document holding ``shape`` and ``items`` (values as the file is to state them)."""
    k = UNITS[unit]
    doc = TDocStd_Document(TCollection_ExtendedString("XmlOcaf"))
    XCAFApp_Application.GetApplication_s().InitDocument(doc)
    shapes = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())
    tool = XCAFDoc_DocumentTool.DimTolTool_s(doc.Main())
    root = shapes.AddShape(shape, False)
    fmap, emap = _map(shape), _map(shape, TopAbs_EDGE)

    def labels(it) -> list[list]:
        sides = [[shapes.AddSubShape(root, fmap.FindKey(i)) for i in side] for side in it.get("faces", ())]
        for j, side in enumerate(it.get("edges", ())):
            while len(sides) <= j:
                sides.append([])
            sides[j] += [shapes.AddSubShape(root, emap.FindKey(i)) for i in side]
        return sides

    def scaled(p):
        return gp_Pnt(*(c / k for c in p))

    def where(it) -> list:
        pts = list(it.get("points", ()))
        for j in range(len(pts), len(it.get("faces", ()))):
            pts.append(centroid(shape, it["faces"][j]))
        return pts

    def presentation(pts):
        b, c = BRep_Builder(), TopoDS_Compound()
        b.MakeCompound(c)
        a, z = pts[0], pts[-1] if len(pts) > 1 else (pts[0][0] + 10, pts[0][1], pts[0][2])
        b.Add(c, BRepBuilderAPI_MakeEdge(scaled(a), scaled(z)).Edge())
        return c, tuple((p + q) / 2 for p, q in zip(a, z, strict=True))

    datums: dict[str, object] = {}
    links = []
    notes = []
    for n, it in enumerate(items, 1):
        kind = it["kind"]
        sides = labels(it)
        if kind in _DIMS:
            lab = tool.AddDimension()
            tool.SetDimension(_seq(sides[0]), _seq(sides[1] if len(sides) > 1 else []), lab)
            obj = XCAFDimTolObjects_DimensionObject()
            obj.SetType(_DIMS[kind])
            if kind != "graphical":
                obj.SetValue(it["nominal"])
                if it.get("upper") is not None:
                    # OCCT keeps the lower deviation as a magnitude and writes lower_bound = -magnitude
                    obj.SetUpperTolValue(it["upper"])
                    obj.SetLowerTolValue(0.0 - it["lower"])
            for j, p in enumerate(it.get("points", ())):
                (obj.SetPoint if j == 0 else obj.SetPoint2)(scaled(p))
            if it.get("presentation") or kind == "graphical":
                pres, text_point = presentation(where(it))
                obj.SetPresentation(pres, TCollection_HAsciiString(f"PMI {n}"))
                obj.SetPointTextAttach(scaled(text_point))
            XCAFDoc_Dimension.Set_s(lab).SetObject(obj)
        elif kind == "datum":
            lab = tool.AddDatum()
            tool.SetDatum(_seq(sides[0]), lab)
            datums[it["name"]] = (lab, sides)
        elif kind == "note":
            notes.append((it["text"], sides))
        else:
            lab = tool.AddGeomTolerance()
            tool.SetGeomTolerance(_seq(sides[0]), lab)
            obj = XCAFDimTolObjects_GeomToleranceObject()
            obj.SetType(_tol_type(kind))
            obj.SetValue(it["zone"] * zone_scale)
            if it.get("presentation"):
                raise ValueError("presentations are for dimensions only (OCCT fails to write one on GD&T alone)")
            XCAFDoc_GeomTolerance.Set_s(lab).SetObject(obj)
            links += [(name, lab) for name in it.get("datums", ())]
    order = {name: i for i, (name, _) in enumerate(links, 1)}
    for name, (lab, _) in datums.items():
        obj = XCAFDimTolObjects_DatumObject()
        obj.SetName(TCollection_HAsciiString(name))
        obj.SetPosition(order.get(name, 1))
        XCAFDoc_Datum.Set_s(lab).SetObject(obj)
    for name, tol in links:
        tool.SetDatumToGeomTol(datums[name][0], tol)
    if notes:
        nt = XCAFDoc_DocumentTool.NotesTool_s(doc.Main())
        for text, sides in notes:
            note = nt.CreateComment(
                TCollection_ExtendedString("flatmat"), TCollection_ExtendedString(""), TCollection_ExtendedString(text)
            )
            for lab in [lab for side in sides for lab in side] or [root]:
                nt.AddNote(note.Label(), lab)
    return doc


def _zones(path: str) -> list[float]:
    reader = STEPCAFControl_Reader()
    reader.SetGDTMode(True)
    reader.ReadFile(path)
    doc = TDocStd_Document(TCollection_ExtendedString("XmlOcaf"))
    XCAFApp_Application.GetApplication_s().InitDocument(doc)
    reader.Transfer(doc)
    seq = LabelSequence()
    XCAFDoc_DocumentTool.DimTolTool_s(doc.Main()).GetGeomToleranceLabels(seq)
    return [XCAFDoc_GeomTolerance.Set_s(seq.Value(i)).GetObject().GetValue() for i in range(1, seq.Length() + 1)]


def _write(doc, path: str, unit: str) -> None:
    STEPCAFControl_Writer()  # initialises the STEP statics
    old_schema, old_unit = Interface_Static.IVal_s("write.step.schema"), Interface_Static.CVal_s("write.step.unit")
    Interface_Static.SetIVal_s("write.step.schema", 5)  # AP242
    Interface_Static.SetCVal_s("write.step.unit", unit)
    try:
        writer = STEPCAFControl_Writer()
        writer.SetDimTolMode(True)
        if not writer.Transfer(doc) or writer.Write(path) != IFSelect_RetDone:
            raise RuntimeError(f"AP242 write failed: {path}")
    finally:
        Interface_Static.SetIVal_s("write.step.schema", old_schema)
        Interface_Static.SetCVal_s("write.step.unit", old_unit)


def write_step_with_pmi(shape, items: list[dict], path: str, unit: str = "MM") -> None:
    """Write ``shape`` with ``items`` (see the module docstring) as AP242 in ``unit``."""
    if any(it["kind"] == "note" for it in items):
        raise ValueError("OCCT's STEP writer drops XCAF notes: use make_doc_with_pmi and read_pmi")
    k = UNITS[unit]
    want = [it["zone"] * k for it in items if it["kind"] not in (*_DIMS, "datum")]
    _write(make_doc_with_pmi(shape, items, unit), path, unit)
    if want and not _close(_zones(path), want):
        _write(make_doc_with_pmi(shape, items, unit, zone_scale=k / 1000.0), path, unit)  # labelled in metres
        if not _close(_zones(path), want):
            raise RuntimeError(f"GD&T zones read back as {_zones(path)}, not {want}")


def _close(got: list[float], want: list[float]) -> bool:
    return len(got) == len(want) and all(math.isclose(g, w, rel_tol=1e-9) for g, w in zip(got, want, strict=True))
