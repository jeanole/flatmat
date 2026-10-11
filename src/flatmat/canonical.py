"""Canonical recognition of B-spline geometry (docs/design.md §5.5).

Exporters often write planes, cylinders, lines and circles as B-splines (some write
*everything* that way). flatmat needs analytic planes and cylinders to recognise
bends, and lines and circles to write clean DXF, so this module converts them back
before recognition:

- edges: ``ShapeAnalysis_CanonicalRecognition`` -> ``Geom_Line`` / ``Geom_Circle``
  with the parameter direction of the original curve;
- faces: ``GeomConvert_SurfToAnaSurf`` per face (never the shell mode, which in
  OCCT 7.8.0/8.0.0 compares the first face with itself) -> plane / cylinder / cone,
  with the normal direction of the original surface;
- the old pcurves are removed and rebuilt on the new geometry.

Works on a copy; the input shape is not modified.
"""

from __future__ import annotations

import math
from typing import Any

from OCP.BRep import BRep_Builder, BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy
from OCP.BRepTools import BRepTools
from OCP.ElCLib import ElCLib
from OCP.Geom import Geom_Circle, Geom_Line
from OCP.GeomAbs import (
    GeomAbs_BezierCurve,
    GeomAbs_BezierSurface,
    GeomAbs_BSplineCurve,
    GeomAbs_BSplineSurface,
    GeomAbs_Cone,
    GeomAbs_Cylinder,
    GeomAbs_OffsetSurface,
    GeomAbs_OtherCurve,
    GeomAbs_OtherSurface,
    GeomAbs_Plane,
    GeomAbs_SurfaceOfExtrusion,
    GeomAbs_SurfaceOfRevolution,
)
from OCP.GeomAdaptor import GeomAdaptor_Surface
from OCP.GeomConvert import GeomConvert_SurfToAnaSurf
from OCP.GeomLProp import GeomLProp_SLProps
from OCP.gp import gp_Ax2, gp_Circ, gp_Lin
from OCP.ShapeAnalysis import ShapeAnalysis_CanonicalRecognition, ShapeAnalysis_Surface
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_FORWARD, TopAbs_REVERSED, TopAbs_VERTEX
from OCP.TopExp import TopExp, TopExp_Explorer
from OCP.TopLoc import TopLoc_Location

from ._occ import EdgeFaceMap, as_edge, as_face, as_vertex
from .diagnostics import Code, Diagnostic, warn
from .io import add_missing_pcurves

_FREEFORM_SURFACES = {
    GeomAbs_BSplineSurface,
    GeomAbs_BezierSurface,
    GeomAbs_SurfaceOfExtrusion,
    GeomAbs_SurfaceOfRevolution,
    GeomAbs_OffsetSurface,
    GeomAbs_OtherSurface,
}
_FREEFORM_CURVES = {GeomAbs_BSplineCurve, GeomAbs_BezierCurve, GeomAbs_OtherCurve}
_TARGETS = {GeomAbs_Plane: "plane", GeomAbs_Cylinder: "cylinder", GeomAbs_Cone: "cone"}


def _faces(shape) -> list:
    out = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        out.append(as_face(ex.Current()))
        ex.Next()
    return out


def _edges(shape) -> list:
    out = []
    ex = TopExp_Explorer(shape, TopAbs_EDGE)
    while ex.More():
        out.append(as_edge(ex.Current()))
        ex.Next()
    return out


def needs_canonicalization(shape) -> bool:
    if any(BRepAdaptor_Surface(f, False).GetType() in _FREEFORM_SURFACES for f in _faces(shape)):
        return True
    return any(
        not BRep_Tool.Degenerated_s(e) and BRepAdaptor_Curve(e).GetType() in _FREEFORM_CURVES for e in _edges(shape)
    )


def _ends(e) -> tuple[Any, Any]:
    """First and last vertex in the parameter direction of the edge's curve."""
    first = last = None
    ex = TopExp_Explorer(e.Oriented(TopAbs_FORWARD), TopAbs_VERTEX)
    while ex.More():
        v = as_vertex(ex.Current())
        if v.Orientation() == TopAbs_FORWARD:
            first = v
        elif v.Orientation() == TopAbs_REVERSED:
            last = v
        ex.Next()
    return first, last


def _drop_pcurves(builder, e, faces) -> None:
    for f in faces:
        if BRep_Tool.IsClosed_s(e, f):
            builder.UpdateEdge(e, None, None, f, BRep_Tool.Tolerance_s(e))
        else:
            builder.UpdateEdge(e, None, f, BRep_Tool.Tolerance_s(e))


def _canonical_edge(e, tol: float):
    """(Geom curve, first, last) for an analytic replacement of e, or None."""
    c = BRepAdaptor_Curve(e)
    a, b = c.FirstParameter(), c.LastParameter()
    v0, v1 = _ends(e)
    if v0 is None or v1 is None:
        return None
    p0, p1 = BRep_Tool.Pnt_s(v0), BRep_Tool.Pnt_s(v1)
    pm = c.Value(0.5 * (a + b))
    rec = ShapeAnalysis_CanonicalRecognition(e)
    lin = gp_Lin()
    if rec.IsLine(tol, lin):
        t0, t1 = ElCLib.Parameter_s(lin, p0), ElCLib.Parameter_s(lin, p1)
        if t1 < t0:
            lin.Reverse()
            t0, t1 = -t0, -t1
        return Geom_Line(lin), t0, t1
    circ = gp_Circ()
    if rec.IsCircle(tol, circ):
        for attempt in range(2):
            t0 = ElCLib.Parameter_s(circ, p0)
            tm = (ElCLib.Parameter_s(circ, pm) - t0) % (2 * math.pi)
            t1 = (ElCLib.Parameter_s(circ, p1) - t0) % (2 * math.pi)
            if v0.IsSame(v1) or t1 < 1e-9:
                t1 = 2 * math.pi
            if tm < t1:
                return Geom_Circle(circ), t0, t0 + t1
            if attempt == 0:  # runs the other way round: flip the circle's axis
                ax = circ.Position()
                circ = gp_Circ(gp_Ax2(ax.Location(), ax.Direction().Reversed(), ax.XDirection()), circ.Radius())
    return None


def _canonical_surface(f, tol: float):
    """(analytic Geom surface with the original normal direction, location, gap) or None."""
    loc = TopLoc_Location()
    surf = BRep_Tool.Surface_s(f, loc)
    u0, u1, v0, v1 = BRepTools.UVBounds_s(f)
    res, conv = None, None
    for bounded in (True, False):
        conv = GeomConvert_SurfToAnaSurf(surf)
        res = conv.ConvertToAnalytical(tol, u0, u1, v0, v1) if bounded else conv.ConvertToAnalytical(tol)
        if res is not None and GeomAdaptor_Surface(res).GetType() in _TARGETS:
            break
        res = None
    if res is None or conv is None:
        return None
    uc, vc = 0.5 * (u0 + u1), 0.5 * (v0 + v1)
    p = surf.Value(uc, vc)
    n0 = GeomLProp_SLProps(surf, uc, vc, 1, 1e-9).Normal()
    uv = ShapeAnalysis_Surface(res).ValueOfUV(p, 1e-7)
    n1 = GeomLProp_SLProps(res, uv.X(), uv.Y(), 1, 1e-9).Normal()
    if n0.Dot(n1) < 0:
        res.UReverse()
    return res, loc, conv.Gap()


def canonicalize(shape, tol: float = 1e-4) -> tuple[Any, Diagnostic | None]:
    """Return a copy of shape with freeform geometry that is really analytic made analytic.

    The second value is a ``SURFACE_CONVERTED`` warning (with counts and the largest
    deviation) when anything was converted, else None. Geometry that is not analytic
    within ``tol`` is left alone (and refused later if it matters).
    """
    if not needs_canonicalization(shape):
        return shape, None
    shape = BRepBuilderAPI_Copy(shape, True).Shape()
    builder = BRep_Builder()
    efm = EdgeFaceMap()
    TopExp.MapShapesAndAncestors_s(shape, TopAbs_EDGE, TopAbs_FACE, efm)

    def faces_of(e) -> list:
        lst = efm.FindFromIndex(efm.FindIndex(e))
        return [as_face(lst.First()), as_face(lst.Last())] if lst.Size() <= 2 else [as_face(x) for x in lst]

    n_edges = n_faces = 0
    worst = 0.0
    for e in _edges(shape):
        if BRep_Tool.Degenerated_s(e) or BRepAdaptor_Curve(e).GetType() not in _FREEFORM_CURVES:
            continue
        found = _canonical_edge(e, tol)
        if found is None:
            continue
        curve, t0, t1 = found
        _drop_pcurves(builder, e, faces_of(e))
        builder.UpdateEdge(e, curve, TopLoc_Location(), max(BRep_Tool.Tolerance_s(e), tol))
        builder.Range(e, t0, t1)
        builder.SameRange(e, False)
        builder.SameParameter(e, False)
        n_edges += 1
    for f in _faces(shape):
        if BRepAdaptor_Surface(f, False).GetType() not in _FREEFORM_SURFACES:
            continue
        found = _canonical_surface(f, tol)
        if found is None:
            continue
        surf, loc, gap = found
        for e in _edges(f):
            _drop_pcurves(builder, e, [f])
        builder.UpdateFace(f, surf, loc, max(BRep_Tool.Tolerance_s(f), gap))
        worst = max(worst, gap)
        n_faces += 1
    if not (n_edges or n_faces):
        return shape, None
    add_missing_pcurves(shape)
    return shape, warn(
        Code.SURFACE_CONVERTED,
        f"{n_faces} faces and {n_edges} edges were B-splines and are analytic within {worst:.2g} mm",
        faces_converted=n_faces,
        edges_converted=n_edges,
        max_gap=worst,
    )


__all__ = ["canonicalize", "needs_canonicalization"]
