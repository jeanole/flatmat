"""The recognition context: OCP faces and edges as recognize() leaves them for unfold().

Internal to flatmat; nothing here is public API.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve2d
from OCP.BRepLProp import BRepLProp_SLProps
from OCP.gp import gp_Pnt2d, gp_Vec
from OCP.IntTools import IntTools_FClass2d
from OCP.TopAbs import TopAbs_IN, TopAbs_REVERSED, TopAbs_VERTEX
from OCP.TopExp import TopExp_Explorer

from ._occ import as_vertex
from .diagnostics import Code, refuse
from .model import Tolerances


@dataclass
class Face:
    idx: int
    shape: Any
    kind: str
    area: float
    normal: gp_Vec  # outward (face orientation applied), at the UV centre
    mirror: float  # +1 if the surface normal equals the outward normal, else -1
    uv: tuple[float, float, float, float]
    ad: Any  # BRepAdaptor_Surface, one per face
    radius: float = 0.0
    axis_p: Any = None  # gp_Pnt for cylinders
    axis_d: Any = None  # gp_Vec for cylinders and cones
    semi_angle: float = 0.0  # cones
    _cls: Any = None

    @property
    def uv_centre(self) -> tuple[float, float]:
        u0, u1, v0, v1 = self.uv
        return 0.5 * (u0 + u1), 0.5 * (v0 + v1)

    def inside(self, u: float, v: float) -> bool:
        # one classifier per face: BRepClass_FaceClassifier rebuilds per call (~9 ms on
        # a face with 1000 holes against 0.8 ms here)
        if self._cls is None:
            self._cls = IntTools_FClass2d(self.shape, 1e-7)
        return self._cls.Perform(gp_Pnt2d(u, v)) == TopAbs_IN


@dataclass(frozen=True)
class BendFace:
    face: int
    bend_id: int
    r_in: float
    concave: bool  # reference skin is the inside of the bend at this face
    angle: float  # rad


@dataclass(frozen=True)
class SharpEdge:
    """A sharp bend edge between two faces of the reference skin."""

    edge: int
    faces: tuple[int, int]
    angle: float  # deflection, rad
    concave: bool  # reference skin is the inside of the bend here
    length: float
    vertices: tuple[int, int]  # vertex map indices
    bend_id: int = 0


@dataclass(frozen=True)
class Context:
    """Everything unfold() needs; kept out of the public model."""

    solid: Any
    fmap: Any
    vmap: Any
    faces: list[Face]
    edges: list[tuple[Any, tuple[int, ...]]]
    adj: dict[int, list[tuple[int, int]]]  # face -> [(neighbour, edge index)] over SMOOTH edges
    skin_a: set[int]
    root: int
    thickness: float
    bend_faces: dict[int, BendFace] = field(default_factory=dict)
    sharp_edges: dict[int, SharpEdge] = field(default_factory=dict)
    features: dict[int, tuple[int, ...]] = field(default_factory=dict)  # feature id -> faces
    feature_through_r: dict[int, float] = field(default_factory=dict)
    feature_on_a: dict[int, bool] = field(default_factory=dict)
    tol: Tolerances = field(default_factory=Tolerances)


def surface_props(face, ad, u: float, v: float) -> tuple[gp_Vec, float]:
    """Outward normal at (u, v) and the face's mirror sign (see :attr:`Face.mirror`)."""
    pr = BRepLProp_SLProps(ad, u, v, 1, 1e-7)
    n = gp_Vec(pr.Normal())
    sn = pr.D1U().Crossed(pr.D1V())
    if face.Orientation() == TopAbs_REVERSED:
        n.Reverse()
    return n, (1.0 if sn.Dot(n) > 0 else -1.0)


def pcurve(e, f) -> BRepAdaptor_Curve2d:
    """BRepAdaptor_Curve2d, refusing instead of crashing when the pcurve is missing.

    On a missing pcurve the adaptor raises on OCP 7.8 and segfaults on 8.0.
    """
    if BRep_Tool.CurveOnSurface_s(e, f, 0.0, 0.0) is None:
        raise refuse(Code.MISSING_PCURVE, "an edge has no pcurve on its face")
    return BRepAdaptor_Curve2d(e, f)


def edge_vertices(e) -> list:
    """The vertices of an edge, first to last (one for a closed edge, possibly twice)."""
    out = []
    ex = TopExp_Explorer(e, TopAbs_VERTEX)
    while ex.More():
        out.append(as_vertex(ex.Current()))
        ex.Next()
    return out


def inner_uv(f: Face) -> tuple[float, float] | None:
    """A (u, v) inside the face, from a coarse grid over its UV box; None if none is."""
    u0, u1, v0, v1 = f.uv
    for su in (0.5, 0.25, 0.75, 0.1, 0.9):
        for sv in (0.5, 0.25, 0.75, 0.1, 0.9):
            u, v = u0 + su * (u1 - u0), v0 + sv * (v1 - v0)
            if f.inside(u, v):
                return u, v
    return None
