import math

import numpy as np
import pytest
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCP.BRepGProp import BRepGProp
from OCP.Geom2d import Geom2d_BSplineCurve
from OCP.GeomAPI import GeomAPI
from OCP.gp import gp_Ax3, gp_Pln, gp_Pnt, gp_Pnt2d
from OCP.GProp import GProp_GProps

try:  # OCP 7.8
    from OCP.TColgp import TColgp_Array1OfPnt2d
    from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
except ImportError:  # OCP 8.0
    from OCP.collections import Array1_double as TColStd_Array1OfReal
    from OCP.collections import Array1_gp_Pnt2d as TColgp_Array1OfPnt2d
    from OCP.collections import Array1_int as TColStd_Array1OfInteger

from flatmat.model import Line, Loop2D, Spline


def _occ(poles, weights, knots, mults, deg):
    P = TColgp_Array1OfPnt2d(1, len(poles))
    W = TColStd_Array1OfReal(1, len(poles))
    K = TColStd_Array1OfReal(1, len(knots))
    M = TColStd_Array1OfInteger(1, len(mults))
    for i, (x, y) in enumerate(poles):
        P.SetValue(i + 1, gp_Pnt2d(x, y))
        W.SetValue(i + 1, weights[i])
    for i, k in enumerate(knots):
        K.SetValue(i + 1, k)
    for i, m in enumerate(mults):
        M.SetValue(i + 1, m)
    return Geom2d_BSplineCurve(P, W, K, M, deg)


@pytest.mark.parametrize("rational", [False, True])
def test_de_boor_matches_occt(rational):
    poles = [(0, 0), (1, 2), (3, 3), (4, 0), (6, 1), (7, 4)]
    weights = [1, 0.7, 1.3, 1, 0.5, 1] if rational else [1] * 6
    knots, mults, deg = [0.0, 0.3, 0.7, 1.0], [4, 1, 1, 4], 3
    ref = _occ(poles, weights, knots, mults, deg)
    s = Spline(deg, tuple(poles), tuple(knots), tuple(mults), tuple(weights) if rational else None)
    for i in range(41):
        u = i / 40
        p = ref.Value(u)
        assert s.point(u) == pytest.approx((p.X(), p.Y()), abs=1e-12)
    r = s.reversed()
    assert r.point(0.0) == pytest.approx(s.point(1.0))
    assert r.point(0.25) == pytest.approx(s.point(0.75))


@pytest.mark.parametrize("rational", [False, True])
def test_area_matches_occt(rational):
    # a closed loop: the spline from (0, 0) to (7, 4), back along two lines
    poles = [(0, 0), (1, 2), (3, 3), (4, 0), (6, 1), (7, 4)]
    weights = [1, 0.7, 1.3, 1, 0.5, 1] if rational else [1] * 6
    knots, mults, deg = [0.0, 0.3, 0.7, 1.0], [4, 1, 1, 4], 3
    s = Spline(deg, tuple(poles), tuple(knots), tuple(mults), tuple(weights) if rational else None)
    loop = Loop2D((s, Line((7, 4), (7, -1)), Line((7, -1), (0, 0))))
    pln = gp_Pln(gp_Ax3())
    curve = GeomAPI.To3d_s(_occ(poles, weights, knots, mults, deg), pln)
    wire = BRepBuilderAPI_MakeWire(
        BRepBuilderAPI_MakeEdge(curve).Edge(),
        BRepBuilderAPI_MakeEdge(gp_Pnt(7, 4, 0), gp_Pnt(7, -1, 0)).Edge(),
        BRepBuilderAPI_MakeEdge(gp_Pnt(7, -1, 0), gp_Pnt(0, 0, 0)).Edge(),
    ).Wire()
    g = GProp_GProps()
    BRepGProp.SurfaceProperties_s(BRepBuilderAPI_MakeFace(pln, wire).Face(), g, 1e-12)
    assert abs(loop.area()) == pytest.approx(g.Mass(), rel=1e-10)
    assert s.start == s.point(0.0) and s.end == pytest.approx(s.point(1.0), abs=1e-15)


def test_area_of_a_rational_circle():
    h = 0.5**0.5
    poles = ((1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (1, 0))
    s = Spline(2, poles, (0, 1, 2, 3, 4), (3, 2, 2, 2, 3), (1, h, 1, h, 1, h, 1, h, 1))
    assert Loop2D((s,)).area() == pytest.approx(math.pi, abs=1e-13)


@pytest.mark.parametrize("rational", [False, True])
def test_vectorised_points_are_point_to_the_last_bit(rational):
    poles = ((0, 0), (1, 2), (3, 3), (4, 0), (6, 1), (7, 4))
    weights = (1, 0.7, 1.3, 1, 0.5, 1) if rational else None
    s = Spline(3, poles, (0.0, 0.3, 0.7, 1.0), (4, 1, 1, 4), weights)
    us = [i / 333 for i in range(334)] + [0.3, 0.7, -0.1, 1.1]  # knots, and outside the range
    assert [tuple(p) for p in s.points(np.asarray(us)).tolist()] == [s.point(u) for u in us]


def test_polyline_is_sampled_once_and_handed_out_as_copies():
    s = Spline(2, ((0, 0), (1, 1), (2, 0)), (0.0, 1.0), (3, 3))
    loop = Loop2D((s, Line((2, 0), (0, 0))))
    a = loop.polyline(0.01)
    a.append((99.0, 99.0))
    assert loop.polyline(0.01) == a[:-1] and loop.polyline(0.01) is not loop.polyline(0.01)
    assert len(loop.polyline(0.001)) > len(a) - 1  # each tolerance on its own
