import pytest
from OCP.Geom2d import Geom2d_BSplineCurve
from OCP.gp import gp_Pnt2d

try:  # OCP 7.8
    from OCP.TColgp import TColgp_Array1OfPnt2d
    from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
except ImportError:  # OCP 8.0
    from OCP.collections import Array1_double as TColStd_Array1OfReal
    from OCP.collections import Array1_gp_Pnt2d as TColgp_Array1OfPnt2d
    from OCP.collections import Array1_int as TColStd_Array1OfInteger

from flatmat.model import Spline


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
