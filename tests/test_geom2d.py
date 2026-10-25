import math

import pytest
from shapely.geometry import Polygon

from flatmat.geom2d import cross, cut_slits, heading, perp, rigid_fit, transform, unit
from flatmat.model import Line, Loop2D


def test_vectors():
    u, n = unit((1.0, 1.0), (4.0, 5.0))
    assert n == 5.0 and u == pytest.approx((0.6, 0.8))
    assert unit((2.0, 3.0), (2.0, 3.0)) == ((1.0, 0.0), 0.0)  # coincident points: no division by zero
    assert perp((1.0, 0.0)) == (-0.0, 1.0) and cross((1.0, 0.0), perp((1.0, 0.0))) > 0
    assert heading((0.0, 0.0), (0.0, 2.0)) == pytest.approx(math.pi / 2)


def test_rigid_fit_maps_both_points_without_mirroring():
    b1, b2 = (1.0, 2.0), (4.0, 6.0)
    a1, a2 = (10.0, 0.0), (10.0, 5.0)  # same distance, turned and moved
    T = rigid_fit(a1, a2, b1, b2)
    assert transform(T, b1) == pytest.approx(a1) and transform(T, b2) == pytest.approx(a2)
    assert T[0] ** 2 + T[1] ** 2 == pytest.approx(1.0)


def test_a_slit_separates_what_it_cuts():
    square = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    slit = Loop2D((Line((5.0, -1.0), (5.0, 11.0)),))
    assert len(cut_slits(square, [slit], 0.01).geoms) == 2
