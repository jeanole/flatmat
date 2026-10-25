"""Folding the flat back up must land the mid-surface t/2 inside the original solid."""

import math

import numpy as np
import pytest
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.gp import gp_Pnt
from OCP.TopAbs import TopAbs_IN, TopAbs_SHELL
from OCP.TopExp import TopExp_Explorer
from parts import box_solid, cases, make_solid, sharp_cases
from shapely.geometry import Point

import flatmat as fm
from flatmat.fold import FoldModel, apply

K = 0.44


def to_solid_frame(sheet, pat, q):
    """Points of the flat's frame (z towards the reference skin's outside) -> the STEP's frame."""
    f = sheet._ctx.faces[pat.meta["root_face"] - 1]
    pos = f.ad.Plane().Position()
    o, x, y = pos.Location(), pos.XDirection(), pos.YDirection()
    n = f.normal.Normalized()
    dx, dy = pat.meta["root_offset"]
    out = []
    for px, py, pz in q:
        u, v, w = (px - dx) / f.mirror, py - dy, pz - sheet.thickness / 2
        out.append(gp_Pnt(*(o.Coord(i) + u * x.Coord(i) + v * y.Coord(i) + w * n.Coord(i) for i in (1, 2, 3))))
    return out


def mid_surface_errors(sheet, pat, step=6.0):
    """Largest |distance to the skin - t/2| over points well inside each folded flange."""
    shell = TopExp_Explorer(sheet.occ_shape, TopAbs_SHELL).Current()
    m = FoldModel(pat)
    T = m.folded()
    t = sheet.thickness

    def errors(anchor):
        inv = np.linalg.inv(T[anchor])
        worst, n = 0.0, 0
        for r, M in T.items():
            inner = m.regions[r].buffer(-1.5 * t)
            if inner.is_empty:
                continue
            x0, y0, x1, y1 = inner.bounds
            g = [
                (x, y) for x in np.arange(x0, x1, step) for y in np.arange(y0, y1, step) if inner.contains(Point(x, y))
            ]
            for p in to_solid_frame(sheet, pat, apply(inv @ M, np.array(g)) if g else []):
                n += 1
                if BRepClass3d_SolidClassifier(sheet.occ_shape, p, 1e-6).State() != TopAbs_IN:
                    return math.inf, n
                d = BRepExtrema_DistShapeShape(BRepBuilderAPI_MakeVertex(p).Vertex(), shell).Value()
                worst = max(worst, abs(d - t / 2))
        return worst, n

    # the region holding the root face is the frame of the flat; try each as that anchor
    return min((errors(a) for a in T), key=lambda e: e[0])


PARTS = {p.name: p for p in cases() + sharp_cases() if p.name != "joggle"}


@pytest.mark.parametrize("name", list(PARTS))
def test_folding_restores_the_part(name):
    sheet = fm.recognize(make_solid(PARTS[name]))
    pat = fm.unfold(sheet, k=K)
    worst, n = mid_surface_errors(sheet, pat)
    assert n > 10 and worst < 1e-6


def test_folding_a_box():
    sheet = fm.recognize(box_solid(80, 60, 30, 2.0, 2.0, 2.0, corner_gap=0.1))
    pat = fm.unfold(sheet, k=K)
    m = FoldModel(pat)
    assert len(m.regions) == 5 and len(m.hinges) == 4 and not m.unseparated
    worst, n = mid_surface_errors(sheet, pat, step=4.0)
    assert n > 100 and worst < 1e-6


def test_the_mid_surface_does_not_depend_on_k():
    sheet = fm.recognize(make_solid(PARTS["U"]))
    for k in (0.3, 0.5):
        m = FoldModel(fm.unfold(sheet, k=k))
        T = m.folded()
        pts = np.vstack([apply(T[r], np.asarray(m.regions[r].exterior.coords)[:, :2]) for r in T])
        ext = sorted(pts.max(axis=0) - pts.min(axis=0))
        # U: flats 30/60/30, r 1.5, t 1.5, 80 wide: legs 60 + 2 (r + t/2) apart, 30 + r + t/2 high
        assert ext == pytest.approx(sorted([64.5, 80.0, 32.25]))


def test_joggle_without_a_flat_between_is_not_separated():
    (joggle,) = (p for p in cases() if p.name == "joggle")
    m = FoldModel(fm.unfold(fm.recognize(make_solid(joggle)), k=K))
    assert not m.hinges  # two tangent bends, no flange between: no hinge, reported as unchecked
