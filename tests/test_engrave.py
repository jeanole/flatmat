"""Engravings: recesses that do not go through the sheet, marked along their rims.

All on the L of the formed tests (t = 2, r = 2, 100 wide; base 50, flange 40). Profile
coordinates: the base runs along +x at y in [-1, 1], the bend's tangent line is at x = 50
and the flange rises from y = 3 at x in [52, 54] (inner skin at 52). The width runs along z.
"""

from __future__ import annotations

import math
import os

import ezdxf
import pytest
from conftest import write_and_load
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakePrism
from OCP.gp import gp_Ax2, gp_Dir, gp_Elips, gp_Pnt, gp_Vec
from OCP.ShapeUpgrade import ShapeUpgrade_UnifySameDomain
from parts import Part, make_solid, write_step

import flatmat as fm
from flatmat.model import Arc

K = 0.4
L = dict(t=2.0, width=100, flats=[50, 40], bends=[(90, 2.0)])
LENGTH = Part("L", **L).flat_length(K)  # 94.398...: base 50, bend 4.398, flange 40


def _load(solid, name, directory):
    path = os.path.join(directory, f"{name}.step")
    write_step(solid, path)
    return fm.load_step(path)


def _cyl(x, y0, y1, z, r):
    """A cylinder along y from y0 to y1."""
    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x, y0, z), gp_Dir(0, 1 if y1 > y0 else -1, 0)), r, abs(y1 - y0))


def _cyl_x(x0, x1, y, z, r):
    """A cylinder along x from x0 to x1 (into the flange)."""
    return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x0, y, z), gp_Dir(1 if x1 > x0 else -1, 0, 0)), r, abs(x1 - x0))


def _unify(shape):
    u = ShapeUpgrade_UnifySameDomain(shape, True, True, False)
    u.Build()
    return u.Shape()


def _ring(solid, cyl):
    """A ring groove: a round pocket cut by ``cyl(r=8)`` with an island of radius 4 left standing."""
    groove, island = cyl(8.0)
    return _unify(BRepAlgoAPI_Fuse(BRepAlgoAPI_Cut(solid, groove.Shape()).Shape(), island.Shape()).Shape())


def _ellipse_pocket(solid):
    """An elliptical pocket 0.4 deep into the base's +y skin, semi-axes 12 (along x) and 6."""
    e = gp_Elips(gp_Ax2(gp_Pnt(25, 0.6, 50), gp_Dir(0, 1, 0), gp_Dir(1, 0, 0)), 12.0, 6.0)
    face = BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(e).Edge()).Wire()).Face()
    return BRepAlgoAPI_Cut(solid, BRepPrimAPI_MakePrism(face, gp_Vec(0, 2, 0)).Shape()).Shape()


def _cases():
    plain = lambda: make_solid(Part("L", **L))  # noqa: E731
    return {
        "pocket_top": make_solid(Part("pocket_top", cuts=[((10, 0.5, 20), (30, 2, 40))], **L)),
        "pocket_bottom": make_solid(Part("pocket_bottom", cuts=[((10, -2, 20), (30, -0.5, 40))], **L)),
        "pocket_flange_out": make_solid(Part("pocket_flange_out", cuts=[((53.5, 13, 20), (55, 23, 40))], **L)),
        "pocket_flange_in": make_solid(Part("pocket_flange_in", cuts=[((51, 13, 20), (52.5, 23, 40))], **L)),
        "blind_hole": BRepAlgoAPI_Cut(plain(), _cyl(20, 2, 0.4, 70, 3).Shape()).Shape(),
        "ring_base": _ring(plain(), lambda r: (_cyl(20, 2, 0.7, 50, r), _cyl(20, 0, 1, 50, 4))),
        "ring_flange_in": _ring(plain(), lambda r: (_cyl_x(51, 52.3, 23, 50, r), _cyl_x(52.5, 52.0, 23, 50, 4))),
        "ellipse": _ellipse_pocket(plain()),
    }


@pytest.fixture(scope="module")
def engraved(step_dir):
    """name -> (model, flat pattern unfolded with engrave=True)."""
    out = {}
    for name, solid in _cases().items():
        model = _load(solid, name, step_dir)
        out[name] = (model, fm.unfold(fm.recognize(model, engrave=True), k=K))
    return out


def _bbox(loop):
    pts = loop.polyline(0.001)
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    return min(xs), min(ys), max(xs), max(ys)


@pytest.mark.parametrize("name", list(_cases()))
def test_refused_without_engrave(engraved, name):
    model, _ = engraved[name]
    with pytest.raises(fm.FlatmatError) as e:
        fm.recognize(model)
    assert e.value.diagnostic.code is fm.Code.EDGE_FEATURE


@pytest.mark.parametrize("name", list(_cases()))
def test_flat_is_the_plain_part(engraved, name):
    """The recess changes nothing about the flat: same size, same area, nothing cut."""
    _, pat = engraved[name]
    assert pat.size == pytest.approx((LENGTH, 100.0), abs=1e-6)
    assert pat.area == pytest.approx(LENGTH * 100.0, rel=1e-9)
    assert pat.meta["area_rel_err"] < 1e-9
    assert not pat.holes and not pat.slits
    assert all(f.kind == "engraving" for f in pat.features)
    (d,) = (d for d in pat.diagnostics if d.code is fm.Code.ENGRAVED_RECESS)
    assert d.severity is fm.Severity.WARNING and d.data["recesses"] == 1


@pytest.mark.parametrize(
    ("name", "depth", "along", "across", "size"),
    [
        # along: distance of the rim from the base's free end, in the flat; across: from z = 0
        ("pocket_top", 0.5, (10, 30), (20, 40), (20, 20)),
        ("pocket_bottom", 0.5, (10, 30), (20, 40), (20, 20)),
        ("pocket_flange_out", 0.5, (50 + 4.398 + 10, 50 + 4.398 + 20), (20, 40), (10, 20)),
        ("pocket_flange_in", 0.5, (50 + 4.398 + 10, 50 + 4.398 + 20), (20, 40), (10, 20)),
        ("blind_hole", 0.6, (17, 23), (67, 73), (6, 6)),
        ("ellipse", 0.4, (13, 37), (44, 56), (24, 12)),
    ],
)
def test_rim_in_the_flat(engraved, name, depth, along, across, size):
    _, pat = engraved[name]
    (f,) = pat.features
    assert f.depth == pytest.approx(depth, abs=1e-6)
    x0, y0, x1, y1 = _bbox(f.outline)
    # the flat is the L laid along x, the base's free end at x = 0 or x = LENGTH, z along y or mirrored
    assert (x1 - x0, y1 - y0) == pytest.approx(size, abs=2e-3)
    assert min((x0, x1), (LENGTH - x1, LENGTH - x0), key=lambda a: abs(a[0] - along[0])) == pytest.approx(
        along, abs=2e-3
    )
    assert min((y0, y1), (100 - y1, 100 - y0), key=lambda a: abs(a[0] - across[0])) == pytest.approx(across, abs=2e-3)
    assert f.outline.area() > 0
    assert f.center == pytest.approx(((x0 + x1) / 2, (y0 + y1) / 2), abs=2e-3)


def test_sides(engraved):
    """``side`` is the skin the recess opens into: top on the reference skin of the flat."""
    for top, bottom in (("pocket_flange_in", "pocket_flange_out"),):
        assert engraved[top][1].features[0].side != engraved[bottom][1].features[0].side
    sheet = fm.recognize(engraved["pocket_flange_in"][0], engrave=True)
    (f,) = sheet.features
    root = next(fr for fr in sheet.faces if fr.id == sheet.reference_root)
    assert f.side == "top" and root.role == "skin_a"


def test_circular_rims(engraved):
    _, pat = engraved["blind_hole"]
    (f,) = pat.features
    (rim,) = f.outline.segments
    assert isinstance(rim, Arc) and rim.is_circle and rim.radius == pytest.approx(3.0)
    assert f.outer_diameter == pytest.approx(6.0)


@pytest.mark.parametrize("name", ["ring_base", "ring_flange_in"])
def test_ring_groove_is_two_rims(engraved, name):
    """The island a ring groove leaves is not engraved over: its rim is marked as well."""
    _, pat = engraved[name]
    outer, inner = sorted(pat.features, key=lambda f: -f.outer_diameter)
    assert (outer.outer_diameter, inner.outer_diameter) == pytest.approx((16.0, 8.0))
    assert outer.depth == inner.depth
    assert outer.center == pytest.approx(inner.center, abs=1e-6)
    for f, r in ((outer, 8.0), (inner, 4.0)):
        (rim,) = f.outline.segments
        assert isinstance(rim, Arc) and rim.is_circle and rim.radius == pytest.approx(r)


def test_ring_on_the_flange_sits_on_the_flange(engraved):
    _, pat = engraved["ring_flange_in"]
    centre = pat.features[0].center
    # 20 mm up the flange from its tangent line, in the middle of the width
    assert min(abs(centre[0] - (50 + 4.398 + 20)), abs(LENGTH - centre[0] - (50 + 4.398 + 20))) < 2e-3
    assert centre[1] == pytest.approx(50.0)


def test_dxf_layer(engraved, tmp_path):
    _, pat = engraved["pocket_top"]
    path = tmp_path / "engraved.dxf"
    fm.write_dxf(pat, path)
    msp = ezdxf.readfile(path).modelspace()
    assert len(msp.query('LINE[layer=="ENGRAVE"]')) == 4
    assert not msp.query('*[layer=="FORM"]')
    assert not any("ENGRAVE" in t.dxf.text for t in msp.query("TEXT"))  # no label per rim
    table = fm.export.bend_table(pat)
    assert [f["kind"] for f in table["features"]] == ["engraving"]


@pytest.mark.parametrize(
    "cuts",
    [
        # a rebate along the base's free end: it reaches the outline
        [((-1, 0.5, 20), (10, 2, 40))],
        # a groove across the bend line: it runs into the bend
        [((40, 0.5, 20), (51.5, 2, 30))],
    ],
)
def test_still_refused(step_dir, cuts):
    model = write_and_load(Part(f"refused_{cuts[0][0][0]}", cuts=cuts, **L), step_dir)
    with pytest.raises(fm.FlatmatError) as e:
        fm.recognize(model, engrave=True)
    assert e.value.diagnostic.code is fm.Code.EDGE_FEATURE


def test_through_cuts_untouched(step_dir):
    """A recess beside a through hole: the hole is cut, the recess marked."""
    part = Part("pocket_and_hole", cuts=[((10, 0.5, 20), (30, 2, 40))], holes=[("flat", 0, 0.4, 70.0, 8.0)], **L)
    pat = fm.unfold(fm.recognize(write_and_load(part, step_dir), engrave=True), k=K)
    (hole,) = pat.holes
    assert -hole.area() == pytest.approx(math.pi * 16, rel=1e-9)
    assert [f.kind for f in pat.features] == ["engraving"]
    assert pat.area == pytest.approx(LENGTH * 100 - math.pi * 16, rel=1e-9)
