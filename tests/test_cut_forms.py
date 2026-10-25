"""Louvres, bridges and lances: formed out of the sheet and cut free along slits."""

import json
import math
from dataclasses import replace

import ezdxf
import pytest
from conftest import write_and_load
from parts import cut_formed_cases
from shapely.geometry import Point, box

import flatmat as fm
from flatmat.fold import FoldModel
from flatmat.model import Arc, Line, Loop2D

K = 0.44
CASES = {p.name: p for p in cut_formed_cases()}


@pytest.fixture(scope="module")
def flats(step_dir):
    out = {}
    for name, p in CASES.items():
        sheet = fm.recognize(write_and_load(p, step_dir))
        out[name] = (p, sheet, fm.unfold(sheet, k=K))
    return out


def _cut_area(spec) -> float:
    kind, *_, gap, _side = spec
    if kind in ("louvre_arc", "lance_formed", "bridge_formed"):
        rim = spec[4]
        disc = Point(0, 0).buffer(rim, quad_segs=4096)
        if kind == "louvre_arc":
            cx, rc = spec[8], spec[9]
            ring = Point(-cx, 0).buffer(rc, quad_segs=4096).difference(Point(-cx, 0).buffer(rc - gap, quad_segs=4096))
            return disc.intersection(ring).area
        width, tip = spec[8], spec[9]
        opening = disc.intersection(box(-rim - 1, -width / 2 - gap, tip + gap, width / 2 + gap))
        return opening.difference(disc.intersection(box(-rim - 1, -width / 2, tip, width / 2))).area
    if kind == "louvre":
        return 2 * spec[4] * gap  # the slit along the diameter
    if kind == "bridge":
        return 2 * spec[4] * gap  # two slits along the span
    length, width = spec[4], spec[5]
    return (length + gap) * (width + 2 * gap) - length * width  # the opening around the flat tab


def _bend_angle(pat, bend_id) -> float:
    bl = next(b for b in pat.bend_lines if b.bend_id == bend_id)
    return math.degrees(math.atan2(bl.p1[1] - bl.p0[1], bl.p1[0] - bl.p0[0])) % 180


@pytest.mark.parametrize("name", list(CASES))
def test_recognised_and_area_exact(flats, name):
    part, sheet, pat = flats[name]
    (spec,) = part.formed
    (f,) = pat.features
    assert f.kind == spec[0].split("_")[0]
    assert sheet.thickness == pytest.approx(part.t, abs=1e-6)
    exact = 1e-9 if spec[0] in ("louvre", "bridge", "lance") else 1e-5  # curved: through STEP, against polygons
    assert pat.area == pytest.approx(part.width * part.flat_length(K) - _cut_area(spec), rel=exact)
    assert pat.meta["area_rel_err"] < 1e-9
    assert f.holes and all(i < len(pat.holes) for i in f.holes)
    assert sum(abs(pat.holes[i].area()) for i in f.holes) == pytest.approx(_cut_area(spec), rel=exact)
    assert len(pat.holes) == len(f.holes)  # nothing else is cut
    assert fm.check(pat).ok  # slits and openings are not judged as holes or slots


def test_louvre(flats):
    part, sheet, pat = flats["louvre_up"]
    (f,) = pat.features
    assert f.height == pytest.approx(5.0) and f.bend_id is None
    assert {fr.role for fr in sheet.faces if fr.id in f.face_ids} == {"formed"}
    # the slit runs along the flat (across the L's bend line); the rim is a half circle and the slit
    assert abs(f.orientation_deg - _bend_angle(pat, pat.bends[0].id)) == pytest.approx(90.0)
    assert f.size == pytest.approx((20.0, 10.5))
    arcs = [s for s in f.outline.segments if isinstance(s, Arc)]
    assert arcs and all(a.radius == pytest.approx(10.0) for a in arcs)
    (slot,) = (pat.holes[i] for i in f.holes)
    assert all(isinstance(s, Line) for s in slot.segments)


def test_bridge(flats):
    part, sheet, pat = flats["bridge"]
    (f,) = pat.features
    assert f.height == pytest.approx(4.0) and len(f.holes) == 2
    assert abs(f.orientation_deg - _bend_angle(pat, pat.bends[0].id)) == pytest.approx(90.0)
    assert f.size == pytest.approx((30.0, 11.0))
    assert [b.angle_deg for b in pat.bends] == pytest.approx([90.0])  # the joggles are not press-brake bends


@pytest.mark.parametrize("name", ["lance", "lance_90"])
def test_lance_is_developed(flats, name):
    part, sheet, pat = flats[name]
    (spec,) = part.formed
    (f,) = pat.features
    length, width, angle, gap = spec[4], spec[5], spec[6], spec[8]
    tab = next(b for b in pat.bends if b.id == f.bend_id)
    assert tab.angle_deg == pytest.approx(angle) and tab.r_in == pytest.approx(2.0)
    assert f.orientation_deg == pytest.approx(_bend_angle(pat, tab.id))
    assert f.size == pytest.approx((width + 2 * gap, length + gap))
    a = math.radians(angle)
    straight = length - a * (2 + K * 2)  # the tab's flat part
    radius = 2.0 if tab.direction == "up" else 4.0  # the reference skin inside the bend, or outside
    assert f.height == pytest.approx(straight * math.sin(a) + radius * (1 - math.cos(a)), abs=1e-6)
    # the tab's flange is measured to the opening, not across it
    (bc,) = (c for c in fm.check(pat).bends if c.bend_id == tab.id)
    assert any(x == pytest.approx(straight + tab.ossb, abs=1e-6) for x in bc.flanges)
    plan = fm.plan_bends(pat)  # the tab is formed with the blank, not on the press brake
    assert plan.feasible and tab.id in plan.unchecked and tab.id not in plan.order


def test_curved_louvre(flats):
    _, sheet, pat = flats["louvre_arc"]
    (f,) = pat.features
    assert f.kind == "louvre" and f.height == pytest.approx(5.0)
    assert sheet._ctx.features[f.id].gap == pytest.approx(0.5, abs=1e-6)
    (slot,) = (pat.holes[i] for i in f.holes)
    arcs = [s for s in slot.segments if isinstance(s, Arc)]
    # the sheet's edge (r 24.5), the hood's edge projected (r 25), and the rim between them: exact arcs
    assert len(arcs) == len(slot.segments)
    assert sorted({round(a.radius, 6) for a in arcs}) == pytest.approx([12.0, 24.5, 25.0])


@pytest.mark.parametrize("name", ["lance_formed", "bridge_formed"])
def test_formed_lance_and_bridge(flats, name):
    _, sheet, pat = flats[name]
    (f,) = pat.features
    assert f.kind == name.split("_")[0] and f.bend_id is None and f.outline is not None
    assert {fr.role for fr in sheet.faces if fr.id in f.face_ids} == {"formed"}  # not developed
    assert len(f.holes) == (1 if f.kind == "lance" else 2)
    assert [b.angle_deg for b in pat.bends] == pytest.approx([90.0])
    # along the chord of the slits: across the tongue (width 10 + two slits), or the pad's length
    assert f.size[0] == pytest.approx(11.0 if f.kind == "lance" else 24.0, abs=0.05)


def test_sides(flats):
    for name in ("louvre_up", "louvre_down"):
        part, sheet, pat = flats[name]
        left_is_a = sheet.bends[0].direction == "up"
        (spec,) = part.formed
        assert pat.features[0].side == ("top" if (spec[-1] > 0) == left_is_a else "bottom")


def test_exports(flats, tmp_path):
    for name in ("louvre_up", "bridge", "lance"):
        _, _, pat = flats[name]
        fm.write_dxf(pat, tmp_path / f"{name}.dxf")
        doc = ezdxf.readfile(tmp_path / f"{name}.dxf")
        assert not doc.audit().has_errors
        texts = [e.dxf.text for e in doc.modelspace().query("TEXT")]
        assert any(t.startswith(name.split("_")[0].upper()) for t in texts)
        form = doc.modelspace().query("*[layer=='FORM']")
        assert bool(form) == (name != "lance")  # a bent lance is developed: no rim to mark
        fm.write_bend_table(pat, tmp_path / f"{name}.json")
        (f,) = json.loads((tmp_path / f"{name}.json").read_text())["features"]
        assert f["kind"] == name.split("_")[0] and f["orientation_deg"] is not None and f["holes"]


def test_slits_without_width(flats, tmp_path):
    """A slit with no width (a zero-area loop) is cut as a line and still separates the sheet."""
    _, _, pat = flats["lance"]
    (f,) = pat.features
    opening = pat.holes[f.holes[0]]
    pts = opening.polyline(0.01)
    # the far side of the opening as a cut line instead of the opening
    x0, y0, x1, y1 = (min(p[0] for p in pts), min(p[1] for p in pts), max(p[0] for p in pts), max(p[1] for p in pts))
    slit = Loop2D((Line((x0, y0), (x1, y0)), Line((x1, y0), (x1, y1)), Line((x1, y1), (x0, y1))))
    pat2 = replace(pat, holes=(), slits=(slit,), features=(replace(f, holes=()),))
    fm.write_dxf(pat2, tmp_path / "s.dxf")
    lines = ezdxf.readfile(tmp_path / "s.dxf").modelspace().query("LINE[layer=='INNER']")
    assert len(lines) == 3
    assert len(FoldModel(pat2).regions) == len(FoldModel(pat).regions)
    fm.check(pat2)
