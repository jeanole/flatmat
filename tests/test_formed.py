"""Embosses and extruded holes: recognised, reported, not cut (emboss) or cut as a pre-hole."""

import json
import math

import ezdxf
import pytest
from conftest import write_and_load
from parts import Part, formed_cases, make_solid

import flatmat as fm
from flatmat.bendmodel import extruded_hole_prehole
from flatmat.model import Arc

K = 0.44
CASES = {p.name: p for p in formed_cases()}


@pytest.fixture(scope="module")
def formed(step_dir):
    return {name: (p, fm.recognize(write_and_load(p, step_dir))) for name, p in CASES.items()}


def _sides(part: Part, sheet: fm.SheetPart) -> dict[float, str]:
    """Expected Feature.side by height: +1 in the part is the left skin; 'top' if that is the reference skin."""
    left_is_a = sheet.bends[0].direction == "up"  # the part's only bend turns left
    return {round(f[5], 6): "top" if (f[-1] > 0) == left_is_a else "bottom" for f in part.formed}


@pytest.mark.parametrize("name", list(CASES))
def test_recognised_and_area_kept(formed, name):
    part, sheet = formed[name]
    assert sheet.thickness == pytest.approx(part.t, abs=1e-6)
    assert [(round(b.angle_deg, 6), round(b.r_in, 6)) for b in sheet.bends] == [(90.0, 2.0)]
    assert len(sheet.features) == len(part.formed)
    sides = _sides(part, sheet)
    for f in sheet.features:
        assert f.side == sides[round(f.height, 6)]
    pat = fm.unfold(sheet, k=K)
    cut = sum(math.pi * f.through_diameter**2 / 4 for f in pat.features)
    # an emboss's rim is not cut and the formed region is laid flat inside its rim
    assert pat.area == pytest.approx(part.width * part.flat_length(K) - cut, rel=1e-9)
    assert pat.meta["area_rel_err"] < 1e-9
    assert all(any(fr.role == "formed" for fr in sheet.faces if fr.id == i) for f in sheet.features for i in f.face_ids)


def test_emboss(formed):
    part, sheet = formed["emboss_up"]
    (f,) = sheet.features
    assert f.kind == "emboss" and f.through_diameter == 0.0
    assert f.height == pytest.approx(5.0) and f.outer_diameter == pytest.approx(20.0)
    pat = fm.unfold(sheet, k=K)
    assert not pat.holes  # a marking, not a cut
    (f,) = pat.features
    (rim,) = f.outline.segments
    assert isinstance(rim, Arc) and rim.is_circle and rim.radius == pytest.approx(10.0)
    assert f.center == pytest.approx(rim.center)
    # in the flat: half way along the 50 mm base from its free end, 50 across the width
    along, across = sorted(f.center, key=lambda c: abs(c - 50.0))[::-1]
    assert across == pytest.approx(50.0)
    assert along == pytest.approx(25.0) or along == pytest.approx(part.flat_length(K) - 25.0)


def test_extruded_hole(formed):
    part, sheet = formed["collar_up"]
    (f,) = sheet.features
    assert f.kind == "extruded_hole"
    assert f.outer_diameter == pytest.approx(2 * (6.0 - 1.0 - 2.0))  # bore
    assert f.height == pytest.approx(4.5) and f.radius == pytest.approx(1.0)
    # the published rule d = D - 2 (H - 0.43 r - 0.72 t) is the neutral fibre with K = 0.5
    D = f.outer_diameter + part.t
    published = D - 2 * (f.height - (2 - math.pi / 2) * f.radius - (1.5 - math.pi / 4) * part.t)
    assert f.through_diameter == pytest.approx(published)
    assert f.through_diameter == pytest.approx(extruded_hole_prehole(6.0, 1.0, math.pi / 2, 4.5 - 2 - 1, 2.0))
    pat = fm.unfold(sheet, k=K)
    (h,) = pat.holes
    (c,) = h.segments
    assert isinstance(c, Arc) and c.is_circle and 2 * c.radius == pytest.approx(f.through_diameter)
    assert pat.features[0].center == pytest.approx(c.center)
    assert pat.features[0].outline.segments[0].radius == pytest.approx(6.0)


def test_pre_hole_override(formed):
    _, sheet = formed["collar_down_two"]
    pat = fm.unfold(sheet, k=K, pre_hole=lambda f: 0.5 * f.outer_diameter)
    assert sorted(2 * h.segments[0].radius for h in pat.holes) == pytest.approx([3.5, 5.0])
    pat = fm.unfold(sheet, k=K, pre_hole=3.0)
    assert [2 * h.segments[0].radius for h in pat.holes] == pytest.approx([3.0, 3.0])


def test_too_tall_collar_is_pierced(step_dir):
    part = Part("tall", t=2.0, width=100, flats=[50, 40], bends=[(90, 2.0)],
                formed=[("collar", 0, 0.5, 50.0, 8.0, 8.0, 1.0, +1)])  # fmt: skip
    sheet = fm.recognize(make_solid(part))
    assert sheet.features[0].through_diameter == 0.0
    assert any(d.code == fm.Code.PRE_HOLE for d in sheet.diagnostics)
    assert not fm.unfold(sheet, k=K).holes


def test_formed_in_exports(formed, tmp_path):
    pat = fm.unfold(formed["emboss_up"][1], k=K)
    fm.write_bend_table(pat, tmp_path / "t.json")
    (f,) = json.loads((tmp_path / "t.json").read_text())["features"]
    assert f["kind"] == "emboss" and f["height"] == pytest.approx(5.0)
    fm.write_dxf(pat, tmp_path / "e.dxf")
    doc = ezdxf.readfile(tmp_path / "e.dxf")
    assert not doc.audit().has_errors
    assert len(doc.modelspace().query("CIRCLE[layer=='FORM']")) == 1
    assert not doc.modelspace().query("CIRCLE[layer=='INNER']")
    assert any(e.dxf.text.startswith("EMBOSS") for e in doc.modelspace().query("TEXT"))
