"""Countersinks and counterbores: cut the through hole, report the feature."""

import json
import math

import ezdxf
import pytest
from conftest import write_and_load
from parts import feature_cases

import flatmat as fm
from flatmat.model import Arc

K = 0.44
CASES = {p.name: p for p in feature_cases()}


@pytest.fixture(scope="module")
def feats(step_dir):
    return {name: (p, fm.recognize(write_and_load(p, step_dir))) for name, p in CASES.items()}


@pytest.mark.parametrize("name", list(CASES))
def test_through_holes_are_cut(feats, name):
    part, sheet = feats[name]
    pat = fm.unfold(sheet, k=K)
    assert len(pat.holes) == len(part.features)
    for h in pat.holes:
        (c,) = h.segments
        assert isinstance(c, Arc) and c.is_circle
        assert 2 * c.radius == pytest.approx(part.features[0][4])
    through = sum(math.pi * f[4] ** 2 / 4 for f in part.features)
    assert pat.area == pytest.approx(part.width * part.flat_length(K) - through, rel=1e-9)


def test_countersinks(feats):
    part, sheet = feats["csk_L"]
    assert sorted(f.side for f in sheet.features) == ["bottom", "top"]
    for f in sheet.features:
        assert f.kind == "countersink"
        assert f.through_diameter == pytest.approx(6.0)
        assert f.outer_diameter == pytest.approx(10.0, abs=1e-6)
        assert f.angle_deg == pytest.approx(90.0)
    pat = fm.unfold(sheet, k=K)
    assert any(d.code == fm.Code.COUNTERSINK_PROJECTED for d in pat.diagnostics)
    centers = sorted(f.center for f in pat.features)
    holes = sorted(h.segments[0].center for h in pat.holes)
    assert centers == pytest.approx(holes)


def test_counterbores(feats):
    _, sheet = feats["cbore_L"]
    assert {f.kind for f in sheet.features} == {"counterbore"}
    for f in sheet.features:
        assert f.through_diameter == pytest.approx(5.0)
        assert f.outer_diameter == pytest.approx(9.0)
        assert f.depth == pytest.approx(2.0, abs=1e-6)


def test_features_in_exports(feats, tmp_path):
    pat = fm.unfold(feats["csk_L"][1], k=K)
    fm.write_bend_table(pat, tmp_path / "t.json")
    t = json.loads((tmp_path / "t.json").read_text())
    assert len(t["features"]) == 2 and t["features"][0]["kind"] == "countersink"
    fm.write_dxf(pat, tmp_path / "f.dxf")
    doc = ezdxf.readfile(tmp_path / "f.dxf")
    texts = [e.dxf.text for e in doc.modelspace().query("TEXT")]
    assert sum(t.startswith("CSK 90°") for t in texts) == 2
