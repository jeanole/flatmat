"""Manufacturability checks: each test part carries one deliberate problem."""

import json

import ezdxf
import pytest
from parts import cases, dfm_cases, make_solid

import flatmat as fm
from flatmat.dfm import bending_force_kn, die_opening

K = 0.44
DFM = {p.name: p for p in dfm_cases()}
EXPECT = {
    "dfm_ok": set(),
    "dfm_short_flange": {"DFM_FLANGE_TOO_SHORT"},
    "dfm_hole_near_bend": {"DFM_HOLE_NEAR_BEND"},
    "dfm_small_hole": {"DFM_HOLE_TOO_SMALL"},
    "dfm_hole_near_edge": {"DFM_HOLE_NEAR_EDGE"},
    "dfm_holes_close": {"DFM_HOLES_TOO_CLOSE"},
    "dfm_no_relief": {"DFM_RELIEF_MISSING"},
    "dfm_relief": set(),
    "dfm_narrow_slot": {"DFM_FEATURE_TOO_SMALL"},
    "dfm_short_hem": {"DFM_HEM_TOO_SHORT"},
}


@pytest.fixture(scope="module")
def flats():
    return {name: fm.unfold(fm.recognize(make_solid(p)), k=K) for name, p in DFM.items()}


@pytest.mark.parametrize("name", list(EXPECT))
def test_each_part_triggers_exactly_its_rule(flats, name):
    report = fm.check(flats[name])
    assert {d.code.value for d in report.problems} == EXPECT[name]
    assert report.ok == (not EXPECT[name])
    for d in report.problems:
        assert "at" in d.data  # every finding is located in the flat


def test_values_of_a_clean_L(flats):
    report = fm.check(flats["dfm_ok"])
    (b,) = report.bends
    assert report.die_v == 12.0  # 6t for t = 2
    assert b.flanges == pytest.approx((54.0, 44.0)) or b.flanges == pytest.approx((44.0, 54.0))  # legs + (r + t)
    assert b.min_flange == pytest.approx(8.4)  # max(0.7 V, 2.5 t + r)
    assert b.force_kn == pytest.approx(1.42 * 400 * 4 * 100 / (1000 * 12))
    assert b.air_radius == pytest.approx(0.16 * 12)


def test_measured_flange_and_required(flats):
    (d,) = fm.check(flats["dfm_short_flange"]).problems
    assert d.data["flange"] == pytest.approx(7.0)  # 3 mm flat + (r + t)
    assert d.data["required"] == pytest.approx(8.4)


def test_strict_makes_errors(flats):
    report = fm.check(flats["dfm_short_flange"], strict=True)
    assert all(d.severity == fm.Severity.ERROR for d in report.problems)


def test_rules_are_configurable(flats):
    assert fm.check(flats["dfm_short_flange"], rules=fm.DfmRules(min_flange_v=0.5, min_flange_t=2.0)).ok
    assert fm.check(flats["dfm_short_flange"], rules=fm.DfmRules(disabled=frozenset({"DFM_FLANGE_TOO_SHORT"}))).ok


def test_die_opening_rule():
    assert die_opening(1.0) == 6.0 and die_opening(2.9999) == pytest.approx(8 * 2.9999) and die_opening(10) == 100.0
    assert bending_force_kn(450, 2, 1000, 16) == pytest.approx(159.75)


def test_tonnage_and_material(flats):
    pat = flats["dfm_ok"]
    assert "DFM_TONNAGE" in fm.check(pat, machine_max_force_kn=10).codes()
    hem = flats["dfm_short_hem"]
    assert "DFM_HEM_MATERIAL" in fm.check(hem, material="aluminium_6061_t6").codes()
    with pytest.raises(ValueError, match="unknown material 'steel304'"):
        fm.check(pat, material="steel304")


def test_grain_direction(flats):
    pat = flats["dfm_ok"]  # r = 2 = 1.0 t: fine across the grain, too small along it for stainless (1.5 t)
    (bl,) = pat.bend_lines
    along = 0.0 if abs(bl.p0[1] - bl.p1[1]) < 1e-9 else 90.0
    assert "DFM_RADIUS_TOO_SMALL" not in fm.check(pat, material="stainless").codes()
    assert "DFM_RADIUS_TOO_SMALL" in fm.check(pat, material="stainless", rolling_direction_deg=along).codes()
    assert "DFM_RADIUS_TOO_SMALL" not in fm.check(pat, material="stainless", rolling_direction_deg=along + 90).codes()


def test_u_channel_and_hole_in_bend():
    by = {p.name: p for p in cases()}
    hat = fm.unfold(fm.recognize(make_solid(by["hat"])), k=K)
    assert "DFM_U_CHANNEL" in fm.check(hat).codes()
    hb = fm.unfold(fm.recognize(make_solid(by["L_hole_in_bend"])), k=K)
    codes = [d.code.value for d in fm.check(hb).problems]
    assert codes.count("DFM_CUT_IN_BEND_ZONE") == 1 and "DFM_RELIEF_MISSING" not in codes


def test_air_bending_radius_is_only_a_note(flats):
    report = fm.check(flats["dfm_ok"], die_v=40)  # r_air = 6.4 vs model 2
    notes = [d for d in report.findings if d.code == fm.Code.DFM_RADIUS_NOT_AIR_BENDABLE]
    assert notes and notes[0].severity == fm.Severity.INFO and report.ok


def test_reports_in_exports(flats, tmp_path):
    pat = flats["dfm_hole_near_bend"]
    report = fm.check(pat)
    fm.write_bend_table(pat, tmp_path / "t.json", dfm=report)
    t = json.loads((tmp_path / "t.json").read_text())
    assert t["dfm"]["findings"][0]["code"] == "DFM_HOLE_NEAR_BEND"
    fm.write_dxf(pat, tmp_path / "d.dxf", dfm=report)
    doc = ezdxf.readfile(tmp_path / "d.dxf")
    assert not doc.audit().has_errors
    assert len(doc.modelspace().query("CIRCLE[layer=='DFM']")) == 1
