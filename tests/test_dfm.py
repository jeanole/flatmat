"""Manufacturability checks: each test part carries one deliberate problem."""

import json

import ezdxf
import pytest
from parts import Part, cases, dfm_cases, make_solid

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


def _perforated():
    """Round holes far apart (left out of the web and slot morphology), two too close to be,
    and a narrow slot (not round: kept in, and reported)."""
    holes = [("flat", 0, f, z, 4.0) for f in (0.2, 0.5) for z in (20.0, 30.0, 40.0, 50.0)]
    holes += [("flat", 0, 0.8, 75.0, 4.0), ("flat", 0, 0.8, 80.5, 4.0)]  # 1.5 mm apart
    cuts = [((10.0, -2.0, 88.0), (30.0, 2.0, 89.0))]
    return Part("dfm_perforated", t=2.0, width=100, flats=[50, 40], bends=[(90, 2.0)], holes=holes, cuts=cuts)


def test_leaving_isolated_holes_out_changes_no_finding(flats, monkeypatch):
    from flatmat import dfm

    parts = dict(flats) | {"dfm_perforated": fm.unfold(fm.recognize(make_solid(_perforated())), k=K)}
    seen = {}
    real = dfm._isolated_holes

    def counting(sheet, rad):
        seen[id(sheet)] = real(sheet, rad)
        return seen[id(sheet)]

    monkeypatch.setattr(dfm, "_isolated_holes", counting)
    with_shortcut = {name: fm.check(p).to_dict() for name, p in parts.items()}
    assert len(next(v for v in seen.values() if len(v) >= 8)) == 8  # the 8 round holes far apart
    monkeypatch.setattr(dfm, "_isolated_holes", lambda sheet, rad: set())
    for name, p in parts.items():
        assert fm.check(p).to_dict() == with_shortcut[name], name
    assert "DFM_FEATURE_TOO_SMALL" in {f["code"] for f in with_shortcut["dfm_perforated"]["findings"]}


@pytest.mark.parametrize("name", list(EXPECT))
def test_disabling_a_rule_removes_only_its_findings(flats, name):
    # a switched-off rule skips its work; that must not change what the other rules find
    # (e.g. a hole near the edge is still not reported again as a narrow web)
    def found(report):
        return [(d.code.value, d.message) for d in report.findings]

    every = found(fm.check(flats[name]))
    for code in sorted(c.value for c in fm.Code if c.value.startswith("DFM_")):
        report = fm.check(flats[name], rules=fm.DfmRules(disabled=frozenset({code})))
        assert found(report) == [f for f in every if f[0] != code], code


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
