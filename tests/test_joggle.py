"""Joggles: two bends that offset the sheet, recognised as a pair, laid flat by a joggle table."""

import csv
import json
import math

import ezdxf
import pytest
from handworked_parts import cases as hand_cases
from parts import Part, cases, make_solid, sharp_cases

import flatmat as fm

K = 0.44
HAND = {c.name: c for c in hand_cases()}
PARTS = {p.name: p for p in cases() + sharp_cases()}


def ba(deg, r, t, k=K):
    return math.radians(deg) * (r + k * t)


def joggle_values(deg, r1, r2, t, flat):
    """Offset and run of a joggle, by hand from its centre line (radii r + t/2)."""
    th, rc = math.radians(deg), r1 + r2 + t
    return rc * (1 - math.cos(th)) + flat * math.sin(th), rc * math.sin(th) + flat * math.cos(th)


@pytest.fixture(scope="module")
def sheets():
    out = {n: fm.recognize(HAND[n].build()) for n in ("joggle_offset", "joggle_z_tight", "joggle_on_flange")}
    out["joggle"] = fm.recognize(make_solid(PARTS["joggle"]))
    return out


# (name, the joggle's bends' angle, r, t, flat between, ids of the paired bends)
JOGGLES = [
    ("joggle_offset", 45.0, 2.0, 2.0, 5.0, (1, 2)),
    ("joggle_z_tight", 90.0, 1.0, 1.5, 2.0, (1, 2)),
    ("joggle_on_flange", 20.0, 4.0, 2.0, 0.0, (2, 3)),
    ("joggle", 30.0, 3.0, 2.0, 0.0, (1, 2)),
]


@pytest.mark.parametrize(("name", "deg", "r", "t", "flat", "ids"), JOGGLES)
def test_recognised(sheets, name, deg, r, t, flat, ids):
    sheet = sheets[name]
    (j,) = sheet.joggles
    offset, run = joggle_values(deg, r, r, t, flat)
    assert j.bends == ids
    assert j.angle_deg == pytest.approx(deg)
    assert (j.offset, j.flat, j.run) == pytest.approx((offset, flat, run), abs=1e-9)
    kinds = {b.id: b.kind for b in sheet.bends}
    assert {i for i, kind in kinds.items() if kind == "joggle"} == set(ids)
    assert {b.direction for b in sheet.bends if b.id in ids} == {"up", "down"}


@pytest.mark.parametrize("name", ["Z", "hat", "U", "L", "sharp_in_Z", "sharp_both_Z"])
def test_not_joggles(name):
    # a Z or hat with a long web is two bends, not an offset; same-way bends, sharp bends: never
    sheet = fm.recognize(make_solid(PARTS[name]))
    assert sheet.joggles == ()
    assert all(b.kind != "joggle" for b in sheet.bends)


def test_the_offset_limit():
    part = HAND["joggle_offset"].build()  # offset 5.29, t 2
    assert fm.recognize(part, joggle_max_offset=0).joggles == ()
    assert fm.recognize(part, joggle_max_offset=5.2).joggles == ()
    assert len(fm.recognize(part, joggle_max_offset=5.3).joggles) == 1
    z = Part("Z_long", t=2.0, width=50, flats=[30, 20, 30], bends=[(60, 2.0), (-60, 2.0)])
    assert fm.recognize(make_solid(z)).joggles == ()  # offset 6 (1 - cos 60) + 20 sin 60 = 20.3 > 5 t
    assert fm.recognize(make_solid(z), joggle_max_offset=20).joggles == ()
    assert len(fm.recognize(make_solid(z), joggle_max_offset=21).joggles) == 1


def test_unequal_angles_are_not_a_joggle():
    p = Part("skew", t=2.0, width=50, flats=[30, 3, 30], bends=[(45, 2.0), (-40, 2.0)])
    assert fm.recognize(make_solid(p)).joggles == ()


@pytest.mark.parametrize(("name", "deg", "r", "t", "flat", "ids"), JOGGLES)
def test_length_and_allowance_from_k(sheets, name, deg, r, t, flat, ids):
    pat = fm.unfold(sheets[name], k=K)
    (j,) = pat.joggles
    _, run = joggle_values(deg, r, r, t, flat)
    assert j.source == "k"
    assert j.length == pytest.approx(2 * ba(deg, r, t) + flat)
    assert j.allowance == pytest.approx(j.length - run)
    assert all(b.k == K for b in pat.bends)


TABLE = fm.JoggleTable.from_rows([(0.0, 0.1), (4.0, 0.5), (8.0, 1.3)])


def table_value(offset):
    return 0.1 + 0.1 * offset if offset <= 4 else 0.5 + 0.2 * (offset - 4)


@pytest.mark.parametrize(("name", "deg", "r", "t", "flat", "ids"), JOGGLES)
def test_a_joggle_table_sets_the_flat_length(sheets, name, deg, r, t, flat, ids):
    sheet = sheets[name]
    offset, run = joggle_values(deg, r, r, t, flat)
    plain, tabled = fm.unfold(sheet, k=K), fm.unfold(sheet, k=K, joggle=TABLE)
    (j,) = tabled.joggles
    assert j.source == "table"
    assert j.allowance == pytest.approx(table_value(offset))
    assert j.length == pytest.approx(run + table_value(offset))
    # the flat grows by what the joggle gains, and nothing else changes
    gain = j.length - plain.joggles[0].length
    assert sorted(tabled.size)[1] == pytest.approx(sorted(plain.size)[1] + gain, abs=1e-9)
    assert tabled.meta["area_rel_err"] < 1e-9
    by_id = {b.id: b for b in tabled.bends}
    assert by_id[ids[0]].ba + flat + by_id[ids[1]].ba == pytest.approx(j.length)
    assert by_id[ids[0]].ba == pytest.approx(by_id[ids[1]].ba)  # equal radii: an equal share
    for b in tabled.bends:
        if b.id not in ids:
            assert b.k == K and b.ba == pytest.approx(ba(b.angle_deg, b.r_in, t))


def test_holes_either_side_move_with_the_table(sheets):
    """joggle_offset: d6 holes at the middle of both 40 mm flats, at z 20 and 40."""
    for model in (None, TABLE):
        pat = fm.unfold(sheets["joggle_offset"], k=K, joggle=model)
        (j,) = pat.joggles
        c1, c2 = (h.segments[0].center for h in pat.holes)
        assert math.dist(c1, c2) == pytest.approx(math.hypot(20 + j.length + 20, 20), abs=1e-9)


def test_a_callable_and_unequal_radii():
    p = Part("jog_r", t=2.0, width=40, flats=[30, 4, 30], bends=[(30, 2.0), (-30, 6.0)])
    sheet = fm.recognize(make_solid(p))
    (j,) = sheet.joggles
    offset, run = joggle_values(30, 2.0, 6.0, 2.0, 4.0)
    assert (j.offset, j.run) == pytest.approx((offset, run))
    pat = fm.unfold(sheet, k=K, joggle=lambda jog, t: 0.25 * t)
    (j,) = pat.joggles
    assert j.allowance == pytest.approx(0.5)
    b1, b2 = pat.bends
    # shared in proportion to the K model's allowances
    assert b1.ba / b2.ba == pytest.approx(ba(30, 2.0, 2.0) / ba(30, 6.0, 2.0))
    assert b1.ba + 4.0 + b2.ba == pytest.approx(run + 0.5)


def test_table_warnings(sheets):
    sheet = sheets["joggle_offset"]  # offset 5.29
    short = fm.JoggleTable.from_rows([(1.0, 0.2), (3.0, 0.4)], thickness=1.5, name="t1.5")
    pat = fm.unfold(sheet, k=K, joggle=short)
    codes = {d.code for d in pat.diagnostics}
    assert fm.Code.JOGGLE_EXTRAPOLATED in codes
    assert pat.joggles[0].allowance == pytest.approx(0.4)  # clamped
    tabled = fm.unfold(sheet, k=K, joggle=fm.JoggleTable.from_rows([(0, 0.3), (9, 0.3)], thickness=1.5))
    assert fm.Code.JOGGLE_TABLE_THICKNESS in {d.code for d in tabled.diagnostics}
    with pytest.raises(fm.FlatmatError) as e:
        fm.unfold(sheet, k=K, joggle=fm.JoggleTable.from_rows([(0, -5.0)]))
    assert e.value.diagnostic.code == fm.Code.JOGGLE_ALLOWANCE
    with pytest.raises(ValueError):
        fm.JoggleTable((1.0, 0.5), (0.1, 0.2))


def test_load_csv(tmp_path):
    p = tmp_path / "jog.csv"
    p.write_text("# 2 mm steel\noffset,allowance\n0,0.1\n4,0.5\n8,1.3\n")
    t = fm.JoggleTable.load(p, thickness=2.0)
    assert (t.offset, t.allowance, t.thickness) == ((0.0, 4.0, 8.0), (0.1, 0.5, 1.3), 2.0)


def test_exports(sheets, tmp_path):
    pat = fm.unfold(sheets["joggle_on_flange"], k=K, joggle=TABLE)
    table = fm.bend_table(pat)
    (j,) = table["joggles"]
    assert j["bends"] == [2, 3] and j["source"] == "table"
    assert j["allowance"] == pytest.approx(table_value(j["offset"]))
    assert [(b["id"], b["kind"], b["joggle"]) for b in table["bends"]] == [(1, "bend", None), (2, "joggle", 1),
                                                                            (3, "joggle", 1)]  # fmt: skip
    json.dumps(table)
    fm.write_bend_table(pat, tmp_path / "t.csv")
    rows = list(csv.DictReader((tmp_path / "t.csv").open()))
    assert [(r["kind"], r["joggle"]) for r in rows] == [("bend", ""), ("joggle", "1"), ("joggle", "1")]
    fm.write_dxf(pat, tmp_path / "f.dxf")
    texts = {e.dxf.text for e in ezdxf.readfile(tmp_path / "f.dxf").modelspace().query("TEXT")}
    assert {"J+20 R4", "J-20 R4"} <= texts and any(t.endswith("90 R2") and t[0] != "J" for t in texts)


def test_dfm_and_bend_order(sheets):
    pat = fm.unfold(sheets["joggle_z_tight"], k=K)  # a 2 mm flat between the bends: no press brake bends it
    report = fm.check(pat)
    (info,) = (d for d in report.findings if d.code == fm.Code.DFM_JOGGLE)
    assert info.severity == fm.Severity.INFO and info.data["offset"] == pytest.approx(5.5)
    assert fm.Code.DFM_FLANGE_TOO_SHORT not in report.codes()
    plan = fm.plan_bends(fm.unfold(sheets["joggle_on_flange"], k=K))
    assert set(plan.unchecked) >= {2, 3} and [s.bend_id for s in plan.steps] == [1]


def test_cli(tmp_path, capsys):
    from parts import write_step

    from flatmat.cli import main

    step, tab, js = tmp_path / "j.step", tmp_path / "jog.csv", tmp_path / "j.json"
    write_step(HAND["joggle_offset"].build(), str(step))
    tab.write_text("offset,allowance\n0,0.1\n4,0.5\n8,1.3\n")
    assert main(["unfold", str(step), "--joggle-table", str(tab), "--table", str(js)]) == 0
    out = capsys.readouterr().out
    assert "joggle 1: bends 1+2" in out and "(table)" in out
    assert json.loads(js.read_text())["joggles"][0]["allowance"] == pytest.approx(table_value(5.292893218813452))
    assert main(["unfold", str(step), "--joggle-max-offset", "0", "--table", str(js)]) == 0
    assert json.loads(js.read_text())["joggles"] == []
    assert main(["unfold", str(step), "--joggle-table", str(tmp_path / "missing.csv")]) == 2
