"""Tool collisions and the bend order."""

import csv
import json

import ezdxf
import pytest
from parts import Part, box_solid, cases, make_solid

import flatmat as fm
from flatmat.tooling import DIES, PUNCHES, Tooling

K = 0.44
BY = {p.name: p for p in cases()}
EXTRA = {
    # a Z whose long first flange hangs into the lower beam if it is bent first
    "z_long_leg": Part("z_long_leg", t=2.0, width=100, flats=[100, 20, 40], bends=[(90, 2.0), (-90, 2.0)]),
    # a U 20 wide with 60 mm legs: the second leg closes onto a straight punch
    "u_deep": Part("u_deep", t=1.5, width=100, flats=[60, 20, 60], bends=[(90, 1.5), (90, 1.5)]),
    # a U with a lip: the lip must not be last
    "u_lip": Part("u_lip", t=1.5, width=100, flats=[12, 50, 60, 50], bends=[(90, 1.5)] * 3),
}


def flat(name):
    p = BY.get(name) or EXTRA[name]
    return fm.unfold(fm.recognize(make_solid(p)), k=K)


@pytest.fixture(scope="module")
def flats():
    return {n: flat(n) for n in ("U", "Z", "hat", "hem_open", "z_long_leg", "u_deep", "u_lip")}


def test_simple_profiles_are_feasible(flats):
    for name in ("U", "Z", "hat"):
        plan = fm.plan_bends(flats[name])
        assert plan.feasible and sorted(plan.order) == [b.id for b in flats[name].bends]
    assert fm.plan_bends(flats["Z"]).flips == 1  # up then down: the part is turned over once


def test_order_matters(flats):
    pat = flats["z_long_leg"]
    by_dir = {b.direction: b.id for b in pat.bends}
    long_first = [b.id for b in pat.bends if b.id != by_dir["down"]] + [by_dir["down"]]
    bad = fm.plan_bends(pat, order=long_first)
    (c,) = bad.collisions
    assert not bad.feasible and c.tool in ("die", "beam") and c.bend_id == long_first[1]
    assert c.done == (long_first[0],) and c.describe().startswith(f"bend {c.bend_id} after bends {long_first[0]}")
    good = fm.plan_bends(pat)
    assert good.feasible and good.order == tuple(reversed(long_first))


def test_deep_u_needs_a_gooseneck(flats):
    pat = flats["u_deep"]
    plan = fm.plan_bends(pat, Tooling(PUNCHES["straight_88"], DIES["V10"]))
    assert not plan.feasible and {c.tool for c in plan.collisions} <= {"punch", "clamp"}
    plan = fm.plan_bends(pat, Tooling(PUNCHES["gooseneck_88"], DIES["V10"]))
    assert plan.feasible
    assert any(s.turned for s in plan.steps)  # one leg goes into the throat


def test_lip_order(flats):
    import itertools

    pat = flats["u_lip"]
    tools = Tooling(PUNCHES["straight_88"], DIES["V10"])
    plan = fm.plan_bends(pat, tools)
    assert plan.feasible and fm.plan_bends(pat, tools, order=plan.order).feasible
    verdicts = {o: fm.plan_bends(pat, tools, order=o).feasible for o in itertools.permutations(plan.order)}
    assert 0 < sum(verdicts.values()) < len(verdicts)  # some orders work, some do not


def test_tall_box_hits_the_punch():
    pat = fm.unfold(fm.recognize(box_solid(80, 60, 70, 2.0, 2.0, 2.0)), k=K)
    plan = fm.plan_bends(pat)
    assert not plan.feasible and plan.collisions and {c.tool for c in plan.collisions} <= {"punch", "clamp"}
    assert fm.plan_bends(fm.unfold(fm.recognize(box_solid(80, 60, 30, 2.0, 2.0, 2.0)), k=K)).feasible


def test_hems_are_made_first_and_not_simulated(flats):
    plan = fm.plan_bends(flats["hem_open"])
    assert plan.unchecked == (1,) and plan.feasible and plan.order == ()


def test_check_with_sequence(flats, tmp_path):
    pat = flats["z_long_leg"]
    report = fm.check(pat, sequence=True)
    assert report.sequence.feasible and "DFM_TOOL_COLLISION" not in report.codes()
    by_dir = {b.direction: b.id for b in pat.bends}
    bad = fm.check(pat, sequence=[by_dir["up"], by_dir["down"]])
    (d,) = (d for d in bad.problems if d.code == fm.Code.DFM_TOOL_COLLISION)
    assert d.data["bend"] == by_dir["down"] and "at" in d.data
    # exports: steps in the bend table and the DXF labels
    fm.write_bend_table(pat, tmp_path / "t.json", dfm=report)
    t = json.loads((tmp_path / "t.json").read_text())
    assert t["dfm"]["sequence"]["feasible"]
    assert sorted(b["step"] for b in t["bends"]) == [1, 2]
    fm.write_bend_table(pat, tmp_path / "t.csv", dfm=report)
    rows = list(csv.DictReader(open(tmp_path / "t.csv")))
    assert {r["step"] for r in rows} == {"1", "2"}
    fm.write_dxf(pat, tmp_path / "s.dxf", dfm=report)
    texts = [e.dxf.text for e in ezdxf.readfile(tmp_path / "s.dxf").modelspace().query("TEXT")]
    assert any(x.startswith("(1) ") for x in texts) and any(x.startswith("(2) ") for x in texts)


def test_cli_sequence(tmp_path, capsys):
    from flatmat.cli import main

    part = EXTRA["z_long_leg"]
    from parts import write_step

    write_step(make_solid(part), str(tmp_path / "z.step"))
    assert main(["unfold", str(tmp_path / "z.step"), "--check", "--sequence", "--punch", "straight_88",
                 "--die", "V12"]) == 0  # fmt: skip
    out = capsys.readouterr().out
    assert "bend order:" in out and "straight_88 / V12" in out
    assert main(["unfold", str(tmp_path / "z.step"), "--check", "--punch", "nope"]) == 2
