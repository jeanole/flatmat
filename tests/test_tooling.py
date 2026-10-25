"""Punch and die schema, catalogue, JSON, and the rules that need tools."""

import json
import math

import pytest
from parts import cases, make_solid
from shapely.geometry import Point

import flatmat as fm
from flatmat.tooling import DIES, PUNCHES, Die, Punch, Tooling, die_opening

K = 0.44
BY = {p.name: p for p in cases()}


@pytest.fixture(scope="module")
def flats():
    return {n: fm.unfold(fm.recognize(make_solid(BY[n])), k=K) for n in ("L", "obtuse", "acute", "U")}


@pytest.mark.parametrize("name", list(PUNCHES))
def test_punch_sections(name):
    p = PUNCHES[name]
    body, clamp = p.sections()
    assert body.is_valid and clamp.is_valid
    assert body.bounds[1] == pytest.approx(0.0, abs=1e-9)  # the tip's lowest point is at the origin
    assert body.contains(Point(0, 1.0)) and not body.contains(Point(0, -0.01))
    # beside the tip, below the clamp: free
    leg = Point(20.0, 20.0 * math.tan(math.radians(45)) - 3)
    assert body.contains(leg) is False
    assert clamp.bounds[1] == pytest.approx(body.bounds[3])


def test_gooseneck_has_a_throat():
    p = PUNCHES["gooseneck_88"]
    body, _ = p.sections()
    straight, _ = Punch("s", height=p.height, width=p.width).sections()
    probe = Point(-p.throat_depth / 2, (p.throat_bottom + p.throat_top) / 2 + 10)
    assert not body.contains(probe)  # the recess behind the tip
    assert body.bounds[0] < -p.throat_depth  # the neck leans back past it
    assert body.area > straight.area


def test_die_geometry():
    d = DIES["V12"]
    block, beam = d.sections()
    assert block.is_valid and block.bounds[3] == pytest.approx(0.0)
    assert d.depth == pytest.approx(6 / math.tan(math.radians(44)))
    sy, sz = d.shoulder_center()
    assert block.exterior.distance(Point(sy, sz)) == pytest.approx(d.shoulder_radius, abs=0.01)  # polygonised arc
    assert not block.contains(Point(0.0, -1.0))  # the V
    assert beam.bounds[2] - beam.bounds[0] == pytest.approx(d.base_width)


def test_catalogue_and_json(tmp_path):
    tl = Tooling.named("gooseneck_88", "V16")
    assert tl.die.v == 16 and tl.punch.throat_depth > 0
    assert Tooling.named("straight_88", "V7.5").die.v == 7.5
    tl.save(tmp_path / "t.json")
    assert Tooling.load(tmp_path / "t.json") == tl
    d = json.loads((tmp_path / "t.json").read_text())
    d["punch"]["typo"] = 1
    with pytest.raises(ValueError, match="typo"):
        Tooling.from_dict(d)
    with pytest.raises(KeyError):
        Tooling.named("nope", "V12")
    assert Tooling.default(2.0).die.v == die_opening(2.0) == 12.0
    assert Tooling.from_dict({"die": {"v": 20}}).die == Die(v=20)


def test_v_and_radius_from_the_tools(flats):
    report = fm.check(flats["L"], tooling=Tooling.named("straight_88", "V16"))
    assert report.die_v == 16 and report.bends[0].air_radius == pytest.approx(0.16 * 16)
    big = Tooling(Punch("r6", tip_radius=6.0), DIES["V16"])
    assert fm.check(flats["L"], tooling=big).bends[0].air_radius == pytest.approx(6.0)
    assert fm.check(flats["L"], tooling=big, die_v=24).die_v == 24  # an explicit V wins


def test_tool_angle(flats):
    tl = Tooling.default(2.0)
    report = fm.check(flats["obtuse"], tooling=tl)  # 135°: opening 45° < 88°
    hits = [d for d in report.problems if d.code == fm.Code.DFM_TOOL_ANGLE]
    assert {d.data["tool"] for d in hits} == {"die", "punch"}
    acute = Tooling(PUNCHES["acute_30"], DIES["V16_30"])
    assert "DFM_TOOL_ANGLE" not in fm.check(flats["obtuse"], tooling=acute).codes()
    assert "DFM_TOOL_ANGLE" not in fm.check(flats["L"], tooling=tl).codes()  # 90° with 88° tools: overbend
    assert "DFM_TOOL_ANGLE" not in fm.check(flats["obtuse"]).codes()  # no tools, no rule


def test_tool_length(flats):
    short = Tooling(Punch(length=50.0), DIES["V12"])  # the L is 100 wide
    (d,) = (d for d in fm.check(flats["L"], tooling=short).problems if d.code == fm.Code.DFM_TOOL_LENGTH)
    assert d.data["tool"] == "punch" and d.data["length"] == pytest.approx(100.0)
