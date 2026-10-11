"""Sharp bends (r = 0 in the model): sharp inside and round outside, or sharp on both sides."""

import math

import pytest
from conftest import write_and_load
from parts import sharp_cases

import flatmat as fm

K = 0.44
CASES = {p.name: p for p in sharp_cases()}


@pytest.fixture(scope="module")
def sharp(step_dir):
    return {name: (p, fm.recognize(write_and_load(p, step_dir))) for name, p in CASES.items()}


@pytest.mark.parametrize("name", list(CASES))
def test_recognised_and_exact(sharp, name):
    part, sheet = sharp[name]
    assert sheet.thickness == pytest.approx(part.t, abs=1e-6)
    want = sorted((round(abs(a), 6), "sharp" if r == 0 else "bend") for a, r in part.bends)
    assert sorted((round(b.angle_deg, 6), b.kind) for b in sheet.bends) == want
    pat = fm.unfold(sheet, k=K)
    assert pat.area == pytest.approx(part.width * part.flat_length(K), rel=1e-9)
    assert sorted(pat.size) == pytest.approx(sorted((part.width, part.flat_length(K))), abs=1e-6)
    assert pat.meta["area_rel_err"] < 1e-9
    assert len(pat.bend_lines) == len(part.bends)
    for bl in pat.bend_lines:
        assert math.dist(bl.p0, bl.p1) == pytest.approx(part.width, abs=1e-6)


@pytest.mark.parametrize("name", list(CASES))
def test_direction_pattern(sharp, name):
    part, sheet = sharp[name]
    dirs = sorted(b.direction for b in sheet.bends)
    signs = sorted("up" if a > 0 else "down" for a, _ in part.bends)
    flipped = sorted("down" if a > 0 else "up" for a, _ in part.bends)
    assert dirs in (signs, flipped)


def test_bend_line_sits_in_the_middle_of_the_strip(sharp):
    part, sheet = sharp["sharp_both_L"]
    pat = fm.unfold(sheet, k=K)
    (bl,) = pat.bend_lines
    # distance of the bend line from the two short edges of the flat: leg + strip/2
    strip = math.pi / 2 * K * part.t
    x = bl.p0[0] if abs(bl.p0[0] - bl.p1[0]) < 1e-9 else bl.p0[1]
    length = part.flat_length(K)
    assert sorted((x, length - x)) == pytest.approx(sorted((50 + strip / 2, 40 + strip / 2)), abs=1e-6)


@pytest.mark.parametrize("name", ["sharp_in_L", "sharp_both_L", "sharp_both_obtuse"])
def test_sharp_radius_from_tooling(sharp, name):
    part, sheet = sharp[name]
    r = 1.6
    pat = fm.unfold(sheet, k=K, sharp_radius=r)
    ((a, _),) = part.bends
    phi = math.radians(abs(a))
    expected = sum(part.flats) + phi * (r + K * part.t) - 2 * r * math.tan(phi / 2)
    assert max(pat.size) == pytest.approx(expected, abs=1e-6)
    (b,) = pat.bends
    assert b.r_in == r and b.ba == pytest.approx(phi * (r + K * part.t))
    assert not any(d.code == fm.Code.SHARP_BEND_ASSUMED_RADIUS for d in pat.diagnostics)


def test_warns_when_radius_assumed(sharp):
    pat = fm.unfold(sharp["sharp_in_L"][1], k=K)
    assert any(d.code == fm.Code.SHARP_BEND_ASSUMED_RADIUS for d in pat.diagnostics)
