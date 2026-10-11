"""Generated parts with exactly known flat lengths, each round-tripped through STEP."""

import math

import pytest
from conftest import write_and_load
from parts import Part, cases

import flatmat as fm
from flatmat.model import Arc, Spline

K = 0.44
NAMES = [p.name for p in cases()]


@pytest.mark.parametrize("name", NAMES)
def test_recognition(generated, name):
    part, sheet = generated[name]
    assert sheet.thickness == pytest.approx(part.t, abs=1e-6)
    assert sheet.thickness_confidence > 0.99
    got = sorted((round(b.angle_deg, 6), round(b.r_in, 6)) for b in sheet.bends)
    want = sorted((round(abs(a), 6), round(r, 6)) for a, r in part.bends)
    assert got == want
    # up/down is relative to the reference skin: the pattern must match up to a global flip
    dirs = sorted(b.direction for b in sheet.bends)
    signs = sorted("up" if a > 0 else "down" for a, _ in part.bends)
    flipped = sorted("down" if a > 0 else "up" for a, _ in part.bends)
    assert dirs in (signs, flipped)


@pytest.mark.parametrize("name", NAMES)
def test_flat_area(generated, name):
    part, sheet = generated[name]
    pat = fm.unfold(sheet, k=K)
    if any(h[0] == "bend" for h in part.holes):
        return  # a hole through a bend develops to a non-circular shape: covered by the identity check
    hole_area = sum(math.pi * d * d / 4 for kind, *_, d in part.holes)
    expected = part.width * part.flat_length(K) - hole_area
    assert pat.area == pytest.approx(expected, rel=1e-9)  # areas are exact for lines and arcs
    assert pat.meta["area_rel_err"] < 1e-9


@pytest.mark.parametrize("name", ["L", "U", "Z", "hat", "obtuse", "acute", "hem_open", "joggle"])
def test_flat_size_is_width_by_length(generated, name):
    part, sheet = generated[name]
    w, h = fm.unfold(sheet, k=K).size
    assert sorted((w, h)) == pytest.approx(sorted((part.width, part.flat_length(K))), abs=1e-6)


def test_holes_are_exact_circles(generated):
    pat = fm.unfold(generated["L_holes"][1], k=K)
    assert len(pat.holes) == 2
    for h in pat.holes:
        (seg,) = h.segments
        assert isinstance(seg, Arc) and seg.is_circle
    assert sorted(h.segments[0].radius for h in pat.holes) == pytest.approx([4.0, 5.0])


def test_hole_through_bend_is_an_exact_spline_and_splits_the_bend_line(generated):
    pat = fm.unfold(generated["L_hole_in_bend"][1], k=K)
    (hole,) = pat.holes
    assert any(isinstance(s, Spline) for s in hole.segments)
    assert len(pat.bend_lines) == 2  # interrupted by the hole
    total = sum(math.dist(b.p0, b.p1) for b in pat.bend_lines)
    assert total == pytest.approx(100.0 - 6.0, abs=0.05)  # width minus the hole (~ its diameter)


def test_bend_values(generated):
    pat = fm.unfold(generated["L"][1], k=K)
    (b,) = pat.bends
    assert b.k == K
    assert b.ba == pytest.approx(math.pi / 2 * (2.0 + K * 2.0))
    assert b.ossb == pytest.approx(4.0)
    assert b.bd == pytest.approx(8.0 - b.ba)


def test_bend_line_spans_the_width(generated):
    pat = fm.unfold(generated["U"][1], k=K)
    assert len(pat.bend_lines) == 2
    for bl in pat.bend_lines:
        assert math.dist(bl.p0, bl.p1) == pytest.approx(80.0, abs=1e-6)


def test_din_k_model(generated):
    part, sheet = generated["Z"]
    pat = fm.unfold(sheet, k=fm.KTable.din6935())
    k = fm.din6935_k(4.0 / 3.0) / 2
    assert all(b.k == pytest.approx(k) for b in pat.bends)
    assert pat.meta["area_rel_err"] < 1e-9


def test_given_thickness_is_checked(step_dir):
    model = write_and_load(Part("t_check", t=2.0, width=50, flats=[30, 30], bends=[(90, 2.0)]), step_dir)
    assert fm.recognize(model, thickness=2.0).thickness == 2.0
    with pytest.raises(fm.FlatmatError) as e:
        fm.recognize(model, thickness=2.5)
    assert e.value.diagnostic.code == fm.Code.VARIABLE_THICKNESS
