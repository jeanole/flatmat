"""Random profiles: the flat length must be exact for any combination."""

import math

from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from parts import Part, make_solid

import flatmat as fm


@st.composite
def profiles(draw):
    n = draw(st.integers(1, 4))
    t = draw(st.sampled_from([0.8, 1.0, 1.5, 2.0, 3.0, 5.0]))
    bends, sharp_outer = [], []
    for i in range(n):
        kind = draw(st.sampled_from(["round", "round", "sharp_in", "sharp_both"]))
        # a sharp outer corner sits t·tan(φ/2) beyond the inner one: keep it inside the flats
        ang = draw(st.floats(15, 120 if kind == "sharp_both" else 165)) * draw(st.sampled_from([1, -1]))
        r = 0.0 if kind != "round" else draw(st.floats(0.5, 8.0))
        if kind == "sharp_both":
            sharp_outer.append(i)
        bends.append((round(ang, 3), round(r, 3)))
    # flats long enough that neighbouring flanges cannot collide in the folded part, and
    # longer than the outer setback t·tan(φ/2) of a sharp corner
    flats = [round(draw(st.floats(30, 60)), 3) for _ in range(n + 1)]
    width = round(draw(st.floats(10, 200)), 3)
    k = draw(st.floats(0.25, 0.5))
    return Part("prop", t, width, flats, bends, sharp_outer=tuple(sharp_outer)), k


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(profiles())
def test_random_profiles_unfold_exactly(case):
    part, k = case
    sheet = fm.recognize(make_solid(part))
    assert math.isclose(sheet.thickness, part.t, abs_tol=1e-6)
    assert len(sheet.bends) == len(part.bends)
    pat = fm.unfold(sheet, k=k)
    assert math.isclose(pat.area, part.width * part.flat_length(k), rel_tol=1e-7)
    w, h = pat.size
    assert math.isclose(max(w, h), max(part.width, part.flat_length(k)), abs_tol=1e-6)
    assert math.isclose(min(w, h), min(part.width, part.flat_length(k)), abs_tol=1e-6)
