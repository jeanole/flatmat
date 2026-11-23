"""Holes across bends and box corners with flats worked out by hand: see handworked_parts.py.

Each part is unfolded at three K-factors and compared with its closed-form flat: area,
size, inner loops and bend allowances. At K 0.5 the flat area times the thickness is the
folded volume (SolidCraft's check), except for a hole drilled radially through a bend,
whose two skins see different holes.
"""

import itertools
import json
import math
import os

import pytest
from handworked_parts import ba, cases, write_step
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps

import flatmat as fm
from flatmat.model import Arc

CASES = {c.name: c for c in cases()}
DATA = os.path.join(os.path.dirname(__file__), "data", "handworked")


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    d = tmp_path_factory.mktemp("handworked")
    out = {}
    for name, c in CASES.items():
        shape = c.build()
        path = os.path.join(d, f"{name}.step")
        write_step(shape, path)
        props = GProp_GProps()
        BRepGProp.VolumeProperties_s(shape, props)
        out[name] = (props.Mass(), fm.recognize(fm.load_step(path)))
    return out


def _developed_skin(c, sheet):
    """For a radial hole: which skin flatmat develops, from its bend face's area."""
    (area,) = (f.area for f in sheet.faces if f.role == "skin_a" and f.kind == "cylinder")
    (skin,) = (s for s in ("inner", "outer") if c.skin_cylinder(s) == pytest.approx(area, rel=1e-4))
    return skin


def _want_area(c, sheet, k):
    want = c.flat_area(k)
    return want[_developed_skin(c, sheet)] if isinstance(want, dict) else want


@pytest.mark.parametrize("name", CASES)
def test_recognition(built, name):
    c, sheet = CASES[name], built[name][1]
    assert sheet.thickness == pytest.approx(c.t, abs=1e-6)
    got = sorted((round(b.angle_deg, 6), round(b.r_in, 6)) for b in sheet.bends)
    assert got == sorted((round(a, 6), round(r, 6)) for a, r in c.bends)


@pytest.mark.parametrize("name", CASES)
@pytest.mark.parametrize("k", [0.33, 0.44, 0.5])
def test_flat_pattern(built, name, k):
    c, sheet = CASES[name], built[name][1]
    pat = fm.unfold(sheet, k=k)
    assert pat.area == pytest.approx(_want_area(c, sheet, k), rel=1e-9)
    assert sorted(pat.size) == pytest.approx(c.flat_size(k), abs=1e-6)
    assert len(pat.holes) == c.holes
    assert sorted(b.ba for b in pat.bends) == pytest.approx(sorted(ba(a, r, c.t, k) for a, r in c.bends), rel=1e-9)


@pytest.mark.parametrize("name", [n for n, c in CASES.items() if c.volume_rel])
def test_flat_area_times_thickness_is_the_volume(built, name):
    c, (volume, sheet) = CASES[name], built[name]
    assert fm.unfold(sheet, k=0.5).area * c.t == pytest.approx(volume, rel=c.volume_rel)


@pytest.mark.parametrize("name", [n for n, c in CASES.items() if c.circle])
def test_a_hole_cut_flat_unfolds_to_its_circle(built, name):
    c, sheet = CASES[name], built[name][1]
    radius, k = c.circle
    (hole,) = fm.unfold(sheet, k=k).holes
    pts = hole.polyline()
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    centre = ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2)
    assert [math.dist(p, centre) for p in pts] == pytest.approx([radius] * len(pts), abs=1e-6)


@pytest.mark.parametrize("name", [n for n, c in CASES.items() if c.hems])
def test_hems(built, name):
    assert sum(b.kind == "hem" for b in built[name][1].bends) == CASES[name].hems


@pytest.mark.parametrize("name", [n for n, c in CASES.items() if c.features])
def test_countersinks_and_counterbores(built, name):
    c, sheet = CASES[name], built[name][1]

    def key(kind, side, through, outer, angle, depth):
        r = lambda x: None if x is None else round(x, 6)  # noqa: E731
        return (kind, side, r(through), r(outer), r(angle), r(depth))

    got = sorted(
        key(f.kind, f.side, f.through_diameter, f.outer_diameter, f.angle_deg, f.depth) for f in sheet.features
    )
    # which skin is "top" is the reference skin's choice: the sides must match up to a flip
    flips = [{1: "top", -1: "bottom"}, {1: "bottom", -1: "top"}]
    assert got in [sorted(key(f[0], m[f[1]], *f[2:]) for f in c.features) for m in flips]


@pytest.mark.parametrize("name", [n for n, c in CASES.items() if c.centres])
@pytest.mark.parametrize("k", [0.33, 0.44, 0.5])
def test_hole_centres(built, name, k):
    """The flat's round holes are where the hand-worked flat puts them: their distances to each
    other, which do not depend on where the flat is placed."""
    c, sheet = CASES[name], built[name][1]
    holes = fm.unfold(sheet, k=k).holes
    assert all(len(h.segments) == 1 and isinstance(h.segments[0], Arc) for h in holes)
    got = sorted(math.dist(a, b) for a, b in itertools.combinations([h.segments[0].center for h in holes], 2))
    want = sorted(math.dist(a, b) for a, b in itertools.combinations(c.centres(k), 2))
    assert got == pytest.approx(want, abs=1e-9)


def test_the_radial_hole_skins_differ(built):
    # the two candidates must be far apart, or the skin check above proves nothing
    c, sheet = CASES["hb_radial"], built["hb_radial"][1]
    want = c.flat_area(0.44)
    assert abs(want["inner"] - want["outer"]) > 1.0
    assert _developed_skin(c, sheet) in want


@pytest.mark.parametrize("name", CASES)
def test_committed_files(name):
    """The STEP files in tests/data/handworked unfold to the flat in their JSON."""
    with open(os.path.join(DATA, f"{name}.json")) as fh:
        want = json.load(fh)
    model = fm.load_step(os.path.join(DATA, f"{name}.step"))
    sheet = fm.recognize(model)
    pat = fm.unfold(sheet, k=want["k"])
    area = want["flat"].get("area_mm2")
    if area is None:
        area = want["flat"]["area_mm2_by_developed_skin"][_developed_skin(CASES[name], sheet)]
    assert pat.thickness == pytest.approx(want["thickness"], abs=1e-6)
    assert pat.area == pytest.approx(area, rel=1e-9)
    assert sorted(pat.size) == pytest.approx(want["flat"]["size"], abs=1e-6)
    assert len(pat.holes) == want["flat"]["holes"]
    assert sorted(b.ba for b in pat.bends) == pytest.approx(sorted(b["ba"] for b in want["bends"]), rel=1e-9)
