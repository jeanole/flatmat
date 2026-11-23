"""SolidCraft's sheet metal reference parts (tests/oracle 67-76): see oracle_parts.py.

Each part is first checked to be SolidCraft's part (its folded volume and area), then
unfolded at SolidCraft's K-factor and at K 0.5 and compared with the flat worked out by
hand. At K 0.5 the neutral fibre is the mid-surface, so the flat area times the
thickness must be the folded volume: SolidCraft's own check, independent of the
hand-worked flat.
"""

import json
import math
import os

import pytest
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from oracle_parts import ba, cases, write_step

import flatmat as fm

CASES = {c.name: c for c in cases()}


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    d = tmp_path_factory.mktemp("oracle")
    out = {}
    for name, c in CASES.items():
        shape = c.build()
        path = os.path.join(d, f"{name}.step")
        write_step(shape, path)
        out[name] = (shape, fm.recognize(fm.load_step(path)), path)
    return out


def _props(shape):
    v, a = GProp_GProps(), GProp_GProps()
    BRepGProp.VolumeProperties_s(shape, v)
    BRepGProp.SurfaceProperties_s(shape, a)
    return v.Mass(), a.Mass()


@pytest.mark.parametrize("name", CASES)
def test_the_part_is_solidcrafts(built, name):
    c = CASES[name]
    volume, area = _props(built[name][0])
    assert volume == pytest.approx(c.volume, rel=1e-9)
    assert area == pytest.approx(c.area, rel=1e-9)


@pytest.mark.parametrize("name", CASES)
def test_recognition(built, name):
    c, sheet = CASES[name], built[name][1]
    assert sheet.thickness == pytest.approx(c.t, abs=1e-6)
    got = sorted((round(b.angle_deg, 6), round(b.r_in, 6)) for b in sheet.bends)
    assert got == sorted((round(a, 6), round(r, 6)) for a, r in c.bends)


@pytest.mark.parametrize("name", CASES)
@pytest.mark.parametrize("k", ["rule", 0.5])
def test_flat_pattern(built, name, k):
    c, sheet = CASES[name], built[name][1]
    k = c.k if k == "rule" else k
    pat = fm.unfold(sheet, k=k)
    assert pat.area == pytest.approx(c.flat_area(k), rel=1e-9)
    assert sorted(pat.size) == pytest.approx(c.flat_size(k), abs=1e-6)
    assert len(pat.holes) == c.holes
    want = sorted(ba(a, r, c.t, k) for a, r in c.bends)
    assert sorted(b.ba for b in pat.bends) == pytest.approx(want, rel=1e-9)


@pytest.mark.parametrize("name", CASES)
def test_flat_area_times_thickness_is_the_volume(built, name):
    c, sheet = CASES[name], built[name][1]
    assert fm.unfold(sheet, k=0.5).area * c.t == pytest.approx(c.volume, rel=1e-9)


@pytest.mark.parametrize("name", CASES)
def test_dxf_and_bend_table(built, name, tmp_path):
    c, sheet = CASES[name], built[name][1]
    pat = fm.unfold(sheet, k=c.k)
    fm.write_dxf(pat, str(tmp_path / "flat.dxf"))
    fm.write_bend_table(pat, str(tmp_path / "flat.json"))
    assert (tmp_path / "flat.dxf").stat().st_size > 0
    assert math.isclose(sum(b.ba for b in pat.bends), sum(ba(a, r, c.t, c.k) for a, r in c.bends), rel_tol=1e-9)


DATA = os.path.join(os.path.dirname(__file__), "data", "solidcraft_oracle")


@pytest.mark.parametrize("name", CASES)
def test_committed_files(name):
    """The STEP files in tests/data/solidcraft_oracle unfold to the flat in their JSON."""
    with open(os.path.join(DATA, f"{name}.json")) as fh:
        want = json.load(fh)
    pat = fm.unfold_step(os.path.join(DATA, f"{name}.step"), k=want["k"])
    assert pat.thickness == pytest.approx(want["thickness"], abs=1e-6)
    assert pat.area == pytest.approx(want["flat"]["area_mm2"], rel=1e-9)
    assert sorted(pat.size) == pytest.approx(want["flat"]["size"], abs=1e-6)
    assert len(pat.holes) == want["flat"]["holes"]
    assert sorted(b.ba for b in pat.bends) == pytest.approx(sorted(b["ba"] for b in want["bends"]), rel=1e-9)
