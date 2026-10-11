import csv
import json

import ezdxf
import pytest

import flatmat as fm


@pytest.fixture(scope="module")
def pat(generated):
    return fm.unfold(generated["L_holes"][1], k=0.44)


@pytest.fixture(scope="module")
def pat_spline(generated):
    return fm.unfold(generated["L_hole_in_bend"][1], k=0.44)


def test_dxf_entities_are_exact(tmp_path, pat):
    p = tmp_path / "a.dxf"
    fm.write_dxf(pat, p)
    doc = ezdxf.readfile(p)
    assert not doc.audit().has_errors
    msp = doc.modelspace()
    assert doc.header["$INSUNITS"] == 4
    assert len(msp.query("CIRCLE[layer=='INNER']")) == 2
    assert len(msp.query("LINE[layer=='OUTER']")) == 4  # collinear pieces at bend tangents are joined
    assert len(msp.query("LWPOLYLINE")) == 0
    (bl,) = msp.query("LINE[layer=='BEND_UP' | layer=='BEND_DOWN']")
    assert bl.dxf.linetype == "CENTER"
    assert msp.query("TEXT")[0].dxf.text in ("+90 R2", "-90 R2")


def test_split_circles_are_joined():
    from flatmat.model import Arc, Loop2D
    from flatmat.unfold import _merge_collinear

    halves = Loop2D((Arc((1, 2), 3, 0.0, 3.141592653589793), Arc((1, 2), 3, 3.141592653589793, 6.283185307179586)))
    (c,) = _merge_collinear(halves).segments
    assert c.is_circle and c.center == (1, 2)


def test_dxf_spline(tmp_path, pat_spline):
    p = tmp_path / "s.dxf"
    fm.write_dxf(pat_spline, p, labels=None)
    doc = ezdxf.readfile(p)
    assert not doc.audit().has_errors
    assert len(doc.modelspace().query("SPLINE[layer=='INNER']")) >= 1
    assert len(doc.modelspace().query("TEXT")) == 0


def test_dxf_deterministic(tmp_path, pat):
    a, b = tmp_path / "a.dxf", tmp_path / "b.dxf"
    fm.write_dxf(pat, a, deterministic=True)
    fm.write_dxf(pat, b, deterministic=True)
    assert a.read_bytes() == b.read_bytes()


def test_bend_table_json(tmp_path, pat):
    p = tmp_path / "t.json"
    fm.write_bend_table(pat, p)
    t = json.loads(p.read_text())
    assert t["schema"] == "flatmat.bendtable/1"
    assert t["units"] == "mm" and t["thickness"] == pytest.approx(2.0)
    (b,) = t["bends"]
    assert {"direction", "angle_deg", "opening_angle_deg", "r_in", "k", "ba", "bd", "ossb", "segments"} <= set(b)
    assert b["opening_angle_deg"] == pytest.approx(90.0)


def test_bend_table_csv(tmp_path, pat_spline):
    p = tmp_path / "t.csv"
    fm.write_bend_table(pat_spline, p)
    rows = list(csv.DictReader(p.open()))
    assert len(rows) == 2  # one bend, two segments around the hole
    assert {r["bend"] for r in rows} == {"1"}


def test_svg(tmp_path, pat):
    p = tmp_path / "a.svg"
    fm.write_svg(pat, p)
    assert p.read_text().startswith("<svg")
