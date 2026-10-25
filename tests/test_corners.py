"""The corner rule: where two bends meet, a relief around the bend zones, and a gap between the folded flanges."""

import pytest
from parts import box_solid

import flatmat as fm

K = 0.44


def box(**kw):
    a = dict(length=80, width=60, height=30, t=2.0, r=2.0, gap=2.0)
    a.update(kw)
    return fm.unfold(fm.recognize(box_solid(**a)), k=K)


@pytest.fixture(scope="module")
def boxes():
    return {
        "ok": box(relief=5.0),  # square reliefs cut into the base's corners
        "tear": box(gap=0.2),  # the short bends end 0.2 mm short of the long ones: no relief
        "no_relief": box(),  # they end 2 mm short, without a relief into the base
        "closed": box(corner_gap=0.1, relief=5.0),  # tabs close the corners to 0.1 mm
        "closed_ok": box(corner_gap=1.0, relief=5.0),
    }


def test_relieved_box_is_clean(boxes):
    assert fm.check(boxes["ok"]).ok


def test_corner_without_relief(boxes):
    report = fm.check(boxes["tear"])
    hits = [d for d in report.problems if d.code == fm.Code.DFM_CORNER_RELIEF]
    assert len(hits) == 4 and {d.code for d in report.problems} == {fm.Code.DFM_CORNER_RELIEF}
    assert all(d.data["fill"] > 0.6 and "other" in d.data and "at" in d.data for d in hits)
    # a bend that stops short of the other one is an ordinary bend end: the relief rule
    assert fm.check(boxes["no_relief"]).codes() == {"DFM_RELIEF_MISSING"}


def test_corner_relief_is_configurable(boxes):
    rules = fm.DfmRules(corner_fill=0.8)
    assert "DFM_CORNER_RELIEF" not in fm.check(boxes["tear"], rules=rules).codes()


def test_closed_corner_gap(boxes):
    hits = [d for d in fm.check(boxes["closed"]).problems if d.code == fm.Code.DFM_CORNER_GAP]
    assert len(hits) == 4
    assert all(d.data["gap"] == pytest.approx(0.1, abs=1e-6) for d in hits)
    assert "DFM_CORNER_GAP" not in fm.check(boxes["closed_ok"]).codes()
    assert "DFM_CORNER_GAP" not in fm.check(boxes["closed"], rules=fm.DfmRules(corner_gap=0.05)).codes()


def test_open_corner_and_profiles_have_no_gap_findings(boxes):
    assert "DFM_CORNER_GAP" not in fm.check(boxes["ok"]).codes()  # 2 mm apart
