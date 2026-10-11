import math

import pytest

import flatmat as fm
from flatmat.bendmodel import KScalar, as_k_model, outside_setback


def test_bend_allowance_90():
    assert fm.bend_allowance(math.pi / 2, 2.0, 2.0, 0.44) == pytest.approx(math.pi / 2 * (2 + 0.88))


def test_bend_deduction_90():
    # OSSB = r + t at 90°, BD = 2(r + t) - BA
    ba = fm.bend_allowance(math.pi / 2, 1.0, 2.0, 0.5)
    assert fm.bend_deduction(math.pi / 2, 1.0, 2.0, 0.5) == pytest.approx(2 * 3.0 - ba)


def test_hem_has_no_setback():
    assert outside_setback(math.pi, 1.0, 1.0) is None
    assert fm.bend_deduction(math.pi, 1.0, 1.0, 0.4) is None


@pytest.mark.parametrize(
    ("r_over_s", "k"),
    [(1.0, 0.65), (5.0, 0.65 + 0.5 * math.log10(5)), (5.0001, 1.0), (10, 1.0), (0.01, 0.0)],
)
def test_din_k(r_over_s, k):
    assert fm.din6935_k(r_over_s) == pytest.approx(k)


@pytest.mark.parametrize(
    ("beta", "r", "s", "v"),
    # worked examples of the DIN 6935:2011 text (via modulusmetal), see docs/research.md §6
    [(90, 6, 4, -8.26), (135, 32, 12, -7.25)],
)
def test_din_v_worked_examples(beta, r, s, v):
    assert fm.din6935_v(beta, r, s) == pytest.approx(v, abs=0.005)


def test_din_v_flat_is_zero():
    assert fm.din6935_v(170, 2, 2) == 0.0


def test_ktable_interpolates_between_correct_rows():
    # the FreeCAD SheetMetal table example; FreeCAD's own lookup always used rows 0/1
    t = fm.KTable.from_rows([(1, 0.38), (3, 0.43), (99, 0.5)])
    assert t(2.0, 1.0, 1.0)[0] == pytest.approx(0.405)
    assert t(51.0, 1.0, 1.0)[0] == pytest.approx(0.465)
    k, d = t(200.0, 1.0, 1.0)
    assert k == 0.5 and d.code == fm.Code.K_EXTRAPOLATED


def test_din_model_warns_below_its_range():
    k, d = fm.KTable.din6935()(0.5, 1.0, 1.0)
    assert d.code == fm.Code.K_OUTSIDE_STANDARD
    k, d = fm.KTable.din6935()(2.0, 1.0, 1.0)
    assert d is None and k == pytest.approx((0.65 + 0.5 * math.log10(2)) / 2)


def test_k_model_coercion():
    assert isinstance(as_k_model(0.4), KScalar)
    with pytest.raises(ValueError):
        as_k_model(1.5)
    assert as_k_model(lambda r, t, a: 0.3)(1, 1, 1)[0] == 0.3
