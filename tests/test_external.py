"""Optional corpora outside the repository.

FLATMAT_BENDFM_DIR    BenDFM example files (<id>.stp, <id>_unfolded.stp, <id>_sequence.json),
                      e.g. github.com/UGent-CVAMO/bendfm example_files. K = 0.45 there.
FLATMAT_PARTS_DIR     any directory of *.step / *.STEP files that must unfold cleanly
                      (private customer parts; never committed).
"""

import glob
import json
import os

import pytest

import flatmat as fm

BENDFM = os.environ.get("FLATMAT_BENDFM_DIR")
PARTS = os.environ.get("FLATMAT_PARTS_DIR")


def _bendfm_ids():
    if not BENDFM:
        return []
    return sorted(int(os.path.basename(p)[:-4]) for p in glob.glob(os.path.join(BENDFM, "*.stp")) if "_" not in p[-14:])


def _top_skin_area(path):
    sheet = fm.recognize(fm.load_step(path))
    return sum(f.area for f in sheet.faces if f.role == "skin_a")


@pytest.mark.skipif(not BENDFM, reason="FLATMAT_BENDFM_DIR not set")
@pytest.mark.parametrize("i", _bendfm_ids() or [None])
def test_bendfm(i):
    seq = json.load(open(os.path.join(BENDFM, f"{i}_sequence.json")))
    pat = fm.unfold_step(os.path.join(BENDFM, f"{i}.stp"), k=0.45)
    assert pat.thickness == pytest.approx(seq[0]["part_thickness"], abs=1e-3)
    assert len(pat.bends) == len(seq) - 1
    theirs = _top_skin_area(os.path.join(BENDFM, f"{i}_unfolded.stp"))
    assert pat.area == pytest.approx(theirs, rel=1e-4)


def _part_files():
    if not PARTS:
        return []
    return sorted(glob.glob(os.path.join(PARTS, "*.step")) + glob.glob(os.path.join(PARTS, "*.STEP")))


@pytest.mark.skipif(not PARTS, reason="FLATMAT_PARTS_DIR not set")
@pytest.mark.parametrize("path", _part_files() or [None], ids=lambda p: os.path.basename(p or "none"))
def test_private_part(path):
    pat = fm.unfold_step(path, k=0.44)
    assert pat.meta["area_rel_err"] < 1e-4
    assert pat.bends
