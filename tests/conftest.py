import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))

import parts  # noqa: E402

import flatmat as fm  # noqa: E402


@pytest.fixture(scope="session")
def step_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("step")


def write_and_load(part: parts.Part, directory) -> fm.StepModel:
    path = os.path.join(directory, f"{part.name}.step")
    parts.write_step(parts.make_solid(part), path)
    return fm.load_step(path)


@pytest.fixture(scope="session")
def generated(step_dir):
    """name -> (Part, SheetPart), each part round-tripped through STEP."""
    out = {}
    for p in parts.cases():
        out[p.name] = (p, fm.recognize(write_and_load(p, step_dir)))
    return out


def _solid_outcome(model, solid, k: float, dxf: str):
    """The SheetPart, the FlatPattern and the deterministic DXF bytes of one solid, or the
    code flatmat refuses it with."""
    try:
        sheet = fm.recognize(model, solid=solid)
    except fm.FlatmatError as e:
        return ("refused at recognize", e.diagnostic.code)
    try:
        pattern = fm.unfold(sheet, k=k)
    except fm.FlatmatError as e:
        return sheet, ("refused at unfold", e.diagnostic.code)
    fm.write_dxf(pattern, dxf, deterministic=True)
    with open(dxf, "rb") as fh:
        return sheet, pattern, fh.read()


def _outcome(path: str, pmi: bool, k: float, dxf: str, heal: str):
    """What flatmat makes of ``path``: the solid count, the load diagnostics and each solid's
    outcome (the default solid, and each by index when there are several)."""
    try:
        model = fm.load_step(path, heal, pmi=pmi)
    except fm.FlatmatError as e:
        return ("refused at load", e.diagnostic.code)
    picks = [None] + (list(range(len(model.solids))) if len(model.solids) > 1 else [])
    return (
        len(model.solids),
        [d.code for d in model.diagnostics],
        [_solid_outcome(model, i, k, dxf) for i in picks],
    )


def assert_pmi_read_identical(path, out_dir, k: float = 0.44, heal: str = "fast") -> None:
    """Reading ``path`` with PMI on changes nothing: solids (in order), diagnostics, SheetPart,
    FlatPattern and the DXF byte for byte, or the same refusal."""
    path = os.fspath(path)
    old = _outcome(path, False, k, os.path.join(out_dir, "old.dxf"), heal)
    new = _outcome(path, True, k, os.path.join(out_dir, "new.dxf"), heal)
    assert old == new
