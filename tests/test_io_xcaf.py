"""Reading STEP through XCAF (``pmi=True``) gives the same geometry as the plain reader."""

import os

import parts
import pytest
from conftest import assert_pmi_read_identical
from OCP.BRep import BRep_Builder
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetError
from OCP.STEPCAFControl import STEPCAFControl_Reader
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer
from OCP.TopoDS import TopoDS_Compound

import flatmat as fm


@pytest.mark.parametrize("part", parts.cases(), ids=lambda p: p.name)
def test_xcaf_read_is_identical(part, step_dir, tmp_path):
    path = os.path.join(step_dir, f"{part.name}.step")
    parts.write_step(parts.make_solid(part), path)
    assert_pmi_read_identical(path, tmp_path, k=0.44)


@pytest.fixture(scope="module")
def any_case_path(step_dir):
    part = parts.cases()[0]
    path = os.path.join(step_dir, f"xcaf_doc_{part.name}.step")
    parts.write_step(parts.make_solid(part), path)
    return path


def test_xcaf_keeps_the_document(any_case_path):
    assert fm.load_step(any_case_path, pmi=True).doc is not None
    assert fm.load_step(any_case_path, pmi=False).doc is None


def test_pmi_is_on_by_default(any_case_path):
    assert fm.load_step(any_case_path).doc is not None


def test_document_is_not_in_repr(any_case_path):
    assert "doc=" not in repr(fm.load_step(any_case_path, pmi=True))


@pytest.mark.parametrize("name", ["L", "hat", "L_hole_in_bend"])
def test_xcaf_read_is_identical_with_full_healing(name, step_dir, tmp_path):
    (part,) = (p for p in parts.cases() if p.name == name)
    path = os.path.join(step_dir, f"xcaf_full_{name}.step")
    parts.write_step(parts.make_solid(part), path)
    assert_pmi_read_identical(path, tmp_path, k=0.44, heal="full")


@pytest.mark.skipif(not fm.io._HAS_FIX_PARAMETERS, reason="OCP 7.8 heals through a process-wide resource file")
def test_fast_healing_reaches_the_xcaf_reader(any_case_path, monkeypatch):
    seen = []
    real = fm.io._apply_fast_80

    def spy(reader):
        seen.append(type(reader).__name__)
        real(reader)

    monkeypatch.setattr(fm.io, "_apply_fast_80", spy)
    model = fm.load_step(any_case_path, pmi=True)
    assert seen == ["STEPControl_Reader"]
    assert model.doc is not None


class _Reader:
    """The real XCAF reader with some methods replaced."""

    def __init__(self, **replaced):
        self._real = STEPCAFControl_Reader()
        self._replaced = replaced

    def __getattr__(self, name):
        return self._replaced.get(name) or getattr(self._real, name)


def _boom(path, heal):
    raise RuntimeError("boom")


@pytest.mark.parametrize(
    "patch, reason",
    [
        (lambda m: m.setattr(fm.io, "_read_xcaf", _boom), "XCAF read failed: RuntimeError('boom')"),
        (lambda m: m.setattr(fm.io, "_read_xcaf", lambda path, heal: ([], object())), "XCAF read found no solid"),
        (
            lambda m: m.setattr(fm.io, "STEPCAFControl_Reader", lambda: _Reader(ReadFile=lambda p: IFSelect_RetError)),
            "XCAF ReadFile failed",
        ),
        (
            lambda m: m.setattr(fm.io, "STEPCAFControl_Reader", lambda: _Reader(Transfer=lambda doc: False)),
            "XCAF Transfer failed",
        ),
    ],
    ids=["raises", "no_solid", "read_file", "transfer"],
)
def test_xcaf_failure_falls_back_with_a_warning(patch, reason, any_case_path, monkeypatch):
    plain = fm.load_step(any_case_path, pmi=False)
    patch(monkeypatch)
    model = fm.load_step(any_case_path, pmi=True)
    assert model.doc is None
    assert _volumes(model) == _volumes(plain)
    (w,) = (d for d in model.diagnostics if d.code is fm.Code.PMI_READ_FAILED)
    assert w.severity is fm.Severity.WARNING
    assert w.message.startswith(reason) and w.message.endswith("; geometry read without PMI")
    assert [d.code for d in model.diagnostics if d is not w] == [d.code for d in plain.diagnostics]


def test_unreadable_file_is_refused_without_a_pmi_warning(tmp_path):
    path = tmp_path / "garbage.step"
    path.write_text("not a STEP file")
    for pmi in (False, True):
        with pytest.raises(fm.FlatmatError) as e:
            fm.load_step(path, pmi=pmi)
        assert e.value.diagnostic.code is fm.Code.READ_FAILED


def _volumes(model):
    out = []
    for s in model.solids:
        g = GProp_GProps()
        BRepGProp.VolumeProperties_s(s, g)
        out.append(round(g.Mass()))
    return out


def test_xcaf_keeps_the_solid_order_of_several_roots(tmp_path):
    """Three roots of one box each, then a root that is a compound of three boxes."""
    writer = STEPControl_Writer()
    for d in (3, 1, 2):
        writer.Transfer(BRepPrimAPI_MakeBox(d, d, d).Shape(), STEPControl_AsIs)
    builder, compound = BRep_Builder(), TopoDS_Compound()
    builder.MakeCompound(compound)
    for d in (5, 4, 6):
        builder.Add(compound, BRepPrimAPI_MakeBox(d, d, d).Shape())
    writer.Transfer(compound, STEPControl_AsIs)
    path = str(tmp_path / "roots.step")
    writer.Write(path)
    assert _volumes(fm.load_step(path, pmi=False)) == [27, 1, 8, 125, 64, 216]
    assert _volumes(fm.load_step(path, pmi=True)) == [27, 1, 8, 125, 64, 216]
