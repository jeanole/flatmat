"""Parts flatmat must refuse, with the right code."""

import os

import pytest
from OCP.BRep import BRep_Builder
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepFilletAPI import BRepFilletAPI_MakeChamfer
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeWedge
from OCP.gp import gp_Pnt
from OCP.TopAbs import TopAbs_EDGE
from OCP.TopExp import TopExp_Explorer
from OCP.TopoDS import TopoDS_Compound
from parts import write_step

import flatmat as fm
from flatmat._occ import as_edge


def _roundtrip(shape, step_dir, name):
    p = os.path.join(step_dir, f"{name}.step")
    write_step(shape, p)
    return fm.load_step(p)


def _code(fn):
    with pytest.raises(fm.FlatmatError) as e:
        fn()
    return e.value.diagnostic.code


def test_two_solids(step_dir):
    comp = TopoDS_Compound()
    b = BRep_Builder()
    b.MakeCompound(comp)
    b.Add(comp, BRepPrimAPI_MakeBox(10, 10, 2).Shape())
    b.Add(comp, BRepPrimAPI_MakeBox(gp_Pnt(20, 0, 0), 10, 10, 2).Shape())
    model = _roundtrip(comp, step_dir, "two")
    assert len(model.solids) == 2
    assert _code(lambda: fm.recognize(model)) == fm.Code.MULTIPLE_SOLIDS
    assert fm.recognize(model, solid=1).thickness == pytest.approx(2.0)


def test_tapered_plate(step_dir):
    # a wedge: thickness goes from 4 to 1 along x
    wedge = BRepPrimAPI_MakeWedge(100, 4, 50, 1).Shape()
    model = _roundtrip(wedge, step_dir, "wedge")
    assert _code(lambda: fm.recognize(model)) in (fm.Code.VARIABLE_THICKNESS, fm.Code.NO_THICKNESS)
    # (it must not be reported as an edge feature: thickness is checked first)


def test_chamfered_edge(step_dir):
    box = BRepPrimAPI_MakeBox(100, 50, 3).Shape()
    ch = BRepFilletAPI_MakeChamfer(box)
    ex = TopExp_Explorer(box, TopAbs_EDGE)
    while ex.More():  # a long edge of a skin, not a corner edge (those are fine)
        e = as_edge(ex.Current())
        if BRepAdaptor_Curve(e).LastParameter() - BRepAdaptor_Curve(e).FirstParameter() > 10:
            break
        ex.Next()
    ch.Add(1.0, e)
    model = _roundtrip(ch.Shape(), step_dir, "chamfer")
    assert _code(lambda: fm.recognize(model)) == fm.Code.EDGE_FEATURE


def test_flat_plate_without_bends(step_dir):
    model = _roundtrip(BRepPrimAPI_MakeBox(120, 40, 2).Shape(), step_dir, "plate")
    pat = fm.unfold(fm.recognize(model))
    assert pat.bends == ()
    assert sorted(pat.size) == pytest.approx([40, 120])


def test_unreadable_file(tmp_path):
    p = tmp_path / "x.step"
    p.write_text("not a step file")
    assert _code(lambda: fm.load_step(p)) == fm.Code.READ_FAILED
