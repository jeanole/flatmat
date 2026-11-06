"""Reading semantic AP242 PMI into PmiItems (plan 7, Task 3)."""

import math
import os

import parts
import pytest
from OCP.BRep import BRep_Builder
from OCP.BRepGProp import BRepGProp
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt, gp_Trsf, gp_Vec
from OCP.GProp import GProp_GProps
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS_Compound
from pmi_parts import centroid, edge_at, face_at, make_doc_with_pmi, write_step_with_pmi

import flatmat as fm
from flatmat.pmi import PmiData, PmiSide, read_pmi, tolerance

# The L_holes case: t = 2, a 50 long base along +x at y in [-1, 1], 100 wide in z, with an
# 8 mm hole at x = 20, z = 30 (axis along y) and a 90° bend up at x = 50.
(L_HOLES,) = (p for p in parts.cases() if p.name == "L_holes")


@pytest.fixture(scope="module")
def bracket():
    """The L_holes solid and the faces the tests refer to, found by geometry."""
    shape = parts.make_solid(L_HOLES)
    return shape, {
        "end": face_at(shape, (0, 0, 50), "plane"),  # the base's cut end, x = 0
        "hole": face_at(shape, (24, 0, 30), "cylinder"),  # the wall of the 8 mm hole
        "skin": face_at(shape, (10, -1, 50), "plane"),  # the base's outer skin, y = -1
        "side": face_at(shape, (10, 0, 0), "plane"),  # the cut side z = 0
        "corner": edge_at(shape, (0, -1, 50)),  # the edge between "end" and "skin"
    }


def _load(tmp_path, shape, items, unit="MM", name="pmi.step"):
    path = str(tmp_path / name)
    write_step_with_pmi(shape, items, path, unit)
    return fm.load_step(path)


def test_linear_dimension_round_trips(tmp_path, bracket):
    shape, f = bracket
    items = [
        dict(
            kind="linear",
            nominal=20.0,
            upper=0.1,
            lower=-0.1,
            faces=([f["end"]], [f["hole"]]),
            points=((0, 0, 30), (16, 0, 30)),
        )
    ]
    data = _load(tmp_path, shape, items).pmi_for(0)
    assert data.status == "semantic"
    assert data.diagnostics == () and data.dropped == 0
    (item,) = data.items
    assert (item.id, item.kind, item.nominal, item.upper, item.lower) == (1, "linear", 20.0, 0.1, -0.1)
    assert item.sides[0].face_ids == (f["end"],) and item.sides[1].face_ids == (f["hole"],)
    assert item.sides[0].points == ((0, 0, 30),) and item.sides[1].points == ((16, 0, 30),)
    assert item.anchor == pytest.approx((8, 0, 30))  # no text point: the midpoint of the side points
    assert item.text == "20 +0.1/-0.1"


def test_inch_file_reads_in_mm(tmp_path, bracket):
    shape, f = bracket
    items = [
        dict(
            kind="linear",
            nominal=1.0,
            upper=0.005,
            lower=-0.005,
            faces=([f["end"]], [f["hole"]]),
            points=((0, 0, 30), (16, 0, 30)),
        )
    ]
    path = str(tmp_path / "inch.step")
    write_step_with_pmi(shape, items, path, unit="INCH")
    with open(path) as fh:
        step = fh.read()
    assert "CONVERSION_BASED_UNIT('INCH'" in step and "POSITIVE_LENGTH_MEASURE(1.)" in step  # the file says 1 in
    model = fm.load_step(path)
    (item,) = model.pmi_for(0).items
    assert item.nominal == pytest.approx(25.4)
    assert item.upper == pytest.approx(0.127) and item.lower == pytest.approx(-0.127)
    assert item.sides[0].face_ids == (f["end"],) and item.sides[1].face_ids == (f["hole"],)
    assert item.sides[1].points[0] == pytest.approx((16, 0, 30))


def test_inch_gdt_zone_reads_in_mm(tmp_path, bracket):
    shape, f = bracket
    items = [
        dict(kind="datum", name="A", faces=([f["skin"]],)),
        dict(kind="position", zone=0.01, datums=("A",), faces=([f["hole"]],)),
    ]
    data = _load(tmp_path, shape, items, unit="INCH").pmi_for(0)
    (tol,) = (i for i in data.items if i.kind == "position")
    assert tol.upper == pytest.approx(0.254)


@pytest.mark.parametrize(
    "upper, lower, expected",
    [(0.1, 0.0, (0.1, 0.0)), (0.0, -0.2, (0.0, -0.2)), (0.2, -0.05, (0.2, -0.05)), (-0.1, 0.1, (None, None))],
    ids=["plus_only", "minus_only", "asymmetric", "inverted"],
)
def test_one_sided_and_inverted_tolerances(tmp_path, bracket, upper, lower, expected):
    shape, f = bracket
    items = [dict(kind="linear", nominal=40.0, upper=upper, lower=lower, faces=([f["end"]], [f["hole"]]))]
    (item,) = _load(tmp_path, shape, items).pmi_for(0).items
    assert (item.nominal, item.upper, item.lower) == (40.0, *expected)
    assert item.text.startswith("40")
    if expected == (None, None):
        assert "tolerance unreadable" in item.text


@pytest.mark.parametrize(
    "nominal, upper, lower, expected",
    [
        (40.0, 0.1, -0.1, (0.1, -0.1)),
        (40.0, math.nan, -0.1, (None, None)),
        (math.nan, 0.1, -0.1, (None, None)),
        (40.0, -0.1, 0.1, (None, None)),
        (40.0, None, None, (None, None)),
    ],
)
def test_tolerance_pairs(nominal, upper, lower, expected):
    assert tolerance(nominal, upper, lower) == expected


def test_anchor_without_points_is_the_centroid_of_the_side_faces(tmp_path, bracket):
    shape, f = bracket
    items = [dict(kind="linear", nominal=20.0, faces=([f["end"]], [f["hole"]]))]
    (item,) = _load(tmp_path, shape, items).pmi_for(0).items
    both = centroid(shape, [f["end"], f["hole"]])
    assert both != pytest.approx(centroid(shape, [f["hole"]]), abs=1)  # the two faces, not the last one
    assert item.anchor == pytest.approx(both)


def test_semantic_and_graphical_copy_is_one_item(tmp_path, bracket):
    shape, f = bracket
    sides = ([f["end"]], [f["hole"]])
    items = [
        dict(kind="linear", nominal=20.0, upper=0.1, lower=-0.1, faces=sides, presentation=True),
        dict(kind="graphical", faces=sides),
    ]
    data = _load(tmp_path, shape, items).pmi_for(0)
    assert data.status == "semantic" and data.diagnostics == ()
    (item,) = data.items
    assert (item.kind, item.nominal) == ("linear", 20.0)
    assert item.anchor is not None  # the presentation's text point


def test_gdt_position_with_datums(tmp_path, bracket):
    shape, f = bracket
    items = [
        dict(kind="datum", name="A", faces=([f["skin"]],)),
        dict(kind="datum", name="B", faces=([f["end"]],)),
        dict(kind="position", zone=0.2, datums=("A", "B"), faces=([f["hole"]],)),
    ]
    data = _load(tmp_path, shape, items).pmi_for(0)
    assert data.status == "semantic"
    (tol,) = (i for i in data.items if i.kind == "position")
    assert (tol.nominal, tol.upper, tol.lower, tol.datums) == (None, 0.2, None, ("A", "B"))
    assert tol.sides == (PmiSide((f["hole"],)),)
    assert tol.anchor == pytest.approx((20, 0, 30), abs=1e-6)  # no points: the centroid of the hole wall
    assert tol.text == "position 0.2 | A | B"
    named = {i.text: i for i in data.items if i.kind == "datum"}
    assert named["A"].sides == (PmiSide((f["skin"],)),) and named["A"].datums == ("A",)
    assert named["B"].sides == (PmiSide((f["end"],)),)
    assert sorted(i.id for i in data.items) == [1, 2, 3]


def test_other_gdt_keeps_its_name(tmp_path, bracket):
    shape, f = bracket
    items = [
        dict(kind="datum", name="A", faces=([f["skin"]],)),
        dict(kind="profile_of_surface", zone=0.5, datums=("A",), faces=([f["end"]],)),
        dict(kind="flatness", zone=0.05, faces=([f["skin"]],)),
    ]
    data = _load(tmp_path, shape, items).pmi_for(0)
    kinds = {i.kind: i for i in data.items}
    assert kinds["other_gdt"].text == "profile of surface 0.5 | A" and kinds["other_gdt"].upper == 0.5
    assert kinds["flatness"].upper == 0.05 and kinds["flatness"].datums == ()


def test_diameter_and_edge_reference(tmp_path, bracket):
    shape, f = bracket
    items = [
        dict(kind="diameter", nominal=8.0, upper=0.1, lower=0.0, faces=([f["hole"]],)),
        dict(kind="linear", nominal=30.0, upper=0.2, lower=-0.2, edges=([f["corner"]],), faces=([], [f["hole"]])),
    ]
    dia, lin = _load(tmp_path, shape, items).pmi_for(0).items
    assert (dia.kind, dia.nominal, dia.upper, dia.lower, dia.text) == ("diameter", 8.0, 0.1, 0.0, "Ø8 +0.1/-0")
    assert dia.sides == (PmiSide((f["hole"],)),)
    assert lin.sides[0].face_ids == tuple(sorted((f["end"], f["skin"])))  # an edge becomes the faces it bounds
    assert lin.sides[1].face_ids == (f["hole"],)


@pytest.fixture(scope="module")
def no_pmi_models(step_dir):
    out = []
    for part in parts.cases():
        path = os.path.join(step_dir, f"pmi_none_{part.name}.step")
        parts.write_step(parts.make_solid(part), path)
        out.append(fm.load_step(path))
    return out


def test_file_without_pmi_has_status_none(no_pmi_models):
    for model in no_pmi_models:
        assert model.pmi_for(0) == PmiData("none", ()), model.path
        assert fm.Code.PMI_READ_FAILED not in [d.code for d in model.diagnostics]


def test_graphical_only_file_warns(tmp_path, bracket):
    shape, f = bracket
    data = _load(tmp_path, shape, [dict(kind="graphical", faces=([f["end"]], [f["hole"]]))]).pmi_for(0)
    assert (data.status, data.items) == ("graphical_only", ())
    (w,) = data.diagnostics
    assert (w.code, w.severity) == (fm.Code.PMI_GRAPHICAL_ONLY, fm.Severity.WARNING)


class _Raises:
    @staticmethod
    def DimTolTool_s(main):
        raise RuntimeError("no GD&T table")


def test_broken_pmi_never_refuses(tmp_path, bracket, monkeypatch):
    shape, f = bracket
    model = _load(tmp_path, shape, [dict(kind="linear", nominal=20.0, faces=([f["end"]], [f["hole"]]))])
    monkeypatch.setattr(fm.pmi, "XCAFDoc_DocumentTool", _Raises)
    data = model.pmi_for(0)
    assert (data.status, data.items) == ("failed", ())
    (w,) = data.diagnostics
    assert w.code is fm.Code.PMI_READ_FAILED and "no GD&T table" in w.message
    assert len(model.solids) == 1 and _volume(model.solids[0]) == pytest.approx(_volume(shape))
    assert fm.recognize(model).thickness == pytest.approx(2.0)


def _volume(s) -> float:
    g = GProp_GProps()
    BRepGProp.VolumeProperties_s(s, g)
    return g.Mass()


def test_no_document_is_none_or_failed(tmp_path, bracket, monkeypatch):
    shape, f = bracket
    path = str(tmp_path / "dim.step")
    write_step_with_pmi(shape, [dict(kind="linear", nominal=20.0, faces=([f["end"]], [f["hole"]]))], path)
    assert fm.load_step(path, pmi=False).pmi_for(0) == PmiData("none", ())

    def boom(path, heal):
        raise RuntimeError("boom")

    monkeypatch.setattr(fm.io, "_read_xcaf", boom)
    model = fm.load_step(path)
    assert model.pmi_for(0) == PmiData("failed", ())  # the load already warned PMI_READ_FAILED, once
    assert [d.code for d in model.diagnostics].count(fm.Code.PMI_READ_FAILED) == 1


def test_pmi_is_cached_per_solid(tmp_path, bracket):
    shape, f = bracket
    model = _load(tmp_path, shape, [dict(kind="linear", nominal=20.0, faces=([f["end"]], [f["hole"]]))])
    assert model.pmi_for(0) is model.pmi_for(0)


def test_items_on_another_solid_are_dropped(tmp_path, bracket):
    shape, f = bracket
    builder, both = BRep_Builder(), TopoDS_Compound()
    builder.MakeCompound(both)
    builder.Add(both, shape)
    builder.Add(both, BRepPrimAPI_MakeBox(gp_Pnt(200, 0, 0), 10, 20, 30).Shape())
    box_top = face_at(both, (205, 20, 15), "plane")  # an id on the compound, not on either solid
    items = [
        dict(kind="linear", nominal=20.0, faces=([face_at(both, (0, 0, 50))], [face_at(both, (24, 0, 30))])),
        dict(kind="linear", nominal=20.0, faces=([box_top], [face_at(both, (205, 0, 15))])),
    ]
    model = _load(tmp_path, both, items)
    assert len(model.solids) == 2
    first, second = model.pmi_for(0), model.pmi_for(1)
    assert (first.status, [i.id for i in first.items], first.dropped) == ("semantic", [1], 1)
    assert first.items[0].sides[0].face_ids == (f["end"],)
    assert ([i.id for i in second.items], second.dropped) == ([2], 1)
    assert second.items[0].sides[0].face_ids == (face_at(model.solids[1], (205, 20, 15)),)


def test_pmi_on_a_placed_part(bracket):
    """PMI defined on a part's own faces still names them when the solid is an instance placed elsewhere."""
    shape, f = bracket
    doc = make_doc_with_pmi(shape, [dict(kind="linear", nominal=20.0, faces=([f["end"]], [f["hole"]]))])
    move = gp_Trsf()
    move.SetTranslation(gp_Vec(100, 0, 0))
    data = read_pmi(doc, shape.Moved(TopLoc_Location(move)))
    assert data.dropped == 0
    (item,) = data.items
    assert item.sides[0].face_ids == (f["end"],) and item.sides[1].face_ids == (f["hole"],)


def test_notes_in_the_document(bracket):
    """OCCT's STEP writer drops notes, so the notes table is filled in memory."""
    shape, f = bracket
    doc = make_doc_with_pmi(
        shape, [dict(kind="note", text="M6 THRU", faces=([f["hole"]],)), dict(kind="note", text="ISO 2768-mK")]
    )
    data = read_pmi(doc, shape)
    thread, general = data.items
    assert (thread.kind, thread.text, thread.sides) == ("note", "M6 THRU", (PmiSide((f["hole"],)),))
    assert thread.anchor == pytest.approx(centroid(shape, [f["hole"]]))
    assert (general.text, general.sides, general.anchor) == ("ISO 2768-mK", (), None)


@pytest.mark.xfail(
    strict=True, raises=AssertionError, reason="OCCT's reader leaves angular tolerances in the file's unit"
)
def test_angle_tolerance_reads_in_degrees(tmp_path, bracket):
    shape, f = bracket
    rad = math.radians
    items = [dict(kind="angle", nominal=rad(90.0), upper=rad(0.5), lower=-rad(0.5), faces=([f["end"]], [f["skin"]]))]
    (item,) = _load(tmp_path, shape, items).pmi_for(0).items  # OCCT writes angles in radians
    assert item.nominal == pytest.approx(90.0)
    assert (item.upper, item.lower) == pytest.approx((0.5, -0.5))
