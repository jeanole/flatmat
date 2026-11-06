"""PMI on the sheet: paths through flanges and bends, threads, general tolerance (plan 7, Task 4)."""

import dataclasses
import math
import os

import parts
import pytest
from OCP.BRep import BRep_Builder
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.gp import gp_Pnt
from OCP.TopoDS import TopoDS_Compound
from pmi_parts import edge_at, face_at, make_doc_with_pmi, write_step_with_pmi

import flatmat as fm
from flatmat.pmi import PmiData, PmiItem, PmiSide, read_pmi
from flatmat.pmi_map import GeneralTolerance, PmiPath, SheetPmi, map_pmi, parse_general_tolerance, parse_thread


def _case(name: str) -> parts.Part:
    (part,) = (p for p in [*parts.cases(), *parts.feature_cases(), *parts.formed_cases()] if p.name == name)
    return part


# ------------------------------------------------------------------------------ parsers


@pytest.mark.parametrize(
    "text, expected",
    [
        ("ISO 2768-m", ("ISO 2768", "m", None)),
        ("ISO 2768-mK", ("ISO 2768", "m", "K")),
        ("Allgemeintoleranzen DIN ISO 2768-f", ("ISO 2768", "f", None)),
        ("DIN 6930-m", ("DIN 6930", "m", None)),
        ("ISO2768 c", ("ISO 2768", "c", None)),
        ("no standard here", None),
    ],
)
def test_parse_general_tolerance(text, expected):
    got = parse_general_tolerance(text)
    assert (got if got is None else (got.standard, got.linear_class, got.geometric_class)) == expected
    if got is not None:
        assert got.text == text


@pytest.mark.parametrize(
    "text, expected",
    [
        ("M6", ("M6", 6.0, None, False, 1)),
        ("M8x1", ("M8x1", 8.0, 1.0, False, 1)),
        ("M6 THRU", ("M6", 6.0, None, True, 1)),
        ("4x M5 durchgehend", ("M5", 5.0, None, True, 4)),
        ("Ø6.6", None),
        ("M", None),
    ],
)
def test_parse_thread(text, expected):
    assert parse_thread(text) == expected


@pytest.mark.parametrize(
    "text, expected",
    [
        ("ISO 2768-1-m", ("ISO 2768", "m", None)),
        ("iso 2768 MK", ("ISO 2768", "m", "K")),
        ("ISO 2768-m Kanten gebrochen", ("ISO 2768", "m", None)),
        ("ISO 2768", None),
        ("ISO 9001", None),
    ],
)
def test_parse_general_tolerance_variants(text, expected):
    got = parse_general_tolerance(text)
    assert (got if got is None else (got.standard, got.linear_class, got.geometric_class)) == expected


@pytest.mark.parametrize(
    "text, expected",
    [
        ("M8 x 1,25", ("M8x1.25", 8.0, 1.25, False, 1)),
        ("20 +0.1/-0.1; M10 through", ("M10", 10.0, None, True, 1)),
        ("4xM5", ("M5", 5.0, None, False, 4)),
        ("6 x M4x0.7 THRU", ("M4x0.7", 4.0, 0.7, True, 6)),
        ("MM", None),
        ("SUM6", None),
        ("PMI 1", None),
    ],
)
def test_parse_thread_variants(text, expected):
    assert parse_thread(text) == expected


# ------------------------------------------------------------------------------ helpers


def _flange_of(sheet, face: int) -> int:
    return next(f.id for f in sheet.flanges if face in f.face_ids)


def _skins_at(sheet, shape, p, q) -> tuple[int, int]:
    """(skin_a face, skin_b face) of the two skins through p and q."""
    a, b = face_at(shape, p, "plane"), face_at(shape, q, "plane")
    role = {f.id: f.role for f in sheet.faces}
    assert {role[a], role[b]} == {"skin_a", "skin_b"}
    return (a, b) if role[a] == "skin_a" else (b, a)


def _bend_between(sheet, fa: int, fb: int) -> int:
    return next(b.id for b in sheet.bends if set(b.flanges) == {fa, fb})


def _mapped(tmp_path, shape, items):
    """``(sheet, [path per item])`` of ``shape`` written with ``items`` and read back."""
    path = os.fspath(tmp_path / "pmi.step")
    write_step_with_pmi(shape, items, path)
    sheet = fm.recognize(fm.load_step(path))
    assert sheet.pmi.status == "semantic"
    return sheet, [p for _, p in sheet.pmi.items]


def _mapped_items(tmp_path, shape, items):
    """The (item, path) pairs, for files whose items come back in another order (datums)."""
    path = os.fspath(tmp_path / "pmi.step")
    write_step_with_pmi(shape, items, path)
    return fm.recognize(fm.load_step(path)).pmi.items


@pytest.fixture(scope="module")
def u_channel():
    """The U case: three flanges and two bends (t = 1.5, 80 wide in z). Flange 1 runs along
    +x at y in [-0.75, 0.75], flange 2 up at x in [31.5, 33], flange 3 back at y in [63.75, 65.25]."""
    shape = parts.make_solid(_case("U"))
    sheet = fm.recognize(shape)
    skins = [
        _skins_at(sheet, shape, (15, -0.75, 40), (15, 0.75, 40)),
        _skins_at(sheet, shape, (31.5, 30, 40), (33, 30, 40)),
        _skins_at(sheet, shape, (15, 63.75, 40), (15, 65.25, 40)),
    ]
    flanges = [_flange_of(sheet, a) for a, _ in skins]
    assert len(set(flanges)) == 3 and len(sheet.bends) == 2
    return shape, sheet, skins, flanges


@pytest.fixture(scope="module")
def bracket():
    """The L_holes case: t = 2, a 50 long base along +x at y in [-1, 1], 100 wide in z, an
    8 mm hole at x = 20, z = 30, a 90° bend up at x = 50 to a flange at x in [52, 54]."""
    shape = parts.make_solid(_case("L_holes"))
    sheet = fm.recognize(shape)
    base_a, base_b = _skins_at(sheet, shape, (10, -1, 50), (10, 1, 50))
    flange_a, _ = _skins_at(sheet, shape, (52, 20, 50), (54, 20, 50))
    return (
        shape,
        sheet,
        {
            "end": face_at(shape, (0, 0, 50), "plane"),  # the base's cut end, x = 0
            "side": face_at(shape, (10, 0, 0), "plane"),  # the cut side z = 0, along the whole part
            "hole": face_at(shape, (24, 0, 30), "cylinder"),  # the 8 mm hole's wall
            "base_a": base_a,
            "base_b": base_b,
            "flange_a": flange_a,
            "base": _flange_of(sheet, base_a),
            "flange": _flange_of(sheet, flange_a),
        },
    )


# ------------------------------------------------------------------------------ paths


def test_same_flange_is_laser(tmp_path, bracket):
    shape, _, f = bracket
    items = [dict(kind="linear", nominal=16.0, upper=0.1, lower=-0.1, faces=([f["end"]], [f["hole"]]))]
    sheet, (path,) = _mapped(tmp_path, shape, items)
    assert path == PmiPath("laser", (f["base"],), ())


def test_skin_b_face_maps_to_its_flange(tmp_path, bracket):
    shape, _, f = bracket
    items = [dict(kind="linear", nominal=2.0, faces=([f["base_b"]], [f["end"]]))]
    _, (path,) = _mapped(tmp_path, shape, items)
    assert path == PmiPath("laser", (f["base"],))


def test_hole_diameter_is_laser_hole(tmp_path, bracket):
    shape, _, f = bracket
    _, (path,) = _mapped(tmp_path, shape, [dict(kind="diameter", nominal=8.0, faces=([f["hole"]],))])
    assert path == PmiPath("laser_hole", (f["base"],), feature_id=None)

    csk = parts.make_solid(_case("csk_L"))  # a 6 mm through hole under a countersink at x = 24, z = 20
    wall = face_at(csk, (27, -1.0, 20), "cylinder")
    sheet, (path,) = _mapped(tmp_path, csk, [dict(kind="diameter", nominal=6.0, faces=([wall],))])
    (feature,) = (x for x in sheet.features if wall in x.face_ids)
    assert feature.kind == "countersink"
    assert path.kind == "laser_hole" and path.feature_id == feature.id
    assert path.flanges == (_flange_of(sheet, _skins_at(sheet, csk, (10, -1.5, 50), (10, 1.5, 50))[0]),)


def test_across_two_bends_is_bend_chain(tmp_path, u_channel):
    shape, sheet, skins, (f1, f2, f3) = u_channel
    outer_1 = face_at(shape, (15, -0.75, 40), "plane")
    outer_3 = face_at(shape, (15, 65.25, 40), "plane")
    _, (path,) = _mapped(tmp_path, shape, [dict(kind="linear", nominal=66.0, faces=([outer_1], [outer_3]))])
    assert path == PmiPath("bend_chain", (f1, f3), (_bend_between(sheet, f1, f2), _bend_between(sheet, f2, f3)))


def test_angle_between_adjacent_flanges(tmp_path, u_channel):
    shape, sheet, skins, (f1, f2, _) = u_channel
    items = [dict(kind="angle", nominal=math.radians(90.0), faces=([skins[0][0]], [skins[1][0]]))]
    _, (path,) = _mapped(tmp_path, shape, items)
    assert path == PmiPath("bend_angle", (f1, f2), (_bend_between(sheet, f1, f2),))


def test_angle_across_two_bends_is_unmapped(tmp_path, u_channel):
    shape, _, skins, _ = u_channel
    items = [dict(kind="angle", nominal=math.radians(180.0), faces=([skins[0][0]], [skins[2][0]]))]
    _, (path,) = _mapped(tmp_path, shape, items)
    assert path == PmiPath("unmapped", reason="angle across 2 bends")


def test_side_on_bend_face_is_unmapped(tmp_path, u_channel):
    shape, sheet, _, _ = u_channel
    end_1 = face_at(shape, (0, 0, 40), "plane")
    bend_face = sheet.bends[0].segments[0].face_id
    _, (path,) = _mapped(tmp_path, shape, [dict(kind="linear", nominal=30.0, faces=([end_1], [bend_face]))])
    assert path == PmiPath("unmapped", reason="on a bend")


def test_perpendicularity_is_bend_angle(tmp_path, u_channel):
    shape, sheet, skins, (f1, f2, _) = u_channel
    items = [
        dict(kind="datum", name="A", faces=([skins[0][0]],)),
        dict(kind="perpendicularity", zone=0.1, datums=("A",), faces=([skins[1][0]],)),
    ]
    by_kind = {item.kind: path for item, path in _mapped_items(tmp_path, shape, items)}
    assert by_kind == {
        "datum": PmiPath("unmapped", reason="datum"),
        "perpendicularity": PmiPath("bend_angle", (f2, f1), (_bend_between(sheet, f1, f2),)),
    }


def test_flatness_on_one_flange(tmp_path, u_channel):
    shape, _, skins, (f1, _, _) = u_channel
    _, (path,) = _mapped(tmp_path, shape, [dict(kind="flatness", zone=0.2, faces=([skins[0][0]],))])
    assert path == PmiPath("flatness", (f1,))


def test_unsupported_gdt_is_unmapped(tmp_path, u_channel):
    shape, _, skins, _ = u_channel
    _, (path,) = _mapped(tmp_path, shape, [dict(kind="profile_of_surface", zone=0.3, faces=([skins[1][0]],))])
    assert path == PmiPath("unmapped", reason="profile of surface is not checked")


def test_position_takes_the_path_to_its_datum(tmp_path, bracket):
    shape, sheet, f = bracket
    items = [
        dict(kind="datum", name="A", faces=([f["flange_a"]],)),
        dict(kind="datum", name="B", faces=([f["end"]],)),
        dict(kind="position", zone=0.2, datums=("B", "A"), faces=([f["hole"]],)),
    ]
    (path,) = (p for item, p in _mapped_items(tmp_path, shape, items) if item.kind == "position")
    assert path == PmiPath("bend_chain", (f["base"], f["flange"]), (_bend_between(sheet, f["base"], f["flange"]),))


# ------------------------------------------------------------------------------ in memory


def _in_memory(shape, items) -> SheetPmi:
    """Notes do not survive OCCT's STEP round trip: map a document built in memory."""
    return map_pmi(fm.recognize(shape), read_pmi(make_doc_with_pmi(shape, items), shape, trust_xcaf=True))


def test_thread_from_hole_callout(bracket):
    shape, _, f = bracket
    pmi = _in_memory(shape, [dict(kind="note", text="M6 THRU", faces=([f["hole"]],))])
    (thread,) = pmi.threads
    assert (thread.designation, thread.nominal_mm, thread.pitch_mm, thread.through) == ("M6", 6.0, None, True)
    assert thread.count == 1
    assert thread.face_ids == (f["hole"],)
    assert pmi.notes == ("M6 THRU",) and pmi.items == ()


def test_general_tolerance_from_note(bracket):
    shape, _, _ = bracket
    pmi = _in_memory(shape, [dict(kind="note", text="ISO 2768-mK")])
    assert pmi.general == GeneralTolerance("ISO 2768", "m", "K", "ISO 2768-mK")
    assert pmi.notes == ("ISO 2768-mK",) and pmi.threads == ()


def test_thread_in_a_dimension_callout(bracket):
    """Threads come from every text, not only notes: here a diameter's callout. The tap holds
    the tapped size, not the laser, so the diameter itself is unmapped."""
    shape, sheet, f = bracket
    text = "M6x1.0 - 6H THRU"
    item = PmiItem(1, "diameter", 6.0, 0.1, 0.0, sides=(PmiSide((f["hole"], f["base_a"])),), text=text)
    pmi = map_pmi(sheet, PmiData("semantic", (item,)))
    (thread,) = pmi.threads
    assert (thread.designation, thread.pitch_mm, thread.through, thread.face_ids) == ("M6x1.0", 1.0, True, (f["hole"],))
    assert pmi.items == ((item, PmiPath("unmapped", reason="a thread")),)


def test_unreadable_tolerance_keeps_its_path(bracket):
    _, sheet, f = bracket
    sides = (PmiSide((f["end"],)), PmiSide((f["hole"],)))
    items = (
        PmiItem(1, "linear", 20.0, None, None, sides=sides, text="20 +0.1/-0.1 (unit unknown)"),
        PmiItem(2, "linear", None, None, None, sides=sides, text="40–40.2 (limits)"),
    )
    pmi = map_pmi(sheet, PmiData("semantic", items))
    assert [p for _, p in pmi.items] == [PmiPath("laser", (f["base"],))] * 2


def test_sheet_thickness_is_unmapped(bracket):
    _, sheet, f = bracket
    item = PmiItem(1, "linear", 2.0, 0.1, -0.1, sides=(PmiSide((f["base_a"],)), PmiSide((f["base_b"],))))
    ((_, path),) = map_pmi(sheet, PmiData("semantic", (item,))).items
    assert path == PmiPath("unmapped", reason="the sheet thickness")


def test_hole_in_a_bend_is_on_the_bend():
    shape = parts.make_solid(_case("L_hole_in_bend"))
    sheet = fm.recognize(shape)
    (wall,) = (x.id for x in sheet.faces if x.role == "thickness" and x.kind == "cylinder")
    item = PmiItem(1, "diameter", 6.0, 0.1, -0.1, sides=(PmiSide((wall,)),))
    ((_, path),) = map_pmi(sheet, PmiData("semantic", (item,))).items
    assert path == PmiPath("unmapped", reason="on a bend")


def test_side_on_a_formed_feature_is_unmapped():
    shape = parts.make_solid(_case("emboss_up"))
    sheet = fm.recognize(shape)
    (emboss,) = sheet.features
    end = face_at(shape, (0, 0, 50), "plane")
    item = PmiItem(1, "linear", 25.0, None, None, sides=(PmiSide((end,)), PmiSide((emboss.face_ids[0],))))
    ((_, path),) = map_pmi(sheet, PmiData("semantic", (item,))).items
    assert path == PmiPath("unmapped", reason="on an emboss")


# ------------------------------------------------------------------------------ files


def test_multi_solid_drops_foreign_items(tmp_path, bracket):
    shape, _, _ = bracket
    builder, both = BRep_Builder(), TopoDS_Compound()
    builder.MakeCompound(both)
    builder.Add(both, shape)
    builder.Add(both, BRepPrimAPI_MakeBox(gp_Pnt(200, 0, 0), 10, 20, 30).Shape())
    items = [dict(kind="linear", nominal=20.0, faces=([face_at(both, (205, 20, 15))], [face_at(both, (205, 0, 15))]))]
    path = os.fspath(tmp_path / "two.step")
    write_step_with_pmi(both, items, path)
    sheet = fm.recognize(fm.load_step(path), solid=0)
    assert (sheet.pmi.status, sheet.pmi.items, sheet.pmi.dropped) == ("semantic", (), 1)


def test_no_pmi_changes_nothing(tmp_path):
    for part in parts.cases():
        path = os.fspath(tmp_path / f"{part.name}.step")
        parts.write_step(parts.make_solid(part), path)
        with_pmi = fm.recognize(fm.load_step(path))
        without = fm.recognize(fm.load_step(path, pmi=False))
        assert with_pmi.pmi == without.pmi == SheetPmi("none"), part.name
        assert with_pmi.diagnostics == without.diagnostics, part.name
        assert with_pmi == without, part.name


def test_pmi_diagnostics_reach_the_sheet(tmp_path, bracket):
    """A graphical-only file: status and its one warning come through recognize."""
    shape, _, f = bracket
    path = os.fspath(tmp_path / "graphical.step")
    write_step_with_pmi(shape, [dict(kind="graphical", faces=([f["end"]], [f["hole"]]))], path)
    sheet = fm.recognize(fm.load_step(path))
    assert sheet.pmi == SheetPmi("graphical_only")
    assert [d.code for d in sheet.diagnostics].count(fm.Code.PMI_GRAPHICAL_ONLY) == 1


def test_mapping_error_never_refuses(tmp_path, bracket, monkeypatch):
    shape, _, f = bracket
    path = os.fspath(tmp_path / "pmi.step")
    write_step_with_pmi(shape, [dict(kind="linear", nominal=20.0, faces=([f["end"]], [f["hole"]]))], path)

    def boom(*a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr("flatmat.pmi_map.map_pmi", boom)
    sheet = fm.recognize(fm.load_step(path))
    assert sheet.pmi.status == "failed"
    assert [d.code for d in sheet.diagnostics] == [fm.Code.PMI_READ_FAILED]
    assert sheet.flanges and fm.unfold(sheet).area > 0


def test_joggle_chain_crosses_both_bends():
    """A joggle's two bends touch with no flange between them: the chain still crosses both."""
    sheet = fm.recognize(parts.make_solid(_case("joggle")))
    first, second = sheet.flanges
    sides = (PmiSide(first.face_ids[:1]), PmiSide(second.face_ids[:1]))
    ((_, path),) = map_pmi(sheet, PmiData("semantic", (PmiItem(1, "linear", 10.0, None, None, sides=sides),))).items
    (b1,) = (b.id for b in sheet.bends if b.flanges == (first.id,))  # each bend touches one flange
    (b2,) = (b.id for b in sheet.bends if b.flanges == (second.id,))
    assert path == PmiPath("bend_chain", (first.id, second.id), (b1, b2))


def test_no_path_crosses_a_lance():
    """A lance's tab is bent in the punching tool, not on the press brake."""
    cut = parts.cut_formed_cases()
    sheet = fm.recognize(parts.make_solid(next(p for p in cut if p.name == "lance")))
    (lance,) = sheet.features
    base = sheet.flanges[0].face_ids[0]
    items = tuple(
        PmiItem(n, kind, 1.0, None, None, sides=(PmiSide((base,)), PmiSide((f.id,))))
        for n, (kind, f) in enumerate(((k, f) for k in ("linear", "angle") for f in sheet.faces), 1)
    )
    paths = [p for _, p in map_pmi(sheet, PmiData("semantic", items)).items]
    assert all(lance.bend_id not in p.bends for p in paths)
    assert PmiPath("unmapped", reason="across a lance") in paths and PmiPath("unmapped", reason="on a lance") in paths


def _paths(sheet, *items) -> list[PmiPath]:
    return [p for _, p in map_pmi(sheet, PmiData("semantic", tuple(items))).items]


def test_countersink_cone_and_counterbore_step_are_unmapped():
    """Only the through-bore is laser-cut: a side on a countersink's cone, or on a
    counterbore's bore or step, is not."""
    csk = parts.make_solid(_case("csk_L"))  # t = 3; a 90° cone from Ø10 at y = 1.5 to Ø6 at y = -0.5, x = 24, z = 20
    sheet = fm.recognize(csk)
    end, cone = face_at(csk, (0, 0, 50), "plane"), face_at(csk, (28, 0.5, 20))
    assert next(x for x in sheet.features if cone in x.face_ids).kind == "countersink"
    linear = PmiItem(1, "linear", 24.0, None, None, sides=(PmiSide((end,)), PmiSide((cone,))))
    diameter = PmiItem(2, "diameter", 10.0, None, None, sides=(PmiSide((cone,)),))
    assert _paths(sheet, linear, diameter) == [PmiPath("unmapped", reason="on a countersink")] * 2

    cbore = parts.make_solid(_case("cbore_L"))  # t = 4; Ø9 bore 2 deep from y = 2 over a Ø5 hole, x = 24, z = 20
    sheet = fm.recognize(cbore)
    bore, step = face_at(cbore, (28.5, 1, 20), "cylinder"), face_at(cbore, (27, 0, 20), "plane")
    items = [PmiItem(1, "diameter", 9.0, None, None, sides=(PmiSide((bore,)),))]
    items.append(
        PmiItem(2, "linear", 2.0, None, None, sides=(PmiSide((step,)), PmiSide((face_at(cbore, (0, 0, 50)),))))
    )
    assert _paths(sheet, *items) == [PmiPath("unmapped", reason="on a counterbore")] * 2


def test_counterbore_through_bore_is_laser_hole():
    cbore = parts.make_solid(_case("cbore_L"))
    sheet = fm.recognize(cbore)
    wall = face_at(cbore, (26.5, -1, 20), "cylinder")  # the Ø5 through-bore under the counterbore
    (feature,) = (x for x in sheet.features if wall in x.face_ids)
    base = _flange_of(sheet, _skins_at(sheet, cbore, (10, -2, 50), (10, 2, 50))[0])
    (path,) = _paths(sheet, PmiItem(1, "diameter", 5.0, 0.1, 0.0, sides=(PmiSide((wall,)),)))
    assert path == PmiPath("laser_hole", (base,), feature_id=feature.id)


def _note_thread(shape, text, walls):
    (thread,) = _in_memory(shape, [dict(kind="note", text=text, faces=(list(walls),))]).threads
    return thread


def test_thread_count_from_the_callout(bracket):
    shape, _, f = bracket
    thread = _note_thread(shape, "4x M5", [f["hole"]])
    assert (thread.designation, thread.count, thread.face_ids) == ("M5", 4, (f["hole"],))
    assert _note_thread(shape, "M6", [f["hole"]]).count == 1


def test_thread_count_from_the_holes_it_points_at():
    """L_holes with three holes in its base: an "M6" pointing at all three walls is three threads."""
    part = dataclasses.replace(_case("L_holes"), holes=[("flat", 0, 0.4, z, 5.0) for z in (20.0, 50.0, 80.0)])
    shape = parts.make_solid(part)
    walls = [face_at(shape, (22.5, 0, z), "cylinder") for z in (20.0, 50.0, 80.0)]
    thread = _note_thread(shape, "M6", walls)
    assert (thread.count, thread.face_ids) == (3, tuple(sorted(walls)))
    assert _note_thread(shape, "2x M6", walls).count == 3  # the larger of the two


def test_coaxial_walls_are_one_hole():
    """A counterbore's bore and its through-bore are two walls of one hole."""
    cbore = parts.make_solid(_case("cbore_L"))
    walls = [face_at(cbore, (28.5, 1, 20), "cylinder"), face_at(cbore, (26.5, -1, 20), "cylinder")]
    assert _note_thread(cbore, "M5", walls).count == 1


def _rim(r: float, y: float, cx: float = 20.0, cz: float = 30.0, deg: float = 45.0):
    """A point on the rim of the L_holes 8 mm hole (centre x = cx, z = cz) on the skin y."""
    a = math.radians(deg)
    return (cx + r * math.cos(a), y, cz + r * math.sin(a))


def test_hole_diameter_through_an_edge(tmp_path, bracket):
    """A diameter on the hole's rim edge reads as {hole wall, skin plane}: it is the hole."""
    shape, _, f = bracket
    rim = edge_at(shape, _rim(4.0, -1.0))
    items = [dict(kind="diameter", nominal=8.0, upper=0.1, lower=0.0, edges=([rim],))]
    sheet, (path,) = _mapped(tmp_path, shape, items)
    ((item, _),) = sheet.pmi.items
    ids = item.sides[0].face_ids  # the edge's two faces: the hole wall and one skin
    assert len(ids) == 2 and f["hole"] in ids and set(ids) - {f["hole"]} <= {f["base_a"], f["base_b"]}
    assert path == PmiPath("laser_hole", (f["base"],))


def test_thickness_through_edges_is_unmapped(tmp_path, bracket):
    """The two edges of the base's cut end: {skin_a, end} against {skin_b, end} is the thickness."""
    shape, _, f = bracket
    items = [
        dict(
            kind="linear",
            nominal=2.0,
            upper=0.1,
            lower=-0.1,
            edges=([edge_at(shape, (0, -1, 50))], [edge_at(shape, (0, 1, 50))]),
        )
    ]
    _, (path,) = _mapped(tmp_path, shape, items)
    assert path == PmiPath("unmapped", reason="the sheet thickness")


def test_side_on_several_flanges_is_unmapped(bracket):
    """The cut side z = 0 runs along the base, the bend and the flange: which flange a
    dimension to it starts from cannot be told."""
    _, sheet, f = bracket
    reason = f"a side spans flanges {min(f['base'], f['flange'])} and {max(f['base'], f['flange'])}"
    side, end = PmiSide((f["side"],)), PmiSide((f["end"],))
    items = (
        PmiItem(1, "linear", 100.0, None, None, sides=(side,)),
        PmiItem(2, "position", None, 0.1, None, sides=(side,)),
        PmiItem(3, "linear", 50.0, None, None, sides=(end, side)),
    )
    assert _paths(sheet, *items) == [PmiPath("unmapped", reason=reason)] * 3
