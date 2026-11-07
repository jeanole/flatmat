"""B-spline geometry that is really analytic: every test part, written as pure NURBS."""

import math
import os

import pytest
from OCP.BRep import BRep_Tool
from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
from OCP.TopAbs import TopAbs_VERTEX
from parts import Part, cases, feature_cases, make_solid, sharp_cases, write_step

import flatmat as fm
from flatmat._occ import as_vertex
from flatmat.canonical import canonicalize, needs_canonicalization
from flatmat.model import Spline

K = 0.44
ALL = {p.name: p for p in cases() + sharp_cases() + feature_cases()}


def _nurbs_model(p: Part, directory) -> fm.StepModel:
    path = os.path.join(directory, f"nurbs_{p.name}.step")
    write_step(BRepBuilderAPI_NurbsConvert(make_solid(p), True).Shape(), path)
    return fm.load_step(path)


@pytest.mark.parametrize("name", list(ALL))
def test_nurbs_part_unfolds_like_the_analytic_one(name, step_dir):
    p = ALL[name]
    pat = fm.unfold(fm.recognize(_nurbs_model(p, step_dir)), k=K)
    assert len(pat.bends) == len(p.bends)
    assert len(pat.features) == len(p.features)
    assert any(d.code == fm.Code.SURFACE_CONVERTED for d in pat.diagnostics)
    if not any(h[0] == "bend" for h in p.holes):
        through = sum(math.pi * d * d / 4 for kind, *_, d in p.holes) + sum(math.pi * f[4] ** 2 / 4 for f in p.features)
        assert pat.area == pytest.approx(p.width * p.flat_length(K) - through, rel=1e-9)
        # lines and circles come back as LINE / ARC, not as splines
        assert not any(isinstance(s, Spline) for L in (pat.outer, *pat.holes) for s in L.segments)


def test_without_canonicalization_nurbs_is_refused(step_dir):
    model = _nurbs_model(ALL["L"], step_dir)
    with pytest.raises(fm.FlatmatError) as e:
        fm.recognize(model, canonical=False)
    assert e.value.diagnostic.code == fm.Code.UNSUPPORTED_SURFACE


def test_analytic_shapes_are_left_alone():
    solid = make_solid(ALL["U"])
    assert not needs_canonicalization(solid)
    same, diag = canonicalize(solid)
    assert same is solid and diag is None


def test_input_shape_is_not_modified(step_dir):
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    from OCP.GeomAbs import GeomAbs_BSplineSurface
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp_Explorer

    from flatmat._occ import as_face

    solid = _nurbs_model(ALL["L"], step_dir).solid
    converted, diag = canonicalize(solid)
    assert diag is not None and converted is not solid
    ex = TopExp_Explorer(solid, TopAbs_FACE)
    assert BRepAdaptor_Surface(as_face(ex.Current()), False).GetType() == GeomAbs_BSplineSurface


def _face_signatures(shape) -> list[list[tuple[float, float, float]]]:
    """Sorted vertex points of every face, in TopExp.MapShapes order.

    That is the order PMI face references index into. Vertices are exact and survive the
    surface conversion untouched (rounded to 1 um), unlike GProp centroids of trimmed faces, which are
    integrated over rebuilt pcurves and can drift (L_hole_in_bend's hole wall).
    """
    from OCP.TopAbs import TopAbs_FACE
    from OCP.TopExp import TopExp

    from flatmat._occ import ShapeMap

    fmap = ShapeMap()
    TopExp.MapShapes_s(shape, TopAbs_FACE, fmap)
    out = []
    for i in range(1, fmap.Extent() + 1):
        vmap = ShapeMap()
        TopExp.MapShapes_s(fmap.FindKey(i), TopAbs_VERTEX, vmap)
        pts = [BRep_Tool.Pnt_s(as_vertex(vmap.FindKey(k))) for k in range(1, vmap.Extent() + 1)]
        out.append(sorted((round(p.X(), 3), round(p.Y(), 3), round(p.Z(), 3)) for p in pts))
    return out


@pytest.mark.parametrize("name", list(ALL))
def test_canonical_keeps_face_order(name, step_dir):
    """Face i of MapShapes(original) is face i of MapShapes(canonicalize(original)).

    PMI face references index the original solid; flatmat's face ids index the canonical copy.
    """
    original = _nurbs_model(ALL[name], step_dir).solid  # still B-splines: loaded without canonicalization
    assert needs_canonicalization(original)
    canon, diag = canonicalize(original)
    assert canon is not original and diag is not None
    before, after = _face_signatures(original), _face_signatures(canon)
    assert len(before) == len(after)
    for i, (a, b) in enumerate(zip(before, after, strict=True)):
        assert a == b, f"face {i} moved: {a} -> {b}"
    # discriminating: swapping any two faces would be caught, because all signatures differ
    assert len({tuple(a) for a in before}) == len(before)
