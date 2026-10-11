"""STEP import (docs/design.md §6.1).

Healing modes:

``fast`` (default)
    The reader's built-in healing with ``FixShellMode`` and ``FixSolidMode`` off.
    ``ShapeFix_Shell`` is quadratic in the face count (4 s of a 4.4 s import on a
    plate with 1000 holes). Switching it off also skips every face- and wire-level
    fix below the shell, so missing pcurves are added here explicitly
    (SolidWorks writes none for cylindrical faces). Solid orientation is checked
    by volume.

    - OCP 8.0: ``XSControl_Reader.SetShapeFixParameters`` with a complete map,
      set *after* ``ReadFile`` (the actor that stores it only exists then).
    - OCP 7.8: the reader takes its sequence from a resource file named ``STEP`` in
      ``$CSF_STEPDefaults``, read once per process at the first STEP read. flatmat
      writes such a file and sets the variable before its first read, unless the
      variable is already set. This affects every STEP read in the process.

``full``
    The reader's default healing.
"""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from typing import Any, Literal

from OCP.BRep import BRep_Tool
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.ShapeFix import ShapeFix_Wire
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID, TopAbs_WIRE
from OCP.TopExp import TopExp_Explorer

from ._occ import OCP_VERSION, as_edge, as_face, as_solid, as_wire
from .diagnostics import Code, Diagnostic, refuse, warn

Heal = Literal["fast", "full"]

_FAST_OVERRIDES = {"FixShape.FixShellMode": "0", "FixShape.FixSolidMode": "0"}
_RESOURCE_78 = """FromSTEP.exec.op : FixShape
FromSTEP.FixShape.Tolerance3d : &Runtime.Tolerance
FromSTEP.FixShape.MaxTolerance3d : &Runtime.MaxTolerance
FromSTEP.FixShape.MinTolerance3d : 1.e-7
FromSTEP.FixShape.FixShellMode : 0
FromSTEP.FixShape.FixSolidMode : 0
ToSTEP.exec.op : SplitCommonVertex,DirectFaces
"""
_HAS_FIX_PARAMETERS = hasattr(STEPControl_Reader, "SetShapeFixParameters")
_state = {"read_once": False, "fast_78": False}


@dataclass(frozen=True)
class StepModel:
    """A STEP file reduced to its solids."""

    solids: tuple[Any, ...]
    path: str
    heal: Heal
    diagnostics: tuple[Diagnostic, ...] = ()

    @property
    def solid(self):
        if len(self.solids) != 1:
            raise refuse(Code.MULTIPLE_SOLIDS, f"{len(self.solids)} solids; pick one with solid=<index>")
        return self.solids[0]


def _prepare_78_fast() -> None:
    if _state["read_once"] or "CSF_STEPDefaults" in os.environ:
        return
    d = tempfile.mkdtemp(prefix="flatmat_occt78_")
    with open(os.path.join(d, "STEP"), "w") as fh:
        fh.write(_RESOURCE_78)
    os.environ["CSF_STEPDefaults"] = d
    _state["fast_78"] = True


def _apply_fast_80(reader) -> None:
    from OCP.DE import DE_ShapeFixParameters
    from OCP.TCollection import TCollection_AsciiString as A

    reader.SetShapeFixParameters(DE_ShapeFixParameters())  # the defaults, as a map
    params = reader.GetShapeFixParameters()
    for k, v in _FAST_OVERRIDES.items():
        params.Bind(A(k), A(v))
    reader.SetShapeFixParameters(params)


def add_missing_pcurves(shape) -> int:
    """Add pcurves the file lacks, in place. Returns the number of wires fixed."""
    fixed = 0
    fe = TopExp_Explorer(shape, TopAbs_FACE)
    while fe.More():
        f = as_face(fe.Current())
        we = TopExp_Explorer(f, TopAbs_WIRE)
        while we.More():
            w = as_wire(we.Current())
            ee = TopExp_Explorer(w, TopAbs_EDGE)
            lacking = False
            while ee.More() and not lacking:
                lacking = BRep_Tool.CurveOnSurface_s(as_edge(ee.Current()), f, 0.0, 0.0) is None
                ee.Next()
            if lacking:
                ShapeFix_Wire(w, f, 1e-7).FixEdgeCurves()
                fixed += 1
            we.Next()
        fe.Next()
    return fixed


def load_step(path: str | os.PathLike, heal: Heal = "fast") -> StepModel:
    path = os.fspath(path)
    diags: list[Diagnostic] = []
    if heal == "fast" and not _HAS_FIX_PARAMETERS:
        _prepare_78_fast()
    if heal == "full" and _state["fast_78"]:
        diags.append(warn(Code.READ_FAILED, "full healing is unavailable on OCP 7.8 after a fast read in this process"))
    reader = STEPControl_Reader()
    status = reader.ReadFile(path)
    _state["read_once"] = True
    if status != IFSelect_RetDone:
        raise refuse(Code.READ_FAILED, f"cannot read STEP file {path}")
    if heal == "fast" and _HAS_FIX_PARAMETERS:
        _apply_fast_80(reader)
    reader.TransferRoots()
    shape = reader.OneShape()
    solids = []
    ex = TopExp_Explorer(shape, TopAbs_SOLID)
    while ex.More():
        solids.append(as_solid(ex.Current()))
        ex.Next()
    if not solids:
        raise refuse(Code.NO_SOLID, f"no solid in {path}")
    if heal == "fast":
        fixed = []
        for s in solids:
            add_missing_pcurves(s)
            g = GProp_GProps()
            BRepGProp.VolumeProperties_s(s, g)
            if g.Mass() < 0:
                s = as_solid(s.Reversed())
                diags.append(warn(Code.INSIDE_OUT_FIXED, "solid was inside out; reversed"))
            fixed.append(s)
        solids = fixed
    return StepModel(tuple(solids), path, heal, tuple(diags))


def ocp_version() -> str:
    return OCP_VERSION
