"""flatmat: sheet metal recognition and unfolding for STEP files.

>>> import flatmat as fm
>>> pat = fm.unfold_step("bracket.step", k=0.44)
>>> fm.write_dxf(pat, "bracket.dxf")
>>> fm.write_bend_table(pat, "bracket.json")
"""

from __future__ import annotations

__version__ = "0.1.0.dev0"

from .bendmodel import KDin6935, KScalar, KTable, bend_allowance, bend_deduction, din6935_k, din6935_v
from .dfm import MATERIALS, DfmReport, DfmRules, Material, check
from .diagnostics import Code, Diagnostic, FlatmatError, Severity
from .export import LayerMap, bend_table, write_bend_table, write_dxf, write_svg
from .io import StepModel, load_step
from .model import Arc, Bend, BendLine, Feature, FlatPattern, Line, Loop2D, SheetPart, Spline, Tolerances
from .recognize import recognize
from .unfold import unfold


def unfold_step(path, k=0.44, *, heal="fast", solid=None, thickness=None, tol=None, sharp_radius=0.0) -> FlatPattern:
    """Read, recognise and unfold in one call."""
    sheet = recognize(load_step(path, heal=heal), solid=solid, thickness=thickness, tol=tol)
    return unfold(sheet, k=k, sharp_radius=sharp_radius)


__all__ = [
    "__version__",
    "Arc",
    "Bend",
    "BendLine",
    "Code",
    "DfmReport",
    "DfmRules",
    "MATERIALS",
    "Material",
    "Diagnostic",
    "Feature",
    "FlatPattern",
    "FlatmatError",
    "KDin6935",
    "KScalar",
    "KTable",
    "LayerMap",
    "Line",
    "Loop2D",
    "Severity",
    "SheetPart",
    "Spline",
    "StepModel",
    "Tolerances",
    "bend_allowance",
    "bend_deduction",
    "bend_table",
    "check",
    "din6935_k",
    "din6935_v",
    "load_step",
    "recognize",
    "unfold",
    "unfold_step",
    "write_bend_table",
    "write_dxf",
    "write_svg",
]
