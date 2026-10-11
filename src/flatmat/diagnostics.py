"""Structured diagnostics: every refusal and warning names its code and the faces involved."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class Severity(str, Enum):
    ERROR = "error"
    WARNING = "warning"


class Code(str, Enum):
    # input
    READ_FAILED = "READ_FAILED"
    NO_SOLID = "NO_SOLID"
    MULTIPLE_SOLIDS = "MULTIPLE_SOLIDS"
    OPEN_SHELL = "OPEN_SHELL"
    NON_MANIFOLD = "NON_MANIFOLD"
    MISSING_PCURVE = "MISSING_PCURVE"
    INSIDE_OUT_FIXED = "INSIDE_OUT_FIXED"
    # thickness and skins
    NO_THICKNESS = "NO_THICKNESS"
    VARIABLE_THICKNESS = "VARIABLE_THICKNESS"
    NOT_TWO_SKINS = "NOT_TWO_SKINS"
    # geometry
    EDGE_FEATURE = "EDGE_FEATURE"
    UNSUPPORTED_SURFACE = "UNSUPPORTED_SURFACE"
    BEND_TOO_LARGE = "BEND_TOO_LARGE"
    # unfolding
    CYCLE_INCONSISTENT = "CYCLE_INCONSISTENT"
    DEVELOP_FAILED = "DEVELOP_FAILED"
    DEVELOP_FOLDBACK = "DEVELOP_FOLDBACK"
    OUTLINE_OPEN = "OUTLINE_OPEN"
    # flat checks
    FLAT_OVERLAP = "FLAT_OVERLAP"
    AREA_MISMATCH = "AREA_MISMATCH"
    # warnings
    K_EXTRAPOLATED = "K_EXTRAPOLATED"
    K_OUTSIDE_STANDARD = "K_OUTSIDE_STANDARD"
    SHARP_BEND_ASSUMED_RADIUS = "SHARP_BEND_ASSUMED_RADIUS"
    SHARP_BEND_UNSUPPORTED = "SHARP_BEND_UNSUPPORTED"
    COUNTERSINK_PROJECTED = "COUNTERSINK_PROJECTED"
    SURFACE_CONVERTED = "SURFACE_CONVERTED"
    FEATURE_NOT_CIRCULAR = "FEATURE_NOT_CIRCULAR"


@dataclass(frozen=True)
class Diagnostic:
    code: Code
    severity: Severity
    message: str
    face_ids: tuple[int, ...] = ()
    edge_ids: tuple[int, ...] = ()
    data: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code.value,
            "severity": self.severity.value,
            "message": self.message,
            "face_ids": list(self.face_ids),
            "edge_ids": list(self.edge_ids),
            "data": self.data,
        }


class FlatmatError(Exception):
    """Raised when a part is refused. ``diagnostic`` says why and where."""

    def __init__(self, diagnostic: Diagnostic):
        super().__init__(f"{diagnostic.code.value}: {diagnostic.message}")
        self.diagnostic = diagnostic


def refuse(code: Code, message: str, *, faces=(), edges=(), **data: Any) -> FlatmatError:
    return FlatmatError(Diagnostic(code, Severity.ERROR, message, tuple(faces), tuple(edges), dict(data)))


def warn(code: Code, message: str, *, faces=(), edges=(), **data: Any) -> Diagnostic:
    return Diagnostic(code, Severity.WARNING, message, tuple(faces), tuple(edges), dict(data))
