"""The OCP names flatmat uses that OCP spells awkwardly, in one place (OCP 8.0).

OCP 8.0 puts the NCollection instantiations in ``OCP.collections`` under long template
names, and its TopoDS downcasts have no ``_s`` suffix (``TopoDS.Face``).
"""

from __future__ import annotations

import OCP
from OCP.collections import (
    IndexedDataMap_TopoDS_Shape_List_TopoDS_Shape_TopTools_ShapeMapHasher as EdgeFaceMap,
)
from OCP.collections import IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher as ShapeMap
from OCP.collections import List_TopoDS_Shape as ShapeList
from OCP.collections import Sequence_TDF_Label as LabelSequence
from OCP.TopoDS import TopoDS

as_face = TopoDS.Face
as_edge = TopoDS.Edge
as_vertex = TopoDS.Vertex
as_wire = TopoDS.Wire
as_shell = TopoDS.Shell
as_solid = TopoDS.Solid

OCP_VERSION: str = getattr(OCP, "__version__", "unknown")

__all__ = [
    "EdgeFaceMap",
    "ShapeMap",
    "ShapeList",
    "LabelSequence",
    "as_face",
    "as_edge",
    "as_vertex",
    "as_wire",
    "as_shell",
    "as_solid",
    "OCP_VERSION",
]
