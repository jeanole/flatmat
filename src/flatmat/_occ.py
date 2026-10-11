"""The few places where OCP 7.8 and 8.0 differ.

Found by the phase-0 spike (spike/api_probe.py, spike/README.md):

- 8.0 moved the NCollection instantiations to ``OCP.collections`` with long
  template names; 7.8 keeps the ``TopTools_*`` typedefs.
- 8.0 dropped the ``_s`` suffix on the TopoDS downcasts (``TopoDS.Face``).
- ``BRepTopAdaptor_FClass2d`` is gone in 8.0; ``IntTools_FClass2d`` exists in both.
- STEP healing is configured differently (see :mod:`flatmat.io`).
"""

from __future__ import annotations

import OCP
from OCP.TopoDS import TopoDS

try:  # OCP 7.8
    from OCP.TopTools import (  # type: ignore[attr-defined]
        TopTools_IndexedDataMapOfShapeListOfShape as EdgeFaceMap,
    )
    from OCP.TopTools import (  # type: ignore[attr-defined]
        TopTools_IndexedMapOfShape as ShapeMap,
    )
    from OCP.TopTools import (  # type: ignore[attr-defined]
        TopTools_ListOfShape as ShapeList,
    )
except ImportError:  # OCP 8.0
    from OCP.collections import (  # type: ignore[no-redef]
        IndexedDataMap_TopoDS_Shape_List_TopoDS_Shape_TopTools_ShapeMapHasher as EdgeFaceMap,
    )
    from OCP.collections import (  # type: ignore[no-redef]
        IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher as ShapeMap,
    )
    from OCP.collections import (  # type: ignore[no-redef]
        List_TopoDS_Shape as ShapeList,
    )

as_face = getattr(TopoDS, "Face_s", None) or TopoDS.Face
as_edge = getattr(TopoDS, "Edge_s", None) or TopoDS.Edge
as_vertex = getattr(TopoDS, "Vertex_s", None) or TopoDS.Vertex
as_wire = getattr(TopoDS, "Wire_s", None) or TopoDS.Wire
as_solid = getattr(TopoDS, "Solid_s", None) or TopoDS.Solid

OCP_VERSION: str = getattr(OCP, "__version__", "unknown")

__all__ = [
    "EdgeFaceMap",
    "ShapeMap",
    "ShapeList",
    "as_face",
    "as_edge",
    "as_vertex",
    "as_wire",
    "as_solid",
    "OCP_VERSION",
]
