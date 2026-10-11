"""The few places where OCP 7.8 and 8.0 differ, found by `api_probe.py`.

- 8.0 moved the NCollection instantiations to `OCP.collections` with long
  template names; 7.8 keeps the `TopTools_*` typedefs.
- 8.0 dropped the `_s` suffix on the TopoDS downcasts (`TopoDS.Face`), while
  7.8 has `TopoDS.Face_s`.
- `BRepTopAdaptor_FClass2d` is gone in 8.0; `BRepClass_FaceClassifier` exists in
  both and is what we use.
- `XSControl_Reader.SetShapeProcessFlags` (healing control) exists only in 8.0.
"""

from OCP.TopoDS import TopoDS

try:  # OCP 7.8
    from OCP.TopTools import (
        TopTools_IndexedDataMapOfShapeListOfShape as EdgeFaceMap,
    )
    from OCP.TopTools import (
        TopTools_IndexedMapOfShape as ShapeMap,
    )
    from OCP.TopTools import (
        TopTools_ListOfShape as ShapeList,
    )
except ImportError:  # OCP 8.0
    from OCP.collections import (
        IndexedDataMap_TopoDS_Shape_List_TopoDS_Shape_TopTools_ShapeMapHasher as EdgeFaceMap,
    )
    from OCP.collections import (
        IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher as ShapeMap,
    )
    from OCP.collections import (
        List_TopoDS_Shape as ShapeList,
    )

as_face = getattr(TopoDS, "Face_s", None) or TopoDS.Face
as_edge = getattr(TopoDS, "Edge_s", None) or TopoDS.Edge
as_vertex = getattr(TopoDS, "Vertex_s", None) or TopoDS.Vertex
as_solid = getattr(TopoDS, "Solid_s", None) or TopoDS.Solid
as_wire = getattr(TopoDS, "Wire_s", None) or TopoDS.Wire

__all__ = ["EdgeFaceMap", "ShapeMap", "ShapeList", "as_face", "as_edge", "as_vertex", "as_solid", "as_wire"]
