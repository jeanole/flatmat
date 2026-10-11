"""Check that every OCP name the design relies on exists in this OCP build."""

import importlib
import sys

NAMES = {
    "OCP.STEPControl": ["STEPControl_Reader", "STEPControl_Writer", "STEPControl_AsIs"],
    "OCP.XSControl": ["XSControl_Reader"],
    "OCP.Interface": ["Interface_Static"],
    "OCP.TopExp": ["TopExp", "TopExp_Explorer"],
    "OCP.TopTools": ["TopTools_IndexedDataMapOfShapeListOfShape", "TopTools_IndexedMapOfShape"],
    "OCP.BRep": ["BRep_Tool"],
    "OCP.BRepTools": ["BRepTools"],
    "OCP.BRepAdaptor": ["BRepAdaptor_Surface", "BRepAdaptor_Curve"],
    "OCP.BRepLProp": ["BRepLProp_SLProps"],
    "OCP.BRepGProp": ["BRepGProp"],
    "OCP.GProp": ["GProp_GProps"],
    "OCP.IntCurvesFace": ["IntCurvesFace_ShapeIntersector"],
    "OCP.BRepTopAdaptor": ["BRepTopAdaptor_FClass2d"],
    "OCP.BRepClass3d": ["BRepClass3d_SolidClassifier"],
    "OCP.BRepExtrema": ["BRepExtrema_DistShapeShape"],
    "OCP.ShapeUpgrade": ["ShapeUpgrade_UnifySameDomain", "ShapeUpgrade_FaceDivide"],
    "OCP.ShapeAnalysis": ["ShapeAnalysis_CanonicalRecognition", "ShapeAnalysis_ShapeTolerance"],
    "OCP.GeomConvert": ["GeomConvert_SurfToAnaSurf"],
    "OCP.GeomLib": ["GeomLib_IsPlanarSurface"],
    "OCP.Geom2dConvert": ["Geom2dConvert"],
    "OCP.BRepCheck": ["BRepCheck_Analyzer"],
    "OCP.BRepLib": ["BRepLib"],
    "OCP.BRepBuilderAPI": ["BRepBuilderAPI_NurbsConvert", "BRepBuilderAPI_MakeFace"],
    "OCP.STEPCAFControl": ["STEPCAFControl_Reader"],
    "OCP.XCAFDoc": ["XCAFDoc_ShapeTool", "XCAFDoc_DocumentTool"],
}
STATICS = [
    ("OCP.TopExp", "TopExp", "MapShapesAndAncestors_s"),
    ("OCP.TopExp", "TopExp", "MapShapes_s"),
    ("OCP.BRep", "BRep_Tool", "CurveOnSurface_s"),
    ("OCP.BRep", "BRep_Tool", "IsClosed_s"),
    ("OCP.BRep", "BRep_Tool", "Parameter_s"),
    ("OCP.BRep", "BRep_Tool", "Tolerance_s"),
    ("OCP.BRepTools", "BRepTools", "UVBounds_s"),
    ("OCP.BRepGProp", "BRepGProp", "SurfaceProperties_s"),
    ("OCP.BRepGProp", "BRepGProp", "VolumeProperties_s"),
    ("OCP.Geom2dConvert", "Geom2dConvert", "CurveToBSplineCurve_s"),
    ("OCP.BRepLib", "BRepLib", "BuildPCurveForEdgeOnPlane_s"),
    ("OCP.BRepLib", "BRepLib", "OrientClosedSolid_s"),
]
METHODS = [
    ("OCP.XSControl", "XSControl_Reader", "SetShapeProcessFlags"),
    ("OCP.ShapeUpgrade", "ShapeUpgrade_UnifySameDomain", "KeepShapes"),
    ("OCP.ShapeUpgrade", "ShapeUpgrade_UnifySameDomain", "History"),
    ("OCP.GeomConvert", "GeomConvert_SurfToAnaSurf", "ConvertToAnalytical"),
    ("OCP.GeomConvert", "GeomConvert_SurfToAnaSurf", "Gap"),
]

missing = []
for mod, names in NAMES.items():
    try:
        m = importlib.import_module(mod)
    except ImportError as e:
        missing.append(f"{mod} (module): {e}")
        continue
    missing += [f"{mod}.{n}" for n in names if not hasattr(m, n)]
for mod, cls, attr in STATICS + METHODS:
    try:
        c = getattr(importlib.import_module(mod), cls)
    except (ImportError, AttributeError):
        continue  # already reported
    if not hasattr(c, attr):
        missing.append(f"{cls}.{attr}")

import OCP  # noqa: E402

ver = getattr(OCP, "__version__", "?")
print(f"python {sys.version.split()[0]}  OCP {ver}")
print("missing:" if missing else "all names present", *missing, sep="\n  ")
