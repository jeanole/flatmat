"""Pictures of the committed sheet metal test parts for the results page.

    pip install matplotlib
    python tests/render_parts.py            # writes docs/parts/<name>.png

One PNG per ``tests/data/*/<name>.step``: the part in 3D (shaded) beside its flat
pattern from flatmat at the K of its JSON, with the bend lines (up red, down blue,
joggle bends orange) and the holes. Rerun it when a part changes; the results page
(``tests/results_page.py``) links to the pictures that exist.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402
from OCP.BRep import BRep_Tool  # noqa: E402
from OCP.BRepMesh import BRepMesh_IncrementalMesh  # noqa: E402
from OCP.TopAbs import TopAbs_FACE  # noqa: E402
from OCP.TopExp import TopExp_Explorer  # noqa: E402
from OCP.TopLoc import TopLoc_Location  # noqa: E402

import flatmat as fm  # noqa: E402
from flatmat._occ import as_face  # noqa: E402

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "..", "docs", "parts")
FOLDERS = ("solidcraft_oracle", "handworked")
COLOURS = {"up": "#c0392b", "down": "#1f77b4", "joggle": "#e67e22"}


def _triangles(shape) -> np.ndarray:
    BRepMesh_IncrementalMesh(shape, 0.05, False, 0.25, True)
    tris = []
    ex = TopExp_Explorer(shape, TopAbs_FACE)
    while ex.More():
        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation_s(as_face(ex.Current()), loc)
        if tri is not None:
            tr = loc.Transformation()
            pts = [tri.Node(i).Transformed(tr) for i in range(1, tri.NbNodes() + 1)]
            for j in range(1, tri.NbTriangles() + 1):
                tris.append([(pts[i - 1].X(), pts[i - 1].Y(), pts[i - 1].Z()) for i in tri.Triangle(j).Get()])
        ex.Next()
    return np.array(tris)


def _draw_3d(ax, shape) -> None:
    tris = _triangles(shape)
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= np.linalg.norm(n, axis=1)[:, None] + 1e-12
    light = np.array([0.35, -0.55, 0.75])
    light /= np.linalg.norm(light)
    shade = 0.35 + 0.65 * np.abs(n @ light)
    ax.add_collection3d(Poly3DCollection(tris, facecolors=np.outer(shade, [0.42, 0.58, 0.80]), edgecolor="none"))
    lo, hi = tris.reshape(-1, 3).min(0), tris.reshape(-1, 3).max(0)
    mid, half = (lo + hi) / 2, (hi - lo).max() / 2
    ax.set_xlim(mid[0] - half, mid[0] + half)
    ax.set_ylim(mid[1] - half, mid[1] + half)
    ax.set_zlim(mid[2] - half, mid[2] + half)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=24, azim=-58)
    ax.set_axis_off()


def _draw_flat(ax, pat: fm.FlatPattern) -> None:
    outer = pat.outer.polyline()
    ax.fill([p[0] for p in outer], [p[1] for p in outer], color="#dbe6f3", lw=0)
    for loop in (pat.outer, *pat.holes):
        pts = loop.polyline()
        pts.append(pts[0])
        if loop is not pat.outer:
            ax.fill([p[0] for p in pts], [p[1] for p in pts], color="white", lw=0)
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color="#1d3557", lw=0.9)
    kinds = {b.id: b.kind for b in pat.bends}
    for bl in pat.bend_lines:
        colour = COLOURS["joggle"] if kinds.get(bl.bend_id) == "joggle" else COLOURS[bl.direction]
        ax.plot([bl.p0[0], bl.p1[0]], [bl.p0[1], bl.p1[1]], ls="--", color=colour, lw=0.9)
    w, h = pat.size
    ax.set_title(f"flat {w:.1f} × {h:.1f} mm", fontsize=8, pad=2)
    ax.set_aspect("equal")
    ax.set_axis_off()


def render(folder: str, name: str, out_dir: str = OUT) -> str:
    base = os.path.join(HERE, "data", folder, name)
    with open(base + ".json") as fh:
        k = json.load(fh)["k"]
    model = fm.load_step(base + ".step", pmi=False)
    pat = fm.unfold(fm.recognize(model), k=k)
    fig = plt.figure(figsize=(5.2, 2.2))
    _draw_3d(fig.add_axes((0.0, 0.0, 0.48, 1.0), projection="3d"), model.solid)
    _draw_flat(fig.add_axes((0.52, 0.04, 0.46, 0.84)), pat)
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{name}.png")
    fig.savefig(path, dpi=100, metadata={"Software": None})  # no timestamp: same part, same bytes
    plt.close(fig)
    return path


if __name__ == "__main__":
    only = set(sys.argv[1:])
    for folder in FOLDERS:
        for f in sorted(os.listdir(os.path.join(HERE, "data", folder))):
            if f.endswith(".json") and (not only or f[:-5] in only):
                print(os.path.relpath(render(folder, f[:-5])))
