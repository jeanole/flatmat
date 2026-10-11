"""DXF, bend table (JSON/CSV) and SVG output (docs/design.md §2)."""

from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass
from typing import Literal

import ezdxf

from . import __version__
from ._occ import OCP_VERSION
from .model import Arc, FlatPattern, Line, Loop2D, Spline

Labels = Literal["delem", "plain"] | None


@dataclass(frozen=True)
class LayerMap:
    outer: str = "OUTER"
    inner: str = "INNER"
    bend_up: str = "BEND_UP"
    bend_down: str = "BEND_DOWN"
    text: str = "TEXT"

    @classmethod
    def default(cls) -> LayerMap:
        return cls()


def _label(b, style: Labels) -> str:
    if style == "delem":
        sign = "+" if b.direction == "up" else "-"
        head = f"H{sign}{b.angle_deg:.0f}" if b.kind == "hem" else f"{sign}{b.angle_deg:.0f}"
        return f"{head} R{b.r_in:g}"
    return f"{b.direction.upper()} {b.angle_deg:.1f}° R{b.r_in:g}"


def write_dxf(
    pattern: FlatPattern,
    path: str | os.PathLike,
    *,
    layers: LayerMap | None = None,
    labels: Labels = "delem",
    version: str = "R2010",
    text_height: float | None = None,
    deterministic: bool = False,
) -> None:
    """Write the flat pattern as DXF with exact LINE/ARC/CIRCLE/SPLINE entities."""
    layers = layers or LayerMap()
    if deterministic:  # must be set before ezdxf.new(): it stamps the creation time
        ezdxf.options.write_fixed_meta_data_for_testing = True
    try:
        _write_dxf(pattern, path, layers, labels, version, text_height)
    finally:
        if deterministic:
            ezdxf.options.write_fixed_meta_data_for_testing = False


def _write_dxf(pattern, path, layers, labels, version, text_height) -> None:
    doc = ezdxf.new(version, setup=True)
    doc.header["$INSUNITS"] = 4  # mm
    for name, color in (
        (layers.outer, 7),
        (layers.inner, 3),
        (layers.bend_up, 1),
        (layers.bend_down, 5),
        (layers.text, 2),
    ):
        if name not in doc.layers:
            doc.layers.add(name, color=color)
    msp = doc.modelspace()
    _add_loop(msp, pattern.outer, layers.outer)
    for h in pattern.holes:
        _add_loop(msp, h, layers.inner)
    w, h = pattern.size
    th = text_height or max(1.0, min(5.0, 0.02 * max(w, h)))
    for bl in pattern.bend_lines:
        layer = layers.bend_up if bl.direction == "up" else layers.bend_down
        msp.add_line(bl.p0, bl.p1, dxfattribs={"layer": layer, "linetype": "CENTER"})
    if labels:
        done = set()
        for bl in pattern.bend_lines:
            if bl.bend_id in done:
                continue
            done.add(bl.bend_id)
            b = next(b for b in pattern.bends if b.id == bl.bend_id)
            mid = ((bl.p0[0] + bl.p1[0]) / 2, (bl.p0[1] + bl.p1[1]) / 2)
            ang = math.degrees(math.atan2(bl.p1[1] - bl.p0[1], bl.p1[0] - bl.p0[0]))
            if ang > 90 or ang <= -90:
                ang -= 180 if ang > 0 else -180
            msp.add_text(
                _label(b, labels), dxfattribs={"layer": layers.text, "height": th, "rotation": ang}
            ).set_placement(mid)
    doc.saveas(os.fspath(path))


def _add_loop(msp, loop: Loop2D, layer: str) -> None:
    att = {"layer": layer}
    for s in loop.segments:
        if isinstance(s, Line):
            msp.add_line(s.p0, s.p1, dxfattribs=att)
        elif isinstance(s, Arc):
            if s.is_circle:
                msp.add_circle(s.center, s.radius, dxfattribs=att)
            else:
                a0, a1 = (s.a0, s.a1) if s.ccw else (s.a1, s.a0)  # DXF arcs run counter-clockwise
                msp.add_arc(s.center, s.radius, math.degrees(a0), math.degrees(a1), dxfattribs=att)
        elif isinstance(s, Spline):
            knots = s.flat_knots
            if s.weights:
                msp.add_rational_spline(s.poles, s.weights, s.degree, knots, dxfattribs=att)
            else:
                msp.add_open_spline(s.poles, s.degree, knots, dxfattribs=att)


# ------------------------------------------------------------------------------ bend table


def bend_table(pattern: FlatPattern) -> dict:
    by_bend: dict[int, list] = {}
    for bl in pattern.bend_lines:
        by_bend.setdefault(bl.bend_id, []).append(bl)
    bends = []
    for b in pattern.bends:
        segs = [
            {"line": [list(bl.p0), list(bl.p1)], "length": math.dist(bl.p0, bl.p1), "face": bl.segment_face}
            for bl in by_bend.get(b.id, [])
        ]
        bends.append(
            {
                "id": b.id,
                "direction": b.direction,
                "angle_deg": b.angle_deg,
                "opening_angle_deg": b.opening_angle_deg,
                "r_in": b.r_in,
                "k": b.k,
                "ba": b.ba,
                "bd": b.bd,
                "ossb": b.ossb,
                "kind": b.kind,
                "segments": segs,
            }
        )
    w, h = pattern.size
    return {
        "schema": "flatmat.bendtable/1",
        "units": "mm",
        "thickness": pattern.thickness,
        "k_model": pattern.k_model,
        "flat_size": [w, h],
        "flat_area": pattern.area,
        "bends": bends,
        "generator": {"flatmat": __version__, "ocp": OCP_VERSION},
        "diagnostics": [d.to_dict() for d in pattern.diagnostics],
    }


def write_bend_table(pattern: FlatPattern, path: str | os.PathLike) -> None:
    """JSON (``.json``) or CSV (``.csv``, one row per bend line segment)."""
    path = os.fspath(path)
    table = bend_table(pattern)
    if path.lower().endswith(".csv"):
        cols = ["bend", "segment", "direction", "angle_deg", "r_in", "k", "ba", "bd", "x0", "y0", "x1", "y1", "length"]
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for b in table["bends"]:
                for n, s in enumerate(b["segments"]):
                    (x0, y0), (x1, y1) = s["line"]
                    w.writerow(
                        [
                            b["id"],
                            n + 1,
                            b["direction"],
                            b["angle_deg"],
                            b["r_in"],
                            b["k"],
                            b["ba"],
                            b["bd"],
                            x0,
                            y0,
                            x1,
                            y1,
                            s["length"],
                        ]  # fmt: skip
                    )
    else:
        with open(path, "w") as fh:
            json.dump(table, fh, indent=1)


# ------------------------------------------------------------------------------ SVG


def write_svg(pattern: FlatPattern, path: str | os.PathLike, *, margin: float = 5.0) -> None:
    """A quick visual check (not for production): outline, holes, bend lines."""
    w, h = pattern.size
    W, H = w + 2 * margin, h + 2 * margin

    def pts(loop: Loop2D) -> str:
        return " ".join(f"{x + margin:.4f},{H - (y + margin):.4f}" for x, y in loop.polyline(0.02))

    sw = max(w, h) / 800
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.3f} {H:.3f}" '
        f'width="{min(1200, W * 3):.0f}" height="{min(1200, W * 3) * H / W:.0f}">',
        f'<rect width="{W:.3f}" height="{H:.3f}" fill="white"/>',
        f'<polygon points="{pts(pattern.outer)}" fill="#f2f2f2" stroke="black" stroke-width="{sw:.4f}"/>',
    ]
    for hole in pattern.holes:
        out.append(f'<polygon points="{pts(hole)}" fill="white" stroke="#2a7" stroke-width="{sw:.4f}"/>')
    for bl in pattern.bend_lines:
        col = "#d22" if bl.direction == "up" else "#22d"
        (x0, y0), (x1, y1) = bl.p0, bl.p1
        out.append(
            f'<line x1="{x0 + margin:.4f}" y1="{H - (y0 + margin):.4f}" x2="{x1 + margin:.4f}" '
            f'y2="{H - (y1 + margin):.4f}" stroke="{col}" stroke-width="{sw:.4f}" '
            f'stroke-dasharray="{6 * sw:.4f} {3 * sw:.4f}"/>'
        )
    out.append("</svg>")
    with open(os.fspath(path), "w") as fh:
        fh.write("\n".join(out))
