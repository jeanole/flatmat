"""DXF, bend table (JSON/CSV) and SVG output (docs/design.md §2)."""

from __future__ import annotations

import csv
import json
import math
import os
from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

import ezdxf

from . import __version__
from ._occ import OCP_VERSION
from .geom2d import heading, midpoint
from .model import Arc, Feature, FlatPattern, Line, Loop2D, Spline

Labels = Literal["delem", "plain"] | None


@dataclass(frozen=True)
class LayerMap:
    outer: str = "OUTER"
    inner: str = "INNER"
    bend_up: str = "BEND_UP"
    bend_down: str = "BEND_DOWN"
    text: str = "TEXT"
    form: str = "FORM"  # rims of embosses and extruded holes: markings, not cuts
    engrave: str = "ENGRAVE"  # rims of engravings: marked by the laser, not cut

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
    dfm=None,
) -> None:
    """Write the flat pattern as DXF with exact LINE/ARC/CIRCLE/SPLINE entities."""
    layers = layers or LayerMap()
    if deterministic:  # must be set before ezdxf.new(): it stamps the creation time
        ezdxf.options.write_fixed_meta_data_for_testing = True
    try:
        _write_dxf(pattern, path, layers, labels, version, text_height, dfm)
    finally:
        if deterministic:
            ezdxf.options.write_fixed_meta_data_for_testing = False


def _write_dxf(pattern, path, layers, labels, version, text_height, dfm=None) -> None:
    doc = ezdxf.new(version, setup=True)
    doc.header["$INSUNITS"] = 4  # mm
    for name, color in (
        (layers.outer, 7),
        (layers.inner, 3),
        (layers.bend_up, 1),
        (layers.bend_down, 5),
        (layers.text, 2),
        (layers.form, 4),
        (layers.engrave, 6),
    ):
        if name not in doc.layers:
            doc.layers.add(name, color=color)
    msp = doc.modelspace()
    _add_loop(msp, pattern.outer, layers.outer)
    for h in pattern.holes:
        _add_loop(msp, h, layers.inner)
    for sl in pattern.slits:  # cut lines without width
        _add_loop(msp, sl, layers.inner)
    for f in pattern.features:
        if f.outline is not None:
            _add_loop(msp, f.outline, layers.engrave if f.kind_info.rim_layer == "engrave" else layers.form)
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
            mid = midpoint(bl.p0, bl.p1)
            ang = math.degrees(heading(bl.p0, bl.p1))
            if ang > 90 or ang <= -90:
                ang -= 180 if ang > 0 else -180
            msp.add_text(
                _step_prefix(dfm, b.id) + _label(b, labels),
                dxfattribs={"layer": layers.text, "height": th, "rotation": ang},
            ).set_placement(mid)
        for f in pattern.features:
            if f.center is None or not f.kind_info.label:  # an engraving's rim says it all (text is many rims)
                continue
            x, y = f.center
            msp.add_text(_feature_label(f), dxfattribs={"layer": layers.text, "height": th * 0.6}).set_placement(
                (x + max(f.through_diameter, f.outer_diameter) / 2 + th * 0.3, y)
            )
    if dfm is not None and dfm.findings:
        # a marker per finding on its own layer, so the shop can switch it off
        if "DFM" not in doc.layers:
            doc.layers.add("DFM", color=6)
        for d in dfm.findings:
            at = d.data.get("at")
            if at is None:
                continue
            msp.add_circle(at, th * 1.5, dxfattribs={"layer": "DFM"})
            msp.add_text(
                d.code.value.removeprefix("DFM_"), dxfattribs={"layer": "DFM", "height": th * 0.6}
            ).set_placement((at[0] + th * 1.8, at[1]))
    doc.saveas(os.fspath(path))


def _steps(dfm) -> dict[int, int]:
    """Bend id -> its step in the report's bend order, if the order is feasible."""
    plan = getattr(dfm, "sequence", None)
    if plan is None or not plan.feasible:
        return {}
    return {b: i + 1 for i, b in enumerate(plan.order)}


def _step_prefix(dfm, bend_id: int) -> str:
    n = _steps(dfm).get(bend_id)
    return f"({n}) " if n is not None else ""


def _feature_label(f: Feature) -> str:
    return _LABELS[f.kind](f, f.side.upper())


def _counterbore_label(f: Feature, side: str) -> str:
    depth = f" T{f.depth:g}" if f.depth else ""
    return f"CBORE Ø{f.outer_diameter:g}{depth} {side}"


def _oriented_label(f: Feature, side: str) -> str:
    size = f"{f.size[0]:.4g}x{f.size[1]:.4g} " if f.size else ""
    angle = f"A{f.orientation_deg:.4g} " if f.orientation_deg is not None else ""
    return f"{f.kind.upper()} {size}{angle}H{f.height:.4g} {side}"


# the DXF text next to each kind of feature (one entry per kind in model.FEATURE_KINDS)
_LABELS: dict[str, Callable[[Feature, str], str]] = {
    "countersink": lambda f, side: f"CSK {f.angle_deg:.0f}° Ø{f.outer_diameter:g} {side}",
    "counterbore": _counterbore_label,
    "emboss": lambda f, side: f"EMBOSS H{f.height:.4g} {side}",
    "extruded_hole": lambda f, side: (
        f"EXTR Ø{f.outer_diameter:.4g} H{f.height:.4g} {side} (PRE Ø{f.through_diameter:.3f})"
    ),
    "louvre": _oriented_label,
    "bridge": _oriented_label,
    "lance": _oriented_label,
    "engraving": lambda f, side: f"ENGRAVE T{f.depth:.3g} {side}",
}


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
    features = [
        {
            "id": f.id,
            "kind": f.kind,
            "side": f.side,
            "center": list(f.center) if f.center else None,
            "through_diameter": f.through_diameter,
            "outer_diameter": f.outer_diameter,
            "angle_deg": f.angle_deg,
            "depth": f.depth,
            "height": f.height,
            "radius": f.radius,
            "rim": [list(p) for p in f.outline.polyline(0.01)] if f.outline is not None else None,
            "orientation_deg": f.orientation_deg,
            "size": list(f.size) if f.size else None,
            "holes": list(f.holes),
            "bend": f.bend_id,
        }
        for f in pattern.features
    ]
    return {
        "schema": "flatmat.bendtable/1",
        "units": "mm",
        "thickness": pattern.thickness,
        "k_model": pattern.k_model,
        "flat_size": [w, h],
        "flat_area": pattern.area,
        "bends": bends,
        "features": features,
        "slits": [[list(p) for p in sl.polyline(0.01)] for sl in pattern.slits],
        "generator": {"flatmat": __version__, "ocp": OCP_VERSION},
        "diagnostics": [d.to_dict() for d in pattern.diagnostics],
    }


def write_bend_table(pattern: FlatPattern, path: str | os.PathLike, *, dfm=None) -> None:
    """JSON (``.json``) or CSV (``.csv``, one row per bend line segment).

    With a :class:`~flatmat.dfm.DfmReport` the JSON gets a ``dfm`` section, and if the
    report has a feasible bend order, every bend its ``step`` (CSV: a ``step`` column).
    """
    path = os.fspath(path)
    table = bend_table(pattern)
    steps = _steps(dfm)
    if dfm is not None:
        table["dfm"] = dfm.to_dict()
        for b in table["bends"]:
            b["step"] = steps.get(b["id"])
    if path.lower().endswith(".csv"):
        cols = ["bend", "segment", "direction", "angle_deg", "r_in", "k", "ba", "bd", "x0", "y0", "x1", "y1", "length"]
        cols += ["step"] if steps else []
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
                        ]
                        + ([steps.get(b["id"])] if steps else [])  # fmt: skip
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
    for f in pattern.features:
        if f.outline is not None:
            out.append(f'<polygon points="{pts(f.outline)}" fill="none" stroke="#a6a" stroke-width="{sw:.4f}"/>')
    for sl in pattern.slits:
        out.append(f'<polyline points="{pts(sl)}" fill="none" stroke="#2a7" stroke-width="{sw:.4f}"/>')
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
