"""Run the spike: generated parts, BenDFM samples, a scaling series, DXF output.

python spike/run.py [--bendfm DIR] [--out DIR] [--scaling]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(__file__))

import ezdxf  # noqa: E402
import gen  # noqa: E402
import unfold as uf  # noqa: E402

K = 0.44


def check_generated(out_dir, results):
    print(f"\n== generated parts (K = {K}), STEP round trip")
    print(
        f"{'part':16} {'faces':>5} {'t':>5} {'bends':>5} {'area rel.err':>13} {'identity':>10} "
        f"{'r/angle/dir':>12} {'bspl err':>9} {'ms':>7}"
    )
    for part in gen.cases():
        solid = gen.make_solid(part)
        path = os.path.join(out_dir, f"{part.name}.step")
        gen.write_step(solid, path)
        t0 = time.perf_counter()
        solid = uf.read_step(path)
        read_ms = (time.perf_counter() - t0) * 1000
        try:
            res = uf.unfold(solid, k=K)
        except uf.Refused as e:
            print(f"{part.name:16} REFUSED {e}")
            results.append({"part": part.name, "refused": str(e)})
            continue
        # analytic expectation: width x flat length, minus developed hole areas (only
        # checkable for flange holes, which stay circles)
        expected = part.width * part.flat_length(K)
        hole_area = sum(math.pi * d * d / 4 for kind, *_, d in part.holes if kind == "flat")
        has_bend_hole = any(h[0] == "bend" for h in part.holes)
        area_err = abs(res.flat_area - (expected - hole_area)) / expected
        ident_err = abs(res.flat_area - res.expected_area) / res.expected_area
        # angle and radius must match; up/down is relative to the reference skin, so the
        # pattern of directions must match the signs of the generator up to a global flip
        got = sorted((round(b.angle_deg, 3), round(b.r_in, 3)) for b in res.bends)
        want = sorted((round(abs(a), 3), round(r, 3)) for a, r in part.bends)
        ups = sorted(b.direction for b in res.bends)
        signs = sorted("up" if a > 0 else "down" for a, _ in part.bends)
        flipped = sorted("down" if a > 0 else "up" for a, _ in part.bends)
        ok_bends = got == want and ups in (signs, flipped)
        thick_ok = abs(res.thickness - part.t) < 1e-3
        flag = "ok" if ok_bends and thick_ok else "MISMATCH"
        area_txt = "n/a" if has_bend_hole else f"{area_err:.1e}"
        total = read_ms + sum(res.timings.values())
        print(
            f"{part.name:16} {res.n_faces:5d} {res.thickness:5.2f} {len(res.bends):5d} {area_txt:>13} "
            f"{ident_err:10.1e} {flag:>12} {res.bspline_max_err:9.1e} {total:7.1f}"
        )
        write_dxf(res, os.path.join(out_dir, f"{part.name}.dxf"))
        results.append(
            {
                "part": part.name,
                "faces": res.n_faces,
                "thickness": res.thickness,
                "bends": [(b.angle_deg, b.r_in, b.direction) for b in res.bends],
                "area_rel_err": None if has_bend_hole else area_err,
                "identity_rel_err": ident_err,
                "bends_ok": ok_bends,
                "thickness_ok": thick_ok,
                "bspline_max_err": res.bspline_max_err,
                "read_ms": read_ms,
                "timings_ms": res.timings,
            }
        )


def skin_area_of_flat(path):
    """Area of the top skin of a BenDFM unfolded part: planar faces facing +Z at max z."""
    solid = uf.read_step(path)
    _, faces = uf.face_table(solid)
    zs = []
    for f in faces:
        if f.kind == "plane" and abs(abs(f.normal.Z()) - 1) < 1e-9:
            u0, u1, v0, v1 = f.uv
            from OCP.BRepAdaptor import BRepAdaptor_Surface

            p = BRepAdaptor_Surface(f.shape, True).Value((u0 + u1) / 2, (v0 + v1) / 2)
            zs.append((p.Z(), f))
    ztop = max(z for z, _ in zs)
    return sum(f.area for z, f in zs if abs(z - ztop) < 1e-6)


def check_bendfm(bdir, out_dir, results):
    print("\n== BenDFM samples (K = 0.45, their convention), area vs their unfolded STEP")
    print(
        f"{'id':>3} {'faces':>5} {'t':>6} {'t json':>6} {'bends':>5} {'json':>4} {'ours':>11} {'theirs':>11} "
        f"{'rel.err':>8} {'ms':>7}"
    )
    for i in range(6):
        folded = os.path.join(bdir, f"{i}.stp")
        if not os.path.exists(folded):
            continue
        seq = json.load(open(os.path.join(bdir, f"{i}_sequence.json")))
        t_json = seq[0]["part_thickness"]
        n_json = len(seq) - 1
        t0 = time.perf_counter()
        try:
            solid = uf.read_step(folded)
            read_ms = (time.perf_counter() - t0) * 1000
            res = uf.unfold(solid, k=0.45)
        except uf.Refused as e:
            print(f"{i:3d} REFUSED {e}")
            results.append({"bendfm": i, "refused": str(e)})
            continue
        theirs = skin_area_of_flat(os.path.join(bdir, f"{i}_unfolded.stp"))
        err = abs(res.flat_area - theirs) / theirs
        total = read_ms + sum(res.timings.values())
        print(
            f"{i:3d} {res.n_faces:5d} {res.thickness:6.3f} {t_json:6.3f} {len(res.bends):5d} {n_json:4d} "
            f"{res.flat_area:11.2f} {theirs:11.2f} {err:8.1e} {total:7.1f}"
        )
        write_dxf(res, os.path.join(out_dir, f"bendfm_{i}.dxf"))
        results.append(
            {
                "bendfm": i,
                "faces": res.n_faces,
                "thickness": res.thickness,
                "t_json": t_json,
                "bends": len(res.bends),
                "bends_json": n_json,
                "ours": res.flat_area,
                "theirs": theirs,
                "rel_err": err,
                "read_ms": read_ms,
                "timings_ms": res.timings,
            }
        )


def check_scaling(out_dir, results):
    print("\n== scaling: L bracket with N holes in the long flange")
    print(
        f"{'holes':>6} {'faces':>6} {'read':>8} "
        + " ".join(
            f"{k:>9}" for k in ("tables", "classify", "skins", "thickness", "bends", "place", "outline", "bendlines")
        )
        + f" {'total':>8}"
    )
    for n in (0, 10, 100, 400, 1000):
        cols = max(1, int(math.sqrt(n * 4)))
        rows = math.ceil(n / cols) if n else 0
        holes = []
        for j in range(n):
            r_, c_ = divmod(j, cols)
            holes.append(("flat", 0, 0.1 + 0.8 * (r_ + 0.5) / max(rows, 1), 5 + 390 * (c_ + 0.5) / cols, 1.5))
        part = gen.Part(f"perf_{n}", t=2.0, width=400, flats=[200, 40], bends=[(90, 2.0)], holes=holes)
        path = os.path.join(out_dir, f"{part.name}.step")
        gen.write_step(gen.make_solid(part), path)
        t0 = time.perf_counter()
        solid = uf.read_step(path)
        read_ms = (time.perf_counter() - t0) * 1000
        res = uf.unfold(solid, k=K)
        tm = res.timings
        total = read_ms + sum(tm.values())
        print(
            f"{n:6d} {res.n_faces:6d} {read_ms:8.1f} "
            + " ".join(
                f"{tm[k]:9.1f}"
                for k in ("tables", "classify", "skins", "thickness", "bends", "place", "outline", "bendlines")
            )
            + f" {total:8.1f}"
        )
        results.append(
            {
                "scaling_holes": n,
                "faces": res.n_faces,
                "read_ms": read_ms,
                "timings_ms": tm,
                "identity_rel_err": abs(res.flat_area - res.expected_area) / res.expected_area,
            }
        )


def write_dxf(res: uf.Result, path):
    doc = ezdxf.new("R2010", setup=True)
    doc.header["$INSUNITS"] = 4
    for name, color in (("OUTER", 7), ("INNER", 3), ("BEND_UP", 1), ("BEND_DOWN", 5), ("TEXT", 2)):
        doc.layers.add(name, color=color)
    msp = doc.modelspace()
    msp.add_lwpolyline(res.outer, close=True, dxfattribs={"layer": "OUTER"})
    for h in res.holes:
        msp.add_lwpolyline(h, close=True, dxfattribs={"layer": "INNER"})
    for b in res.bends:
        layer = "BEND_UP" if b.direction == "up" else "BEND_DOWN"
        for p0, p1 in b.segments:
            msp.add_line(p0, p1, dxfattribs={"layer": layer, "linetype": "CENTER"})
        if b.segments:
            (x0, y0), (x1, y1) = b.segments[0]
            sign = "+" if b.direction == "up" else "-"
            msp.add_text(
                f"{sign}{b.angle_deg:.0f} R{b.r_in:g}", dxfattribs={"layer": "TEXT", "height": 2.5}
            ).set_placement(((x0 + x1) / 2, (y0 + y1) / 2))
    ezdxf.options.write_fixed_meta_data_for_testing = True
    doc.saveas(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bendfm", help="directory with BenDFM example files 0..5")
    ap.add_argument("--out", default=tempfile.mkdtemp(prefix="flatmat_spike_"))
    ap.add_argument("--scaling", action="store_true")
    ap.add_argument("--json", help="write all results to this file")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    import OCP

    print(f"python {platform.python_version()}  OCP {getattr(OCP, '__version__', '?')}  out={a.out}")
    results = []
    check_generated(a.out, results)
    if a.bendfm:
        check_bendfm(a.bendfm, a.out, results)
    if a.scaling:
        check_scaling(a.out, results)
    if a.json:
        json.dump(
            {"python": platform.python_version(), "ocp": getattr(OCP, "__version__", "?"), "results": results},
            open(a.json, "w"),
            indent=1,
            default=str,
        )


if __name__ == "__main__":
    main()
