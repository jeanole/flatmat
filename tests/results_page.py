"""The test results page: ``pytest tests --results-page docs/test-results.md``.

At the end of the session it writes a markdown page with the outcome of each test file
and one row per committed sheet metal test part (``tests/data/*/<name>.step`` with its
hand-worked ``<name>.json``): the part unfolded at the JSON's K against the flat in the
JSON, and how many of the session's tests that name the part passed. The page is a
report: the tests themselves decide pass or fail.
"""

from __future__ import annotations

import datetime as dt
import json
import os
import platform
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field

import pytest

import flatmat as fm

HERE = os.path.dirname(__file__)
FOLDERS = [
    ("solidcraft_oracle", "SolidCraft reference parts", "test_oracle.py"),
    ("handworked", "Hand-worked parts: holes across bends, boxes, hems, joggles, countersinks", "test_handworked.py"),
]
TOL_AREA, TOL_SIZE = 1e-9, 1e-6  # the tests' own tolerances: relative area, size in mm


@dataclass
class Outcomes:
    """What the session's tests reported: per test file, and per part name in a test's id."""

    files: dict[str, dict[str, int]] = field(default_factory=lambda: defaultdict(lambda: defaultdict(int)))
    parts: dict[str, dict[str, int]] = field(default_factory=lambda: defaultdict(lambda: defaultdict(int)))
    failed: list[str] = field(default_factory=list)

    def add(self, report: pytest.TestReport, names: set[str]) -> None:
        if report.when != "call" and not (report.failed or report.skipped):
            return  # setup and teardown count only when they fail or skip the test
        outcome = report.outcome
        path = report.nodeid.split("::")[0]
        self.files[os.path.basename(path)][outcome] += 1
        params = re.search(r"\[(.*)\]$", report.nodeid)
        if params:
            for name in set(params.group(1).split("-")) & names:
                self.parts[name][outcome] += 1
        if outcome == "failed":
            self.failed.append(report.nodeid)


def _names(folder: str) -> list[str]:
    return sorted(f[:-5] for f in os.listdir(os.path.join(HERE, "data", folder)) if f.endswith(".json"))


def part_names() -> set[str]:
    return {n for folder, _, _ in FOLDERS for n in _names(folder)}


def _descriptions(folder: str) -> dict[str, str]:
    """``| `name` | description | ...`` rows of the folder's README."""
    out = {}
    with open(os.path.join(HERE, "data", folder, "README.md")) as fh:
        for line in fh:
            m = re.match(r"\|\s*`([^`]+)`\s*\|\s*([^|]+)\|", line)
            if m:
                out[m.group(1)] = m.group(2).strip()
    return out


def _skin_area(sheet, want: dict[str, float]) -> float:
    """A flat whose area depends on the skin developed (a radial hole): the developed skin's,
    told by its bend face's area as test_handworked does (inner skin: the smaller face)."""
    (cyl,) = (f.area for f in sheet.faces if f.role == "skin_a" and f.kind == "cylinder")
    others = [f.area for f in sheet.faces if f.role == "skin_b" and f.kind == "cylinder"]
    inner = not others or cyl < others[0]
    return want["inner" if inner else "outer"]


def _row(folder: str, name: str) -> dict:
    with open(os.path.join(HERE, "data", folder, f"{name}.json")) as fh:
        want = json.load(fh)
    flat = want["flat"]
    try:
        sheet = fm.recognize(fm.load_step(os.path.join(HERE, "data", folder, f"{name}.step")))
        pat = fm.unfold(sheet, k=want["k"])
    except fm.FlatmatError as e:
        return {"name": name, "error": f"{e.diagnostic.code.value}: {e.diagnostic.message}", "want": want}
    area = flat.get("area_mm2")
    if area is None:
        area = _skin_area(sheet, flat["area_mm2_by_developed_skin"])
    size = sorted(pat.size)
    return {
        "name": name,
        "want": want,
        "t": pat.thickness,
        "bends": len(pat.bends),
        "area": (area, pat.area, abs(pat.area - area) / area),
        "size": (flat["size"], size, max(abs(a - b) for a, b in zip(size, flat["size"], strict=True))),
        "holes": (flat["holes"], len(pat.holes)),
    }


def _git(*args: str) -> str:
    try:
        return subprocess.run(["git", *args], cwd=HERE, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def _status(row: dict, tests: dict[str, int]) -> str:
    bad = tests.get("failed", 0) or "error" in row
    if not bad:
        bad = row["area"][2] > TOL_AREA or row["size"][2] > TOL_SIZE or row["holes"][0] != row["holes"][1]
    return "❌" if bad else "✅"


def render(outcomes: Outcomes) -> str:
    sha = os.environ.get("GITHUB_SHA") or _git("rev-parse", "HEAD")
    now = dt.datetime.now(dt.UTC).strftime("%Y-%m-%d %H:%M UTC")
    totals: dict[str, int] = defaultdict(int)
    for counts in outcomes.files.values():
        for k, v in counts.items():
            totals[k] += v
    lines = [
        "# Test results",
        "",
        "Written by `pytest tests --results-page docs/test-results.md` (in CI on every push to `main`).",
        "Do not edit by hand.",
        "",
        f"- **Run:** {now}, commit `{sha[:10] or 'unknown'}`",
        f"- **Environment:** Python {platform.python_version()}, OCP {fm.io.ocp_version()}, "
        f"{platform.system()} {platform.machine()}",
        f"- **Tests:** {totals['passed']} passed, {totals['failed']} failed, {totals['skipped']} skipped",
        "",
        "## Test files",
        "",
        "| File | Passed | Failed | Skipped |",
        "|---|---:|---:|---:|",
    ]
    for f in sorted(outcomes.files):
        c = outcomes.files[f]
        lines.append(f"| `{f}` | {c['passed']} | {c['failed']} | {c['skipped']} |")
    if outcomes.failed:
        lines += ["", "**Failed:**", ""] + [f"- `{n}`" for n in outcomes.failed[:50]]
    for folder, title, test_file in FOLDERS:
        desc = _descriptions(folder)
        names = _names(folder)
        lines += [
            "",
            f"## {title}",
            "",
            f"STEP files and hand-worked flats in [`tests/data/{folder}`](../tests/data/{folder}/README.md), "
            f"checked by `tests/{test_file}`. Each part unfolded at its K: flat area and size against the "
            f"hand-worked flat (tolerance {TOL_AREA:g} relative, {TOL_SIZE:g} mm). *Tests*: the session's tests "
            "that name the part, passed / run.",
            "",
            "| | Part | What it is | t | Bends | K | Flat area: hand-worked | flatmat | Rel. error | "
            "Flat size (mm) | Size error (mm) | Holes | Tests |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|",
        ]
        for name in names:
            row = _row(folder, name)
            tests = outcomes.parts.get(name, {})
            run = sum(tests.values())
            t_col = f"{tests.get('passed', 0)} / {run}" if run else "–"
            what = desc.get(name, "").replace("|", "/")
            if "error" in row:
                k = row["want"]["k"]
                lines.append(f"| ❌ | `{name}` | {what} | | | {k:g} | | {row['error']} | | | | | {t_col} |")
                continue
            (wa, ga, ea), (_, gs, es), (wh, gh) = row["area"], row["size"], row["holes"]
            lines.append(
                f"| {_status(row, tests)} | `{name}` | {what} | {row['t']:g} | {row['bends']} | {row['want']['k']:g} | "
                f"{wa:.4f} | {ga:.4f} | {ea:.1e} | {gs[0]:.3f} × {gs[1]:.3f} | {es:.1e} | "
                f"{gh}" + ("" if gh == wh else f" (want {wh})") + f" | {t_col} |"
            )
    return "\n".join(lines) + "\n"


def write(path: str, outcomes: Outcomes) -> None:
    text = render(outcomes)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")  # the Actions run page shows it too
    if summary:
        with open(summary, "a", encoding="utf-8") as fh:
            fh.write(text)
