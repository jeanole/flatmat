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


def _readme(folder: str) -> dict[str, list[str]]:
    """``| `name` | cell | cell |`` rows of the folder's README: what the part is, and the flat by hand."""
    out = {}
    with open(os.path.join(HERE, "data", folder, "README.md"), encoding="utf-8") as fh:
        for line in fh:
            m = re.match(r"\|\s*`([^`]+)`\s*\|(.*)\|\s*$", line)
            if m:
                out[m.group(1)] = [c.strip() for c in m.group(2).split("|")]
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


def _passed(row: dict, tests: dict[str, int]) -> bool:
    if tests.get("failed", 0) or "error" in row:
        return False
    return row["area"][2] <= TOL_AREA and row["size"][2] <= TOL_SIZE and row["holes"][0] == row["holes"][1]


def _source(folder: str, want: dict) -> str:
    m = re.search(r"tests/oracle/(\S+)", want.get("source", ""))
    if folder == "solidcraft_oracle" and m:
        case = m.group(1)
        url = f"https://github.com/tangency-labs/solidcraft/tree/main/tests/oracle/{case}"
        return f"[SolidCraft {case.split('-')[0]}]({url}) (analytic); flat by hand"
    return "by hand"


INTRO = f"""Each sheet metal test part is a STEP file with its flat pattern worked out by hand
(`<name>.json` beside it). flatmat unfolds the part at the K the JSON gives and is compared
with that flat: flat area (relative tolerance {TOL_AREA:g}), flat size ({TOL_SIZE:g} mm) and the
number of holes. The tests that name a part check more: thickness, bends, each bend's
allowance θ·(R + K·t), hole positions, features, other K-factors and, where it holds,
flat area × t = volume at K 0.5.

Where the expected values come from (the **Source** column):

- **SolidCraft (analytic)**: a part of SolidCraft's reference tests
  ([`tests/oracle/`](https://github.com/tangency-labs/solidcraft/tree/main/tests/oracle)),
  rebuilt here. Its folded volume and area match SolidCraft's hand-worked values to 1e-9, so it
  is the same part; its flat is worked out by hand here (SolidCraft checks only the folded part).
- **by hand**: closed-form geometry, the formula in the **Method** column (BA = θ·(R + K·t)).

No expected value comes from a commercial CAD product. The pictures show each part in 3D and
flatmat's flat pattern (bend lines: up red, down blue, joggle orange); they are drawn by
`python tests/render_parts.py`, not by the test run."""


def render(outcomes: Outcomes) -> str:
    sha = os.environ.get("GITHUB_SHA") or _git("rev-parse", "HEAD")
    now = dt.datetime.now(dt.UTC).strftime("%Y-%m-%d %H:%M UTC")
    totals: dict[str, int] = defaultdict(int)
    for counts in outcomes.files.values():
        for k, v in counts.items():
            totals[k] += v
    sections, n_pass, n_parts, sources = [], 0, 0, []
    for folder, title, test_file in FOLDERS:
        readme = _readme(folder)
        rows = []
        for name in _names(folder):
            row = _row(folder, name)
            tests = outcomes.parts.get(name, {})
            run = sum(tests.values())
            ok = _passed(row, tests)
            n_parts += 1
            n_pass += ok
            cells = readme.get(name, [])
            what = cells[0] if cells else ""
            method = cells[1] if len(cells) > 1 else ""
            pic = f'<img src="parts/{name}.png" width="360" alt="{name}">' if os.path.exists(
                os.path.join(HERE, "..", "docs", "parts", f"{name}.png")) else ""  # fmt: skip
            part = f"{pic}<br>`{name}`" if pic else f"`{name}`"
            result = "✅ pass" if ok else "❌ **fail**"
            test_cell = f"{tests.get('passed', 0)} / {run}" if run else "–"
            if "error" in row:
                rows.append(f"| {part} | {result} | {_source(folder, row['want'])} | {row['error']} | | {test_cell} | "
                            f"{method} | {what} |")  # fmt: skip
                continue
            (_, ga, ea), (_, gs, es), (wh, gh) = row["area"], row["size"], row["holes"]
            flat = f"{gs[0]:.3f} × {gs[1]:.3f} mm<br>area {ga:.3f} mm²<br>t {row['t']:g}, K {row['want']['k']:g}, "
            flat += f"{row['bends']} bend{'' if row['bends'] == 1 else 's'}, {gh} hole{'' if gh == 1 else 's'}" + (
                "" if gh == wh else f" (hand: {wh})"
            )
            err = f"area {ea:.1e}<br>size {es:.1e} mm"
            rows.append(f"| {part} | {result} | {_source(folder, row['want'])} | {flat} | {err} | {test_cell} | "
                        f"{method} | {what} |")  # fmt: skip
        sources.append(f"{len(rows)} {title.split(':')[0].lower()}")
        sections += [
            "",
            f"## {title}",
            "",
            f"Files: [`tests/data/{folder}`](../tests/data/{folder}/README.md); tests: `tests/{test_file}`.",
            "",
            "| Part | Result | Source | Flat (flatmat) | Error vs hand | Tests | Method: flat by hand | What it is |",
            "|---|---|---|---|---|---:|---|---|",
            *rows,
        ]
    lines = [
        "# Test results",
        "",
        "Generated by `pytest tests --results-page docs/test-results.md` (in CI on every push to `main`).",
        "Do not edit by hand.",
        "",
        INTRO,
        "",
        f"**{n_pass} / {n_parts} parts pass**; sources: {', '.join(sources)}. "
        f"All tests: {totals['passed']} passed, {totals['failed']} failed, {totals['skipped']} skipped.",
        "",
        f"Run {now}, commit `{sha[:10] or 'unknown'}`, Python {platform.python_version()}, "
        f"OCP {fm.io.ocp_version()}, {platform.system()} {platform.machine()}.",
        *sections,
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
        lines += ["", "**Failed tests:**", ""] + [f"- `{n}`" for n in outcomes.failed[:50]]
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
