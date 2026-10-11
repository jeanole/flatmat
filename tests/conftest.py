import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))

import parts  # noqa: E402

import flatmat as fm  # noqa: E402


@pytest.fixture(scope="session")
def step_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("step")


def write_and_load(part: parts.Part, directory) -> fm.StepModel:
    path = os.path.join(directory, f"{part.name}.step")
    parts.write_step(parts.make_solid(part), path)
    return fm.load_step(path)


@pytest.fixture(scope="session")
def generated(step_dir):
    """name -> (Part, SheetPart), each part round-tripped through STEP."""
    out = {}
    for p in parts.cases():
        out[p.name] = (p, fm.recognize(write_and_load(p, step_dir)))
    return out
