"""Feature kinds: the table every stage asks, and the places that print each kind."""

from typing import get_args

import flatmat as fm
from flatmat import cli, export
from flatmat.model import FeatureKindName


def test_every_kind_is_described_everywhere():
    # a new kind must be added to the type, the table and both label lookups
    kinds = set(get_args(FeatureKindName))
    assert set(fm.FEATURE_KINDS) == kinds
    assert set(export._LABELS) == kinds
    assert set(cli._DETAILS) == kinds


def test_kind_properties_hang_together():
    for name, k in fm.FEATURE_KINDS.items():
        assert not (k.owns_holes and k.flat in ("hole", "rim")), name  # holes it owns are cut by it
        assert not k.tool_bent or k.formed, name
        assert not k.oriented or k.flat == "slits", name
