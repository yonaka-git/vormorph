"""領域の定義（TOML）の読み込みと検査。gmsh は使わない。"""

import tomllib
from pathlib import Path

import numpy as np
import pytest

from vormesh.domain import Domain, DomainError

EXAMPLES = sorted((Path(__file__).parents[2] / "examples" / "domains").glob("*.toml"))

RECT = {
    "title": "rect",
    "cell_size": 1.0,
    "outer": [
        {"kind": "upstream", "line": [[0.0, 0.0], [0.0, -5.0]]},
        {"kind": "wall", "line": [[0.0, -5.0], [10.0, -5.0]]},
        {"kind": "downstream", "line": [[10.0, -5.0], [10.0, 0.0]]},
        {"kind": "wall", "line": [[10.0, 0.0], [0.0, 0.0]]},
    ],
}


def test_examples_are_valid():
    assert EXAMPLES, "examples/domains/*.toml not found"
    for path in EXAMPLES:
        domain = Domain.load(path)
        assert domain.cell_size > 0.0


def test_arcs_are_split_into_quarter_circles_or_less():
    bend = Domain.load(next(p for p in EXAMPLES if p.stem == "bend"))
    arcs = [s for s in bend.outer if s.arc_center is not None]
    assert len(arcs) == 2
    for seg in arcs:
        r = np.hypot(*(seg.points - seg.arc_center).T)
        assert np.allclose(r, r[0])
        assert len(seg.points) == 3  # 180° → 90° × 2


def test_clockwise_outer_is_rejected():
    data = {**RECT, "outer": RECT["outer"]}
    Domain.from_dict(data)  # 反時計回りなら通る
    reversed_outer = [{"kind": s["kind"], "line": s["line"][::-1]} for s in RECT["outer"][::-1]]
    with pytest.raises(DomainError, match="counterclockwise"):
        Domain.from_dict({**RECT, "outer": reversed_outer})


def test_gap_between_segments_is_rejected():
    outer = [dict(s) for s in RECT["outer"]]
    outer[1] = {"kind": "wall", "line": [[0.0, -5.0], [9.0, -5.0]]}
    with pytest.raises(DomainError, match="gap"):
        Domain.from_dict({**RECT, "outer": outer})


@pytest.mark.parametrize(
    "patch",
    [
        {"cell_sise": 1.0},
        {"lloyd": {"iteration": 10}},
        {"outer": [{"kind": "inlet", "line": [[0.0, 0.0], [1.0, 0.0]]}]},
    ],
)
def test_typos_are_rejected(patch):
    with pytest.raises(DomainError):
        Domain.from_dict({**RECT, **patch})


def test_example_files_parse_as_plain_toml():
    for path in EXAMPLES:
        with open(path, "rb") as f:
            assert "outer" in tomllib.load(f)
