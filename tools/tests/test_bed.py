"""河床の定義（領域の定義の [bed]）、セルの値のファイル、tools/make_bed.py のテスト。gmsh は使わない。"""

import math
from pathlib import Path

import numpy as np
import pytest

import make_bed
from vormesh import cellsfile, meshfile
from vormesh.bed import BedError, CenterlineBed, PlaneBed, parse_bed
from vormesh.domain import Domain, DomainError

ROOT = Path(__file__).parents[2]
BEND = ROOT / "examples" / "domains" / "bend.toml"
STRAIGHT = ROOT / "examples" / "domains" / "straight.toml"
STRAIGHT_MESH = ROOT / "tests" / "data" / "straight_coarse.mesh.toml"
STRAIGHT_CELLS = ROOT / "tests" / "data" / "straight_coarse.cells.toml"


def test_examples_define_their_beds():
    bend = Domain.load(BEND).bed
    assert isinstance(bend, CenterlineBed)
    # 40 m + π 30 m + 40 m
    assert bend.length == pytest.approx(80.0 + 30.0 * math.pi, rel=1e-12)
    assert isinstance(Domain.load(STRAIGHT).bed, PlaneBed)


def test_centerline_distance_is_exact_on_lines_and_arcs():
    bed = Domain.load(BEND).bed
    radius = 30.0
    # 上流の直線区間（x = 30、y = −40 → 0）の中心線上と、横にずれた点
    assert bed.distance(np.array([[30.0, -40.0], [30.0, -10.0], [33.0, -10.0]])) == pytest.approx(
        [0.0, 30.0, 30.0]
    )
    # 湾曲部：角度 φ（x 軸から反時計回り）の点は、半径によらず s = 40 + 30 φ
    for phi in (0.3, math.pi / 2, 2.9):
        for r in (26.0, 30.0, 34.0):
            point = np.array([[r * math.cos(phi), r * math.sin(phi)]])
            assert bed.distance(point)[0] == pytest.approx(40.0 + radius * phi, abs=1e-12)
    # 下流の直線区間（x = −30、y = 0 → −40）
    end = 40.0 + radius * math.pi
    assert bed.distance(np.array([[-28.0, -25.0]]))[0] == pytest.approx(end + 25.0, abs=1e-12)
    # 河床高：下流端で z_downstream、上流端で z_downstream + slope × 全長
    z = bed.elevation(np.array([30.0, -30.0]), np.array([-40.0, -40.0]))
    assert z == pytest.approx([0.002 * bed.length, 0.0], abs=1e-14)


def test_clockwise_arcs_and_points_outside_the_arc():
    bed = parse_bed(
        {
            "kind": "centerline",
            "slope": 0.01,
            "z_downstream": 1.0,
            "centerline": [{"arc": {"center": [0.0, 0.0], "radius": 10.0, "start": 90.0, "end": 0.0}}],
        }
    )
    # 時計回り（90° → 0°）：45° の点は s = 10 π/4
    point = np.array([[7.0 * math.cos(math.pi / 4), 7.0 * math.sin(math.pi / 4)]])
    assert bed.distance(point)[0] == pytest.approx(10.0 * math.pi / 4, abs=1e-12)
    # 円弧の範囲の外の点は、近い方の端に射影する
    near_start = np.array([[-1.0, 10.0]])  # 90° の端（始点）の少し先
    near_end = np.array([[10.0, -1.0]])  # 0° の端（終点）の少し先
    assert bed.distance(near_start)[0] == pytest.approx(0.0)
    assert bed.distance(near_end)[0] == pytest.approx(10.0 * math.pi / 2)


def test_plane_bed():
    bed = parse_bed({"kind": "plane", "z0": 0.5, "gradient": [-0.002, 0.01]})
    assert bed.elevation(np.array([10.0]), np.array([2.0])) == pytest.approx([0.5 - 0.02 + 0.02])


@pytest.mark.parametrize(
    ("spec", "message"),
    [
        ({"kind": "cone"}, "kind must be one of"),
        ({"kind": "plane", "z0": 0.0}, "gradient"),
        ({"kind": "plane", "z0": 0.0, "gradient": [0.0, 0.0], "slope": 1.0}, "unknown keys"),
        ({"kind": "centerline", "slope": 0.001, "z_downstream": 0.0}, "centerline"),
        (
            {
                "kind": "centerline",
                "slope": 0.001,
                "z_downstream": 0.0,
                "centerline": [{"line": [[0.0, 0.0], [1.0, 0.0]]}, {"line": [[2.0, 0.0], [3.0, 0.0]]}],
            },
            "gap",
        ),
        (
            {
                "kind": "centerline",
                "slope": 0.001,
                "z_downstream": 0.0,
                "centerline": [{"line": [[0.0, 0.0], [1.0, 0.0]], "arc": {}}],
            },
            "exactly one",
        ),
    ],
)
def test_invalid_bed_definitions_are_rejected(spec, message):
    with pytest.raises(BedError, match=message):
        parse_bed(spec)


def test_domain_reports_bed_errors():
    data = {
        "title": "rect",
        "cell_size": 1.0,
        "outer": [
            {"kind": "upstream", "line": [[0.0, 0.0], [0.0, -5.0]]},
            {"kind": "wall", "line": [[0.0, -5.0], [10.0, -5.0]]},
            {"kind": "downstream", "line": [[10.0, -5.0], [10.0, 0.0]]},
            {"kind": "wall", "line": [[10.0, 0.0], [0.0, 0.0]]},
        ],
        "bed": {"kind": "plane", "z0": 0.0},
    }
    with pytest.raises(DomainError, match=r"\[bed\]"):
        Domain.from_dict(data)
    del data["bed"]
    assert Domain.from_dict(data).bed is None


def test_cells_file_round_trip(tmp_path: Path):
    x, y = np.array([0.1, 1.0 / 3.0]), np.array([-2.5, 7.0])
    zb, z_ne = np.array([1e-17, 0.25]), np.array([-1.0, 0.2])
    path = tmp_path / "a.cells.toml"
    cellsfile.write(path, 'odd "title"', x, y, zb, {"generator": "test"}, z_ne=z_ne)
    values, header = cellsfile.read(path)
    assert header["title"] == 'odd "title"' and header["meta"]["generator"] == "test"
    for key, expected in (("x", x), ("y", y), ("zb", zb), ("z_ne", z_ne)):
        assert np.array_equal(values[key], expected), key
    cellsfile.write(path, "no z_ne", x, y, zb, {})
    assert "z_ne" not in cellsfile.read(path)[0]
    with pytest.raises(ValueError, match="same length"):
        cellsfile.write(path, "bad", x, y, zb[:1], {})


def test_make_bed_writes_the_fixture(tmp_path: Path):
    output = tmp_path / "straight.cells.toml"
    assert make_bed.main([str(STRAIGHT), str(STRAIGHT_MESH), "-o", str(output)]) == 0
    values, header = cellsfile.read(output)
    mesh, mesh_header = meshfile.read(STRAIGHT_MESH)
    assert header["title"] == mesh_header["title"] == "straight"
    assert np.array_equal(values["x"], mesh.x_cell) and np.array_equal(values["y"], mesh.y_cell)
    assert values["zb"] == pytest.approx(0.12 - 0.002 * mesh.x_cell, abs=1e-15)
    # tests/data/straight_coarse.cells.toml は、このコマンドで作ったもの（パスの meta だけ違う）
    fixture, _ = cellsfile.read(STRAIGHT_CELLS)
    for key in ("x", "y", "zb"):
        assert np.array_equal(fixture[key], values[key]), key


def test_make_bed_needs_a_bed_definition(tmp_path: Path, capsys):
    domain = tmp_path / "no_bed.toml"
    text = STRAIGHT.read_text(encoding="utf-8")
    domain.write_text(text[: text.index("[bed]")], encoding="utf-8")
    assert make_bed.main([str(domain), str(STRAIGHT_MESH), "-o", str(tmp_path / "x.toml")]) == 1
    assert "[bed] is required" in capsys.readouterr().err
