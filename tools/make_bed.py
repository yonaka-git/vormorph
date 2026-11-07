"""格子ファイルの母点で河床高を求め、セルの値のファイルを作る。

使い方::

    uv run --project tools python tools/make_bed.py examples/domains/bend.toml \\
        output/mesh/bend.mesh.toml -o output/mesh/bend.cells.toml

領域の定義（`make_mesh.py` と同じ TOML）の `[bed]`（中心線に沿った勾配、または平面）から、格子
ファイルの各母点の河床高を求め、格子と同じセルの順に書き出す。vormorph の設定では
`[mesh] kind = "file"` の `cells` にこのファイルを指定する（ADR-0020）。格子を作り直さずに、
河床の定義だけを変えて書き直せる。詳細は docs/design/mesh-generation.md。
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from vormesh import Domain, cellsfile, meshfile
from vormesh.bed import CenterlineBed
from vormesh.domain import DomainError


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("domain", type=Path, help="domain definition with [bed] (TOML)")
    parser.add_argument("mesh", type=Path, help="mesh file made by make_mesh.py (TOML)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="cell values file to write (TOML)")
    args = parser.parse_args(argv)

    try:
        domain = Domain.load(args.domain)
        mesh, header = meshfile.read(args.mesh)
    except (DomainError, OSError, ValueError, KeyError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    if domain.bed is None:
        print(f"error: {args.domain}: [bed] is required to make the bed elevations", file=sys.stderr)
        return 1

    zb = domain.bed.elevation(mesh.x_cell, mesh.y_cell)
    kind = "centerline" if isinstance(domain.bed, CenterlineBed) else "plane"
    meta = {
        "generator": "tools/make_bed.py",
        "domain": args.domain.as_posix(),
        "mesh": args.mesh.as_posix(),
        "bed": kind,
    }
    title = str(header.get("title", domain.title))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cellsfile.write(args.output, title, mesh.x_cell, mesh.y_cell, zb, meta)

    print(f"{args.output}: {title}, {len(zb)} cells, bed '{kind}'")
    print(f"zb             : {zb.min():.4f} .. {zb.max():.4f} m")
    if isinstance(domain.bed, CenterlineBed):
        print(f"centerline     : {domain.bed.length:.3f} m, slope {domain.bed.slope}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
