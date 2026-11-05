"""計算領域の定義（TOML）から vormorph の格子ファイルを作る。

使い方::

    uv run --project tools python tools/make_mesh.py examples/domains/bend.toml \\
        -o output/mesh/bend.mesh.toml --plot output/mesh/bend.png

gmsh で境界の分割と初期の母点を作り、領域で切り取ったボロノイ図を Lloyd 法で重心ボロノイ分割に
近づけ、セルをヒルベルト曲線の順に並べて書き出す。詳細は docs/design/mesh-generation.md。
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
from collections.abc import Sequence
from pathlib import Path

from vormesh import Domain, generate, meshfile
from vormesh.domain import DomainError
from vormesh.hilbert import DEFAULT_ORDER
from vormesh.voronoi import MeshError


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("domain", type=Path, help="domain definition (TOML)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="mesh file to write (TOML)")
    parser.add_argument("--plot", type=Path, default=None, help="also draw the mesh (PNG etc.)")
    parser.add_argument("--cell-size", type=float, default=None, help="override cell_size [m]")
    parser.add_argument("--iterations", type=int, default=None, help="override lloyd.iterations")
    parser.add_argument("--verbose", action="store_true", help="show gmsh messages")
    args = parser.parse_args(argv)

    try:
        domain = Domain.load(args.domain)
        if args.cell_size is not None:
            domain = dataclasses.replace(domain, cell_size=args.cell_size)
        if args.iterations is not None:
            lloyd = dataclasses.replace(domain.lloyd, iterations=args.iterations)
            domain = dataclasses.replace(domain, lloyd=lloyd)
        result = generate(domain, verbose=args.verbose)
    except (DomainError, MeshError, OSError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 1

    lloyd = result.lloyd
    last = lloyd.displacement[-1] if lloyd.displacement else float("nan")
    last_max = lloyd.max_displacement[-1] if lloyd.max_displacement else float("nan")
    meta = {
        "generator": "tools/make_mesh.py",
        "domain": args.domain.as_posix(),
        "cell_size": domain.cell_size,
        "lloyd_iterations": lloyd.iterations,
        "lloyd_rms_displacement": float(f"{last:.3e}"),
        "hilbert_order": DEFAULT_ORDER,
        "domain_area": result.quality.domain_area,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    meshfile.write(args.output, result.mesh, domain.title, "hilbert", meta)

    print(f"{args.output}: {domain.title}")
    print(
        f"lloyd          : {lloyd.iterations} iterations, displacement / h = {last:.1e} (rms), "
        f"{last_max:.1e} (max); {lloyd.n_removed} crowded generators removed first"
    )
    for line in result.quality.lines():
        print(line)
    if args.plot is not None:
        from vormesh.plot import plot_mesh

        args.plot.parent.mkdir(parents=True, exist_ok=True)
        plot_mesh(result.mesh, args.plot, title=f"{domain.title}: {result.mesh.n_cells} cells")
        print(f"wrote {args.plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
