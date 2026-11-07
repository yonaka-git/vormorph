# tools

vormorph の格子生成と、計算結果（VTK）の後処理に使う Python のツールです。環境は [uv](https://docs.astral.sh/uv/) で管理します。

```bash
# 依存パッケージ（gmsh, numpy, scipy, shapely, matplotlib, meshio）を tools/.venv に入れる
uv sync --project tools

# 格子生成：領域の定義（TOML）から格子ファイルと図を作る
uv run --project tools python tools/make_mesh.py examples/domains/bend.toml \
    -o output/mesh/bend.mesh.toml --plot output/mesh/bend.png

# 縦断図：セル中心を水路中心線の近くで切り出し、河床高と水位を描く
uv run --project tools python tools/plot_profile.py output/sill/sill_*.vtu -o profile.png

# 検証ケースの図（docs/verification.md）：検証テストに VTK を書き出させてから描く
VORMORPH_VERIFICATION_OUTPUT=output/verification \
    cargo test --release --test verification v06 -- --nocapture
uv run --project tools python tools/plot_verification.py v06 output/verification \
    -o docs/images/v06_uniform_flow.png

# 移流と時間積分の限界クーラン数（ADR-0016）
uv run --project tools python tools/advection_stability.py

# テスト
uv run --project tools pytest tools/tests
```

| ファイル | 内容 |
|---|---|
| `make_mesh.py` | 格子生成のコマンド。`--cell-size` と `--iterations` で領域の定義の値を上書きできる |
| `vormesh/domain.py` | 領域の定義の読み込みと、gmsh による境界の分割・初期の母点 |
| `vormesh/voronoi.py` | 境界の母点の対、Lloyd 法、多角形格子の組み立てと検査 |
| `vormesh/hilbert.py` | ヒルベルト曲線による並べ替え |
| `vormesh/meshfile.py` | 格子ファイル（vormorph mesh 形式 v1、TOML）の書き出しと読み込み |
| `vormesh/plot.py` | 格子の図（セルを番号で塗り分け、母点を番号順に結ぶ） |
| `plot_profile.py` | `.vtu` のセルデータ（`x_cell`, `y_cell`, `zb`, `eta`, `h`）から縦断図を描く |
| `plot_verification.py` | 検証ケースの図。ケースの条件は `tests/verification.rs` の定数と合わせる |
| `bar_theory.py` | 交互砂州の線形安定解析（V09 の比較対象） |
| `advection_stability.py` | 移流と時間積分の1次元の線形安定解析（限界クーラン数。ADR-0016 の根拠） |
| `tests/` | pytest のテスト。`test_gmsh.py` は gmsh を読み込めない環境では飛ばす（`VORMORPH_REQUIRE_GMSH=1` のときは失敗させる） |

手順と格子ファイルの形式は [docs/design/mesh-generation.md](../docs/design/mesh-generation.md) にあります。格子ファイルの形式は Rust 側の `src/mesh/file.rs` と、セルデータの名前は `src/io.rs` の出力の規約と合わせています。変えるときは両方を直してください。

## gmsh について

- PyPI の gmsh は、OpenGL と X11 の共有ライブラリを必要とします。Ubuntu では `sudo apt install libglu1-mesa libxft2` で入ります。ないと `import gmsh` が `OSError: libGLU.so.1: cannot open shared object file` で失敗します。
- gmsh を読み込めないと、`tests/test_gmsh.py` は飛ばされ、pytest は成功のまま終わります。環境変数 `VORMORPH_REQUIRE_GMSH=1` を設定すると、飛ばさずに失敗させます。CI と Codespaces ではこの値を設定しているので、共有ライブラリが足りなければテストが失敗して気づけます。手元で確かめるときは `VORMORPH_REQUIRE_GMSH=1 uv run --project tools pytest tools/tests` のように実行します。
- gmsh は GPL-2.0-or-later（外部ライブラリとのリンクの例外つき）で配布されています。このリポジトリは gmsh を含まず、配布もしません。
