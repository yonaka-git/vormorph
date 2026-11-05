# tools

vormorph の計算結果（VTK）を後処理する Python スクリプトです。環境は [uv](https://docs.astral.sh/uv/) で管理します。

```bash
# 依存パッケージ（numpy, matplotlib, meshio）を tools/.venv に入れる
uv sync --project tools

# 縦断図：セル中心を水路中心線の近くで切り出し、河床高と水位を描く
uv run --project tools python tools/plot_profile.py output/sill/sill_*.vtu -o profile.png
```

| スクリプト | 内容 |
|---|---|
| `plot_profile.py` | `.vtu` のセルデータ（`x_cell`, `y_cell`, `zb`, `eta`, `h`）から縦断図を描く |

セルデータの名前は `src/io.rs` の出力の規約に合わせています。名前を変えるときは両方を直してください。

## 格子生成（任意）：jigsawpy

外部格子（JIGSAW の `.msh`）を作るときに使います。PyPI では配布されていないため、`pyproject.toml` の依存には入れていません。使う場合は、C++ コンパイラと cmake を用意して、公式リポジトリからビルドします。

```bash
git clone --depth 1 https://github.com/dengwirda/jigsaw-python.git ../jigsaw-python
tools/.venv/bin/python ../jigsaw-python/build.py                    # C++ の本体を cmake でビルド
uv pip install --python tools/.venv/bin/python ../jigsaw-python     # tools/.venv に入れる
```

- `uv sync` は依存に書かれていないパッケージを消します。jigsawpy を残すには `uv sync --inexact` を使います（`uv run` は消しません）。
- conda を使う場合は、conda-forge の `jigsawpy` でも入ります。
- JIGSAW は独自のライセンスです。研究・個人での利用は無料ですが、商用のシステムに組み込むには作者との取り決めが必要です。このリポジトリには JIGSAW のコードを含めません。
