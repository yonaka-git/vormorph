# GitHub Codespaces での実行

GitHub Codespaces を使うと、手元に何も入れずに、ブラウザだけでビルド・テスト・格子生成・計算ができます。開発環境は `.devcontainer/` に定義しています。

> **v0.1.0 の時点で動くもの**：ビルド、テスト、設定ファイルの検証（`check`）、格子生成ツールです。計算（`run`）は、ソルバーを実装すると動くようになります。

## 環境の中身

| 項目 | 内容 |
|---|---|
| ベースイメージ | `mcr.microsoft.com/devcontainers/rust:1-trixie`（Debian 13） |
| Rust | `rust-toolchain.toml` に従う（stable、rustfmt、clippy） |
| Python | uv 0.13.7。`tools/` の環境は `uv.lock` から作る（Python 本体も uv が用意する） |
| gmsh 用の共有ライブラリ | `libglu1-mesa`, `libxft2`, `libxcursor1`, `libxinerama1`, `libgomp1` |
| その他 | `tmux`（長い計算用） |
| VS Code の拡張機能 | rust-analyzer, Even Better TOML, Python, Ruff |

codespace を作ると、最後に次のコマンドが自動で実行されます（`postCreateCommand`）。初回は数分かかります。

```bash
rustup toolchain install                 # rust-toolchain.toml のツールチェーン
cargo build --release --locked           # リリースビルドを済ませておく
uv sync --project tools --locked         # tools/.venv を作る
```

## 起動する

1. README の「Open in GitHub Codespaces」のボタンを押すか、リポジトリのページで **Code → Codespaces → Create codespace on main** を選びます。
2. ターミナルに出る `postCreateCommand` のログが終わるまで待ちます。
3. 動作を確かめます。

```bash
cargo test
cargo run --release -- check examples/sill.toml
uv run --project tools pytest tools/tests
```

### マシンの種類

既定は 2 コアです。今のソルバーは逐次計算（並列化は将来 rayon で行う予定）なので、計算には 2 コアで足ります。ビルドを速くしたいときや、大きな格子でメモリが足りないときは、4 コア以上を選びます（**Code → Codespaces → … → Change machine type**、または `gh codespace create -m`）。コア数が多いほど、無料枠を早く使い切ります（[費用](#費用無料枠)を参照）。

## 長い計算をするときの注意

### アイドルタイムアウト

Codespaces は、一定時間操作がないと自動で停止します。停止すると、実行中の計算も止まります。

- 既定は 30 分で、5〜240 分の範囲で変えられます（[github.com/settings/codespaces](https://github.com/settings/codespaces) の「Default idle timeout」）。
- キーボードやマウスの操作に加えて、**ターミナルへの出力**でもタイマーがリセットされます。進み具合をターミナルに出し続ける計算なら、操作しなくても止まりません。
- 一方で、ターミナルに何も出さないバックグラウンドの処理は、操作がなければアイドルとみなされます。
- アイドルタイムアウトの設定とは別に、codespace を続けて動かせるのは **最長 12 時間** です。12 時間を超える計算は、ケースを分けて実行してください。

### 実行のしかた

tmux の中で計算を実行し、ログをファイルにも残す方法を勧めます。ブラウザとの接続が切れても計算は続き、つなぎ直したあとに `tmux attach` で戻れます。

```bash
mkdir -p output
tmux new -s sill
cargo run --release -- run examples/sill.toml 2>&1 | tee output/sill.log
```

- 計算中は、tmux の画面をターミナルに表示したままにします。tmux から切り離す（`Ctrl-b d`）と、ターミナルへの出力がなくなるので、アイドルタイムアウトで止まる可能性があります。
- 戻るときは `tmux attach -t sill` を実行します。
- `gh codespace ssh` で手元のターミナルから接続して、同じように tmux を使うこともできます。

### 停止したときに残るもの

`/workspaces/vormorph` の下のファイル（`output/` の計算結果や `target/` のビルド結果）は、codespace を停止しても残ります。停止の前に書き出した VTK ファイルは、再開したあとも使えます。

## 結果を手元に持ってくる

ParaView での可視化は手元の PC で行います。計算結果は次のどちらかの方法で取り出します。

**ブラウザの VS Code から**：エクスプローラーで `output/sill` を右クリックし、「Download...」を選びます。

**GitHub CLI から**：手元のターミナルで実行します。ファイルが多いときは、先に codespace の中で `tar czf output/sill.tar.gz -C output sill` のようにまとめておくと速く取り出せます。

```bash
gh codespace list                                      # codespace の名前を調べる
gh codespace cp -r -c <名前> remote:/workspaces/vormorph/output/sill ./output/
```

Python での後処理（図の作成）は codespace の中でもできます。できた PNG は VS Code のエディターで開いて確認できます。

```bash
uv run --project tools python tools/plot_profile.py output/sill/sill_*.vtu -o output/profile.png
```

## 費用（無料枠）

個人アカウントには、毎月の無料枠があります（2026 年 11 月時点。最新の値は [GitHub Codespaces の課金の説明](https://docs.github.com/en/billing/concepts/product-billing/github-codespaces)で確認してください）。

| プラン | 計算時間（コア時間） | ストレージ |
|---|---|---|
| GitHub Free | 120 コア時間／月（2 コアのマシンで 60 時間） | 15 GB／月 |
| GitHub Pro | 180 コア時間／月（2 コアのマシンで 90 時間） | 20 GB／月 |

- 計算時間は、codespace が動いている間だけ数えます。使い終わったら停止します（`gh codespace stop`、または **Code → Codespaces → … → Stop codespace**）。
- ストレージは、停止している codespace にもかかります。`target/` と `output/` は大きくなりやすいので、不要になったら `cargo clean` や不要な結果の削除で減らします。使わない codespace は削除します。

## 開発環境を変えるとき

- `.devcontainer/` を変えたら、コマンドパレットの **Codespaces: Rebuild Container** で作り直します。
- uv の版は `.devcontainer/Dockerfile` で固定しています。上げるときは、手元の uv と合わせます。
- gmsh の版を上げて必要な共有ライブラリが変わったら、`Dockerfile` の `apt-get install` の一覧を直します。必要なライブラリは `readelf -d tools/.venv/lib/libgmsh.so.*` で確かめられます。
