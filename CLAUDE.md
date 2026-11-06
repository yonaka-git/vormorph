# CLAUDE.md

このファイルは、このリポジトリで作業する Claude Code への指示書です。作業を始める前に必ず読んでください。

## プロジェクトの目的

- vormorph は、ボロノイ格子と半陰解法による平面二次元河床変動モデルです（Rust）。
- 作者の個人 GitHub で公開し、転職活動のポートフォリオとして使います。コードの品質と同じくらい、README と設計文書の読みやすさを重視します。
- 現在は v0.1.0（雛形）の段階です。数値計算の本体は `todo!()` のままです。

## 設計決定（変えるときは ADR を追加する）

| 項目 | 決定事項 |
|---|---|
| 支配方程式 | 水深平均の浅水流方程式、掃流砂・浮遊砂、Exner 式 |
| 格子 | ボロノイ格子（母点を結ぶ線が境界面に直交する性質を使う） |
| 変数配置 | C 格子（スタガード）：セル中心に η, zb, hC、境界面に法線流速 u_n |
| 時間積分 | Casulli 型 θ 法の半陰解法（θ = 0.55〜0.6）。水位の連立方程式は対称正定値 → Newton（wet/dry の体積の折れ点）＋行列フリー PCG |
| 移流 | 陽的な保存型風上差分（Stelling–Duinmeijer / Kramer–Stelling 系）。対象河道に跳水・射流（堰・床止め）があるため。将来、加速流ではエネルギー水頭保存形に切り替え |
| 摩擦 | マニング式、陰的に扱う |
| 流速の再構成 | Perot 法（セル中心の流速ベクトル） |
| 掃流砂 | 境界面で法線成分を評価。MPM／芦田・道上を切り替え可能。横断勾配補正は直交性を使った2点差分。非侵食層（床止め）で流出量を頭打ち |
| 浮遊砂 | hC の移流拡散（連続式と同じ θ 平均流量を使用）。沈降は陰的に扱い、浅いセルでは浮上を止める |
| 河床 | Exner 式、MORFAC、安息角補正。河床更新時は湛水セルの水位を保持 |
| 境界条件 | 上流：流量（通水能で配分）、下流：水位固定／水深一定、壁 |
| 格子生成 | gmsh を使った自作ツール（`tools/make_mesh.py`）。境界は鏡映対の母点、内部は重心ボロノイ分割。セルはヒルベルト曲線の順（ADR-0005） |
| 言語 | Rust（stable）。外部依存は最小限（serde, toml 程度から始める） |
| 出力 | ParaView 用 VTK |

離散式と1ステップの手順は `docs/design/overview.md`、格子生成の手順と格子ファイルの形式は `docs/design/mesh-generation.md`、判断の理由は `docs/design/adr/` にあります。

## リポジトリの構成

| パス | 内容 |
|---|---|
| `src/config.rs` | TOML 設定（物理定数の受け渡し口） |
| `src/mesh/` | ボロノイ格子（MPAS 形式）、多角形の並びからの組み立て、直線水路の生成、格子ファイルの読み込み |
| `src/state.rs` | 状態変数 |
| `src/hydro.rs` | 半陰解法、Perot 再構成、移流 |
| `src/linsolve.rs` | 行列フリー PCG |
| `src/sediment/` | 掃流砂・浮遊砂・Exner 式、流砂量式 |
| `src/model.rs` | 時間ループ、質量収支 |
| `src/io.rs` | VTK 出力 |
| `tests/verification.rs` | 検証ケース（`docs/verification.md` と 1 対 1） |
| `examples/` | 例題の設定ファイル。`examples/domains/` は格子生成ツールの入力 |
| `tools/` | Python の格子生成（`make_mesh.py`, `vormesh/`）と後処理（uv で管理）。テストは `tools/tests/` |
| `references/` | 文献 PDF（コミットしない） |
| `.devcontainer/` | GitHub Codespaces の開発環境（使い方は `docs/codespaces.md`） |

## コーディング規約

- 浮動小数点数は `f64` だけを使う。
- `unsafe` は使わない（`Cargo.toml` の `unsafe_code = "forbid"` で禁止している）。
- 数値計算のループはセル単位で書く。各セルは自分の値だけを書き込む（将来 rayon で並列化するため）。辺のループで両側のセルへ足し込む書き方（scatter）は避ける。
- 物理定数（重力加速度・密度・動粘性係数など）はハードコードせず、`config` の構造体から渡す。既定値は `src/config.rs` にだけ置く。
- 格子の規約を守る。
  - 境界面の法線は `cells_on_edge[e][0]` から外向き
  - `edge_sign_on_cell` は、そのセルから見て外向きなら +1（MPAS-Ocean の `edgeSignOnCell` とは逆符号）
  - 境界の外側を表すセル番号は `NONE = usize::MAX`
  - 接続配列は MPAS 形式（`cells_on_edge`, `edges_on_cell`, `dc_edge`, `dv_edge`, `area_cell` など）
- 公開する項目には `///` のドキュメントを書く（`missing_docs` を警告にしている）。モジュールの冒頭の `//!` には、役割・変数の規約・参照する文献を書く。
- ドキュメントとコメントは日本語、識別子は英語、実行時のメッセージ（エラーなど）は英語にする。
- ライブラリのコードでは `panic!` / `unwrap` を避けて `Result` を返す。ただし雛形の段階の `todo!()` は除く。
- 依存クレートを追加するときは、理由を ADR かプルリクエストの説明に書く。

## 検証の方針

**検証テストが通らない変更はマージしない。** 変更のたびに、次がすべて通ることを確かめる。

```bash
cargo build
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test
uv run --project tools pytest tools/tests   # tools/ を変えたとき
```

- 検証ケースを実装したら、`tests/verification.rs` の `#[ignore]` を外し、`docs/verification.md` と README の状態欄（✅ / 🚧 / 📋）を更新する。
- 離散式を変えたら `docs/design/overview.md` を、設計判断を変えたら ADR を、同じコミットで更新する。
- README に、未実装の機能を実装済みのように書かない。
- 変更点は `CHANGELOG.md` の該当する版に追記する。

## コミットメッセージ

Conventional Commits 形式で書く。

```text
<type>(<scope>): <summary>

type : feat, fix, docs, test, refactor, perf, chore, ci, build
scope: config, mesh, hydro, linsolve, sediment, model, io, tools など（任意）
例   : feat(hydro): implement Perot velocity reconstruction
```

## GitHub の運用

作業ブランチの push、PR の作成とマージは、作者の確認を取らずに次の手順で行ってよい。

1. main から作業ブランチを作る（`feat/...`, `fix/...`, `docs/...` など）。main には直接コミットしない。
2. コミットのたびに「検証の方針」のコマンドをすべて通す。
3. push して PR を作る。
   - タイトルは Conventional Commits 形式（英語）。本文は日本語で、変更点・検証の結果（検証ケースの値を含む）・設計上の判断を書く。
   - 本文に生成ツールの表記（「Generated with Claude Code」など）は入れない。AI との協働は、コミットの `Co-Authored-By` で示す。
4. CI がすべて通ったら、merge commit でマージしてブランチを消す（`gh pr merge --merge --delete-branch`）。`docs/verification.md` がコミットのハッシュを参照するので、squash と rebase は使わない。検証の結果にハッシュを書き足すコミットは、同じ PR に入れる。
5. CI が失敗したら、ブランチで直して push し直す。失敗したままマージしない。直せないときは作者に報告する。

## 恒久ルール

1. main への直接の push、force push などの履歴の書き換え、リポジトリの作成・設定の変更・削除、タグと GitHub Releases の作成、issue の作成は、必ず作者の確認を取ってから行う（作業ブランチと PR は「GitHub の運用」の手順で行ってよい）。
2. メールの送信、他人が作った issue・PR・Discussion へのコメント、掲示板などへの書き込みは一切しない。
3. 作者の勤務先の業務データ・コード・社名を含めない。例題・テストデータは、すべて合成データか公開データにする。
4. 参考文献の PDF はリポジトリにコミットしない（著作権のため）。`references/` に置き、`.gitignore` で除外する。一覧は `docs/references.md`。

## 開発環境のメモ

- 作者の環境は WSL2 上の Ubuntu で、リポジトリは Windows 側のドライブ（`/mnt/c/...`）にある。ファイルアクセスが遅いので、ビルドに時間がかかっても異常ではない。
- Python の後処理と格子生成は `tools/` で uv を使う。例：`uv run --project tools python tools/plot_profile.py output/sill/*.vtu`
- PyPI の gmsh は `libGLU.so.1` と `libXft.so.2` を必要とする（Ubuntu では `sudo apt install libglu1-mesa libxft2`）。ないと `import gmsh` が `OSError` になり、`tools/tests/test_gmsh.py` は飛ばされる。ただし `VORMORPH_REQUIRE_GMSH=1`（CI と Codespaces で設定）のときは失敗させる。
- 格子ファイルの形式（`tools/vormesh/meshfile.py` と `src/mesh/file.rs`）を変えるときは、両方と `docs/design/mesh-generation.md` を同じコミットで直し、`tests/data/island_coarse.mesh.toml` を作り直す（コマンドは `tests/mesh_file.rs` の冒頭）。
- 計算は GitHub Codespaces でも実行できる（`.devcontainer/`、手順は `docs/codespaces.md`）。アイドルタイムアウト（ターミナルへの出力でリセット）と最長 12 時間の制限がある。gmsh が必要とする共有ライブラリを変えたら `.devcontainer/Dockerfile` も直す。
- CI（`.github/workflows/ci.yml`）は ubuntu と windows で fmt・clippy・test を、ubuntu で MSRV（1.85）の check と `cargo doc`、`tools/` の pytest を実行する。
