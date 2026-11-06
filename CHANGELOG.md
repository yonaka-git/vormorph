# Changelog

このプロジェクトの主な変更を記録します。

形式は [Keep a Changelog](https://keepachangelog.com/ja/1.1.0/) に、版番号は [Semantic Versioning](https://semver.org/lang/ja/) に従います。

## [0.1.0] - Unreleased

### Added

- リポジトリの雛形：Cargo プロジェクト、`rust-toolchain.toml`（stable 固定）、`rustfmt.toml`、`.editorconfig`、`.gitattributes`、`.gitignore`
- モジュールの雛形：`config`、`mesh`（`channel`、`voronoi`）、`state`、`hydro`、`linsolve`、`sediment`（`formulas`）、`model`、`io`。数値計算の本体は `todo!()`
- TOML 設定の読み込みと検証（`Config::load`、未知のキーの検出、値の範囲の検査）
- コマンドライン：`vormorph check <config.toml>`（設定の検証）。`run` は未実装のエラーを返す
- 格子の接続配列（CSR 形式の `Jagged`）と、法線の符号の規約（`Mesh::derive_cell_connectivity`）
- 例題 `examples/sill.toml`（床止め下流の局所洗掘）
- 検証ケース 11 件の一覧（`tests/verification.rs`、未実装は `#[ignore]`）と `docs/verification.md`
- 設計文書 `docs/design/overview.md`（支配方程式・離散式・1ステップの手順）と ADR 0001〜0004
- 文献リスト `docs/references.md`
- Python の後処理の雛形 `tools/plot_profile.py`（uv で管理）
- 格子生成ツール `tools/make_mesh.py`（パッケージ `tools/vormesh/`）：gmsh による境界の分割と初期母点、境界に関して対称な母点の対、Lloyd 法による重心ボロノイ分割、ヒルベルト曲線の順への並べ替え、格子の検査、格子ファイル（TOML）の書き出し、格子の図。pytest のテスト
- 格子生成ツールの入力例 `examples/domains/`（180° 湾曲水路、島のある池）
- 格子ファイルの読み込み `mesh::file`（形式・版・番号の範囲の検査）と、多角形の並び `mesh::voronoi::PolygonMesh`
- 多角形の並びから MPAS 形式の接続配列と幾何量を組み立てる `mesh::voronoi::build` と `compute_geometry`。辺の共有（T 字の接合）、境界辺の一覧との照合、セルの向き、境界辺に対する母点の位置を検査し、見つかった問題を `MeshError` で列挙して返す。`mesh::file::read` はこの誤りを `MeshFileError::Mesh` で返す
- 格子の幾何の検査 `Mesh::check_geometry`（総面積、内部の辺の直交誤差、セルの閉合誤差）
- 直線水路の格子生成 `mesh::channel::build`（正方形配置と六角形配置）。水路の矩形を近傍の母点との垂直二等分線で切り取ってセルを求め、別々に求めた頂点を許容誤差で同一視する。側壁は最も外側の行の母点と鏡映点の二等分線に一致する
- 検証ケース V01（格子幾何）を実装。面積・直交誤差・閉合誤差がいずれも $10^{-12}$ 以下であることを、正方形配置と六角形配置で確かめる。結果は `docs/verification.md` に記録
- Perot 法によるセル中心の流速の再構成 `hydro::reconstruct_cell_velocity`
- 検証ケース V02（Perot 再構成）を実装。一様流で、直線水路の2つの配置と格子生成ツールの島の格子の全セル（境界セルを含む）が $10^{-12}$ m/s 以内で一致することを確かめる
- VTK 出力 `io::VtkWriter`：時刻ごとの VTK XML UnstructuredGrid（`.vtu`、ASCII）と ParaView Data（`.pvd`）。各 `.vtu` に時刻（`TimeValue`）も書く。浮動小数点数は読み戻すと同じ値になる最短の10進表記で書き、NaN・無限大や配列の長さの不一致は書き出さずに誤りにする。濃度 `c` の出力に使う `dry_depth` を `VtkWriter::new` の引数に追加。読み戻して確かめるテスト `tests/vtk_output.rs`（meshio と `tools/plot_profile.py` でも読めることを確認済み）
- 格子生成ツールの出力を Rust 側で読み、面積・辺の共有・直交性・境界の母点の位置・番号の局所性を確かめるテスト `tests/mesh_file.rs`（データは `tests/data/island_coarse.mesh.toml`）。Rust 側で組み立てた接続配列と幾何量（境界辺の数、オイラー標数、`check_geometry` の結果）も確かめる
- 設計文書 `docs/design/mesh-generation.md` と ADR-0005（格子生成は gmsh を使った自作ツールで行い、セルをヒルベルト曲線の順に並べる）
- CI：fmt・clippy・test（ubuntu・windows）、MSRV の check、`cargo doc`、`tools/` の pytest
- 環境変数 `VORMORPH_REQUIRE_GMSH=1` で、gmsh を読み込めないときに `tools/tests/test_gmsh.py` を飛ばさずに失敗させる。CI と Codespaces で設定し、共有ライブラリの不足を見逃さないようにする
- GitHub Codespaces の開発環境 `.devcontainer/`（Rust、uv、gmsh 用の共有ライブラリ、tmux、SSH サーバー。作成時にリリースビルドと `uv sync` を実行）と、使い方の文書 `docs/codespaces.md`（アイドルタイムアウト、結果の取り出し、無料枠）。README に「Open in GitHub Codespaces」のボタン
- 行列フリーの Jacobi 前処理付き共役勾配法 `linsolve::pcg`。右辺が 0 のとき、作業配列の長さが違うとき、正定値でない行列（対角が負、pᵀAp ≤ 0）、反復の打ち切りを単体テストで確かめる
- 半陰解法の1ステップ `hydro::HydroSolver::step`（壁で囲んだ水域）：辺の水深 `hydro::edge_depths`（風上の水位と高い方の河床高、`dry_depth` 未満の辺は乾いた辺）、陰的な摩擦（|u| は Perot 流速の接線成分を含む）、行列フリーの係数行列 `hydro::WaterLevelOperator`、修正量の形の Newton 法と PCG、θ 平均流量、乾いたセルの水位を河床高に揃える処理。移流は次の段階で実装し、それまでは移流項を省く
- ADR-0006（wet/dry を辺の水深で判定し、水位の連立方程式を修正量の形の Newton 法で解く）。`docs/design/overview.md` に Newton 法の初期値と単調収束の根拠を追記
- ライセンス：MIT OR Apache-2.0
