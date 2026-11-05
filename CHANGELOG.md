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
- 設計文書 `docs/design/overview.md`（支配方程式・離散式・1ステップの手順）と ADR 4 件
- 文献リスト `docs/references.md`
- Python の後処理の雛形 `tools/plot_profile.py`（uv で管理）
- CI：fmt・clippy・test（ubuntu・windows）、MSRV の check、`cargo doc`
- ライセンス：MIT OR Apache-2.0
