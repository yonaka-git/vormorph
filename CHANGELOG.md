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
- 時間ループ `model::Model`：`with_state`（与えた格子と状態から組み立てる。乾いたセルの水位と壁の流速を整え、配列の長さと境界条件の対応を検査して `ModelError::Setup` を返す）、`step`（時刻は開始時刻 + ステップ数 × Δt）、`run`、`write_output`、水の収支（上流端・下流端の流量は θ 平均流量から）
- `hydro.manning_n = 0`（摩擦なし）を許す。`[sediment]` があるときは正の値を要求する
- 検証ケース V03（静水保持）と V04（閉鎖水域の振動）を実装。V03 は流速・水位の変化が厳密に 0、V04 は体積の相対誤差が $10^{-13}$ 程度。V04 の解析解との比較は、移流がないため振動が減衰するので、移流の実装後に基準を決める。結果は `docs/verification.md` に記録
- 環境変数 `VORMORPH_VERIFICATION_OUTPUT` で、検証ケースの VTK を書き出す
- `cargo test` のビルドを最適化する（`[profile.test] opt-level = 2`。検証テストの実行時間を短くするため）
- 開境界 `hydro::HydroSolver::step`：上流端の流量境界（断面で共通の水位から求めた通水能で配分）、下流端の水位境界と水深一定の境界（水位を境界辺の上に置き、母点から辺までの距離で勾配を評価）
- ADR-0007（開境界の水位を境界辺の上に置き、流量は断面で共通の水位から求めた通水能で配分する）。`docs/design/overview.md` の第6章（境界条件）と第7章（辺の水深を風上からとることによるクーラン数の制約）を更新
- 検証ケース V06（等流）を実装。3つの格子で中央部の水深の相対誤差が $10^{-9}$ 程度（合格基準を $10^{-4}$ に確定）。図 `docs/images/v06_uniform_flow.png` と、検証ケースの図を描く `tools/plot_verification.py`（pytest のテストつき）
- 運動量の移流 `hydro::advection`（保存型の1次風上差分）：セルの運動量収支を流入側の辺から求め、時刻 n の流量で連続式を陽的に進めた水深の和で割って辺へ配る。局所クーラン数が 1 を超える辺では増分を縮める（安全策）。`StepReport` に局所クーラン数の最大値 `max_courant` と、増分を縮めた辺の数 `limited_edges` を追加
- 設定 `[hydro] advection = "upwind" | "none"`（既定は `"upwind"`。`"none"` は移流項を省く検証用）
- ADR-0008（移流は、セルの運動量収支を、陽的に進めた水深の和で割って辺へ配る）。`docs/design/overview.md` の第4.4節（移流の離散式）、第5章、第7章、第8章を更新
- 検証ケース V04 に、Thacker の解析解との比較の合格基準を追加（水塊の重心の誤差が 1周期後 0.12 m・3周期後 0.25 m 以下、水深の L1 誤差が 1周期後 20%・3周期後 40% 以下）。図 `docs/images/v04_thacker.png` と `tools/plot_verification.py v04`。VTK を書き出すときは、移流なしの計算と重心の時系列（`centroid.csv`）も書く。検証の図の色を、検証済みのカテゴリ配色にそろえる
- 検証ケース V05（ダム破壊）を実装。SWASHES と同じ条件の Stoker（下流が湛水）と Ritter（下流が乾燥）の解析解と比べ、格子幅 0.1〜0.0125 m で L1 誤差の平均の収束次数 0.65〜0.70、Stoker の衝撃波の位置の誤差が格子幅の 1.2 倍以下。図 `docs/images/v05_dam_break.png` と `tools/plot_verification.py v05`
- 検証ケース V07（跳水）を実装。SWASHES と同じ条件のこぶを越える流れ（頂部で射流に変わり、下流側の斜面で跳水する）を正方形配置で解き、跳水の位置の誤差が格子幅の 0.6 倍、跳水の直前と直後の水深の誤差が格子幅 0.05 m で 0.5% 以内。図 `docs/images/v07_hydraulic_jump.png` と `tools/plot_verification.py v07`。六角形配置で跳水の位置に横断方向の流れが生じる問題を、未決事項として記録
- 移流を入れた後の V03・V06 の結果を `docs/verification.md` に記録（どちらも同じ合格基準で通る）
- 移流の「自分の流速」を、セル中心の流速の法線成分から辺の流速に変更（ADR-0009。ADR-0008 を一部置き換え）。辺の流速を上流の流速へ局所クーラン数の割合で近づける形になり、六角形配置で、セル中心の流速に現れない格子規模の辺の流速の成分が跳水にたまって横断方向の流れが生じる問題を解消。単体テストを、1次元の定常跳水の釣り合い（Bélanger の式）と、六角形配置の格子規模の乱れの減衰に置き換え。V07 の跳水の直前・直後の水深を確かめる位置を、跳水から 0.4 m に変更
- V07 の合格基準に六角形配置（0.05 m）を加える。跳水の位置の誤差が格子幅の 0.5 倍、跳水の直前と直後の水深の誤差が 0.6% 以内で定常。図 `docs/images/v07_hydraulic_jump.png` を差し替え
- 移流を ADR-0009 の形にした後の V03〜V06 の結果を `docs/verification.md` に記録（どれも同じ合格基準で通る）
- 格子生成ツールを、平面直角座標のような原点から遠い座標でも使えるように修正：重心をセルの頂点の平均を原点にして計算し、境界の母点・Lloyd 法・組み立て・検査を、領域の近くに原点を移した局所座標で行う。$(10^5, 10^5)$ m 平行移動した領域で、以前は重心が桁落ちし、格子幅 1.5 m で Lloyd 法が収束せず、0.5 m で生成が失敗していた。平行移動のテストを追加（`tools/tests/test_voronoi.py`, `test_gmsh.py`）。`docs/design/mesh-generation.md` に「座標の扱い」を追加
- 格子の幾何量が座標の原点によらないこと（$(10^5, 10^5)$ m 平行移動して組み立てても、座標の丸めの範囲で一致する）と、一次の流速場の Perot 法の再構成の誤差が壁ぎわのセルでも小さいことを確かめるテストを追加（`tests/mesh_file.rs`、`hydro` の単体テスト）。`docs/design/mesh-generation.md` に、壁ぎわのセルの母点と重心のずれ、中心を母点のままにする理由、外周の凸な角に接するセルの制限を記録
- 設定からモデルを組み立てる `model::Model::new`（`mesh.kind = "channel"`）：直線水路の格子を作り、`[bed]` の河床高と床止め（かさ上げ・非侵食層）を母点で評価する（端にちょうど乗る母点を丸め誤差で取りこぼさないよう、格子幅の $10^{-9}$ 倍の許容差）。初期状態は乾いた河床（η = z_b、u = 0）。`mesh.kind = "file"`（格子ファイルに河床高がない）と、未実装の流砂（`[sediment]`）は `ModelError::Setup` で拒否する。格子の誤り `ModelError::Mesh` を追加
- 可変の時間刻み（ADR-0010）：`[time]` に `dt`（固定）か `courant` と `dt_max`（可変）のどちらか一方を書く。可変のときは、時刻 n の流量でセルの流出量の比 max_i Σ_流出 L q / V（`hydro::HydroSolver::max_outflow_rate`）を求め、Δt = min(C / 比, `dt_max`) とし、出力時刻・`morphology_start`・終了時刻にちょうど止める。`config::TimeStep`、`ModelError::TimeStep`。`StepReport` に時間刻み `dt` と実際のクーラン数 `outflow_courant` を追加。`docs/design/overview.md` の第5章と第7章を更新
- `model::Model::run_with`（出力のたびに直前の出力からの記録 `model::Progress` を渡す）、`Model::finished`、`Model::time_step`
- コマンドライン `vormorph run <config.toml>`：計算を実行し、出力時刻ごとに時刻・ステップ数・時間刻み・クーラン数・反復回数・水の収支の誤差を表示する。バイナリを実行するテスト `tests/cli.rs`
- 例題 `examples/sill.toml` を可変の時間刻み（`courant = 0.5`、`dt_max = 0.2`）に変更。固定の 0.2 s では、床止めの直下の射流で移流の局所クーラン数が 1.75 になり、安全策が毎ステップ約 860 辺で働いていた。乾いた状態から定常になるまで約 500 s かかるので、`morphology_start` を 300 s から 600 s に延ばした
- 掃流砂量式 `sediment::formulas`：Meyer-Peter & Müller 式、芦田・道上式（有効掃流力 τ*e = τ*）、長谷川の式による勾配補正の係数 γ。係数は原著で確認。手計算の値と比べる単体テスト
- 辺の掃流砂 `sediment::SedimentSolver::compute_bedload`（ADR-0011）：流れを進めた後の流速と、流れの計算と同じ辺の水深・流速の大きさ（`hydro::edge_speed` を公開）で掃流力を求め、勾配補正は全方向に2点差分でかける。上流端（平衡給砂）と下流端（自由流出）の辺は、内側のセルの河床勾配（最小二乗。`sediment::bed_gradient`）を外挿して同じ式で求める。非侵食層より上の土砂量で流出を頭打ちにする。`BedloadReport`
- 河床の更新 `sediment::SedimentSolver::update_bed`（ADR-0012）：Exner 式と MORFAC、安息角補正（緩和係数 0.5 の Jacobi 型の反復。辺の配列に書いてセルで集め、体積を厳密に保存。非侵食層より上の量で頭打ち）、水位の調整（水深が `dry_depth` 以上なら水位、未満なら水深を保つ）。`BedUpdateReport`
- `model::Model` に河床変動を組み込む：`morphology_start` の後に掃流砂と河床を更新し、土砂の収支（`SedimentBalance`）と、水位を保ったことによる水の体積の変化（`MassBalance::bed_exchange`）を記録する。`Model::sediment`、`Model::last_sediment_report`、`model::SedimentStepReport`。`Progress` に河床の変化の最大値と安息角補正の反復回数を追加。`[sediment]` を拒否していたのをやめ、未実装の浮遊砂（`[sediment.suspended]`）だけを拒否する
- `vormorph run` の表示に、1ステップの河床の変化の最大値と土砂の収支の誤差を追加
- 移動床（`[sediment]`）と下流端の水深一定の境界の組み合わせを、設定の誤りにする（水位が河床とともに上下し、下流端から河床低下が進むため。ADR-0011）。`examples/sill.toml` の下流端を水位固定（1.07 m）に変更
- VTK 出力に、掃流砂量ベクトル `qb`（辺の q_b を Perot 法と同じ式でセル中心へ再構成）と非侵食層の印 `fixed_bed` を追加（移動床のとき）。`VtkWriter::write_step` に辺の掃流砂量の引数を追加
- 検証ケース V12（平衡河床の保持）を追加・実装。等流に平衡給砂を与え、正方形配置と格子生成ツールの不規則な格子（凸な角を含む）で河床の変化が丸め誤差の水準（$10^{-14}$ m 以下）であることを確かめる。格子生成ツールの入力例 `examples/domains/straight.toml` と、テスト用の格子 `tests/data/straight_coarse.mesh.toml`
- `docs/design/overview.md` の第4.6節・第4.8節・第6章・第7章を更新。ADR-0011、ADR-0012
- 上流端の給砂に「一定の給砂量」（`sediment_supply = "constant"`）を追加（ADR-0013）。河床変動を始めた最初のステップの平衡給砂量の合計を保ち、各ステップの平衡給砂量に比例して辺に配る。上流端の流れが速くなる条件では上流の区間が給砂不足で下がり、流路が分化しうるので、`examples/sill.toml` は平衡給砂のまま。V12 に一定の給砂量の六角形配置 2 m のケースを追加
- ライセンス：MIT OR Apache-2.0
