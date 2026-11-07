//! # vormorph
//!
//! ボロノイ格子・半陰解法による平面二次元河床変動モデル。
//!
//! **現在は v0.1.0 の開発中です。** 流れの半陰解法（wet/dry、流量・水位境界、陰的な摩擦、
//! 保存型風上差分の移流）、固定または可変の時間刻み、掃流砂と河床変動（Exner 式、MORFAC、
//! 非侵食層、安息角補正、二次流による掃流砂の向きの補正）は実装済みです。設定から直線水路を組み立てる [`model::Model::new`] か、
//! 格子と状態を与える [`model::Model::with_state`] から使えます。浮遊砂は未実装で、該当する
//! 関数の本体は `todo!()` になっています（`[sediment.suspended]` のある設定は、モデルの
//! 組み立ての時点でエラーにします）。離散式と1ステップの手順は
//! `docs/design/overview.md`、設計判断は `docs/design/adr/` を参照してください。
//!
//! ## モジュール構成
//!
//! | モジュール | 役割 |
//! |---|---|
//! | [`config`] | TOML 設定の読み込みと検証。物理定数の受け渡し口 |
//! | [`mesh`] | ボロノイ格子（MPAS 形式の接続配列）、直線水路の格子生成、格子ファイルの読み込み |
//! | [`state`] | 状態変数（セル中心の η, z_b, hC, I、境界面の法線流速 u_n） |
//! | [`hydro`] | Casulli 型 θ 法の半陰解法、保存型風上差分の移流、Perot 法による流速の再構成 |
//! | [`linsolve`] | 行列フリーの前処理付き共役勾配法（PCG） |
//! | [`sediment`] | 掃流砂・Exner 式・安息角補正・二次流（浮遊砂は未実装） |
//! | [`transport`] | セル中心のスカラーの移流拡散（二次流の強さ。浮遊砂と共用する予定） |
//! | [`model`] | 時間ループと質量収支 |
//! | [`io`] | ParaView 用の VTK 出力 |
//!
//! ## 格子と符号の規約（全モジュール共通）
//!
//! - 境界面 `e` の法線 `n_e` は `cells_on_edge[e][0]` から外向き（`[0]` から `[1]` へ向く）。
//! - `edge_sign_on_cell[c][k]` は、辺 `edges_on_cell[c][k]` の法線がセル `c` から見て外向きなら
//!   `+1.0`、内向きなら `-1.0`。
//! - 領域の外側を表すセル番号は [`mesh::NONE`]（`usize::MAX`）。境界辺では
//!   `cells_on_edge[e][1] == NONE` で、法線は領域の外を向く。
//! - 接続配列と幾何量の名前は MPAS 形式（`cellsOnEdge` → `cells_on_edge` のように snake_case）。
//! - 物理量はすべて `f64`、SI 単位（角度のみ設定ファイルで度）。
//!
//! ## 参考文献
//!
//! 文献の一覧は `docs/references.md` にあります。各モジュールの冒頭に、そのモジュールが
//! 依拠する文献を記載しています。

pub mod config;
pub mod hydro;
pub mod io;
pub mod linsolve;
pub mod mesh;
pub mod model;
pub mod sediment;
pub mod state;
pub mod transport;
