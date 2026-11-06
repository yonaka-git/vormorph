//! 検証ケース。`docs/verification.md` の表と 1 対 1 に対応する。
//!
//! 未実装のケースは `#[ignore = "..."]` を付けて一覧として残している。実装したら
//! `#[ignore]` を外し、`docs/verification.md` と README の状態欄を更新する。
//!
//! ```text
//! cargo test                              # 実装済みのケースだけ
//! cargo test -- --ignored                 # 未実装のケースも実行（現状は todo!() で失敗する）
//! cargo test --test verification v01      # 1ケースだけ
//! ```

use std::path::Path;
use vormorph::config::{ChannelLayout, Config};
use vormorph::hydro;
use vormorph::mesh::channel::{self, ChannelSpec};
use vormorph::mesh::{Mesh, file};

const LAYOUTS: [ChannelLayout; 2] = [ChannelLayout::Rectangular, ChannelLayout::Hexagonal];

/// 検証用の直線水路（延長 200 m、幅 20 m）。
fn channel_mesh(layout: ChannelLayout, cell_size: f64) -> Mesh {
    channel::build(&ChannelSpec {
        length: 200.0,
        width: 20.0,
        cell_size,
        layout,
    })
    .expect("the verification channel must be a valid mesh")
}

/// 格子生成ツール（`tools/make_mesh.py`）で作った、島のある池の格子。
///
/// 内部は重心ボロノイ分割に近いが母点はセルの重心と一致せず、境界は鏡映対の母点で作った
/// 不規則な格子（`tests/mesh_file.rs` と同じデータ）。
fn island_mesh() -> Mesh {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/island_coarse.mesh.toml");
    file::read(&path).expect("the island fixture must be a valid mesh")
}

/// V01 格子幾何：総面積、直交誤差、各セルが閉じていること。
///
/// 結果の値は `cargo test --test verification v01 -- --nocapture` で表示される
/// （`docs/verification.md` の記録に使う）。
#[test]
fn v01_mesh_geometry() {
    for layout in LAYOUTS {
        let mesh = channel_mesh(layout, 2.0);
        let report = mesh.check_geometry();
        let exact_area = 200.0 * 20.0;
        let area_error = (report.total_area - exact_area).abs() / exact_area;
        println!(
            "V01 {layout:?}: {} cells, {} edges, area error {area_error:.1e}, \
             orthogonality {:.1e}, closure {:.1e}",
            mesh.n_cells(),
            mesh.n_edges(),
            report.max_orthogonality_error,
            report.max_closure_error
        );
        assert!(
            area_error <= 1e-12,
            "{layout:?}: total area {} != {exact_area}",
            report.total_area
        );
        assert!(
            report.max_orthogonality_error <= 1e-12,
            "{layout:?}: orthogonality error {:e}",
            report.max_orthogonality_error
        );
        assert!(
            report.max_closure_error <= 1e-12,
            "{layout:?}: closure error {:e}",
            report.max_closure_error
        );
    }
}

/// V02 Perot 再構成：一様流で厳密に一致すること（境界セルを含む）。
///
/// 直線水路の2つの配置に加え、格子生成ツールで作った不規則な格子でも確かめる。
/// 結果の値は `cargo test --test verification v02 -- --nocapture` で表示される。
#[test]
fn v02_perot_uniform_flow() {
    let (ux, uy) = (0.7, -0.3);
    let meshes = [
        (
            "channel, rectangular",
            channel_mesh(ChannelLayout::Rectangular, 2.0),
        ),
        (
            "channel, hexagonal",
            channel_mesh(ChannelLayout::Hexagonal, 2.0),
        ),
        ("island (tools/make_mesh.py)", island_mesh()),
    ];
    for (name, mesh) in &meshes {
        let u_n: Vec<f64> = mesh
            .normal_edge
            .iter()
            .map(|n| ux * n[0] + uy * n[1])
            .collect();
        let mut u_cell = vec![[0.0; 2]; mesh.n_cells()];
        hydro::reconstruct_cell_velocity(mesh, &u_n, &mut u_cell);
        let max_error = u_cell
            .iter()
            .map(|u| (u[0] - ux).abs().max((u[1] - uy).abs()))
            .fold(0.0, f64::max);
        println!(
            "V02 {name}: {} cells, max error {max_error:.1e} m/s",
            mesh.n_cells()
        );
        for (c, u) in u_cell.iter().enumerate() {
            assert!(
                (u[0] - ux).abs() <= 1e-12 && (u[1] - uy).abs() <= 1e-12,
                "{name}: cell {c}: reconstructed {u:?}, expected [{ux}, {uy}]"
            );
        }
    }
}

/// V03 静水保持（干出する島を含む）：流速が発生しないこと（well-balanced）。
///
/// 四辺を壁とした矩形水域の中央に、水面から突き出るガウス型の島を置く。水位一定・
/// 流速 0 から 1000 ステップ進め、max |u_n| ≤ 1e-12 m/s、湿潤セルの |η − η₀| ≤ 1e-12 m、
/// 島の上のセルが乾いたままであることを確かめる。
#[test]
#[ignore = "未実装: model::Model（半陰解法）"]
fn v03_lake_at_rest_with_island() {
    todo!("Model::with_state で島のある静水を与えて時間発展")
}

/// V04 閉鎖水域の振動：質量保存、負の水深が出ないこと。
///
/// 放物面の容器内の自由振動（Thacker, 1981 の解析解）を、wet/dry の境界が動く条件で解く。
/// 水の体積の相対誤差 ≤ 1e-10（反復解法の収束判定値に応じて設定）、全時刻で h ≥ 0。
#[test]
#[ignore = "未実装: model::Model（半陰解法、wet/dry）"]
fn v04_closed_basin_oscillation() {
    todo!("放物面容器の振動と質量収支")
}

/// V05 ダム破壊：Ritter（下流が乾燥）／ Stoker（下流が湛水）の解析解との比較。
///
/// 水位・流速分布の L1 誤差と、格子を細かくしたときの収束次数を確かめる。
#[test]
#[ignore = "未実装: model::Model（保存型移流）"]
fn v05_dam_break() {
    todo!("Ritter / Stoker の解析解との比較")
}

/// V06 等流：マニング式の等流水深に一致すること。
///
/// 直線水路（壁は滑り）で、中央部の水深が h_n = (n q / √I)^(3/5) と相対誤差 0.1% 以内。
#[test]
#[ignore = "未実装: model::Model（摩擦・流量境界）"]
fn v06_uniform_flow() {
    todo!("等流水深との比較")
}

/// V07 跳水：共役水深と跳水位置。
///
/// 射流で流入し、下流端の水位で跳水を起こす。共役水深比が Bélanger の式
/// h₂/h₁ = (√(1 + 8 Fr₁²) − 1) / 2 と 2% 以内、跳水位置が運動量の釣り合いと整合すること。
#[test]
#[ignore = "未実装: model::Model（保存型移流）"]
fn v07_hydraulic_jump() {
    todo!("共役水深と跳水位置")
}

/// V08 トレンチ移動：浮遊砂の検証。
///
/// 水路を横断するトレンチの移動速度と形状を既往の実験・計算と比較する
/// （比較対象は実装時に選定し、docs/references.md に追加する）。
#[test]
#[ignore = "未実装: sediment（浮遊砂）"]
fn v08_trench_migration() {
    todo!("トレンチの移動速度と形状")
}

/// V09 交互砂州：波長と波高。
///
/// 直線水路で交互砂州を発達させ、波長・波高を既往の実験や線形安定解析と比較する。
#[test]
#[ignore = "未実装: sediment（掃流砂・勾配補正）"]
fn v09_alternate_bars() {
    todo!("交互砂州の波長と波高")
}

/// V10 湾曲水路の河床変動：外岸側の洗掘。
///
/// 湾曲部で外岸側が洗掘され内岸側に堆積すること、その横断形状を既往の実験と比較する。
/// 二次流による流砂の向きの補正の扱いを含めて、実装時に条件を決める。
#[test]
#[ignore = "未実装: sediment、外部格子の読み込み"]
fn v10_bend_scour() {
    todo!("湾曲水路の外岸洗掘")
}

/// V11 床止め下流の局所洗掘：洗掘深と格子依存性（`examples/sill.toml`）。
///
/// 格子幅を 3 段階に変え、最大洗掘深とその位置が格子幅に対して収束することを確かめる。
#[test]
#[ignore = "未実装: model::Model（流れ・流砂・非侵食層）"]
fn v11_sill_local_scour() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/sill.toml");
    let config = Config::load(path).expect("examples/sill.toml must be valid");
    let _ = config;
    todo!("格子幅を変えた計算と洗掘深の比較")
}
