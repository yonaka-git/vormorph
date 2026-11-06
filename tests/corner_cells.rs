//! 凸な角に接するセルの掃流砂（`docs/design/mesh-generation.md` の第8節）。
//!
//! 格子生成ツールの格子では、外周の凸な角（90°）に接するセルの母点が重心から約 0.2h ずれ、
//! 流速が流れの向きに変わる流れ（∂u/∂x ≠ 0）で Perot 法の誤差が内部の 6〜8 倍になる。
//! この誤差が掃流砂を通して角のセルの河床に偏りを生まないかを確かめる。
//!
//! 格子生成ツールの直線水路（延長 60 m、幅 10 m、格子幅 2 m、`tests/data/straight_coarse.mesh.toml`）
//! で、下流端の水位を等流水位より下げ、下流端へ向かって加速する流れ（低下背水）を作る。流れは
//! 1次元なので、河床の変化は x だけの関数のはず。同じ条件を正方形配置の直線水路（格子幅 0.5 m、
//! 横断方向の流速は丸め誤差で 0）で解いた縦断形を基準にし、各セルの河床の変化との差を、
//! 内部・直線の壁ぎわ・凸な角のセルに分けて比べる。
//!
//! 結果（2026-11-06）：河床の変化の最大値 0.111 m に対し、基準との差の最大値は、内部 1.27 mm、
//! 直線の壁ぎわ 1.35 mm、凸な角 1.37 mm で、角のセルに偏りはなかった（差は主に格子幅の違い）。

use std::path::Path;
use vormorph::config::{ChannelLayout, Config, SedimentConfig};
use vormorph::mesh::channel::{self, ChannelSpec};
use vormorph::mesh::{EdgeKind, Mesh, NONE, file};
use vormorph::model::Model;
use vormorph::state::State;

const LENGTH: f64 = 60.0;
const WIDTH: f64 = 10.0;
const SLOPE: f64 = 0.002;
const MANNING_N: f64 = 0.025;
const UNIT_DISCHARGE: f64 = 2.0;
/// 下流端の水位を等流水位から下げる量 \[m\]。
const DRAWDOWN: f64 = 0.3;
/// 固定床の助走の時間と、河床変動を計算する時間（流れの時間）\[s\]。
const SPIN_UP: f64 = 300.0;
const DURATION: f64 = 300.0;

fn sill_sediment() -> SedimentConfig {
    Config::load(concat!(env!("CARGO_MANIFEST_DIR"), "/examples/sill.toml"))
        .expect("examples/sill.toml must be valid")
        .sediment
        .expect("examples/sill.toml has [sediment]")
}

/// 低下背水の流れで河床変動を解き、各セルの（x, y, 河床の変化）を返す。
fn drawdown_run(mesh: Mesh) -> (Mesh, Vec<f64>) {
    let h_n = (MANNING_N * UNIT_DISCHARGE / SLOPE.sqrt()).powf(0.6);
    let text = format!(
        r#"
title = "corner"
[mesh]
kind = "file"
path = "unused"
[time]
courant = 0.5
dt_max = 0.5
end = {:?}
morphology_start = {SPIN_UP:?}
[hydro]
theta = 0.55
manning_n = {MANNING_N:?}
[boundary.upstream]
kind = "discharge"
discharge = {:?}
[boundary.downstream]
kind = "water_level"
level = {:?}
[output]
dir = {:?}
interval = 1.0e9
"#,
        SPIN_UP + DURATION,
        UNIT_DISCHARGE * WIDTH,
        h_n - DRAWDOWN,
        Path::new(env!("CARGO_TARGET_TMPDIR")).join("corner_cells")
    );
    let mut config: Config = text.parse().expect("valid config");
    config.sediment = Some(sill_sediment());
    let mut state = State::new(mesh.n_cells(), mesh.n_edges());
    for c in 0..mesh.n_cells() {
        state.zb[c] = SLOPE * (LENGTH - mesh.x_cell[c]);
        state.eta[c] = state.zb[c] + h_n;
    }
    for e in 0..mesh.n_edges() {
        state.u_n[e] = UNIT_DISCHARGE / h_n * mesh.normal_edge[e][0];
    }
    let zb0 = state.zb.clone();
    let mut model = Model::with_state(config, mesh, state).expect("valid setup");
    while !model.finished() {
        model.step().expect("the step must succeed");
    }
    let change = model
        .state()
        .zb
        .iter()
        .zip(&zb0)
        .map(|(z, z0)| z - z0)
        .collect();
    (model.mesh().clone(), change)
}

/// セルの分類：0 = 内部、1 = 直線の壁ぎわ、2 = 凸な角に接する（外周の角の点を頂点に持つ）、
/// 3 = それ以外の上下流端に接するセル。
fn category(mesh: &Mesh, c: usize) -> usize {
    let at_corner = mesh.vertices_on_cell[c].iter().any(|&v| {
        let (x, y) = (mesh.x_vertex[v], mesh.y_vertex[v]);
        let near = |a: f64, b: f64| (a - b).abs() < 1e-9;
        (near(x, 0.0) || near(x, LENGTH)) && (near(y, 0.0) || near(y, WIDTH))
    });
    if at_corner {
        return 2;
    }
    let kinds: Vec<EdgeKind> = mesh.edges_on_cell[c]
        .iter()
        .zip(&mesh.cells_on_cell[c])
        .filter(|&(_, &k)| k == NONE)
        .map(|(&e, _)| mesh.edge_kind[e])
        .collect();
    if kinds.contains(&EdgeKind::Wall) {
        1
    } else if kinds.is_empty() {
        0
    } else {
        3
    }
}

#[test]
fn convex_corner_cells_follow_the_cross_section_in_accelerating_flow() {
    // 1次元の基準（正方形配置、格子幅 0.5 m）
    let reference_mesh = channel::build(&ChannelSpec {
        length: LENGTH,
        width: WIDTH,
        cell_size: 0.5,
        layout: ChannelLayout::Rectangular,
    })
    .expect("valid channel");
    let (reference_mesh, reference_change) = drawdown_run(reference_mesh);
    let mut reference: Vec<(f64, f64)> = (0..reference_mesh.n_cells())
        .filter(|&c| reference_mesh.y_cell[c] < 0.5)
        .map(|c| (reference_mesh.x_cell[c], reference_change[c]))
        .collect();
    reference.sort_by(|a, b| a.0.total_cmp(&b.0));
    let interpolate = |x: f64| {
        let k = reference
            .partition_point(|p| p.0 < x)
            .clamp(1, reference.len() - 1);
        let ((x0, z0), (x1, z1)) = (reference[k - 1], reference[k]);
        z0 + (z1 - z0) * ((x - x0) / (x1 - x0)).clamp(0.0, 1.0)
    };

    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/straight_coarse.mesh.toml");
    let (mesh, change) = drawdown_run(file::read(&path).expect("valid fixture"));
    let names = [
        "interior",
        "straight wall",
        "convex corner",
        "inlet or outlet",
    ];
    let mut worst = [(0.0_f64, 0usize); 4];
    for (c, dz) in change.iter().enumerate() {
        let k = category(&mesh, c);
        let deviation = (dz - interpolate(mesh.x_cell[c])).abs();
        if deviation > worst[k].0 {
            worst[k] = (deviation, c);
        }
    }
    let scale = reference.iter().fold(0.0_f64, |m, p| m.max(p.1.abs()));
    for (name, (deviation, c)) in names.iter().zip(worst) {
        println!(
            "{name}: max |dz_b - 1D reference| {deviation:.2e} m at ({:.2}, {:.2}) \
             (largest bed change in the reference {scale:.3} m)",
            mesh.x_cell[c], mesh.y_cell[c]
        );
    }
    // 角のセルの差は、内部・直線の壁ぎわのセルの差（主に格子幅 2 m と 0.5 m の離散化の差）の
    // 2 倍以下
    let others = worst[0].0.max(worst[1].0);
    assert!(worst[2].0 > 0.0, "the mesh must have convex corner cells");
    assert!(
        worst[2].0 <= 2.0 * others,
        "convex corner cells deviate by {:e} m, others by {others:e} m",
        worst[2].0
    );
    assert!(
        others <= 0.03 * scale,
        "the 2D mesh must follow the 1D reference"
    );
}
