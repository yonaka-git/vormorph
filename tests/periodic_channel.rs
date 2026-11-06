//! 周期境界の継ぎ目が、計算に跡を残さないことを確かめる（ADR-0014）。
//!
//! 周期 L の水路と、同じ初期状態を2回並べた周期 2L の水路を、移動床で同じ時間刻みで進める。
//! 2L の水路の継ぎ目は x = 0 だけで、x = L は内部の辺になる。継ぎ目の扱い（体積力、Perot 法の
//! 辺の中点、河床勾配の最小二乗、安息角補正の高低差）が正しければ、2L の水路の解は L の水路の
//! 解を2回並べたものに丸め誤差の範囲で一致する。正方形配置と六角形配置の両方で確かめる。
//!
//! 周期的な水路のセルは行ごと（y の小さい行から、各行を x の順）に並び、どの行も n_x 個なので、
//! 2L の水路のセル (i, j) は L の水路のセル (i mod n_x, j) に対応する。

use std::f64::consts::PI;
use vormorph::config::{ChannelLayout, Config};
use vormorph::mesh::Mesh;
use vormorph::mesh::channel::{self, ChannelSpec};
use vormorph::model::Model;
use vormorph::state::State;

const LENGTH: f64 = 20.0;
const WIDTH: f64 = 4.0;
const SLOPE: f64 = 0.002;
const MANNING_N: f64 = 0.025;
const UNIT_DISCHARGE: f64 = 0.5;

/// 周期境界と移動床（`examples/sill.toml` の流砂の条件、最初から河床変動）の設定。
fn config(name: &str) -> Config {
    let dir = std::env::temp_dir().join(format!("vormorph-{name}-{}", std::process::id()));
    let mut config: Config = format!(
        r#"
title = "{name}"
[mesh]
kind = "channel"
length = {LENGTH:?}
width = {WIDTH:?}
cell_size = 1.0
[bed]
slope = {SLOPE:?}
[time]
dt = 0.2
end = 40.0
[hydro]
theta = 0.55
manning_n = {MANNING_N:?}
[boundary.upstream]
kind = "periodic"
discharge = {:?}
[boundary.downstream]
kind = "periodic"
[output]
dir = {dir:?}
interval = 1000.0
"#,
        UNIT_DISCHARGE * WIDTH
    )
    .parse()
    .expect("valid periodic config");
    let sill = Config::load(concat!(env!("CARGO_MANIFEST_DIR"), "/examples/sill.toml"))
        .expect("examples/sill.toml must be valid");
    config.sediment = sill.sediment;
    config
}

/// 周期 `length` の水路で、周期 L の起伏 z_b = 0.05 cos(2πx/L) cos(πy/W) + 0.02 sin(4πx/L)
/// の上に、水位一定・等流の流速を与えた状態から `steps` ステップ進める。
fn run(layout: ChannelLayout, length: f64, steps: usize) -> (Mesh, State) {
    let mesh = channel::build_periodic(&ChannelSpec {
        length,
        width: WIDTH,
        cell_size: 1.0,
        layout,
    })
    .expect("valid periodic channel");
    let h_n = (MANNING_N * UNIT_DISCHARGE / SLOPE.sqrt()).powf(0.6);
    let mut state = State::new(mesh.n_cells(), mesh.n_edges());
    for c in 0..mesh.n_cells() {
        let (x, y) = (mesh.x_cell[c], mesh.y_cell[c]);
        state.zb[c] = 0.05 * (2.0 * PI * x / LENGTH).cos() * (PI * y / WIDTH).cos()
            + 0.02 * (4.0 * PI * x / LENGTH).sin();
        state.eta[c] = h_n;
    }
    for (u, n) in state.u_n.iter_mut().zip(&mesh.normal_edge) {
        *u = UNIT_DISCHARGE / h_n * n[0];
    }
    let mut model = Model::with_state(config("periodic-seam"), mesh, state).expect("valid setup");
    for _ in 0..steps {
        model.step().expect("the step must succeed");
    }
    let water = model.water_balance();
    assert!(
        water.error().abs() <= 1e-12 * water.initial_volume,
        "{layout:?}, L = {length}: {water:?}"
    );
    let sediment = model.sediment_balance();
    assert!(
        sediment.error().abs() <= 1e-12,
        "{layout:?}, L = {length}: {sediment:?}"
    );
    (model.mesh().clone(), model.state().clone())
}

#[test]
fn doubled_periodic_channel_repeats_the_single_period_solution() {
    for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
        let steps = 200;
        let (mesh, single) = run(layout, LENGTH, steps);
        let (mesh2, double) = run(layout, 2.0 * LENGTH, steps);
        let rows = if layout == ChannelLayout::Rectangular {
            4
        } else {
            5
        };
        let nx = mesh.n_cells() / rows;
        assert_eq!(mesh2.n_cells(), 2 * mesh.n_cells(), "{layout:?}");

        let (mut max_eta, mut max_zb, mut max_u, mut bed_change) =
            (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
        for c2 in 0..mesh2.n_cells() {
            let (j, i) = (c2 / (2 * nx), c2 % (2 * nx));
            let c = j * nx + i % nx;
            let shift = if i < nx { 0.0 } else { LENGTH };
            assert!(
                (mesh2.x_cell[c2] - shift - mesh.x_cell[c]).abs() < 1e-12
                    && (mesh2.y_cell[c2] - mesh.y_cell[c]).abs() < 1e-12,
                "{layout:?}: cell {c2} does not correspond to cell {c}"
            );
            max_eta = max_eta.max((double.eta[c2] - single.eta[c]).abs());
            max_zb = max_zb.max((double.zb[c2] - single.zb[c]).abs());
            for k in 0..2 {
                max_u = max_u.max((double.u_cell[c2][k] - single.u_cell[c][k]).abs());
            }
            let (x, y) = (mesh.x_cell[c], mesh.y_cell[c]);
            let zb0 = 0.05 * (2.0 * PI * x / LENGTH).cos() * (PI * y / WIDTH).cos()
                + 0.02 * (4.0 * PI * x / LENGTH).sin();
            bed_change = bed_change.max((single.zb[c] - zb0).abs());
        }
        println!(
            "{layout:?}: bed change {bed_change:.2e} m, max difference eta {max_eta:.1e} m, \
             zb {max_zb:.1e} m, u {max_u:.1e} m/s"
        );
        // 起伏の上の流れで河床が実際に変わっていること
        assert!(bed_change > 1e-4, "{layout:?}: bed change {bed_change:e}");
        assert!(max_eta <= 1e-10, "{layout:?}: eta differs by {max_eta:e}");
        assert!(max_zb <= 1e-10, "{layout:?}: zb differs by {max_zb:e}");
        assert!(max_u <= 1e-10, "{layout:?}: velocity differs by {max_u:e}");
    }
}
