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
//!
//! 流れのケース（V03 以降）は計算量が多いので、値を記録するときはリリースビルドで実行する。
//! 環境変数 `VORMORPH_VERIFICATION_OUTPUT` にディレクトリを指定すると、ケースごとの VTK を
//! `<dir>/<ケース名>/` に書き出す（図の作成に使う。`docs/verification.md`）。
//!
//! ```text
//! VORMORPH_VERIFICATION_OUTPUT=output/verification \
//!     cargo test --release --test verification v04 -- --nocapture
//! ```

use std::f64::consts::PI;
use std::path::{Path, PathBuf};
use vormorph::config::{AdvectionScheme, ChannelLayout, Config};
use vormorph::hydro;
use vormorph::mesh::channel::{self, ChannelSpec};
use vormorph::mesh::{Mesh, NONE, file};
use vormorph::model::Model;
use vormorph::state::State;

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

/// 長方形の水域の格子（四辺は壁として扱う）。
fn basin_mesh(layout: ChannelLayout, length: f64, width: f64, cell_size: f64) -> Mesh {
    channel::build(&ChannelSpec {
        length,
        width,
        cell_size,
        layout,
    })
    .expect("the verification basin must be a valid mesh")
}

/// 流れの検証ケースの設定。格子と初期状態は [`Model::with_state`] に別に渡す。
///
/// `boundary` は `[boundary.upstream]` と `[boundary.downstream]` の TOML。出力は、
/// `VORMORPH_VERIFICATION_OUTPUT` が指定されたときだけ `output_interval` ごとに書く。
fn flow_config(
    case: &str,
    dt: f64,
    end: f64,
    manning_n: f64,
    boundary: &str,
    output_interval: f64,
) -> Config {
    let text = format!(
        r#"
title = "{case}"

[mesh]
kind = "file"
path = "unused: the mesh is given to Model::with_state"

[time]
dt = {dt:?}
end = {end:?}

[hydro]
theta = 0.55
manning_n = {manning_n:?}

{boundary}

[output]
dir = "unused"
interval = {output_interval:?}
"#
    );
    let mut config: Config = text.parse().expect("the verification config must be valid");
    match output_root() {
        Some(dir) => config.output.dir = PathBuf::from(dir).join(case),
        None => {
            // 終了時刻より長い間隔にして、何も書かないようにする
            config.output.dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(case);
            config.output.interval = 10.0 * end;
        }
    }
    config
}

/// VTK の出力先（`VORMORPH_VERIFICATION_OUTPUT`）。指定がなければ `None`。
fn output_root() -> Option<std::ffi::OsString> {
    std::env::var_os("VORMORPH_VERIFICATION_OUTPUT")
}

const CLOSED: &str = r#"
[boundary.upstream]
kind = "wall"

[boundary.downstream]
kind = "wall"
"#;

/// 状態の水の体積 V_i = A_i max(η_i − z_{b,i}, 0) \[m³\]（セルごと）。
fn cell_volumes(mesh: &Mesh, state: &State) -> Vec<f64> {
    (0..mesh.n_cells())
        .map(|c| mesh.area_cell[c] * state.depth(c))
        .collect()
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
/// 四辺を壁とした矩形水域（40 m × 20 m、格子幅 1 m）の中央に、水面から突き出るガウス型の
/// 島を置き、河床に起伏も与える。水位一定・流速 0 から 1000 ステップ進め、max |u_n| ≤ 1e-12 m/s、
/// 湿潤セルの |η − η₀| ≤ 1e-12 m、島の上のセルが乾いたまま（η = z_b）であることを確かめる。
/// 格子生成ツールの島の格子（穴あり、不規則なセル）でも、河床に島と起伏を与えて確かめる。
#[test]
fn v03_lake_at_rest_with_island() {
    let eta0 = 1.0;
    let ripple = |x: f64, y: f64| 0.2 * (0.3 * x).sin() * (0.2 * y).cos();
    let gaussian = |x: f64, y: f64, (x0, y0): (f64, f64), sigma: f64| {
        2.0 * (-((x - x0).powi(2) + (y - y0).powi(2)) / (2.0 * sigma * sigma)).exp()
    };
    let cases: [(&str, Mesh, (f64, f64), f64); 3] = [
        (
            "basin, rectangular",
            basin_mesh(ChannelLayout::Rectangular, 40.0, 20.0, 1.0),
            (20.0, 10.0),
            3.0,
        ),
        (
            "basin, hexagonal",
            basin_mesh(ChannelLayout::Hexagonal, 40.0, 20.0, 1.0),
            (20.0, 10.0),
            3.0,
        ),
        (
            "island (tools/make_mesh.py)",
            island_mesh(),
            (12.0, 15.0),
            4.0,
        ),
    ];
    for (name, mesh, center, sigma) in cases {
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        for c in 0..mesh.n_cells() {
            let (x, y) = (mesh.x_cell[c], mesh.y_cell[c]);
            state.zb[c] = gaussian(x, y, center, sigma) + ripple(x, y);
            state.eta[c] = eta0;
        }
        let dry: Vec<bool> = (0..mesh.n_cells()).map(|c| state.zb[c] >= eta0).collect();
        let n_dry = dry.iter().filter(|&&d| d).count();
        assert!(
            n_dry > 0 && n_dry < mesh.n_cells(),
            "{name}: the island must emerge"
        );

        let config = flow_config("v03", 1.0, 1000.0, 0.03, CLOSED, 100.0);
        let mut model = Model::with_state(config, mesh, state).expect("valid setup");
        let mut newton = 0;
        for _ in 0..1000 {
            newton += model
                .step()
                .expect("the step must succeed")
                .newton_iterations;
        }

        let state = model.state();
        let max_u = state.u_n.iter().fold(0.0_f64, |m, u| m.max(u.abs()));
        let max_deta = (0..dry.len())
            .filter(|&c| !dry[c])
            .fold(0.0_f64, |m, c| m.max((state.eta[c] - eta0).abs()));
        println!(
            "V03 {name}: {} cells ({n_dry} dry), max |u_n| {max_u:.1e} m/s, \
             max |eta - eta0| {max_deta:.1e} m, Newton iterations {newton}",
            dry.len()
        );
        assert!(max_u <= 1e-12, "{name}: max |u_n| = {max_u:e}");
        assert!(max_deta <= 1e-12, "{name}: max |eta - eta0| = {max_deta:e}");
        for (c, &is_dry) in dry.iter().enumerate() {
            if is_dry {
                assert_eq!(state.eta[c], state.zb[c], "{name}: cell {c} must stay dry");
            }
        }
    }
}

/// V04 の放物面容器（Thacker, 1981。SWASHES の "planar surface in a paraboloid"）の諸元。
mod thacker {
    /// 容器の半径（静水の汀線の半径）a \[m\]。
    pub const A: f64 = 1.0;
    /// 中心の静水深 h₀ \[m\]。
    pub const H0: f64 = 0.1;
    /// 水塊の中心の振れ幅 η \[m\]（水塊は半径 a の放物面のまま、中心の周りを回る）。
    pub const SHIFT: f64 = 0.5;
    /// 計算領域の一辺 \[m\]（中心は (L/2, L/2)）。
    pub const L: f64 = 4.0;
    /// 重力加速度 \[m/s²\]（設定の既定値）。
    pub const G: f64 = 9.81;

    /// 角振動数 ω = √(2 g h₀) / a \[1/s\]。
    pub fn omega() -> f64 {
        (2.0 * G * H0).sqrt() / A
    }

    /// 河床高 z_b = −h₀ (1 − r²/a²) \[m\]。
    pub fn bed(x: f64, y: f64) -> f64 {
        let (x, y) = (x - 0.5 * L, y - 0.5 * L);
        -H0 * (1.0 - (x * x + y * y) / (A * A))
    }

    /// 水面 ζ = (η h₀ / a²)(2 x cos ωt + 2 y sin ωt − η) \[m\]（水のある範囲で有効）。
    pub fn surface(x: f64, y: f64, t: f64) -> f64 {
        let (x, y) = (x - 0.5 * L, y - 0.5 * L);
        let wt = omega() * t;
        SHIFT * H0 / (A * A) * (2.0 * x * wt.cos() + 2.0 * y * wt.sin() - SHIFT)
    }

    /// 水深 h = max(ζ − z_b, 0) \[m\]。
    pub fn depth(x: f64, y: f64, t: f64) -> f64 {
        (surface(x, y, t) - bed(x, y)).max(0.0)
    }

    /// 水のある範囲の流速 (u, v) = η ω (−sin ωt, cos ωt) \[m/s\]（空間的に一様）。
    pub fn velocity(t: f64) -> [f64; 2] {
        let wt = omega() * t;
        [-SHIFT * omega() * wt.sin(), SHIFT * omega() * wt.cos()]
    }

    /// 水塊の重心 (L/2 + η cos ωt, L/2 + η sin ωt) \[m\]。
    pub fn centroid(t: f64) -> [f64; 2] {
        let wt = omega() * t;
        [0.5 * L + SHIFT * wt.cos(), 0.5 * L + SHIFT * wt.sin()]
    }
}

/// V04 の1回の計算の結果。
struct ThackerRun {
    /// 水の体積の相対誤差。
    volume_error: f64,
    /// 全ステップ・全セルの連続式の残差の最大値 \[m\]。
    max_continuity: f64,
    /// Newton 法の収束判定値 \[m\]。
    newton_tol: f64,
    /// Newton 法の反復回数の最大値。
    max_newton: usize,
    /// 1周期ごとの（水深の L1 誤差の解析解の体積に対する比、水塊の重心の誤差 \[m\]）。
    errors: Vec<(f64, f64)>,
}

/// V04 の条件で `layout` の格子を `periods` 周期解く。
///
/// VTK を書き出すときは、ステップごとの水塊の重心と水深の L1 誤差を `centroid.csv` にも書く
/// （`tools/plot_verification.py v04` の入力）。
fn thacker_run(layout: ChannelLayout, advection: AdvectionScheme, case: &str) -> ThackerRun {
    const CELL_SIZE: f64 = 0.05;
    const STEPS_PER_PERIOD: u64 = 200;
    const PERIODS: u64 = 3;
    let period = 2.0 * PI / thacker::omega();
    let dt = period / STEPS_PER_PERIOD as f64;

    let mesh = basin_mesh(layout, thacker::L, thacker::L, CELL_SIZE);
    let mut state = State::new(mesh.n_cells(), mesh.n_edges());
    for c in 0..mesh.n_cells() {
        let (x, y) = (mesh.x_cell[c], mesh.y_cell[c]);
        state.zb[c] = thacker::bed(x, y);
        state.eta[c] = thacker::surface(x, y, 0.0).max(state.zb[c]);
    }
    let [u0, v0] = thacker::velocity(0.0);
    for e in 0..mesh.n_edges() {
        let [c1, c2] = mesh.cells_on_edge[e];
        if c2 != NONE && state.depth(c1) > 0.0 && state.depth(c2) > 0.0 {
            let [nx, ny] = mesh.normal_edge[e];
            state.u_n[e] = u0 * nx + v0 * ny;
        }
    }

    let mut config = flow_config(case, dt, PERIODS as f64 * period, 0.0, CLOSED, period / 8.0);
    config.hydro.advection = advection;
    let newton_tol = config.hydro.solver.newton_tol;
    let record = output_root()
        .is_some()
        .then(|| config.output.dir.join("centroid.csv"));
    let mut model = Model::with_state(config, mesh, state).expect("valid setup");
    if record.is_some() {
        model
            .write_output()
            .expect("the initial state must be written");
    }

    let mut csv = String::from("time,centroid_x,centroid_y,l1_error\n");
    let mut result = ThackerRun {
        volume_error: 0.0,
        max_continuity: 0.0,
        newton_tol,
        max_newton: 0,
        errors: Vec::new(),
    };
    for step in 1..=PERIODS * STEPS_PER_PERIOD {
        let before = cell_volumes(model.mesh(), model.state());
        let report = model.step().expect("the step must succeed");
        result.max_newton = result.max_newton.max(report.newton_iterations);
        let (mesh, state, q) = (model.mesh(), model.state(), model.hydro().theta_flux());
        let after = cell_volumes(mesh, state);
        for c in 0..mesh.n_cells() {
            let outflow: f64 = mesh.edges_on_cell[c]
                .iter()
                .zip(&mesh.edge_sign_on_cell[c])
                .map(|(&e, &s)| s * mesh.dv_edge[e] * q[e])
                .sum();
            let r = (after[c] - before[c] + dt * outflow).abs() / mesh.area_cell[c];
            result.max_continuity = result.max_continuity.max(r);
        }

        let t = state.time;
        let (mut l1, mut total, mut volume) = (0.0, 0.0, 0.0);
        let (mut cx, mut cy) = (0.0, 0.0);
        for (c, &v) in after.iter().enumerate() {
            let (x, y) = (mesh.x_cell[c], mesh.y_cell[c]);
            let exact = thacker::depth(x, y, t);
            l1 += mesh.area_cell[c] * (state.depth(c) - exact).abs();
            total += mesh.area_cell[c] * exact;
            volume += v;
            cx += v * x;
            cy += v * y;
        }
        let (l1, cx, cy) = (l1 / total, cx / volume, cy / volume);
        csv.push_str(&format!("{t:?},{cx:?},{cy:?},{l1:?}\n"));
        if step % (STEPS_PER_PERIOD / 2) == 0 {
            let [ex, ey] = thacker::centroid(t);
            let centroid_error = (cx - ex).hypot(cy - ey);
            println!(
                "  t = {:.2} T: L1 error of h {:.2}%, centroid error {:.4} m",
                t / period,
                100.0 * l1,
                centroid_error
            );
            if step % STEPS_PER_PERIOD == 0 {
                result.errors.push((l1, centroid_error));
            }
        }
    }
    let balance = model.water_balance();
    result.volume_error = balance.error().abs() / balance.initial_volume;
    println!(
        "  volume error {:.1e}, max continuity residual {:.1e} m (tolerance {newton_tol:e}), \
         max Newton iterations {}",
        result.volume_error, result.max_continuity, result.max_newton
    );
    if let Some(path) = record {
        std::fs::write(&path, csv).expect("the centroid history must be written");
    }
    result
}

/// V04 閉鎖水域の振動：質量保存、負の水深が出ないこと、Thacker の解析解との比較。
///
/// 放物面の容器内で、平面の水面が回りながら振動する Thacker (1981) の解析解の条件で、
/// 汀線が動く流れを3周期解く（摩擦なし、格子幅 0.05 m、時間刻み T/200）。次を確かめる。
///
/// - 水の体積の相対誤差 ≤ 1e-10
/// - 全ステップで、乾いたセルを含む全セルの連続式の残差 |V_i^{n+1} − V_i^n + Δt Σ_j s L q^θ| / A_i
///   が Newton 法の収束判定値以下（負の体積を避けるために水を作り出していない）
/// - 水塊の重心の誤差が 1周期後に 0.12 m 以下、3周期後に 0.25 m 以下
/// - 水深の L1 誤差（解析解の体積に対する比）が 1周期後に 20% 以下、3周期後に 40% 以下
///
/// 解析解との基準は、移流がない場合（1周期後に重心の誤差 0.27 m、L1 誤差 46%）を検出できる
/// 水準で、1次精度の誤差（格子幅を半分にすると約 0.6 倍）に余裕を持たせて決めた
/// （`docs/verification.md`）。VTK を書き出すときは、図のために移流なしの計算も行う。
#[test]
fn v04_closed_basin_oscillation() {
    for layout in LAYOUTS {
        let case = format!("v04_{}", format!("{layout:?}").to_lowercase());
        println!("V04 {layout:?}:");
        let run = thacker_run(layout, AdvectionScheme::Upwind, &case);
        assert!(
            run.volume_error <= 1e-10,
            "{layout:?}: volume error {:e}",
            run.volume_error
        );
        assert!(
            run.max_continuity <= run.newton_tol,
            "{layout:?}: continuity residual {:e} m",
            run.max_continuity
        );
        let limits = [(1, 0.20, 0.12), (3, 0.40, 0.25)];
        for (periods, max_l1, max_centroid) in limits {
            let (l1, centroid) = run.errors[periods - 1];
            assert!(
                l1 <= max_l1,
                "{layout:?}: L1 error {l1} after {periods} periods"
            );
            assert!(
                centroid <= max_centroid,
                "{layout:?}: centroid error {centroid} m after {periods} periods"
            );
        }
    }
    if output_root().is_some() {
        println!("V04 Rectangular, without advection (for the figure):");
        thacker_run(
            ChannelLayout::Rectangular,
            AdvectionScheme::Disabled,
            "v04_rectangular_none",
        );
    }
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
/// 直線水路（延長 200 m、幅 20 m、勾配 1/500、側壁は滑り）に流量 40 m³/s を与え、下流端の
/// 水位を等流水位に固定する（`examples/sill.toml` と同じ流量・勾配・粗度）。一様な水深・
/// 流速 0 から 1200 s 進め、中央部（50 m ≤ x ≤ 150 m）の水深と h_n = (n q / √I)^(3/5) の
/// 相対誤差の最大値が 1e-4 以下であることを確かめる。格子は 2 m の正方形配置・六角形配置と、
/// 上下流端に短い境界辺ができる 1 m の六角形配置（sill と同じ）。時間刻みは、辺の水深を
/// 風上側からとることによるクーラン数の制約 U Δt / Δx ≤ 1 に対して、約 0.47 になるようにする。
#[test]
fn v06_uniform_flow() {
    const LENGTH: f64 = 200.0;
    const WIDTH: f64 = 20.0;
    const SLOPE: f64 = 0.002;
    const DISCHARGE: f64 = 40.0;
    const MANNING_N: f64 = 0.025;
    const END: f64 = 1200.0;
    let q = DISCHARGE / WIDTH;
    let h_n = (MANNING_N * q / SLOPE.sqrt()).powf(0.6);
    let boundary = format!(
        r#"
[boundary.upstream]
kind = "discharge"
discharge = {DISCHARGE:?}

[boundary.downstream]
kind = "water_level"
level = {h_n:?}
"#
    );

    let cases = [
        ("rectangular, 2 m", ChannelLayout::Rectangular, 2.0),
        ("hexagonal, 2 m", ChannelLayout::Hexagonal, 2.0),
        ("hexagonal, 1 m", ChannelLayout::Hexagonal, 1.0),
    ];
    for (name, layout, cell_size) in cases {
        let mesh = basin_mesh(layout, LENGTH, WIDTH, cell_size);
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        for c in 0..mesh.n_cells() {
            state.zb[c] = SLOPE * (LENGTH - mesh.x_cell[c]);
            state.eta[c] = state.zb[c] + h_n;
        }
        let case = format!("v06_{}_{cell_size}m", format!("{layout:?}").to_lowercase());
        let dt = 0.25 * cell_size;
        let config = flow_config(&case, dt, END, MANNING_N, &boundary, 300.0);
        let mut model = Model::with_state(config, mesh, state).expect("valid setup");
        if output_root().is_some() {
            model
                .write_output()
                .expect("the initial state must be written");
        }
        let mut pcg = 0;
        while model.state().time < END {
            let report = model.step().expect("the step must succeed");
            pcg = pcg.max(report.pcg_iterations);
        }

        let (mesh, state) = (model.mesh(), model.state());
        let relative_error = |c: usize| (state.depth(c) - h_n) / h_n;
        let central = (0..mesh.n_cells()).filter(|&c| (50.0..=150.0).contains(&mesh.x_cell[c]));
        let max_error = central.fold(0.0_f64, |m, c| m.max(relative_error(c).abs()));
        let overall = (0..mesh.n_cells()).fold(0.0_f64, |m, c| m.max(relative_error(c).abs()));
        let max_v = state.u_cell.iter().fold(0.0_f64, |m, u| m.max(u[1].abs()));
        let balance = model.water_balance();
        println!(
            "V06 {name}: {} cells, dt {dt} s, h_n = {h_n:.6} m, max relative error {max_error:.1e} \
             (central), {overall:.1e} (all cells), max |v| {max_v:.1e} m/s, \
             water balance error {:.1e} of inflow, max PCG iterations {pcg}",
            mesh.n_cells(),
            balance.error() / balance.inflow
        );
        assert!(
            max_error <= 1e-4,
            "{name}: relative depth error {max_error:e}"
        );
    }
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
