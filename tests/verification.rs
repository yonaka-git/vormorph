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
use vormorph::mesh::{EdgeKind, Mesh, NONE, file};
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

/// V05 のダム破壊（SWASHES の "dam break on a wet / dry domain without friction"）の諸元と解析解。
mod dam_break {
    /// 水路の延長 \[m\]。
    pub const L: f64 = 10.0;
    /// 水路の幅 \[m\]（1次元の流れを、側壁が滑り条件の細い水路で解く）。
    pub const WIDTH: f64 = 0.2;
    /// ダムの位置 \[m\]。
    pub const X0: f64 = 5.0;
    /// 上流側の水深 \[m\]。
    pub const HL: f64 = 0.005;
    /// 比べる時刻 \[s\]。
    pub const T: f64 = 6.0;
    /// 重力加速度 \[m/s²\]（設定の既定値）。
    pub const G: f64 = 9.81;

    /// Stoker の解の中間の状態の波速 c_m = √(g h_m)。下流側の水深 `hr` > 0 で、
    /// −8 g h_r c_m² (√(g h_l) − c_m)² + (c_m² − g h_r)² (c_m² + g h_r) = 0 を二分法で解く。
    pub fn stoker_cm(hr: f64) -> f64 {
        let c0 = (G * HL).sqrt();
        let f = |c: f64| {
            -8.0 * G * hr * c * c * (c0 - c).powi(2) + (c * c - G * hr).powi(2) * (c * c + G * hr)
        };
        let (mut lo, mut hi) = ((G * hr).sqrt(), c0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if f(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// 衝撃波の速さ 2 c_m² (√(g h_l) − c_m) / (c_m² − g h_r) \[m/s\]（Stoker）。
    pub fn shock_speed(hr: f64) -> f64 {
        let (c0, cm) = ((G * HL).sqrt(), stoker_cm(hr));
        2.0 * cm * cm * (c0 - cm) / (cm * cm - G * hr)
    }

    /// 時刻 `t` > 0 の水深 \[m\]。`hr` = 0 なら Ritter、`hr` > 0 なら Stoker の解。
    pub fn depth(x: f64, t: f64, hr: f64) -> f64 {
        let c0 = (G * HL).sqrt();
        let xi = (x - X0) / t;
        let rarefaction = 4.0 / (9.0 * G) * (c0 - 0.5 * xi).powi(2);
        if xi <= -c0 {
            HL
        } else if hr == 0.0 {
            if xi <= 2.0 * c0 { rarefaction } else { 0.0 }
        } else {
            let cm = stoker_cm(hr);
            if xi <= 2.0 * c0 - 3.0 * cm {
                rarefaction
            } else if xi <= shock_speed(hr) {
                cm * cm / G
            } else {
                hr
            }
        }
    }
}

/// V05 の1回の計算の結果。
struct DamBreakRun {
    /// 水深の L1 誤差の、解析解の体積に対する比。
    l1_error: f64,
    /// Stoker：質量から求めた衝撃波の位置 \[m\]。Ritter：水深が 1e-4 m 以上の最も下流のセルの位置 \[m\]。
    front: f64,
    /// 横断方向の流速の最大値 \[m/s\]。
    max_v: f64,
    /// 局所クーラン数の最大値。
    max_courant: f64,
}

/// V05 の条件で、下流側の水深 `hr` のダム破壊を時刻 T まで解く。
fn dam_break_run(layout: ChannelLayout, cell_size: f64, hr: f64, case: &str) -> DamBreakRun {
    use dam_break::{G, HL, L, T, WIDTH, X0};
    // 時間刻みは、最大の流速 2√(g h_l)（Ritter の先端）に対するクーラン数が 0.5 以下で、
    // T をちょうど割り切る値
    let steps = (T * 2.0 * (G * HL).sqrt() / (0.5 * cell_size)).ceil();
    let dt = T / steps;
    let mesh = basin_mesh(layout, L, WIDTH, cell_size);
    let mut state = State::new(mesh.n_cells(), mesh.n_edges());
    for c in 0..mesh.n_cells() {
        state.eta[c] = if mesh.x_cell[c] < X0 { HL } else { hr };
    }
    let mut config = flow_config(case, dt, T, 0.0, CLOSED, T);
    // 水深は mm の単位なので、wet/dry の閾値を水深より十分小さくする
    config.hydro.dry_depth = 1.0e-6;
    let mut model = Model::with_state(config, mesh, state).expect("valid setup");
    if output_root().is_some() {
        model
            .write_output()
            .expect("the initial state must be written");
    }
    let mut max_courant = 0.0_f64;
    for _ in 0..steps as u64 {
        let report = model.step().expect("the step must succeed");
        max_courant = max_courant.max(report.max_courant);
    }
    let (mesh, state) = (model.mesh(), model.state());
    assert!((state.time - T).abs() <= 1e-9 * T, "time {}", state.time);

    let (mut l1, mut total) = (0.0, 0.0);
    for c in 0..mesh.n_cells() {
        let exact = dam_break::depth(mesh.x_cell[c], T, hr);
        l1 += mesh.area_cell[c] * (state.depth(c) - exact).abs();
        total += mesh.area_cell[c] * exact;
    }
    let front = if hr > 0.0 {
        // 中間の一様な区間の中ほどの x_m より下流の水の量 ∫ (h − h_r) dx は、解析解では
        // (h_m − h_r)(x_C − x_m)。これから衝撃波の位置 x_C を求める
        let (c0, cm) = ((G * HL).sqrt(), dam_break::stoker_cm(hr));
        let x_b = X0 + T * (2.0 * c0 - 3.0 * cm);
        let x_c = X0 + T * dam_break::shock_speed(hr);
        let x_m = 0.5 * (x_b + x_c);
        let excess: f64 = (0..mesh.n_cells())
            .filter(|&c| mesh.x_cell[c] > x_m)
            .map(|c| mesh.area_cell[c] * (state.depth(c) - hr))
            .sum();
        x_m + excess / (WIDTH * (cm * cm / G - hr))
    } else {
        (0..mesh.n_cells())
            .filter(|&c| state.depth(c) >= 1.0e-4)
            .map(|c| mesh.x_cell[c])
            .fold(f64::NEG_INFINITY, f64::max)
    };
    let max_v = state.u_cell.iter().fold(0.0_f64, |m, u| m.max(u[1].abs()));
    DamBreakRun {
        l1_error: l1 / total,
        front,
        max_v,
        max_courant,
    }
}

/// V05 ダム破壊：Stoker（下流が湛水）／ Ritter（下流が乾燥）の解析解との比較。
///
/// SWASHES（Delestre et al., 2013）の §4.1.1, §4.1.2 と同じ条件（延長 10 m、x₀ = 5 m、
/// h_l = 5 mm、h_r = 1 mm／0、摩擦なし、t = 6 s）を、幅 0.2 m の水路で解く。格子幅を
/// 0.1, 0.05, 0.025, 0.0125 m と変え、両方の配置で次を確かめる。
///
/// - 水深の L1 誤差（解析解の体積に対する比）が、格子を細かくするたびに減り、0.1 m から
///   0.0125 m までの平均の収束次数が 0.5 以上（不連続を含む解での1次精度のスキームの下限）
/// - 0.0125 m での L1 誤差が、Stoker で 0.5% 以下、Ritter で 0.7% 以下
/// - Stoker：質量から求めた衝撃波の位置の誤差が、どの格子でも格子幅の 2 倍以下
/// - Ritter：先端（水深 1e-4 m の位置）の遅れが、格子を細かくするたびに減る
#[test]
fn v05_dam_break() {
    use dam_break::{G, HL, T, X0};
    const CELL_SIZES: [f64; 4] = [0.1, 0.05, 0.025, 0.0125];
    let mut csv = String::from("case,layout,cell_size,l1_error,front\n");
    for (name, hr, max_l1) in [("stoker", 0.001, 0.005), ("ritter", 0.0, 0.007)] {
        let exact_front = if hr > 0.0 {
            X0 + T * dam_break::shock_speed(hr)
        } else {
            // 解析解の水深が 1e-4 m になる位置
            let c0 = (G * HL).sqrt();
            X0 + T * 2.0 * (c0 - (1.0e-4 * 9.0 * G / 4.0).sqrt())
        };
        for layout in LAYOUTS {
            let layout_name = format!("{layout:?}").to_lowercase();
            println!("V05 {name}, {layout:?}: front (exact) {exact_front:.3} m");
            let mut runs = Vec::new();
            for cell_size in CELL_SIZES {
                let case = format!("v05_{name}_{layout_name}_{cell_size}m");
                let run = dam_break_run(layout, cell_size, hr, &case);
                let order = runs
                    .last()
                    .map(|p: &DamBreakRun| (p.l1_error / run.l1_error).log2());
                println!(
                    "  dx {cell_size} m: L1 error {:.2}%, order {}, front {:.3} m \
                     (error {:+.3} m = {:+.1} dx), max |v| {:.1e} m/s, max Courant {:.2}",
                    100.0 * run.l1_error,
                    order.map_or("-".to_owned(), |p| format!("{p:.2}")),
                    run.front,
                    run.front - exact_front,
                    (run.front - exact_front) / cell_size,
                    run.max_v,
                    run.max_courant
                );
                csv.push_str(&format!(
                    "{name},{layout_name},{cell_size:?},{:?},{:?}\n",
                    run.l1_error, run.front
                ));
                if hr > 0.0 {
                    assert!(
                        (run.front - exact_front).abs() <= 2.0 * cell_size,
                        "{name}, {layout:?}, dx {cell_size}: shock at {} m, exact {exact_front} m",
                        run.front
                    );
                }
                runs.push(run);
            }
            for (pair, sizes) in runs.windows(2).zip(CELL_SIZES.windows(2)) {
                assert!(
                    pair[1].l1_error < pair[0].l1_error,
                    "{name}, {layout:?}: L1 error grows from dx {} to {}",
                    sizes[0],
                    sizes[1]
                );
                if hr == 0.0 {
                    assert!(
                        pair[1].front > pair[0].front,
                        "{name}, {layout:?}: the front lag grows from dx {} to {}",
                        sizes[0],
                        sizes[1]
                    );
                }
            }
            let (first, last) = (&runs[0], &runs[runs.len() - 1]);
            let mean_order =
                (first.l1_error / last.l1_error).log2() / (CELL_SIZES.len() - 1) as f64;
            println!("  mean order {mean_order:.2}");
            assert!(
                mean_order >= 0.5,
                "{name}, {layout:?}: mean order {mean_order}"
            );
            assert!(
                last.l1_error <= max_l1,
                "{name}, {layout:?}: L1 error {} on the finest grid",
                last.l1_error
            );
        }
    }
    if let Some(root) = output_root() {
        let path = PathBuf::from(root).join("v05_convergence.csv");
        std::fs::write(path, csv).expect("the convergence table must be written");
    }
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

/// V07 のこぶを越える跳水（SWASHES の "transcritical flow with shock"）の諸元と解析解。
mod bump_jump {
    /// 水路の延長 \[m\]。
    pub const L: f64 = 25.0;
    /// 水路の幅 \[m\]（1次元の流れを、側壁が滑り条件の細い水路で解く）。
    pub const WIDTH: f64 = 0.2;
    /// 単位幅流量 \[m²/s\]。
    pub const Q: f64 = 0.18;
    /// 下流端の水深（= 水位。下流端の河床高は 0）\[m\]。
    pub const H_DOWNSTREAM: f64 = 0.33;
    /// こぶの頂部の高さ \[m\]。
    pub const Z_MAX: f64 = 0.2;
    /// 重力加速度 \[m/s²\]（設定の既定値）。
    pub const G: f64 = 9.81;

    /// 河床高 z(x) = 0.2 − 0.05 (x − 10)²（8 m < x < 12 m）、それ以外は 0 \[m\]。
    pub fn bed(x: f64) -> f64 {
        if 8.0 < x && x < 12.0 {
            Z_MAX - 0.05 * (x - 10.0).powi(2)
        } else {
            0.0
        }
    }

    /// 限界水深 h_c = (q² / g)^(1/3) \[m\]。
    pub fn critical_depth() -> f64 {
        (Q * Q / G).cbrt()
    }

    /// 比エネルギー E(h) = h + q² / (2 g h²) \[m\]。
    fn specific_energy(h: f64) -> f64 {
        h + Q * Q / (2.0 * G * h * h)
    }

    /// 全水頭 `head` の流れの、河床高 `z` での水深 \[m\]。`supercritical` なら射流側
    /// （h < h_c）、そうでなければ常流側の根。比エネルギーが足りなければ `None`。
    pub fn depth_for_head(head: f64, z: f64, supercritical: bool) -> Option<f64> {
        let (hc, e) = (critical_depth(), head - z);
        if e < specific_energy(hc) {
            return None;
        }
        // 射流側では E は h について減少、常流側では増加する
        let (mut lo, mut hi) = if supercritical { (1e-6, hc) } else { (hc, e) };
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            let above = specific_energy(mid) > e;
            if above == supercritical {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(0.5 * (lo + hi))
    }

    /// 上流側の全水頭（こぶの頂部で限界流）\[m\]。
    pub fn upstream_head() -> f64 {
        Z_MAX + 1.5 * critical_depth()
    }

    /// 下流側の全水頭 \[m\]。
    pub fn downstream_head() -> f64 {
        H_DOWNSTREAM + Q * Q / (2.0 * G * H_DOWNSTREAM * H_DOWNSTREAM)
    }

    /// 跳水の位置 x_s \[m\]：射流側の水深 h₁ と下流側の常流の水深 h₂ が
    /// q² (1/h₁ − 1/h₂) + g (h₁² − h₂²) / 2 = 0（Bélanger の式）を満たす位置。
    pub fn jump_position() -> f64 {
        let balance = |x: f64| {
            let z = bed(x);
            let h1 = depth_for_head(upstream_head(), z, true)?;
            let h2 = depth_for_head(downstream_head(), z, false)?;
            Some(Q * Q * (1.0 / h1 - 1.0 / h2) + 0.5 * G * (h1 * h1 - h2 * h2))
        };
        // 下流側の常流の根がある範囲で、頂部側（運動量が射流の方が大きい）から二分法
        let (mut lo, mut hi) = (10.0, 12.0);
        while balance(lo).is_none() {
            lo += 1e-3;
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if balance(mid).is_some_and(|b| b > 0.0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// 跳水の位置での共役水深 (h₁, h₂) \[m\]。
    pub fn conjugate_depths() -> (f64, f64) {
        let z = bed(jump_position());
        let h1 = depth_for_head(upstream_head(), z, true).expect("h1 exists at the jump");
        let h2 = depth_for_head(downstream_head(), z, false).expect("h2 exists at the jump");
        (h1, h2)
    }

    /// 解析解の水深 \[m\]。跳水より上流は上流側の全水頭（頂部より上流は常流、下流は射流）、
    /// 跳水より下流は下流側の全水頭の常流。
    pub fn depth(x: f64) -> f64 {
        let z = bed(x);
        let h = if x < jump_position() {
            depth_for_head(upstream_head(), z, x > 10.0)
        } else {
            depth_for_head(downstream_head(), z, false)
        };
        h.expect("the analytic solution exists everywhere")
    }
}

/// V07 の1回の計算の結果。
struct JumpRun {
    /// 跳水の位置（x 方向に平均した水深が (h₁ + h₂)/2 を上回る位置）\[m\]。
    jump: f64,
    /// 跳水の直前（解析解の跳水の位置から 0.3 m 上流）の水深の相対誤差。
    before: f64,
    /// 跳水の直後（同じく 0.3 m 下流）の水深の相対誤差。
    after: f64,
    /// 跳水の位置（解析解・計算値）から格子幅の 3 倍以上離れた範囲の水深の相対誤差の最大値と、
    /// その位置 \[m\]。
    max_error: (f64, f64),
    /// 下流端の流出量と上流端の流入量の比 − 1。
    discharge_error: f64,
    /// 横断方向の流速の最大値 \[m/s\]。
    max_v: f64,
}

/// V07 の条件で、`end` 秒まで解く。
fn bump_jump_run(layout: ChannelLayout, cell_size: f64, end: f64, case: &str) -> JumpRun {
    use bump_jump::{H_DOWNSTREAM, L, Q, WIDTH};
    let boundary = format!(
        r#"
[boundary.upstream]
kind = "discharge"
discharge = {:?}

[boundary.downstream]
kind = "water_level"
level = {H_DOWNSTREAM:?}
"#,
        Q * WIDTH
    );
    // 跳水の直前の流速 q/h₁ ≈ 2.4 m/s に対するクーラン数が 0.5 以下
    let dt = 0.5 * cell_size / 2.5;
    let mesh = basin_mesh(layout, L, WIDTH, cell_size);
    let mut state = State::new(mesh.n_cells(), mesh.n_edges());
    for c in 0..mesh.n_cells() {
        state.zb[c] = bump_jump::bed(mesh.x_cell[c]);
        state.eta[c] = H_DOWNSTREAM;
    }
    let config = flow_config(case, dt, end, 0.0, &boundary, end / 4.0);
    let mut model = Model::with_state(config, mesh, state).expect("valid setup");
    if output_root().is_some() {
        model
            .write_output()
            .expect("the initial state must be written");
    }
    while model.state().time < end - 1e-9 * dt {
        model.step().expect("the step must succeed");
    }
    let (mesh, state) = (model.mesh(), model.state());

    // x 方向に格子幅ごとのビンで平均した水深
    let bins = (L / cell_size).round() as usize;
    let (mut x_sum, mut h_sum, mut count) = (vec![0.0; bins], vec![0.0; bins], vec![0.0; bins]);
    for c in 0..mesh.n_cells() {
        let k = ((mesh.x_cell[c] / cell_size) as usize).min(bins - 1);
        x_sum[k] += mesh.x_cell[c];
        h_sum[k] += state.depth(c);
        count[k] += 1.0;
    }
    let profile: Vec<(f64, f64)> = (0..bins)
        .filter(|&k| count[k] > 0.0)
        .map(|k| (x_sum[k] / count[k], h_sum[k] / count[k]))
        .collect();
    let error = |x: f64, h: f64| h / bump_jump::depth(x) - 1.0;
    let error_near = |x0: f64| {
        let &(x, h) = profile
            .iter()
            .min_by(|a, b| (a.0 - x0).abs().total_cmp(&(b.0 - x0).abs()))
            .expect("the profile is not empty");
        error(x, h)
    };

    let x_s = bump_jump::jump_position();
    let (h1, h2) = bump_jump::conjugate_depths();
    let threshold = 0.5 * (h1 + h2);
    let jump = profile
        .windows(2)
        .find(|w| w[0].0 > 10.5 && w[0].1 < threshold && w[1].1 >= threshold)
        .map(|w| w[0].0 + (threshold - w[0].1) / (w[1].1 - w[0].1) * (w[1].0 - w[0].0))
        .expect("the jump must be found downstream of the crest");
    let max_error = profile
        .iter()
        .filter(|(x, _)| (x - x_s).abs() >= 3.0 * cell_size && (x - jump).abs() >= 3.0 * cell_size)
        .map(|&(x, h)| (error(x, h).abs(), x))
        .fold((0.0, f64::NAN), |m, e| if e.0 > m.0 { e } else { m });
    let (inflow, outflow) = mesh
        .edge_kind
        .iter()
        .zip(model.hydro().theta_flux())
        .zip(&mesh.dv_edge)
        .fold((0.0, 0.0), |(i, o), ((&kind, &q), &l)| match kind {
            EdgeKind::Upstream => (i - l * q, o),
            EdgeKind::Downstream => (i, o + l * q),
            EdgeKind::Interior | EdgeKind::Wall => (i, o),
        });
    JumpRun {
        jump,
        before: error_near(x_s - 0.3),
        after: error_near(x_s + 0.3),
        max_error,
        discharge_error: outflow / inflow - 1.0,
        max_v: state.u_cell.iter().fold(0.0_f64, |m, u| m.max(u[1].abs())),
    }
}

/// V07 跳水：こぶを越える流れの跳水の位置と、跳水の前後の水深。
///
/// SWASHES（Delestre et al., 2013）の §3.1.5 と同じ条件（延長 25 m、こぶの高さ 0.2 m、
/// q = 0.18 m²/s、下流端の水深 0.33 m、摩擦なし）で、こぶの頂部で常流から射流に変わり、
/// 下流側で跳水して常流に戻る流れを、幅 0.2 m の水路で 300 s 解く。格子は正方形配置の
/// 0.1 m と 0.05 m。次を確かめる。
///
/// - 定常になっている：下流端の流出量と上流端の流入量の差が 0.1% 以下
/// - 跳水の位置（水深が共役水深の平均を上回る位置）の誤差が格子幅の 2 倍以下
/// - 跳水の直前と直後（解析解の跳水の位置から 0.3 m 上流・下流）の水深の誤差が 2% 以下
///   （解析解の共役水深の比 h₂/h₁ = 3.41 は Bélanger の式を満たす）
/// - 跳水の近く（格子幅の 3 倍以内）を除く水深の誤差の最大値が 4% 以下で、格子を細かくすると
///   減る
///
/// 六角形配置では、跳水の位置で側壁沿いに横断方向の流れが生じて定常にならないので、
/// 合格基準には含めない（`docs/verification.md`）。
#[test]
fn v07_hydraulic_jump() {
    let x_s = bump_jump::jump_position();
    let (h1, h2) = bump_jump::conjugate_depths();
    println!(
        "V07: jump at {x_s:.4} m, conjugate depths {h1:.4} m / {h2:.4} m (ratio {:.3}), \
         h_c = {:.4} m",
        h2 / h1,
        bump_jump::critical_depth()
    );
    let mut previous_max = f64::INFINITY;
    for cell_size in [0.1, 0.05] {
        let layout = ChannelLayout::Rectangular;
        let case = format!("v07_{}_{cell_size}m", format!("{layout:?}").to_lowercase());
        let run = bump_jump_run(layout, cell_size, 300.0, &case);
        println!(
            "V07 {layout:?}, {cell_size} m: jump at {:.3} m (error {:+.3} m = {:+.1} dx), \
             depth error {:+.2}% before / {:+.2}% after the jump, \
             max {:.2}% at x = {:.2} m, outflow / inflow − 1 = {:.1e}, max |v| {:.1e} m/s",
            run.jump,
            run.jump - x_s,
            (run.jump - x_s) / cell_size,
            100.0 * run.before,
            100.0 * run.after,
            100.0 * run.max_error.0,
            run.max_error.1,
            run.discharge_error,
            run.max_v
        );
        assert!(
            run.discharge_error.abs() <= 1e-3,
            "{cell_size} m: not steady ({:e})",
            run.discharge_error
        );
        assert!(
            (run.jump - x_s).abs() <= 2.0 * cell_size,
            "{cell_size} m: jump at {} m, exact {x_s} m",
            run.jump
        );
        assert!(
            run.before.abs() <= 0.02 && run.after.abs() <= 0.02,
            "{cell_size} m: depth error {} before / {} after the jump",
            run.before,
            run.after
        );
        assert!(
            run.max_error.0 <= 0.04 && run.max_error.0 < previous_max,
            "{cell_size} m: max depth error {} at {} m",
            run.max_error.0,
            run.max_error.1
        );
        previous_max = run.max_error.0;
    }
    if output_root().is_some() {
        // 図のために六角形配置も計算する（合格基準には含めない）
        let run = bump_jump_run(ChannelLayout::Hexagonal, 0.05, 300.0, "v07_hexagonal_0.05m");
        println!(
            "V07 Hexagonal, 0.05 m (for the figure): jump at {:.3} m, depth error {:+.2}% before / \
             {:+.2}% after the jump, outflow / inflow − 1 = {:.1e}, max |v| {:.1e} m/s",
            run.jump,
            100.0 * run.before,
            100.0 * run.after,
            run.discharge_error,
            run.max_v
        );
    }
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
