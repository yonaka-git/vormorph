//! vormorph のコマンドライン。
//!
//! ```text
//! vormorph check <config.toml>   設定ファイルを読み込んで検証し、要約を表示する
//! vormorph run <config.toml>     計算を実行する（出力時刻ごとに進み具合を表示する）
//! vormorph --help | --version
//! ```
//!
//! 依存を増やさないため、引数は標準ライブラリだけで解析する。出力先（`output.dir`）は、
//! 実行したディレクトリからの相対パスとして扱う。

use std::process::ExitCode;
use std::time::Instant;
use vormorph::config::{Config, TimeStep};
use vormorph::model::{Model, Progress};

const USAGE: &str = "\
vormorph - planar 2D morphodynamic model on Voronoi C-grids

Usage:
  vormorph check <config.toml>   validate a configuration file and print a summary
  vormorph run <config.toml>     run a simulation and write VTK output
  vormorph --help | --version";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    match args.as_slice() {
        ["check", path] => check(path),
        ["run", path] => run(path),
        ["-h" | "--help"] => {
            println!("{USAGE}");
            ExitCode::SUCCESS
        }
        ["-V" | "--version"] => {
            println!("vormorph {}", env!("CARGO_PKG_VERSION"));
            ExitCode::SUCCESS
        }
        _ => {
            eprintln!("{USAGE}");
            ExitCode::from(2)
        }
    }
}

/// 設定ファイルを検証して要約を表示する。
fn check(path: &str) -> ExitCode {
    match Config::load(path) {
        Ok(config) => {
            println!("{path}: OK\n{}", config.summary());
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("error: {path}: {e}");
            ExitCode::FAILURE
        }
    }
}

/// 計算を実行する。出力時刻ごとに1行の進み具合を表示する。
fn run(path: &str) -> ExitCode {
    let config = match Config::load(path) {
        Ok(config) => config,
        Err(e) => {
            eprintln!("error: {path}: {e}");
            return ExitCode::FAILURE;
        }
    };
    let mut model = match Model::new(config) {
        Ok(model) => model,
        Err(e) => {
            eprintln!("error: {path}: {e}");
            return ExitCode::FAILURE;
        }
    };
    let (mesh, config) = (model.mesh(), model.config());
    println!(
        "running '{}': {} cells, {} edges, until t = {} s, output to {}",
        config.title,
        mesh.n_cells(),
        mesh.n_edges(),
        config.time.end,
        config.output.dir.display()
    );
    let started = Instant::now();
    let result = model.run_with(|model, progress| println!("{}", progress_line(model, progress)));
    match result {
        Ok(()) => {
            println!(
                "finished: t = {} s, {} steps, {:.1} s wall time",
                model.state().time,
                model.steps(),
                started.elapsed().as_secs_f64()
            );
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("error: at t = {} s: {e}", model.state().time);
            ExitCode::FAILURE
        }
    }
}

/// 進み具合の1行：時刻、ステップ数、時間刻み、クーラン数、反復回数、水の収支の誤差。
fn progress_line(model: &Model, progress: &Progress) -> String {
    let dt = match model.time_step() {
        TimeStep::Fixed { dt } => format!("dt {dt} s"),
        TimeStep::Adaptive { .. } => {
            format!("dt {:.4}..{:.4} s", progress.dt_min, progress.dt_max)
        }
    };
    let water = model.water_balance();
    let scale = (water.initial_volume + water.inflow).max(f64::MIN_POSITIVE);
    let mut line = format!(
        "t = {:>10.2} s | {:>7} steps | {dt} | Courant {:.3} | Newton <= {} | PCG {:.1}/step | water balance {:.1e}",
        model.state().time,
        progress.steps,
        progress.max_outflow_courant,
        progress.max_newton_iterations,
        progress.pcg_iterations as f64 / progress.steps.max(1) as f64,
        water.error() / scale,
    );
    if progress.limited_edges > 0 {
        line.push_str(&format!(
            " | advection limited on {} edges",
            progress.limited_edges
        ));
    }
    line
}
