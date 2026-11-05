//! vormorph のコマンドライン。
//!
//! ```text
//! vormorph check <config.toml>   設定ファイルを読み込んで検証し、要約を表示する
//! vormorph run <config.toml>     計算を実行する（v0.1.0 では未実装）
//! vormorph --help | --version
//! ```
//!
//! 依存を増やさないため、引数は標準ライブラリだけで解析する。

use std::process::ExitCode;
use vormorph::config::Config;

const USAGE: &str = "\
vormorph - planar 2D morphodynamic model on Voronoi C-grids

Usage:
  vormorph check <config.toml>   validate a configuration file and print a summary
  vormorph run <config.toml>     run a simulation (not implemented yet in v0.1.0)
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

/// 計算を実行する。v0.1.0 では設定の検証までで止める。
fn run(path: &str) -> ExitCode {
    let config = match Config::load(path) {
        Ok(config) => config,
        Err(e) => {
            eprintln!("error: {path}: {e}");
            return ExitCode::FAILURE;
        }
    };
    // 数値計算の本体は未実装。実装後は次のように呼ぶ。
    //   let mut model = vormorph::model::Model::new(config)?;
    //   model.run()?;
    eprintln!(
        "error: case '{}' was loaded, but the solver is not implemented yet (v0.1.0 is a scaffold)",
        config.title
    );
    ExitCode::FAILURE
}
