//! コマンドライン（`vormorph check` / `vormorph run`）のテスト。
//!
//! ビルドしたバイナリを小さな設定で実行し、終了コード・表示・出力ファイルを確かめる。

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const SILL: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/sill.toml");

fn vormorph(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vormorph"))
        .args(args)
        .output()
        .expect("the vormorph binary must start")
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// テストごとの作業ディレクトリ（前回の結果は消す）。
fn work_dir(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    if dir.exists() {
        std::fs::remove_dir_all(&dir).expect("the old work directory must be removable");
    }
    std::fs::create_dir_all(&dir).expect("the work directory must be creatable");
    dir
}

/// 固定床の小さな直線水路（乾いた状態から 2 s、出力は 1 s ごと）。
fn small_case(dir: &Path) -> PathBuf {
    let config = format!(
        r#"
title = "small"

[mesh]
kind = "channel"
length = 10.0
width = 2.0
cell_size = 0.5
layout = "hexagonal"

[bed]
slope = 0.01

[[structures]]
kind = "sill"
x_start = 4.0
x_end = 5.0
height = 0.3

[time]
courant = 0.5
dt_max = 0.1
end = 2.0

[hydro]
theta = 0.55
manning_n = 0.03

[boundary.upstream]
kind = "discharge"
discharge = 0.5

[boundary.downstream]
kind = "water_level"
level = 0.2

[output]
dir = {:?}
interval = 1.0
"#,
        dir.join("out")
    );
    let path = dir.join("small.toml");
    std::fs::write(&path, config).expect("the config must be writable");
    path
}

#[test]
fn check_prints_the_summary_of_the_example() {
    let out = vormorph(&["check", SILL]);
    let stdout = text(&out.stdout);
    assert!(out.status.success(), "{}", text(&out.stderr));
    assert!(stdout.contains("OK"), "{stdout}");
    assert!(stdout.contains("adaptive dt (Courant 0.8"), "{stdout}");
}

#[test]
fn run_writes_vtk_and_reports_progress() {
    let dir = work_dir("cli-run");
    let config = small_case(&dir);
    let out = vormorph(&["run", config.to_str().expect("utf-8 path")]);
    let stdout = text(&out.stdout);
    assert!(out.status.success(), "{}\n{stdout}", text(&out.stderr));
    assert!(stdout.contains("running 'small'"), "{stdout}");
    assert_eq!(stdout.matches("| water balance").count(), 2, "{stdout}");
    assert!(stdout.contains("finished: t = 2 s"), "{stdout}");

    let pvd = std::fs::read_to_string(dir.join("out/small.pvd")).expect("the .pvd must exist");
    assert_eq!(pvd.matches("<DataSet").count(), 3, "{pvd}");
    for k in 0..3 {
        assert!(dir.join(format!("out/small_{k:06}.vtu")).exists());
    }
}

#[test]
fn run_reports_bed_changes_on_a_movable_bed() {
    let dir = work_dir("cli-movable");
    let config = small_case(&dir);
    let text_config = std::fs::read_to_string(&config).expect("the config must be readable");
    let movable = text_config.replace(
        "[boundary.upstream]",
        "[sediment]\ngrain_size = 0.002\nmorfac = 10.0\nangle_of_repose = 32.0\n\n\
         [sediment.bedload]\nformula = \"mpm\"\ncritical_shields = 0.047\n\n\
         [boundary.upstream]",
    );
    std::fs::write(&config, movable).expect("the config must be writable");
    let out = vormorph(&["run", config.to_str().expect("utf-8 path")]);
    let stdout = text(&out.stdout);
    assert!(out.status.success(), "{}\n{stdout}", text(&out.stderr));
    assert_eq!(stdout.matches("| bed change <=").count(), 2, "{stdout}");
    let vtu =
        std::fs::read_to_string(dir.join("out/small_000002.vtu")).expect("the .vtu must exist");
    assert!(vtu.contains("Name=\"qb\""));
    assert!(vtu.contains("Name=\"fixed_bed\""));
}

#[test]
fn run_rejects_suspended_sediment_until_it_is_implemented() {
    let dir = work_dir("cli-suspended");
    let text_config = std::fs::read_to_string(SILL)
        .expect("examples/sill.toml must be readable")
        .replacen("# [sediment.suspended]", "[sediment.suspended]", 1)
        .replacen("# entrainment =", "entrainment =", 1)
        .replacen("# min_depth =", "min_depth =", 1)
        // 浮遊砂の節の行（二次流の節の同名の行より前にある）
        .replacen("# diffusivity_coef =", "diffusivity_coef =", 1);
    let config = dir.join("suspended.toml");
    std::fs::write(&config, text_config).expect("the config must be writable");
    let out = vormorph(&["run", config.to_str().expect("utf-8 path")]);
    assert!(!out.status.success());
    let stderr = text(&out.stderr);
    assert!(
        stderr.contains("suspended sediment is not implemented yet"),
        "{stderr}"
    );
}

#[test]
fn unknown_commands_print_the_usage() {
    let out = vormorph(&["simulate"]);
    assert_eq!(out.status.code(), Some(2));
    assert!(text(&out.stderr).contains("Usage:"));
}
