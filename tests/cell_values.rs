//! 格子生成ツールの格子ファイルと、セルの値のファイル（河床高）から、設定でモデルを組み立てる
//! （ADR-0020）。
//!
//! `tests/data/straight_coarse.cells.toml` は、`tests/data/straight_coarse.mesh.toml` の母点で
//! `examples/domains/straight.toml` の `[bed]`（勾配 1/500、x = 0 で 0.12 m の平面）を評価したもの。
//!
//! ```text
//! uv run --project tools python tools/make_bed.py examples/domains/straight.toml \
//!     tests/data/straight_coarse.mesh.toml -o tests/data/straight_coarse.cells.toml
//! ```
//!
//! 格子ファイルを作り直したら、このファイルも作り直す（母点の座標で照合するので、作り直さないと
//! 読み込みの誤りになる）。

use std::path::{Path, PathBuf};
use vormorph::config::Config;
use vormorph::mesh::cells::CellValuesFile;
use vormorph::model::{Model, ModelError};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// 格子ファイルの直線水路（延長 60 m、幅 10 m）に、上流端から 20 m³/s を流し、下流端の水位を
/// 等流水位（1.07 m）にする設定。
fn config(name: &str, mesh: &Path, cells: &Path) -> Config {
    let dir = std::env::temp_dir().join(format!("vormorph-{name}-{}", std::process::id()));
    format!(
        r#"
title = "{name}"
[mesh]
kind = "file"
path = {mesh:?}
cells = {cells:?}
[time]
courant = 0.5
dt_max = 0.5
end = 60.0
[hydro]
theta = 0.55
manning_n = 0.025
[boundary.upstream]
kind = "discharge"
discharge = 20.0
[boundary.downstream]
kind = "water_level"
level = 1.07
[output]
dir = {dir:?}
interval = 1000.0
"#
    )
    .parse()
    .expect("valid config")
}

/// 河床高はセルの値のファイルの値（母点で評価した平面）で、初期状態は乾いた河床。流れを進めても
/// 水の収支が合う。
#[test]
fn model_is_built_from_the_mesh_and_cell_values_files() {
    let cells = data("straight_coarse.cells.toml");
    let config = config("file-mesh", &data("straight_coarse.mesh.toml"), &cells);
    assert!(config.summary().to_string().contains("cell values"));
    let mut model = Model::new(config).expect("the model must build from the files");
    let (mesh, state) = (model.mesh(), model.state());
    assert_eq!(mesh.n_cells(), 173);
    let values = CellValuesFile::load(&cells).expect("the fixture must be valid");
    assert_eq!(state.zb, values.cells.zb);
    for c in 0..mesh.n_cells() {
        let plane = 0.12 - 0.002 * mesh.x_cell[c];
        assert!((state.zb[c] - plane).abs() < 1e-15, "cell {c}");
        assert_eq!(state.eta[c], state.zb[c], "the initial bed must be dry");
        assert_eq!(state.z_ne[c], f64::NEG_INFINITY);
    }
    while !model.finished() {
        model.step().expect("the step must succeed");
    }
    let water = model.water_balance();
    assert!(water.inflow > 0.0);
    assert!(water.error().abs() <= 1e-10 * water.inflow, "{water:?}");
}

/// 別の格子のセルの値、母点のずれたセルの値、ないファイルは、組み立ての時点で誤りになる。
#[test]
fn cell_values_of_another_mesh_are_rejected() {
    let cells = data("straight_coarse.cells.toml");
    let island = data("island_coarse.mesh.toml");
    let err = Model::new(config("other-mesh", &island, &cells)).unwrap_err();
    assert!(
        matches!(&err, ModelError::CellValues { .. })
            && err.to_string().contains("cells but the mesh has"),
        "{err}"
    );

    // 母点を1つずらしたファイル（格子を作り直したのにセルの値を作り直さなかった場合）
    let text = std::fs::read_to_string(&cells).expect("the fixture must be readable");
    let first = text
        .find("x = [\n  ")
        .map(|i| i + "x = [\n  ".len())
        .expect("cells.x");
    let end = first + text[first..].find(',').expect("a value");
    let moved = format!("{}{}{}", &text[..first], "0.9", &text[end..]);
    let dir = std::env::temp_dir().join(format!("vormorph-cells-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("a temporary directory");
    let path = dir.join("moved.cells.toml");
    std::fs::write(&path, moved).expect("the file must be writable");
    let mesh = data("straight_coarse.mesh.toml");
    let err = Model::new(config("moved", &mesh, &path)).unwrap_err();
    assert!(err.to_string().contains("cell 0 is at (0.9,"), "{err}");

    let err = Model::new(config("missing", &mesh, &dir.join("missing.cells.toml"))).unwrap_err();
    assert!(
        matches!(&err, ModelError::CellValues { path, .. } if path.ends_with("missing.cells.toml")),
        "{err}"
    );
    std::fs::remove_dir_all(&dir).expect("clean up");
}
