//! VTK 出力（`io::VtkWriter`）が、格子と状態を読み戻せる形で書くことを確かめる。
//!
//! 出力は `target/tmp/vtk_output/` に残るので、ParaView や meshio で開いて確かめられる。
//!
//! ```text
//! uv run --project tools python -c \
//!     "import meshio; print(meshio.read('target/tmp/vtk_output/fields/case_000000.vtu'))"
//! ```

use std::fs;
use std::io;
use std::path::PathBuf;
use vormorph::config::ChannelLayout;
use vormorph::hydro;
use vormorph::io::VtkWriter;
use vormorph::mesh::Mesh;
use vormorph::mesh::channel::{self, ChannelSpec};
use vormorph::state::State;

const DRY_DEPTH: f64 = 1.0e-3;

/// テストごとの空の出力先。
fn output_dir(test: &str) -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
        .join("vtk_output")
        .join(test);
    let _ = fs::remove_dir_all(&dir);
    dir
}

/// 六角形配置の小さな水路（3 行、17 セル）と、セルごとに値の異なる状態。
///
/// 最後のセルは干出させ（η = z_b）、`h` と `c` が 0 になることを確かめる。
fn case() -> (Mesh, State, Vec<f64>) {
    let mesh = channel::build(&ChannelSpec {
        length: 6.0,
        width: 3.0,
        cell_size: 1.0,
        layout: ChannelLayout::Hexagonal,
    })
    .expect("valid channel");
    let n = mesh.n_cells();
    let mut state = State::new(n, mesh.n_edges());
    state.time = 12.5;
    for c in 0..n {
        state.zb[c] = 0.2 + 0.05 * mesh.y_cell[c];
        state.eta[c] = 1.0 + 0.01 * mesh.x_cell[c];
        state.hc[c] = 0.002 * state.depth(c);
    }
    state.eta[n - 1] = state.zb[n - 1];
    state.hc[n - 1] = 0.0;
    state.u_n = mesh
        .normal_edge
        .iter()
        .map(|n| 0.7 * n[0] - 0.3 * n[1])
        .collect();
    hydro::reconstruct_cell_velocity(&mesh, &state.u_n, &mut state.u_cell);
    let zb_initial = (0..n).map(|c| state.zb[c] - 1.0e-3 * c as f64).collect();
    (mesh, state, zb_initial)
}

/// `from` 以降で最初に現れる DataArray の中身。
fn array_after(text: &str, from: &str) -> Vec<f64> {
    let start = text
        .find(from)
        .unwrap_or_else(|| panic!("{from} not found"));
    let rest = &text[start..];
    let rest = &rest[rest.find("<DataArray").expect("DataArray")..];
    let body = &rest[rest.find('>').expect("end of tag") + 1..];
    let body = &body[..body.find("</DataArray>").expect("closing tag")];
    body.split_whitespace()
        .map(|t| t.parse().unwrap_or_else(|_| panic!("not a number: {t}")))
        .collect()
}

/// `Name="name"` の DataArray の中身。
fn array(text: &str, name: &str) -> Vec<f64> {
    // `Name="x"` を含むタグの先頭から探す
    let tag = format!("Name=\"{name}\"");
    let start = text
        .find(&tag)
        .unwrap_or_else(|| panic!("DataArray {name} not found"));
    let tag_start = text[..start].rfind("<DataArray").expect("DataArray");
    array_after(&text[tag_start..], "<DataArray")
}

fn as_indices(values: Vec<f64>) -> Vec<usize> {
    values.into_iter().map(|v| v as usize).collect()
}

#[test]
fn cells_points_and_fields_are_written() {
    let (mesh, state, zb_initial) = case();
    let n = mesh.n_cells();
    assert_eq!(n, 17);
    let mut writer = VtkWriter::new(output_dir("fields"), "case", DRY_DEPTH);
    let path = writer
        .write_step(&mesh, &state, &zb_initial)
        .expect("write .vtu");
    assert!(path.ends_with("case_000000.vtu"), "{}", path.display());
    let text = fs::read_to_string(&path).expect("read back");

    assert!(text.contains(&format!(
        "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{n}\">",
        mesh.n_vertices()
    )));
    assert_eq!(array(&text, "TimeValue"), vec![12.5]);

    // 点：ボロノイ頂点（z = 0）
    let points = array_after(&text, "<Points>");
    let expected: Vec<f64> = (0..mesh.n_vertices())
        .flat_map(|v| [mesh.x_vertex[v], mesh.y_vertex[v], 0.0])
        .collect();
    assert_eq!(points, expected);

    // セル：頂点の並び、終わりの位置、型（多角形 = 7）
    let connectivity: Vec<usize> = mesh.vertices_on_cell.rows().flatten().copied().collect();
    assert_eq!(as_indices(array(&text, "connectivity")), connectivity);
    let offsets: Vec<usize> = mesh
        .vertices_on_cell
        .rows()
        .scan(0, |end, row| {
            *end += row.len();
            Some(*end)
        })
        .collect();
    assert_eq!(as_indices(array(&text, "offsets")), offsets);
    assert_eq!(as_indices(array(&text, "types")), vec![7; n]);

    // セルデータ：最短の10進表記で書くので、読み戻した値は元の値と一致する
    let h: Vec<f64> = (0..n).map(|c| state.depth(c)).collect();
    let conc: Vec<f64> = (0..n).map(|c| state.concentration(c, DRY_DEPTH)).collect();
    let dzb: Vec<f64> = (0..n).map(|c| state.zb[c] - zb_initial[c]).collect();
    assert_eq!(array(&text, "eta"), state.eta);
    assert_eq!(array(&text, "zb"), state.zb);
    assert_eq!(array(&text, "h"), h);
    assert_eq!(array(&text, "c"), conc);
    assert_eq!(array(&text, "dzb"), dzb);
    assert_eq!(array(&text, "x_cell"), mesh.x_cell);
    assert_eq!(array(&text, "y_cell"), mesh.y_cell);
    let velocity: Vec<f64> = state
        .u_cell
        .iter()
        .flat_map(|u| [u[0], u[1], 0.0])
        .collect();
    assert_eq!(array(&text, "velocity"), velocity);

    // 干出させたセル
    assert_eq!((h[n - 1], conc[n - 1]), (0.0, 0.0));
    assert!(conc[0] > 0.0);
}

#[test]
fn collection_lists_every_step() {
    let (mesh, mut state, zb_initial) = case();
    let dir = output_dir("collection");
    let mut writer = VtkWriter::new(&dir, "case", DRY_DEPTH);
    for time in [0.0, 30.0, 60.0] {
        state.time = time;
        writer
            .write_step(&mesh, &state, &zb_initial)
            .expect("write .vtu");
    }
    let path = writer.write_collection().expect("write .pvd");
    assert_eq!(path, dir.join("case.pvd"));
    let text = fs::read_to_string(&path).expect("read back");

    let datasets: Vec<&str> = text.lines().filter(|l| l.contains("<DataSet")).collect();
    assert_eq!(datasets.len(), 3);
    for (k, (line, time)) in datasets.iter().zip(["0", "30", "60"]).enumerate() {
        let file = format!("case_{k:06}.vtu");
        assert!(line.contains(&format!("timestep=\"{time}\"")), "{line}");
        assert!(line.contains(&format!("file=\"{file}\"")), "{line}");
        assert!(dir.join(&file).is_file());
    }
    assert_eq!(writer.entries().len(), 3);
}

#[test]
fn invalid_states_are_rejected_without_writing() {
    let (mesh, mut state, zb_initial) = case();
    let dir = output_dir("invalid");
    let mut writer = VtkWriter::new(&dir, "case", DRY_DEPTH);

    let err = writer
        .write_step(&mesh, &state, &zb_initial[1..])
        .expect_err("zb_initial is too short");
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    assert!(err.to_string().contains("zb_initial"), "{err}");

    state.u_cell[3] = [f64::NAN, 0.0];
    let err = writer
        .write_step(&mesh, &state, &zb_initial)
        .expect_err("NaN velocity");
    assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    assert!(err.to_string().contains("velocity[6]"), "{err}");

    assert!(writer.entries().is_empty());
    assert!(!dir.exists(), "nothing must be written");
}
