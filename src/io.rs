//! ParaView 用の VTK 出力。
//!
//! - 出力時刻ごとに VTK XML UnstructuredGrid（`.vtu`）を1ファイル書く。点はボロノイ頂点
//!   （z = 0）、セルは `VTK_POLYGON`（型番号 7）で、頂点の並びは `vertices_on_cell` のまま
//!   （反時計回り）。
//! - 時系列は ParaView Data（`.pvd`）にまとめ、各 `.vtu` と時刻を対応付ける。ParaView で
//!   `.pvd` を開くと時間方向に再生できる。各 `.vtu` にも時刻を `TimeValue`（FieldData）として
//!   書くので、1ファイルだけ開いても時刻が分かる。
//! - 依存を増やさないため、当面は ASCII 形式で書く。ファイルが大きくなったら base64 の
//!   バイナリ形式に切り替える。
//!
//! ## 数値の書き方
//!
//! 浮動小数点数は、読み戻すと同じ値になる最短の10進表記で書く（Rust の `{}` / `{:e}`）。
//! 絶対値が 10⁻⁴ 未満または 10¹⁵ 以上の値は指数表記にする。NaN と無限大は ASCII 形式の
//! VTK の読み込み側が解釈できないので、書き出さずに誤り（[`io::ErrorKind::InvalidData`]）にする。
//!
//! ## セルデータの名前
//!
//! `tools/plot_profile.py` もこの名前で読むので、変えるときは両方を直す。
//!
//! | 名前 | 内容 | 単位 |
//! |---|---|---|
//! | `eta` | 水位 η | m |
//! | `zb` | 河床高 z_b | m |
//! | `h` | 水深 h | m |
//! | `c` | 水深平均の浮遊砂濃度 C（水深が `dry_depth` 未満のセルは 0） | – |
//! | `dzb` | 初期からの河床変動量 | m |
//! | `velocity` | Perot 法で再構成した流速（3成分、z 成分は 0） | m/s |
//! | `x_cell`, `y_cell` | 母点の座標 | m |
//! | `qb` | 掃流砂量ベクトル（移動床のときだけ。3成分、z 成分は 0） | m²/s |
//! | `fixed_bed` | 非侵食層のあるセルなら 1、なければ 0（移動床のときだけ） | – |
//!
//! 辺の法線流速 u_n は出力せず、セル中心へ再構成した値（[`State::u_cell`]）を出す。掃流砂量も
//! 同じく、辺の法線成分 q_{b,j} から Perot 法と同じ式で再構成したセル中心のベクトルを出す
//! （一様な掃流砂量なら厳密）。MORFAC は掛けない（流れの時間あたりの量）。

use crate::hydro::reconstruct_cell_velocity;
use crate::mesh::Mesh;
use crate::state::State;
use std::borrow::Cow;
use std::fs;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

/// `.vtu` の時系列と `.pvd` を書き出す。
#[derive(Debug, Clone)]
pub struct VtkWriter {
    /// 出力先のディレクトリ
    dir: PathBuf,
    /// ファイル名の接頭辞（`<prefix>_000000.vtu`, `<prefix>.pvd`）
    prefix: String,
    /// これより浅いセルの濃度を 0 として出力する水深 \[m\]
    dry_depth: f64,
    /// これまでに書いた（時刻, ファイル名）
    entries: Vec<(f64, String)>,
}

impl VtkWriter {
    /// 出力先と接頭辞を決める。ディレクトリは最初の書き出しで作る。
    ///
    /// `dry_depth` は濃度 `c` を出力するときの wet/dry の閾値（[`State::concentration`]）。
    pub fn new(dir: impl Into<PathBuf>, prefix: impl Into<String>, dry_depth: f64) -> Self {
        Self {
            dir: dir.into(),
            prefix: prefix.into(),
            dry_depth,
            entries: Vec::new(),
        }
    }

    /// これまでに書いた（時刻 \[s\], ファイル名）。
    pub fn entries(&self) -> &[(f64, String)] {
        &self.entries
    }

    /// 1時刻分の `.vtu` を書き、書いたファイルのパスを返す。
    ///
    /// `zb_initial` は初期の河床高（`dzb` の計算に使う）。`bedload` は辺の掃流砂量 q_{b,j}
    /// \[m²/s\]（移動床のとき。`Some` なら `qb` と `fixed_bed` も書く）。配列の長さが格子と
    /// 合わない場合は [`io::ErrorKind::InvalidInput`]、出力する値に NaN や無限大がある場合は
    /// [`io::ErrorKind::InvalidData`] の誤りを返し、ファイルを書かない。
    pub fn write_step(
        &mut self,
        mesh: &Mesh,
        state: &State,
        zb_initial: &[f64],
        bedload: Option<&[f64]>,
    ) -> io::Result<PathBuf> {
        check_lengths(mesh, state, zb_initial, bedload)?;
        let fields = CellFields::new(mesh, state, zb_initial, bedload, self.dry_depth);
        fields.check_finite(mesh, state.time)?;

        fs::create_dir_all(&self.dir)?;
        let name = format!("{}_{:06}.vtu", self.prefix, self.entries.len());
        let path = self.dir.join(&name);
        write_file(&path, |w| write_vtu(w, mesh, state.time, &fields))?;
        self.entries.push((state.time, name));
        Ok(path)
    }

    /// これまでに書いた全時刻の `.pvd` を書き、そのパスを返す。
    ///
    /// 出力のたびに呼んで書き直すと、計算の途中でも ParaView で開ける。
    pub fn write_collection(&self) -> io::Result<PathBuf> {
        fs::create_dir_all(&self.dir)?;
        let path = self.dir.join(format!("{}.pvd", self.prefix));
        write_file(&path, |w| {
            writeln!(w, "<?xml version=\"1.0\"?>")?;
            writeln!(
                w,
                "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"
            )?;
            writeln!(w, "  <Collection>")?;
            for (time, name) in &self.entries {
                writeln!(
                    w,
                    "    <DataSet timestep=\"{}\" group=\"\" part=\"0\" file=\"{}\"/>",
                    number(*time),
                    escape(name)
                )?;
            }
            writeln!(w, "  </Collection>")?;
            writeln!(w, "</VTKFile>")
        })?;
        Ok(path)
    }
}

/// 出力するセルデータ。状態と格子にある配列は借り、計算で求める配列だけを持つ。
struct CellFields<'a> {
    /// スカラー（名前、値）
    scalars: Vec<(&'static str, Cow<'a, [f64]>)>,
    /// 平面のベクトル（名前、値）。最初は流速（`velocity`）
    vectors: Vec<(&'static str, Cow<'a, [[f64; 2]]>)>,
}

impl<'a> CellFields<'a> {
    fn new(
        mesh: &'a Mesh,
        state: &'a State,
        zb_initial: &[f64],
        bedload: Option<&[f64]>,
        dry_depth: f64,
    ) -> Self {
        let n = mesh.n_cells();
        let h = (0..n).map(|c| state.depth(c)).collect();
        let c = (0..n).map(|c| state.concentration(c, dry_depth)).collect();
        let dzb = state
            .zb
            .iter()
            .zip(zb_initial)
            .map(|(z, z0)| z - z0)
            .collect();
        let mut scalars = vec![
            ("eta", Cow::Borrowed(state.eta.as_slice())),
            ("zb", Cow::Borrowed(state.zb.as_slice())),
            ("h", Cow::Owned(h)),
            ("c", Cow::Owned(c)),
            ("dzb", Cow::Owned(dzb)),
            ("x_cell", Cow::Borrowed(mesh.x_cell.as_slice())),
            ("y_cell", Cow::Borrowed(mesh.y_cell.as_slice())),
        ];
        let mut vectors = vec![("velocity", Cow::Borrowed(state.u_cell.as_slice()))];
        if let Some(bedload) = bedload {
            let mut qb = vec![[0.0; 2]; n];
            reconstruct_cell_velocity(mesh, bedload, &mut qb);
            vectors.push(("qb", Cow::Owned(qb)));
            let fixed = state
                .z_ne
                .iter()
                .map(|z| if z.is_finite() { 1.0 } else { 0.0 })
                .collect();
            scalars.push(("fixed_bed", Cow::Owned(fixed)));
        }
        Self { scalars, vectors }
    }

    /// 書き出す値がすべて有限であることを確かめる。
    fn check_finite(&self, mesh: &Mesh, time: f64) -> io::Result<()> {
        check_finite("time", [time])?;
        check_finite("x_vertex", mesh.x_vertex.iter().copied())?;
        check_finite("y_vertex", mesh.y_vertex.iter().copied())?;
        for (name, values) in &self.scalars {
            check_finite(name, values.iter().copied())?;
        }
        for (name, values) in &self.vectors {
            check_finite(name, values.iter().flatten().copied())?;
        }
        Ok(())
    }
}

/// `.vtu` の本文を書く。
fn write_vtu(w: &mut impl Write, mesh: &Mesh, time: f64, fields: &CellFields) -> io::Result<()> {
    writeln!(w, "<?xml version=\"1.0\"?>")?;
    writeln!(
        w,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
    )?;
    writeln!(w, "  <UnstructuredGrid>")?;
    writeln!(w, "    <FieldData>")?;
    writeln!(
        w,
        "      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">"
    )?;
    writeln!(w, "        {}", number(time))?;
    writeln!(w, "      </DataArray>")?;
    writeln!(w, "    </FieldData>")?;
    writeln!(
        w,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        mesh.n_vertices(),
        mesh.n_cells()
    )?;

    writeln!(w, "      <Points>")?;
    begin_array(w, "Float64", None, 3)?;
    for (x, y) in mesh.x_vertex.iter().zip(&mesh.y_vertex) {
        writeln!(w, "          {} {} 0", number(*x), number(*y))?;
    }
    end_array(w)?;
    writeln!(w, "      </Points>")?;

    writeln!(w, "      <Cells>")?;
    begin_array(w, "Int64", Some("connectivity"), 1)?;
    for verts in mesh.vertices_on_cell.rows() {
        write!(w, "         ")?;
        for v in verts {
            write!(w, " {v}")?;
        }
        writeln!(w)?;
    }
    end_array(w)?;
    begin_array(w, "Int64", Some("offsets"), 1)?;
    let mut offset = 0;
    for verts in mesh.vertices_on_cell.rows() {
        offset += verts.len();
        writeln!(w, "          {offset}")?;
    }
    end_array(w)?;
    begin_array(w, "UInt8", Some("types"), 1)?;
    for _ in 0..mesh.n_cells() {
        writeln!(w, "          {VTK_POLYGON}")?;
    }
    end_array(w)?;
    writeln!(w, "      </Cells>")?;

    writeln!(w, "      <CellData Scalars=\"h\" Vectors=\"velocity\">")?;
    for (name, values) in &fields.scalars {
        begin_array(w, "Float64", Some(name), 1)?;
        for v in values.iter() {
            writeln!(w, "          {}", number(*v))?;
        }
        end_array(w)?;
    }
    for (name, values) in &fields.vectors {
        begin_array(w, "Float64", Some(name), 3)?;
        for [u, v] in values.iter() {
            writeln!(w, "          {} {} 0", number(*u), number(*v))?;
        }
        end_array(w)?;
    }
    writeln!(w, "      </CellData>")?;

    writeln!(w, "    </Piece>")?;
    writeln!(w, "  </UnstructuredGrid>")?;
    writeln!(w, "</VTKFile>")
}

/// VTK のセルの型番号（多角形）。
const VTK_POLYGON: u8 = 7;

fn begin_array(
    w: &mut impl Write,
    ty: &str,
    name: Option<&str>,
    components: usize,
) -> io::Result<()> {
    write!(w, "        <DataArray type=\"{ty}\"")?;
    if let Some(name) = name {
        write!(w, " Name=\"{name}\"")?;
    }
    writeln!(w, " NumberOfComponents=\"{components}\" format=\"ascii\">")
}

fn end_array(w: &mut impl Write) -> io::Result<()> {
    writeln!(w, "        </DataArray>")
}

/// 浮動小数点数を、読み戻すと同じ値になる最短の10進表記にする。
///
/// 絶対値が 10⁻⁴ 未満（0 を除く）または 10¹⁵ 以上なら指数表記にする。`{}` は指数表記を
/// 使わないので、極端な値で桁数が増えすぎないようにするためである。
fn number(v: f64) -> String {
    if v == 0.0 || (1.0e-4..1.0e15).contains(&v.abs()) {
        format!("{v}")
    } else {
        format!("{v:e}")
    }
}

/// XML の属性値に入れる文字列のエスケープ。
fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(ch),
        }
    }
    out
}

fn check_lengths(
    mesh: &Mesh,
    state: &State,
    zb_initial: &[f64],
    bedload: Option<&[f64]>,
) -> io::Result<()> {
    let (n, n_edges) = (mesh.n_cells(), mesh.n_edges());
    let lengths = [
        ("state.eta", state.eta.len(), n, "cells"),
        ("state.zb", state.zb.len(), n, "cells"),
        ("state.z_ne", state.z_ne.len(), n, "cells"),
        ("state.hc", state.hc.len(), n, "cells"),
        ("state.u_cell", state.u_cell.len(), n, "cells"),
        ("zb_initial", zb_initial.len(), n, "cells"),
        (
            "bedload",
            bedload.map_or(n_edges, <[f64]>::len),
            n_edges,
            "edges",
        ),
    ];
    match lengths
        .iter()
        .find(|&&(_, len, expected, _)| len != expected)
    {
        Some((name, len, expected, what)) => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("{name} has {len} values but the mesh has {expected} {what}"),
        )),
        None => Ok(()),
    }
}

fn check_finite(name: &str, values: impl IntoIterator<Item = f64>) -> io::Result<()> {
    match values.into_iter().enumerate().find(|(_, v)| !v.is_finite()) {
        Some((k, v)) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{name}[{k}] = {v} cannot be written to an ASCII VTK file"),
        )),
        None => Ok(()),
    }
}

/// ファイルを作り、`body` で書いて閉じる。
fn write_file(
    path: &Path,
    body: impl FnOnce(&mut BufWriter<fs::File>) -> io::Result<()>,
) -> io::Result<()> {
    let mut w = BufWriter::new(fs::File::create(path)?);
    body(&mut w)?;
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_round_trip_and_stay_short() {
        for v in [
            0.0, -0.0, 1.0, -2.5, 0.1, 1.0e-4, 123.456, 1.0e-17, 6.02e23, -3.0e-300,
        ] {
            let text = number(v);
            assert_eq!(text.parse::<f64>(), Ok(v), "{text}");
            assert!(text.len() <= 24, "{text}");
        }
        assert_eq!(number(0.5), "0.5");
        assert_eq!(number(1.0e-17), "1e-17");
    }

    #[test]
    fn attribute_values_are_escaped() {
        assert_eq!(
            escape(r#"a&b<"c">'d'"#),
            "a&amp;b&lt;&quot;c&quot;&gt;&apos;d&apos;"
        );
    }
}
