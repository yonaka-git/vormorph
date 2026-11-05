//! ParaView 用の VTK 出力。
//!
//! - 出力時刻ごとに VTK XML UnstructuredGrid（`.vtu`）を1ファイル書く。点はボロノイ頂点、
//!   セルは `VTK_POLYGON`（型番号 7）。
//! - 時系列は ParaView Data（`.pvd`）にまとめ、各 `.vtu` と時刻を対応付ける。ParaView で
//!   `.pvd` を開くと時間方向に再生できる。
//! - 依存を増やさないため、当面は ASCII 形式で書く。ファイルが大きくなったら base64 の
//!   バイナリ形式に切り替える。
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
//! | `c` | 水深平均の浮遊砂濃度 C | – |
//! | `dzb` | 初期からの河床変動量 | m |
//! | `velocity` | Perot 法で再構成した流速（3成分、z 成分は 0） | m/s |
//! | `x_cell`, `y_cell` | 母点の座標 | m |
//!
//! 辺の法線流速 u_n は出力せず、セル中心へ再構成した値を出す。

use crate::mesh::Mesh;
use crate::state::State;
use std::io;
use std::path::PathBuf;

/// `.vtu` の時系列と `.pvd` を書き出す。
#[derive(Debug, Clone)]
#[expect(dead_code, reason = "v0.1.0 の雛形。VTK 出力の実装で使う")]
pub struct VtkWriter {
    /// 出力先のディレクトリ
    dir: PathBuf,
    /// ファイル名の接頭辞（`<prefix>_000000.vtu`, `<prefix>.pvd`）
    prefix: String,
    /// これまでに書いた（時刻, ファイル名）
    entries: Vec<(f64, String)>,
}

impl VtkWriter {
    /// 出力先と接頭辞を決める。ディレクトリは最初の書き出しで作る。
    pub fn new(dir: impl Into<PathBuf>, prefix: impl Into<String>) -> Self {
        Self {
            dir: dir.into(),
            prefix: prefix.into(),
            entries: Vec::new(),
        }
    }

    /// 1時刻分の `.vtu` を書き、書いたファイルのパスを返す。
    ///
    /// `zb_initial` は初期の河床高（`dzb` の計算に使う）。
    pub fn write_step(
        &mut self,
        mesh: &Mesh,
        state: &State,
        zb_initial: &[f64],
    ) -> io::Result<PathBuf> {
        let _ = (mesh, state, zb_initial);
        todo!("VTK XML UnstructuredGrid の書き出し")
    }

    /// これまでに書いた全時刻の `.pvd` を書き、そのパスを返す。
    pub fn write_collection(&self) -> io::Result<PathBuf> {
        todo!("ParaView Data (.pvd) の書き出し")
    }
}
