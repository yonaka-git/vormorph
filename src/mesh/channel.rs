//! 直線矩形水路の格子生成（検証ケースとデモ用）。
//!
//! x 軸を流下方向（上流端 x = 0、下流端 x = L）、y 軸を横断方向（側壁 y = 0 と y = W）にとる。
//!
//! - [`ChannelLayout::Rectangular`]：母点を正方格子に置く。セルは正方形で、直交誤差は 0。
//! - [`ChannelLayout::Hexagonal`]：母点を三角格子に置く。内部のセルは正六角形で、側壁と
//!   上下流端に接するセルは、ボロノイ図を水路の矩形で切り取った多角形になる。
//!
//! 境界辺の種類は、x = 0 が [`EdgeKind::Upstream`](super::EdgeKind::Upstream)、
//! x = L が [`EdgeKind::Downstream`](super::EdgeKind::Downstream)、
//! y = 0 と y = W が [`EdgeKind::Wall`](super::EdgeKind::Wall)。
//!
//! 生成した格子は [`Mesh::derive_cell_connectivity`] と
//! [`voronoi::compute_geometry`](super::voronoi::compute_geometry) を通して、
//! 符号と幾何量をほかの格子と同じ手順で計算する。

use super::Mesh;
use crate::config::{ChannelLayout, MeshConfig};

/// 直線水路の諸元。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelSpec {
    /// 水路延長 L \[m\]。
    pub length: f64,
    /// 水路幅 W \[m\]。
    pub width: f64,
    /// 代表格子幅（母点の間隔）\[m\]。
    pub cell_size: f64,
    /// 母点の配置。
    pub layout: ChannelLayout,
}

impl ChannelSpec {
    /// 設定の `[mesh]` から作る。`kind = "channel"` 以外なら `None`。
    pub fn from_config(mesh: &MeshConfig) -> Option<Self> {
        match *mesh {
            MeshConfig::Channel {
                length,
                width,
                cell_size,
                layout,
            } => Some(Self {
                length,
                width,
                cell_size,
                layout,
            }),
            MeshConfig::File { .. } => None,
        }
    }
}

/// 直線水路の格子を生成する。
///
/// 横断方向のセル数は W / `cell_size` を丸めた整数とし、母点の間隔はそれに合わせて
/// わずかに調整する（縦断方向も同様）。
pub fn build(spec: &ChannelSpec) -> Mesh {
    todo!("直線水路の格子生成 ({:?})", spec.layout)
}
