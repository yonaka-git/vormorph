//! 状態変数（C 格子のスタガード配置）。
//!
//! | フィールド | 位置 | 記号 | 意味 | 単位 |
//! |---|---|---|---|---|
//! | `eta` | セル | η | 水位 | m |
//! | `zb` | セル | z_b | 河床高 | m |
//! | `z_ne` | セル | z_ne | 非侵食層の上面の高さ（床止めなど。なければ −∞） | m |
//! | `hc` | セル | hC | 水深と水深平均濃度の積（浮遊砂の保存量） | m |
//! | `u_n` | 辺 | u_n | 法線流速（`cells_on_edge[e][0]` から `[1]` へ向かう向きを正） | m/s |
//! | `u_cell` | セル | **u**_i | Perot 法で再構成した流速ベクトル（流砂の評価と出力に使う） | m/s |
//!
//! 水深 h = max(η − z_b, 0) は保存せず、[`State::depth`] で都度計算する。
//! 辺の水深 H_j は風上側の水位で決まるため、[`crate::hydro`] が計算する。

/// モデルの状態。
#[derive(Debug, Clone, PartialEq)]
pub struct State {
    /// 時刻 t \[s\]（流れの時間）。
    pub time: f64,
    /// 水位 η \[m\]（セル）。
    pub eta: Vec<f64>,
    /// 河床高 z_b \[m\]（セル）。
    pub zb: Vec<f64>,
    /// 非侵食層の上面の高さ z_ne \[m\]（セル）。非侵食層がなければ `f64::NEG_INFINITY`。
    pub z_ne: Vec<f64>,
    /// 水深と水深平均濃度の積 hC \[m\]（セル）。
    pub hc: Vec<f64>,
    /// 法線流速 u_n \[m/s\]（辺）。
    pub u_n: Vec<f64>,
    /// セル中心の流速ベクトル \[m/s\]（Perot 法による再構成）。
    pub u_cell: Vec<[f64; 2]>,
}

impl State {
    /// すべての変数を 0（非侵食層は −∞）で初期化する。
    pub fn new(n_cells: usize, n_edges: usize) -> Self {
        Self {
            time: 0.0,
            eta: vec![0.0; n_cells],
            zb: vec![0.0; n_cells],
            z_ne: vec![f64::NEG_INFINITY; n_cells],
            hc: vec![0.0; n_cells],
            u_n: vec![0.0; n_edges],
            u_cell: vec![[0.0; 2]; n_cells],
        }
    }

    /// セル `c` の水深 h = max(η − z_b, 0) \[m\]。
    pub fn depth(&self, c: usize) -> f64 {
        (self.eta[c] - self.zb[c]).max(0.0)
    }

    /// セル `c` の水深平均濃度 C = hC / h \[-\]。水深が `dry_depth` 未満なら 0。
    pub fn concentration(&self, c: usize, dry_depth: f64) -> f64 {
        let h = self.depth(c);
        if h < dry_depth { 0.0 } else { self.hc[c] / h }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn depth_is_clipped_at_zero() {
        let mut s = State::new(2, 1);
        s.eta = vec![1.5, 0.2];
        s.zb = vec![1.0, 0.5];
        s.hc = vec![0.01, 0.0];
        assert_eq!(s.depth(0), 0.5);
        assert_eq!(s.depth(1), 0.0);
        assert!((s.concentration(0, 1e-3) - 0.02).abs() < 1e-15);
        assert_eq!(s.concentration(1, 1e-3), 0.0);
    }
}
