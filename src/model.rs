//! 時間ループと質量収支。
//!
//! ## 1ステップの流れ
//!
//! 1. [`HydroSolver::step`] で η と u_n を時刻 n+1 へ進める。
//! 2. t ≥ `time.morphology_start` かつ移動床なら、
//!    [`SedimentSolver::compute_bedload`]、[`SedimentSolver::step_suspended`]、
//!    [`SedimentSolver::update_bed`] の順に流砂と河床を更新する。
//! 3. 水と土砂の収支（[`MassBalance`], [`SedimentBalance`]）を更新する。
//! 4. 出力時刻なら [`VtkWriter::write_step`] で書き出す。
//!
//! ## 質量収支
//!
//! 水の体積は、流入・流出の累積と、河床更新のときに湛水セルの水位を保ったことによる
//! 体積の変化で説明できなければならない。
//!
//! ```text
//! V(t) − V(0) = Σ (Q_in − Q_out) Δt + ΔV_bed
//! ```
//!
//! 半陰解法では、Newton 反復と PCG の収束判定値の範囲でこれが成り立つ。収支の誤差は
//! 検証ケース（閉鎖水域の振動など）で確かめる。

use crate::config::{Config, ConfigError};
use crate::hydro::{HydroError, HydroSolver};
use crate::io::VtkWriter;
use crate::mesh::Mesh;
use crate::sediment::SedimentSolver;
use crate::state::State;
use std::error::Error;
use std::fmt;

/// 水の収支（累積）\[m³\]。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct MassBalance {
    /// 初期の水の体積。
    pub initial_volume: f64,
    /// 現在の水の体積。
    pub volume: f64,
    /// 上流端からの累積流入量。
    pub inflow: f64,
    /// 下流端からの累積流出量。
    pub outflow: f64,
    /// 河床更新のときに湛水セルの水位を保ったことによる、体積の累積変化。
    pub bed_exchange: f64,
}

impl MassBalance {
    /// 収支の誤差 V − V₀ − (流入 − 流出) − ΔV_bed \[m³\]。
    pub fn error(&self) -> f64 {
        self.volume - self.initial_volume - (self.inflow - self.outflow) - self.bed_exchange
    }
}

/// 土砂の収支（累積、空隙を含まない実体積）\[m³\]。
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SedimentBalance {
    /// 初期からの河床の土砂量の変化（(1 − λ) Σ A_i Δz_{b,i}）。
    pub bed_change: f64,
    /// 現在の浮遊砂の量（Σ A_i (hC)_i）。
    pub suspended: f64,
    /// 初期の浮遊砂の量。
    pub initial_suspended: f64,
    /// 上流端からの累積給砂量。
    pub inflow: f64,
    /// 下流端からの累積流出量。
    pub outflow: f64,
}

impl SedimentBalance {
    /// 収支の誤差 \[m³\]。
    pub fn error(&self) -> f64 {
        self.bed_change + (self.suspended - self.initial_suspended) - (self.inflow - self.outflow)
    }
}

/// モデルの組み立てと実行の失敗。
#[derive(Debug)]
pub enum ModelError {
    /// 設定が正しくない。
    Config(ConfigError),
    /// 流れの計算が失敗した。
    Hydro(HydroError),
    /// 入出力の失敗。
    Io(std::io::Error),
}

impl fmt::Display for ModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Config(e) => write!(f, "configuration error: {e}"),
            Self::Hydro(e) => write!(f, "hydrodynamic solver error: {e}"),
            Self::Io(e) => write!(f, "I/O error: {e}"),
        }
    }
}

impl Error for ModelError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Config(e) => Some(e),
            Self::Hydro(e) => Some(e),
            Self::Io(e) => Some(e),
        }
    }
}

impl From<ConfigError> for ModelError {
    fn from(e: ConfigError) -> Self {
        Self::Config(e)
    }
}

impl From<HydroError> for ModelError {
    fn from(e: HydroError) -> Self {
        Self::Hydro(e)
    }
}

impl From<std::io::Error> for ModelError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

/// 河床変動モデル。格子・状態・各ソルバー・出力をまとめ、時間ループを回す。
#[derive(Debug)]
#[expect(dead_code, reason = "v0.1.0 の雛形。時間ループの実装で使う")]
pub struct Model {
    config: Config,
    mesh: Mesh,
    state: State,
    /// 初期の河床高（出力の dzb に使う）
    zb_initial: Vec<f64>,
    hydro: HydroSolver,
    sediment: Option<SedimentSolver>,
    writer: VtkWriter,
    water: MassBalance,
    sediment_balance: SedimentBalance,
    steps: u64,
}

impl Model {
    /// 設定から格子と初期条件を作り、モデルを組み立てる。
    pub fn new(config: Config) -> Result<Self, ModelError> {
        let _ = config;
        todo!("格子の生成・初期条件の設定")
    }

    /// 格子と初期状態を外から与えて組み立てる（検証ケースで任意の地形を与えるため）。
    pub fn with_state(config: Config, mesh: Mesh, state: State) -> Result<Self, ModelError> {
        let _ = (config, mesh, state);
        todo!("与えられた格子と状態からの組み立て")
    }

    /// 格子。
    pub fn mesh(&self) -> &Mesh {
        &self.mesh
    }

    /// 現在の状態。
    pub fn state(&self) -> &State {
        &self.state
    }

    /// 水の収支。
    pub fn water_balance(&self) -> &MassBalance {
        &self.water
    }

    /// 土砂の収支。
    pub fn sediment_balance(&self) -> &SedimentBalance {
        &self.sediment_balance
    }

    /// 1ステップ進める。
    pub fn step(&mut self) -> Result<(), ModelError> {
        todo!("1ステップ（流れ → 流砂 → 河床 → 収支 → 出力）")
    }

    /// 設定の終了時刻まで進める。
    pub fn run(&mut self) -> Result<(), ModelError> {
        todo!("時間ループ")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_balance_error_accounts_for_all_terms() {
        let b = MassBalance {
            initial_volume: 100.0,
            volume: 112.0,
            inflow: 20.0,
            outflow: 5.0,
            bed_exchange: -3.0,
        };
        assert_eq!(b.error(), 0.0);
    }
}
