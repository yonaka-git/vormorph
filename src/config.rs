//! TOML 設定ファイルの読み込みと検証。
//!
//! 設定ファイルの例は `examples/sill.toml`。ファイルからは [`Config::load`]、文字列からは
//! `str::parse::<Config>()`（[`FromStr`] の実装）で読み込む。どちらも読み込んだ直後に
//! [`Config::validate`] で値の範囲を確かめる。未知のキーはエラーにする（綴りの誤りで設定が
//! 黙って既定値になるのを防ぐため）。
//!
//! ## 物理定数の扱い
//!
//! 重力加速度・密度・動粘性係数などの物理定数は、数値計算のモジュールに直接書かず、必ず
//! このモジュールの構造体を通して渡す。既定値を持つ項目の既定値は、このファイルの
//! `default_*` 関数と `Default` 実装にだけ置く。
//!
//! ## 単位
//!
//! すべて SI 単位（m, s, kg）。角度（水中安息角）のみ度で与える。
//!
//! ## 設定の構成
//!
//! | キー | 型 | 内容 |
//! |---|---|---|
//! | `title` | [`Config::title`] | 計算ケースの名前（出力ファイル名の接頭辞） |
//! | `[mesh]` | [`MeshConfig`] | 格子。`kind = "channel"` で直線水路を生成、`"file"` で外部格子 |
//! | `[bed]` | [`BedConfig`] | 初期河床（直線水路のとき必須。周期境界では、流れを駆動する勾配にも使う） |
//! | `[[structures]]` | [`StructureConfig`] | 床止めなどの構造物（任意、複数可） |
//! | `[physics]` | [`PhysicsConfig`] | 物理定数（省略時は既定値） |
//! | `[time]` | [`TimeConfig`] | 時間刻み・終了時刻・河床変動の開始時刻 |
//! | `[hydro]` | [`HydroConfig`] | θ、粗度係数、wet/dry の閾値、移流、時間積分、反復解法 |
//! | `[sediment]` | [`SedimentConfig`] | 粒径・空隙率・流砂量式・MORFAC・二次流など（省略時は固定床） |
//! | `[boundary]` | [`BoundaryConfig`] | 上流端・下流端の境界条件（それ以外の境界は壁。両端を `periodic` にすると周期境界） |
//! | `[output]` | [`OutputConfig`] | 出力先と出力間隔 |

use serde::Deserialize;
use std::error::Error;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// 計算ケース全体の設定。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    /// 計算ケースの名前。出力ファイル名の接頭辞にも使う。
    pub title: String,
    /// 格子。
    pub mesh: MeshConfig,
    /// 初期河床。`mesh.kind = "channel"` のとき必須。
    pub bed: Option<BedConfig>,
    /// 床止めなどの構造物。
    #[serde(default)]
    pub structures: Vec<StructureConfig>,
    /// 物理定数。
    #[serde(default)]
    pub physics: PhysicsConfig,
    /// 時間積分。
    pub time: TimeConfig,
    /// 流れの計算。
    pub hydro: HydroConfig,
    /// 流砂と河床変動。`None`（セクションを省略）のときは固定床で流れだけを解く。
    pub sediment: Option<SedimentConfig>,
    /// 境界条件。
    pub boundary: BoundaryConfig,
    /// 出力。
    pub output: OutputConfig,
}

/// 格子の作り方。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum MeshConfig {
    /// 直線矩形水路の格子を内部で生成する（[`crate::mesh::channel`]）。
    ///
    /// x 軸を流下方向（上流端 x = 0、下流端 x = `length`）、y 軸を横断方向
    /// （側壁 y = 0 と y = `width`）にとる。
    Channel {
        /// 水路延長 L \[m\]。
        length: f64,
        /// 水路幅 W \[m\]。
        width: f64,
        /// 代表格子幅（母点の間隔）\[m\]。
        cell_size: f64,
        /// 母点の配置。
        #[serde(default)]
        layout: ChannelLayout,
    },
    /// 格子生成ツール（`tools/make_mesh.py`）が書き出した格子ファイルを読み込む
    /// （[`crate::mesh::file`]）。セルの河床高は、`tools/make_bed.py` が書き出すセルの値の
    /// ファイル（[`crate::mesh::cells`]、ADR-0020）で与える。
    ///
    /// パスは、実行したディレクトリからの相対パス（`output.dir` と同じ）。
    File {
        /// 格子ファイルのパス。
        path: PathBuf,
        /// セルの値のファイルのパス（河床高、非侵食層）。[`crate::model::Model::new`] で
        /// 組み立てるときに必要（[`crate::model::Model::with_state`] では使わない）。
        #[serde(default)]
        cells: Option<PathBuf>,
    },
}

/// 直線水路の母点の配置。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ChannelLayout {
    /// 母点を三角格子状に置く。内部のセルはほぼ正六角形になる（[`crate::mesh::channel`]）。
    #[default]
    Hexagonal,
    /// 母点を正方格子状に置く。セルはほぼ正方形の長方形になる。
    Rectangular,
}

/// 初期河床（直線水路用）。
///
/// 河床高は z_b(x) = `z_downstream` + `slope` × (L − x)。床止めがある場合は
/// [`StructureConfig::Sill`] の説明のとおり、床止めより上流をかさ上げする。
///
/// 周期境界（[`UpstreamBoundary::Periodic`]）では、河床高と水位を勾配を除いた値で持つので、
/// 初期の河床高は一様に `z_downstream` になり、`slope` は運動量式の体積力 g I n_x の勾配として
/// 使う（ADR-0014）。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BedConfig {
    /// 縦断勾配 I \[-\]。下流に向かって下がる向きを正とする。
    pub slope: f64,
    /// 下流端（x = L）の河床高 \[m\]。
    #[serde(default)]
    pub z_downstream: f64,
}

/// 構造物。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum StructureConfig {
    /// 床止め。
    ///
    /// x ≤ `x_end` の初期河床を `height` だけかさ上げし、`x_start` ≤ x ≤ `x_end` の帯を
    /// 非侵食層とする（堆砂が天端まで進んだ床止めを想定）。下流端 x = `x_end` で、河床に
    /// 高さ `height` の段差ができる。
    Sill {
        /// 床止めの上流端の位置 \[m\]。
        x_start: f64,
        /// 床止めの下流端の位置 \[m\]。
        x_end: f64,
        /// 落差（天端高 − 直下流の初期河床高）\[m\]。
        height: f64,
    },
}

/// 物理定数。省略した項目は既定値になる。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct PhysicsConfig {
    /// 重力加速度 g \[m/s²\]。
    pub gravity: f64,
    /// 水の密度 ρ \[kg/m³\]。
    pub water_density: f64,
    /// 水の動粘性係数 ν \[m²/s\]。
    pub kinematic_viscosity: f64,
    /// カルマン定数 κ \[-\]（二次流の適応長に使う）。
    pub von_karman: f64,
}

impl Default for PhysicsConfig {
    fn default() -> Self {
        Self {
            gravity: 9.81,
            water_density: 1000.0,
            kinematic_viscosity: 1.0e-6,
            von_karman: 0.41,
        }
    }
}

/// 時間積分の設定。時刻はすべて流れの時間（MORFAC を掛ける前）で表す。
///
/// 時間刻みは、固定（`dt`）か可変（`courant` と `dt_max`）のどちらか一方で与える
/// （[`TimeConfig::step_control`]）。可変のときは、各ステップの初めに、セルから1ステップで
/// 流出する体積がセルの体積の `courant` 倍以下になる最大の値を、`dt_max` で頭打ちにして使う
/// （ADR-0010）。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TimeConfig {
    /// 固定の時間刻み Δt \[s\]。移流の CFL 条件で上限が決まる（重力波には制約されない）。
    /// `courant` と同時には書けない。
    pub dt: Option<f64>,
    /// 可変の時間刻みの目標のクーラン数 C（0 < C ≤ 1）\[-\]。セルの流出量の比
    /// Δt Σ_流出 L_j q_j / V_i の最大値がこの値になるように Δt を決める。`dt` と同時には書けない。
    pub courant: Option<f64>,
    /// 可変の時間刻みの上限 \[s\]。`courant` を書くときは必須（流れのない状態から始めると、
    /// クーラン数の条件からは時間刻みが決まらないため）。
    pub dt_max: Option<f64>,
    /// 計算の終了時刻 \[s\]。
    pub end: f64,
    /// 河床変動の計算を始める時刻 \[s\]。それまでは固定床で流れだけを解く（助走計算）。
    #[serde(default)]
    pub morphology_start: f64,
}

/// 時間刻みの決め方（[`TimeConfig::step_control`] が設定から作る）。
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TimeStep {
    /// 固定の時間刻み。
    Fixed {
        /// 時間刻み Δt \[s\]。
        dt: f64,
    },
    /// クーラン数で決める可変の時間刻み。
    Adaptive {
        /// 目標のクーラン数（セルの流出量の比の最大値）\[-\]。
        courant: f64,
        /// 時間刻みの上限 \[s\]。
        dt_max: f64,
    },
}

impl TimeStep {
    /// 時間刻みの上限 \[s\]（固定ならその値）。時刻の比較の許容差の目安に使う。
    pub fn max_dt(self) -> f64 {
        match self {
            Self::Fixed { dt } => dt,
            Self::Adaptive { dt_max, .. } => dt_max,
        }
    }
}

impl TimeConfig {
    /// 時間刻みの決め方。`dt` だけなら固定、`courant` と `dt_max` なら可変。それ以外の
    /// 組み合わせ（両方ある、どちらもない、`courant` に `dt_max` がない、`dt` に `dt_max` が
    /// ある）は `None`（[`Config::validate`] がエラーにする）。
    pub fn step_control(&self) -> Option<TimeStep> {
        match (self.dt, self.courant, self.dt_max) {
            (Some(dt), None, None) => Some(TimeStep::Fixed { dt }),
            (None, Some(courant), Some(dt_max)) => Some(TimeStep::Adaptive { courant, dt_max }),
            _ => None,
        }
    }
}

/// 流れの計算の設定。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HydroConfig {
    /// θ 法の重み θ（0.5 ≤ θ ≤ 1、推奨 0.55〜0.6）。
    pub theta: f64,
    /// マニングの粗度係数 n [s/m^(1/3)]。0 なら摩擦なし（解析解と比べる検証用。流砂を
    /// 計算するときは掃流力の評価に使うので正の値が必要）。
    pub manning_n: f64,
    /// wet/dry を判定する水深の閾値 \[m\]。
    #[serde(default = "default_dry_depth")]
    pub dry_depth: f64,
    /// 運動量の移流の扱い（省略時は2次精度の保存型の差分）。
    #[serde(default)]
    pub advection: AdvectionScheme,
    /// 時間積分の方法（省略時は予測子・修正子）。2次の移流には `"predictor_corrector"` が必要
    /// （ADR-0016、ADR-0017）。
    #[serde(default)]
    pub time_integration: TimeIntegration,
    /// 水位の連立方程式の反復解法。
    #[serde(default)]
    pub solver: SolverConfig,
}

fn default_dry_depth() -> f64 {
    1.0e-3
}

/// 運動量の移流の扱い。設定ファイルでは `"second_order"`（既定）、`"upwind"`、`"none"` のいずれか。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
pub enum AdvectionScheme {
    /// 保存型の1次風上差分（[`crate::hydro::advection`]。ADR-0008, ADR-0009）。
    #[serde(rename = "upwind")]
    Upwind,
    /// 保存型の2次精度の差分（[`crate::hydro::advection_second_order`]。ADR-0016）。
    /// セルの流速の勾配（最小二乗）で面の値を再構成する。`time_integration = "predictor_corrector"`
    /// が必要。
    #[default]
    #[serde(rename = "second_order")]
    SecondOrder,
    /// 移流項を省く（F u_j = u_j）。検証で移流の効果を比べるためのもので、急変流には使えない。
    #[serde(rename = "none")]
    Disabled,
}

/// 流れの時間積分の方法。設定ファイルでは `"predictor_corrector"`（既定）または `"euler"`。
///
/// どちらも、水位勾配は θ 法、摩擦は陰的に扱う。違いは陽的な部分（移流と、連続式の風上の
/// 辺の水深）の扱い（[`crate::hydro::HydroSolver::step`]）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
pub enum TimeIntegration {
    /// 時刻 n の陽的な部分で1回解く（移流は前進 Euler）。1次の移流か、移流なしと組み合わせる。
    #[serde(rename = "euler")]
    Euler,
    /// 予測子・修正子：時刻 n で1段解いた予測子 U* と時刻 n の陽的な部分を平均し（Heun 法）、
    /// 時刻 n から水位を解き直す。計算量は2倍だが、2次の移流でも安定で、定常解は Δt によらず、
    /// 重力波の減衰は θ 法と同じ（ADR-0017）。
    #[default]
    #[serde(rename = "predictor_corrector")]
    PredictorCorrector,
}

/// 水位の連立方程式（Newton 法＋PCG）の設定。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct SolverConfig {
    /// Newton 法の収束判定値：体積残差をセル面積で割った値の最大値 \[m\]。
    pub newton_tol: f64,
    /// Newton 法の最大反復回数。
    pub newton_max_iter: usize,
    /// PCG の収束判定値：相対残差 ‖r‖₂ / ‖b‖₂。
    pub pcg_rel_tol: f64,
    /// PCG の最大反復回数。
    pub pcg_max_iter: usize,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            newton_tol: 1.0e-10,
            newton_max_iter: 50,
            pcg_rel_tol: 1.0e-12,
            pcg_max_iter: 1000,
        }
    }
}

/// 流砂と河床変動の設定（均一粒径）。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SedimentConfig {
    /// 粒径 d \[m\]。
    pub grain_size: f64,
    /// 土粒子の密度 ρ_s \[kg/m³\]。
    #[serde(default = "default_sediment_density")]
    pub density: f64,
    /// 河床の空隙率 λ \[-\]。
    #[serde(default = "default_porosity")]
    pub porosity: f64,
    /// 地形変化の加速係数 MORFAC \[-\]。
    #[serde(default = "default_morfac")]
    pub morfac: f64,
    /// 水中安息角 φ \[°\]。河床勾配が tan φ を超える辺では土砂を再配分する。
    pub angle_of_repose: f64,
    /// 掃流砂。
    pub bedload: BedloadConfig,
    /// 浮遊砂。`None`（セクションを省略）のときは浮遊砂を考慮しない。
    pub suspended: Option<SuspendedConfig>,
    /// 二次流による掃流砂の向きの補正。`None`（セクションを省略）のときは補正しない
    /// （掃流砂は水深平均流速の向き）。
    pub secondary_flow: Option<SecondaryFlowConfig>,
}

fn default_sediment_density() -> f64 {
    2650.0
}

fn default_porosity() -> f64 {
    0.4
}

fn default_morfac() -> f64 {
    1.0
}

/// 掃流砂の設定。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BedloadConfig {
    /// 掃流砂量式。
    pub formula: BedloadFormula,
    /// 無次元限界掃流力 τ*_c \[-\]。
    pub critical_shields: f64,
    /// 勾配補正（長谷川の式）に使う静止摩擦係数 μ_s \[-\]。
    #[serde(default = "default_static_friction")]
    pub static_friction: f64,
    /// 勾配補正（長谷川の式）に使う動摩擦係数 μ_k \[-\]。
    #[serde(default = "default_kinetic_friction")]
    pub kinetic_friction: f64,
}

fn default_static_friction() -> f64 {
    1.0
}

fn default_kinetic_friction() -> f64 {
    0.45
}

/// 掃流砂量式。設定ファイルでは `"mpm"` または `"ashida_michiue"`。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BedloadFormula {
    /// Meyer-Peter & Müller (1948)。
    Mpm,
    /// 芦田・道上 (1972)。
    AshidaMichiue,
}

/// 浮遊砂の設定。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SuspendedConfig {
    /// 浮上量式。
    pub entrainment: EntrainmentFormula,
    /// この水深未満のセルでは浮上を止める \[m\]。
    pub min_depth: f64,
    /// 水平拡散係数 ε = α u_* h の係数 α \[-\]。
    pub diffusivity_coef: f64,
}

/// 二次流（らせん流）の設定（[`crate::sediment::secondary_flow`]、ADR-0019）。
///
/// 二次流の強さ I を求め、掃流砂の向きを内岸側へ tan δ = N* I / |u| だけ回す。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SecondaryFlowConfig {
    /// I の求め方（省略時は輸送方程式型）。
    #[serde(default)]
    pub kind: SecondaryFlowKind,
    /// 掃流砂の偏角の係数 N* \[-\]（Engelund, 1974 の 7）。
    #[serde(default = "default_n_star")]
    pub n_star: f64,
    /// I の水平拡散係数 D = α_D u_* h の係数 α_D \[-\]（輸送方程式型だけで使う）。
    #[serde(default = "default_secondary_flow_diffusivity")]
    pub diffusivity_coef: f64,
}

fn default_n_star() -> f64 {
    7.0
}

fn default_secondary_flow_diffusivity() -> f64 {
    0.15
}

/// 二次流の強さ I の求め方。設定ファイルでは `"transport"`（既定）または `"equilibrium"`。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SecondaryFlowKind {
    /// 輸送方程式型（非平衡型）：遠心加速度を生成項とし、適応長 L_I で平衡値 I_e へ緩和する
    /// I の移流拡散方程式を解く。湾曲部の下流でゆっくり減衰する二次流も表せる。
    #[default]
    Transport,
    /// 平衡型：各時刻の流れの平衡値 I = I_e を使う（輸送方程式型の L_I → 0 の極限。比較用）。
    Equilibrium,
}

/// 浮遊砂の浮上量式。設定ファイルでは `"itakura_kishi"` または `"garcia_parker"`。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EntrainmentFormula {
    /// 板倉・岸の式。
    ItakuraKishi,
    /// Garcia & Parker (1991)。
    GarciaParker,
}

/// 境界条件。上流端・下流端以外の境界は壁（法線流速 0、流砂の出入りなし）とする。
///
/// 上流端と下流端をともに `kind = "periodic"` にすると、直線水路の上流端と下流端をつないだ
/// 周期境界になる（[`UpstreamBoundary::Periodic`]）。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BoundaryConfig {
    /// 上流端（直線水路では x = 0）。
    pub upstream: UpstreamBoundary,
    /// 下流端（直線水路では x = L）。
    pub downstream: DownstreamBoundary,
}

impl BoundaryConfig {
    /// 周期境界（上流端が `periodic`）なら `true`。[`Config::validate`] は、下流端も `periodic`
    /// であることを確かめる。
    pub fn is_periodic(&self) -> bool {
        matches!(self.upstream, UpstreamBoundary::Periodic { .. })
    }
}

/// 上流端の境界条件。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum UpstreamBoundary {
    /// 流量を与える。各境界辺へは通水能 K_j = L_j H_j^(5/3) / n に比例して配分する。H_j は、
    /// 上流端の湿ったセルの水位の平均（断面で共通の水位）から内側のセルの河床高を引いた値。
    Discharge {
        /// 流量 Q \[m³/s\]。
        discharge: f64,
        /// 給砂の条件。
        #[serde(default)]
        sediment_supply: SedimentSupply,
    },
    /// 壁（閉鎖水域の検証用）。
    Wall,
    /// 周期境界：上流端を下流端とつなぐ（直線水路の格子だけ。下流端も `periodic` にする）。
    ///
    /// 河床高と水位は `[bed] slope` の勾配 I を除いた値で持ち、勾配は運動量式の体積力
    /// g I n_x として与える（ADR-0014）。給砂と流出はなく、水と土砂の量は保たれる。流量は
    /// 水の量と摩擦で決まり、河床が変わると変わる。
    Periodic {
        /// 初期の等流の流量 Q \[m³/s\]（[`crate::model::Model::new`] が初期状態を作るのに使う）。
        discharge: f64,
    },
}

/// 上流端の給砂の条件（ADR-0011、ADR-0013）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SedimentSupply {
    /// 平衡給砂：各ステップの上流端の流れから求めた掃流砂量（流砂量式の値）を供給する。
    /// 給砂量は上流端の流れに追従するので、長い時間では河床勾配が初期の状態に戻るとは
    /// 限らない。
    #[default]
    Equilibrium,
    /// 一定の給砂量：河床変動を始めた最初のステップの平衡給砂量の合計を保ち、各ステップの
    /// 平衡給砂量に比例して辺に配る。上流端の流れが時間とともに速くなる条件（構造物の水面低下が
    /// 上流端に届く場合など）では、上流の区間が給砂不足で下がり、流路が分化しうる（ADR-0013）。
    Constant,
    /// 給砂しない。
    None,
}

/// 下流端の境界条件。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum DownstreamBoundary {
    /// 水位を固定する（境界辺の上の水位）。
    WaterLevel {
        /// 水位 η \[m\]。
        level: f64,
    },
    /// 水深を一定に保つ（境界セルの河床高に `depth` を加えた値を、境界辺の上の水位とする）。
    ///
    /// 移動床（`[sediment]`）とは組み合わせられない。水位が下流端の河床とともに上下するので、
    /// 河床が下がると水面勾配が増えてさらに掘れ、上がると堰上げでさらに堆積する（正の
    /// フィードバック）。等流の平衡河床でも、下流端から河床低下が進む（ADR-0011）。
    Depth {
        /// 水深 h \[m\]。
        depth: f64,
    },
    /// 壁（閉鎖水域の検証用）。
    Wall,
    /// 周期境界（上流端も `periodic` にする。[`UpstreamBoundary::Periodic`]）。
    Periodic,
}

/// 出力の設定。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OutputConfig {
    /// 出力先のディレクトリ。
    pub dir: PathBuf,
    /// 出力間隔 \[s\]。
    pub interval: f64,
}

/// 設定の読み込み・検証の失敗。
#[derive(Debug)]
pub enum ConfigError {
    /// ファイルを読めない。
    Io {
        /// 読もうとしたファイル。
        path: PathBuf,
        /// 元のエラー。
        source: std::io::Error,
    },
    /// TOML として解釈できない、未知のキーがある、または型が合わない。
    Parse(toml::de::Error),
    /// 値が範囲外など、意味として正しくない。項目ごとのメッセージを持つ。
    Invalid(Vec<String>),
}

impl fmt::Display for ConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "cannot read {}: {source}", path.display()),
            Self::Parse(e) => write!(f, "invalid TOML: {e}"),
            Self::Invalid(messages) => {
                write!(f, "invalid configuration:")?;
                for m in messages {
                    write!(f, "\n  - {m}")?;
                }
                Ok(())
            }
        }
    }
}

impl Error for ConfigError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            Self::Parse(e) => Some(e),
            Self::Invalid(_) => None,
        }
    }
}

impl FromStr for Config {
    type Err = ConfigError;

    /// TOML の文字列から読み込み、[`Config::validate`] で検証する。
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let config: Self = toml::from_str(s).map_err(ConfigError::Parse)?;
        config.validate()?;
        Ok(config)
    }
}

impl Config {
    /// ファイルから読み込み、[`Config::validate`] で検証する。
    pub fn load(path: impl AsRef<Path>) -> Result<Self, ConfigError> {
        let path = path.as_ref();
        let text = fs::read_to_string(path).map_err(|source| ConfigError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        text.parse()
    }

    /// 値の範囲と、項目どうしの整合を確かめる。問題をすべて集めて返す。
    pub fn validate(&self) -> Result<(), ConfigError> {
        let mut errors = Vec::new();
        let mut require = |ok: bool, message: &str| {
            if !ok {
                errors.push(message.to_owned());
            }
        };

        require(!self.title.trim().is_empty(), "title must not be empty");

        match &self.mesh {
            MeshConfig::Channel {
                length,
                width,
                cell_size,
                ..
            } => {
                require(positive(*length), "mesh.length must be > 0");
                require(positive(*width), "mesh.width must be > 0");
                require(
                    positive(*cell_size) && *cell_size <= width.min(*length),
                    "mesh.cell_size must be > 0 and not larger than the channel",
                );
                require(
                    self.bed.is_some(),
                    "[bed] is required when mesh.kind = \"channel\"",
                );
                for structure in &self.structures {
                    match structure {
                        StructureConfig::Sill {
                            x_start,
                            x_end,
                            height,
                        } => {
                            require(
                                0.0 <= *x_start && x_start < x_end && *x_end <= *length,
                                "structures: a sill must satisfy 0 <= x_start < x_end <= mesh.length",
                            );
                            require(
                                height.is_finite() && *height >= 0.0,
                                "structures: sill height must be >= 0",
                            );
                        }
                    }
                }
            }
            MeshConfig::File { path, cells } => {
                require(!path.as_os_str().is_empty(), "mesh.path must not be empty");
                require(
                    cells.as_ref().is_none_or(|c| !c.as_os_str().is_empty()),
                    "mesh.cells must not be empty",
                );
                require(
                    self.structures.is_empty(),
                    "structures are supported only with mesh.kind = \"channel\" for now",
                );
            }
        }

        if let Some(bed) = &self.bed {
            require(bed.slope.is_finite(), "bed.slope must be finite");
            require(
                bed.z_downstream.is_finite(),
                "bed.z_downstream must be finite",
            );
        }

        let physics = &self.physics;
        require(positive(physics.gravity), "physics.gravity must be > 0");
        require(
            positive(physics.water_density),
            "physics.water_density must be > 0",
        );
        require(
            positive(physics.kinematic_viscosity),
            "physics.kinematic_viscosity must be > 0",
        );
        require(
            positive(physics.von_karman) && physics.von_karman < 1.0,
            "physics.von_karman must be in (0, 1)",
        );

        let time = &self.time;
        require(positive(time.end), "time.end must be > 0");
        match (time.dt, time.courant, time.dt_max) {
            (Some(_), Some(_), _) => require(
                false,
                "time: give either dt (fixed time step) or courant (adaptive time step), not both",
            ),
            (None, None, _) => require(
                false,
                "time: dt (fixed time step) or courant (adaptive time step) is required",
            ),
            (Some(dt), None, dt_max) => {
                require(
                    positive(dt) && dt <= time.end,
                    "time.dt must be > 0 and <= time.end",
                );
                require(
                    dt_max.is_none(),
                    "time.dt_max is used only with time.courant (adaptive time step)",
                );
            }
            (None, Some(courant), dt_max) => {
                require(
                    positive(courant) && courant <= 1.0,
                    "time.courant must be in (0, 1]",
                );
                match dt_max {
                    Some(dt_max) => require(
                        positive(dt_max) && dt_max <= time.end,
                        "time.dt_max must be > 0 and <= time.end",
                    ),
                    None => require(
                        false,
                        "time.dt_max is required with time.courant (the Courant condition gives no bound while nothing flows)",
                    ),
                }
            }
        }
        require(
            time.morphology_start.is_finite() && time.morphology_start >= 0.0,
            "time.morphology_start must be >= 0",
        );

        let hydro = &self.hydro;
        require(
            (0.5..=1.0).contains(&hydro.theta),
            "hydro.theta must be in [0.5, 1.0]",
        );
        require(
            hydro.manning_n.is_finite() && hydro.manning_n >= 0.0,
            "hydro.manning_n must be >= 0",
        );
        require(positive(hydro.dry_depth), "hydro.dry_depth must be > 0");
        require(
            hydro.advection != AdvectionScheme::SecondOrder
                || hydro.time_integration == TimeIntegration::PredictorCorrector,
            "hydro.advection = \"second_order\" needs hydro.time_integration = \"predictor_corrector\" (it is unstable with \"euler\")",
        );
        require(
            positive(hydro.solver.newton_tol) && hydro.solver.newton_max_iter > 0,
            "hydro.solver: newton_tol and newton_max_iter must be > 0",
        );
        require(
            positive(hydro.solver.pcg_rel_tol) && hydro.solver.pcg_max_iter > 0,
            "hydro.solver: pcg_rel_tol and pcg_max_iter must be > 0",
        );

        if let Some(sediment) = &self.sediment {
            require(
                positive(sediment.grain_size),
                "sediment.grain_size must be > 0",
            );
            require(
                sediment.density.is_finite() && sediment.density > physics.water_density,
                "sediment.density must be larger than physics.water_density",
            );
            require(
                (0.0..1.0).contains(&sediment.porosity),
                "sediment.porosity must be in [0, 1)",
            );
            require(positive(sediment.morfac), "sediment.morfac must be > 0");
            require(
                positive(hydro.manning_n),
                "hydro.manning_n must be > 0 when [sediment] is given (it sets the bed shear stress)",
            );
            require(
                sediment.angle_of_repose > 0.0 && sediment.angle_of_repose < 90.0,
                "sediment.angle_of_repose must be in (0, 90) degrees",
            );
            let bedload = &sediment.bedload;
            require(
                positive(bedload.critical_shields),
                "sediment.bedload.critical_shields must be > 0",
            );
            require(
                positive(bedload.static_friction) && positive(bedload.kinetic_friction),
                "sediment.bedload: static_friction and kinetic_friction must be > 0",
            );
            if let Some(suspended) = &sediment.suspended {
                require(
                    positive(suspended.min_depth),
                    "sediment.suspended.min_depth must be > 0",
                );
                require(
                    suspended.diffusivity_coef.is_finite() && suspended.diffusivity_coef >= 0.0,
                    "sediment.suspended.diffusivity_coef must be >= 0",
                );
            }
            if let Some(secondary) = &sediment.secondary_flow {
                require(
                    secondary.n_star.is_finite() && secondary.n_star >= 0.0,
                    "sediment.secondary_flow.n_star must be >= 0",
                );
                require(
                    secondary.diffusivity_coef.is_finite() && secondary.diffusivity_coef >= 0.0,
                    "sediment.secondary_flow.diffusivity_coef must be >= 0",
                );
            }
        }

        match &self.boundary.upstream {
            UpstreamBoundary::Discharge { discharge, .. }
            | UpstreamBoundary::Periodic { discharge } => require(
                discharge.is_finite() && *discharge >= 0.0,
                "boundary.upstream.discharge must be >= 0",
            ),
            UpstreamBoundary::Wall => {}
        }
        let downstream_periodic = matches!(self.boundary.downstream, DownstreamBoundary::Periodic);
        require(
            self.boundary.is_periodic() == downstream_periodic,
            "boundary: kind = \"periodic\" must be given at both the upstream and the downstream end",
        );
        if self.boundary.is_periodic() {
            require(
                self.bed.is_some(),
                "[bed] is required with periodic boundaries (bed.slope drives the flow as a body force)",
            );
            require(
                self.structures.is_empty(),
                "structures cannot be used with periodic boundaries",
            );
        }
        match &self.boundary.downstream {
            DownstreamBoundary::WaterLevel { level } => require(
                level.is_finite(),
                "boundary.downstream.level must be finite",
            ),
            DownstreamBoundary::Depth { depth } => {
                require(positive(*depth), "boundary.downstream.depth must be > 0");
                require(
                    self.sediment.is_none(),
                    "boundary.downstream: kind = \"depth\" cannot be used with [sediment] (the water level follows the outlet bed, so the outlet erodes or aggrades without limit); use kind = \"water_level\"",
                );
            }
            DownstreamBoundary::Wall | DownstreamBoundary::Periodic => {}
        }

        require(
            positive(self.output.interval),
            "output.interval must be > 0",
        );

        if errors.is_empty() {
            Ok(())
        } else {
            Err(ConfigError::Invalid(errors))
        }
    }

    /// 設定の要約（`vormorph check` で表示する）。
    pub fn summary(&self) -> impl fmt::Display {
        Summary(self)
    }
}

/// 有限の正の値なら true（NaN と無限大は false）。
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

/// [`Config::summary`] の表示。
struct Summary<'a>(&'a Config);

impl ChannelLayout {
    /// 設定ファイルでの表記。
    fn name(self) -> &'static str {
        match self {
            Self::Hexagonal => "hexagonal",
            Self::Rectangular => "rectangular",
        }
    }
}

impl BedloadFormula {
    /// 設定ファイルでの表記。
    fn name(self) -> &'static str {
        match self {
            Self::Mpm => "mpm",
            Self::AshidaMichiue => "ashida_michiue",
        }
    }
}

impl AdvectionScheme {
    /// 設定ファイルでの表記。
    fn name(self) -> &'static str {
        match self {
            Self::Upwind => "upwind",
            Self::SecondOrder => "second_order",
            Self::Disabled => "none",
        }
    }
}

impl TimeIntegration {
    /// 設定ファイルでの表記。
    fn name(self) -> &'static str {
        match self {
            Self::Euler => "euler",
            Self::PredictorCorrector => "predictor_corrector",
        }
    }
}

impl SecondaryFlowKind {
    /// 設定ファイルでの表記。
    fn name(self) -> &'static str {
        match self {
            Self::Transport => "transport",
            Self::Equilibrium => "equilibrium",
        }
    }
}

impl SedimentSupply {
    /// 設定ファイルでの表記。
    fn name(self) -> &'static str {
        match self {
            Self::Equilibrium => "equilibrium",
            Self::Constant => "constant",
            Self::None => "none",
        }
    }
}

impl fmt::Display for Summary<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let c = self.0;
        writeln!(f, "title      : {}", c.title)?;
        match &c.mesh {
            MeshConfig::Channel {
                length,
                width,
                cell_size,
                layout,
            } => writeln!(
                f,
                "mesh       : channel {length} m x {width} m, cell size {cell_size} m, {}",
                layout.name()
            )?,
            MeshConfig::File { path, cells } => match cells {
                Some(cells) => writeln!(
                    f,
                    "mesh       : file {}, cell values {}",
                    path.display(),
                    cells.display()
                )?,
                None => writeln!(f, "mesh       : file {}", path.display())?,
            },
        }
        if let Some(bed) = &c.bed {
            writeln!(
                f,
                "bed        : slope {}, z_downstream {} m",
                bed.slope, bed.z_downstream
            )?;
        }
        for s in &c.structures {
            match s {
                StructureConfig::Sill {
                    x_start,
                    x_end,
                    height,
                } => writeln!(
                    f,
                    "structure  : sill x = {x_start}..{x_end} m, drop {height} m"
                )?,
            }
        }
        match c.time.step_control() {
            Some(TimeStep::Fixed { dt }) => write!(f, "time       : dt {dt} s")?,
            Some(TimeStep::Adaptive { courant, dt_max }) => write!(
                f,
                "time       : adaptive dt (Courant {courant}, dt <= {dt_max} s)"
            )?,
            None => write!(f, "time       : invalid time step")?,
        }
        writeln!(
            f,
            ", end {} s, morphology from {} s",
            c.time.end, c.time.morphology_start
        )?;
        writeln!(
            f,
            "hydro      : theta {}, Manning n {}, advection {}, time integration {}",
            c.hydro.theta,
            c.hydro.manning_n,
            c.hydro.advection.name(),
            c.hydro.time_integration.name()
        )?;
        match &c.sediment {
            Some(s) => {
                write!(
                    f,
                    "sediment   : d {} m, bedload {}, MORFAC {}, suspended {}",
                    s.grain_size,
                    s.bedload.formula.name(),
                    s.morfac,
                    if s.suspended.is_some() { "on" } else { "off" }
                )?;
                match &s.secondary_flow {
                    Some(sf) => {
                        writeln!(f, ", secondary flow {} (N* {})", sf.kind.name(), sf.n_star)?
                    }
                    None => writeln!(f, ", secondary flow off")?,
                }
            }
            None => writeln!(f, "sediment   : fixed bed")?,
        }
        match &c.boundary.upstream {
            UpstreamBoundary::Discharge {
                discharge,
                sediment_supply,
            } => writeln!(
                f,
                "upstream   : discharge {discharge} m^3/s, sediment supply {}",
                sediment_supply.name()
            )?,
            UpstreamBoundary::Wall => writeln!(f, "upstream   : wall")?,
            UpstreamBoundary::Periodic { discharge } => writeln!(
                f,
                "upstream   : periodic (connected to the downstream end), initial discharge {discharge} m^3/s"
            )?,
        }
        match &c.boundary.downstream {
            DownstreamBoundary::WaterLevel { level } => {
                writeln!(f, "downstream : water level {level} m")?
            }
            DownstreamBoundary::Depth { depth } => writeln!(f, "downstream : depth {depth} m")?,
            DownstreamBoundary::Wall => writeln!(f, "downstream : wall")?,
            DownstreamBoundary::Periodic => writeln!(f, "downstream : periodic")?,
        }
        write!(
            f,
            "output     : {} every {} s",
            c.output.dir.display(),
            c.output.interval
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SILL: &str = include_str!("../examples/sill.toml");

    #[test]
    fn sill_example_is_valid() {
        let config: Config = SILL.parse().expect("examples/sill.toml must be valid");
        assert!(matches!(
            config.mesh,
            MeshConfig::Channel { length, width, .. } if length == 200.0 && width == 20.0
        ));
        assert_eq!(config.bed.as_ref().map(|b| b.slope), Some(0.002));
        assert!(matches!(
            config.boundary.upstream,
            UpstreamBoundary::Discharge { discharge, .. } if discharge == 40.0
        ));
        assert!(matches!(
            config.structures.as_slice(),
            [StructureConfig::Sill { height, .. }] if *height == 1.0
        ));
        assert_eq!(config.sediment.as_ref().map(|s| s.grain_size), Some(0.002));
    }

    #[test]
    fn unknown_key_is_rejected() {
        let text = SILL.replacen("theta =", "theta_typo =", 1);
        assert!(matches!(text.parse::<Config>(), Err(ConfigError::Parse(_))));
    }

    #[test]
    fn out_of_range_theta_is_rejected() {
        let text = SILL.replacen("theta = 0.55", "theta = 0.4", 1);
        match text.parse::<Config>() {
            Err(ConfigError::Invalid(messages)) => {
                assert!(messages.iter().any(|m| m.contains("hydro.theta")));
            }
            other => panic!("expected ConfigError::Invalid, got {other:?}"),
        }
    }

    #[test]
    fn frictionless_flow_is_allowed_only_without_sediment() {
        let text = SILL.replacen("manning_n = 0.025", "manning_n = 0.0", 1);
        match text.parse::<Config>() {
            Err(ConfigError::Invalid(messages)) => {
                assert!(messages.iter().any(|m| m.contains("hydro.manning_n")));
            }
            other => panic!("expected ConfigError::Invalid, got {other:?}"),
        }
        let start = text
            .find("\n[sediment]\n")
            .expect("sill.toml has [sediment]");
        let end = text
            .find("\n[boundary.upstream]\n")
            .expect("sill.toml has [boundary.upstream]");
        let fixed_bed = format!("{}{}", &text[..start], &text[end..]);
        let config: Config = fixed_bed.parse().expect("n = 0 is valid on a fixed bed");
        assert_eq!(config.hydro.manning_n, 0.0);
        assert!(config.sediment.is_none());

        let negative = fixed_bed.replacen("manning_n = 0.0", "manning_n = -0.01", 1);
        assert!(matches!(
            negative.parse::<Config>(),
            Err(ConfigError::Invalid(_))
        ));
    }

    #[test]
    fn advection_can_be_disabled_and_rejects_unknown_schemes() {
        const SECOND: &str = "advection = \"second_order\"";
        let config: Config = SILL.parse().expect("valid");
        assert_eq!(config.hydro.advection, AdvectionScheme::SecondOrder);
        assert_eq!(
            config.hydro.time_integration,
            TimeIntegration::PredictorCorrector
        );
        assert!(
            config
                .summary()
                .to_string()
                .contains("advection second_order, time integration predictor_corrector")
        );
        let text = SILL.replacen(SECOND, "advection = \"none\"", 1);
        let config: Config = text.parse().expect("advection = \"none\" is valid");
        assert_eq!(config.hydro.advection, AdvectionScheme::Disabled);
        let text = SILL.replacen(SECOND, "advection = \"upwind\"", 1);
        let config: Config = text
            .parse()
            .expect("upwind with predictor_corrector is valid");
        assert_eq!(config.hydro.advection, AdvectionScheme::Upwind);
        let omitted = SILL.replacen(SECOND, "", 1);
        assert_eq!(
            omitted.parse::<Config>().expect("valid").hydro.advection,
            AdvectionScheme::SecondOrder
        );
        let unknown = SILL.replacen(SECOND, "advection = \"elm\"", 1);
        assert!(matches!(
            unknown.parse::<Config>(),
            Err(ConfigError::Parse(_))
        ));
    }

    /// 時間積分は省略すると予測子・修正子。2次の移流は1段（euler）では不安定なので、組み合わせを
    /// 誤りにする。1段は1次の移流となら使える。
    #[test]
    fn second_order_advection_needs_predictor_corrector() {
        const PREDICTOR_CORRECTOR: &str = "time_integration = \"predictor_corrector\"";
        let omitted = SILL.replacen(PREDICTOR_CORRECTOR, "", 1);
        assert_eq!(
            omitted
                .parse::<Config>()
                .expect("valid")
                .hydro
                .time_integration,
            TimeIntegration::PredictorCorrector
        );
        let euler = SILL.replacen(PREDICTOR_CORRECTOR, "time_integration = \"euler\"", 1);
        match euler.parse::<Config>() {
            Err(ConfigError::Invalid(messages)) => {
                assert!(
                    messages
                        .iter()
                        .any(|m| m.contains("hydro.time_integration"))
                );
            }
            other => panic!("expected ConfigError::Invalid, got {other:?}"),
        }
        let upwind = euler.replacen("advection = \"second_order\"", "advection = \"upwind\"", 1);
        let config: Config = upwind.parse().expect("upwind with euler is valid");
        assert_eq!(config.hydro.time_integration, TimeIntegration::Euler);
        assert!(
            config
                .summary()
                .to_string()
                .contains("advection upwind, time integration euler")
        );
        let unknown = SILL.replacen(PREDICTOR_CORRECTOR, "time_integration = \"rk4\"", 1);
        assert!(matches!(
            unknown.parse::<Config>(),
            Err(ConfigError::Parse(_))
        ));
    }

    /// `[time]` の時間刻みの行を `lines` に置き換えた設定。
    fn with_time_step(lines: &str) -> Result<Config, ConfigError> {
        let start = SILL.find("\n[time]\n").expect("sill.toml has [time]");
        let end = SILL[start..]
            .find("end =")
            .map(|i| start + i)
            .expect("sill.toml has time.end");
        format!("{}[time]\n{lines}\n{}", &SILL[..start], &SILL[end..]).parse()
    }

    fn invalid_messages(result: Result<Config, ConfigError>) -> Vec<String> {
        match result {
            Err(ConfigError::Invalid(messages)) => messages,
            other => panic!("expected ConfigError::Invalid, got {other:?}"),
        }
    }

    #[test]
    fn time_step_is_either_fixed_or_adaptive() {
        let fixed = with_time_step("dt = 0.1").expect("a fixed dt is valid");
        assert_eq!(fixed.time.step_control(), Some(TimeStep::Fixed { dt: 0.1 }));
        let adaptive =
            with_time_step("courant = 0.5\ndt_max = 0.2").expect("an adaptive dt is valid");
        assert_eq!(
            adaptive.time.step_control(),
            Some(TimeStep::Adaptive {
                courant: 0.5,
                dt_max: 0.2
            })
        );
        assert!(
            adaptive
                .summary()
                .to_string()
                .contains("adaptive dt (Courant 0.5, dt <= 0.2 s)")
        );

        for (lines, expected) in [
            ("dt = 0.1\ncourant = 0.5\ndt_max = 0.2", "not both"),
            ("", "is required"),
            ("courant = 0.5", "time.dt_max is required"),
            ("dt = 0.1\ndt_max = 0.2", "used only with time.courant"),
            ("courant = 1.5\ndt_max = 0.2", "time.courant"),
            ("courant = 0.5\ndt_max = 0.0", "time.dt_max"),
            ("dt = -0.1", "time.dt"),
        ] {
            let messages = invalid_messages(with_time_step(lines));
            assert!(
                messages.iter().any(|m| m.contains(expected)),
                "{lines:?}: {messages:?}"
            );
        }
    }

    #[test]
    fn depth_boundary_is_rejected_with_a_movable_bed() {
        let text = SILL.replacen(
            "kind = \"water_level\"\nlevel = 1.07",
            "kind = \"depth\"\ndepth = 1.07",
            1,
        );
        assert_ne!(text, SILL, "sill.toml has a water level boundary");
        let messages = invalid_messages(text.parse());
        assert!(
            messages
                .iter()
                .any(|m| m.contains("cannot be used with [sediment]")),
            "{messages:?}"
        );
        let start = text
            .find("\n[sediment]\n")
            .expect("sill.toml has [sediment]");
        let end = text
            .find("\n[boundary.upstream]\n")
            .expect("sill.toml has [boundary.upstream]");
        let fixed_bed = format!("{}{}", &text[..start], &text[end..]);
        assert!(
            fixed_bed.parse::<Config>().is_ok(),
            "a fixed bed may use it"
        );
    }

    #[test]
    fn physics_section_can_be_omitted() {
        let start = SILL.find("\n[physics]\n").expect("sill.toml has [physics]");
        let end = SILL.find("\n[time]\n").expect("sill.toml has [time]");
        let text = format!("{}{}", &SILL[..start], &SILL[end..]);
        let config: Config = text.parse().expect("valid without [physics]");
        assert_eq!(config.physics, PhysicsConfig::default());
    }

    /// `[sediment.secondary_flow]` は省略できる（補正なし）。節だけ書けば輸送方程式型・N* = 7・
    /// α_D = 0.15 になる。
    #[test]
    fn secondary_flow_section_has_defaults() {
        let config: Config = SILL.parse().expect("valid");
        let sediment = config.sediment.as_ref().expect("sill.toml has [sediment]");
        assert_eq!(sediment.secondary_flow, None);
        assert!(config.summary().to_string().contains("secondary flow off"));
        assert_eq!(config.physics.von_karman, 0.41);

        let with = |lines: &str| {
            SILL.replacen(
                "\n[boundary.upstream]\n",
                &format!("\n[sediment.secondary_flow]\n{lines}\n\n[boundary.upstream]\n"),
                1,
            )
            .parse::<Config>()
        };
        let config = with("").expect("an empty section is valid");
        let secondary = config
            .sediment
            .as_ref()
            .and_then(|s| s.secondary_flow.clone())
            .expect("secondary flow is on");
        assert_eq!(
            secondary,
            SecondaryFlowConfig {
                kind: SecondaryFlowKind::Transport,
                n_star: 7.0,
                diffusivity_coef: 0.15,
            }
        );
        assert!(
            config
                .summary()
                .to_string()
                .contains("secondary flow transport (N* 7)")
        );
        let equilibrium = with("kind = \"equilibrium\"\nn_star = 11.0").expect("valid");
        let secondary = equilibrium
            .sediment
            .and_then(|s| s.secondary_flow)
            .expect("secondary flow is on");
        assert_eq!(secondary.kind, SecondaryFlowKind::Equilibrium);
        assert_eq!(secondary.n_star, 11.0);

        let messages = invalid_messages(with("n_star = -1.0\ndiffusivity_coef = -0.1"));
        assert!(
            messages.iter().any(|m| m.contains("n_star")),
            "{messages:?}"
        );
        assert!(
            messages.iter().any(|m| m.contains("diffusivity_coef")),
            "{messages:?}"
        );
        assert!(matches!(
            with("kind = \"engelund\""),
            Err(ConfigError::Parse(_))
        ));
        assert!(matches!(with("nstar = 7.0"), Err(ConfigError::Parse(_))));
        let kappa = SILL.replacen("von_karman = 0.41", "von_karman = 0.0", 1);
        assert_ne!(kappa, SILL, "sill.toml has von_karman");
        let messages = invalid_messages(kappa.parse());
        assert!(
            messages.iter().any(|m| m.contains("physics.von_karman")),
            "{messages:?}"
        );
    }

    #[test]
    fn summary_uses_config_spellings() {
        let config: Config = SILL.parse().expect("examples/sill.toml must be valid");
        let summary = config.summary().to_string();
        assert!(summary.contains("hexagonal"), "{summary}");
        assert!(summary.contains("bedload mpm"), "{summary}");
        assert!(summary.contains("advection second_order"), "{summary}");
        assert!(
            summary.contains("upstream   : discharge 40 m^3/s"),
            "{summary}"
        );
        assert!(
            summary.contains("downstream : water level 1.07 m"),
            "{summary}"
        );
        assert!(!summary.contains('{'), "Debug output leaked: {summary}");
    }
}
