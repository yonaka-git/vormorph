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
//! | `[bed]` | [`BedConfig`] | 初期河床（直線水路のとき必須） |
//! | `[[structures]]` | [`StructureConfig`] | 床止めなどの構造物（任意、複数可） |
//! | `[physics]` | [`PhysicsConfig`] | 物理定数（省略時は既定値） |
//! | `[time]` | [`TimeConfig`] | 時間刻み・終了時刻・河床変動の開始時刻 |
//! | `[hydro]` | [`HydroConfig`] | θ、粗度係数、wet/dry の閾値、反復解法 |
//! | `[sediment]` | [`SedimentConfig`] | 粒径・空隙率・流砂量式・MORFAC など（省略時は固定床） |
//! | `[boundary]` | [`BoundaryConfig`] | 上流端・下流端の境界条件（それ以外の境界は壁） |
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
    /// （[`crate::mesh::file`]）。
    File {
        /// 格子ファイルのパス。
        path: PathBuf,
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
}

impl Default for PhysicsConfig {
    fn default() -> Self {
        Self {
            gravity: 9.81,
            water_density: 1000.0,
            kinematic_viscosity: 1.0e-6,
        }
    }
}

/// 時間積分の設定。時刻はすべて流れの時間（MORFAC を掛ける前）で表す。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TimeConfig {
    /// 時間刻み Δt \[s\]。移流の CFL 条件で上限が決まる（重力波には制約されない）。
    pub dt: f64,
    /// 計算の終了時刻 \[s\]。
    pub end: f64,
    /// 河床変動の計算を始める時刻 \[s\]。それまでは固定床で流れだけを解く（助走計算）。
    #[serde(default)]
    pub morphology_start: f64,
}

/// 流れの計算の設定。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HydroConfig {
    /// θ 法の重み θ（0.5 ≤ θ ≤ 1、推奨 0.55〜0.6）。
    pub theta: f64,
    /// マニングの粗度係数 n [s/m^(1/3)]。
    pub manning_n: f64,
    /// wet/dry を判定する水深の閾値 \[m\]。
    #[serde(default = "default_dry_depth")]
    pub dry_depth: f64,
    /// 水位の連立方程式の反復解法。
    #[serde(default)]
    pub solver: SolverConfig,
}

fn default_dry_depth() -> f64 {
    1.0e-3
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
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BoundaryConfig {
    /// 上流端（直線水路では x = 0）。
    pub upstream: UpstreamBoundary,
    /// 下流端（直線水路では x = L）。
    pub downstream: DownstreamBoundary,
}

/// 上流端の境界条件。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum UpstreamBoundary {
    /// 流量を与える。各境界辺へは通水能 K_j = L_j H_j^(5/3) / n に比例して配分する。
    Discharge {
        /// 流量 Q \[m³/s\]。
        discharge: f64,
        /// 給砂の条件。
        #[serde(default)]
        sediment_supply: SedimentSupply,
    },
    /// 壁（閉鎖水域の検証用）。
    Wall,
}

/// 上流端の給砂の条件。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SedimentSupply {
    /// 平衡給砂：境界セルの流砂量と同じ量を供給する。
    #[default]
    Equilibrium,
    /// 給砂しない。
    None,
}

/// 下流端の境界条件。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum DownstreamBoundary {
    /// 水位を固定する。
    WaterLevel {
        /// 水位 η \[m\]。
        level: f64,
    },
    /// 水深を一定に保つ（境界セルの河床高に `depth` を加えた値を水位とする）。
    Depth {
        /// 水深 h \[m\]。
        depth: f64,
    },
    /// 壁（閉鎖水域の検証用）。
    Wall,
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
            MeshConfig::File { path } => {
                require(!path.as_os_str().is_empty(), "mesh.path must not be empty");
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

        let time = &self.time;
        require(positive(time.dt), "time.dt must be > 0");
        require(
            positive(time.end) && time.dt <= time.end,
            "time.end must be > 0 and >= time.dt",
        );
        require(
            time.morphology_start.is_finite() && time.morphology_start >= 0.0,
            "time.morphology_start must be >= 0",
        );

        let hydro = &self.hydro;
        require(
            (0.5..=1.0).contains(&hydro.theta),
            "hydro.theta must be in [0.5, 1.0]",
        );
        require(positive(hydro.manning_n), "hydro.manning_n must be > 0");
        require(positive(hydro.dry_depth), "hydro.dry_depth must be > 0");
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
        }

        match &self.boundary.upstream {
            UpstreamBoundary::Discharge { discharge, .. } => require(
                discharge.is_finite() && *discharge >= 0.0,
                "boundary.upstream.discharge must be >= 0",
            ),
            UpstreamBoundary::Wall => {}
        }
        match &self.boundary.downstream {
            DownstreamBoundary::WaterLevel { level } => require(
                level.is_finite(),
                "boundary.downstream.level must be finite",
            ),
            DownstreamBoundary::Depth { depth } => {
                require(positive(*depth), "boundary.downstream.depth must be > 0")
            }
            DownstreamBoundary::Wall => {}
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

impl SedimentSupply {
    /// 設定ファイルでの表記。
    fn name(self) -> &'static str {
        match self {
            Self::Equilibrium => "equilibrium",
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
            MeshConfig::File { path } => writeln!(f, "mesh       : file {}", path.display())?,
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
        writeln!(
            f,
            "time       : dt {} s, end {} s, morphology from {} s",
            c.time.dt, c.time.end, c.time.morphology_start
        )?;
        writeln!(
            f,
            "hydro      : theta {}, Manning n {}",
            c.hydro.theta, c.hydro.manning_n
        )?;
        match &c.sediment {
            Some(s) => writeln!(
                f,
                "sediment   : d {} m, bedload {}, MORFAC {}, suspended {}",
                s.grain_size,
                s.bedload.formula.name(),
                s.morfac,
                if s.suspended.is_some() { "on" } else { "off" }
            )?,
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
        }
        match &c.boundary.downstream {
            DownstreamBoundary::WaterLevel { level } => {
                writeln!(f, "downstream : water level {level} m")?
            }
            DownstreamBoundary::Depth { depth } => writeln!(f, "downstream : depth {depth} m")?,
            DownstreamBoundary::Wall => writeln!(f, "downstream : wall")?,
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
    fn physics_section_can_be_omitted() {
        let start = SILL.find("[physics]").expect("sill.toml has [physics]");
        let end = SILL.find("[time]").expect("sill.toml has [time]");
        let text = format!("{}{}", &SILL[..start], &SILL[end..]);
        let config: Config = text.parse().expect("valid without [physics]");
        assert_eq!(config.physics, PhysicsConfig::default());
    }

    #[test]
    fn summary_uses_config_spellings() {
        let config: Config = SILL.parse().expect("examples/sill.toml must be valid");
        let summary = config.summary().to_string();
        assert!(summary.contains("hexagonal"), "{summary}");
        assert!(summary.contains("bedload mpm"), "{summary}");
        assert!(
            summary.contains("upstream   : discharge 40 m^3/s"),
            "{summary}"
        );
        assert!(summary.contains("downstream : depth 1.07 m"), "{summary}");
        assert!(!summary.contains('{'), "Debug output leaked: {summary}");
    }
}
