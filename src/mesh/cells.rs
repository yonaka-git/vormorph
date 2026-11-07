//! セルの値のファイル（vormorph cells 形式 v1、TOML）の読み込み。
//!
//! 格子ファイル（[`super::file`]）にはセルの河床高がないので、河床高などのセルの値は別のファイルで
//! 与える（ADR-0020）。`tools/make_bed.py` が、格子ファイルの母点で河床の定義を評価して書き出す。
//! 形式の仕様は `docs/design/mesh-generation.md` の第7章にある。
//!
//! ```toml
//! format = "vormorph-cells"
//! version = 1
//! title = "bend180"           # 格子の名前（照合には使わない）
//!
//! [meta]                      # 任意の付加情報。読み込み側は解釈しない
//!
//! [cells]                     # 格子ファイルと同じセルの順
//! x = [ ... ]                 # 母点の座標 [m]（格子との照合に使う）
//! y = [ ... ]
//! zb = [ ... ]                # 河床高 [m]
//! z_ne = [ ... ]              # 非侵食層の上面の高さ [m]（任意。なければ非侵食層なし）
//! ```
//!
//! - 読み込み時に、形式と版、配列の長さ、値が有限であること、z_ne ≤ z_b を検査する
//!   （[`CellValuesFile::validate`]）。
//! - 格子と組み合わせるときに、セル数と母点の座標が一致することを確かめる
//!   （[`CellValuesFile::check_mesh`]。許容差は母点の広がりの 10⁻⁹ 倍）。格子を作り直したのに
//!   セルの値のファイルを作り直さなかった、という取り違えを防ぐため。

use super::Mesh;
use serde::Deserialize;
use std::error::Error;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// `format` の値。
pub const FORMAT: &str = "vormorph-cells";

/// 読み込める形式の版。
pub const VERSION: u32 = 1;

/// 母点の座標の照合の許容差（母点の広がりに対する比）。
const COORDINATE_TOLERANCE: f64 = 1.0e-9;

/// セルの値のファイルの内容。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CellValuesFile {
    /// 形式の名前（[`FORMAT`]）。
    pub format: String,
    /// 形式の版（[`VERSION`]）。
    pub version: u32,
    /// 格子の名前。
    pub title: String,
    /// 任意の付加情報（生成条件など）。
    #[serde(default)]
    pub meta: toml::Table,
    /// セルの値。
    pub cells: CellValuesSection,
}

/// `[cells]`。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CellValuesSection {
    /// 母点の x 座標 \[m\]。
    pub x: Vec<f64>,
    /// 母点の y 座標 \[m\]。
    pub y: Vec<f64>,
    /// 河床高 \[m\]。
    pub zb: Vec<f64>,
    /// 非侵食層の上面の高さ \[m\]（なければ非侵食層なし）。
    #[serde(default)]
    pub z_ne: Option<Vec<f64>>,
}

/// セルの値のファイルの読み込みの誤り。
#[derive(Debug)]
pub enum CellValuesError {
    /// ファイルを読めない。
    Io {
        /// ファイルのパス。
        path: PathBuf,
        /// 元の誤り。
        source: std::io::Error,
    },
    /// TOML として解釈できない、または項目が合わない。
    Parse(toml::de::Error),
    /// 値が規約を満たさない、または格子と合わない（すべての問題を列挙する）。
    Invalid(Vec<String>),
}

impl fmt::Display for CellValuesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "cannot read {}: {source}", path.display()),
            Self::Parse(err) => write!(f, "invalid cell values file: {err}"),
            Self::Invalid(problems) => {
                write!(f, "invalid cell values file:")?;
                for p in problems {
                    write!(f, "\n  - {p}")?;
                }
                Ok(())
            }
        }
    }
}

impl Error for CellValuesError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            Self::Parse(err) => Some(err),
            Self::Invalid(_) => None,
        }
    }
}

impl FromStr for CellValuesFile {
    type Err = CellValuesError;

    /// TOML の文字列から読み込み、[`CellValuesFile::validate`] で検査する。
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let file: Self = toml::from_str(s).map_err(CellValuesError::Parse)?;
        file.validate()?;
        Ok(file)
    }
}

impl CellValuesFile {
    /// ファイルから読み込む。
    pub fn load(path: &Path) -> Result<Self, CellValuesError> {
        let text = fs::read_to_string(path).map_err(|source| CellValuesError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        text.parse()
    }

    /// セル数。
    pub fn n_cells(&self) -> usize {
        self.cells.x.len()
    }

    /// 形式・版、配列の長さ、値が有限であること、z_ne ≤ z_b を検査する。
    pub fn validate(&self) -> Result<(), CellValuesError> {
        let mut problems = Vec::new();
        if self.format != FORMAT {
            problems.push(format!(
                "format must be \"{FORMAT}\", got \"{}\"",
                self.format
            ));
        }
        if self.version != VERSION {
            problems.push(format!(
                "version {} is not supported (expected {VERSION})",
                self.version
            ));
        }
        let n = self.n_cells();
        let cells = &self.cells;
        let mut arrays = vec![("y", &cells.y), ("zb", &cells.zb)];
        if let Some(z_ne) = &cells.z_ne {
            arrays.push(("z_ne", z_ne));
        }
        for (name, values) in &arrays {
            if values.len() != n {
                problems.push(format!(
                    "cells.{name} has {} values but cells.x has {n}",
                    values.len()
                ));
            }
        }
        arrays.push(("x", &cells.x));
        for (name, values) in &arrays {
            if let Some(k) = values.iter().position(|v| !v.is_finite()) {
                problems.push(format!("cells.{name}[{k}] = {} is not finite", values[k]));
            }
        }
        if let Some(z_ne) = &cells.z_ne {
            if let Some(k) =
                (0..n.min(z_ne.len()).min(cells.zb.len())).find(|&k| z_ne[k] > cells.zb[k])
            {
                problems.push(format!(
                    "cells.z_ne[{k}] = {} is above cells.zb[{k}] = {}",
                    z_ne[k], cells.zb[k]
                ));
            }
        }
        if problems.is_empty() {
            Ok(())
        } else {
            Err(CellValuesError::Invalid(problems))
        }
    }

    /// セル数と母点の座標が格子 `mesh` と一致することを確かめる。
    ///
    /// 許容差は、格子の母点の広がり（x と y の範囲の大きい方、1 m 未満なら 1 m）の 10⁻⁹ 倍。
    /// 一致しない母点は、最初の1つを報告する。
    pub fn check_mesh(&self, mesh: &Mesh) -> Result<(), CellValuesError> {
        let n = mesh.n_cells();
        if self.n_cells() != n {
            return Err(CellValuesError::Invalid(vec![format!(
                "the file has {} cells but the mesh has {n}",
                self.n_cells()
            )]));
        }
        let span = |v: &[f64]| {
            let (lo, hi) = v
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &x| {
                    (lo.min(x), hi.max(x))
                });
            if n == 0 { 0.0 } else { hi - lo }
        };
        let tolerance = COORDINATE_TOLERANCE * span(&mesh.x_cell).max(span(&mesh.y_cell)).max(1.0);
        let mismatch = (0..n).find(|&c| {
            (self.cells.x[c] - mesh.x_cell[c]).abs() > tolerance
                || (self.cells.y[c] - mesh.y_cell[c]).abs() > tolerance
        });
        match mismatch {
            None => Ok(()),
            Some(c) => Err(CellValuesError::Invalid(vec![format!(
                "cell {c} is at ({}, {}) in the file but at ({}, {}) in the mesh \
                 (was the mesh regenerated without the cell values?)",
                self.cells.x[c], self.cells.y[c], mesh.x_cell[c], mesh.y_cell[c]
            )])),
        }
    }

    /// 非侵食層の上面の高さ（ファイルになければ、すべて −∞）。
    pub fn non_erodible(&self) -> Vec<f64> {
        self.cells
            .z_ne
            .clone()
            .unwrap_or_else(|| vec![f64::NEG_INFINITY; self.n_cells()])
    }
}

/// セルの値のファイルを読み込み、格子 `mesh` と照合する。
pub fn read(path: &Path, mesh: &Mesh) -> Result<CellValuesFile, CellValuesError> {
    let file = CellValuesFile::load(path)?;
    file.check_mesh(mesh)?;
    Ok(file)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::ChannelLayout;
    use crate::mesh::channel::{self, ChannelSpec};

    /// 2セル（母点 (0.5, 0.5), (1.5, 0.5)）の値。
    const TWO_CELLS: &str = r#"
format = "vormorph-cells"
version = 1
title = "two cells"

[meta]
generator = "hand-written"

[cells]
x = [0.5, 1.5]
y = [0.5, 0.5]
zb = [0.3, 0.2]
"#;

    fn two_squares() -> Mesh {
        channel::build(&ChannelSpec {
            length: 2.0,
            width: 1.0,
            cell_size: 1.0,
            layout: ChannelLayout::Rectangular,
        })
        .expect("valid channel")
    }

    #[test]
    fn values_are_read_and_matched_with_the_mesh() {
        let file: CellValuesFile = TWO_CELLS.parse().expect("valid cell values");
        assert_eq!(file.title, "two cells");
        assert_eq!(file.cells.zb, [0.3, 0.2]);
        assert_eq!(file.non_erodible(), [f64::NEG_INFINITY; 2]);
        let mesh = two_squares();
        file.check_mesh(&mesh).expect("the generators match");

        let with_layer = TWO_CELLS.replace("zb = [0.3, 0.2]", "zb = [0.3, 0.2]\nz_ne = [0.1, 0.2]");
        let file: CellValuesFile = with_layer.parse().expect("valid cell values");
        assert_eq!(file.non_erodible(), [0.1, 0.2]);
    }

    #[test]
    fn inconsistent_values_are_reported() {
        let cases = [
            ("vormorph-cells\"", "other\"", "format must be"),
            ("version = 1", "version = 2", "version 2"),
            ("zb = [0.3, 0.2]", "zb = [0.3]", "cells.zb has 1 values"),
            ("zb = [0.3, 0.2]", "zb = [0.3, nan]", "cells.zb[1] = NaN"),
            (
                "zb = [0.3, 0.2]",
                "zb = [0.3, 0.2]\nz_ne = [0.1, 0.25]",
                "z_ne[1] = 0.25 is above",
            ),
        ];
        for (from, to, expected) in cases {
            let text = TWO_CELLS.replacen(from, to, 1);
            match text.parse::<CellValuesFile>() {
                Err(CellValuesError::Invalid(problems)) => assert!(
                    problems.iter().any(|p| p.contains(expected)),
                    "{to}: {problems:?}"
                ),
                other => panic!("{to}: expected Invalid, got {other:?}"),
            }
        }
        let typo = TWO_CELLS.replace("zb =", "z_b =");
        assert!(matches!(
            typo.parse::<CellValuesFile>(),
            Err(CellValuesError::Parse(_))
        ));
    }

    #[test]
    fn a_different_mesh_is_rejected() {
        let file: CellValuesFile = TWO_CELLS.parse().expect("valid cell values");
        let moved: CellValuesFile = TWO_CELLS
            .replace("x = [0.5, 1.5]", "x = [0.5, 1.6]")
            .parse()
            .expect("valid cell values");
        let err = moved.check_mesh(&two_squares()).unwrap_err();
        assert!(err.to_string().contains("cell 1 is at (1.6, 0.5)"), "{err}");
        let larger = channel::build(&ChannelSpec {
            length: 3.0,
            width: 1.0,
            cell_size: 1.0,
            layout: ChannelLayout::Rectangular,
        })
        .expect("valid channel");
        let err = file.check_mesh(&larger).unwrap_err();
        assert!(
            err.to_string().contains("2 cells but the mesh has 3"),
            "{err}"
        );
    }
}
