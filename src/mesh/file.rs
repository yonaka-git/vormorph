//! 格子ファイル（vormorph mesh 形式 v1、TOML）の読み込み。
//!
//! 格子生成ツール（`tools/make_mesh.py`）が書き出すファイルを読み、[`PolygonMesh`] にする。
//! 形式の仕様は `docs/design/mesh-generation.md` にある。概略は次のとおり。
//!
//! ```toml
//! format = "vormorph-mesh"
//! version = 1
//! title = "bend180"
//! cell_order = "hilbert"      # セルの並び（hilbert / unspecified）
//!
//! [meta]                      # 任意の付加情報。読み込み側は解釈しない
//!
//! [cells]
//! x = [ ... ]                 # 母点の座標 [m]
//! y = [ ... ]
//! vertices = [[0, 1, 2, 3, 4, 5], ...]   # セルの頂点番号（反時計回り）
//!
//! [vertices]
//! x = [ ... ]                 # ボロノイ頂点の座標 [m]
//! y = [ ... ]
//!
//! [boundary]                  # 境界辺の両端の頂点番号（セルから見て反時計回り）
//! wall = [[a, b], ...]
//! upstream = [[a, b], ...]
//! downstream = [[a, b], ...]
//! ```
//!
//! 読み込み時には、配列の長さと番号の範囲だけを検査する。辺の共有や直交性などの幾何の検査は、
//! 接続配列を組み立てる [`voronoi::build`] と
//! [`Mesh::check_geometry`](super::Mesh::check_geometry) で行う。

use super::voronoi::{self, BoundaryEdge, PolygonMesh};
use super::{EdgeKind, Jagged, Mesh};
use serde::Deserialize;
use std::error::Error;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// `format` の値。
pub const FORMAT: &str = "vormorph-mesh";

/// 読み込める形式の版。
pub const VERSION: u32 = 1;

/// セルの並び順。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CellOrder {
    /// 母点のヒルベルト曲線上の位置の順。
    Hilbert,
    /// 特に定めない。
    Unspecified,
}

/// 格子ファイルの内容。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MeshFile {
    /// 形式の名前（[`FORMAT`]）。
    pub format: String,
    /// 形式の版（[`VERSION`]）。
    pub version: u32,
    /// 格子の名前。
    pub title: String,
    /// セルの並び順。
    pub cell_order: CellOrder,
    /// 任意の付加情報（生成条件など）。
    #[serde(default)]
    pub meta: toml::Table,
    /// セル。
    pub cells: CellsSection,
    /// ボロノイ頂点。
    pub vertices: VerticesSection,
    /// 境界辺。
    #[serde(default)]
    pub boundary: BoundarySection,
}

/// `[cells]`。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CellsSection {
    /// 母点の x 座標 \[m\]。
    pub x: Vec<f64>,
    /// 母点の y 座標 \[m\]。
    pub y: Vec<f64>,
    /// セルの頂点番号（反時計回り）。
    pub vertices: Vec<Vec<usize>>,
}

/// `[vertices]`。
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VerticesSection {
    /// x 座標 \[m\]。
    pub x: Vec<f64>,
    /// y 座標 \[m\]。
    pub y: Vec<f64>,
}

/// `[boundary]`。種類ごとの境界辺の両端の頂点番号。
#[derive(Debug, Clone, PartialEq, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BoundarySection {
    /// 壁。
    #[serde(default)]
    pub wall: Vec<[usize; 2]>,
    /// 上流端。
    #[serde(default)]
    pub upstream: Vec<[usize; 2]>,
    /// 下流端。
    #[serde(default)]
    pub downstream: Vec<[usize; 2]>,
}

/// 格子ファイルの読み込みの誤り。
#[derive(Debug)]
pub enum MeshFileError {
    /// ファイルを読めない。
    Io {
        /// ファイルのパス。
        path: PathBuf,
        /// 元の誤り。
        source: std::io::Error,
    },
    /// TOML として解釈できない、または項目が合わない。
    Parse(toml::de::Error),
    /// 値が規約を満たさない（すべての問題を列挙する）。
    Invalid(Vec<String>),
}

impl fmt::Display for MeshFileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io { path, source } => write!(f, "cannot read {}: {source}", path.display()),
            Self::Parse(err) => write!(f, "invalid mesh file: {err}"),
            Self::Invalid(problems) => {
                write!(f, "invalid mesh file:")?;
                for p in problems {
                    write!(f, "\n  - {p}")?;
                }
                Ok(())
            }
        }
    }
}

impl Error for MeshFileError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            Self::Parse(err) => Some(err),
            Self::Invalid(_) => None,
        }
    }
}

impl FromStr for MeshFile {
    type Err = MeshFileError;

    /// TOML の文字列から読み込み、[`MeshFile::validate`] で検査する。
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let file: Self = toml::from_str(s).map_err(MeshFileError::Parse)?;
        file.validate()?;
        Ok(file)
    }
}

impl MeshFile {
    /// ファイルから読み込む。
    pub fn load(path: &Path) -> Result<Self, MeshFileError> {
        let text = fs::read_to_string(path).map_err(|source| MeshFileError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        text.parse()
    }

    /// 形式・版、配列の長さ、番号の範囲を検査する。
    pub fn validate(&self) -> Result<(), MeshFileError> {
        let mut problems = Vec::new();
        let mut require = |ok: bool, msg: String| {
            if !ok {
                problems.push(msg);
            }
        };
        require(
            self.format == FORMAT,
            format!("format must be \"{FORMAT}\", got \"{}\"", self.format),
        );
        require(
            self.version == VERSION,
            format!(
                "version {} is not supported (expected {VERSION})",
                self.version
            ),
        );

        let n_cells = self.cells.x.len();
        let n_vertices = self.vertices.x.len();
        require(n_cells > 0, "the mesh has no cells".to_string());
        require(
            self.cells.y.len() == n_cells && self.cells.vertices.len() == n_cells,
            format!(
                "cells.x, cells.y and cells.vertices must have the same length ({}, {}, {})",
                n_cells,
                self.cells.y.len(),
                self.cells.vertices.len()
            ),
        );
        require(
            self.vertices.y.len() == n_vertices,
            format!(
                "vertices.x and vertices.y must have the same length ({}, {})",
                n_vertices,
                self.vertices.y.len()
            ),
        );
        let coords = self.cells.x.iter().chain(&self.cells.y);
        let coords = coords.chain(&self.vertices.x).chain(&self.vertices.y);
        require(
            coords.into_iter().all(|v| v.is_finite()),
            "coordinates must be finite".to_string(),
        );
        for (c, verts) in self.cells.vertices.iter().enumerate() {
            require(
                verts.len() >= 3,
                format!(
                    "cell {c} has {} vertices (at least 3 are needed)",
                    verts.len()
                ),
            );
            require(
                verts.iter().all(|&v| v < n_vertices),
                format!("cell {c} refers to a vertex out of range (n_vertices = {n_vertices})"),
            );
            let mut sorted = verts.clone();
            sorted.sort_unstable();
            sorted.dedup();
            require(
                sorted.len() == verts.len(),
                format!("cell {c} lists a vertex twice"),
            );
        }
        for (name, pairs) in self.boundary.by_kind() {
            for (k, &[a, b]) in pairs.iter().enumerate() {
                require(
                    a < n_vertices && b < n_vertices && a != b,
                    format!("boundary.{name}[{k}] = [{a}, {b}] is not a pair of distinct vertices"),
                );
            }
        }

        if problems.is_empty() {
            Ok(())
        } else {
            Err(MeshFileError::Invalid(problems))
        }
    }

    /// 境界辺の数（種類の合計）。
    pub fn n_boundary_edges(&self) -> usize {
        self.boundary.by_kind().iter().map(|(_, p)| p.len()).sum()
    }

    /// 多角形の並びに変換する。
    pub fn to_polygons(&self) -> PolygonMesh {
        let boundary_edges = self
            .boundary
            .by_kind()
            .into_iter()
            .flat_map(|(name, pairs)| {
                let kind = match name {
                    "wall" => EdgeKind::Wall,
                    "upstream" => EdgeKind::Upstream,
                    _ => EdgeKind::Downstream,
                };
                pairs
                    .iter()
                    .map(move |&vertices| BoundaryEdge { vertices, kind })
            })
            .collect();
        PolygonMesh {
            x_cell: self.cells.x.clone(),
            y_cell: self.cells.y.clone(),
            x_vertex: self.vertices.x.clone(),
            y_vertex: self.vertices.y.clone(),
            vertices_on_cell: Jagged::from_rows(self.cells.vertices.iter().cloned()),
            boundary_edges,
        }
    }
}

/// 格子ファイルを読み込み、接続配列と幾何量を組み立てる。
pub fn read(path: &Path) -> Result<Mesh, MeshFileError> {
    let file = MeshFile::load(path)?;
    Ok(voronoi::build(&file.to_polygons()))
}

impl BoundarySection {
    /// 種類の名前と境界辺の組。
    pub fn by_kind(&self) -> [(&'static str, &[[usize; 2]]); 3] {
        [
            ("wall", &self.wall),
            ("upstream", &self.upstream),
            ("downstream", &self.downstream),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 2つの正方形のセル（x = 0..1 と 1..2、y = 0..1）。左端が上流、右端が下流。
    const TWO_SQUARES: &str = r#"
format = "vormorph-mesh"
version = 1
title = "two squares"
cell_order = "unspecified"

[meta]
note = "hand-written test case"

[cells]
x = [0.5, 1.5]
y = [0.5, 0.5]
vertices = [[0, 1, 2, 3], [1, 4, 5, 2]]

[vertices]
x = [0.0, 1.0, 1.0, 0.0, 2.0, 2.0]
y = [0.0, 0.0, 1.0, 1.0, 0.0, 1.0]

[boundary]
wall = [[0, 1], [2, 3], [1, 4], [5, 2]]
upstream = [[3, 0]]
downstream = [[4, 5]]
"#;

    #[test]
    fn two_squares_are_read() {
        let file: MeshFile = TWO_SQUARES.parse().expect("valid mesh file");
        assert_eq!(file.cell_order, CellOrder::Unspecified);
        assert_eq!(file.n_boundary_edges(), 6);
        let poly = file.to_polygons();
        assert_eq!(poly.n_cells(), 2);
        assert_eq!(poly.n_vertices(), 6);
        assert_eq!(&poly.vertices_on_cell[1], &[1, 4, 5, 2]);
        let upstream: Vec<_> = poly
            .boundary_edges
            .iter()
            .filter(|e| e.kind == EdgeKind::Upstream)
            .collect();
        assert_eq!(upstream.len(), 1);
        assert_eq!(upstream[0].vertices, [3, 0]);
    }

    #[test]
    fn out_of_range_vertex_is_reported() {
        let text = TWO_SQUARES.replace("[1, 4, 5, 2]]", "[1, 4, 6, 2]]");
        match text.parse::<MeshFile>() {
            Err(MeshFileError::Invalid(problems)) => {
                assert!(
                    problems.iter().any(|p| p.contains("cell 1")),
                    "{problems:?}"
                );
            }
            other => panic!("expected Invalid, got {other:?}"),
        }
    }

    #[test]
    fn unknown_format_and_keys_are_rejected() {
        let text = TWO_SQUARES.replace("vormorph-mesh", "other-mesh");
        assert!(matches!(
            text.parse::<MeshFile>(),
            Err(MeshFileError::Invalid(_))
        ));
        let text = TWO_SQUARES.replace("wall =", "walls =");
        assert!(matches!(
            text.parse::<MeshFile>(),
            Err(MeshFileError::Parse(_))
        ));
    }
}
