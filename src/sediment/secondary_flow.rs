//! 二次流（らせん流）の強さ I と、掃流砂の向きの補正（ADR-0019）。
//!
//! 離散式の詳細は `docs/design/overview.md` の第2.5節と第4.10節。
//!
//! ## 二次流の強さ
//!
//! 湾曲した流れでは、遠心力で水面付近の水が外岸へ、河床付近の水が内岸へ向かう二次流が生じる。
//! その強さ I \[m/s\] を、遠心加速度を生成項とする輸送方程式で求める（`kind = "transport"`。
//! D-Flow FM の式 5.13〜5.20 と同じ形）。
//!
//! ```text
//! ∂(hI)/∂t + ∇·(q I) = ∇·(h D ∇I) + h (I_e − I) / T_I
//! I_e = h a_n / |u|,   a_n = (u × (u·∇)u)_z / |u|
//! T_I = L_I / |u|,     L_I = (1 − 2α) h / (2κ² α),   α = √g n / (κ h^{1/6})
//! D   = α_D u_* h,     u_* = √g n |u| / h^{1/6}
//! ```
//!
//! - a_n は遠心加速度（向心加速度の大きさに、左へ曲がる流れで正の符号をつけた値）で、
//!   流線の曲率を 1/R とすると a_n = |u|² / R、I_e = h |u| / R。曲率を求めて |u|³ で割るのでは
//!   なく、セル中心の流速（Perot 法）とその勾配（最小二乗。[`CellGradient`]）から直接求める。
//! - 緩和の項は h |u| (I_e − I) / L_I = (h² a_n − h |u| I) / L_I なので、|u| → 0 でも特異にならない。
//! - α ≥ 1/2（ごく浅い流れ。n = 0.025 なら h ≲ 3 mm）では L_I = 0（平衡）とする。
//! - `kind = "equilibrium"` では、各時刻の流れの平衡値 I = I_e を使う（L_I → 0 の極限）。
//!
//! ## 離散式
//!
//! 流れを進めた後（予測子・修正子では修正子の後）に1ステップに1回、共用の移流拡散
//! （[`crate::transport`]）で解く。移流は連続式と同じ θ 平均流量の1次風上差分、拡散は辺の2点差分
//! （(hD)_j = α_D u_{*j} H_j²。H_j と |u|_j は掃流力と同じ辺の水深と速さ）、緩和は陰的に扱う。
//!
//! ```text
//! I^{n+1} = (L_I I* + Δt h a_n) / (L_I + Δt |u|)    （I* は移流と拡散だけで進めた値）
//! ```
//!
//! h、|u|、a_n、L_I は時刻 n+1 の流れの値。
//!
//! - **汀線では a_n = 0**：セルか隣のセルに乾いた内部の辺（H_j = 0）があれば、流速の勾配を
//!   使わない（2次の移流と同じ規則。汀線のセルの Perot 流速は乾いた辺を流速 0 として数えるので
//!   小さく出る）。開境界のセルは 0 にしない（二次流は流れに戻らないので、2次の移流のような
//!   正のフィードバックがない）。
//! - 乾いたセル（h < `dry_depth`）では I = 0。
//! - 境界：上流端から流れ込む水の I は 0、下流端は自由流出（逆流のときは内側のセルの値）、
//!   壁と開境界の拡散の流束は 0。周期境界の継ぎ目の辺は内部の辺と同じ。
//!
//! ## 掃流砂の向き
//!
//! 河床付近の流速の向きを、水深平均流速から反時計回りに δ だけ回す（左へ曲がる流れで内岸側）。
//!
//! ```text
//! tan δ = N* I / |u|
//! q_{b,j} = q_b [ (u_j − tan δ u_t) / (|u| √(1 + tan² δ)) − γ ∂z_b/∂n ]
//! ```
//!
//! u_t = **u**·(**k** × **n**_j) は辺の接線成分（[`hydro::edge_tangential_velocity`]）、|u| は辺の
//! 速さ（[`hydro::edge_speed`]）、I は両側のセルの平均（境界辺は内側のセル）。河床付近の流速を
//! **u** + tan δ (**k** × **u**) とした向きの法線成分で、tan δ が大きくても大きさは 1 を超えない。
//! Nays2DH（Engelund モデル）は河床付近の速さを流線方向の成分で近似するので、分母の
//! √(1 + tan² δ) がない。十分発達した湾曲（I = h |u| / R）で横断方向の掃流砂が 0 になる河床の
//! 横断勾配は sin δ = γ ∂z_b/∂n（δ が小さければ N* h / R = γ ∂z_b/∂n）。
//!
//! ## 参考文献
//!
//! - Engelund (1974): 湾曲部の二次流と掃流砂の偏角（N* = 7）
//! - Deltares (2026) D-Flow FM Technical Reference Manual §5.6: 二次流の強さの輸送方程式
//! - iRIC Project (2023) Nays2DH マニュアル II.5.3〜II.5.4: 河床付近の流速と掃流砂の向き

use crate::config::{SecondaryFlowConfig, SecondaryFlowKind};
use crate::hydro;
use crate::mesh::gradient::CellGradient;
use crate::mesh::{EdgeKind, Mesh, NONE};
use crate::state::State;
use crate::transport::{ScalarTransport, TransportTerms};

/// 二次流の計算に使うパラメータ。物理定数は設定から受け取る。
#[derive(Debug, Clone, PartialEq)]
pub struct SecondaryFlowParams {
    /// 二次流の設定（求め方、N*、拡散係数の係数 α_D）。
    pub config: SecondaryFlowConfig,
    /// 重力加速度 g \[m/s²\]。
    pub gravity: f64,
    /// マニングの粗度係数 n。
    pub manning_n: f64,
    /// カルマン定数 κ。
    pub von_karman: f64,
    /// wet/dry を判定する水深の閾値 \[m\]。
    pub dry_depth: f64,
    /// 周期境界で河床高と水位から除いた縦断勾配 I \[-\]（辺の水深の計算に使う。ADR-0015）。
    pub background_slope: f64,
}

impl SecondaryFlowParams {
    /// 適応長 L_I = (1 − 2α) h / (2κ² α)、α = √g n / (κ h^{1/6}) \[m\]。
    ///
    /// α ≥ 1/2（ごく浅い流れ）では 0（平衡）。水深が 0 以下なら 0。
    pub fn adaptation_length(&self, depth: f64) -> f64 {
        if depth <= 0.0 {
            return 0.0;
        }
        let kappa = self.von_karman;
        let alpha = self.gravity.sqrt() * self.manning_n / (kappa * depth.powf(1.0 / 6.0));
        if alpha >= 0.5 {
            0.0
        } else {
            (1.0 - 2.0 * alpha) * depth / (2.0 * kappa * kappa * alpha)
        }
    }

    /// 辺の拡散の係数 (hD)_j = α_D u_* H²、u_* = √g n |u| / H^{1/6} \[m³/s\]。
    pub fn diffusion(&self, depth: f64, speed: f64) -> f64 {
        if depth <= 0.0 {
            return 0.0;
        }
        let u_star = self.gravity.sqrt() * self.manning_n * speed / depth.powf(1.0 / 6.0);
        self.config.diffusivity_coef * u_star * depth * depth
    }
}

/// 流速 `u` と、その勾配 `gradient`（∂u_a/∂x_b を `[a][b]`）から求める遠心加速度
/// a_n = (u × (u·∇)u)_z / |u| \[m/s²\]（左へ曲がる流れで正）。|u| = 0 なら 0。
pub fn centripetal_acceleration(u: [f64; 2], gradient: [[f64; 2]; 2]) -> f64 {
    let speed = u[0].hypot(u[1]);
    if speed == 0.0 {
        return 0.0;
    }
    curvature_term(u, gradient) / speed
}

/// (u × (u·∇)u)_z = u² ∂v/∂x + uv (∂v/∂y − ∂u/∂x) − v² ∂u/∂y \[m²/s³\]。
fn curvature_term(u: [f64; 2], g: [[f64; 2]; 2]) -> f64 {
    let [ux, uy] = u;
    let ax = ux * g[0][0] + uy * g[0][1];
    let ay = ux * g[1][0] + uy * g[1][1];
    ux * ay - uy * ax
}

/// 辺 `e` での掃流砂の向きの法線成分 (u_j − tan δ u_t) / (|u| √(1 + tan² δ))、
/// tan δ = `n_star` I_j / |u|（`speed` = |u| > 0）。
///
/// I_j は両側のセルの I の平均（境界辺は内側のセルの値）。I = 0 なら u_j / |u| と一致する。
pub fn bedload_direction(mesh: &Mesh, state: &State, e: usize, speed: f64, n_star: f64) -> f64 {
    let [c1, c2] = mesh.cells_on_edge[e];
    let intensity = if c2 == NONE {
        state.spiral[c1]
    } else {
        0.5 * (state.spiral[c1] + state.spiral[c2])
    };
    let tan_delta = n_star * intensity / speed;
    let tangential = hydro::edge_tangential_velocity(mesh, state, e);
    (state.u_n[e] - tan_delta * tangential) / (speed * (1.0 + tan_delta * tan_delta).sqrt())
}

/// 二次流の強さ I のソルバー。作業配列を持ち回す。
#[derive(Debug, Clone)]
pub struct SecondaryFlow {
    params: SecondaryFlowParams,
    /// セル中心の流速の勾配の重み（最小二乗）
    gradient: CellGradient,
    /// 移流拡散の作業配列
    transport: ScalarTransport,
    /// 辺の水深 H_j \[m\]
    edge_depth: Vec<f64>,
    /// 境界辺の外側の水位（すべて −∞。[`hydro::edge_depths`] に渡す）
    outside_level: Vec<f64>,
    /// セルが濡れていて、内部の辺がすべて濡れているか
    wet: Vec<bool>,
    /// 辺の拡散の係数 (hD)_j \[m³/s\]
    diffusion: Vec<f64>,
    /// 境界辺から流れ込む水の I \[m/s\]
    inflow: Vec<f64>,
    /// 緩和の速さ |u| / L_I \[1/s\]
    rate: Vec<f64>,
    /// 平衡値 I_e \[m/s\]
    target: Vec<f64>,
}

impl SecondaryFlow {
    /// 格子の大きさに合わせて作業配列を確保し、最小二乗の重みを計算する。
    pub fn new(mesh: &Mesh, params: SecondaryFlowParams) -> Self {
        let (n_cells, n_edges) = (mesh.n_cells(), mesh.n_edges());
        Self {
            params,
            gradient: CellGradient::new(mesh),
            transport: ScalarTransport::new(n_cells),
            edge_depth: vec![0.0; n_edges],
            outside_level: vec![f64::NEG_INFINITY; n_edges],
            wet: vec![false; n_cells],
            diffusion: vec![0.0; n_edges],
            inflow: vec![0.0; n_edges],
            rate: vec![0.0; n_cells],
            target: vec![0.0; n_cells],
        }
    }

    /// パラメータ。
    pub fn params(&self) -> &SecondaryFlowParams {
        &self.params
    }

    /// `state.spiral`（時刻 n の I）を `dt` だけ進める。
    ///
    /// `state` は流れを進めた後（時刻 n+1）の状態、`theta_flux` はその流れの計算の θ 平均の
    /// 単位幅流量（[`hydro::HydroSolver::theta_flux`]）。
    pub fn step(&mut self, mesh: &Mesh, state: &mut State, theta_flux: &[f64], dt: f64) {
        let Self {
            params,
            gradient,
            transport,
            edge_depth,
            outside_level,
            wet,
            diffusion,
            inflow,
            rate,
            target,
        } = self;
        let dry_depth = params.dry_depth;
        hydro::edge_depths(
            mesh,
            state,
            outside_level,
            dry_depth,
            params.background_slope,
            edge_depth,
        );
        for (c, w) in wet.iter_mut().enumerate() {
            *w = state.depth(c) >= dry_depth
                && mesh.edges_on_cell[c]
                    .iter()
                    .zip(&mesh.cells_on_cell[c])
                    .all(|(&e, &k)| k == NONE || edge_depth[e] > 0.0);
        }

        // セル：平衡値 I_e = h a_n / |u| と、緩和の速さ |u| / L_I
        for (c, (k, i_e)) in rate.iter_mut().zip(target.iter_mut()).enumerate() {
            let depth = state.depth(c);
            let u = state.u_cell[c];
            let speed = u[0].hypot(u[1]);
            let reliable = wet[c] && mesh.cells_on_cell[c].iter().all(|&n| n == NONE || wet[n]);
            *i_e = if reliable && speed > 0.0 {
                let g = gradient.of_vector(mesh, &state.u_cell, c);
                depth * curvature_term(u, g) / (speed * speed)
            } else {
                0.0
            };
            let length = params.adaptation_length(depth);
            *k = if length > 0.0 {
                speed / length
            } else if speed > 0.0 {
                f64::INFINITY
            } else {
                0.0
            };
        }

        if params.config.kind == SecondaryFlowKind::Equilibrium {
            for (c, &i_e) in target.iter().enumerate() {
                state.spiral[c] = if state.depth(c) >= dry_depth {
                    i_e
                } else {
                    0.0
                };
            }
            return;
        }

        // 辺：拡散の係数と、境界から流れ込む水の I
        for (e, (d, i_in)) in diffusion.iter_mut().zip(inflow.iter_mut()).enumerate() {
            let interior = mesh.cells_on_edge[e][1] != NONE;
            *d = if interior && edge_depth[e] > 0.0 {
                params.diffusion(edge_depth[e], hydro::edge_speed(mesh, state, e))
            } else {
                0.0
            };
            *i_in = match mesh.edge_kind[e] {
                EdgeKind::Downstream => state.spiral[mesh.cells_on_edge[e][0]],
                EdgeKind::Upstream | EdgeKind::Wall | EdgeKind::Interior => 0.0,
            };
        }
        let terms = TransportTerms {
            flux: theta_flux,
            diffusion,
            inflow,
            rate,
            target,
        };
        let mut spiral = std::mem::take(&mut state.spiral);
        transport.step(mesh, state, &terms, dt, dry_depth, &mut spiral);
        state.spiral = spiral;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::ChannelLayout;
    use crate::mesh::channel::{self, ChannelSpec};

    fn params(kind: SecondaryFlowKind) -> SecondaryFlowParams {
        SecondaryFlowParams {
            config: SecondaryFlowConfig {
                kind,
                n_star: 7.0,
                diffusivity_coef: 0.15,
            },
            gravity: 9.81,
            manning_n: 0.025,
            von_karman: 0.41,
            dry_depth: 1.0e-3,
            background_slope: 0.0,
        }
    }

    fn spec(layout: ChannelLayout) -> ChannelSpec {
        ChannelSpec {
            length: 20.0,
            width: 8.0,
            cell_size: 1.0,
            layout,
        }
    }

    /// 境界辺を持つセルからの距離（隣のセルをたどった数。境界辺を持つセルは 0）。
    fn rings_from_boundary(mesh: &Mesh) -> Vec<usize> {
        let mut ring: Vec<usize> = (0..mesh.n_cells())
            .map(|c| {
                if mesh.edges_on_cell[c]
                    .iter()
                    .any(|&e| mesh.is_boundary_edge(e))
                {
                    0
                } else {
                    usize::MAX
                }
            })
            .collect();
        for level in 0..mesh.n_cells() {
            let mut changed = false;
            for c in 0..mesh.n_cells() {
                if ring[c] == usize::MAX
                    && mesh.cells_on_cell[c]
                        .iter()
                        .any(|&k| k != NONE && ring[k] == level)
                {
                    ring[c] = level + 1;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        ring
    }

    /// 点 (x0, y0) のまわりの角速度 ω の剛体回転（反時計回りが正）。水深 `h` は一様。
    /// 壁の辺にも同じ流速を与える（Perot 法を全セルで厳密にするため）。
    fn solid_body_rotation(mesh: &Mesh, center: [f64; 2], omega: f64, h: f64) -> State {
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        state.eta.fill(h);
        for e in 0..mesh.n_edges() {
            let [x, y] = [mesh.x_edge[e], mesh.y_edge[e]];
            let u = [-omega * (y - center[1]), omega * (x - center[0])];
            let [nx, ny] = mesh.normal_edge[e];
            state.u_n[e] = u[0] * nx + u[1] * ny;
        }
        hydro::reconstruct_cell_velocity(mesh, &state.u_n, &mut state.u_cell);
        state
    }

    /// 剛体回転の流れでは a_n = |u|² / r、I_e = h |u| / r = h ω（一様）。正多角形のセルでは
    /// Perot 法と最小二乗の勾配が1次式の流速場で厳密なので、丸め誤差まで一致する。反時計回り
    /// （左へ曲がる流れ）で正。汀線（乾いたセルの2つ隣まで）では a_n = 0。
    #[test]
    fn equilibrium_intensity_of_solid_body_rotation() {
        let (h, omega) = (0.5, 0.05);
        for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
            let mesh = channel::build(&spec(layout)).expect("valid channel");
            let rings = rings_from_boundary(&mesh);
            for (center, sign) in [([10.0, -30.0], 1.0), ([10.0, 40.0], -1.0)] {
                let mut state = solid_body_rotation(&mesh, center, sign * omega, h);
                let mut solver = SecondaryFlow::new(&mesh, params(SecondaryFlowKind::Equilibrium));
                let flux = vec![0.0; mesh.n_edges()];
                solver.step(&mesh, &mut state, &flux, 0.1);
                for (c, &ring) in rings.iter().enumerate() {
                    // 正方形配置は全セルが正方形。六角形配置は、境界の切れたセルの Perot 流速が
                    // 厳密でないので、その隣のセルまで除く
                    if layout == ChannelLayout::Hexagonal && ring < 2 {
                        continue;
                    }
                    let i = state.spiral[c];
                    assert!(
                        (i - sign * h * omega).abs() < 1e-12,
                        "{layout:?} cell {c}: {i}"
                    );
                    let a = centripetal_acceleration(
                        state.u_cell[c],
                        solver.gradient.of_vector(&mesh, &state.u_cell, c),
                    );
                    let r = (mesh.x_cell[c] - center[0]).hypot(mesh.y_cell[c] - center[1]);
                    assert!((a - sign * omega * omega * r).abs() < 1e-12, "{a}");
                }
            }

            // 1つのセルを乾かすと、そのセルと、隣（乾いた辺を持つ）、その隣（勾配に乾いた辺の
            // セルの流速を使う）は a_n = 0
            let mut state = solid_body_rotation(&mesh, [10.0, -30.0], omega, h);
            let dry = (0..mesh.n_cells())
                .find(|&c| rings[c] >= 3)
                .expect("an interior cell");
            state.zb[dry] = 1.0;
            state.eta[dry] = 1.0;
            let mut solver = SecondaryFlow::new(&mesh, params(SecondaryFlowKind::Equilibrium));
            solver.step(&mesh, &mut state, &vec![0.0; mesh.n_edges()], 0.1);
            let neighbors = || mesh.cells_on_cell[dry].iter().copied();
            let near: Vec<usize> = neighbors()
                .flat_map(|k| mesh.cells_on_cell[k].iter().copied())
                .chain(neighbors())
                .filter(|&k| k != NONE)
                .collect();
            for &c in &near {
                assert_eq!(state.spiral[c], 0.0, "{layout:?} cell {c}");
            }
            let far = (0..mesh.n_cells())
                .find(|c| rings[*c] >= 2 && !near.contains(c))
                .expect("a far cell");
            assert!(
                (state.spiral[far] - h * omega).abs() < 1e-12,
                "{layout:?} cell {far} (dry {dry}): {}",
                state.spiral[far]
            );
        }
    }

    /// 適応長：L_I = (1 − 2α) h / (2κ² α)（手計算の値）。ごく浅い（α ≥ 1/2）なら 0。
    #[test]
    fn adaptation_length_matches_the_formula() {
        let p = params(SecondaryFlowKind::Transport);
        // h = 1 m：α = √9.81 × 0.025 / 0.41 = 0.190981、L_I = 9.6256 m
        assert!((p.adaptation_length(1.0) - 9.625571).abs() < 1e-6);
        // h = 0.1 m：α = 0.190981 × 10^{1/6} = 0.280336、L_I = (1 − 2α) 0.1 / (2κ² α)
        let alpha = 0.190981 * 10.0_f64.powf(1.0 / 6.0);
        let expected = (1.0 - 2.0 * alpha) * 0.1 / (2.0 * 0.41 * 0.41 * alpha);
        assert!((p.adaptation_length(0.1) - expected).abs() < 1e-5 * expected);
        // α = 1/2 となる水深 (2√g n / κ)^6 = 3.1 mm より浅いと 0
        let limit = (2.0 * 9.81_f64.sqrt() * 0.025 / 0.41).powi(6);
        assert!((limit - 3.105e-3).abs() < 1e-6, "{limit}");
        assert!(p.adaptation_length(1.01 * limit) > 0.0);
        assert_eq!(p.adaptation_length(0.99 * limit), 0.0);
        assert_eq!(p.adaptation_length(0.0), 0.0);
    }

    /// 直線水路の一様流（周期境界）では I_e = 0 で、一様な I は (1 + Δt |u| / L_I)^{−n} で
    /// 減衰する（陰的な緩和。移流と拡散は一様な場を変えない）。平衡型なら直ちに 0。
    #[test]
    fn intensity_relaxes_in_a_straight_channel() {
        let (h, u, dt, steps) = (1.0, 1.5, 0.5, 20);
        for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
            let mesh = channel::build_periodic(&spec(layout)).expect("valid periodic channel");
            let mut state = State::new(mesh.n_cells(), mesh.n_edges());
            state.eta.fill(h);
            for e in 0..mesh.n_edges() {
                if mesh.edge_kind[e] != EdgeKind::Wall {
                    state.u_n[e] = u * mesh.normal_edge[e][0];
                }
            }
            hydro::reconstruct_cell_velocity(&mesh, &state.u_n, &mut state.u_cell);
            let flux: Vec<f64> = state.u_n.iter().map(|v| h * v).collect();
            let i0 = 0.02;
            state.spiral.fill(i0);
            let mut solver = SecondaryFlow::new(&mesh, params(SecondaryFlowKind::Transport));
            for _ in 0..steps {
                solver.step(&mesh, &mut state, &flux, dt);
            }
            let length = solver.params().adaptation_length(h);
            let expected = i0 * (1.0 + dt * u / length).powi(-steps);
            for (c, &i) in state.spiral.iter().enumerate() {
                assert!(
                    (i - expected).abs() < 1e-14,
                    "{layout:?} cell {c}: {i} vs {expected}"
                );
            }
            // 平衡型なら直ちに I_e = 0（Perot 流速の丸め誤差の分だけ残る）
            let mut equilibrium = SecondaryFlow::new(&mesh, params(SecondaryFlowKind::Equilibrium));
            equilibrium.step(&mesh, &mut state, &flux, dt);
            assert!(state.spiral.iter().all(|&i| i.abs() < 1e-15));
        }
    }

    /// 掃流砂の向き：I > 0（左へ曲がる流れ）なら、流速の向きから反時計回りに δ 回る。
    /// 法線成分は cos δ n_x + sin δ n_y（x 方向の流れ）、I = 0 なら u_j / |u| と一致する。
    #[test]
    fn bedload_direction_turns_towards_the_inner_bank() {
        let u = 1.2;
        for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
            let mesh = channel::build_periodic(&spec(layout)).expect("valid periodic channel");
            let mut state = State::new(mesh.n_cells(), mesh.n_edges());
            state.eta.fill(1.0);
            for e in 0..mesh.n_edges() {
                if mesh.edge_kind[e] != EdgeKind::Wall {
                    state.u_n[e] = u * mesh.normal_edge[e][0];
                }
            }
            hydro::reconstruct_cell_velocity(&mesh, &state.u_n, &mut state.u_cell);
            for e in (0..mesh.n_edges()).filter(|&e| mesh.edge_kind[e] != EdgeKind::Wall) {
                let speed = hydro::edge_speed(&mesh, &state, e);
                let plain = bedload_direction(&mesh, &state, e, speed, 7.0);
                assert_eq!(plain, state.u_n[e] / speed);
            }
            state.spiral.fill(0.03);
            let delta = (7.0 * 0.03 / u).atan();
            for e in (0..mesh.n_edges()).filter(|&e| mesh.edge_kind[e] != EdgeKind::Wall) {
                let speed = hydro::edge_speed(&mesh, &state, e);
                let turned = bedload_direction(&mesh, &state, e, speed, 7.0);
                let [nx, ny] = mesh.normal_edge[e];
                let expected = delta.cos() * nx + delta.sin() * ny;
                assert!((turned - expected).abs() < 1e-14, "{layout:?} edge {e}");
            }
        }
    }
}
