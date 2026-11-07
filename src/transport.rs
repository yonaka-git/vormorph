//! セル中心のスカラーの移流拡散。
//!
//! 水深平均の量 φ を、連続式と同じ θ 平均の単位幅流量で運ぶ。二次流の強さ I
//! （[`crate::sediment::secondary_flow`]）と、浮遊砂の濃度 C（未実装。V08）で共用する。
//! 離散式の詳細は `docs/design/overview.md` の第4.10節、判断の理由は ADR-0019。
//!
//! ## 支配方程式
//!
//! ```text
//! ∂(hφ)/∂t + ∇·(q φ) = ∇·(h D ∇φ) + h k (φ_e − φ)
//! ```
//!
//! k は緩和の速さ \[1/s\]、φ_e は緩和の目標値（生成項と消失項を合わせて書いたもの）。
//!
//! ## 離散式（セル i）
//!
//! ```text
//! A_i h_i^{n+1} φ*_i = A_i h_i^n φ_i^n − Δt Σ_j s_ij L_j q_j φ_{up(j)}^n
//!                      + Δt Σ_j (L_j K_j / δ_j) (φ_{k(i,j)}^n − φ_i^n)
//! φ_i^{n+1} = (φ*_i + Δt k_i φ_{e,i}) / (1 + Δt k_i)
//! h_i^n     = h_i^{n+1} + (Δt / A_i) Σ_j s_ij L_j q_j
//! ```
//!
//! - q_j は連続式と同じ θ 平均の単位幅流量（[`crate::hydro::HydroSolver::theta_flux`]）。
//!   up(j) は q_j の向きで決まる風上のセル（1次の風上差分）。
//! - 時刻 n の水深 h^n は、時刻 n+1 の水深と q_j から連続式を逆にたどって求める。流れの計算は
//!   V^{n+1} = V^n − Δt Σ s L q を Newton 法の残差の範囲で満たすので、h^n はその範囲で
//!   実際の値と一致し、しかも上の式と丸め誤差まで整合する。φ が一様なら φ* = φ になり、移流で
//!   一様な状態が崩れない（時刻 n の水深を別に持たなくてよい）。
//! - K_j = (hD)_j \[m³/s\] は辺の拡散の係数で、呼ぶ側が決める。拡散は陽的な2点差分で、境界辺
//!   （壁・開境界）では 0。周期境界の継ぎ目の辺は内部の辺と同じ。
//! - 境界辺から流れ込む水の φ は、呼ぶ側が辺ごとに与える（`inflow`）。流れ出る辺は内側のセルの値。
//! - 緩和は陰的に扱う（k_i = ∞ なら φ^{n+1} = φ_e）。
//! - 水深 h^{n+1} が `dry_depth` 未満のセルでは φ^{n+1} = 0 とする。
//! - φ* が φ^n の凸結合になる（新しい極値を作らない）条件は
//!   Δt Σ_j (s_ij L_j q_j の流出分 + L_j K_j / δ_j) ≤ A_i h_i^n。移流の部分は時間刻みの条件
//!   （overview の第7章）と同じで、拡散の部分は通常これよりずっと緩い。
//!
//! 計算はセルのループで、各セルは自分の値だけを書く。

use crate::mesh::{Mesh, NONE};
use crate::state::State;

/// 1ステップの移流拡散に使う値。配列の長さは、`flux`・`diffusion`・`inflow` が辺の数、
/// `rate`・`target` がセル数。
#[derive(Debug, Clone, Copy)]
pub struct TransportTerms<'a> {
    /// θ 平均の単位幅流量 q_j \[m²/s\]（向きは u_n と同じ）。
    pub flux: &'a [f64],
    /// 辺の拡散の係数 K_j = (hD)_j \[m³/s\]（境界辺の値は使わない）。
    pub diffusion: &'a [f64],
    /// 境界辺から流れ込む水の φ（内部の辺の値は使わない）。
    pub inflow: &'a [f64],
    /// 緩和の速さ k_i \[1/s\]（0 なら緩和なし、+∞ なら φ = φ_e）。
    pub rate: &'a [f64],
    /// 緩和の目標値 φ_{e,i}。
    pub target: &'a [f64],
}

/// スカラーの移流拡散の作業配列。
#[derive(Debug, Clone)]
pub struct ScalarTransport {
    /// 時刻 n の φ
    previous: Vec<f64>,
}

impl ScalarTransport {
    /// セル数 `n_cells` に合わせて確保する。
    pub fn new(n_cells: usize) -> Self {
        Self {
            previous: vec![0.0; n_cells],
        }
    }

    /// `phi`（時刻 n の φ）を `dt` だけ進めて時刻 n+1 の値にする。
    ///
    /// `state` は流れを進めた後（時刻 n+1）の状態で、水深 h^{n+1} に使う。`terms.flux` は、
    /// その流れの計算の θ 平均流量でなければならない（モジュールの説明の離散式）。
    pub fn step(
        &mut self,
        mesh: &Mesh,
        state: &State,
        terms: &TransportTerms,
        dt: f64,
        dry_depth: f64,
        phi: &mut [f64],
    ) {
        debug_assert_eq!(
            phi.len(),
            mesh.n_cells(),
            "phi must have one value per cell"
        );
        debug_assert_eq!(terms.flux.len(), mesh.n_edges());
        debug_assert_eq!(terms.diffusion.len(), mesh.n_edges());
        debug_assert_eq!(terms.inflow.len(), mesh.n_edges());
        debug_assert_eq!(terms.rate.len(), mesh.n_cells());
        debug_assert_eq!(terms.target.len(), mesh.n_cells());
        let previous = &mut self.previous;
        previous.copy_from_slice(phi);
        for (c, value) in phi.iter_mut().enumerate() {
            let depth = state.depth(c);
            if depth < dry_depth {
                *value = 0.0;
                continue;
            }
            let area = mesh.area_cell[c];
            let own = previous[c];
            // 正味の流出量、運ばれる量、拡散の量（いずれも体積 × φ / 時間）
            let (mut net, mut carried, mut spread) = (0.0, 0.0, 0.0);
            let rows = mesh.edges_on_cell[c]
                .iter()
                .zip(&mesh.edge_sign_on_cell[c])
                .zip(&mesh.cells_on_cell[c]);
            for ((&e, &s), &k) in rows {
                let out = s * mesh.dv_edge[e] * terms.flux[e];
                let upwind = if out > 0.0 {
                    own
                } else if k != NONE {
                    previous[k]
                } else {
                    terms.inflow[e]
                };
                net += out;
                carried += out * upwind;
                if k != NONE {
                    spread += mesh.dv_edge[e] * terms.diffusion[e] / mesh.dc_edge[e]
                        * (previous[k] - own);
                }
            }
            let old_depth = (depth + dt * net / area).max(0.0);
            let advected = (old_depth * own - dt * (carried - spread) / area) / depth;
            let (k, target) = (terms.rate[c], terms.target[c]);
            *value = if k.is_infinite() {
                target
            } else {
                (advected + dt * k * target) / (1.0 + dt * k)
            };
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{AdvectionScheme, SedimentSupply, TimeIntegration};
    use crate::config::{
        BoundaryConfig, ChannelLayout, DownstreamBoundary, SolverConfig, UpstreamBoundary,
    };
    use crate::hydro::{self, HydroParams, HydroSolver};
    use crate::mesh::EdgeKind;
    use crate::mesh::channel::{self, ChannelSpec};

    fn channel(length: f64, width: f64, cell_size: f64, layout: ChannelLayout) -> Mesh {
        channel::build(&ChannelSpec {
            length,
            width,
            cell_size,
            layout,
        })
        .expect("valid channel")
    }

    /// 水深 `h`、流速 `u`（x 方向）の一様流（水平な河床）。
    fn uniform_flow(mesh: &Mesh, h: f64, u: f64) -> State {
        let mut state = State::new(mesh.n_cells(), mesh.n_edges());
        state.eta.fill(h);
        for e in 0..mesh.n_edges() {
            if mesh.edge_kind[e] != EdgeKind::Wall {
                state.u_n[e] = u * mesh.normal_edge[e][0];
            }
        }
        hydro::reconstruct_cell_velocity(mesh, &state.u_n, &mut state.u_cell);
        state
    }

    /// こぶのある河床（一部は乾く）の上で流れを進めても、φ が一様なら一様のまま（丸め誤差まで）。
    /// 流れの計算の θ 平均流量と連続式から逆算した h^n が整合するため。
    #[test]
    fn uniform_scalar_stays_uniform_over_an_unsteady_flow() {
        for layout in [ChannelLayout::Rectangular, ChannelLayout::Hexagonal] {
            let mesh = channel(30.0, 6.0, 1.0, layout);
            let mut state = State::new(mesh.n_cells(), mesh.n_edges());
            for c in 0..mesh.n_cells() {
                let (x, y) = (mesh.x_cell[c], mesh.y_cell[c]);
                let bump = 0.5 * (-((x - 12.0).powi(2) + (y - 2.0).powi(2)) / 8.0).exp();
                state.zb[c] = 0.01 * (30.0 - x) + bump;
                state.eta[c] = 0.45;
            }
            let boundary = BoundaryConfig {
                upstream: UpstreamBoundary::Discharge {
                    discharge: 3.0,
                    sediment_supply: SedimentSupply::None,
                },
                downstream: DownstreamBoundary::WaterLevel { level: 0.45 },
            };
            let params = HydroParams {
                gravity: 9.81,
                theta: 0.55,
                manning_n: 0.025,
                dry_depth: 1.0e-3,
                advection: AdvectionScheme::SecondOrder,
                time_integration: TimeIntegration::PredictorCorrector,
                solver: SolverConfig::default(),
                background_slope: 0.0,
            };
            let mut hydro = HydroSolver::new(&mesh, params);
            let mut transport = ScalarTransport::new(mesh.n_cells());
            let value = 0.37;
            let mut phi: Vec<f64> = (0..mesh.n_cells())
                .map(|c| if state.depth(c) >= 1.0e-3 { value } else { 0.0 })
                .collect();
            // 境界から入る水も同じ値、拡散あり、緩和なし
            let diffusion = vec![0.05; mesh.n_edges()];
            let inflow = vec![value; mesh.n_edges()];
            let (rate, target) = (vec![0.0; mesh.n_cells()], vec![0.0; mesh.n_cells()]);
            let mut dry_cells_seen = false;
            for _ in 0..60 {
                // 乾いていたセルが濡れると、そこは φ = 0 から始まるので、濡れたままのセルだけ比べる
                let wet_before: Vec<bool> = (0..mesh.n_cells())
                    .map(|c| state.depth(c) >= 1.0e-3)
                    .collect();
                dry_cells_seen |= wet_before.iter().any(|w| !w);
                hydro
                    .step(&mesh, &mut state, &boundary, 0.05)
                    .expect("flow step");
                let terms = TransportTerms {
                    flux: hydro.theta_flux(),
                    diffusion: &diffusion,
                    inflow: &inflow,
                    rate: &rate,
                    target: &target,
                };
                transport.step(&mesh, &state, &terms, 0.05, 1.0e-3, &mut phi);
                for c in 0..mesh.n_cells() {
                    let wet_cell = wet_before[c] && state.depth(c) >= 1.0e-3;
                    let neighbors_wet = mesh.cells_on_cell[c]
                        .iter()
                        .all(|&k| k == NONE || wet_before[k]);
                    if wet_cell && neighbors_wet {
                        assert!(
                            (phi[c] - value).abs() < 1e-11,
                            "{layout:?} cell {c}: {}",
                            phi[c]
                        );
                    }
                }
                // 次のステップのために、乾いていたセルにも同じ値を入れる（一様な場から始める）
                for (c, p) in phi.iter_mut().enumerate() {
                    if state.depth(c) >= 1.0e-3 {
                        *p = value;
                    }
                }
            }
            assert!(dry_cells_seen, "{layout:?}: the bump must emerge");
        }
    }

    /// 一様流の中では、緩和は陰的な Euler 法の減衰 (1 + Δt k)^{−n} になる（移流と拡散は効かない）。
    #[test]
    fn relaxation_is_implicit_euler() {
        let mesh = channel::build_periodic(&ChannelSpec {
            length: 10.0,
            width: 4.0,
            cell_size: 1.0,
            layout: ChannelLayout::Hexagonal,
        })
        .expect("valid periodic channel");
        let state = uniform_flow(&mesh, 0.8, 1.2);
        let flux: Vec<f64> = state.u_n.iter().map(|u| 0.8 * u).collect();
        let diffusion = vec![0.1; mesh.n_edges()];
        let inflow = vec![0.0; mesh.n_edges()];
        let (k, phi_e, dt) = (0.3, 0.05, 0.2);
        let rate = vec![k; mesh.n_cells()];
        let target = vec![phi_e; mesh.n_cells()];
        let terms = TransportTerms {
            flux: &flux,
            diffusion: &diffusion,
            inflow: &inflow,
            rate: &rate,
            target: &target,
        };
        let mut transport = ScalarTransport::new(mesh.n_cells());
        let mut phi = vec![1.0; mesh.n_cells()];
        for _ in 0..10 {
            transport.step(&mesh, &state, &terms, dt, 1.0e-3, &mut phi);
        }
        let expected = phi_e + (1.0 - phi_e) * (1.0 + dt * k).powi(-10);
        for (c, p) in phi.iter().enumerate() {
            assert!((p - expected).abs() < 1e-14, "cell {c}: {p} vs {expected}");
        }
        // 緩和の速さが無限大なら、目標値そのもの
        let rate = vec![f64::INFINITY; mesh.n_cells()];
        let terms = TransportTerms {
            rate: &rate,
            ..terms
        };
        transport.step(&mesh, &state, &terms, dt, 1.0e-3, &mut phi);
        assert!(phi.iter().all(|&p| p == phi_e));
    }

    /// 上流端から φ₀ の水が流れ込む一様流の定常状態は、e^{−k x / u}（適応長 u / k）に近づき、
    /// 1次の風上差分なので、誤差は格子幅にほぼ比例して減る。
    #[test]
    fn steady_decay_from_the_inlet_converges_at_first_order() {
        let (length, h, u, k, phi0) = (40.0, 0.5, 1.0, 0.2, 1.0);
        let adaptation = u / k; // 5 m
        let errors: Vec<f64> = [1.0, 0.5]
            .iter()
            .map(|&dx| {
                let mesh = channel(length, 2.0, dx, ChannelLayout::Rectangular);
                let state = uniform_flow(&mesh, h, u);
                let flux: Vec<f64> = state.u_n.iter().map(|v| h * v).collect();
                let diffusion = vec![0.0; mesh.n_edges()];
                let inflow = vec![phi0; mesh.n_edges()];
                let rate = vec![k; mesh.n_cells()];
                let target = vec![0.0; mesh.n_cells()];
                let terms = TransportTerms {
                    flux: &flux,
                    diffusion: &diffusion,
                    inflow: &inflow,
                    rate: &rate,
                    target: &target,
                };
                let mut transport = ScalarTransport::new(mesh.n_cells());
                let mut phi = vec![0.0; mesh.n_cells()];
                let dt = 0.5 * dx / u;
                for _ in 0..(4.0 * length / (u * dt)) as usize {
                    transport.step(&mesh, &state, &terms, dt, 1.0e-3, &mut phi);
                }
                // 上流端から 10 m（適応長の2倍）までの最大誤差
                (0..mesh.n_cells())
                    .filter(|&c| mesh.x_cell[c] <= 10.0)
                    .map(|c| (phi[c] - phi0 * (-mesh.x_cell[c] / adaptation).exp()).abs())
                    .fold(0.0, f64::max)
            })
            .collect();
        let ratio = errors[1] / errors[0];
        // 誤差は 0.072、0.042（比 0.59。最初のセルの半セルのずれが残る）
        assert!(errors[0] < 0.08, "{errors:?}");
        assert!((0.5..0.65).contains(&ratio), "{errors:?}, ratio {ratio}");
    }
}
