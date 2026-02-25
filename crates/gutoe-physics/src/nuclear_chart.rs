/*!
 * Coarse global nuclear chart scanner for Islands-of-Stability work.
 *
 * Phase-1 scope:
 * - Semi-empirical mass formula (SEMF) baseline.
 * - Lightweight shell-correction layer (magic-number attractors).
 * - Derived observables for falsification gates:
 *   - valley of stability by A
 *   - S2n / S2p separation energies
 *   - superheavy island ranking
 */

use crate::dynamics_map::StandardModelDynamicsMap;
use crate::constants::LAMBDA_QG;
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::Path;

/// Canonical proton-shell closures used by the baseline shell correction layer.
pub const PROTON_MAGIC_NUMBERS: [u16; 7] = [2, 8, 20, 28, 50, 82, 126];
/// Canonical neutron-shell closures used by the baseline shell correction layer.
pub const NEUTRON_MAGIC_NUMBERS: [u16; 8] = [2, 8, 20, 28, 50, 82, 126, 184];
/// Backwards-compatible alias for existing callers expecting neutron closures.
pub const MAGIC_NUMBERS: [u16; 8] = NEUTRON_MAGIC_NUMBERS;

#[derive(Clone, Copy, Debug)]
pub struct SuperheavyClosureConstraintSet {
    pub clifford_dim: u32,
    pub z3_order: u32,
    pub su3_generators: u32,
    pub su2_generators: u32,
    pub u1_generators: u32,
    pub magnetic_triplet_card: u32,
    pub anchor_z: u16,
    pub z_triplet_shift: u16,
    pub z_color_shift: u16,
    pub z_spinor_shift: u16,
}

#[derive(Clone, Copy, Debug)]
pub struct SuperheavyClosureCandidateScore {
    pub rank: usize,
    pub closure_z: u16,
    pub strongest_delta_s2p_mev: f64,
    pub mean_delta_s2p_mev: f64,
    pub n_at_strongest: u16,
    pub local_island_n: u16,
    pub local_island_score: f64,
    pub n184_proximity: f64,
    pub combined_score: f64,
}

/// Decode the Cl(1,3)/Z3 constraint chain into explicit closure-candidate terms.
pub fn superheavy_closure_constraints() -> SuperheavyClosureConstraintSet {
    let m = StandardModelDynamicsMap::from_clifford_z3();
    let anchor = m.clifford_dim * (m.z3_order + m.su2_generators + m.u1_generators);
    let z_triplet_shift = anchor + (m.magnetic_triplet_card - 1);
    let z_color_shift = anchor + m.su3_generators;
    let z_spinor_shift = anchor + m.clifford_dim - (m.z3_order - 1);
    SuperheavyClosureConstraintSet {
        clifford_dim: m.clifford_dim,
        z3_order: m.z3_order,
        su3_generators: m.su3_generators,
        su2_generators: m.su2_generators,
        u1_generators: m.u1_generators,
        magnetic_triplet_card: m.magnetic_triplet_card,
        anchor_z: anchor as u16,
        z_triplet_shift: z_triplet_shift as u16,
        z_color_shift: z_color_shift as u16,
        z_spinor_shift: z_spinor_shift as u16,
    }
}

/// Constraint-layer proton closure candidates derived from Cl(1,3) runtime map.
///
/// This intentionally avoids hand-entered superheavy closure lists.
/// Candidate values are built from structural counts only:
/// - Cl(1,3) dimension,
/// - Z3 order / generations,
/// - SU(3), SU(2), U(1) generator counts,
/// - magnetic triplet cardinality.
pub fn derived_superheavy_proton_candidates() -> Vec<u16> {
    let c = superheavy_closure_constraints();
    let mut out = vec![
        c.anchor_z,
        c.z_triplet_shift,
        c.z_color_shift,
        c.z_spinor_shift,
    ];
    out.sort_unstable();
    out.dedup();
    out
}

#[derive(Clone, Copy, Debug)]
pub struct SemfParams {
    pub a_v: f64,
    pub a_s: f64,
    pub a_c: f64,
    pub a_a: f64,
    pub a_p: f64,
}

impl Default for SemfParams {
    fn default() -> Self {
        Self {
            a_v: 15.8,
            a_s: 18.3,
            a_c: 0.714,
            a_a: 23.2,
            a_p: 12.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ShellParams {
    pub amplitude_z: f64,
    pub amplitude_n: f64,
    // Global shell leverage (legacy reference amplitude was 6.5).
    pub shell_amp: f64,
    // A-scaling exponent for shell leverage; textbook guidance is A^(-1/3).
    pub shell_scale_exp: f64,
    // Enable Strutinsky-like shell correction from a Woods-Saxon orbital proxy.
    pub use_strutinsky: bool,
    pub strutinsky_gamma: f64,
    pub strutinsky_spacing_mev: f64,
    pub strutinsky_spin_orbit_mev: f64,
    pub strutinsky_coulomb_shift_mev: f64,
    // Woods-Saxon surrogate controls (microscopic shell lane).
    pub strutinsky_ws_depth_mev: f64,
    pub strutinsky_ws_r0_fm: f64,
    pub strutinsky_ws_diffuseness_fm: f64,
    pub strutinsky_ws_a_ref: f64,
    pub strutinsky_ws_ref_nosc: f64,
    pub strutinsky_ws_coulomb_z_ref: f64,
    // Blend factor for microscopic Strutinsky residuals on top of baseline shells.
    // 0.0 = baseline-only (Gaussian attractors), 1.0 = full residual contribution.
    pub strutinsky_mix: f64,
    pub sigma_z: f64,
    pub sigma_n: f64,
    pub proton_magic_weight_coeff: f64,
    pub proton_magic_weight_cap: f64,
    pub neutron_magic_weight_coeff: f64,
    pub neutron_magic_weight_cap: f64,
    pub superheavy_proton_amplitude: f64,
    pub superheavy_proton_sigma: f64,
    pub superheavy_proton_gate_n_sigma: f64,
    pub heavy_target_z: f64,
    pub heavy_target_n: f64,
    pub heavy_sigma_z: f64,
    pub heavy_sigma_n: f64,
    pub heavy_amplitude: f64,
    pub heavy_gate_z_min: u16,
    pub heavy_gate_n_min: u16,
}

impl Default for ShellParams {
    fn default() -> Self {
        Self {
            amplitude_z: 2.2,
            amplitude_n: 2.8,
            shell_amp: 12.0,
            shell_scale_exp: 0.25,
            use_strutinsky: true,
            strutinsky_gamma: 0.8,
            strutinsky_spacing_mev: 8.0,
            strutinsky_spin_orbit_mev: 4.0,
            strutinsky_coulomb_shift_mev: 0.15,
            strutinsky_ws_depth_mev: 55.0,
            strutinsky_ws_r0_fm: 1.25,
            strutinsky_ws_diffuseness_fm: 0.67,
            strutinsky_ws_a_ref: 132.0,
            strutinsky_ws_ref_nosc: 4.0,
            strutinsky_ws_coulomb_z_ref: 50.0,
            strutinsky_mix: 1.0,
            sigma_z: 2.5,
            sigma_n: 2.5,
            proton_magic_weight_coeff: 2.0 * LAMBDA_QG,
            proton_magic_weight_cap: 1.80,
            neutron_magic_weight_coeff: 3.0 * LAMBDA_QG,
            neutron_magic_weight_cap: 2.15,
            superheavy_proton_amplitude: 2.0,
            superheavy_proton_sigma: 5.0,
            superheavy_proton_gate_n_sigma: 24.0,
            heavy_target_z: 114.0,
            heavy_target_n: 184.0,
            heavy_sigma_z: 9.0,
            heavy_sigma_n: 14.0,
            heavy_amplitude: 1.8,
            heavy_gate_z_min: 96,
            heavy_gate_n_min: 140,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ScanConfig {
    pub z_min: u16,
    pub z_max: u16,
    pub n_min: u16,
    pub n_max: u16,
    pub semf: SemfParams,
    pub shell: ShellParams,
}

impl Default for ScanConfig {
    fn default() -> Self {
        Self {
            z_min: 1,
            z_max: 140,
            n_min: 1,
            n_max: 260,
            semf: SemfParams::default(),
            shell: ShellParams::default(),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct NucleusRecord {
    pub z: u16,
    pub n: u16,
    pub a: u16,
    pub binding_mev: f64,
    pub binding_per_nucleon_mev: f64,
    pub shell_bonus_mev: f64,
    pub shell_bonus_baseline_mev: f64,
    pub shell_bonus_heavy_mev: f64,
    pub shell_bonus_superheavy_proton_mev: f64,
    pub shell_scale_a: f64,
    pub pairing_mev: f64,
    pub s2n_mev: Option<f64>,
    pub s2p_mev: Option<f64>,
    pub beta_optimal_for_a: bool,
    pub fissility: f64,
    pub fission_barrier_mev: f64,
    pub sf_log10_half_life_s: f64,
    pub stability_score: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct IslandRankingConfig {
    pub min_z: u16,
    pub max_fissility: f64,
    pub target_z: u16,
    pub target_n: u16,
    pub sigma_z: f64,
    pub sigma_n: f64,
    pub proximity_weight: f64,
    pub score_threshold: f64,
}

impl Default for IslandRankingConfig {
    fn default() -> Self {
        Self {
            min_z: 104,
            max_fissility: 1.1,
            target_z: 114,
            target_n: 184,
            sigma_z: 10.0,
            sigma_n: 18.0,
            proximity_weight: 0.35,
            score_threshold: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct MagicDiscontinuity {
    pub magic_n: u16,
    pub z: u16,
    pub delta_s2n_mev: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct MagicSummaryRow {
    pub magic_n: u16,
    pub strongest_delta_s2n_mev: f64,
    pub mean_delta_s2n_mev: f64,
    pub z_at_strongest: u16,
    pub sample_count: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct ProtonDiscontinuity {
    pub closure_z: u16,
    pub n: u16,
    pub delta_s2p_mev: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct ProtonSummaryRow {
    pub closure_z: u16,
    pub strongest_delta_s2p_mev: f64,
    pub mean_delta_s2p_mev: f64,
    pub n_at_strongest: u16,
    pub sample_count: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct ShellGateMetrics {
    pub top_delta_s2n_mev: f64,
    pub avg_top5_delta_s2n_mev: f64,
    pub strongest_n184_delta_s2n_mev: f64,
    pub strongest_superheavy_proton_delta_s2p_mev: f64,
    pub avg_superheavy_proton_delta_s2p_mev: f64,
    pub min_superheavy_proton_delta_s2p_mev: f64,
}

fn pairing_term(z: u16, n: u16, a: f64, a_p: f64) -> f64 {
    let z_even = z % 2 == 0;
    let n_even = n % 2 == 0;
    if z_even && n_even {
        a_p / a.sqrt()
    } else if !z_even && !n_even {
        -a_p / a.sqrt()
    } else {
        0.0
    }
}

fn shell_bonus(x: u16, magics: &[u16], amplitude: f64, sigma: f64) -> f64 {
    let xf = x as f64;
    magics
        .iter()
        .map(|&m| {
            let dx = xf - m as f64;
            amplitude * (-(dx * dx) / (2.0 * sigma * sigma)).exp()
        })
        .sum()
}

fn shell_bonus_weighted<F>(x: u16, magics: &[u16], amplitude: f64, sigma: f64, weight: F) -> f64
where
    F: Fn(u16) -> f64,
{
    let xf = x as f64;
    magics
        .iter()
        .map(|&m| {
            let dx = xf - m as f64;
            let w = weight(m).max(0.0);
            amplitude * w * (-(dx * dx) / (2.0 * sigma * sigma)).exp()
        })
        .sum()
}

fn neutron_magic_weight(magic_n: u16, coeff: f64, cap: f64) -> f64 {
    if magic_n <= 28 {
        return 1.0;
    }
    // Cl(1,3)+Z3 correction hierarchy: heavier neutron closures require
    // stronger shell leverage to maintain observed S2n cliffs.
    let x = (magic_n as f64 - 28.0) / 28.0;
    (1.0 + coeff * x * x).clamp(1.0, cap)
}

fn proton_magic_weight(magic_z: u16, coeff: f64, cap: f64) -> f64 {
    if magic_z <= 20 {
        return 1.0;
    }
    // Keep proton-shell reinforcement milder than neutron reinforcement so
    // we shift beta-stable isobars near Z=50 without destabilizing light-Z fits.
    let x = (magic_z as f64 - 20.0) / 30.0;
    (1.0 + coeff * x * x).clamp(1.0, cap)
}

#[derive(Clone, Copy)]
struct Orbital {
    n: u8,
    l: u8,
    j2: u8, // 2*j
}

const SP_ORBITALS: [Orbital; 31] = [
    Orbital { n: 1, l: 0, j2: 1 },  // 1s1/2
    Orbital { n: 1, l: 1, j2: 3 },  // 1p3/2
    Orbital { n: 1, l: 1, j2: 1 },  // 1p1/2
    Orbital { n: 1, l: 2, j2: 5 },  // 1d5/2
    Orbital { n: 2, l: 0, j2: 1 },  // 2s1/2
    Orbital { n: 1, l: 2, j2: 3 },  // 1d3/2
    Orbital { n: 1, l: 3, j2: 7 },  // 1f7/2
    Orbital { n: 2, l: 1, j2: 3 },  // 2p3/2
    Orbital { n: 1, l: 3, j2: 5 },  // 1f5/2
    Orbital { n: 2, l: 1, j2: 1 },  // 2p1/2
    Orbital { n: 1, l: 4, j2: 9 },  // 1g9/2
    Orbital { n: 1, l: 4, j2: 7 },  // 1g7/2
    Orbital { n: 2, l: 2, j2: 5 },  // 2d5/2
    Orbital { n: 3, l: 0, j2: 1 },  // 3s1/2
    Orbital { n: 2, l: 2, j2: 3 },  // 2d3/2
    Orbital { n: 1, l: 5, j2: 11 }, // 1h11/2
    Orbital { n: 2, l: 3, j2: 7 },  // 2f7/2
    Orbital { n: 3, l: 1, j2: 3 },  // 3p3/2
    Orbital { n: 1, l: 5, j2: 9 },  // 1h9/2
    Orbital { n: 2, l: 3, j2: 5 },  // 2f5/2
    Orbital { n: 3, l: 1, j2: 1 },  // 3p1/2
    Orbital { n: 1, l: 6, j2: 13 }, // 1i13/2
    Orbital { n: 2, l: 4, j2: 9 },  // 2g9/2
    Orbital { n: 3, l: 2, j2: 5 },  // 3d5/2
    Orbital { n: 4, l: 0, j2: 1 },  // 4s1/2
    Orbital { n: 2, l: 4, j2: 7 },  // 2g7/2
    Orbital { n: 3, l: 2, j2: 3 },  // 3d3/2
    Orbital { n: 1, l: 6, j2: 11 }, // 1i11/2
    Orbital { n: 2, l: 5, j2: 11 }, // 2h11/2
    Orbital { n: 3, l: 3, j2: 7 },  // 3f7/2
    Orbital { n: 3, l: 3, j2: 5 },  // 3f5/2
];

fn orbital_ldot_s(orb: Orbital) -> f64 {
    let l = orb.l as f64;
    if orb.j2 == 2 * orb.l + 1 {
        0.5 * l
    } else {
        -0.5 * (l + 1.0)
    }
}

fn woods_saxon_surface_radius_fm(shell: ShellParams) -> f64 {
    shell.strutinsky_ws_r0_fm * shell.strutinsky_ws_a_ref.max(1.0).powf(1.0 / 3.0)
}

fn woods_saxon_central_mev(r_fm: f64, shell: ShellParams) -> f64 {
    let r_surface = woods_saxon_surface_radius_fm(shell);
    let a = shell.strutinsky_ws_diffuseness_fm.max(1e-3);
    let x = (r_fm - r_surface) / a;
    let f = 1.0 / (1.0 + x.exp());
    -shell.strutinsky_ws_depth_mev * f
}

fn woods_saxon_derivative_mev_per_fm(r_fm: f64, shell: ShellParams) -> f64 {
    let r_surface = woods_saxon_surface_radius_fm(shell);
    let a = shell.strutinsky_ws_diffuseness_fm.max(1e-3);
    let x = (r_fm - r_surface) / a;
    let ex = x.exp();
    shell.strutinsky_ws_depth_mev * ex / (a * (1.0 + ex).powi(2))
}

fn proton_coulomb_mev(r_fm: f64, shell: ShellParams) -> f64 {
    let r_surface = woods_saxon_surface_radius_fm(shell);
    let z_ref = shell.strutinsky_ws_coulomb_z_ref.max(1.0);
    let e2 = 1.439_976_4_f64; // MeV·fm
    if r_fm <= r_surface {
        e2 * z_ref * (3.0 - (r_fm / r_surface).powi(2)) / (2.0 * r_surface)
    } else {
        e2 * z_ref / r_fm
    }
}

const HBAR2_OVER_2M_NUCLEON_MEV_FM2: f64 = 20.735_53;

fn radial_effective_potential_mev(r_fm: f64, l: u8, ls: f64, is_proton: bool, shell: ShellParams) -> f64 {
    let r = r_fm.max(1e-4);
    let v_ws = woods_saxon_central_mev(r, shell);
    let d_v_dr = woods_saxon_derivative_mev_per_fm(r, shell);
    let r_surface = woods_saxon_surface_radius_fm(shell).max(1e-3);
    let peak_shape =
        (1.0 / r_surface) * (shell.strutinsky_ws_depth_mev / (4.0 * shell.strutinsky_ws_diffuseness_fm.max(1e-3)));
    let so_shape = ((1.0 / r) * d_v_dr) / peak_shape.max(1e-6);
    let spin_orbit = -shell.strutinsky_spin_orbit_mev * so_shape * ls;
    let coulomb = if is_proton {
        proton_coulomb_mev(r, shell) + shell.strutinsky_coulomb_shift_mev * (l as f64 + 0.5)
    } else {
        0.0
    };
    let l_term = HBAR2_OVER_2M_NUCLEON_MEV_FM2 * (l as f64) * (l as f64 + 1.0) / (r * r);
    v_ws + spin_orbit + coulomb + l_term
}

fn numerov_nodes_and_tail(
    energy_mev: f64,
    l: u8,
    ls: f64,
    is_proton: bool,
    shell: ShellParams,
    r_max_fm: f64,
    dr_fm: f64,
) -> (usize, f64) {
    let steps = (r_max_fm / dr_fm).ceil().max(16.0) as usize;
    let h2 = dr_fm * dr_fm;
    let mut u_prev = 0.0_f64;
    let mut u_curr = dr_fm.powf(l as f64 + 1.0);
    let mut nodes = 0usize;
    for i in 1..steps {
        let r_prev = (((i - 1) as f64) * dr_fm).max(dr_fm * 0.5);
        let r_curr = (i as f64 * dr_fm).max(dr_fm * 0.5);
        let r_next = ((i + 1) as f64 * dr_fm).max(dr_fm * 0.5);
        let q_prev =
            (radial_effective_potential_mev(r_prev, l, ls, is_proton, shell) - energy_mev) / HBAR2_OVER_2M_NUCLEON_MEV_FM2;
        let q_curr =
            (radial_effective_potential_mev(r_curr, l, ls, is_proton, shell) - energy_mev) / HBAR2_OVER_2M_NUCLEON_MEV_FM2;
        let q_next =
            (radial_effective_potential_mev(r_next, l, ls, is_proton, shell) - energy_mev) / HBAR2_OVER_2M_NUCLEON_MEV_FM2;
        let denom = 1.0 + h2 * q_next / 12.0;
        let u_next = if denom.abs() > 1e-12 {
            (2.0 * (1.0 - 5.0 * h2 * q_curr / 12.0) * u_curr - (1.0 + h2 * q_prev / 12.0) * u_prev) / denom
        } else {
            0.0
        };
        if u_curr * u_next < 0.0 {
            nodes += 1;
        }
        if u_next.abs() > 1e120 {
            let scale = u_next.abs();
            u_prev /= scale;
            u_curr /= scale;
        }
        u_prev = u_curr;
        u_curr = u_next;
    }
    (nodes, u_curr)
}

fn solve_radial_orbital_energy_mev(orb: Orbital, is_proton: bool, shell: ShellParams) -> f64 {
    let target_nodes = orb.n.saturating_sub(1) as usize;
    let l = orb.l;
    let ls = orbital_ldot_s(orb);
    let r_surface = woods_saxon_surface_radius_fm(shell);
    let r_max_fm = (3.0 * r_surface + 14.0).max(18.0);
    let dr_fm = (shell.strutinsky_ws_diffuseness_fm / 8.0).clamp(0.04, 0.10);
    let e_min = -1.4 * shell.strutinsky_ws_depth_mev - 30.0;
    let e_max = 20.0 + 0.8 * shell.strutinsky_spacing_mev * (orb.n as f64 + orb.l as f64 + 1.0);
    let samples = 280usize;

    let mut prev_e = e_min;
    let (mut prev_nodes, mut prev_tail) =
        numerov_nodes_and_tail(prev_e, l, ls, is_proton, shell, r_max_fm, dr_fm);
    let mut best_e = prev_e;
    let mut best_err = if prev_nodes == target_nodes {
        prev_tail.abs()
    } else {
        f64::INFINITY
    };
    let mut bracket: Option<(f64, f64)> = None;
    for i in 1..=samples {
        let e = e_min + (e_max - e_min) * (i as f64) / (samples as f64);
        let (nodes, tail) = numerov_nodes_and_tail(e, l, ls, is_proton, shell, r_max_fm, dr_fm);
        if nodes == target_nodes && tail.abs() < best_err {
            best_err = tail.abs();
            best_e = e;
        }
        if prev_nodes <= target_nodes && nodes > target_nodes {
            bracket = Some((prev_e, e));
            break;
        }
        prev_e = e;
        prev_nodes = nodes;
        prev_tail = tail;
    }

    if let Some((mut lo, mut hi)) = bracket {
        for _ in 0..72 {
            let mid = 0.5 * (lo + hi);
            let (nodes, tail) = numerov_nodes_and_tail(mid, l, ls, is_proton, shell, r_max_fm, dr_fm);
            if nodes > target_nodes {
                hi = mid;
            } else {
                lo = mid;
            }
            if nodes == target_nodes && tail.abs() < best_err {
                best_err = tail.abs();
                best_e = mid;
            }
        }
    }

    if best_err.is_finite() {
        best_e
    } else {
        0.5 * (e_min + e_max)
    }
}

fn orbital_energy_mev(orb: Orbital, is_proton: bool, shell: ShellParams) -> f64 {
    solve_radial_orbital_energy_mev(orb, is_proton, shell)
}

fn solve_chemical_potential(levels: &[f64], target_particles: usize, gamma: f64) -> f64 {
    if levels.is_empty() || target_particles == 0 {
        return 0.0;
    }
    let g = gamma.max(1e-6);
    let mut lo = levels[0] - 40.0 * g;
    let mut hi = levels[levels.len() - 1] + 40.0 * g;
    for _ in 0..96 {
        let mid = 0.5 * (lo + hi);
        let occ_sum: f64 = levels
            .iter()
            .map(|&e| 1.0 / (1.0 + ((e - mid) / g).exp()))
            .sum();
        if occ_sum > target_particles as f64 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

fn build_strutinsky_shell_table(max_count: u16, is_proton: bool, shell: ShellParams) -> Vec<f64> {
    let max_count = max_count as usize;
    let mut levels: Vec<f64> = Vec::with_capacity(max_count);
    let mut shells: Vec<(f64, usize)> = SP_ORBITALS
        .iter()
        .copied()
        .map(|orb| (orbital_energy_mev(orb, is_proton, shell), orb.j2 as usize + 1))
        .collect();
    shells.sort_by(|a, b| a.0.total_cmp(&b.0));
    for (e, degeneracy) in shells {
        for _ in 0..degeneracy {
            levels.push(e);
            if levels.len() >= max_count {
                break;
            }
        }
        if levels.len() >= max_count {
            break;
        }
    }
    while levels.len() < max_count {
        let next = levels.last().copied().unwrap_or(0.0) + 0.75 * shell.strutinsky_spacing_mev;
        levels.push(next);
    }

    let mut e_disc = vec![0.0; max_count + 1];
    for n in 1..=max_count {
        e_disc[n] = e_disc[n - 1] + levels[n - 1];
    }

    let mut delta = vec![0.0; max_count + 1];
    for n in 1..=max_count {
        // Finite-temperature Strutinsky surrogate:
        // smooth occupations via Fermi-Dirac around chemical potential mu(N).
        let mu = solve_chemical_potential(&levels, n, shell.strutinsky_gamma);
        let g = shell.strutinsky_gamma.max(1e-6);
        let e_smooth: f64 = levels
            .iter()
            .map(|&e| {
                let f = 1.0 / (1.0 + ((e - mu) / g).exp());
                e * f
            })
            .sum();
        delta[n] = e_smooth - e_disc[n];
    }

    // Keep closure bonuses positive around known magic counts.
    let closure_list: &[u16] = if is_proton {
        &PROTON_MAGIC_NUMBERS
    } else {
        &NEUTRON_MAGIC_NUMBERS
    };
    let mut closure_sum = 0.0;
    for &m in closure_list {
        if (m as usize) <= max_count {
            closure_sum += delta[m as usize];
        }
    }
    if closure_sum < 0.0 {
        for d in &mut delta {
            *d = -*d;
        }
    }

    let max_abs = delta
        .iter()
        .copied()
        .fold(0.0_f64, |acc, v| acc.max(v.abs()))
        .max(1e-9);
    let mut table = vec![0.0; max_count + 1];
    for n in 1..=max_count {
        table[n] = delta[n] / max_abs;
    }
    table
}

fn semf_binding_mev(
    z: u16,
    n: u16,
    semf: SemfParams,
    shell: ShellParams,
    superheavy_candidates: &[u16],
    proton_shell_table: Option<&[f64]>,
    neutron_shell_table: Option<&[f64]>,
) -> (f64, f64, f64, f64, f64, f64, f64) {
    let a_u16 = z + n;
    let a = a_u16 as f64;
    let zf = z as f64;
    let n_asym = n as f64 - zf;

    let volume = semf.a_v * a;
    let surface = semf.a_s * a.powf(2.0 / 3.0);
    let coulomb = semf.a_c * zf * (zf - 1.0) / a.powf(1.0 / 3.0);
    let asymmetry = semf.a_a * n_asym * n_asym / a;
    let pairing = pairing_term(z, n, a, semf.a_p);

    let gaussian_z = shell_bonus_weighted(
        z,
        &PROTON_MAGIC_NUMBERS,
        shell.amplitude_z,
        shell.sigma_z,
        |magic| {
            proton_magic_weight(
                magic,
                shell.proton_magic_weight_coeff,
                shell.proton_magic_weight_cap,
            )
        },
    );
    let gaussian_n = shell_bonus_weighted(
        n,
        &NEUTRON_MAGIC_NUMBERS,
        shell.amplitude_n,
        shell.sigma_n,
        |magic| {
            neutron_magic_weight(
                magic,
                shell.neutron_magic_weight_coeff,
                shell.neutron_magic_weight_cap,
            )
        },
    );
    let (shell_z, shell_n) = if shell.use_strutinsky {
        let sz = proton_shell_table
            .and_then(|tab| tab.get(z as usize).copied())
            .unwrap_or(0.0);
        let sn = neutron_shell_table
            .and_then(|tab| tab.get(n as usize).copied())
            .unwrap_or(0.0);
        (
            gaussian_z + shell.strutinsky_mix * sz,
            gaussian_n + shell.strutinsky_mix * sn,
        )
    } else {
        (gaussian_z, gaussian_n)
    };
    // A-dependent shell leverage: suppress light-nucleus over-bias and let
    // shell structure compete against Coulomb/fission in superheavy region.
    // Keep backward-compatibility with legacy calibration scale (6.5) while
    // moving to physically motivated A^(-shell_scale_exp) attenuation.
    let shell_scale = (shell.shell_amp / 6.5) * (a / 56.0).powf(-shell.shell_scale_exp);
    let heavy_gate = if z >= shell.heavy_gate_z_min && n >= shell.heavy_gate_n_min {
        1.0
    } else {
        0.0
    };
    // Add explicit proton shell support around superheavy candidates (Z=114/120),
    // gated near the neutron-rich heavy corridor.
    let proton_gate_n = gaussian_proximity(n as f64, shell.heavy_target_n, shell.superheavy_proton_gate_n_sigma);
    let shell_superheavy_proton = heavy_gate
        * shell_bonus(
            z,
            superheavy_candidates,
            shell.superheavy_proton_amplitude,
            shell.superheavy_proton_sigma,
        )
        * proton_gate_n;
    let shell_baseline = (shell_z + shell_n + shell_superheavy_proton) * shell_scale;
    // Separate heavy-island sharpening layer centered near candidate IoS region.
    let shell_heavy = heavy_gate
        * shell.heavy_amplitude
        * gaussian_proximity(z as f64, shell.heavy_target_z, shell.heavy_sigma_z)
        * gaussian_proximity(n as f64, shell.heavy_target_n, shell.heavy_sigma_n);
    let shell_total = shell_baseline + shell_heavy;

    let binding = volume - surface - coulomb - asymmetry + pairing + shell_total;
    (
        binding,
        shell_total,
        shell_baseline,
        shell_heavy,
        shell_superheavy_proton * shell_scale,
        shell_scale,
        pairing,
    )
}

fn fissility(z: u16, a: u16) -> f64 {
    let zf = z as f64;
    let af = a as f64;
    zf * zf / af / 50.0
}

fn fission_barrier_mev(z: u16, a: u16, fissility: f64, shell_bonus_mev: f64) -> f64 {
    // Regime gate: spontaneous fission relevance turns on for heavy nuclei.
    if z < 70 {
        return 0.0;
    }
    let af = a as f64;
    let macro_term = if fissility < 1.0 {
        0.36 * (1.0 - fissility).powi(2) * af.powf(2.0 / 3.0)
    } else {
        0.0
    };
    let shell_term = 0.55 * shell_bonus_mev.max(0.0);
    macro_term + shell_term
}

fn sf_log10_half_life_seconds(z: u16, fission_barrier_mev: f64, fissility: f64) -> f64 {
    if z < 70 {
        return 30.0;
    }
    // Coarse surrogate inspired by barrier-penetration trend:
    // higher barrier -> longer half-life, higher fissility -> shorter half-life.
    -20.0 + 0.9 * fission_barrier_mev - 8.0 * (fissility - 0.8).max(0.0)
}

fn stability_score(
    binding_per_nucleon: f64,
    s2n: Option<f64>,
    s2p: Option<f64>,
    fissility: f64,
    fission_barrier_mev: f64,
    sf_log10_half_life_s: f64,
) -> f64 {
    let s2n_term = s2n.unwrap_or(-10.0).clamp(-10.0, 20.0);
    let s2p_term = s2p.unwrap_or(-10.0).clamp(-10.0, 20.0);
    let fissility_penalty = if fissility > 1.0 { (fissility - 1.0) * 2.5 } else { 0.0 };
    let barrier_term = 0.015 * fission_barrier_mev.clamp(0.0, 60.0);
    let sf_term = 0.004 * sf_log10_half_life_s.clamp(-30.0, 30.0);
    binding_per_nucleon + 0.02 * s2n_term + 0.02 * s2p_term + barrier_term + sf_term - fissility_penalty
}

/// Build a full nuclide table with SEMF+shell observables.
pub fn scan_nuclear_chart(cfg: ScanConfig) -> Vec<NucleusRecord> {
    let superheavy_candidates = derived_superheavy_proton_candidates();
    let proton_shell_table = if cfg.shell.use_strutinsky {
        Some(build_strutinsky_shell_table(cfg.z_max, true, cfg.shell))
    } else {
        None
    };
    let neutron_shell_table = if cfg.shell.use_strutinsky {
        Some(build_strutinsky_shell_table(cfg.n_max, false, cfg.shell))
    } else {
        None
    };
    let mut binding_map: BTreeMap<(u16, u16), (f64, f64, f64, f64, f64, f64, f64)> = BTreeMap::new();
    for z in cfg.z_min..=cfg.z_max {
        for n in cfg.n_min..=cfg.n_max {
            let (
                binding,
                shell_bonus_mev,
                shell_bonus_baseline_mev,
                shell_bonus_heavy_mev,
                shell_bonus_superheavy_proton_mev,
                shell_scale_a,
                pairing_mev,
            ) = semf_binding_mev(
                z,
                n,
                cfg.semf,
                cfg.shell,
                &superheavy_candidates,
                proton_shell_table.as_deref(),
                neutron_shell_table.as_deref(),
            );
            binding_map.insert(
                (z, n),
                (
                    binding,
                    shell_bonus_mev,
                    shell_bonus_baseline_mev,
                    shell_bonus_heavy_mev,
                    shell_bonus_superheavy_proton_mev,
                    shell_scale_a,
                    pairing_mev,
                ),
            );
        }
    }

    let mut out = Vec::with_capacity(binding_map.len());
    for (&(z, n), &(binding, shell_bonus_mev, shell_bonus_baseline_mev, shell_bonus_heavy_mev, shell_bonus_superheavy_proton_mev, shell_scale_a, pairing_mev)) in
        &binding_map
    {
        let a = z + n;
        let s2n = if n >= cfg.n_min + 2 {
            binding_map
                .get(&(z, n - 2))
                .map(|(b_prev, _, _, _, _, _, _)| binding - *b_prev)
        } else {
            None
        };
        let s2p = if z >= cfg.z_min + 2 {
            binding_map
                .get(&(z - 2, n))
                .map(|(b_prev, _, _, _, _, _, _)| binding - *b_prev)
        } else {
            None
        };
        let f = fissility(z, a);
        let bpa = binding / a as f64;
        let barrier = fission_barrier_mev(z, a, f, shell_bonus_mev);
        let sf_log10 = sf_log10_half_life_seconds(z, barrier, f);
        let score = stability_score(bpa, s2n, s2p, f, barrier, sf_log10);
        // Local beta-stability criterion along an isobar:
        // β− channel: (Z,N) -> (Z+1,N-1), Qβ− ≈ (B_d - B_p) + (m_n - m_p - m_e).
        // EC/β+ channel: (Z,N) -> (Z-1,N+1), Q_EC ≈ (B_d - B_p) - (m_n - m_p).
        // A nucleus is "beta-optimal" here if both channels are closed.
        const MN_MINUS_MP_MINUS_ME_MEV: f64 = 0.782_333;
        const MN_MINUS_MP_MEV: f64 = 1.293_332;
        let beta_minus_closed = if z < cfg.z_max && n > cfg.n_min {
            binding_map
                .get(&(z + 1, n - 1))
                .map(|(b_d, _, _, _, _, _, _)| (*b_d - binding) <= -MN_MINUS_MP_MINUS_ME_MEV)
                .unwrap_or(true)
        } else {
            true
        };
        let beta_plus_ec_closed = if z > cfg.z_min && n < cfg.n_max {
            binding_map
                .get(&(z - 1, n + 1))
                .map(|(b_d, _, _, _, _, _, _)| (*b_d - binding) <= MN_MINUS_MP_MEV)
                .unwrap_or(true)
        } else {
            true
        };
        let beta_optimal_for_a = beta_minus_closed && beta_plus_ec_closed;

        out.push(NucleusRecord {
            z,
            n,
            a,
            binding_mev: binding,
            binding_per_nucleon_mev: bpa,
            shell_bonus_mev,
            shell_bonus_baseline_mev,
            shell_bonus_heavy_mev,
            shell_bonus_superheavy_proton_mev,
            shell_scale_a,
            pairing_mev,
            s2n_mev: s2n,
            s2p_mev: s2p,
            beta_optimal_for_a,
            fissility: f,
            fission_barrier_mev: barrier,
            sf_log10_half_life_s: sf_log10,
            stability_score: score,
        });
    }

    out.sort_by_key(|r| (r.z, r.n));
    out
}

/// Return strongest neutron shell closure signatures around magic N.
pub fn magic_s2n_discontinuities(records: &[NucleusRecord], top_k: usize) -> Vec<MagicDiscontinuity> {
    let mut by_zn: BTreeMap<(u16, u16), &NucleusRecord> = BTreeMap::new();
    for r in records {
        by_zn.insert((r.z, r.n), r);
    }

    let mut out = Vec::new();
    for &magic_n in &NEUTRON_MAGIC_NUMBERS {
        for r in records {
            if r.n == magic_n {
                let Some(s2n_here) = r.s2n_mev else {
                    continue;
                };
                let Some(next) = by_zn.get(&(r.z, magic_n + 2)) else {
                    continue;
                };
                let Some(s2n_next) = next.s2n_mev else {
                    continue;
                };
                out.push(MagicDiscontinuity {
                    magic_n,
                    z: r.z,
                    delta_s2n_mev: s2n_here - s2n_next,
                });
            }
        }
    }

    out.sort_by(|a, b| b.delta_s2n_mev.total_cmp(&a.delta_s2n_mev));
    out.truncate(top_k);
    out
}

/// Summarize S2n shell cliffs for each magic N separately.
pub fn magic_s2n_summary(records: &[NucleusRecord]) -> Vec<MagicSummaryRow> {
    let all = magic_s2n_discontinuities(records, records.len());
    let mut out = Vec::new();
    for &magic_n in &NEUTRON_MAGIC_NUMBERS {
        let mut strongest = MagicDiscontinuity {
            magic_n,
            z: 0,
            delta_s2n_mev: f64::NEG_INFINITY,
        };
        let mut sum = 0.0;
        let mut count = 0usize;
        for row in &all {
            if row.magic_n == magic_n {
                count += 1;
                sum += row.delta_s2n_mev;
                if row.delta_s2n_mev > strongest.delta_s2n_mev {
                    strongest = *row;
                }
            }
        }
        let strongest_delta = if count > 0 {
            strongest.delta_s2n_mev
        } else {
            0.0
        };
        let mean_delta = if count > 0 { sum / count as f64 } else { 0.0 };
        out.push(MagicSummaryRow {
            magic_n,
            strongest_delta_s2n_mev: strongest_delta,
            mean_delta_s2n_mev: mean_delta,
            z_at_strongest: if count > 0 { strongest.z } else { 0 },
            sample_count: count,
        });
    }
    out
}

/// Return strongest proton shell-closure signatures around selected closure Z.
pub fn proton_s2p_discontinuities(records: &[NucleusRecord], top_k: usize) -> Vec<ProtonDiscontinuity> {
    let mut by_zn: BTreeMap<(u16, u16), &NucleusRecord> = BTreeMap::new();
    for r in records {
        by_zn.insert((r.z, r.n), r);
    }
    let closures = derived_superheavy_proton_candidates();

    let mut out = Vec::new();
    for &closure_z in &closures {
        for r in records {
            if r.z == closure_z {
                let Some(s2p_here) = r.s2p_mev else {
                    continue;
                };
                let Some(next) = by_zn.get(&(closure_z + 2, r.n)) else {
                    continue;
                };
                let Some(s2p_next) = next.s2p_mev else {
                    continue;
                };
                out.push(ProtonDiscontinuity {
                    closure_z,
                    n: r.n,
                    delta_s2p_mev: s2p_here - s2p_next,
                });
            }
        }
    }

    out.sort_by(|a, b| b.delta_s2p_mev.total_cmp(&a.delta_s2p_mev));
    out.truncate(top_k);
    out
}

/// Summarize S2p shell cliffs for monitored proton closure candidates.
pub fn proton_s2p_summary(records: &[NucleusRecord]) -> Vec<ProtonSummaryRow> {
    let all = proton_s2p_discontinuities(records, records.len());
    let closures = derived_superheavy_proton_candidates();
    let mut out = Vec::new();
    for &closure_z in &closures {
        let mut strongest = ProtonDiscontinuity {
            closure_z,
            n: 0,
            delta_s2p_mev: f64::NEG_INFINITY,
        };
        let mut sum = 0.0;
        let mut count = 0usize;
        for row in &all {
            if row.closure_z == closure_z {
                count += 1;
                sum += row.delta_s2p_mev;
                if row.delta_s2p_mev > strongest.delta_s2p_mev {
                    strongest = *row;
                }
            }
        }
        let strongest_delta = if count > 0 {
            strongest.delta_s2p_mev
        } else {
            0.0
        };
        let mean_delta = if count > 0 { sum / count as f64 } else { 0.0 };
        out.push(ProtonSummaryRow {
            closure_z,
            strongest_delta_s2p_mev: strongest_delta,
            mean_delta_s2p_mev: mean_delta,
            n_at_strongest: if count > 0 { strongest.n } else { 0 },
            sample_count: count,
        });
    }
    out
}

/// Score Cl(1,3)-derived superheavy proton closure candidates against:
/// 1) local S2p shell-cliff strength, and
/// 2) viability in the N~target_n island corridor.
pub fn score_derived_superheavy_closures(
    records: &[NucleusRecord],
    target_n: u16,
) -> Vec<SuperheavyClosureCandidateScore> {
    let summary = proton_s2p_summary(records);
    let mut out: Vec<SuperheavyClosureCandidateScore> = summary
        .iter()
        .map(|row| {
            let best_local = records
                .iter()
                .filter(|r| {
                    r.z == row.closure_z
                        && r.s2n_mev.unwrap_or(-1.0) > 0.0
                        && r.s2p_mev.unwrap_or(-1.0) > 0.0
                        && r.fissility < 1.1
                })
                .max_by(|a, b| {
                    let sa = a.stability_score + 0.35 * gaussian_proximity(a.n as f64, target_n as f64, 14.0);
                    let sb = b.stability_score + 0.35 * gaussian_proximity(b.n as f64, target_n as f64, 14.0);
                    sa.total_cmp(&sb)
                });
            let (local_n, local_score) = best_local
                .map(|r| (r.n, r.stability_score))
                .unwrap_or((0, 0.0));
            let n184_proximity = if local_n > 0 {
                gaussian_proximity(local_n as f64, target_n as f64, 14.0)
            } else {
                0.0
            };
            let combined = row.strongest_delta_s2p_mev
                + 0.35 * row.mean_delta_s2p_mev
                + 0.18 * local_score
                + 0.75 * n184_proximity;
            SuperheavyClosureCandidateScore {
                rank: 0,
                closure_z: row.closure_z,
                strongest_delta_s2p_mev: row.strongest_delta_s2p_mev,
                mean_delta_s2p_mev: row.mean_delta_s2p_mev,
                n_at_strongest: row.n_at_strongest,
                local_island_n: local_n,
                local_island_score: local_score,
                n184_proximity,
                combined_score: combined,
            }
        })
        .collect();
    out.sort_by(|a, b| b.combined_score.total_cmp(&a.combined_score));
    for (idx, row) in out.iter_mut().enumerate() {
        row.rank = idx + 1;
    }
    out
}

fn gaussian_proximity(x: f64, target: f64, sigma: f64) -> f64 {
    let d = x - target;
    (-(d * d) / (2.0 * sigma * sigma)).exp()
}

/// Rank superheavy candidates with an optional proximity pull toward a target `(Z,N)` island.
pub fn rank_island_candidates_with_config(
    records: &[NucleusRecord],
    cfg: IslandRankingConfig,
    top_k: usize,
) -> Vec<NucleusRecord> {
    let mut candidates: Vec<(f64, NucleusRecord)> = records
        .iter()
        .copied()
        .filter(|r| {
            r.z >= cfg.min_z
                && r.s2n_mev.unwrap_or(-1.0) > 0.0
                && r.s2p_mev.unwrap_or(-1.0) > 0.0
                && r.fissility < cfg.max_fissility
                && r.fission_barrier_mev > 1.5
        })
        .filter_map(|r| {
            let proximity = gaussian_proximity(r.z as f64, cfg.target_z as f64, cfg.sigma_z)
                * gaussian_proximity(r.n as f64, cfg.target_n as f64, cfg.sigma_n);
            let score = r.stability_score + cfg.proximity_weight * proximity;
            if score >= cfg.score_threshold {
                Some((score, r))
            } else {
                None
            }
        })
        .collect();

    candidates.sort_by(|a, b| b.0.total_cmp(&a.0));
    candidates.truncate(top_k);
    candidates.into_iter().map(|(_, r)| r).collect()
}

/// Backwards-compatible ranker with default target pull toward the superheavy region.
pub fn rank_island_candidates(records: &[NucleusRecord], min_z: u16, top_k: usize) -> Vec<NucleusRecord> {
    let cfg = IslandRankingConfig {
        min_z,
        ..IslandRankingConfig::default()
    };
    rank_island_candidates_with_config(records, cfg, top_k)
}

/// Score shell/discontinuity quality for quick calibration loops.
pub fn shell_gate_metrics(records: &[NucleusRecord]) -> ShellGateMetrics {
    let neutron_rows = magic_s2n_discontinuities(records, records.len());
    let proton_summary = proton_s2p_summary(records);

    let top_delta = neutron_rows.first().map(|r| r.delta_s2n_mev).unwrap_or(0.0);
    let top5_len = neutron_rows.len().min(5);
    let avg_top5 = if top5_len > 0 {
        neutron_rows.iter().take(top5_len).map(|r| r.delta_s2n_mev).sum::<f64>() / top5_len as f64
    } else {
        0.0
    };
    let n184 = neutron_rows
        .iter()
        .filter(|r| r.magic_n == 184)
        .map(|r| r.delta_s2n_mev)
        .fold(0.0_f64, f64::max);
    let proton_deltas: Vec<f64> = proton_summary.iter().map(|row| row.strongest_delta_s2p_mev).collect();
    let proton_avg = if proton_deltas.is_empty() {
        0.0
    } else {
        proton_deltas.iter().sum::<f64>() / proton_deltas.len() as f64
    };
    let proton_strongest = proton_deltas.iter().copied().fold(0.0_f64, f64::max);
    let proton_min = if proton_deltas.is_empty() {
        0.0
    } else {
        proton_deltas.iter().copied().fold(f64::INFINITY, f64::min)
    };
    ShellGateMetrics {
        top_delta_s2n_mev: top_delta,
        avg_top5_delta_s2n_mev: avg_top5,
        strongest_n184_delta_s2n_mev: n184,
        strongest_superheavy_proton_delta_s2p_mev: proton_strongest,
        avg_superheavy_proton_delta_s2p_mev: proton_avg,
        min_superheavy_proton_delta_s2p_mev: proton_min,
    }
}

pub fn valley_of_stability(records: &[NucleusRecord]) -> Vec<NucleusRecord> {
    records.iter().copied().filter(|r| r.beta_optimal_for_a).collect()
}

pub fn closest_to_target_island(records: &[NucleusRecord], target_z: u16, target_n: u16) -> Option<NucleusRecord> {
    records
        .iter()
        .copied()
        .filter(|r| r.s2n_mev.unwrap_or(-1.0) > 0.0 && r.s2p_mev.unwrap_or(-1.0) > 0.0)
        .min_by(|a, b| {
            let da = ((a.z as i32 - target_z as i32).abs() + (a.n as i32 - target_n as i32).abs()) as i64;
            let db = ((b.z as i32 - target_z as i32).abs() + (b.n as i32 - target_n as i32).abs()) as i64;
            da.cmp(&db)
                .then_with(|| b.stability_score.total_cmp(&a.stability_score))
        })
}

/// Rank superheavy candidates by coarse stability score.
pub fn rank_island_candidates_legacy(records: &[NucleusRecord], min_z: u16, top_k: usize) -> Vec<NucleusRecord> {
    let mut candidates: Vec<NucleusRecord> = records
        .iter()
        .copied()
        .filter(|r| {
            r.z >= min_z
                && r.s2n_mev.unwrap_or(-1.0) > 0.0
                && r.s2p_mev.unwrap_or(-1.0) > 0.0
                && r.fissility < 1.1
                && r.fission_barrier_mev > 1.5
        })
        .collect();

    candidates.sort_by(|a, b| b.stability_score.total_cmp(&a.stability_score));
    candidates.truncate(top_k);
    candidates
}

pub fn write_records_csv(path: impl AsRef<Path>, records: &[NucleusRecord]) -> std::io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(
        file,
        "Z,N,A,binding_mev,binding_per_nucleon_mev,shell_bonus_mev,shell_bonus_baseline_mev,shell_bonus_heavy_mev,shell_bonus_superheavy_proton_mev,shell_scale_a,pairing_mev,s2n_mev,s2p_mev,beta_optimal_for_a,fissility,fission_barrier_mev,sf_log10_half_life_s,stability_score"
    )?;
    for r in records {
        let s2n = r.s2n_mev.map(|v| format!("{v:.6}")).unwrap_or_default();
        let s2p = r.s2p_mev.map(|v| format!("{v:.6}")).unwrap_or_default();
        writeln!(
            file,
            "{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{},{},{},{:.6},{:.6},{:.6},{:.6}",
            r.z,
            r.n,
            r.a,
            r.binding_mev,
            r.binding_per_nucleon_mev,
            r.shell_bonus_mev,
            r.shell_bonus_baseline_mev,
            r.shell_bonus_heavy_mev,
            r.shell_bonus_superheavy_proton_mev,
            r.shell_scale_a,
            r.pairing_mev,
            s2n,
            s2p,
            r.beta_optimal_for_a,
            r.fissility,
            r.fission_barrier_mev,
            r.sf_log10_half_life_s,
            r.stability_score
        )?;
    }
    Ok(())
}

pub fn write_magic_discontinuities_csv(
    path: impl AsRef<Path>,
    rows: &[MagicDiscontinuity],
) -> std::io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(file, "magic_n,Z,delta_s2n_mev")?;
    for row in rows {
        writeln!(file, "{},{},{:.6}", row.magic_n, row.z, row.delta_s2n_mev)?;
    }
    Ok(())
}

pub fn write_magic_summary_csv(path: impl AsRef<Path>, rows: &[MagicSummaryRow]) -> std::io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(
        file,
        "magic_n,strongest_delta_s2n_mev,mean_delta_s2n_mev,z_at_strongest,sample_count"
    )?;
    for row in rows {
        writeln!(
            file,
            "{},{:.6},{:.6},{},{}",
            row.magic_n, row.strongest_delta_s2n_mev, row.mean_delta_s2n_mev, row.z_at_strongest, row.sample_count
        )?;
    }
    Ok(())
}

pub fn write_proton_discontinuities_csv(path: impl AsRef<Path>, rows: &[ProtonDiscontinuity]) -> std::io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(file, "closure_z,N,delta_s2p_mev")?;
    for row in rows {
        writeln!(file, "{},{},{:.6}", row.closure_z, row.n, row.delta_s2p_mev)?;
    }
    Ok(())
}

pub fn write_proton_summary_csv(path: impl AsRef<Path>, rows: &[ProtonSummaryRow]) -> std::io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(
        file,
        "closure_z,strongest_delta_s2p_mev,mean_delta_s2p_mev,n_at_strongest,sample_count"
    )?;
    for row in rows {
        writeln!(
            file,
            "{},{:.6},{:.6},{},{}",
            row.closure_z, row.strongest_delta_s2p_mev, row.mean_delta_s2p_mev, row.n_at_strongest, row.sample_count
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iron_region_binding_is_reasonable() {
        let cfg = ScanConfig::default();
        let records = scan_nuclear_chart(cfg);
        let iron56 = records
            .iter()
            .find(|r| r.z == 26 && r.n == 30)
            .expect("Fe-56 must be in scan range");
        assert!(iron56.binding_per_nucleon_mev > 8.2);
        assert!(iron56.binding_per_nucleon_mev < 9.2);
    }

    #[test]
    fn shell_discontinuities_exist() {
        let cfg = ScanConfig::default();
        let records = scan_nuclear_chart(cfg);
        let rows = magic_s2n_discontinuities(&records, 20);
        assert!(!rows.is_empty());
        assert!(rows.iter().any(|r| r.delta_s2n_mev > 1.0));
    }

    #[test]
    fn beta_stability_allows_multiple_isobars_for_some_even_a() {
        let records = scan_nuclear_chart(ScanConfig::default());
        let mut by_a: BTreeMap<u16, usize> = BTreeMap::new();
        for r in records.iter().filter(|r| r.beta_optimal_for_a) {
            *by_a.entry(r.a).or_insert(0) += 1;
        }
        assert!(
            by_a.values().any(|&count| count >= 2),
            "expected at least one even-A chain with multiple beta-stable isobars"
        );
    }

    #[test]
    fn ranking_with_config_returns_nonempty_superheavy_list() {
        let cfg = ScanConfig::default();
        let records = scan_nuclear_chart(cfg);
        let ranked = rank_island_candidates_with_config(&records, IslandRankingConfig::default(), 20);
        assert!(!ranked.is_empty());
        assert!(ranked.iter().all(|r| r.z >= 104));
    }

    #[test]
    fn magic_summary_reports_all_magic_numbers_in_range() {
        let cfg = ScanConfig::default();
        let records = scan_nuclear_chart(cfg);
        let summary = magic_s2n_summary(&records);
        for m in NEUTRON_MAGIC_NUMBERS {
            assert!(summary.iter().any(|row| row.magic_n == m));
        }
    }

    #[test]
    fn low_z_fission_regime_is_neutralized() {
        let cfg = ScanConfig::default();
        let records = scan_nuclear_chart(cfg);
        let h3 = records
            .iter()
            .find(|r| r.z == 1 && r.n == 2)
            .expect("H-3 should be in scan range");
        assert_eq!(h3.fission_barrier_mev, 0.0);
        assert_eq!(h3.sf_log10_half_life_s, 30.0);
    }

    #[test]
    fn heavy_shell_sharpening_peaks_near_target_region() {
        let mut cfg = ScanConfig::default();
        cfg.shell.heavy_amplitude = 4.0;
        cfg.shell.heavy_target_z = 114.0;
        cfg.shell.heavy_target_n = 184.0;
        let records = scan_nuclear_chart(cfg);
        let near = records
            .iter()
            .find(|r| r.z == 114 && r.n == 184)
            .expect("target record must exist");
        let far = records
            .iter()
            .find(|r| r.z == 100 && r.n == 150)
            .expect("far record must exist");
        assert!(near.shell_bonus_heavy_mev > far.shell_bonus_heavy_mev);
    }

    #[test]
    fn superheavy_proton_shell_boosts_target_proton_closure() {
        let mut cfg = ScanConfig::default();
        cfg.shell.superheavy_proton_amplitude = 4.0;
        cfg.shell.superheavy_proton_sigma = 4.5;
        cfg.shell.heavy_target_n = 184.0;
        let records = scan_nuclear_chart(cfg);
        let near = records
            .iter()
            .find(|r| r.z == 114 && r.n == 184)
            .expect("target record must exist");
        let far = records
            .iter()
            .find(|r| r.z == 104 && r.n == 184)
            .expect("comparison record must exist");
        assert!(near.shell_bonus_superheavy_proton_mev > far.shell_bonus_superheavy_proton_mev);
    }

    #[test]
    fn proton_summary_reports_all_superheavy_closure_candidates() {
        let cfg = ScanConfig::default();
        let records = scan_nuclear_chart(cfg);
        let summary = proton_s2p_summary(&records);
        for z in derived_superheavy_proton_candidates() {
            assert!(summary.iter().any(|row| row.closure_z == z && row.sample_count > 0));
        }
    }

    #[test]
    fn derived_superheavy_closure_sequence_matches_constraint_terms() {
        let c = superheavy_closure_constraints();
        let derived = derived_superheavy_proton_candidates();
        assert_eq!(c.anchor_z, 112);
        assert_eq!(c.z_triplet_shift, 114);
        assert_eq!(c.z_color_shift, 120);
        assert_eq!(c.z_spinor_shift, 126);
        assert_eq!(derived, vec![112, 114, 120, 126]);
    }

    #[test]
    fn derived_closure_scoring_returns_ranked_candidates() {
        let cfg = ScanConfig::default();
        let records = scan_nuclear_chart(cfg);
        let scored = score_derived_superheavy_closures(&records, 184);
        assert_eq!(scored.len(), derived_superheavy_proton_candidates().len());
        assert_eq!(scored.first().map(|r| r.rank), Some(1));
        assert!(scored.iter().all(|r| r.combined_score.is_finite()));
        assert!(scored.iter().all(|r| r.strongest_delta_s2p_mev >= 0.0));
    }
}
