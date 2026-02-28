//! Structural nuclear lane for GRAND-106/107/108.
//!
//! This module exposes a deterministic, no-env-override coefficient chain:
//! Cl(1,3) counts -> NN potential proxy + shell controls + scan config.

use crate::constants::{LAMBDA_QG, VISIBLE_STATE_COUNT_STRUCTURAL};
use crate::dynamics_map::StandardModelDynamicsMap;
use crate::nuclear_chart::{
    derived_superheavy_proton_candidates, ScanConfig, SemfParams, ShellParams,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NnPotentialParams {
    pub attractive_depth_mev: f64,
    pub repulsive_core_mev: f64,
    pub range_fm: f64,
    pub core_radius_fm: f64,
    pub spin_orbit_mev: f64,
}

#[derive(Clone, Debug)]
pub struct StructuralNuclearModel {
    pub map: StandardModelDynamicsMap,
    pub semf: SemfParams,
    pub shell: ShellParams,
    pub nn: NnPotentialParams,
}

/// Derive a no-fit nuclear model from shared Cl(1,3) primitives.
pub fn derive_structural_nuclear_model() -> StructuralNuclearModel {
    let map = StandardModelDynamicsMap::from_clifford_z3();

    let clifford = map.clifford_dim as f64;
    let su3 = map.su3_generators as f64;
    let su2 = map.su2_generators as f64;
    let u1 = map.u1_generators as f64;
    let gauge_total = map.total_gauge_generators as f64;

    // Ticket GRAND-106: NN potential proxy from structural depth/scale terms.
    // Ticket GRAND-107/108: same coefficients feed shell + mass-lane scan.
    let semf = SemfParams {
        a_v: clifford - 2.0 * LAMBDA_QG,
        a_s: clifford + su2 - su3 * LAMBDA_QG,
        a_c: su3 * LAMBDA_QG,
        a_a: clifford + gauge_total / 2.0 + su3 / 8.0,
        a_p: gauge_total,
    };

    let candidates = derived_superheavy_proton_candidates();
    let heavy_target_z = candidates
        .iter()
        .copied()
        .find(|&z| z == 114)
        .unwrap_or(114) as f64;

    let shell = ShellParams {
        amplitude_z: su2 - LAMBDA_QG,
        amplitude_n: su3 / 2.0 - LAMBDA_QG,
        shell_amp: clifford - su2 + gauge_total / 2.0,
        shell_scale_exp: 1.0 / (su2 + u1),
        use_strutinsky: true,
        strutinsky_gamma: 2.0 / (su2 + u1),
        strutinsky_spacing_mev: su3,
        strutinsky_spin_orbit_mev: su2 + u1,
        strutinsky_coulomb_shift_mev: 2.0 * LAMBDA_QG,
        strutinsky_ws_depth_mev: 3.0 * clifford + gauge_total / 2.0,
        strutinsky_ws_r0_fm: (gauge_total + su2) / gauge_total,
        strutinsky_ws_diffuseness_fm: 2.0 / 3.0,
        strutinsky_ws_a_ref: VISIBLE_STATE_COUNT_STRUCTURAL * gauge_total,
        strutinsky_ws_ref_nosc: 4.0,
        strutinsky_ws_coulomb_z_ref: 5.0 * (4.0 + 6.0),
        strutinsky_mix: 1.0,
        sigma_z: 5.0 / 2.0,
        sigma_n: 5.0 / 2.0,
        // Structural no-fit lane disables empirical heavy-closure boost;
        // attenuation handles high-shell damping directly.
        proton_magic_weight_coeff: 0.0,
        neutron_magic_weight_coeff: 0.0,
        proton_magic_weight_cap: 1.80,
        neutron_magic_weight_cap: 2.15,
        closure_index_attenuation: 1.0 / 4.0,
        superheavy_proton_amplitude: 2.0,
        superheavy_proton_sigma: 5.0,
        superheavy_proton_gate_n_sigma: 2.0 * gauge_total,
        heavy_target_z,
        heavy_target_n: clifford * VISIBLE_STATE_COUNT_STRUCTURAL + su3,
        heavy_sigma_z: (map.z3_order * map.z3_order) as f64,
        heavy_sigma_n: (map.clifford_dim - 2) as f64,
        heavy_amplitude: 3.0 / 2.0,
        heavy_gate_z_min: (6 * map.clifford_dim) as u16,
        heavy_gate_n_min: (9 * map.clifford_dim) as u16,
        z50_isovector_valley_amplitude: 0.0,
        z50_isovector_beta_coeff: 0.0,
    };

    let nn = NnPotentialParams {
        attractive_depth_mev: shell.strutinsky_ws_depth_mev,
        repulsive_core_mev: semf.a_a,
        range_fm: shell.strutinsky_ws_r0_fm,
        core_radius_fm: shell.strutinsky_ws_diffuseness_fm,
        spin_orbit_mev: shell.strutinsky_spin_orbit_mev,
    };

    StructuralNuclearModel {
        map,
        semf,
        shell,
        nn,
    }
}

/// Structural scan config for periodic-table binding-energy lane (Z <= 118).
pub fn structural_scan_config_z118() -> ScanConfig {
    let model = derive_structural_nuclear_model();
    ScanConfig {
        z_min: 1,
        z_max: 118,
        n_min: 1,
        n_max: 260,
        semf: model.semf,
        shell: model.shell,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn structural_coefficients_are_deterministic() {
        let m = derive_structural_nuclear_model();
        assert!((m.semf.a_v - 95.0 / 6.0).abs() < 1e-12);
        assert!((m.semf.a_s - 55.0 / 3.0).abs() < 1e-12);
        assert!((m.semf.a_c - 2.0 / 3.0).abs() < 1e-12);
        assert!((m.semf.a_a - 23.0).abs() < 1e-12);
        assert!((m.semf.a_p - 12.0).abs() < 1e-12);

        assert!((m.shell.shell_scale_exp - 0.25).abs() < 1e-12);
        assert!((m.shell.strutinsky_ws_depth_mev - 54.0).abs() < 1e-12);
        assert!((m.shell.strutinsky_ws_r0_fm - 1.25).abs() < 1e-12);
        assert!((m.shell.strutinsky_ws_a_ref - 132.0).abs() < 1e-12);

        assert!((m.nn.attractive_depth_mev - 54.0).abs() < 1e-12);
        assert!((m.nn.repulsive_core_mev - 23.0).abs() < 1e-12);
        assert!(m.nn.attractive_depth_mev > 0.0);
        assert!(m.nn.range_fm > 0.0);
    }

    #[test]
    fn structural_targets_match_shared_superheavy_lane() {
        let m = derive_structural_nuclear_model();
        assert!((m.shell.heavy_target_z - 114.0).abs() < 1e-12);
        assert!((m.shell.heavy_target_n - 184.0).abs() < 1e-12);

        let cfg = structural_scan_config_z118();
        assert_eq!(cfg.z_max, 118);
        assert_eq!(cfg.n_max, 260);
    }
}
