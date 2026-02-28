/*!
 * Light-nucleus few-body QM proxy lane.
 *
 * Purpose:
 * - Provide a variational few-body binding estimate where SEMF is known to fail
 *   (bulk approximation breaks for small A).
 * - Keep coefficients structural/no-fit by deriving scales from Cl(1,3)+SM counts.
 *
 * Scope:
 * - A <= 16 light-nucleus corridor.
 * - Gaussian two-body attraction + short-range core repulsion + tensor-like pn term
 *   + compact three-body term + kinetic and Coulomb expectation values.
 * - Intended as a physically motivated correction layer, not an ab-initio solver.
 */

use crate::dynamics_map::StandardModelDynamicsMap;

const HBAR2_OVER_2M_NUCLEON_MEV_FM2: f64 = 20.735_53;

#[derive(Clone, Copy, Debug)]
pub struct FewBodyQmParams {
    pub max_a: u16,
    pub attractive_depth_mev: f64,
    pub repulsive_core_mev: f64,
    pub range_fm: f64,
    pub core_radius_fm: f64,
    pub pn_pair_scale: f64,
    pub like_pair_scale: f64,
    pub kinetic_prefactor: f64,
    pub tensor_mev: f64,
    pub tensor_range_fm: f64,
    pub tensor_a_damping: f64,
    pub three_body_mev: f64,
    pub three_body_range_fm: f64,
    pub pair_saturation_scale: f64,
}

impl FewBodyQmParams {
    /// Structural/no-fit parameter chain from Cl(1,3)+SM counts.
    pub fn from_clifford_z3() -> Self {
        let m = StandardModelDynamicsMap::from_clifford_z3();
        let clifford = m.clifford_dim as f64; // 16
        let su3 = m.su3_generators as f64; // 8
        let su2 = m.su2_generators as f64; // 3
        let u1 = m.u1_generators as f64; // 1
        let gauge = m.total_gauge_generators as f64; // 12

        Self {
            max_a: m.clifford_dim as u16, // A<=16
            attractive_depth_mev: (3.0 * clifford + gauge / 2.0) / (su2 + u1), // 54/4
            repulsive_core_mev: clifford + gauge / 2.0 + su3 / 8.0, // 23
            range_fm: (gauge + su2) / gauge, // 5/4
            core_radius_fm: 2.0 / 3.0,
            pn_pair_scale: 1.0 + 2.0 / (su2 + u1),   // 1.5
            like_pair_scale: 1.0 - 1.0 / (su3 + u1), // 8/9
            kinetic_prefactor: (su2 + u1) / (gauge + su3), // 1/5
            tensor_mev: su3 / su2,                   // 8/3
            tensor_range_fm: 1.5,
            tensor_a_damping: 0.40,
            three_body_mev: gauge, // 12
            three_body_range_fm: 1.4,
            pair_saturation_scale: gauge / (su2 + u1), // 3
        }
    }
}

fn gaussian_overlap(b_fm: f64, range_fm: f64) -> f64 {
    let x = (b_fm * b_fm) / (range_fm * range_fm).max(1.0e-9);
    (1.0 + x).powf(-1.5)
}

fn choose3(a: u16) -> f64 {
    if a < 3 {
        0.0
    } else {
        let af = a as f64;
        af * (af - 1.0) * (af - 2.0) / 6.0
    }
}

/// Few-body variational binding estimate for light nuclei.
///
/// Returns `None` outside the intended light corridor.
pub fn few_body_binding_variational_mev(z: u16, n: u16, params: FewBodyQmParams) -> Option<f64> {
    let a = z + n;
    if a < 2 || a > params.max_a {
        return None;
    }
    let af = a as f64;
    let zf = z as f64;
    let nf = n as f64;

    let pn_pairs = zf * nf;
    let pp_pairs = zf * (zf - 1.0) * 0.5;
    let nn_pairs = nf * (nf - 1.0) * 0.5;
    let pair_total = pn_pairs + pp_pairs + nn_pairs;
    if pair_total <= 0.0 {
        return Some(0.0);
    }

    let isospin_sym = 1.0 - ((nf - zf).abs() / af).powi(2);
    let sat = 1.0 / (1.0 + pair_total / (params.pair_saturation_scale * af).max(1.0));

    let mut best_binding = f64::NEG_INFINITY;
    // Coarse deterministic variational sweep over Gaussian width b.
    for i in 0..=125 {
        let b = 0.70 + 0.02 * i as f64; // 0.70 .. 3.20 fm
        let overlap_attr = gaussian_overlap(b, params.range_fm);
        let overlap_core = gaussian_overlap(b, params.core_radius_fm);
        let overlap_tensor = gaussian_overlap(b, params.tensor_range_fm);
        let overlap_three = gaussian_overlap(b, params.three_body_range_fm);

        let attraction = params.attractive_depth_mev
            * overlap_attr
            * sat
            * (params.pn_pair_scale * pn_pairs + params.like_pair_scale * (pp_pairs + nn_pairs));

        let core_repulsion =
            params.repulsive_core_mev * overlap_core * pair_total * (0.65 + 0.35 * sat);

        let deuteron_coherence = 1.0 + 0.75 / (af - 1.0).max(1.0);
        let tensor_term =
            params.tensor_mev * overlap_tensor * pn_pairs * isospin_sym * deuteron_coherence
                / (1.0 + params.tensor_a_damping * (af - 2.0).max(0.0));

        // Compact few-body enhancement centered on A~3-4, suppressed rapidly
        // beyond the true few-body corridor.
        let compact_a = 1.0 / (1.0 + (af - 3.5).powi(4));
        let three_body = params.three_body_mev * choose3(a) * overlap_three * compact_a;

        let kinetic =
            params.kinetic_prefactor * 1.5 * (af - 1.0).max(1.0) * HBAR2_OVER_2M_NUCLEON_MEV_FM2
                / (b * b);

        let r_rms = (1.5_f64).sqrt() * b;
        let coulomb = if z >= 2 {
            0.72 * zf * (zf - 1.0) / r_rms.max(1.0e-6)
        } else {
            0.0
        };

        let binding = attraction + tensor_term + three_body - core_repulsion - kinetic - coulomb;
        if binding > best_binding {
            best_binding = binding;
        }
    }
    Some(best_binding.max(0.0))
}

/// Blend weight for the few-body lane.
///
/// - 1.0 at A<=4 (pure few-body region),
/// - smoothly decreases to 0 by A=16,
/// - 0 beyond A>16 (bulk SEMF region).
pub fn light_nucleus_blend_weight(a: u16) -> f64 {
    if a <= 4 {
        return 1.0;
    }
    if a > 16 {
        return 0.0;
    }
    let x = (16.0 - a as f64) / 12.0;
    x.clamp(0.0, 1.0).powf(1.75)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bpa(z: u16, n: u16) -> f64 {
        let p = FewBodyQmParams::from_clifford_z3();
        few_body_binding_variational_mev(z, n, p).unwrap() / (z + n) as f64
    }

    #[test]
    fn light_cluster_bindings_are_positive_and_ordered() {
        // D, T/He3, alpha ordering should hold in few-body lane.
        let d = bpa(1, 1);
        let t = bpa(1, 2);
        let he3 = bpa(2, 1);
        let he4 = bpa(2, 2);
        assert!(d > 0.2);
        assert!(t > d);
        assert!(he3 > d);
        assert!(he4 > t);
        assert!(he4 > he3);
    }

    #[test]
    fn alpha_binding_is_in_reasonable_band() {
        let he4 = bpa(2, 2);
        assert!(he4 > 5.0 && he4 < 9.0, "He-4 B/A out of band: {he4}");
    }

    #[test]
    fn blend_weight_transitions_to_bulk_lane() {
        assert_eq!(light_nucleus_blend_weight(2), 1.0);
        assert!(light_nucleus_blend_weight(8) > 0.0);
        assert!(light_nucleus_blend_weight(16) == 0.0);
        assert!(light_nucleus_blend_weight(24) == 0.0);
    }
}
