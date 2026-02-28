/*!
 * GUTOE Physics - Derived Reionization Optical Depth
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-343 / GRAND-355:
 *   Derive tau_reio from structure-formation timing and background expansion
 *   using shared Cl(1,3) primitives (no CMB-spectrum fitting).
 */

use crate::bbn::ETA10_REF;
use crate::constants::{
    BIVECTOR_TIMELIKE_SPACELIKE_COUNT, BIVECTOR_TOTAL_COUNT, C, DARK_GEOMETRIC_AMPLIFICATION,
    GRADE1_STATE_COUNT_STRUCTURAL,
};
use crate::microphysics::{MicrophysicsAssumptions, M_PROTON, SIGMA_T, X_H_PRIMORDIAL};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReionizationDerived {
    pub z_reion_structural: f64,
    pub tau_reio: f64,
    pub suppression_e2tau: f64,
}

fn e2_of_z(z: f64, a: MicrophysicsAssumptions) -> f64 {
    let x = 1.0 + z;
    a.omega_r0 * x.powi(4) + a.omega_m0 * x.powi(3) + a.omega_k0 * x.powi(2) + a.omega_lambda0
}

fn hubble_s_inv(z: f64, a: MicrophysicsAssumptions) -> f64 {
    let meter_per_mpc = 3.085_677_581_491_367e22;
    let h0 = a.h0_km_s_mpc * 1_000.0 / meter_per_mpc;
    let e2 = e2_of_z(z, a);
    if e2 <= 0.0 {
        return f64::NAN;
    }
    h0 * e2.sqrt()
}

fn critical_density0(a: MicrophysicsAssumptions) -> f64 {
    let h0 = hubble_s_inv(0.0, a);
    3.0 * h0 * h0 / (8.0 * std::f64::consts::PI * crate::constants::G)
}

/// Structural reionization redshift proxy from shared Cl(1,3) counts.
///
/// Base factor:
///   (12 - 4) * (10/9) = 8.888...
/// where
///   12 = non-grade1 channel count,
///   4 = grade1 count,
///   10 = grade1+grade2,
///   9 = grade2+SU(2) triplet.
///
/// This is then scaled by the baryogenesis-derived eta10 ratio to keep the
/// timing coupled to the upstream matter asymmetry lane.
pub fn structural_reionization_redshift(eta10: f64) -> f64 {
    if eta10 <= 0.0 {
        return f64::NAN;
    }
    let base = (DARK_GEOMETRIC_AMPLIFICATION - GRADE1_STATE_COUNT_STRUCTURAL)
        * ((GRADE1_STATE_COUNT_STRUCTURAL + BIVECTOR_TOTAL_COUNT)
            / (BIVECTOR_TOTAL_COUNT + BIVECTOR_TIMELIKE_SPACELIKE_COUNT));
    let eta_factor = (eta10 / ETA10_REF).powf(1.0 / 3.0);
    base * eta_factor
}

/// Integrate optical depth from z=0 to structural z_reion, with fully ionized
/// hydrogen and a first-order helium electron correction.
pub fn derive_tau_reio(
    a: MicrophysicsAssumptions,
    eta10: f64,
) -> Result<ReionizationDerived, String> {
    let z_reion = structural_reionization_redshift(eta10);
    if !z_reion.is_finite() || z_reion <= 0.0 {
        return Err(format!("invalid structural z_reion: {z_reion}"));
    }

    let rho_c0 = critical_density0(a);
    let n_b0 = a.omega_b0 * rho_c0 / M_PROTON;
    let n_h0 = X_H_PRIMORDIAL * n_b0;

    // First-order helium electron correction from primordial Yp ~ 0.245.
    let helium_factor = 1.0 + 0.245 / (4.0 * X_H_PRIMORDIAL);

    let n = 8192usize;
    let dz = z_reion / n as f64;
    let mut tau = 0.0;
    for i in 0..n {
        let z = (i as f64 + 0.5) * dz;
        let h = hubble_s_inv(z, a);
        if !(h > 0.0) {
            continue;
        }
        let n_e = helium_factor * n_h0 * (1.0 + z).powi(3);
        let dtaudz = C * SIGMA_T * n_e / ((1.0 + z) * h);
        tau += dtaudz * dz;
    }

    Ok(ReionizationDerived {
        z_reion_structural: z_reion,
        tau_reio: tau,
        suppression_e2tau: (-2.0 * tau).exp(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn baseline() -> MicrophysicsAssumptions {
        MicrophysicsAssumptions {
            h0_km_s_mpc: 68.0,
            omega_b0: 0.0493,
            omega_m0: 0.318,
            omega_r0: 9.0e-5,
            omega_k0: 0.0,
            omega_lambda0: 1.0 - 0.318 - 9.0e-5,
            eta10: 6.0,
        }
    }

    #[test]
    fn structural_reionization_redshift_is_physical() {
        let z = structural_reionization_redshift(6.0);
        assert!(z > 6.0 && z < 12.0, "z_reion={z}");
    }

    #[test]
    fn derived_tau_reio_is_physical() {
        let d = derive_tau_reio(baseline(), 6.0).expect("derive tau");
        assert!(d.z_reion_structural > 6.0 && d.z_reion_structural < 12.0);
        assert!(d.tau_reio > 0.02 && d.tau_reio < 0.12, "tau={}", d.tau_reio);
        assert!(d.suppression_e2tau > 0.7 && d.suppression_e2tau < 1.0);
    }
}
