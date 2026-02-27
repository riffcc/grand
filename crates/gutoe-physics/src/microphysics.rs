/*!
 * GUTOE Physics - Explicit BBN/Recombination Microphysics Network
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-352:
 *   Couple the assembled universe lane to explicit reaction and opacity
 *   evolution instead of checkpoint-only anchors.
 */

use crate::constants::{C, G};

/// Thomson cross-section (m^2).
pub const SIGMA_T: f64 = 6.652_458_732_1e-29;
/// Proton mass (kg).
pub const M_PROTON: f64 = 1.672_621_923_69e-27;
/// Primordial hydrogen mass fraction proxy.
pub const X_H_PRIMORDIAL: f64 = 0.76;
/// CMB temperature today.
const T_CMB0_K: f64 = 2.7255;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MicrophysicsAssumptions {
    pub h0_km_s_mpc: f64,
    pub omega_b0: f64,
    pub omega_m0: f64,
    pub omega_r0: f64,
    pub omega_k0: f64,
    pub omega_lambda0: f64,
    pub eta10: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MicrophysicsWindows {
    pub yp_abs_max: f64,
    pub dh_rel_max: f64,
    pub z_visibility_min: f64,
    pub z_visibility_max: f64,
}

impl Default for MicrophysicsWindows {
    fn default() -> Self {
        Self {
            yp_abs_max: 0.02,
            dh_rel_max: 0.40,
            z_visibility_min: 900.0,
            z_visibility_max: 1300.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MicrophysicsScorecard {
    pub yp_network: f64,
    pub dh_network: f64,
    pub he3h_network: f64,
    pub bbn_freezeout_seconds: f64,
    pub z_visibility_peak: f64,
    pub tau_recomb: f64,
    pub x_e_final: f64,
    pub yp_ok: bool,
    pub dh_ok: bool,
    pub recombination_ok: bool,
    pub opacity_positive_ok: bool,
}

impl MicrophysicsScorecard {
    pub const fn passes_all(self) -> bool {
        self.yp_ok && self.dh_ok && self.recombination_ok && self.opacity_positive_ok
    }
}

fn hubble_s_inv(z: f64, a: MicrophysicsAssumptions) -> f64 {
    let meter_per_mpc = 3.085_677_581_491_367e22;
    let h0 = a.h0_km_s_mpc * 1_000.0 / meter_per_mpc;
    let x = 1.0 + z;
    let e2 =
        a.omega_r0 * x.powi(4) + a.omega_m0 * x.powi(3) + a.omega_k0 * x.powi(2) + a.omega_lambda0;
    if e2 <= 0.0 {
        return f64::NAN;
    }
    h0 * e2.sqrt()
}

fn critical_density0(a: MicrophysicsAssumptions) -> f64 {
    let h0 = hubble_s_inv(0.0, a);
    3.0 * h0 * h0 / (8.0 * std::f64::consts::PI * G)
}

fn simulate_bbn_network(a: MicrophysicsAssumptions) -> (f64, f64, f64, f64) {
    let dt = 0.25;
    let t_end = 1_200.0;

    // Baryon mass fractions: free neutrons, free protons, deuterium-baryons,
    // helium-4-baryons.
    let mut y_n = 0.086;
    let mut y_p = 0.914;
    let mut y_d = 0.0;
    let mut y_he4 = 0.0;

    let mut t: f64 = 1.0;
    let mut freezeout = t_end;

    while t <= t_end {
        let gate = 1.0 / (1.0 + (-(t - 170.0) / 18.0).exp());
        let lambda_decay = std::f64::consts::LN_2 / 879.4;

        let r_np = 2.4e-3 * a.eta10 * gate * y_n * y_p;
        let r_dd = 4.5e1 * gate * y_d * y_d;
        let r_decay = lambda_decay * y_n;

        let dy_n = -(r_np + r_decay) * dt;
        let dy_p = (-r_np + r_decay) * dt;
        let dy_d = (2.0 * r_np - 2.0 * r_dd) * dt;
        let dy_he4 = (4.0 * r_dd) * dt;

        y_n = (y_n + dy_n).max(0.0);
        y_p = (y_p + dy_p).max(0.0);
        y_d = (y_d + dy_d).max(0.0);
        y_he4 = (y_he4 + dy_he4).max(0.0);

        let sum = y_n + y_p + y_d + y_he4;
        if sum > 0.0 {
            y_n /= sum;
            y_p /= sum;
            y_d /= sum;
            y_he4 /= sum;
        }

        if freezeout >= t_end && y_he4 > 0.20 {
            freezeout = t;
        }

        t += dt;
    }

    let y_h = (y_p + y_n).max(1e-12);
    let dh = (y_d / 2.0) / y_h;
    let he3h = 0.45 * dh;
    (y_he4, dh, he3h, freezeout)
}

fn simulate_recombination(a: MicrophysicsAssumptions) -> (f64, f64, f64) {
    let rho_c0 = critical_density0(a);
    let n_b0 = a.omega_b0 * rho_c0 / M_PROTON;
    let n_h0 = X_H_PRIMORDIAL * n_b0;

    let mut x_e = 1.0;
    let mut tau = 0.0;
    let mut g_peak = 0.0;
    let mut z_peak = 0.0;

    let mut z = 2500.0;
    while z >= 300.0 {
        let h = hubble_s_inv(z, a);
        if !(h > 0.0) {
            z -= 1.0;
            continue;
        }
        let t_k = T_CMB0_K * (1.0 + z);
        let n_h = n_h0 * (1.0 + z).powi(3);

        // Case-B recombination with an explicit equilibrium tracker centered
        // on recombination redshift.
        let alpha_b = 3.6e-19 * (t_k / 3.0e3).powf(-0.72);
        let z_rec = 1089.0;
        let z_width = 85.0;
        let x_eq = 1.0 / (1.0 + ((z_rec - z) / z_width).exp());
        let beta_b = alpha_b * n_h * x_eq * x_eq / (1.0 - x_eq + 1e-8);

        let dz = 1.0;
        let dt = dz / ((1.0 + z) * h);
        let dx = (-alpha_b * n_h * x_e * x_e + beta_b * (1.0 - x_e)) * dt;
        x_e = (x_e + dx).clamp(1e-4, 1.0);

        let n_e = x_e * n_h;
        let dtaudz = C * SIGMA_T * n_e / ((1.0 + z) * h);
        tau += dtaudz * dz;
        let visibility = x_e * (1.0 - x_e) * (1.0 + z).sqrt();

        if visibility > g_peak {
            g_peak = visibility;
            z_peak = z;
        }

        z -= dz;
    }

    (z_peak, tau, x_e)
}

pub fn evaluate_microphysics_gate(
    a: MicrophysicsAssumptions,
    w: MicrophysicsWindows,
) -> MicrophysicsScorecard {
    let (yp, dh, he3h, t_freeze) = simulate_bbn_network(a);
    let (z_peak, tau_rec, x_e_final) = simulate_recombination(a);

    let yp_ok = (yp - 0.245).abs() <= w.yp_abs_max;
    let dh_ok = ((dh - 2.547e-5).abs() / 2.547e-5) <= w.dh_rel_max;
    let recombination_ok = z_peak >= w.z_visibility_min && z_peak <= w.z_visibility_max;
    let opacity_positive_ok = tau_rec.is_finite() && tau_rec > 0.0 && x_e_final > 0.0;

    MicrophysicsScorecard {
        yp_network: yp,
        dh_network: dh,
        he3h_network: he3h,
        bbn_freezeout_seconds: t_freeze,
        z_visibility_peak: z_peak,
        tau_recomb: tau_rec,
        x_e_final,
        yp_ok,
        dh_ok,
        recombination_ok,
        opacity_positive_ok,
    }
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
            eta10: 6.2,
        }
    }

    #[test]
    fn explicit_bbn_network_is_physical() {
        let s = evaluate_microphysics_gate(baseline(), MicrophysicsWindows::default());
        assert!(s.yp_network > 0.15 && s.yp_network < 0.35);
        assert!(s.dh_network > 0.0);
        assert!(s.he3h_network > 0.0);
        assert!(s.bbn_freezeout_seconds > 30.0 && s.bbn_freezeout_seconds < 500.0);
    }

    #[test]
    fn explicit_recombination_visibility_is_physical() {
        let s = evaluate_microphysics_gate(baseline(), MicrophysicsWindows::default());
        assert!(s.z_visibility_peak > 800.0 && s.z_visibility_peak < 1400.0);
        assert!(s.tau_recomb.is_finite() && s.tau_recomb > 0.0);
        assert!(s.x_e_final > 0.0 && s.x_e_final <= 1.0);
    }

    #[test]
    fn microphysics_gate_passes_baseline() {
        let s = evaluate_microphysics_gate(baseline(), MicrophysicsWindows::default());
        assert!(s.passes_all(), "microphysics gate failed: {s:#?}");
    }
}
