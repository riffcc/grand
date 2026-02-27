/*!
 * GUTOE Physics - Perturbation Transfer + BAO/CMB Checks
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-351:
 *   Add a first-principles transfer lane that maps derived background
 *   cosmology + inflation seeds onto observable BAO/CMB transfer quantities.
 */

use std::f64::consts::PI;

/// Speed of light in km/s.
pub const C_KM_S: f64 = 299_792.458;
/// Photon+radiation split factor for `N_eff=3.046`:
/// `Ω_r / Ω_γ ≈ 1 + 0.2271 N_eff ≈ 1.6813`.
pub const OMEGA_R_OVER_OMEGA_GAMMA: f64 = 1.6813;

/// Observational anchors used by the transfer gate.
pub const RS_DRAG_OBS_MPC: f64 = 147.1;
pub const THETA_STAR_OBS_RAD: f64 = 0.01041;
pub const L_PEAK1_OBS: f64 = 220.0;
pub const L_PEAK2_OBS: f64 = 537.0;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransferAssumptions {
    pub h0_km_s_mpc: f64,
    pub omega_b0: f64,
    pub omega_m0: f64,
    pub omega_r0: f64,
    pub omega_k0: f64,
    pub omega_lambda0: f64,
    pub n_s: f64,
    pub a_s: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransferWindows {
    pub rs_rel_max: f64,
    pub theta_star_rel_max: f64,
    pub l1_rel_max: f64,
    pub l2_rel_max: f64,
}

impl Default for TransferWindows {
    fn default() -> Self {
        Self {
            rs_rel_max: 0.08,
            theta_star_rel_max: 0.08,
            l1_rel_max: 0.10,
            l2_rel_max: 0.12,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransferScorecard {
    pub h: f64,
    pub omega_b_h2: f64,
    pub omega_m_h2: f64,
    pub z_drag: f64,
    pub z_recomb: f64,
    pub rs_drag_mpc: f64,
    pub dm_recomb_mpc: f64,
    pub theta_star_rad: f64,
    pub acoustic_scale_la: f64,
    pub l_peak1: f64,
    pub l_peak2: f64,
    pub growth_z0: f64,
    pub growth_z1: f64,
    pub pk_pivot_z0: f64,
    pub pk_pivot_z1: f64,
    pub rs_rel_error: f64,
    pub theta_star_rel_error: f64,
    pub l1_rel_error: f64,
    pub l2_rel_error: f64,
    pub rs_ok: bool,
    pub theta_star_ok: bool,
    pub l1_ok: bool,
    pub l2_ok: bool,
    pub transfer_positive_ok: bool,
}

impl TransferScorecard {
    pub const fn passes_all(self) -> bool {
        self.rs_ok && self.theta_star_ok && self.l1_ok && self.l2_ok && self.transfer_positive_ok
    }
}

fn e2_of_z(z: f64, a: TransferAssumptions) -> f64 {
    let x = 1.0 + z;
    a.omega_r0 * x.powi(4) + a.omega_m0 * x.powi(3) + a.omega_k0 * x.powi(2) + a.omega_lambda0
}

fn hubble_km_s_mpc(z: f64, a: TransferAssumptions) -> f64 {
    let e2 = e2_of_z(z, a);
    if e2 <= 0.0 {
        return f64::NAN;
    }
    a.h0_km_s_mpc * e2.sqrt()
}

fn omega_gamma0(a: TransferAssumptions) -> f64 {
    a.omega_r0 / OMEGA_R_OVER_OMEGA_GAMMA
}

fn baryon_loading_r(z: f64, a: TransferAssumptions) -> f64 {
    let og = omega_gamma0(a);
    if og <= 0.0 {
        return f64::NAN;
    }
    (3.0 * a.omega_b0) / (4.0 * og * (1.0 + z))
}

fn sound_speed_km_s(z: f64, a: TransferAssumptions) -> f64 {
    let r = baryon_loading_r(z, a);
    if !(r.is_finite()) || r <= -1.0 {
        return f64::NAN;
    }
    C_KM_S / (3.0 * (1.0 + r)).sqrt()
}

fn integrate_midpoint_logz<F>(z0: f64, z1: f64, n: usize, mut f: F) -> f64
where
    F: FnMut(f64) -> f64,
{
    if !(z1 > z0) || n == 0 {
        return 0.0;
    }
    let x0 = (1.0 + z0).ln();
    let x1 = (1.0 + z1).ln();
    let dx = (x1 - x0) / n as f64;
    let mut acc = 0.0;
    for i in 0..n {
        let x = x0 + (i as f64 + 0.5) * dx;
        let z = x.exp() - 1.0;
        let weight = (1.0 + z) * dx; // dz
        acc += f(z) * weight;
    }
    acc
}

fn comoving_distance_mpc(z: f64, a: TransferAssumptions) -> f64 {
    integrate_midpoint_logz(0.0, z, 8_192, |zp| {
        let h = hubble_km_s_mpc(zp, a);
        if h <= 0.0 {
            return 0.0;
        }
        C_KM_S / h
    })
}

fn drag_redshift(a: TransferAssumptions) -> f64 {
    let h = a.h0_km_s_mpc / 100.0;
    let omh2 = a.omega_m0 * h * h;
    let obh2 = a.omega_b0 * h * h;

    let b1 = 0.313 * omh2.powf(-0.419) * (1.0 + 0.607 * omh2.powf(0.674));
    let b2 = 0.238 * omh2.powf(0.223);

    1291.0 * omh2.powf(0.251) / (1.0 + 0.659 * omh2.powf(0.828)) * (1.0 + b1 * obh2.powf(b2))
}

fn sound_horizon_mpc(z_drag: f64, a: TransferAssumptions) -> f64 {
    integrate_midpoint_logz(z_drag, 1.0e7, 16_384, |z| {
        let h = hubble_km_s_mpc(z, a);
        if h <= 0.0 {
            return 0.0;
        }
        let cs = sound_speed_km_s(z, a);
        if cs <= 0.0 || !cs.is_finite() {
            return 0.0;
        }
        cs / h
    })
}

/// BBKS transfer function with Sugiyama shape parameter.
fn bbks_transfer(k_h_mpc: f64, gamma_eff: f64) -> f64 {
    if !(k_h_mpc > 0.0) || !(gamma_eff > 0.0) {
        return f64::NAN;
    }
    let q = k_h_mpc / gamma_eff;
    let ln_term = (1.0 + 2.34 * q).ln() / (2.34 * q);
    let poly = 1.0 + 3.89 * q + (16.1 * q).powi(2) + (5.46 * q).powi(3) + (6.71 * q).powi(4);
    ln_term / poly.powf(0.25)
}

/// Carroll, Press, Turner growth suppression factor.
fn growth_suppression_g(z: f64, a: TransferAssumptions) -> f64 {
    let e2 = e2_of_z(z, a);
    if e2 <= 0.0 {
        return f64::NAN;
    }
    let x = 1.0 + z;
    let omega_m_z = a.omega_m0 * x.powi(3) / e2;
    let omega_l_z = a.omega_lambda0 / e2;
    let denom =
        omega_m_z.powf(4.0 / 7.0) - omega_l_z + (1.0 + omega_m_z / 2.0) * (1.0 + omega_l_z / 70.0);
    if denom <= 0.0 {
        return f64::NAN;
    }
    5.0 * omega_m_z / (2.0 * denom)
}

fn normalized_growth_d(z: f64, a: TransferAssumptions) -> f64 {
    let g0 = growth_suppression_g(0.0, a);
    let gz = growth_suppression_g(z, a);
    if !(g0 > 0.0) || !(gz > 0.0) {
        return f64::NAN;
    }
    (gz / (1.0 + z)) / g0
}

/// Dimensionless linear power proxy at wavenumber k (1/Mpc) and redshift z.
fn linear_power_proxy(k_mpc_inv: f64, z: f64, a: TransferAssumptions) -> f64 {
    if !(k_mpc_inv > 0.0) {
        return f64::NAN;
    }
    let h = a.h0_km_s_mpc / 100.0;
    let gamma_eff = a.omega_m0 * h * (-a.omega_b0 * (1.0 + (2.0 * h).sqrt() / a.omega_m0)).exp();
    let t = bbks_transfer(k_mpc_inv / h, gamma_eff);
    let d = normalized_growth_d(z, a);
    if !(t > 0.0) || !(d > 0.0) {
        return f64::NAN;
    }
    let k0 = 0.05; // Mpc^-1 pivot
    a.a_s * (k_mpc_inv / k0).powf(a.n_s - 1.0) * t * t * d * d
}

pub fn evaluate_transfer_gate(a: TransferAssumptions, w: TransferWindows) -> TransferScorecard {
    let h = a.h0_km_s_mpc / 100.0;
    let omega_b_h2 = a.omega_b0 * h * h;
    let omega_m_h2 = a.omega_m0 * h * h;

    let z_drag = drag_redshift(a);
    let z_recomb = 1089.0;

    let rs_drag_mpc = sound_horizon_mpc(z_drag, a);
    let chi_recomb_mpc = comoving_distance_mpc(z_recomb, a);
    // CMB acoustic angle uses comoving angular-diameter distance D_M.
    let dm_recomb_mpc = chi_recomb_mpc;
    let theta_star_rad = rs_drag_mpc / dm_recomb_mpc;
    let acoustic_scale_la = PI / theta_star_rad;

    // First-order phase shift from baryon loading at recombination.
    let r_star = baryon_loading_r(z_recomb, a);
    let peak_shift = (0.73 + 0.04 * r_star).clamp(0.65, 0.85);
    let l_peak1 = acoustic_scale_la * peak_shift;
    let l_peak2 = l_peak1 * 2.44;

    let growth_z0 = normalized_growth_d(0.0, a);
    let growth_z1 = normalized_growth_d(1.0, a);
    let pk_pivot_z0 = linear_power_proxy(0.05, 0.0, a);
    let pk_pivot_z1 = linear_power_proxy(0.05, 1.0, a);

    let rs_rel_error = (rs_drag_mpc - RS_DRAG_OBS_MPC).abs() / RS_DRAG_OBS_MPC;
    let theta_star_rel_error = (theta_star_rad - THETA_STAR_OBS_RAD).abs() / THETA_STAR_OBS_RAD;
    let l1_rel_error = (l_peak1 - L_PEAK1_OBS).abs() / L_PEAK1_OBS;
    let l2_rel_error = (l_peak2 - L_PEAK2_OBS).abs() / L_PEAK2_OBS;

    let rs_ok = rs_rel_error <= w.rs_rel_max;
    let theta_star_ok = theta_star_rel_error <= w.theta_star_rel_max;
    let l1_ok = l1_rel_error <= w.l1_rel_max;
    let l2_ok = l2_rel_error <= w.l2_rel_max;
    let transfer_positive_ok = growth_z0.is_finite()
        && growth_z1.is_finite()
        && growth_z0 > growth_z1
        && pk_pivot_z0.is_finite()
        && pk_pivot_z1.is_finite()
        && pk_pivot_z0 > pk_pivot_z1
        && pk_pivot_z1 > 0.0;

    TransferScorecard {
        h,
        omega_b_h2,
        omega_m_h2,
        z_drag,
        z_recomb,
        rs_drag_mpc,
        dm_recomb_mpc,
        theta_star_rad,
        acoustic_scale_la,
        l_peak1,
        l_peak2,
        growth_z0,
        growth_z1,
        pk_pivot_z0,
        pk_pivot_z1,
        rs_rel_error,
        theta_star_rel_error,
        l1_rel_error,
        l2_rel_error,
        rs_ok,
        theta_star_ok,
        l1_ok,
        l2_ok,
        transfer_positive_ok,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn baseline() -> TransferAssumptions {
        TransferAssumptions {
            h0_km_s_mpc: 67.8,
            omega_b0: 0.0493,
            omega_m0: 0.318,
            omega_r0: 9.0e-5,
            omega_k0: 0.0,
            omega_lambda0: 1.0 - 0.318 - 9.0e-5,
            n_s: 0.965,
            a_s: 2.1e-9,
        }
    }

    #[test]
    fn transfer_observables_are_physical() {
        let s = evaluate_transfer_gate(baseline(), TransferWindows::default());
        assert!(s.z_drag > 900.0 && s.z_drag < 1300.0);
        assert!(s.rs_drag_mpc.is_finite() && s.rs_drag_mpc > 100.0);
        assert!(s.theta_star_rad.is_finite() && s.theta_star_rad > 0.0);
        assert!(s.l_peak1 > 100.0 && s.l_peak2 > s.l_peak1);
        assert!(s.transfer_positive_ok);
    }

    #[test]
    fn transfer_gate_passes_baseline() {
        let s = evaluate_transfer_gate(baseline(), TransferWindows::default());
        assert!(s.passes_all(), "transfer gate failed: {s:#?}");
    }
}
