/*!
 * GUTOE Physics - Microphysics-Derived CMB Damping Envelope
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-344 / GRAND-355:
 *   Derive a Silk/visibility damping envelope from in-framework microphysics
 *   (no CMB fit knobs), then project it into multipole space.
 */

use crate::cmb_class::ClassTtPoint;
use crate::constants::{C, G};
use crate::microphysics::{MicrophysicsAssumptions, M_PROTON, SIGMA_T, X_H_PRIMORDIAL};

/// Speed of light in km/s.
const C_KM_S: f64 = 299_792.458;
/// Radiation to photon ratio for Neff ~ 3.046.
const OMEGA_R_OVER_OMEGA_GAMMA: f64 = 1.6813;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SilkDampingDerived {
    pub z_star: f64,
    pub sigma_z: f64,
    pub diffusion_length_mpc: f64,
    pub visibility_width_mpc: f64,
    pub k_diff_mpc_inv: f64,
    pub ell_diff: f64,
    pub ell_vis: f64,
    pub d_m_star_mpc: f64,
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

fn hubble_km_s_mpc(z: f64, a: MicrophysicsAssumptions) -> f64 {
    let e2 = e2_of_z(z, a);
    if e2 <= 0.0 {
        return f64::NAN;
    }
    a.h0_km_s_mpc * e2.sqrt()
}

fn critical_density0(a: MicrophysicsAssumptions) -> f64 {
    let h0 = hubble_s_inv(0.0, a);
    3.0 * h0 * h0 / (8.0 * std::f64::consts::PI * G)
}

fn omega_gamma0(a: MicrophysicsAssumptions) -> f64 {
    a.omega_r0 / OMEGA_R_OVER_OMEGA_GAMMA
}

fn baryon_loading_r(z: f64, a: MicrophysicsAssumptions) -> f64 {
    let og = omega_gamma0(a);
    if og <= 0.0 {
        return f64::NAN;
    }
    (3.0 * a.omega_b0) / (4.0 * og * (1.0 + z))
}

fn comoving_distance_mpc(z: f64, a: MicrophysicsAssumptions) -> f64 {
    if z <= 0.0 {
        return 0.0;
    }
    let n = 8_192usize;
    let x0 = (1.0f64).ln();
    let x1 = (1.0 + z).ln();
    let dx = (x1 - x0) / n as f64;
    let mut acc = 0.0;
    for i in 0..n {
        let x = x0 + (i as f64 + 0.5) * dx;
        let zp = x.exp() - 1.0;
        let dz = (1.0 + zp) * dx;
        let h = hubble_km_s_mpc(zp, a);
        if h > 0.0 {
            acc += (C_KM_S / h) * dz;
        }
    }
    acc
}

/// Derive Silk diffusion and visibility-width scales directly from the
/// in-framework recombination microphysics.
pub fn derive_silk_damping_envelope(a: MicrophysicsAssumptions) -> Result<SilkDampingDerived, String> {
    let rho_c0 = critical_density0(a);
    let n_b0 = a.omega_b0 * rho_c0 / M_PROTON;
    let n_h0 = X_H_PRIMORDIAL * n_b0;

    let mut x_e = 1.0;
    let mut z = 2_500.0;
    let dz = 1.0;

    let mut samples: Vec<(f64, f64, f64, f64)> = Vec::new(); // (z, visibility, n_e, h)
    let mut z_star = 0.0;
    let mut g_star = -1.0_f64;

    while z >= 300.0 {
        let h = hubble_s_inv(z, a);
        if !(h > 0.0) {
            z -= dz;
            continue;
        }
        let t_k = 2.7255 * (1.0 + z);
        let n_h = n_h0 * (1.0 + z).powi(3);

        let alpha_b = 3.6e-19 * (t_k / 3.0e3).powf(-0.72);
        let z_rec = 1089.0;
        let z_width = 85.0;
        let x_eq = 1.0 / (1.0 + ((z_rec - z) / z_width).exp());
        let beta_b = alpha_b * n_h * x_eq * x_eq / (1.0 - x_eq + 1e-8);

        let dt = dz / ((1.0 + z) * h);
        let dx = (-alpha_b * n_h * x_e * x_e + beta_b * (1.0 - x_e)) * dt;
        x_e = (x_e + dx).clamp(1e-4, 1.0);

        let n_e = x_e * n_h;
        let g_vis = x_e * (1.0 - x_e) * (1.0 + z).sqrt();
        if g_vis > g_star {
            g_star = g_vis;
            z_star = z;
        }
        samples.push((z, g_vis, n_e, h));

        z -= dz;
    }

    if !(g_star > 0.0) {
        return Err("visibility integral is non-positive".to_string());
    }

    // Visibility width around decoupling peak: ignore tiny tails that are
    // numerically noisy and physically irrelevant to line-of-sight weighting.
    let vis_cut = 0.01 * g_star;
    let mut sum_g = 0.0;
    let mut sum_gz2 = 0.0;
    for (zz, gg, _, _) in &samples {
        if *gg >= vis_cut {
            sum_g += *gg;
            sum_gz2 += *gg * *zz * *zz;
        }
    }
    if !(sum_g > 0.0) {
        return Err("visibility support is empty after thresholding".to_string());
    }

    // Keep z* as the actual visibility peak, not the weighted mean.
    let var_z = (sum_gz2 / sum_g) - z_star * z_star;
    let sigma_z = var_z.max(0.0).sqrt();

    // Diffusion integral (k_D^{-2}) in SI units, only up to decoupling.
    // Integrating past decoupling overestimates damping in TT.
    let mut inv_kd2_m2 = 0.0;
    for (zz, _, n_e, h) in &samples {
        if *zz < z_star {
            continue;
        }
        let r = baryon_loading_r(*zz, a);
        if *n_e > 0.0 && r.is_finite() && r > -1.0 {
            let bracket = (r * r + (16.0 / 15.0) * (1.0 + r)) / ((1.0 + r) * (1.0 + r));
            let pref = (C / *h) / (6.0 * *n_e * SIGMA_T);
            inv_kd2_m2 += pref * bracket * dz;
        }
    }
    if !(inv_kd2_m2 > 0.0) {
        return Err("diffusion integral is non-positive".to_string());
    }

    let mpc_m = 3.085_677_581_491_367e22;
    let diffusion_length_mpc = inv_kd2_m2.sqrt() / mpc_m;
    let k_diff_mpc_inv = 1.0 / diffusion_length_mpc;

    let d_m_star_mpc = comoving_distance_mpc(z_star, a);
    let ell_diff = k_diff_mpc_inv * d_m_star_mpc;

    let h_star = hubble_km_s_mpc(z_star, a);
    if !(h_star > 0.0) {
        return Err("H(z*) is not positive".to_string());
    }
    let dchi_dz_star = C_KM_S / h_star;
    let visibility_width_mpc = sigma_z * dchi_dz_star;
    let ell_vis = if visibility_width_mpc > 0.0 {
        d_m_star_mpc / visibility_width_mpc
    } else {
        f64::INFINITY
    };

    Ok(SilkDampingDerived {
        z_star,
        sigma_z,
        diffusion_length_mpc,
        visibility_width_mpc,
        k_diff_mpc_inv,
        ell_diff,
        ell_vis,
        d_m_star_mpc,
    })
}

/// Apply a microphysics-derived damping envelope to an existing TT spectrum.
///
/// This is intentionally parameter-free: all scale inputs come from
/// `derive_silk_damping_envelope`.
pub fn apply_microphysics_damping(
    class_tt: &[ClassTtPoint],
    d: SilkDampingDerived,
) -> Vec<ClassTtPoint> {
    class_tt
        .iter()
        .map(|p| {
            let ell = p.ell as f64;
            let e_diff = if d.ell_diff.is_finite() && d.ell_diff > 0.0 {
                (-(ell / d.ell_diff).powi(2)).exp()
            } else {
                1.0
            };
            let e_vis = if d.ell_vis.is_finite() && d.ell_vis > 0.0 {
                (-(ell / d.ell_vis).powi(2)).exp()
            } else {
                1.0
            };
            ClassTtPoint {
                ell: p.ell,
                d_ell_tt_uk2: p.d_ell_tt_uk2 * e_diff * e_vis,
            }
        })
        .collect()
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
    fn derived_scales_are_physical() {
        let d = derive_silk_damping_envelope(baseline()).expect("derive scales");
        assert!(d.z_star > 850.0 && d.z_star < 1300.0, "z*= {}", d.z_star);
        assert!(d.sigma_z > 10.0 && d.sigma_z < 250.0, "sigma_z= {}", d.sigma_z);
        assert!(d.diffusion_length_mpc > 0.1 && d.diffusion_length_mpc < 80.0);
        assert!(d.ell_diff > 100.0 && d.ell_diff < 6000.0);
        assert!(d.ell_vis > 100.0 && d.ell_vis < 8000.0);
    }

    #[test]
    fn damping_is_multiplicative_and_positive() {
        let d = SilkDampingDerived {
            z_star: 1060.0,
            sigma_z: 100.0,
            diffusion_length_mpc: 5.0,
            visibility_width_mpc: 8.0,
            k_diff_mpc_inv: 0.2,
            ell_diff: 1500.0,
            ell_vis: 2200.0,
            d_m_star_mpc: 14_000.0,
        };
        let src = vec![
            ClassTtPoint {
                ell: 50,
                d_ell_tt_uk2: 1_000.0,
            },
            ClassTtPoint {
                ell: 1000,
                d_ell_tt_uk2: 2_000.0,
            },
        ];
        let out = apply_microphysics_damping(&src, d);
        assert_eq!(out.len(), src.len());
        assert!(out[0].d_ell_tt_uk2 > 0.0);
        assert!(out[1].d_ell_tt_uk2 > 0.0);
        assert!(out[1].d_ell_tt_uk2 < src[1].d_ell_tt_uk2);
    }
}
