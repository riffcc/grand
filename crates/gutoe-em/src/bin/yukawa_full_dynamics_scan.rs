//! Full-dynamics Yukawa UV scan.
//!
//! Upgrades beyond the previous one-loop lane:
//! - coupled running of (alpha_s, alpha_em)
//! - 2-loop QCD beta for alpha_s
//! - flavor-sensitive mass flow (QCD + QED + Yukawa self term)
//! - threshold-aware active flavor counting
//!
//! Outputs the same L_g / S_g decomposition over a UV scan.

use gutoe_em::alpha::{z3_extract_params, z3_harmonic_masses};
use gutoe_em::weak::electroweak_vev_from_fermi;
use serde::Serialize;
use std::f64::consts::PI;
use std::fs;

// Reference inputs (MeV) and scales (GeV).
const MU_REF_MEV: f64 = 2.16;
const MD_REF_MEV: f64 = 4.67;
const MS_REF_MEV: f64 = 93.0;
const MC_REF_MEV: f64 = 1270.0;
const MB_REF_MEV: f64 = 4180.0;
const MT_REF_MEV: f64 = 172_760.0;

const MU_REF_GEV: f64 = 2.0;
const MD_REF_GEV: f64 = 2.0;
const MS_REF_GEV: f64 = 2.0;
const MC_REF_GEV: f64 = 1.27;
const MB_REF_GEV: f64 = 4.18;
const MT_REF_GEV: f64 = 172.76;

const MZ_GEV: f64 = 91.1876;
const ALPHA_S_MZ: f64 = 0.118;
const ALPHA_EM_MZ: f64 = 1.0 / 127.95;
const G_F: f64 = 1.166_378_7e-5;

const Q_UP: f64 = 2.0 / 3.0;
const Q_DOWN: f64 = -1.0 / 3.0;

#[derive(Debug, Clone, Serialize)]
struct Z3Fit {
    m_scale: f64,
    s: f64,
    s2: f64,
    delta_deg: f64,
    masses_pred: [f64; 3],
    rms_rel: f64,
}

#[derive(Debug, Clone, Serialize)]
struct FixedSFit {
    s_fixed: f64,
    m_scale: f64,
    delta_deg: f64,
    masses_pred: [f64; 3],
    rms_rel: f64,
}

#[derive(Debug, Clone, Serialize)]
struct LinearFit {
    slope: f64,
    intercept: f64,
    rmse: f64,
    r2: f64,
    y_pred: [f64; 3],
}

#[derive(Debug, Clone, Serialize)]
struct YukawaResponse {
    dmdv_mev_per_gev: [f64; 6], // [u,d,s,c,b,t]
    y_eff: [f64; 6],            // y_eff = sqrt(2) * d(m_GeV)/dv
    d_lgdv_mev_per_gev: [f64; 3],
    d_sgdv_per_gev: [f64; 3],
    eps_rel: f64,
}

#[derive(Debug, Clone, Serialize)]
struct ScaleRow {
    mu_gev: f64,
    alpha_s_mu: f64,
    alpha_em_mu: f64,
    masses_mev: [f64; 6], // [u,d,s,c,b,t]
    lg: [f64; 3],
    sg: [f64; 3],
    lg_z3_free: Z3Fit,
    lg_z3_fixed_s2_2: FixedSFit,
    lg_z3_fixed_s2_3: FixedSFit,
    sg_linear: LinearFit,
    response: YukawaResponse,
}

#[derive(Debug, Clone, Serialize)]
struct Report {
    alpha_s_mz: f64,
    alpha_em_mz: f64,
    vev_gev: f64,
    scan_mu_gev: Vec<f64>,
    rows: Vec<ScaleRow>,
    summary: String,
}

fn nf_active(mu: f64) -> i32 {
    if mu >= MT_REF_GEV {
        6
    } else if mu >= MB_REF_GEV {
        5
    } else if mu >= MC_REF_GEV {
        4
    } else {
        3
    }
}

fn beta0_qcd(nf: i32) -> f64 {
    11.0 - 2.0 * nf as f64 / 3.0
}

fn beta1_qcd(nf: i32) -> f64 {
    102.0 - 38.0 * nf as f64 / 3.0
}

fn sum_nc_q2_active(mu: f64) -> f64 {
    // Active charged leptons.
    let mut s = 0.0;
    if mu >= 0.000_511 {
        s += 1.0; // e
    }
    if mu >= 0.105_658 {
        s += 1.0; // mu
    }
    if mu >= 1.776_86 {
        s += 1.0; // tau
    }

    // Active quarks (N_c Q^2 summed per flavor).
    if mu >= 0.002 {
        s += 3.0 * (2.0 / 3.0) * (2.0 / 3.0); // u
    }
    if mu >= 0.005 {
        s += 3.0 * (1.0 / 3.0) * (1.0 / 3.0); // d
    }
    if mu >= 0.095 {
        s += 3.0 * (1.0 / 3.0) * (1.0 / 3.0); // s
    }
    if mu >= MC_REF_GEV {
        s += 3.0 * (2.0 / 3.0) * (2.0 / 3.0); // c
    }
    if mu >= MB_REF_GEV {
        s += 3.0 * (1.0 / 3.0) * (1.0 / 3.0); // b
    }
    if mu >= MT_REF_GEV {
        s += 3.0 * (2.0 / 3.0) * (2.0 / 3.0); // t
    }
    s
}

fn d_alpha_s_dt(alpha_s: f64, mu: f64) -> f64 {
    let nf = nf_active(mu);
    let b0 = beta0_qcd(nf);
    let b1 = beta1_qcd(nf);
    // t = ln(mu)
    -(b0 / (2.0 * PI)) * alpha_s * alpha_s - (b1 / (4.0 * PI * PI)) * alpha_s * alpha_s * alpha_s
}

fn d_alpha_em_dt(alpha_em: f64, mu: f64) -> f64 {
    let b = (4.0 / 3.0) * sum_nc_q2_active(mu);
    (b / (2.0 * PI)) * alpha_em * alpha_em
}

fn rk4_couplings_step(alpha_s: f64, alpha_em: f64, t: f64, h: f64) -> (f64, f64) {
    let mu = t.exp();

    let k1s = d_alpha_s_dt(alpha_s, mu);
    let k1e = d_alpha_em_dt(alpha_em, mu);

    let as2 = alpha_s + 0.5 * h * k1s;
    let ae2 = alpha_em + 0.5 * h * k1e;
    let mu2 = (t + 0.5 * h).exp();
    let k2s = d_alpha_s_dt(as2, mu2);
    let k2e = d_alpha_em_dt(ae2, mu2);

    let as3 = alpha_s + 0.5 * h * k2s;
    let ae3 = alpha_em + 0.5 * h * k2e;
    let k3s = d_alpha_s_dt(as3, mu2);
    let k3e = d_alpha_em_dt(ae3, mu2);

    let as4 = alpha_s + h * k3s;
    let ae4 = alpha_em + h * k3e;
    let mu4 = (t + h).exp();
    let k4s = d_alpha_s_dt(as4, mu4);
    let k4e = d_alpha_em_dt(ae4, mu4);

    let next_as = alpha_s + (h / 6.0) * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
    let next_ae = alpha_em + (h / 6.0) * (k1e + 2.0 * k2e + 2.0 * k3e + k4e);

    (next_as.max(1e-8), next_ae.max(1e-10))
}

fn couplings_at(mu_target: f64) -> (f64, f64) {
    if (mu_target - MZ_GEV).abs() < 1e-15 {
        return (ALPHA_S_MZ, ALPHA_EM_MZ);
    }
    let t0 = MZ_GEV.ln();
    let t1 = mu_target.ln();
    let span = (t1 - t0).abs();
    let n_steps = ((span * 1200.0).ceil() as usize).max(50);
    let h = (t1 - t0) / n_steps as f64;

    let mut as_ = ALPHA_S_MZ;
    let mut ae = ALPHA_EM_MZ;
    let mut t = t0;
    for _ in 0..n_steps {
        (as_, ae) = rk4_couplings_step(as_, ae, t, h);
        t += h;
    }
    (as_, ae)
}

fn dlnm_dt(m_mev: f64, charge: f64, alpha_s: f64, alpha_em: f64, mu: f64, vev_gev: f64) -> f64 {
    let nf = nf_active(mu);
    let a = alpha_s / PI;
    let gamma1 = 101.0 / 24.0 - 5.0 * nf as f64 / 36.0;

    let gamma_qcd = a + gamma1 * a * a;
    let gamma_qed = 3.0 * charge * charge * alpha_em / (4.0 * PI);

    let m_gev = m_mev / 1000.0;
    let y = (2.0_f64).sqrt() * m_gev / vev_gev;
    // Simple one-loop self-Yukawa reinforcement term.
    let gamma_yuk = 3.0 * y * y / (32.0 * PI * PI);

    -gamma_qcd - gamma_qed + gamma_yuk
}

fn run_mass_full(m_ref_mev: f64, mu_ref: f64, mu_target: f64, charge: f64, vev_gev: f64) -> f64 {
    if (mu_ref - mu_target).abs() < 1e-15 {
        return m_ref_mev;
    }

    let (as_ref, ae_ref) = couplings_at(mu_ref);
    let t0 = mu_ref.ln();
    let t1 = mu_target.ln();
    let span = (t1 - t0).abs();
    let n_steps = ((span * 1400.0).ceil() as usize).max(80);
    let h = (t1 - t0) / n_steps as f64;

    let mut t = t0;
    let mut m = m_ref_mev;
    let mut as_ = as_ref;
    let mut ae = ae_ref;

    for _ in 0..n_steps {
        let mu = t.exp();

        let f_m = |mm: f64, ass: f64, aee: f64, muu: f64| -> f64 {
            mm * dlnm_dt(mm, charge, ass, aee, muu, vev_gev)
        };

        // k1
        let k1m = f_m(m, as_, ae, mu);
        let k1s = d_alpha_s_dt(as_, mu);
        let k1e = d_alpha_em_dt(ae, mu);

        // k2
        let m2 = m + 0.5 * h * k1m;
        let as2 = as_ + 0.5 * h * k1s;
        let ae2 = ae + 0.5 * h * k1e;
        let mu2 = (t + 0.5 * h).exp();
        let k2m = f_m(m2, as2, ae2, mu2);
        let k2s = d_alpha_s_dt(as2, mu2);
        let k2e = d_alpha_em_dt(ae2, mu2);

        // k3
        let m3 = m + 0.5 * h * k2m;
        let as3 = as_ + 0.5 * h * k2s;
        let ae3 = ae + 0.5 * h * k2e;
        let k3m = f_m(m3, as3, ae3, mu2);
        let k3s = d_alpha_s_dt(as3, mu2);
        let k3e = d_alpha_em_dt(ae3, mu2);

        // k4
        let m4 = m + h * k3m;
        let as4 = as_ + h * k3s;
        let ae4 = ae + h * k3e;
        let mu4 = (t + h).exp();
        let k4m = f_m(m4, as4, ae4, mu4);
        let k4s = d_alpha_s_dt(as4, mu4);
        let k4e = d_alpha_em_dt(ae4, mu4);

        m += (h / 6.0) * (k1m + 2.0 * k2m + 2.0 * k3m + k4m);
        as_ += (h / 6.0) * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
        ae += (h / 6.0) * (k1e + 2.0 * k2e + 2.0 * k3e + k4e);

        m = m.max(1e-12);
        as_ = as_.max(1e-8);
        ae = ae.max(1e-10);
        t += h;
    }

    m
}

fn masses_at_mu_with_vev(mu: f64, vev_gev: f64) -> [f64; 6] {
    let m_u = run_mass_full(MU_REF_MEV, MU_REF_GEV, mu, Q_UP, vev_gev);
    let m_d = run_mass_full(MD_REF_MEV, MD_REF_GEV, mu, Q_DOWN, vev_gev);
    let m_s = run_mass_full(MS_REF_MEV, MS_REF_GEV, mu, Q_DOWN, vev_gev);
    let m_c = run_mass_full(MC_REF_MEV, MC_REF_GEV, mu, Q_UP, vev_gev);
    let m_b = run_mass_full(MB_REF_MEV, MB_REF_GEV, mu, Q_DOWN, vev_gev);
    let m_t = run_mass_full(MT_REF_MEV, MT_REF_GEV, mu, Q_UP, vev_gev);
    [m_u, m_d, m_s, m_c, m_b, m_t]
}

fn z3_fit(masses: [f64; 3]) -> Z3Fit {
    let (m_scale, s, delta) = z3_extract_params(masses);
    let pred = z3_harmonic_masses(m_scale, s, delta);
    let rms_rel = {
        let e0 = ((pred[0] - masses[0]) / masses[0]).powi(2);
        let e1 = ((pred[1] - masses[1]) / masses[1]).powi(2);
        let e2 = ((pred[2] - masses[2]) / masses[2]).powi(2);
        ((e0 + e1 + e2) / 3.0).sqrt()
    };
    Z3Fit {
        m_scale,
        s,
        s2: s * s,
        delta_deg: delta.to_degrees(),
        masses_pred: pred,
        rms_rel,
    }
}

fn fit_fixed_s(masses: [f64; 3], s_fixed: f64) -> FixedSFit {
    let a = [masses[0].sqrt(), masses[1].sqrt(), masses[2].sqrt()];
    let mut best_obj = f64::INFINITY;
    let mut best_m = 0.0;
    let mut best_delta = 0.0;
    let n = 30_000usize;

    for i in 0..n {
        let delta = -PI + (2.0 * PI) * (i as f64) / (n as f64);
        let b0 = 1.0 + s_fixed * (delta).cos();
        let b1 = 1.0 + s_fixed * (delta + 2.0 * PI / 3.0).cos();
        let b2 = 1.0 + s_fixed * (delta + 4.0 * PI / 3.0).cos();
        let denom = b0 * b0 + b1 * b1 + b2 * b2;
        if denom <= 1e-14 {
            continue;
        }
        let m = (a[0] * b0 + a[1] * b1 + a[2] * b2) / denom;
        if m <= 0.0 {
            continue;
        }
        let r0 = a[0] - m * b0;
        let r1 = a[1] - m * b1;
        let r2 = a[2] - m * b2;
        let obj = r0 * r0 + r1 * r1 + r2 * r2;
        if obj < best_obj {
            best_obj = obj;
            best_m = m;
            best_delta = delta;
        }
    }

    let pred = z3_harmonic_masses(best_m, s_fixed, best_delta);
    let rms_rel = {
        let e0 = ((pred[0] - masses[0]) / masses[0]).powi(2);
        let e1 = ((pred[1] - masses[1]) / masses[1]).powi(2);
        let e2 = ((pred[2] - masses[2]) / masses[2]).powi(2);
        ((e0 + e1 + e2) / 3.0).sqrt()
    };

    FixedSFit {
        s_fixed,
        m_scale: best_m,
        delta_deg: best_delta.to_degrees(),
        masses_pred: pred,
        rms_rel,
    }
}

fn linear_fit_sg(sg: [f64; 3]) -> LinearFit {
    let x = [1.0_f64, 2.0, 3.0];
    let y = sg;
    let mx = (x[0] + x[1] + x[2]) / 3.0;
    let my = (y[0] + y[1] + y[2]) / 3.0;
    let num = (x[0] - mx) * (y[0] - my) + (x[1] - mx) * (y[1] - my) + (x[2] - mx) * (y[2] - my);
    let den = (x[0] - mx).powi(2) + (x[1] - mx).powi(2) + (x[2] - mx).powi(2);
    let slope = num / den;
    let intercept = my - slope * mx;

    let y_pred = [intercept + slope * x[0], intercept + slope * x[1], intercept + slope * x[2]];
    let mse = ((y[0] - y_pred[0]).powi(2) + (y[1] - y_pred[1]).powi(2) + (y[2] - y_pred[2]).powi(2)) / 3.0;
    let rmse = mse.sqrt();
    let ss_tot = (y[0] - my).powi(2) + (y[1] - my).powi(2) + (y[2] - my).powi(2);
    let ss_res = (y[0] - y_pred[0]).powi(2) + (y[1] - y_pred[1]).powi(2) + (y[2] - y_pred[2]).powi(2);
    let r2 = if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { 1.0 };

    LinearFit {
        slope,
        intercept,
        rmse,
        r2,
        y_pred,
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_YUKAWA_FULL_DYN_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/yukawa_full_dynamics_scan.txt");
    let csv_path = format!("{out_dir}/yukawa_full_dynamics_scan.csv");
    let json_path = format!("{out_dir}/yukawa_full_dynamics_scan.json");

    let vev_gev = electroweak_vev_from_fermi(G_F);
    let eps_rel = 1.0e-3;
    let scan_mu = vec![MT_REF_GEV, 1e4, 1e8, 1e12, 1e16, 1e18, 1e19];

    let mut rows = Vec::new();
    for &mu in &scan_mu {
        let (as_mu, ae_mu) = couplings_at(mu);

        let masses_mev = masses_at_mu_with_vev(mu, vev_gev);
        let [m_u, m_d, m_s, m_c, m_b, m_t] = masses_mev;

        let lg = [
            (m_u * m_d).sqrt(),
            (m_c * m_s).sqrt(),
            (m_t * m_b).sqrt(),
        ];
        let sg = [
            0.5 * (m_u / m_d).ln(),
            0.5 * (m_c / m_s).ln(),
            0.5 * (m_t / m_b).ln(),
        ];

        // Finite-difference Yukawa response around current vev.
        let dv = vev_gev * eps_rel;
        let masses_plus = masses_at_mu_with_vev(mu, vev_gev + dv);
        let masses_minus = masses_at_mu_with_vev(mu, vev_gev - dv);

        let mut dmdv = [0.0_f64; 6];
        let mut y_eff = [0.0_f64; 6];
        for i in 0..6 {
            dmdv[i] = (masses_plus[i] - masses_minus[i]) / (2.0 * dv); // MeV / GeV
            // y_eff = sqrt(2) * d(m_GeV)/dv = sqrt(2) * d(m_MeV)/dv / 1000
            y_eff[i] = (2.0_f64).sqrt() * dmdv[i] / 1000.0;
        }

        let lg_plus = [
            (masses_plus[0] * masses_plus[1]).sqrt(),
            (masses_plus[3] * masses_plus[2]).sqrt(),
            (masses_plus[5] * masses_plus[4]).sqrt(),
        ];
        let lg_minus = [
            (masses_minus[0] * masses_minus[1]).sqrt(),
            (masses_minus[3] * masses_minus[2]).sqrt(),
            (masses_minus[5] * masses_minus[4]).sqrt(),
        ];
        let sg_plus = [
            0.5 * (masses_plus[0] / masses_plus[1]).ln(),
            0.5 * (masses_plus[3] / masses_plus[2]).ln(),
            0.5 * (masses_plus[5] / masses_plus[4]).ln(),
        ];
        let sg_minus = [
            0.5 * (masses_minus[0] / masses_minus[1]).ln(),
            0.5 * (masses_minus[3] / masses_minus[2]).ln(),
            0.5 * (masses_minus[5] / masses_minus[4]).ln(),
        ];
        let d_lgdv = [
            (lg_plus[0] - lg_minus[0]) / (2.0 * dv),
            (lg_plus[1] - lg_minus[1]) / (2.0 * dv),
            (lg_plus[2] - lg_minus[2]) / (2.0 * dv),
        ];
        let d_sgdv = [
            (sg_plus[0] - sg_minus[0]) / (2.0 * dv),
            (sg_plus[1] - sg_minus[1]) / (2.0 * dv),
            (sg_plus[2] - sg_minus[2]) / (2.0 * dv),
        ];

        rows.push(ScaleRow {
            mu_gev: mu,
            alpha_s_mu: as_mu,
            alpha_em_mu: ae_mu,
            masses_mev,
            lg,
            sg,
            lg_z3_free: z3_fit(lg),
            lg_z3_fixed_s2_2: fit_fixed_s(lg, (2.0_f64).sqrt()),
            lg_z3_fixed_s2_3: fit_fixed_s(lg, (3.0_f64).sqrt()),
            sg_linear: linear_fit_sg(sg),
            response: YukawaResponse {
                dmdv_mev_per_gev: dmdv,
                y_eff,
                d_lgdv_mev_per_gev: d_lgdv,
                d_sgdv_per_gev: d_sgdv,
                eps_rel,
            },
        });
    }

    let s2_vals: Vec<f64> = rows.iter().map(|r| r.lg_z3_free.s2).collect();
    let mu_last = *scan_mu.last().unwrap_or(&MT_REF_GEV);
    let summary = format!(
        "full-dynamics scan complete: s2(Lg) at mt={:.6}, at {:.1e}={:.6}; S1 sign at mt={}",
        s2_vals[0],
        mu_last,
        s2_vals[s2_vals.len() - 1],
        if rows[0].sg[0] > 0.0 { "+" } else { "-" }
    );

    let report = Report {
        alpha_s_mz: ALPHA_S_MZ,
        alpha_em_mz: ALPHA_EM_MZ,
        vev_gev,
        scan_mu_gev: scan_mu,
        rows,
        summary,
    };

    // TXT
    let mut txt = String::new();
    txt.push_str("[yukawa_full_dynamics_scan]\n");
    txt.push_str(&format!("alpha_s_mz = {:.9}\n", report.alpha_s_mz));
    txt.push_str(&format!("alpha_em_mz = {:.12}\n", report.alpha_em_mz));
    txt.push_str(&format!("vev_gev = {:.9}\n\n", report.vev_gev));

    for row in &report.rows {
        txt.push_str(&format!("[mu = {:.6e} GeV]\n", row.mu_gev));
        txt.push_str(&format!(
            "alpha_s={:.9} alpha_em={:.12}\n",
            row.alpha_s_mu, row.alpha_em_mu
        ));
        txt.push_str(&format!(
            "masses_mev [u,d,s,c,b,t] = [{:.6}, {:.6}, {:.6}, {:.6}, {:.6}, {:.6}]\n",
            row.masses_mev[0],
            row.masses_mev[1],
            row.masses_mev[2],
            row.masses_mev[3],
            row.masses_mev[4],
            row.masses_mev[5]
        ));
        txt.push_str(&format!("L_g = [{:.6}, {:.6}, {:.6}]\n", row.lg[0], row.lg[1], row.lg[2]));
        txt.push_str(&format!("S_g = [{:.6}, {:.6}, {:.6}]\n", row.sg[0], row.sg[1], row.sg[2]));
        txt.push_str(&format!(
            "L_free: s2={:.9} delta={:.6} rms={:.3e}\n",
            row.lg_z3_free.s2, row.lg_z3_free.delta_deg, row.lg_z3_free.rms_rel
        ));
        txt.push_str(&format!(
            "L_fixed(s2=2): rms={:.6}\n",
            row.lg_z3_fixed_s2_2.rms_rel
        ));
        txt.push_str(&format!(
            "L_fixed(s2=3): rms={:.6}\n",
            row.lg_z3_fixed_s2_3.rms_rel
        ));
        txt.push_str(&format!(
            "S_linear: slope={:.9} intercept={:.9} rmse={:.9} r2={:.9}\n\n",
            row.sg_linear.slope,
            row.sg_linear.intercept,
            row.sg_linear.rmse,
            row.sg_linear.r2
        ));
        txt.push_str(&format!(
            "dm/dv [u,d,s,c,b,t] (MeV/GeV) = [{:.6}, {:.6}, {:.6}, {:.6}, {:.6}, {:.6}]\n",
            row.response.dmdv_mev_per_gev[0],
            row.response.dmdv_mev_per_gev[1],
            row.response.dmdv_mev_per_gev[2],
            row.response.dmdv_mev_per_gev[3],
            row.response.dmdv_mev_per_gev[4],
            row.response.dmdv_mev_per_gev[5]
        ));
        txt.push_str(&format!(
            "y_eff [u,d,s,c,b,t] = [{:.9}, {:.9}, {:.9}, {:.9}, {:.9}, {:.9}]\n",
            row.response.y_eff[0],
            row.response.y_eff[1],
            row.response.y_eff[2],
            row.response.y_eff[3],
            row.response.y_eff[4],
            row.response.y_eff[5]
        ));
        txt.push_str(&format!(
            "dL/dv [L1,L2,L3] = [{:.6}, {:.6}, {:.6}]  dS/dv [S1,S2,S3] = [{:.6e}, {:.6e}, {:.6e}] (eps={:.1e})\n\n",
            row.response.d_lgdv_mev_per_gev[0],
            row.response.d_lgdv_mev_per_gev[1],
            row.response.d_lgdv_mev_per_gev[2],
            row.response.d_sgdv_per_gev[0],
            row.response.d_sgdv_per_gev[1],
            row.response.d_sgdv_per_gev[2],
            row.response.eps_rel
        ));
    }
    txt.push_str(&format!("summary = {}\n", report.summary));
    fs::write(&txt_path, txt).expect("write txt");

    // CSV
    let mut csv = String::new();
    csv.push_str("mu_gev,alpha_s_mu,alpha_em_mu,m_u,m_d,m_s,m_c,m_b,m_t,L1,L2,L3,S1,S2,S3,L_free_s2,L_fixed2_rms,L_fixed3_rms,S_lin_rmse,S_lin_r2,dm_dv_u,dm_dv_d,dm_dv_s,dm_dv_c,dm_dv_b,dm_dv_t,y_eff_u,y_eff_d,y_eff_s,y_eff_c,y_eff_b,y_eff_t,dLdv_1,dLdv_2,dLdv_3,dSdv_1,dSdv_2,dSdv_3,eps_rel\n");
    for r in &report.rows {
        csv.push_str(&format!(
            "{:.6e},{:.9},{:.12},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9e},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.12},{:.12},{:.12},{:.12},{:.12},{:.12},{:.9},{:.9},{:.9},{:.12e},{:.12e},{:.12e},{:.1e}\n",
            r.mu_gev,
            r.alpha_s_mu,
            r.alpha_em_mu,
            r.masses_mev[0], r.masses_mev[1], r.masses_mev[2], r.masses_mev[3], r.masses_mev[4], r.masses_mev[5],
            r.lg[0], r.lg[1], r.lg[2],
            r.sg[0], r.sg[1], r.sg[2],
            r.lg_z3_free.s2,
            r.lg_z3_fixed_s2_2.rms_rel,
            r.lg_z3_fixed_s2_3.rms_rel,
            r.sg_linear.rmse,
            r.sg_linear.r2,
            r.response.dmdv_mev_per_gev[0], r.response.dmdv_mev_per_gev[1], r.response.dmdv_mev_per_gev[2],
            r.response.dmdv_mev_per_gev[3], r.response.dmdv_mev_per_gev[4], r.response.dmdv_mev_per_gev[5],
            r.response.y_eff[0], r.response.y_eff[1], r.response.y_eff[2],
            r.response.y_eff[3], r.response.y_eff[4], r.response.y_eff[5],
            r.response.d_lgdv_mev_per_gev[0], r.response.d_lgdv_mev_per_gev[1], r.response.d_lgdv_mev_per_gev[2],
            r.response.d_sgdv_per_gev[0], r.response.d_sgdv_per_gev[1], r.response.d_sgdv_per_gev[2],
            r.response.eps_rel
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    fs::write(&json_path, serde_json::to_string_pretty(&report).expect("serialize"))
        .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {csv_path}");
    println!("wrote {json_path}");
    println!("{}", report.summary);
}
