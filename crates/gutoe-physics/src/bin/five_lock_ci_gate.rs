//! Integrated five-lock CI gate.
//!
//! Locks bundled:
//! 1) single-RGE lock (α₁, α₂, α₃ from one M_Z structural anchor set)
//! 2) neutrino triple lock (Δm²₂₁, Δm²₃₂, ratio + hierarchy/character checks)
//! 3) charged-lepton hierarchy lock (μ/e and τ/μ ratios)
//! 4) dual g-2 lock (electron Schwinger + muon unresolved-gap candidate)
//! 5) BBN three-isotope lock (Yp, D/H, ³He/H)

use gutoe_em::alpha::{
    lepton_masses_from_electron_structural_alpha, ALPHA_INVERSE_PHYSICAL,
    ALPHA_INVERSE_STRUCTURAL,
};
use gutoe_em::{
    electron_mass_from_proton_anchor, neutrino_dirac_majorana_prediction,
    neutrino_hierarchy_prediction, triangulate_neutrino_from_splittings,
};
use gutoe_physics::constants::ALPHA;
use gutoe_physics::{evaluate_bbn_gate, BbnWindows, StandardModelDynamicsMap};
use serde_json::json;
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::Write;
use std::process;

const MZ_GEV: f64 = 91.1876;
const B1: f64 = 41.0 / 10.0;
const B2: f64 = -19.0 / 6.0;
const B3: f64 = -7.0;
const RGE_SPREAD_INV_MAX: f64 = 4.5;

const SOLAR_DM21_TARGET_EV2: f64 = 7.53e-5;
const ATMOSPHERIC_DM32_TARGET_EV2: f64 = 2.453e-3;
const NEUTRINO_SPLIT_REL_MAX: f64 = 1.0e-9;
const NEUTRINO_SUM_CAP_EV: f64 = 0.12;

const MU_OVER_E_RATIO_OBS: f64 = 206.768_283_0;
const TAU_OVER_MU_RATIO_OBS: f64 = 16.816_708_0;
const LEPTON_RATIO_REL_MAX: f64 = 5.0e-3;

const A_E_EXP: f64 = 0.001_159_652_180_59;
const A_MU_EXP: f64 = 0.001_165_920_59;
const A_MU_SM: f64 = 0.001_165_918_10;
const A_E_REL_MAX: f64 = 2.0e-3;
const DELTA_A_MU_REL_MAX: f64 = 1.0e-3;

fn rel_err(observed: f64, target: f64) -> f64 {
    if target.abs() < 1.0e-30 {
        0.0
    } else {
        (observed - target) / target
    }
}

fn alpha_inv_running(alpha_inv_mz: f64, b: f64, mu_gev: f64) -> f64 {
    alpha_inv_mz - (b / (2.0 * PI)) * (mu_gev / MZ_GEV).ln()
}

fn unification_spread_scan(alpha1_mz: f64, alpha2_mz: f64, alpha3_mz: f64) -> (f64, f64, [f64; 3]) {
    let a1_inv_mz = 1.0 / alpha1_mz;
    let a2_inv_mz = 1.0 / alpha2_mz;
    let a3_inv_mz = 1.0 / alpha3_mz;

    let mut best_mu = MZ_GEV;
    let mut best_spread = f64::INFINITY;
    let mut best_vals = [a1_inv_mz, a2_inv_mz, a3_inv_mz];

    for i in 0..=2000 {
        let log10_mu = 2.0 + (17.0 * i as f64) / 2000.0;
        let mu = 10f64.powf(log10_mu);
        let a1 = alpha_inv_running(a1_inv_mz, B1, mu);
        let a2 = alpha_inv_running(a2_inv_mz, B2, mu);
        let a3 = alpha_inv_running(a3_inv_mz, B3, mu);
        let spread = a1.max(a2).max(a3) - a1.min(a2).min(a3);
        if spread < best_spread {
            best_spread = spread;
            best_mu = mu;
            best_vals = [a1, a2, a3];
        }
    }
    (best_mu, best_spread, best_vals)
}

fn main() {
    let out_dir =
        std::env::var("GUTOE_FIVE_LOCK_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/five_lock_ci_gate.json");

    let sm = StandardModelDynamicsMap::from_clifford_z3();
    let sin2_mz = sm.sin2_theta_w_at_mz();
    let cos2_mz = 1.0 - sin2_mz;
    let alpha_em_mz_anchor = ALPHA;
    let alpha1_mz = (5.0 / 3.0) * alpha_em_mz_anchor / cos2_mz;
    let alpha2_mz = alpha_em_mz_anchor / sin2_mz;
    let alpha3_mz = (16.0 / ALPHA_INVERSE_STRUCTURAL) * (67.0 / 66.0);
    let (best_mu_gev, best_spread_inv, best_vals) =
        unification_spread_scan(alpha1_mz, alpha2_mz, alpha3_mz);
    let rge_ok = best_spread_inv <= RGE_SPREAD_INV_MAX;

    let tri = triangulate_neutrino_from_splittings(SOLAR_DM21_TARGET_EV2, ATMOSPHERIC_DM32_TARGET_EV2);
    let dm21_rel = rel_err(tri.dm21_ev2, SOLAR_DM21_TARGET_EV2);
    let dm32_rel = rel_err(tri.dm32_ev2, ATMOSPHERIC_DM32_TARGET_EV2);
    let ratio_target = ATMOSPHERIC_DM32_TARGET_EV2 / SOLAR_DM21_TARGET_EV2;
    let ratio_rel = rel_err(tri.ratio_fit, ratio_target);
    let neutrino_sum_ev = tri.m1_ev + tri.m2_ev + tri.m3_ev;
    let hierarchy = neutrino_hierarchy_prediction();
    let mass_character = neutrino_dirac_majorana_prediction();
    let hierarchy_ok = hierarchy == "normal";
    let mass_character_ok = mass_character == "dirac";
    let split_ok =
        dm21_rel.abs() <= NEUTRINO_SPLIT_REL_MAX
            && dm32_rel.abs() <= NEUTRINO_SPLIT_REL_MAX
            && ratio_rel.abs() <= NEUTRINO_SPLIT_REL_MAX;
    let sum_ok = neutrino_sum_ev <= NEUTRINO_SUM_CAP_EV;
    let neutrino_ok = split_ok && sum_ok && hierarchy_ok && mass_character_ok;

    let me_anchor_mev = electron_mass_from_proton_anchor();
    let [me_pred, mmu_pred, mtau_pred] = lepton_masses_from_electron_structural_alpha(me_anchor_mev);
    let mu_over_e_pred = mmu_pred / me_pred;
    let tau_over_mu_pred = mtau_pred / mmu_pred;
    let mu_over_e_rel = rel_err(mu_over_e_pred, MU_OVER_E_RATIO_OBS);
    let tau_over_mu_rel = rel_err(tau_over_mu_pred, TAU_OVER_MU_RATIO_OBS);
    let lepton_ok =
        mu_over_e_rel.abs() <= LEPTON_RATIO_REL_MAX && tau_over_mu_rel.abs() <= LEPTON_RATIO_REL_MAX;

    let alpha_phys = 1.0 / ALPHA_INVERSE_PHYSICAL;
    let a_e_schwinger = alpha_phys / (2.0 * PI);
    let a_e_rel = rel_err(a_e_schwinger, A_E_EXP);
    let delta_a_mu_ref = A_MU_EXP - A_MU_SM;
    let denom =
        (sm.total_gauge_generators as f64) * ((sm.clifford_dim - sm.magnetic_triplet_card) as f64);
    let delta_a_mu_candidate = alpha_phys.powi(3) / denom;
    let delta_a_mu_rel = rel_err(delta_a_mu_candidate, delta_a_mu_ref);
    let g2_dual_ok = a_e_rel.abs() <= A_E_REL_MAX && delta_a_mu_rel.abs() <= DELTA_A_MU_REL_MAX;

    let bbn = evaluate_bbn_gate(BbnWindows::default());
    let bbn_three_ok = bbn.yp_ok && bbn.dh_ok && bbn.he3_ok;

    let overall_pass = rge_ok && neutrino_ok && lepton_ok && g2_dual_ok && bbn_three_ok;

    let payload = json!({
      "overall_pass": overall_pass,
      "locks": {
        "rge_lock": {
          "pass": rge_ok,
          "spread_inv_max": RGE_SPREAD_INV_MAX,
          "best_spread_inv": best_spread_inv,
          "best_mu_gev": best_mu_gev,
          "inputs_mz": {
            "alpha_em_anchor": alpha_em_mz_anchor,
            "sin2_theta_w_mz": sin2_mz,
            "alpha1_mz": alpha1_mz,
            "alpha2_mz": alpha2_mz,
            "alpha3_mz_structural": alpha3_mz
          },
          "alpha_inv_at_best": {
            "alpha1_inv": best_vals[0],
            "alpha2_inv": best_vals[1],
            "alpha3_inv": best_vals[2]
          }
        },
        "neutrino_triple_lock": {
          "pass": neutrino_ok,
          "windows": {
            "split_rel_max": NEUTRINO_SPLIT_REL_MAX,
            "sum_cap_ev": NEUTRINO_SUM_CAP_EV
          },
          "targets": {
            "dm21_ev2": SOLAR_DM21_TARGET_EV2,
            "dm32_ev2": ATMOSPHERIC_DM32_TARGET_EV2,
            "ratio": ratio_target
          },
          "triangulated": {
            "p_ratio": tri.p_triangulated,
            "kappa_geo": tri.kappa_geo,
            "m1_ev": tri.m1_ev,
            "m2_ev": tri.m2_ev,
            "m3_ev": tri.m3_ev,
            "sum_ev": neutrino_sum_ev,
            "dm21_ev2": tri.dm21_ev2,
            "dm32_ev2": tri.dm32_ev2,
            "ratio": tri.ratio_fit
          },
          "residuals": {
            "dm21_rel": dm21_rel,
            "dm32_rel": dm32_rel,
            "ratio_rel": ratio_rel
          },
          "checks": {
            "split_ok": split_ok,
            "sum_ok": sum_ok,
            "hierarchy_ok": hierarchy_ok,
            "mass_character_ok": mass_character_ok,
            "hierarchy": hierarchy,
            "mass_character": mass_character
          }
        },
        "charged_lepton_hierarchy_lock": {
          "pass": lepton_ok,
          "ratio_rel_max": LEPTON_RATIO_REL_MAX,
          "predicted": {
            "me_mev": me_pred,
            "mmu_mev": mmu_pred,
            "mtau_mev": mtau_pred,
            "mu_over_e": mu_over_e_pred,
            "tau_over_mu": tau_over_mu_pred
          },
          "observed": {
            "mu_over_e": MU_OVER_E_RATIO_OBS,
            "tau_over_mu": TAU_OVER_MU_RATIO_OBS
          },
          "residuals": {
            "mu_over_e_rel": mu_over_e_rel,
            "tau_over_mu_rel": tau_over_mu_rel
          }
        },
        "g2_dual_lock": {
          "pass": g2_dual_ok,
          "windows": {
            "a_e_rel_max": A_E_REL_MAX,
            "delta_a_mu_rel_max": DELTA_A_MU_REL_MAX
          },
          "electron": {
            "a_e_exp": A_E_EXP,
            "a_e_schwinger": a_e_schwinger,
            "a_e_rel_err": a_e_rel
          },
          "muon": {
            "delta_a_mu_ref": delta_a_mu_ref,
            "delta_a_mu_candidate": delta_a_mu_candidate,
            "delta_a_mu_rel_err": delta_a_mu_rel,
            "denominator": denom
          }
        },
        "bbn_three_isotope_lock": {
          "pass": bbn_three_ok,
          "eta10": bbn.eta10,
          "predicted": {
            "yp": bbn.yp_pred,
            "dh": bbn.dh_pred,
            "he3h": bbn.he3h_pred,
            "li7h": bbn.li7h_pred
          },
          "residuals": {
            "yp_delta": bbn.yp_delta,
            "dh_rel_err": bbn.dh_rel_error,
            "he3_rel_err": bbn.he3_rel_error,
            "li_tension_ratio": bbn.li_tension_ratio
          },
          "checks": {
            "yp_ok": bbn.yp_ok,
            "dh_ok": bbn.dh_ok,
            "he3_ok": bbn.he3_ok,
            "li_tension_ok": bbn.li_tension_ok
          }
        }
      }
    });

    let mut file = File::create(&json_path).expect("create five-lock gate json");
    writeln!(
        file,
        "{}",
        serde_json::to_string_pretty(&payload).expect("serialize five-lock gate json")
    )
    .expect("write five-lock gate json");

    println!(
        "five_lock_ci_gate: pass={} rge={} neutrino={} lepton={} g2={} bbn3={}",
        overall_pass, rge_ok, neutrino_ok, lepton_ok, g2_dual_ok, bbn_three_ok
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: rge_ok={} neutrino_ok={} lepton_ok={} g2_dual_ok={} bbn_three_ok={}",
            rge_ok, neutrino_ok, lepton_ok, g2_dual_ok, bbn_three_ok
        );
        process::exit(2);
    }
}
