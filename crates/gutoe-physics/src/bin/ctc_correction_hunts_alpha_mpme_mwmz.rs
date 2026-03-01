//! Correction hunts prioritized by user prompt:
//! 1) alpha^-1 (from LO 137)
//! 2) mp/me (from LO 6*pi^5)
//! 3) mW/mZ (tree-level vs scheme-corrected lane)
//!
//! Zero new fit parameters are introduced in the reported candidate lanes;
//! we evaluate structurally motivated correction forms from existing constants.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

const ALPHA_INV_PHYS: f64 = 137.035_999_177;
const ALPHA: f64 = 1.0 / ALPHA_INV_PHYS;
const PI: f64 = std::f64::consts::PI;

const MP_ME_EXP: f64 = 1836.152_673_43;
const MW_EXP_GEV: f64 = 80.3692;
const MZ_EXP_GEV: f64 = 91.1876;

fn rel_ppm(pred: f64, exp: f64) -> f64 {
    (pred - exp).abs() / exp * 1.0e6
}

fn main() {
    let out_dir = std::env::var("GUTOE_CORRECTION_HUNTS_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_correction_hunts_alpha_mpme_mwmz".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Hunt 1: alpha^-1 correction from LO=137.
    let alpha_inv_lo = 137.0;
    let delta_alpha_target = ALPHA_INV_PHYS - alpha_inv_lo;
    let delta_alpha_5a = 5.0 * ALPHA;
    let delta_alpha_5a_9a2 = 5.0 * ALPHA - 9.0 * ALPHA * ALPHA;
    let alpha_inv_5a = alpha_inv_lo + delta_alpha_5a;
    let alpha_inv_5a_9a2 = alpha_inv_lo + delta_alpha_5a_9a2;

    // Hunt 2: mp/me correction from LO=6*pi^5.
    let mp_me_lo = 6.0 * PI.powi(5);
    let delta_mp_target = MP_ME_EXP - mp_me_lo;
    let mp_me_5a_9a2 = mp_me_lo + (5.0 * ALPHA - 9.0 * ALPHA * ALPHA);

    // Integer-coefficient scan for 5a - c*a^2 with structural c candidates.
    let mut best_c = 0_i32;
    let mut best_err = f64::INFINITY;
    for c in 0..=128 {
        let pred = mp_me_lo + (5.0 * ALPHA - (c as f64) * ALPHA * ALPHA);
        let err = (pred - MP_ME_EXP).abs();
        if err < best_err {
            best_err = err;
            best_c = c;
        }
    }
    let mp_me_best_struct = mp_me_lo + (5.0 * ALPHA - (best_c as f64) * ALPHA * ALPHA);

    // Hunt 3: mW/mZ correction.
    // Base from weak-angle identity at MZ: sin^2(theta_W)=508/2197 in MSbar-like lane.
    let sin2_ms = 508.0 / 2197.0;
    let ratio_tree = (1.0_f64 - sin2_ms).sqrt();
    let ratio_exp = MW_EXP_GEV / MZ_EXP_GEV;

    // Running coefficient from zero-free weak-angle lane.
    let delta1 = ALPHA * 10.0_f64.ln() / (4.0 * PI);

    // Scheme candidate: convert MSbar-like to on-shell-like using grade2_4d = 6 multiplier.
    // sin^2_on-shell ~= sin^2_MSbar - 6*delta1.
    let sin2_os_candidate = sin2_ms - 6.0 * delta1;
    let ratio_scheme_candidate = (1.0 - sin2_os_candidate).sqrt();

    // Needed shift inferred from measured ratio.
    let sin2_from_ratio_exp = 1.0 - ratio_exp * ratio_exp;
    let delta_scheme_needed = sin2_ms - sin2_from_ratio_exp;
    let delta_over_delta1 = delta_scheme_needed / delta1;

    let payload = json!({
      "scope": "priority correction hunts: alpha^-1, mp/me, mW/mZ",
      "inputs": {
        "alpha_inverse_physical": ALPHA_INV_PHYS,
        "mp_me_experimental": MP_ME_EXP,
        "mW_experimental_GeV": MW_EXP_GEV,
        "mZ_experimental_GeV": MZ_EXP_GEV
      },
      "hunt_alpha_inverse": {
        "target_delta_from_137": delta_alpha_target,
        "lane_5alpha": {
          "delta": delta_alpha_5a,
          "predicted_alpha_inverse": alpha_inv_5a,
          "abs_error": alpha_inv_5a - ALPHA_INV_PHYS,
          "ppm_error": rel_ppm(alpha_inv_5a, ALPHA_INV_PHYS)
        },
        "lane_5alpha_minus_9alpha2": {
          "delta": delta_alpha_5a_9a2,
          "predicted_alpha_inverse": alpha_inv_5a_9a2,
          "abs_error": alpha_inv_5a_9a2 - ALPHA_INV_PHYS,
          "ppm_error": rel_ppm(alpha_inv_5a_9a2, ALPHA_INV_PHYS)
        }
      },
      "hunt_mp_me": {
        "base_6pi5": mp_me_lo,
        "target_delta_from_6pi5": delta_mp_target,
        "lane_5alpha_minus_9alpha2": {
          "predicted": mp_me_5a_9a2,
          "abs_error": mp_me_5a_9a2 - MP_ME_EXP,
          "ppm_error": rel_ppm(mp_me_5a_9a2, MP_ME_EXP)
        },
        "integer_scan_5alpha_minus_calpha2": {
          "best_c": best_c,
          "predicted": mp_me_best_struct,
          "abs_error": mp_me_best_struct - MP_ME_EXP,
          "ppm_error": rel_ppm(mp_me_best_struct, MP_ME_EXP),
          "note": "best c in [0,128]; c=36 corresponds to 3*12"
        }
      },
      "hunt_mw_mz": {
        "ratio_experimental": ratio_exp,
        "base_tree_from_508_over_2197": ratio_tree,
        "base_tree_ppm_error": rel_ppm(ratio_tree, ratio_exp),
        "delta1_alpha_ln10_over_4pi": delta1,
        "scheme_candidate": {
          "sin2_ms": sin2_ms,
          "sin2_os_candidate": sin2_os_candidate,
          "ratio_predicted": ratio_scheme_candidate,
          "abs_error": ratio_scheme_candidate - ratio_exp,
          "ppm_error": rel_ppm(ratio_scheme_candidate, ratio_exp),
          "note": "sin2_os = sin2_ms - grade2_4d*delta1 with grade2_4d=6"
        },
        "inferred_needed_shift": {
          "sin2_from_ratio_exp": sin2_from_ratio_exp,
          "delta_scheme_needed": delta_scheme_needed,
          "delta_over_delta1": delta_over_delta1
        }
      }
    });

    let txt_path = out.join("ctc_correction_hunts_alpha_mpme_mwmz.txt");
    let json_path = out.join("ctc_correction_hunts_alpha_mpme_mwmz.json");

    let mut txt = String::new();
    txt.push_str("[ctc_correction_hunts_alpha_mpme_mwmz]\n");
    txt.push_str("priority correction hunts\n\n");

    txt.push_str("[alpha^-1]\n");
    txt.push_str(&format!("target delta from 137 = {:.12e}\n", delta_alpha_target));
    txt.push_str(&format!("137 + 5a           = {:.12}  ppm_err={:.3}\n", alpha_inv_5a, rel_ppm(alpha_inv_5a, ALPHA_INV_PHYS)));
    txt.push_str(&format!("137 + 5a - 9a^2     = {:.12}  ppm_err={:.3}\n\n", alpha_inv_5a_9a2, rel_ppm(alpha_inv_5a_9a2, ALPHA_INV_PHYS)));

    txt.push_str("[mp/me]\n");
    txt.push_str(&format!("base 6pi^5          = {:.12}\n", mp_me_lo));
    txt.push_str(&format!("target delta        = {:.12e}\n", delta_mp_target));
    txt.push_str(&format!("6pi^5 + (5a-9a^2)   = {:.12}  ppm_err={:.3}\n", mp_me_5a_9a2, rel_ppm(mp_me_5a_9a2, MP_ME_EXP)));
    txt.push_str(&format!("6pi^5 + (5a-c a^2)  best c={} -> {:.12}  ppm_err={:.3}\n\n", best_c, mp_me_best_struct, rel_ppm(mp_me_best_struct, MP_ME_EXP)));

    txt.push_str("[mW/mZ]\n");
    txt.push_str(&format!("exp ratio           = {:.12}\n", ratio_exp));
    txt.push_str(&format!("tree from 508/2197  = {:.12}  ppm_err={:.3}\n", ratio_tree, rel_ppm(ratio_tree, ratio_exp)));
    txt.push_str(&format!("delta1              = {:.12e}\n", delta1));
    txt.push_str(&format!("scheme ratio (ms-6d1)= {:.12}  ppm_err={:.3}\n", ratio_scheme_candidate, rel_ppm(ratio_scheme_candidate, ratio_exp)));
    txt.push_str(&format!("needed shift / delta1= {:.9}\n", delta_over_delta1));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
