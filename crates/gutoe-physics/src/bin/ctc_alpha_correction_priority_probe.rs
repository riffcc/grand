//! Priority probe: structural interpretation of alpha correction lane.
//!
//! Targets from current correction hunts:
//! - alpha^-1 ~ 137 + 5a - 9a^2
//! - mp/me    ~ 6pi^5 + 5a - 36a^2
//!
//! Goal: test structural candidate decompositions for coefficients 5, 9, 36
//! from small Cl(1,3)/Z3 counts.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

const ALPHA_INV_PHYS: f64 = 137.035_999_177;
const ALPHA: f64 = 1.0 / ALPHA_INV_PHYS;
const MP_ME_EXP: f64 = 1836.152_673_43;
const PI: f64 = std::f64::consts::PI;

const LO_ALPHA_INV: f64 = 137.0;

// Structural integers from current framework lanes.
const Z3_ORDER: i32 = 3;
const GRADE_LEVELS: i32 = 5; // grades 0..4
const GRADE1_DIM: i32 = 4;
const GRADE2_DIM: i32 = 6;
const FIBER_DIM: i32 = 12;
const TOTAL_BASIS: i32 = 16;
const NON_VOID_DIM: i32 = 13;

fn ppm(pred: f64, truth: f64) -> f64 {
    (pred - truth).abs() / truth * 1e6
}

fn lane_alpha(a1: f64, a2: f64) -> f64 {
    LO_ALPHA_INV + a1 * ALPHA - a2 * ALPHA * ALPHA
}

fn lane_mp(lo_mp_me: f64, a1: f64, a2: f64) -> f64 {
    lo_mp_me + a1 * ALPHA - a2 * ALPHA * ALPHA
}

fn main() {
    let out_dir = std::env::var("GUTOE_ALPHA_PRIORITY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_alpha_correction_priority_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);
    let lo_mp_me = 6.0 * PI.powi(5);

    // Implied coefficients if we fix linear term a1 = 5.
    let a1 = GRADE_LEVELS as f64;
    let b_alpha_exact = (LO_ALPHA_INV + a1 * ALPHA - ALPHA_INV_PHYS) / (ALPHA * ALPHA);
    let b_mp_exact = (lo_mp_me + a1 * ALPHA - MP_ME_EXP) / (ALPHA * ALPHA);

    // Structural candidate families.
    let b9_z3_sq = (Z3_ORDER * Z3_ORDER) as f64; // 9
    let b9_grade_mix = (GRADE2_DIM + Z3_ORDER) as f64; // 9
    let b36_g1_z3_sq = (GRADE1_DIM * Z3_ORDER * Z3_ORDER) as f64; // 36
    let b36_g2_sq = (GRADE2_DIM * GRADE2_DIM) as f64; // 36

    let alpha_pred_9 = lane_alpha(a1, b9_z3_sq);
    let alpha_pred_9_alt = lane_alpha(a1, b9_grade_mix);

    let mp_pred_36 = lane_mp(lo_mp_me, a1, b36_g1_z3_sq);
    let mp_pred_36_alt = lane_mp(lo_mp_me, a1, b36_g2_sq);

    let payload = json!({
      "scope": "alpha correction priority structural probe",
      "constants": {
        "alpha_inverse_physical": ALPHA_INV_PHYS,
        "alpha": ALPHA,
        "mp_me_physical": MP_ME_EXP,
        "framework_counts": {
          "z3_order": Z3_ORDER,
          "grade_levels": GRADE_LEVELS,
          "grade1_dim": GRADE1_DIM,
          "grade2_dim": GRADE2_DIM,
          "fiber_dim": FIBER_DIM,
          "total_basis": TOTAL_BASIS,
          "non_void_dim": NON_VOID_DIM
        }
      },
      "implied_coefficients_with_a1_5": {
        "b_alpha_exact": b_alpha_exact,
        "b_mp_exact": b_mp_exact,
        "b_mp_over_b_alpha": b_mp_exact / b_alpha_exact
      },
      "candidate_coefficients": {
        "b9": {
          "z3_sq": b9_z3_sq,
          "grade2_plus_z3": b9_grade_mix
        },
        "b36": {
          "grade1_times_z3_sq": b36_g1_z3_sq,
          "grade2_sq": b36_g2_sq
        }
      },
      "predictions": {
        "alpha_lane": {
          "formula": "137 + 5a - b a^2",
          "b=9_z3_sq": {
            "pred": alpha_pred_9,
            "ppm_error": ppm(alpha_pred_9, ALPHA_INV_PHYS)
          },
          "b=9_grade2_plus_z3": {
            "pred": alpha_pred_9_alt,
            "ppm_error": ppm(alpha_pred_9_alt, ALPHA_INV_PHYS)
          }
        },
        "mp_me_lane": {
          "formula": "6pi^5 + 5a - B a^2",
          "B=36_grade1_times_z3_sq": {
            "pred": mp_pred_36,
            "ppm_error": ppm(mp_pred_36, MP_ME_EXP)
          },
          "B=36_grade2_sq": {
            "pred": mp_pred_36_alt,
            "ppm_error": ppm(mp_pred_36_alt, MP_ME_EXP)
          }
        }
      }
    });

    let txt_path = out.join("ctc_alpha_correction_priority_probe.txt");
    let json_path = out.join("ctc_alpha_correction_priority_probe.json");

    let mut txt = String::new();
    txt.push_str("[ctc_alpha_correction_priority_probe]\n");
    txt.push_str("priority structural probe for alpha correction lane\n\n");

    txt.push_str("[implied coefficients with a1=5]\n");
    txt.push_str(&format!("b_alpha_exact = {:.9}\n", b_alpha_exact));
    txt.push_str(&format!("b_mp_exact    = {:.9}\n", b_mp_exact));
    txt.push_str(&format!("ratio b_mp/b_alpha = {:.9}\n\n", b_mp_exact / b_alpha_exact));

    txt.push_str("[candidate decompositions]\n");
    txt.push_str("9  = z3^2 = 3^2\n");
    txt.push_str("9  = grade2 + z3 = 6 + 3\n");
    txt.push_str("36 = grade1 * z3^2 = 4*9\n");
    txt.push_str("36 = grade2^2 = 6^2\n\n");

    txt.push_str("[alpha lane]\n");
    txt.push_str(&format!("137 + 5a - 9a^2 = {:.12}  ppm_err={:.3}\n", alpha_pred_9, ppm(alpha_pred_9, ALPHA_INV_PHYS)));

    txt.push_str("\n[mp/me lane]\n");
    txt.push_str(&format!("6pi^5 + 5a - 36a^2 = {:.12}  ppm_err={:.3}\n", mp_pred_36, ppm(mp_pred_36, MP_ME_EXP)));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
