//! Origin-energy closure probe for multi-sheet fan-in.
//!
//! Goal: quantify what effective fan-in gain is required to match a target
//! origin energy from a seed contribution under finite causal depth.

use gutoe_physics::constants::C;
use serde_json::json;
use std::f64::consts::PI;
use std::fs;
use std::path::PathBuf;

const FC_VOID: f64 = 3.0 / 16.0;
const REAR_FACE_FACTOR: f64 = 1.0 / 10.0;
const V_EWSB_GEV: f64 = 245.3;
const GEV_TO_J: f64 = 1.602_176_634e-10;
const HBARC_GEV_M: f64 = 0.197_326_980_4e-15;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_string(name: &str, default: &str) -> String {
    std::env::var(name).unwrap_or_else(|_| default.to_string())
}

fn wall_tension_front_j_m2(delta_theta: f64, thickness_m: f64) -> f64 {
    let l_nat = thickness_m / HBARC_GEV_M;
    let sigma_gev3 = FC_VOID * V_EWSB_GEV.powi(2) * delta_theta.powi(2) / (2.0 * l_nat);
    let gev3_to_j_m2 = GEV_TO_J / HBARC_GEV_M.powi(2);
    sigma_gev3 * gev3_to_j_m2
}

fn derived_kappa_j_per_m_s() -> f64 {
    let thickness_m = HBARC_GEV_M / V_EWSB_GEV;
    let sigma_front = wall_tension_front_j_m2(PI, thickness_m);
    let sigma_rear = REAR_FACE_FACTOR * sigma_front;
    2.0 * PI * C * sigma_rear / FC_VOID
}

fn threshold_j(radius_m: f64, period_s: f64, kappa: f64) -> f64 {
    kappa * FC_VOID * radius_m.abs() * period_s.abs()
}

fn log_sum_geo(b: f64, k: f64) -> f64 {
    if b <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if (b - 1.0).abs() < 1e-14 {
        return (k + 1.0).ln();
    }
    if b > 1.0 {
        let a = (k + 1.0) * b.ln();
        if a > 50.0 {
            // (b^(k+1)-1)/(b-1) ~ b^(k+1)/(b-1)
            return a - (b - 1.0).ln();
        }
        let num = b.powf(k + 1.0) - 1.0;
        return (num / (b - 1.0)).ln();
    }
    // 0 < b < 1
    let num = 1.0 - b.powf(k + 1.0);
    (num / (1.0 - b)).ln()
}

fn reaches_target(seed: f64, b: f64, k: f64, target: f64) -> bool {
    if seed <= 0.0 || target <= 0.0 || b <= 0.0 {
        return false;
    }
    let lhs = seed.ln() + log_sum_geo(b, k);
    lhs >= target.ln()
}

fn find_bmin(seed: f64, k: f64, target: f64) -> Option<f64> {
    if !reaches_target(seed, 2.0, k, target) {
        return None;
    }
    let mut lo = 1e-18_f64;
    let mut hi = 2.0_f64;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if reaches_target(seed, mid, k, target) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(hi)
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_ORIGIN_CLOSURE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_origin_energy_closure_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Quark-scale defaults from prior lane.
    let radius_m = env_f64("GUTOE_CTC_ORIGIN_RADIUS_M", 8.044_312_286_995_516e-19).max(0.0);
    let period_s = env_f64("GUTOE_CTC_ORIGIN_PERIOD_S", 5.366_854_127e-27).max(1e-30);
    let target_j = env_f64("GUTOE_CTC_ORIGIN_TARGET_J", 1e69).max(1e-300);
    let horizon_s = env_f64("GUTOE_CTC_ORIGIN_HORIZON_S", 4.354_0e17).max(period_s);

    // Effective gain decomposition knobs.
    // Mode `even_subalgebra` hard-wires the structural split:
    // G_eff = 3 * (3/16) * (2/3) * (8/3) = 1.
    let mode_raw = env_string("GUTOE_CTC_ORIGIN_MODE", "legacy");
    let mode = mode_raw.to_ascii_lowercase();
    let (branching, merge_fraction, eta, infra_gain, mode_note) = match mode.as_str() {
        "even_subalgebra" | "z3_void_grade_split" => (
            3.0,
            FC_VOID,
            2.0 / 3.0,
            8.0 / 3.0,
            "structural override: 3*(3/16)*(2/3)*(8/3)",
        ),
        "canonical_even_subalgebra" => (
            2.0,
            1.0,
            0.5,
            1.0,
            "canonical-even override: 2*1*(1/2)*1",
        ),
        _ => (
            env_f64("GUTOE_CTC_ORIGIN_BRANCHING", 2.0).max(0.0),
            env_f64("GUTOE_CTC_ORIGIN_MERGE_FRACTION", FC_VOID).clamp(0.0, 1.0),
            env_f64("GUTOE_CTC_ORIGIN_ETA", 0.999_999_999_999).max(0.0),
            env_f64("GUTOE_CTC_ORIGIN_INFRA_GAIN", 1.0).max(0.0),
            "legacy/env knobs",
        ),
    };

    let kappa = env_f64("GUTOE_CTC_ORIGIN_KAPPA_J_PER_M_S", derived_kappa_j_per_m_s()).max(0.0);
    let seed_j = threshold_j(radius_m, period_s, kappa);

    let kmax = (horizon_s / period_s).floor().max(1.0);
    let ratio = target_j / seed_j;

    let b_eff = branching * merge_fraction * eta * infra_gain;
    let finite_reaches = reaches_target(seed_j, b_eff, kmax, target_j);

    let bmin = find_bmin(seed_j, kmax, target_j);
    // Infinite closed geometric lane needs B = 1 - epsilon with epsilon = seed/target.
    let eps_closed = seed_j / target_j;
    let b_closed_infinite = 1.0 - eps_closed; // numerically rounds to 1.0 in f64 for tiny eps.

    let payload = json!({
        "inputs": {
            "mode": mode,
            "mode_note": mode_note,
            "radius_m": radius_m,
            "period_s": period_s,
            "horizon_s": horizon_s,
            "target_j": target_j,
            "branching": branching,
            "merge_fraction": merge_fraction,
            "eta": eta,
            "infra_gain": infra_gain
        },
        "derived": {
            "kappa_j_per_m_s": kappa,
            "seed_j": seed_j,
            "target_over_seed": ratio,
            "kmax": kmax,
            "b_eff": b_eff,
            "finite_horizon_reaches_target": finite_reaches,
            "bmin_for_finite_horizon": bmin,
            "b_closed_infinite": b_closed_infinite,
            "epsilon_closed_infinite": eps_closed,
            "b_closed_infinite_precision_note": "for tiny epsilon, f64 rounds 1-epsilon to 1.0"
        }
    });

    let txt_path = out.join("ctc_origin_energy_closure_probe.txt");
    let json_path = out.join("ctc_origin_energy_closure_probe.json");

    let mut txt = String::new();
    txt.push_str("[ctc_origin_energy_closure_probe]\n");
    txt.push_str(&format!(
        "mode={} ({})\n",
        mode_raw, mode_note
    ));
    txt.push_str(&format!(
        "seed_j={:.12e}, target_j={:.12e}, ratio={:.12e}\n",
        seed_j, target_j, ratio
    ));
    txt.push_str(&format!(
        "period_s={:.12e}, horizon_s={:.12e}, kmax={:.12e}\n",
        period_s, horizon_s, kmax
    ));
    txt.push_str(&format!(
        "b_eff={:.12e} (branching*merge*eta*infra)\n",
        b_eff
    ));
    txt.push_str(&format!(
        "finite_horizon_reaches_target={}\n",
        finite_reaches
    ));
    txt.push_str(&format!(
        "bmin_for_finite_horizon={:?}\n",
        bmin
    ));
    txt.push_str(&format!(
        "b_closed_infinite={:.18e}, epsilon_closed_infinite={:.18e}\n",
        b_closed_infinite, eps_closed
    ));
    txt.push_str("note: b_closed_infinite may round to 1.0 in f64 when epsilon is extremely small\n");
    txt.push_str("\n[interpretation]\n");
    txt.push_str("finite-horizon model: target is reachable iff b_eff >= bmin_for_finite_horizon\n");
    txt.push_str("infinite closed geometric model (B<1): requires epsilon = seed/target\n");

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
