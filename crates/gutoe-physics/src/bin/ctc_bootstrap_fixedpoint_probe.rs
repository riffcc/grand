//! CTC bootstrap fixed-point probe.
//!
//! Model:
//!   E_past = eta * E_future - loss
//! Fixed point:
//!   E* = eta * E* - loss
//!
//! We test if a threshold-targeted bootstrap point exists:
//!   E* >= E_threshold
//! for both closed-cycle assumptions and open-flux extensions.

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

fn closed_fixedpoint(eta: f64, loss_j: f64) -> Option<f64> {
    // E = eta*E - loss
    let den = 1.0 - eta;
    if den.abs() < 1e-18 {
        if loss_j.abs() < 1e-18 {
            // Degenerate continuum of fixed points (ideal circulation)
            None
        } else {
            // No finite fixed point
            None
        }
    } else {
        Some(-loss_j / den)
    }
}

fn open_fixedpoint(eta: f64, inflow_j: f64, loss_j: f64) -> Option<f64> {
    // E = eta*E + inflow - loss
    let den = 1.0 - eta;
    if den.abs() < 1e-18 {
        None
    } else {
        Some((inflow_j - loss_j) / den)
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_BOOTSTRAP_FIXEDPOINT_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_bootstrap_fixedpoint_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let radius_m = env_f64("GUTOE_CTC_BOOTSTRAP_RADIUS_M", 8.044_312_286_995_516e-19).max(0.0);
    let period_s = env_f64("GUTOE_CTC_BOOTSTRAP_PERIOD_S", 5.366_854_127e-27).max(1e-30);
    let eta = env_f64("GUTOE_CTC_BOOTSTRAP_ETA", 0.999_999_999_999);
    let loss_j = env_f64("GUTOE_CTC_BOOTSTRAP_LOSS_J", 1e-12).max(0.0);
    let inflow_j = env_f64("GUTOE_CTC_BOOTSTRAP_INFLOW_J", 0.0).max(0.0);
    let kappa = env_f64("GUTOE_CTC_BOOTSTRAP_KAPPA_J_PER_M_S", derived_kappa_j_per_m_s()).max(0.0);

    let thr = threshold_j(radius_m, period_s, kappa);

    let closed_fp = closed_fixedpoint(eta, loss_j);
    let closed_ideal_continuum = (1.0 - eta).abs() < 1e-18 && loss_j.abs() < 1e-18;
    let closed_finite_feasible = closed_fp.is_some_and(|e| e.is_finite() && e >= thr && e >= 0.0);
    let closed_feasible = closed_ideal_continuum || closed_finite_feasible;
    let closed_reason = if closed_ideal_continuum {
        "ideal lossless continuum (eta=1, loss=0)"
    } else if (1.0 - eta).abs() < 1e-18 && loss_j > 0.0 {
        "no fixed point (eta=1 with positive loss)"
    } else if eta <= 1.0 && loss_j >= 0.0 {
        "under eta<=1 and loss>=0, finite positive closed fixed point cannot exceed threshold unless loss=0 and eta=1"
    } else {
        "eta>1 amplification regime (not closed-cycle conservative)"
    };

    let open_fp = open_fixedpoint(eta, inflow_j, loss_j);
    let open_feasible = open_fp.is_some_and(|e| e.is_finite() && e >= thr && e >= 0.0);

    let payload = json!({
        "inputs": {
            "radius_m": radius_m,
            "period_s": period_s,
            "eta": eta,
            "loss_j": loss_j,
            "inflow_j": inflow_j,
            "kappa_j_per_m_s": kappa
        },
        "threshold_j": thr,
        "closed_cycle": {
            "fixed_point_j": closed_fp,
            "ideal_continuum": closed_ideal_continuum,
            "finite_feasible_at_threshold": closed_finite_feasible,
            "feasible_at_threshold": closed_feasible,
            "reason": closed_reason
        },
        "open_flux": {
            "fixed_point_j": open_fp,
            "feasible_at_threshold": open_feasible
        }
    });

    let txt_path = out.join("ctc_bootstrap_fixedpoint_probe.txt");
    let json_path = out.join("ctc_bootstrap_fixedpoint_probe.json");

    let mut txt = String::new();
    txt.push_str("[ctc_bootstrap_fixedpoint_probe]\n");
    txt.push_str("closed map: E_past = eta*E_future - loss\n");
    txt.push_str("open map:   E_past = eta*E_future + inflow - loss\n");
    txt.push_str(&format!(
        "threshold_j={:.12e}, eta={:.12e}, loss_j={:.12e}, inflow_j={:.12e}\n",
        thr, eta, loss_j, inflow_j
    ));
    txt.push_str("\n[closed_cycle]\n");
    txt.push_str(&format!(
        "fixed_point_j={:?}, ideal_continuum={}, finite_feasible_at_threshold={}, feasible_at_threshold={}, reason={}\n",
        closed_fp, closed_ideal_continuum, closed_finite_feasible, closed_feasible, closed_reason
    ));
    txt.push_str("\n[open_flux]\n");
    txt.push_str(&format!(
        "fixed_point_j={:?}, feasible_at_threshold={}\n",
        open_fp, open_feasible
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
