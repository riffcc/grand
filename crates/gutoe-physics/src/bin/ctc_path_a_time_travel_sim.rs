//! Path-A time-travel simulation with branch-shift quanta.
//!
//! Worldline model:
//! - Local timelike travel start -> node and node -> goal at speed beta*c.
//! - n local loops at node, each with proper/cover time T.
//! - Effective coordinate branch shift: q_eff*T per loop.
//!
//! Effective arrival:
//!   t_eff = dt_in + dt_out + n*(1-q_eff)*T
//!
//! Budget-closed lane:
//!   q_eff = budget_j / threshold_j
//!   threshold_j = kappa * (3/16) * |R| * |T|
//!
//! q_eff=1: break-even.
//! q_eff>1: possible pre-departure coordinate arrival for large enough n.
//! q_eff<1: no time gain.

use serde_json::json;
use std::f64::consts::PI;
use std::fs;
use std::path::PathBuf;

use gutoe_physics::constants::C;

const FC_VOID: f64 = 3.0 / 16.0;
const REAR_FACE_FACTOR: f64 = 1.0 / 10.0;
const V_EWSB_GEV: f64 = 245.3;
const GEV_TO_J: f64 = 1.602_176_634e-10;
const HBARC_GEV_M: f64 = 0.197_326_980_4e-15;

fn abs(x: f64) -> f64 {
    x.abs()
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

fn threshold_units(radius_m: f64, period_s: f64) -> f64 {
    FC_VOID * abs(radius_m) * abs(period_s)
}

fn segment_dt(dx_m: f64, beta: f64) -> f64 {
    abs(dx_m) / (beta * C)
}

fn segment_proper(dt_s: f64, beta: f64) -> f64 {
    dt_s * (1.0 - beta * beta).sqrt()
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_PATH_A_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_path_a_time_travel_sim".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let x_start = std::env::var("GUTOE_CTC_PATH_A_START_X_M")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(-1_000.0);
    let x_goal = std::env::var("GUTOE_CTC_PATH_A_GOAL_X_M")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1_000.0);
    let x_node = std::env::var("GUTOE_CTC_PATH_A_NODE_X_M")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0);
    let beta = std::env::var("GUTOE_CTC_PATH_A_BETA")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.8)
        .clamp(1e-6, 0.999_999);
    let t_loop = std::env::var("GUTOE_CTC_PATH_A_LOOP_PERIOD_S")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(50.0)
        .max(1e-30);
    let n_loops_max = std::env::var("GUTOE_CTC_PATH_A_N_MAX")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(120);

    // Structural creation geometry used to close q from budget.
    let radius_m = std::env::var("GUTOE_CTC_PATH_A_RADIUS_M")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(20.0)
        .max(0.0);
    let kappa_j_per_m_s = std::env::var("GUTOE_CTC_PATH_A_KAPPA_J_PER_M_S")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or_else(derived_kappa_j_per_m_s)
        .max(0.0);
    let threshold_j = kappa_j_per_m_s * threshold_units(radius_m, t_loop);
    let budget_j = std::env::var("GUTOE_CTC_PATH_A_BUDGET_J")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(threshold_j)
        .max(0.0);
    let q_eff = if threshold_j > 0.0 {
        budget_j / threshold_j
    } else {
        1.0
    };

    // Optional integer sweep retained as a diagnostic lane.
    let q_max = std::env::var("GUTOE_CTC_PATH_A_Q_MAX")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(4)
        .max(1);

    let dt_in = segment_dt(x_node - x_start, beta);
    let dt_out = segment_dt(x_goal - x_node, beta);
    let d_tau_in = segment_proper(dt_in, beta);
    let d_tau_out = segment_proper(dt_out, beta);

    let baseline_light = abs(x_goal - x_start) / C;
    let baseline_cover = dt_in + dt_out;
    let baseline_proper = d_tau_in + d_tau_out;

    let mut derived_first_predeparture: Option<usize> = None;
    let mut derived_best_eff = f64::INFINITY;
    let mut derived_best_n = 0usize;
    let mut derived_best_cover = 0.0_f64;
    let mut derived_best_proper = 0.0_f64;

    for n in 0..=n_loops_max {
        let nf = n as f64;
        let cover = dt_in + nf * t_loop + dt_out;
        let proper = d_tau_in + nf * t_loop + d_tau_out;
        let eff = dt_in + dt_out + nf * (1.0 - q_eff) * t_loop;

        if eff < derived_best_eff {
            derived_best_eff = eff;
            derived_best_n = n;
            derived_best_cover = cover;
            derived_best_proper = proper;
        }
        if derived_first_predeparture.is_none() && eff < 0.0 {
            derived_first_predeparture = Some(n);
        }
    }
    let derived_effective_superluminal = derived_best_eff < baseline_light;
    let derived_pre_departure = derived_best_eff < 0.0;
    let derived_n_required_predeparture = if q_eff > 1.0 {
        let denom = (q_eff - 1.0) * t_loop;
        Some(((dt_in + dt_out) / denom).ceil())
    } else {
        None
    };

    let mut rows = Vec::new();

    for q in 1..=q_max {
        let qf = q as f64;
        let mut first_predeparture: Option<usize> = None;
        let mut best_eff = f64::INFINITY;
        let mut best_n = 0usize;
        let mut best_cover = 0.0_f64;
        let mut best_proper = 0.0_f64;

        for n in 0..=n_loops_max {
            let nf = n as f64;
            let cover = dt_in + nf * t_loop + dt_out;
            let proper = d_tau_in + nf * t_loop + d_tau_out;
            let eff = dt_in + dt_out + nf * (1.0 - qf) * t_loop;

            if eff < best_eff {
                best_eff = eff;
                best_n = n;
                best_cover = cover;
                best_proper = proper;
            }
            if first_predeparture.is_none() && eff < 0.0 {
                first_predeparture = Some(n);
            }
        }

        let effective_superluminal = best_eff < baseline_light;
        let pre_departure = best_eff < 0.0;
        let n_required_predeparture = if qf > 1.0 {
            let denom = (qf - 1.0) * t_loop;
            Some(((dt_in + dt_out) / denom).ceil())
        } else {
            None
        };

        rows.push(json!({
            "q_shift_quanta_per_loop": q,
            "n_required_predeparture": n_required_predeparture,
            "first_predeparture_n": first_predeparture,
            "best_n": best_n,
            "best_effective_arrival_s": best_eff,
            "best_cover_arrival_s": best_cover,
            "best_local_proper_time_s": best_proper,
            "effective_superluminal": effective_superluminal,
            "pre_departure": pre_departure
        }));
    }

    let payload = json!({
        "inputs": {
            "x_start_m": x_start,
            "x_goal_m": x_goal,
            "x_node_m": x_node,
            "beta": beta,
            "radius_m": radius_m,
            "loop_period_s": t_loop,
            "n_loops_max": n_loops_max,
            "q_max": q_max,
            "budget_j": budget_j,
            "kappa_j_per_m_s": kappa_j_per_m_s
        },
        "baselines": {
            "baseline_light_time_s": baseline_light,
            "baseline_cover_time_s": baseline_cover,
            "baseline_local_proper_time_s": baseline_proper,
            "dt_in_s": dt_in,
            "dt_out_s": dt_out,
            "d_tau_in_s": d_tau_in,
            "d_tau_out_s": d_tau_out
        },
        "derived_q_from_budget": {
            "formula_q_eff": "q_eff = budget_j / (kappa*(3/16)*|R|*|T|)",
            "threshold_j": threshold_j,
            "q_eff": q_eff,
            "n_required_predeparture": derived_n_required_predeparture,
            "first_predeparture_n": derived_first_predeparture,
            "best_n": derived_best_n,
            "best_effective_arrival_s": derived_best_eff,
            "best_cover_arrival_s": derived_best_cover,
            "best_local_proper_time_s": derived_best_proper,
            "effective_superluminal": derived_effective_superluminal,
            "pre_departure": derived_pre_departure
        },
        "q_sweep": rows,
        "timelike_local_segments": true,
        "formula": "t_eff = dt_in + dt_out + n*(1-q_eff)*T"
    });

    let txt_path = out.join("ctc_path_a_time_travel_sim.txt");
    let json_path = out.join("ctc_path_a_time_travel_sim.json");

    let mut txt = String::new();
    txt.push_str("[ctc_path_a_time_travel_sim]\n");
    txt.push_str("formula = dt_in + dt_out + n*(1-q_eff)*T\n");
    txt.push_str(&format!(
        "baseline_light={:.12e}s, baseline_cover={:.12e}s, baseline_proper={:.12e}s\n",
        baseline_light, baseline_cover, baseline_proper
    ));
    txt.push_str("\n[derived_q_from_budget]\n");
    txt.push_str(&format!(
        "radius={:.6e}m, T={:.6e}s, kappa={:.6e}J/(m*s), budget={:.6e}J, threshold={:.6e}J, q_eff={:.12e}\n",
        radius_m, t_loop, kappa_j_per_m_s, budget_j, threshold_j, q_eff
    ));
    txt.push_str(&format!(
        "n_required_pre={:?}, first_pre_n={:?}, best_n={}, best_eff={:.12e}s, best_cover={:.12e}s, best_proper={:.12e}s, superluminal={}, pre_departure={}\n",
        derived_n_required_predeparture,
        derived_first_predeparture,
        derived_best_n,
        derived_best_eff,
        derived_best_cover,
        derived_best_proper,
        derived_effective_superluminal,
        derived_pre_departure
    ));
    txt.push_str("\n[q_sweep]\n");
    for r in payload["q_sweep"].as_array().expect("array") {
        txt.push_str(&format!(
            "q={}: n_required_pre={:?}, first_pre_n={:?}, best_n={}, best_eff={:.12e}s, best_cover={:.12e}s, best_proper={:.12e}s, superluminal={}, pre_departure={}\n",
            r["q_shift_quanta_per_loop"].as_u64().unwrap_or(0),
            r["n_required_predeparture"].as_f64(),
            r["first_predeparture_n"].as_u64(),
            r["best_n"].as_u64().unwrap_or(0),
            r["best_effective_arrival_s"].as_f64().unwrap_or(f64::NAN),
            r["best_cover_arrival_s"].as_f64().unwrap_or(f64::NAN),
            r["best_local_proper_time_s"].as_f64().unwrap_or(f64::NAN),
            r["effective_superluminal"].as_bool().unwrap_or(false),
            r["pre_departure"].as_bool().unwrap_or(false)
        ));
    }

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
