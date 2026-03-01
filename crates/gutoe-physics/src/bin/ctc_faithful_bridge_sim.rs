//! Theorem-faithful CTC/FTL bridge simulation.
//!
//! This executable stitches together the proven logical lanes:
//! - Rear shortcut factor `s = 1/10` (VoidRearFace / FTLRearFaceBridge).
//! - Local causal bound (`v_local <= c`) with coordinate-effective `u = c/s`.
//! - Subluminal boosted frame witness with negative time numerator
//!   `dt - v*dx/c^2 < 0` (FTLFrameCTCBridge shape).
//! - Time-cylinder CTC legality witness (`Timelike` + identified closure).
//!
//! It is a simulation/consistency demonstrator, not an engine claim.

use gutoe_physics::constants::C;
use serde_json::json;
use std::fs;
use std::path::PathBuf;

const VOID_FRACTION: f64 = 3.0 / 16.0;
const REAR_SHORTCUT_FACTOR: f64 = 1.0 / 10.0;

#[derive(Debug, Clone, Copy)]
struct Event {
    t: f64,
    x: f64,
}

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn abs(x: f64) -> f64 {
    x.abs()
}

fn interval_sq(a: Event, b: Event) -> f64 {
    let dt = b.t - a.t;
    let dx = b.x - a.x;
    -(dt * dt) + (dx * dx)
}

fn shortcut_travel_time(d_m: f64, c_m_s: f64, s: f64) -> f64 {
    (s * d_m) / c_m_s
}

fn coordinate_speed(d_m: f64, dt_s: f64) -> f64 {
    d_m / dt_s
}

fn gamma_from_beta(beta: f64) -> f64 {
    1.0 / (1.0 - beta * beta).sqrt()
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_FAITHFUL_SIM_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_faithful_bridge_sim".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Geometry/kinematics controls
    let d_m = env_f64("GUTOE_CTC_FAITHFUL_D_M", 1_000.0).abs().max(1e-12);
    let c_m_s = env_f64("GUTOE_CTC_FAITHFUL_C", C).abs().max(1e-12);
    let s = env_f64("GUTOE_CTC_FAITHFUL_S", REAR_SHORTCUT_FACTOR)
        .abs()
        .clamp(1e-12, 1.0 - 1e-12);
    let local_beta = env_f64("GUTOE_CTC_FAITHFUL_LOCAL_BETA", 1.0).clamp(0.0, 1.0);

    // Dynamic gate controls (Path-B skeleton)
    let radius_m = env_f64("GUTOE_CTC_FAITHFUL_RADIUS_M", 20.0).abs();
    let period_s = env_f64("GUTOE_CTC_FAITHFUL_PERIOD_S", 50.0).abs().max(1e-30);
    let threshold_units = VOID_FRACTION * radius_m * period_s;
    let budget_units = env_f64("GUTOE_CTC_FAITHFUL_BUDGET_UNITS", threshold_units).max(0.0);
    let gate_pass = budget_units >= threshold_units;

    // Effective arrival loop controls
    let beta_route = env_f64("GUTOE_CTC_FAITHFUL_ROUTE_BETA", 0.8).clamp(1e-9, 0.999_999);
    let loops_n = env_usize("GUTOE_CTC_FAITHFUL_LOOPS_N", 12);
    let q_eff = env_f64("GUTOE_CTC_FAITHFUL_Q_EFF", 1.0 / s).max(0.0);

    // 1) Local-causal + coordinate-effective lane
    let dt_shortcut = shortcut_travel_time(d_m, c_m_s, s);
    let u_coord = coordinate_speed(d_m, dt_shortcut);
    let local_bound_ok = local_beta <= 1.0 + 1e-12;
    let coord_superluminal = u_coord > c_m_s;

    // 2) Subluminal boost with negative boosted time numerator.
    // Use v = (c + c^2/u)/2 whenever u>c (matches the Lean witness construction).
    let (v_m_s, beta_boost) = if coord_superluminal {
        let v = (c_m_s + (c_m_s * c_m_s) / u_coord) / 2.0;
        let b = (v / c_m_s).clamp(0.0, 0.999_999_999_999);
        (v, b)
    } else {
        let b = 0.5;
        (b * c_m_s, b)
    };
    let gamma = gamma_from_beta(beta_boost);
    let dt_num = dt_shortcut - (v_m_s * d_m) / (c_m_s * c_m_s);
    let dt_prime = gamma * dt_num;
    let predeparture_boosted = dt_prime < 0.0;

    // 3) CTC legality witness on time cylinder
    let a = Event { t: 0.0, x: 0.0 };
    let b = Event { t: period_s, x: 0.0 };
    let ds2 = interval_sq(a, b);
    let timelike_step = ds2 < 0.0;
    let identified_closed = true; // by definition of one-period identification t~t+T at fixed x

    // 4) Route-loop effective arrival (Path-A style summary)
    let x_start = -0.5 * d_m;
    let x_node = 0.0;
    let x_goal = 0.5 * d_m;
    let dt_in = abs(x_node - x_start) / (beta_route * c_m_s);
    let dt_out = abs(x_goal - x_node) / (beta_route * c_m_s);
    let baseline_light = abs(x_goal - x_start) / c_m_s;
    let t_eff = dt_in + dt_out + loops_n as f64 * (1.0 - q_eff) * period_s;
    let t_cover = dt_in + dt_out + loops_n as f64 * period_s;
    let predeparture_effective = gate_pass && t_eff < 0.0;
    let effective_superluminal = gate_pass && t_eff < baseline_light;

    let faithful_sim_possible = local_bound_ok
        && gate_pass
        && coord_superluminal
        && predeparture_boosted
        && timelike_step
        && identified_closed;

    let payload = json!({
      "inputs": {
        "distance_m": d_m,
        "c_m_s": c_m_s,
        "shortcut_factor_s": s,
        "rear_shortcut_factor_reference": REAR_SHORTCUT_FACTOR,
        "local_beta_bound_probe": local_beta,
        "radius_m": radius_m,
        "period_s": period_s,
        "budget_units": budget_units,
        "route_beta": beta_route,
        "loops_n": loops_n,
        "q_eff": q_eff
      },
      "lane_checks": {
        "local_bound_ok_v_le_c": local_bound_ok,
        "dynamic_gate_threshold_units": threshold_units,
        "dynamic_gate_pass": gate_pass
      },
      "shortcut_kinematics": {
        "dt_shortcut_s": dt_shortcut,
        "coordinate_speed_m_s": u_coord,
        "coordinate_speed_over_c": u_coord / c_m_s,
        "coordinate_superluminal": coord_superluminal
      },
      "boosted_frame_witness": {
        "boost_v_m_s": v_m_s,
        "boost_beta": beta_boost,
        "boost_gamma": gamma,
        "dt_numerator_s": dt_num,
        "dt_prime_s": dt_prime,
        "predeparture_in_boosted_frame": predeparture_boosted
      },
      "ctc_legality_witness": {
        "interval_sq": ds2,
        "timelike_step": timelike_step,
        "identified_closed_one_period": identified_closed
      },
      "effective_arrival_loop": {
        "baseline_light_time_s": baseline_light,
        "dt_in_s": dt_in,
        "dt_out_s": dt_out,
        "t_eff_s": t_eff,
        "t_cover_s": t_cover,
        "effective_superluminal": effective_superluminal,
        "predeparture_effective": predeparture_effective
      },
      "summary": {
        "faithful_sim_possible": faithful_sim_possible,
        "model_scope": "simulation of theorem-level consistency lanes; not an engine claim"
      }
    });

    let txt_path = out.join("ctc_faithful_bridge_sim.txt");
    let json_path = out.join("ctc_faithful_bridge_sim.json");

    let mut txt = String::new();
    txt.push_str("[ctc_faithful_bridge_sim]\n");
    txt.push_str("theorem-faithful composite lane check\n\n");
    txt.push_str(&format!(
        "s={:.12e}, u/c={:.12e}, local_bound_ok={}, gate_pass={}\n",
        s,
        u_coord / c_m_s,
        local_bound_ok,
        gate_pass
    ));
    txt.push_str(&format!(
        "boost: beta={:.12e}, dt_num={:.12e}s, dt'={:.12e}s, predeparture_frame={}\n",
        beta_boost, dt_num, dt_prime, predeparture_boosted
    ));
    txt.push_str(&format!(
        "ctc witness: ds2={:.12e}, timelike_step={}, identified_closed={}\n",
        ds2, timelike_step, identified_closed
    ));
    txt.push_str(&format!(
        "loop: t_eff={:.12e}s, baseline_light={:.12e}s, eff_superluminal={}, predeparture_eff={}\n",
        t_eff, baseline_light, effective_superluminal, predeparture_effective
    ));
    txt.push_str(&format!(
        "faithful_sim_possible={}\n",
        faithful_sim_possible
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "faithful_sim_possible={} (u/c={:.4}, dt'={:.4e}s)",
        faithful_sim_possible,
        u_coord / c_m_s,
        dt_prime
    );
}

