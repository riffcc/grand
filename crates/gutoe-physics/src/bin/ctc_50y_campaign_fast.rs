//! Accelerated 50-year CTC/transport campaign simulator.
//!
//! Purpose:
//! - Encode a long-horizon program (theory -> prototypes -> operations) in one executable.
//! - Keep faith with current theorem lanes:
//!   - rear shortcut factor `s = 1/10`
//!   - budget threshold scale `~ (3/16)|R||T|`
//!   - local timelike routing (`beta < 1`)
//! - Simulate infrastructure, device fleet, and transport of simulated beings.
//!
//! Scope:
//! - Systems simulation and planning aid.
//! - Not a physical engine claim.

use gutoe_physics::constants::C;
use serde_json::json;
use std::fs;
use std::path::PathBuf;

const FC_VOID: f64 = 3.0 / 16.0;
const REAR_SHORTCUT_FACTOR: f64 = 1.0 / 10.0;
const V_EWSB_GEV: f64 = 245.3;
const GEV_TO_J: f64 = 1.602_176_634e-10;
const HBARC_GEV_M: f64 = 0.197_326_980_4e-15;
const SECONDS_PER_YEAR: f64 = 365.25 * 24.0 * 3600.0;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_u64(name: &str, default: u64) -> u64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(default)
}

fn clamp01(x: f64) -> f64 {
    x.clamp(0.0, 1.0)
}

fn wall_tension_front_j_m2(delta_theta: f64, thickness_m: f64) -> f64 {
    let l_nat = thickness_m / HBARC_GEV_M;
    let sigma_gev3 = FC_VOID * V_EWSB_GEV.powi(2) * delta_theta.powi(2) / (2.0 * l_nat);
    let gev3_to_j_m2 = GEV_TO_J / HBARC_GEV_M.powi(2);
    sigma_gev3 * gev3_to_j_m2
}

fn derived_kappa_j_per_m_s() -> f64 {
    let thickness_m = HBARC_GEV_M / V_EWSB_GEV;
    let sigma_front = wall_tension_front_j_m2(std::f64::consts::PI, thickness_m);
    let sigma_rear = REAR_SHORTCUT_FACTOR * sigma_front;
    2.0 * std::f64::consts::PI * C * sigma_rear / FC_VOID
}

fn structural_threshold_j(kappa_j_per_m_s: f64, radius_m: f64, period_s: f64) -> f64 {
    kappa_j_per_m_s * FC_VOID * radius_m.abs() * period_s.abs()
}

fn phase_name(year: u64) -> &'static str {
    match year {
        0..=5 => "phase_1_theory_and_metrology",
        6..=15 => "phase_2_micro_patch_prototypes",
        16..=30 => "phase_3_mesoscale_device_and_safety",
        31..=40 => "phase_4_operational_trials",
        _ => "phase_5_networked_operations",
    }
}

fn phase_multipliers(year: u64) -> (f64, f64, f64) {
    // (fabrication_efficiency, mission_load, risk_penalty)
    match year {
        0..=5 => (0.05, 0.00, 0.70),
        6..=15 => (0.20, 0.05, 0.45),
        16..=30 => (0.45, 0.20, 0.28),
        31..=40 => (0.70, 0.60, 0.16),
        _ => (0.85, 1.00, 0.08),
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_50Y_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_50y_campaign_fast".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let years = env_u64("GUTOE_CTC_50Y_YEARS", 50).max(1);
    let c_m_s = env_f64("GUTOE_CTC_50Y_C_M_S", C).abs().max(1e-12);
    let beta_local = env_f64("GUTOE_CTC_50Y_LOCAL_BETA", 0.80).clamp(1e-6, 0.999_999);
    let s = env_f64("GUTOE_CTC_50Y_SHORTCUT_FACTOR", REAR_SHORTCUT_FACTOR)
        .abs()
        .clamp(1e-9, 1.0 - 1e-9);

    // Infrastructure and program controls
    let mut power_w = env_f64("GUTOE_CTC_50Y_POWER0_W", 5.0e8).max(1.0);
    let base_power_growth = env_f64("GUTOE_CTC_50Y_POWER_GROWTH", 0.18).clamp(-0.5, 1.0);
    let budget_frac = env_f64("GUTOE_CTC_50Y_BUDGET_FRAC", 0.08).clamp(0.0, 1.0);
    let maint_frac = env_f64("GUTOE_CTC_50Y_MAINT_FRAC", 0.015).clamp(0.0, 1.0);
    let reinvest_frac = env_f64("GUTOE_CTC_50Y_REINVEST_FRAC", 0.50).clamp(0.0, 1.0);

    // Device geometry lane
    let period_s = env_f64("GUTOE_CTC_50Y_PERIOD_S", 5.366_854_127e-27)
        .abs()
        .max(1e-30);
    let mut radius_m = env_f64("GUTOE_CTC_50Y_RADIUS0_M", 8.044_312_286_995_516e-19)
        .abs()
        .max(1e-30);
    let radius_growth_per_year = env_f64("GUTOE_CTC_50Y_RADIUS_GROWTH", 0.07).clamp(0.0, 2.0);
    let loops_per_mission = env_f64("GUTOE_CTC_50Y_LOOPS_PER_MISSION", 120.0).max(1.0);
    let retry_depth_cap = env_u64("GUTOE_CTC_50Y_RETRY_DEPTH_CAP", 12).max(1) as f64;

    // Simulated beings lane
    let mut population = env_f64("GUTOE_CTC_50Y_POP0", 10_000.0).max(1.0);
    let pop_growth = env_f64("GUTOE_CTC_50Y_POP_GROWTH", 0.10).clamp(-0.5, 2.0);
    let mission_unit = env_f64("GUTOE_CTC_50Y_MISSION_UNIT", 50.0).max(1.0);

    let kappa = env_f64("GUTOE_CTC_50Y_KAPPA_J_PER_M_S", derived_kappa_j_per_m_s()).max(0.0);
    let u_coord = c_m_s / s;
    let local_bound_ok = beta_local < 1.0;
    let coordinate_superluminal = u_coord > c_m_s;

    let mut doors_active = 0.0_f64;
    let mut infra_index = 0.05_f64;
    let mut safety_index = 0.45_f64;
    let mut transported_total = 0.0_f64;
    let mut transported_eventual_total = 0.0_f64;
    let mut transported_asymptotic_total = 0.0_f64;
    let mut missions_total = 0.0_f64;
    let mut missions_eventual_success_total = 0.0_f64;
    let mut missions_asymptotic_success_total = 0.0_f64;
    let mut predeparture_missions = 0.0_f64;
    let mut mission_failures_first_pass = 0.0_f64;
    let mut energy_invested_total_j = 0.0_f64;
    let mut energy_maint_total_j = 0.0_f64;
    let mut peak_power_w = power_w;

    let mut rows = Vec::new();

    for year in 0..years {
        let (fab_eff, mission_load, risk_penalty) = phase_multipliers(year);
        let phase = phase_name(year);

        let annual_energy_j = power_w * SECONDS_PER_YEAR;
        let program_budget_j = annual_energy_j * budget_frac;
        let threshold_j = structural_threshold_j(kappa, radius_m, period_s).max(1e-30);
        // Fabrication bottleneck: budget alone does not define throughput.
        let ideal_new_doors = ((program_budget_j / threshold_j) * fab_eff).max(0.0);
        let fab_cap = (25.0 + 140.0 * infra_index + 22.0 * year as f64) * (0.5 + 0.5 * fab_eff);
        let new_doors = ideal_new_doors.min(fab_cap.max(0.0));
        doors_active += new_doors;

        // Effective loop gain uses per-door budget, not whole-program budget.
        let budget_per_door_j = if doors_active > 0.0 {
            program_budget_j / doors_active
        } else {
            0.0
        };
        let q_eff = budget_per_door_j / threshold_j;

        // Maintenance
        let maint_j = doors_active * threshold_j * maint_frac;
        energy_maint_total_j += maint_j;
        let net_ops_budget_j = (program_budget_j - maint_j).max(0.0);
        energy_invested_total_j += program_budget_j;

        // Mission envelope
        let mission_capacity = doors_active * mission_load * (4.0 + 40.0 * infra_index);
        let mission_demand = (population / mission_unit).max(1.0);
        let missions = mission_capacity.min(mission_demand).max(0.0);
        missions_total += missions;

        let mission_success_prob = clamp01(
            0.45 + 0.45 * safety_index + 0.15 * fab_eff + 0.10 * infra_index
                - risk_penalty
                - 0.10 * (1.0 - local_bound_ok as u8 as f64),
        );
        let missions_ok = missions * mission_success_prob;
        let missions_fail = missions - missions_ok;
        mission_failures_first_pass += missions_fail;

        // Path-A style effective-arrival criterion using q_eff and loop count.
        // Baseline access+egress model is represented by normalized constant 1.0 loop period.
        let t_eff_norm = 1.0 + loops_per_mission * (1.0 - q_eff);
        let is_predeparture = t_eff_norm < 0.0;
        let predeparture_enabled = is_predeparture && year >= 31 && safety_index >= 0.75;
        let pre_missions = if predeparture_enabled { missions_ok } else { 0.0 };
        predeparture_missions += pre_missions;

        // Retro-retry closure metrics.
        let retry_depth = if predeparture_enabled {
            (1.0 + retry_depth_cap * safety_index * infra_index).floor().max(1.0)
        } else {
            1.0
        };
        let eventual_success_prob = 1.0 - (1.0 - mission_success_prob).powf(retry_depth);
        let asymptotic_success_prob = if predeparture_enabled && mission_success_prob > 0.0 {
            1.0
        } else {
            mission_success_prob
        };
        let missions_eventual_success = missions * eventual_success_prob;
        let missions_asymptotic_success = missions * asymptotic_success_prob;
        missions_eventual_success_total += missions_eventual_success;
        missions_asymptotic_success_total += missions_asymptotic_success;

        // Simulated beings transported
        let beings_transported = missions_ok * mission_unit;
        let beings_transported_eventual = missions_eventual_success * mission_unit;
        let beings_transported_asymptotic = missions_asymptotic_success * mission_unit;
        transported_total += beings_transported;
        transported_eventual_total += beings_transported_eventual;
        transported_asymptotic_total += beings_transported_asymptotic;

        // Learning and growth dynamics
        let learning_gain = 0.020 * missions_ok.ln_1p();
        let failure_drag = 0.004 * missions_fail.ln_1p();
        safety_index = clamp01(safety_index + 0.03 * fab_eff + learning_gain - failure_drag);
        infra_index = clamp01(infra_index + 0.025 * fab_eff + 0.010 * missions_ok.ln_1p());

        // Reinvestment closes the loop in finite program terms.
        let reinvest_j = net_ops_budget_j * reinvest_frac;
        let power_gain = reinvest_j / SECONDS_PER_YEAR * 0.000_000_04;
        power_w = (power_w * (1.0 + base_power_growth * (0.25 + 0.75 * infra_index)) + power_gain).max(1.0);
        peak_power_w = peak_power_w.max(power_w);

        radius_m *= 1.0 + radius_growth_per_year * (0.2 + 0.8 * fab_eff);
        population *= 1.0 + pop_growth * (0.2 + 0.8 * safety_index);

        rows.push(json!({
            "year": year,
            "phase": phase,
            "power_w": power_w,
            "annual_energy_j": annual_energy_j,
            "program_budget_j": program_budget_j,
            "threshold_j_per_device": threshold_j,
            "q_eff": q_eff,
            "radius_m": radius_m,
            "period_s": period_s,
            "doors_new": new_doors,
            "doors_active": doors_active,
            "missions": missions,
            "missions_success_prob": mission_success_prob,
            "missions_success": missions_ok,
            "missions_failure": missions_fail,
            "retry_depth": retry_depth,
            "eventual_success_prob": eventual_success_prob,
            "asymptotic_success_prob": asymptotic_success_prob,
            "missions_eventual_success": missions_eventual_success,
            "missions_asymptotic_success": missions_asymptotic_success,
            "predeparture_mode": is_predeparture,
            "predeparture_enabled": predeparture_enabled,
            "predeparture_missions": pre_missions,
            "beings_transported": beings_transported,
            "beings_transported_eventual": beings_transported_eventual,
            "beings_transported_asymptotic": beings_transported_asymptotic,
            "population": population,
            "infra_index": infra_index,
            "safety_index": safety_index
        }));
    }

    let avg_success = if missions_total > 0.0 {
        (missions_total - mission_failures_first_pass) / missions_total
    } else {
        0.0
    };
    let eventual_success_rate = if missions_total > 0.0 {
        missions_eventual_success_total / missions_total
    } else {
        0.0
    };
    let asymptotic_success_rate = if missions_total > 0.0 {
        missions_asymptotic_success_total / missions_total
    } else {
        0.0
    };
    let predeparture_fraction = if missions_total > 0.0 {
        predeparture_missions / missions_total
    } else {
        0.0
    };

    let summary = json!({
      "inputs": {
        "years": years,
        "local_beta": beta_local,
        "shortcut_factor_s": s,
        "void_fraction": FC_VOID,
        "rear_shortcut_factor_reference": REAR_SHORTCUT_FACTOR,
        "kappa_j_per_m_s": kappa
      },
      "theorem_faithful_checks": {
        "local_bound_ok": local_bound_ok,
        "coordinate_superluminal": coordinate_superluminal,
        "coordinate_speed_over_c": u_coord / c_m_s
      },
      "campaign_totals": {
        "missions_total": missions_total,
        "mission_failures_first_pass_total": mission_failures_first_pass,
        "mission_success_rate_first_pass": avg_success,
        "mission_success_rate_eventual": eventual_success_rate,
        "mission_success_rate_asymptotic": asymptotic_success_rate,
        "missions_eventual_success_total": missions_eventual_success_total,
        "missions_asymptotic_success_total": missions_asymptotic_success_total,
        "predeparture_missions_total": predeparture_missions,
        "predeparture_fraction": predeparture_fraction,
        "transported_sim_beings_total_first_pass": transported_total,
        "transported_sim_beings_total_eventual": transported_eventual_total,
        "transported_sim_beings_total_asymptotic": transported_asymptotic_total,
        "doors_active_final": doors_active,
        "infra_index_final": infra_index,
        "safety_index_final": safety_index,
        "energy_invested_total_j": energy_invested_total_j,
        "energy_maintenance_total_j": energy_maint_total_j,
        "peak_power_w": peak_power_w
      },
      "rows": rows,
      "scope": "program-level accelerated simulation; not an engine claim"
    });

    let txt_path = out.join("ctc_50y_campaign_fast.txt");
    let json_path = out.join("ctc_50y_campaign_fast.json");

    let mut txt = String::new();
    txt.push_str("[ctc_50y_campaign_fast]\n");
    txt.push_str("accelerated long-horizon campaign simulation\n\n");
    txt.push_str(&format!("years = {}\n", years));
    txt.push_str(&format!("local_bound_ok = {}\n", local_bound_ok));
    txt.push_str(&format!("coordinate_superluminal = {}\n", coordinate_superluminal));
    txt.push_str(&format!("coordinate_speed_over_c = {:.6}\n", u_coord / c_m_s));
    txt.push_str(&format!("missions_total = {:.3}\n", missions_total));
    txt.push_str(&format!("mission_success_rate_first_pass = {:.6}\n", avg_success));
    txt.push_str(&format!(
        "mission_success_rate_eventual = {:.6}\n",
        eventual_success_rate
    ));
    txt.push_str(&format!(
        "mission_success_rate_asymptotic = {:.6}\n",
        asymptotic_success_rate
    ));
    txt.push_str(&format!("predeparture_fraction = {:.6}\n", predeparture_fraction));
    txt.push_str(&format!(
        "transported_sim_beings_total_first_pass = {:.3}\n",
        transported_total
    ));
    txt.push_str(&format!(
        "transported_sim_beings_total_eventual = {:.3}\n",
        transported_eventual_total
    ));
    txt.push_str(&format!(
        "transported_sim_beings_total_asymptotic = {:.3}\n",
        transported_asymptotic_total
    ));
    txt.push_str(&format!("doors_active_final = {:.3}\n", doors_active));
    txt.push_str(&format!("infra_index_final = {:.6}\n", infra_index));
    txt.push_str(&format!("safety_index_final = {:.6}\n", safety_index));
    txt.push_str(&format!("energy_invested_total_j = {:.6e}\n", energy_invested_total_j));
    txt.push_str(&format!(
        "energy_maintenance_total_j = {:.6e}\n",
        energy_maint_total_j
    ));
    txt.push_str(&format!("peak_power_w = {:.6e}\n", peak_power_w));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&summary).expect("json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "campaign complete: missions={:.1}, first_pass={:.1}, eventual={:.1}, predep_frac={:.4}",
        missions_total, transported_total, transported_eventual_total, predeparture_fraction
    );
}
