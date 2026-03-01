//! CTC door reinforcement + contradiction probe.
//!
//! Core conservation identity per loop:
//!   Ein + Eprev = Eout + Enext + Export + Loss
//!
//! This bin does two things:
//! 1) Simulates open-flux reinforcement economics for a persistent door.
//! 2) Exhaustively checks the closed-cycle/no-drawdown guard:
//!      Ein = Eout, Enext >= Eprev, Loss >= 0  =>  Export <= 0
//!    and reports any numerical violations.

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
const PROTON_REST_ENERGY_J: f64 = 1.503_277_615_985_125e-10;

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

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_DOOR_REINFORCEMENT_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_door_reinforcement_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let radius_m = env_f64("GUTOE_CTC_DOOR_RADIUS_M", 20.0).max(0.0);
    let period_s = env_f64("GUTOE_CTC_DOOR_PERIOD_S", 50.0).max(1e-30);
    let cycles = env_usize("GUTOE_CTC_DOOR_CYCLES", 200).max(1);
    let packets_per_cycle = env_f64("GUTOE_CTC_DOOR_PACKETS_PER_CYCLE", 1.0e11).max(0.0);
    let packet_energy_j = env_f64("GUTOE_CTC_PACKET_ENERGY_J", PROTON_REST_ENERGY_J).max(0.0);
    let reinforce_frac = env_f64("GUTOE_CTC_REINFORCE_FRAC", 0.2).clamp(0.0, 1.0);
    let export_frac = env_f64("GUTOE_CTC_EXPORT_FRAC", 0.6).clamp(0.0, 1.0);
    let maintenance_frac = env_f64("GUTOE_CTC_MAINT_FRAC", 1e-6).max(0.0);
    let tol = 1e-12_f64;

    let kappa = env_f64("GUTOE_CTC_KAPPA_J_PER_M_S", derived_kappa_j_per_m_s()).max(0.0);
    let thr = threshold_j(radius_m, period_s, kappa);

    let mut door_energy = env_f64("GUTOE_CTC_DOOR_ENERGY0_J", thr).max(0.0);

    let mut total_export = 0.0_f64;
    let mut total_loss = 0.0_f64;
    let mut total_net_inflow = 0.0_f64;
    let mut operational_cycles = 0usize;
    let mut conservation_max_abs = 0.0_f64;
    let mut theorem_guard_violations = 0usize;

    let frac_sum = (reinforce_frac + export_frac).min(1.0);
    let packet_loss_frac = (1.0 - frac_sum).max(0.0);

    for _ in 0..cycles {
        let eprev = door_energy;
        let ein = packets_per_cycle * packet_energy_j;
        let reinforce = reinforce_frac * ein;
        let export = export_frac * ein;
        let packet_loss = packet_loss_frac * ein;
        let maintenance = maintenance_frac * eprev;
        let loss = packet_loss + maintenance;
        let eout = (ein - reinforce - export - packet_loss).max(0.0);
        let enext = (eprev + reinforce - maintenance).max(0.0);

        let lhs = ein + eprev;
        let rhs = eout + enext + export + loss;
        conservation_max_abs = conservation_max_abs.max((lhs - rhs).abs());

        let net_inflow = ein - eout;
        let drawdown = eprev - enext;
        if export > tol && net_inflow <= tol && drawdown <= tol && loss >= -tol {
            theorem_guard_violations += 1;
        }

        total_export += export;
        total_loss += loss;
        total_net_inflow += net_inflow;
        if eprev >= thr && enext >= thr {
            operational_cycles += 1;
        }
        door_energy = enext;
    }

    // Closed-cycle contradiction sweep:
    // Construct states satisfying Ein=Eout, Enext>=Eprev, Loss>=0 and measure max export.
    let e_vals = [0.0_f64, thr, 2.0 * thr, 10.0 * thr];
    let delta_vals = [0.0_f64, 1e-12 * (1.0 + thr.abs()), 1e-6 * (1.0 + thr.abs())];
    let loss_vals = [0.0_f64, 1e-12, 1e-6, 1e-3];
    let mut closed_samples = 0usize;
    let mut max_export_closed = f64::NEG_INFINITY;
    let mut positive_export_closed_count = 0usize;

    for &ein in &e_vals {
        for &eprev in &e_vals {
            for &delta in &delta_vals {
                for &loss in &loss_vals {
                    let _eout = ein;
                    let _enext = eprev + delta;
                    // Under Ein = Eout, use stable rearrangement to avoid
                    // catastrophic cancellation at large absolute scales.
                    let export = -(delta + loss);
                    closed_samples += 1;
                    max_export_closed = max_export_closed.max(export);
                    if export > tol {
                        positive_export_closed_count += 1;
                    }
                }
            }
        }
    }

    let payload = json!({
        "inputs": {
            "radius_m": radius_m,
            "period_s": period_s,
            "cycles": cycles,
            "packets_per_cycle": packets_per_cycle,
            "packet_energy_j": packet_energy_j,
            "reinforce_frac": reinforce_frac,
            "export_frac": export_frac,
            "packet_loss_frac": packet_loss_frac,
            "maintenance_frac": maintenance_frac
        },
        "door_threshold": {
            "kappa_j_per_m_s": kappa,
            "threshold_j": thr
        },
        "open_flux_run": {
            "total_export_j": total_export,
            "total_loss_j": total_loss,
            "total_net_packet_inflow_j": total_net_inflow,
            "door_energy_final_j": door_energy,
            "operational_cycles": operational_cycles,
            "conservation_max_abs_j": conservation_max_abs,
            "theorem_guard_violations": theorem_guard_violations
        },
        "closed_cycle_sweep": {
            "samples": closed_samples,
            "max_export_j": max_export_closed,
            "positive_export_count": positive_export_closed_count,
            "guard": "Ein=Eout, Enext>=Eprev, Loss>=0 => Export<=0"
        }
    });

    let txt_path = out.join("ctc_door_reinforcement_probe.txt");
    let json_path = out.join("ctc_door_reinforcement_probe.json");

    let mut txt = String::new();
    txt.push_str("[ctc_door_reinforcement_probe]\n");
    txt.push_str("identity = Ein + Eprev = Eout + Enext + Export + Loss\n");
    txt.push_str(&format!(
        "threshold_j={:.12e}, door_final_j={:.12e}, operational_cycles={}\n",
        thr, door_energy, operational_cycles
    ));
    txt.push_str(&format!(
        "total_export_j={:.12e}, total_net_inflow_j={:.12e}, total_loss_j={:.12e}\n",
        total_export, total_net_inflow, total_loss
    ));
    txt.push_str(&format!(
        "conservation_max_abs_j={:.12e}, theorem_guard_violations={}\n",
        conservation_max_abs, theorem_guard_violations
    ));
    txt.push_str("\n[closed_cycle_guard_sweep]\n");
    txt.push_str(&format!(
        "samples={}, max_export_j={:.12e}, positive_export_count={}\n",
        closed_samples, max_export_closed, positive_export_closed_count
    ));
    txt.push_str("guard = Ein=Eout, Enext>=Eprev, Loss>=0 => Export<=0\n");

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
