//! Closed-loop temporal-loan probe.
//!
//! Demonstrates a two-phase closed packet cycle:
//! - Phase A: positive local export via door drawdown.
//! - Phase B: repayment (negative export) restoring door state.
//!
//! Per phase conservation identity:
//!   Ein + Eprev = Eout + Enext + Export + Loss
//! with closed packet condition `Ein = Eout`.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn residual(ein: f64, eprev: f64, eout: f64, enext: f64, export: f64, loss: f64) -> f64 {
    (ein + eprev) - (eout + enext + export + loss)
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_CLOSED_LOAN_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_closed_loop_temporal_loan_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let e0 = env_f64("GUTOE_CTC_CLOSED_LOAN_E0_J", 10.0).max(0.0);
    let draw = env_f64("GUTOE_CTC_CLOSED_LOAN_DRAWDOWN_J", 1.0).max(0.0).min(e0);
    let loss_a = env_f64("GUTOE_CTC_CLOSED_LOAN_LOSS_A_J", 0.0).max(0.0);
    let loss_b = env_f64("GUTOE_CTC_CLOSED_LOAN_LOSS_B_J", 0.0).max(0.0);

    // Closed packet per phase.
    let ein_a = 0.0;
    let eout_a = 0.0;
    let eprev_a = e0;
    let enext_a = e0 - draw;
    let export_a = draw - loss_a;

    let ein_b = 0.0;
    let eout_b = 0.0;
    let eprev_b = enext_a;
    let enext_b = e0;
    let export_b = (eprev_b - enext_b) - loss_b; // closed-cycle formula

    let res_a = residual(ein_a, eprev_a, eout_a, enext_a, export_a, loss_a);
    let res_b = residual(ein_b, eprev_b, eout_b, enext_b, export_b, loss_b);

    let closed_packet_a = (ein_a - eout_a).abs() < 1e-12;
    let closed_packet_b = (ein_b - eout_b).abs() < 1e-12;

    let phase_a_positive = export_a > 0.0;
    let phase_b_negative = export_b < 0.0;

    let door_restored = (enext_b - e0).abs() < 1e-12;
    let net_export = export_a + export_b;
    let net_loss = loss_a + loss_b;
    let net_export_plus_loss = net_export + net_loss;

    let payload = json!({
      "inputs": {
        "e0_j": e0,
        "drawdown_j": draw,
        "loss_a_j": loss_a,
        "loss_b_j": loss_b
      },
      "phase_a": {
        "ein_j": ein_a,
        "eprev_j": eprev_a,
        "eout_j": eout_a,
        "enext_j": enext_a,
        "export_j": export_a,
        "loss_j": loss_a,
        "closed_packet": closed_packet_a,
        "positive_export_window": phase_a_positive,
        "conservation_residual_j": res_a
      },
      "phase_b": {
        "ein_j": ein_b,
        "eprev_j": eprev_b,
        "eout_j": eout_b,
        "enext_j": enext_b,
        "export_j": export_b,
        "loss_j": loss_b,
        "closed_packet": closed_packet_b,
        "repayment_negative_export": phase_b_negative,
        "conservation_residual_j": res_b
      },
      "cycle_summary": {
        "door_restored": door_restored,
        "net_export_j": net_export,
        "net_loss_j": net_loss,
        "net_export_plus_loss_j": net_export_plus_loss,
        "closed_loop_temporal_loan_pattern": phase_a_positive && phase_b_negative && door_restored
      },
      "interpretation": "phase-A positive local export can coexist with global closure; losses set net cycle export floor"
    });

    let txt_path = out.join("ctc_closed_loop_temporal_loan_probe.txt");
    let json_path = out.join("ctc_closed_loop_temporal_loan_probe.json");

    let mut txt = String::new();
    txt.push_str("[ctc_closed_loop_temporal_loan_probe]\n");
    txt.push_str("identity: Ein + Eprev = Eout + Enext + Export + Loss\n\n");
    txt.push_str(&format!(
        "phaseA: export={:.12e}J, closed_packet={}, residual={:.12e}\n",
        export_a, closed_packet_a, res_a
    ));
    txt.push_str(&format!(
        "phaseB: export={:.12e}J, closed_packet={}, residual={:.12e}\n",
        export_b, closed_packet_b, res_b
    ));
    txt.push_str(&format!(
        "cycle: door_restored={}, net_export={:.12e}J, net_loss={:.12e}J, net_export_plus_loss={:.12e}J\n",
        door_restored, net_export, net_loss, net_export_plus_loss
    ));
    txt.push_str(&format!(
        "closed_loop_temporal_loan_pattern={}\n",
        phase_a_positive && phase_b_negative && door_restored
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}

