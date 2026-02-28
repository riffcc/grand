//! Quantitative baryogenesis CI gate (GRAND-348).

use gutoe_physics::{evaluate_baryogenesis_gate, BaryogenesisWindows};
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn main() {
    let windows = BaryogenesisWindows::default();
    let score = evaluate_baryogenesis_gate(windows);

    let out_dir =
        std::env::var("GUTOE_BARYO_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/baryogenesis_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"eta_rel_error_max\": {:.12},\n  \"score\": {{\"jarlskog_ckm_direct\": {:.12e}, \"jarlskog_ckm_texture\": {:.12e}, \"pmns_theta23_alpha2_c\": {:.12e}, \"leptogenesis_pmns_gain\": {:.12e}, \"leptogenesis_pmns_scalar\": {:.12e}, \"leptogenesis_multiplier\": {:.12e}, \"eta_predicted\": {:.12e}, \"eta_observed\": {:.12e}, \"eta_rel_error\": {:.12}, \"cp_violation_ok\": {}, \"baryon_violation_channel_ok\": {}, \"nonequilibrium_ok\": {}, \"eta_window_ok\": {}, \"sakharov_ok\": {}, \"passes_all\": {}}}\n}}",
        score.passes_all(),
        windows.eta_rel_error_max,
        score.jarlskog_ckm_direct,
        score.jarlskog_ckm_texture,
        score.pmns_theta23_alpha2_c,
        score.leptogenesis_pmns_gain,
        score.leptogenesis_pmns_scalar,
        score.leptogenesis_multiplier,
        score.eta_predicted,
        score.eta_observed,
        score.eta_rel_error,
        score.cp_violation_ok,
        score.baryon_violation_channel_ok,
        score.nonequilibrium_ok,
        score.eta_window_ok,
        score.sakharov_ok(),
        score.passes_all()
    )
    .expect("write gate json");

    println!(
        "baryogenesis gate: pass={} (η_pred={:.6e}, η_obs={:.6e}, rel_err={:.4})",
        score.passes_all(),
        score.eta_predicted,
        score.eta_observed,
        score.eta_rel_error
    );
    println!("wrote {json_path}");

    if !score.passes_all() {
        eprintln!(
            "FAIL: cp={} baryon_channel={} nonequilibrium={} eta_window={}",
            score.cp_violation_ok,
            score.baryon_violation_channel_ok,
            score.nonequilibrium_ok,
            score.eta_window_ok
        );
        process::exit(2);
    }
}
