//! Big Bang Nucleosynthesis CI gate for GRAND-349.

use gutoe_physics::{
    evaluate_bbn_gate, lithium7_be7_additional_destruction_multiplier,
    lithium7_be7_dark_suppression_required, lithium7_be7_destruction_enhancement_required,
    lithium7_residual_post_bbn_depletion_factor, BbnWindows,
    LITHIUM7_BE7_DESTRUCTION_ENHANCEMENT_STRUCTURAL,
};
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn main() {
    let windows = BbnWindows::default();
    let score = evaluate_bbn_gate(windows);
    let li7_be7_required_suppression = lithium7_be7_dark_suppression_required(score.eta10);
    let li7_be7_required_enhancement = lithium7_be7_destruction_enhancement_required(score.eta10);
    let li7_be7_additional_multiplier = lithium7_be7_additional_destruction_multiplier(score.eta10);
    let li7_residual_depletion_factor = lithium7_residual_post_bbn_depletion_factor(score.eta10);
    let li7_residual_depletion_percent = 100.0 * (1.0 - li7_residual_depletion_factor);

    let out_dir =
        std::env::var("GUTOE_BBN_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/bbn_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"windows\": {{\"yp_abs_max\": {:.12}, \"dh_rel_max\": {:.12}, \"he3_rel_max\": {:.12}, \"li_tension_ratio_min\": {:.12}, \"li_tension_ratio_max\": {:.12}}},\n  \"score\": {{\"eta10\": {:.12}, \"yp_pred\": {:.12e}, \"dh_pred\": {:.12e}, \"he3h_pred\": {:.12e}, \"li7h_pred\": {:.12e}, \"yp_delta\": {:.12e}, \"dh_rel_error\": {:.12}, \"he3_rel_error\": {:.12}, \"li_tension_ratio\": {:.12}, \"li7_be7_dark_suppression_required\": {:.12}, \"li7_be7_destruction_enhancement_structural\": {:.12}, \"li7_be7_destruction_enhancement_required\": {:.12}, \"li7_be7_additional_destruction_multiplier\": {:.12}, \"li7_residual_post_bbn_depletion_factor\": {:.12}, \"li7_residual_post_bbn_depletion_percent\": {:.12}, \"yp_ok\": {}, \"dh_ok\": {}, \"he3_ok\": {}, \"li_tension_ok\": {}, \"passes_primary\": {}, \"passes_all\": {}}}\n}}",
        score.passes_all(),
        windows.yp_abs_max,
        windows.dh_rel_max,
        windows.he3_rel_max,
        windows.li_tension_ratio_min,
        windows.li_tension_ratio_max,
        score.eta10,
        score.yp_pred,
        score.dh_pred,
        score.he3h_pred,
        score.li7h_pred,
        score.yp_delta,
        score.dh_rel_error,
        score.he3_rel_error,
        score.li_tension_ratio,
        li7_be7_required_suppression,
        LITHIUM7_BE7_DESTRUCTION_ENHANCEMENT_STRUCTURAL,
        li7_be7_required_enhancement,
        li7_be7_additional_multiplier,
        li7_residual_depletion_factor,
        li7_residual_depletion_percent,
        score.yp_ok,
        score.dh_ok,
        score.he3_ok,
        score.li_tension_ok,
        score.passes_primary(),
        score.passes_all(),
    )
    .expect("write gate json");

    println!(
        "BBN gate: pass={} (η10={:.4}, YpΔ={:.3e}, D/H err={:.4}, 3He/H err={:.4}, Li tension ratio={:.3}, Be7 extra destruction x{:.3})",
        score.passes_all(),
        score.eta10,
        score.yp_delta,
        score.dh_rel_error,
        score.he3_rel_error,
        score.li_tension_ratio,
        li7_be7_additional_multiplier
    );
    println!("wrote {json_path}");

    if !score.passes_all() {
        eprintln!(
            "FAIL: yp_ok={} dh_ok={} he3_ok={} li_tension_ok={}",
            score.yp_ok, score.dh_ok, score.he3_ok, score.li_tension_ok
        );
        process::exit(2);
    }
}
