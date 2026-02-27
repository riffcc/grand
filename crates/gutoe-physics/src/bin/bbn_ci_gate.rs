//! Big Bang Nucleosynthesis CI gate for GRAND-349.

use gutoe_physics::{evaluate_bbn_gate, BbnWindows};
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn main() {
    let windows = BbnWindows::default();
    let score = evaluate_bbn_gate(windows);

    let out_dir =
        std::env::var("GUTOE_BBN_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/bbn_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"windows\": {{\"yp_abs_max\": {:.12}, \"dh_rel_max\": {:.12}, \"he3_rel_max\": {:.12}, \"li_tension_ratio_min\": {:.12}, \"li_tension_ratio_max\": {:.12}}},\n  \"score\": {{\"eta10\": {:.12}, \"yp_pred\": {:.12e}, \"dh_pred\": {:.12e}, \"he3h_pred\": {:.12e}, \"li7h_pred\": {:.12e}, \"yp_delta\": {:.12e}, \"dh_rel_error\": {:.12}, \"he3_rel_error\": {:.12}, \"li_tension_ratio\": {:.12}, \"yp_ok\": {}, \"dh_ok\": {}, \"he3_ok\": {}, \"li_tension_ok\": {}, \"passes_primary\": {}, \"passes_all\": {}}}\n}}",
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
        score.yp_ok,
        score.dh_ok,
        score.he3_ok,
        score.li_tension_ok,
        score.passes_primary(),
        score.passes_all(),
    )
    .expect("write gate json");

    println!(
        "BBN gate: pass={} (η10={:.4}, YpΔ={:.3e}, D/H err={:.4}, 3He/H err={:.4}, Li tension ratio={:.3})",
        score.passes_all(),
        score.eta10,
        score.yp_delta,
        score.dh_rel_error,
        score.he3_rel_error,
        score.li_tension_ratio
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
