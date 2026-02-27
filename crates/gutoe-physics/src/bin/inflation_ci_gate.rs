//! Structural inflation CI gate (GRAND-347).

use gutoe_physics::{evaluate_inflation_gate, InflationWindows};
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn main() {
    let windows = InflationWindows::default();
    let score = evaluate_inflation_gate(windows);

    let out_dir =
        std::env::var("GUTOE_INFLATION_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/inflation_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"windows\": {{\"n_efolds_min\": {:.6}, \"n_efolds_max\": {:.6}, \"n_s_center\": {:.6}, \"n_s_tol\": {:.6}, \"r_max\": {:.6}}},\n  \"score\": {{\"n_efolds\": {:.12}, \"epsilon\": {:.12e}, \"eta\": {:.12e}, \"n_s\": {:.12}, \"r\": {:.12}, \"n_end\": {:.12}, \"expansion_factor\": {:.12e}, \"n_efolds_ok\": {}, \"n_s_ok\": {}, \"r_ok\": {}, \"graceful_exit_ok\": {}, \"passes_all\": {}}}\n}}",
        score.passes_all(),
        windows.n_efolds_min,
        windows.n_efolds_max,
        windows.ns_center,
        windows.ns_tol,
        windows.r_max,
        score.n_efolds,
        score.epsilon,
        score.eta,
        score.n_s,
        score.r,
        score.n_end,
        score.expansion_factor,
        score.n_efolds_ok,
        score.n_s_ok,
        score.r_ok,
        score.graceful_exit_ok,
        score.passes_all()
    )
    .expect("write gate json");

    println!(
        "inflation gate: pass={} (N={:.3}, n_s={:.6}, r={:.6}, exit={})",
        score.passes_all(),
        score.n_efolds,
        score.n_s,
        score.r,
        score.graceful_exit_ok
    );
    println!("wrote {json_path}");

    if !score.passes_all() {
        eprintln!(
            "FAIL: N_ok={} n_s_ok={} r_ok={} graceful_exit_ok={}",
            score.n_efolds_ok, score.n_s_ok, score.r_ok, score.graceful_exit_ok
        );
        process::exit(2);
    }
}
