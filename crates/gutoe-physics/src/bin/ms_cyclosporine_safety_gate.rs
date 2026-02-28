//! Cyclosporine safety gate for MS translational lane.
//!
//! Runs PK bridge uncertainty and checks exposure probabilities against
//! configurable safety windows.

use gutoe_physics::{
    default_cyclosporine_pk_bridge_input, default_cyclosporine_safety_gate_input,
    evaluate_cyclosporine_safety_gate, simulate_cyclosporine_pk_bridge,
    summarize_cyclosporine_pk_bridge,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;
use std::process;

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_u64(key: &str, default: u64) -> u64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(default)
}

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn main() {
    let mut input = default_cyclosporine_pk_bridge_input();
    input.site_target_nanomolar = env_f64("GUTOE_MS_SITE_TARGET_NM", input.site_target_nanomolar);
    input.molecular_weight_g_mol = env_f64("GUTOE_CYCLOSPORINE_MW_G_MOL", input.molecular_weight_g_mol);
    input.blood_to_site_gain_median =
        env_f64("GUTOE_MS_PK_GAIN_MEDIAN", input.blood_to_site_gain_median);
    input.blood_to_site_gain_gsd = env_f64("GUTOE_MS_PK_GAIN_GSD", input.blood_to_site_gain_gsd);
    input.samples = env_usize("GUTOE_MS_PK_SAMPLES", input.samples);
    input.seed = env_u64("GUTOE_MS_PK_SEED", input.seed);

    let mut gate = default_cyclosporine_safety_gate_input();
    gate.windows.target_zone_low_ng_ml =
        env_f64("GUTOE_MS_TARGET_ZONE_LOW_NG_ML", gate.windows.target_zone_low_ng_ml);
    gate.windows.target_zone_high_ng_ml =
        env_f64("GUTOE_MS_TARGET_ZONE_HIGH_NG_ML", gate.windows.target_zone_high_ng_ml);
    gate.windows.renal_caution_ng_ml =
        env_f64("GUTOE_MS_RENAL_CAUTION_NG_ML", gate.windows.renal_caution_ng_ml);
    gate.windows.renal_high_ng_ml =
        env_f64("GUTOE_MS_RENAL_HIGH_NG_ML", gate.windows.renal_high_ng_ml);
    gate.windows.neuro_caution_ng_ml =
        env_f64("GUTOE_MS_NEURO_CAUTION_NG_ML", gate.windows.neuro_caution_ng_ml);
    gate.max_prob_above_renal_caution = env_f64(
        "GUTOE_MS_MAX_P_ABOVE_RENAL_CAUTION",
        gate.max_prob_above_renal_caution,
    );
    gate.max_prob_above_renal_high =
        env_f64("GUTOE_MS_MAX_P_ABOVE_RENAL_HIGH", gate.max_prob_above_renal_high);
    gate.max_prob_above_neuro_caution = env_f64(
        "GUTOE_MS_MAX_P_ABOVE_NEURO_CAUTION",
        gate.max_prob_above_neuro_caution,
    );
    gate.min_prob_in_target_zone =
        env_f64("GUTOE_MS_MIN_P_IN_TARGET_ZONE", gate.min_prob_in_target_zone);

    let ensemble = simulate_cyclosporine_pk_bridge(input);
    let summary = summarize_cyclosporine_pk_bridge(&ensemble);
    let score = evaluate_cyclosporine_safety_gate(&ensemble, gate);

    let out_dir = std::env::var("GUTOE_MS_CYCLOSPORINE_SAFETY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ms_cyclosporine_safety_gate".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("ms_cyclosporine_safety_gate.txt");
    let json_path = out.join("ms_cyclosporine_safety_gate.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ms_cyclosporine_safety_gate]").expect("write");
    writeln!(txt, "overall_pass = {}", score.overall_pass).expect("write");
    writeln!(txt, "p50_ng_mL = {:.9}", summary.p50_ng_ml).expect("write");
    writeln!(txt, "p95_ng_mL = {:.9}", summary.p95_ng_ml).expect("write");
    writeln!(txt, "prob_in_target_zone = {:.9}", score.prob_in_target_zone).expect("write");
    writeln!(txt, "prob_above_renal_caution = {:.9}", score.prob_above_renal_caution).expect("write");
    writeln!(txt, "prob_above_renal_high = {:.9}", score.prob_above_renal_high).expect("write");
    writeln!(txt, "prob_above_neuro_caution = {:.9}", score.prob_above_neuro_caution).expect("write");

    let payload = json!({
        "overall_pass": score.overall_pass,
        "meta": {
            "lane": "ms_cyclosporine_safety_gate",
            "note": "simulation safety gate, not clinical decision support"
        },
        "pk_input": {
            "site_target_nM": input.site_target_nanomolar,
            "molecular_weight_g_mol": input.molecular_weight_g_mol,
            "gain_median": input.blood_to_site_gain_median,
            "gain_gsd": input.blood_to_site_gain_gsd,
            "samples": input.samples,
            "seed": input.seed
        },
        "pk_summary": {
            "p05_ng_mL": summary.p05_ng_ml,
            "p50_ng_mL": summary.p50_ng_ml,
            "p95_ng_mL": summary.p95_ng_ml,
            "mean_ng_mL": summary.mean_ng_ml
        },
        "windows_ng_mL": {
            "target_zone_low": gate.windows.target_zone_low_ng_ml,
            "target_zone_high": gate.windows.target_zone_high_ng_ml,
            "renal_caution": gate.windows.renal_caution_ng_ml,
            "renal_high": gate.windows.renal_high_ng_ml,
            "neuro_caution": gate.windows.neuro_caution_ng_ml
        },
        "thresholds": {
            "max_prob_above_renal_caution": gate.max_prob_above_renal_caution,
            "max_prob_above_renal_high": gate.max_prob_above_renal_high,
            "max_prob_above_neuro_caution": gate.max_prob_above_neuro_caution,
            "min_prob_in_target_zone": gate.min_prob_in_target_zone
        },
        "probabilities": {
            "in_target_zone": score.prob_in_target_zone,
            "above_target_zone": score.prob_above_target_zone,
            "above_renal_caution": score.prob_above_renal_caution,
            "above_renal_high": score.prob_above_renal_high,
            "above_neuro_caution": score.prob_above_neuro_caution
        },
        "gate": {
            "target_zone_ok": score.target_zone_ok,
            "renal_caution_ok": score.renal_caution_ok,
            "renal_high_ok": score.renal_high_ok,
            "neuro_caution_ok": score.neuro_caution_ok
        }
    });
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("serialize"))
        .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "ms_cyclosporine_safety_gate: pass={} p_in_zone={:.3} p_renal_caution={:.3} p_renal_high={:.3}",
        score.overall_pass,
        score.prob_in_target_zone,
        score.prob_above_renal_caution,
        score.prob_above_renal_high
    );

    if !score.overall_pass {
        eprintln!(
            "FAIL: target_ok={} renal_caution_ok={} renal_high_ok={} neuro_ok={}",
            score.target_zone_ok, score.renal_caution_ok, score.renal_high_ok, score.neuro_caution_ok
        );
        process::exit(2);
    }
}
