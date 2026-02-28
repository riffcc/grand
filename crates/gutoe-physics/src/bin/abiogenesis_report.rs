//! Abiogenesis threshold report: theorem-style Kauffman closure lane.

use gutoe_physics::{evaluate_abiogenesis_gate, AbiogenesisWindows};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

const DEFAULT_TEMP_K: f64 = 298.15;

fn main() {
    let temperature_k = std::env::var("GUTOE_ABIOGENESIS_TEMP_K")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(DEFAULT_TEMP_K);
    let windows = AbiogenesisWindows::default();
    let score = evaluate_abiogenesis_gate(windows, temperature_k);

    let out_dir = std::env::var("GUTOE_ABIOGENESIS_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/abiogenesis".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);
    let txt_path = out.join("abiogenesis_report.txt");
    let json_path = out.join("abiogenesis_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[abiogenesis_kauffman_lane]").expect("write");
    writeln!(txt, "temperature_k = {:.6}", temperature_k).expect("write");
    writeln!(txt, "threshold = {:.12}", score.closure.threshold).expect("write");
    writeln!(txt, "N = {:.12}", score.closure.monomer_count).expect("write");
    writeln!(txt, "p_min = {:.12}", score.closure.catalytic_probability).expect("write");
    writeln!(txt, "N_times_p = {:.12}", score.closure.n_times_p).expect("write");
    writeln!(txt, "closure_excess = {:.12}", score.closure.closure_excess).expect("write");
    writeln!(
        txt,
        "N_times_p_lower_3sigma = {:.12}",
        score.inevitability.n_times_p_lower_3sigma
    )
    .expect("write");
    writeln!(
        txt,
        "robust_margin = {:.12}",
        score.inevitability.robust_margin
    )
    .expect("write");
    writeln!(
        txt,
        "thermal_chirality_bias = {:.12e}",
        score.inevitability.thermal_chirality_bias
    )
    .expect("write");
    writeln!(txt, "prebiotic_ok = {}", score.prebiotic_ok).expect("write");
    writeln!(txt, "closure_ok = {}", score.closure_ok).expect("write");
    writeln!(txt, "inevitability_ok = {}", score.inevitability_ok).expect("write");
    writeln!(txt, "overall_pass = {}", score.passes_all()).expect("write");

    let payload = json!({
        "temperature_k": temperature_k,
        "windows": {
            "closure_threshold": windows.closure_threshold,
            "catalytic_probability_min": windows.catalytic_probability_min,
            "robust_margin_min": windows.robust_margin_min
        },
        "prebiotic": {
            "feedstock_species": score.prebiotic.feedstock_species,
            "amino_acid_pool_left": score.prebiotic.amino_acid_pool_left,
            "nucleotide_pool": score.prebiotic.nucleotide_pool,
            "peptide_channels": score.prebiotic.peptide_channels,
            "nucleotide_synthesis_channels": score.prebiotic.nucleotide_synthesis_channels,
            "phosphodiester_channels": score.prebiotic.phosphodiester_channels,
            "k_peptide": score.prebiotic.k_peptide,
            "k_nucleotide": score.prebiotic.k_nucleotide,
            "k_phosphodiester": score.prebiotic.k_phosphodiester,
            "catalytic_probability_lower_bound": score.prebiotic.catalytic_probability_lower_bound
        },
        "closure": {
            "monomer_count": score.closure.monomer_count,
            "catalytic_probability": score.closure.catalytic_probability,
            "n_times_p": score.closure.n_times_p,
            "threshold": score.closure.threshold,
            "closure_excess": score.closure.closure_excess
        },
        "inevitability": {
            "pved_delta_e_ev": score.inevitability.pved_delta_e_ev,
            "thermal_chirality_bias": score.inevitability.thermal_chirality_bias,
            "alpha_rel_uncertainty": score.inevitability.alpha_rel_uncertainty,
            "micro_rel_uncertainty": score.inevitability.micro_rel_uncertainty,
            "network_rel_uncertainty": score.inevitability.network_rel_uncertainty,
            "total_rel_uncertainty": score.inevitability.total_rel_uncertainty,
            "n_times_p_lower_3sigma": score.inevitability.n_times_p_lower_3sigma,
            "robust_margin": score.inevitability.robust_margin
        },
        "gate": {
            "prebiotic_ok": score.prebiotic_ok,
            "closure_ok": score.closure_ok,
            "inevitability_ok": score.inevitability_ok,
            "passes_all": score.passes_all()
        }
    });
    let mut json_file = File::create(&json_path).expect("create json");
    writeln!(
        json_file,
        "{}",
        serde_json::to_string_pretty(&payload).expect("serialize")
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "Abiogenesis lane: N*p={:.6}, lower_3sigma={:.6}, margin={:.6}, pass={}",
        score.closure.n_times_p,
        score.inevitability.n_times_p_lower_3sigma,
        score.inevitability.robust_margin,
        score.passes_all()
    );
}
