//! Multi-cannabinoid panel report (CB1-focused neuron lane).

use gutoe_physics::{default_cannabinoid_specs, evaluate_cannabinoid_panel, NeuronCouplingInput};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let temperature_k = env_f64("GUTOE_CANNABINOID_TEMP_K", 310.15);
    let coupling = NeuronCouplingInput {
        intrinsic_efficacy: 0.55,
        max_release_inhibition_fraction: env_f64("GUTOE_CANNABINOID_MAX_RELEASE_INHIBITION", 0.75),
        max_firing_suppression_fraction: env_f64("GUTOE_CANNABINOID_MAX_FIRING_SUPPRESSION", 0.45),
        hill_coefficient: env_f64("GUTOE_CANNABINOID_HILL", 1.0),
        baseline_release_probability: env_f64("GUTOE_CANNABINOID_BASELINE_RELEASE_P", 0.35),
        baseline_firing_rate_hz: env_f64("GUTOE_CANNABINOID_BASELINE_FIRING_HZ", 8.0),
    };

    let specs = default_cannabinoid_specs();
    let mut rows = evaluate_cannabinoid_panel(&specs, temperature_k, coupling);
    rows.sort_by(|a, b| {
        a.ki_cb1_nanomolar
            .partial_cmp(&b.ki_cb1_nanomolar)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let out_dir = std::env::var("GUTOE_CANNABINOID_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/cannabinoid_panel".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("cannabinoid_panel_report.txt");
    let csv_path = out.join("cannabinoid_panel_report.csv");
    let json_path = out.join("cannabinoid_panel_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[cannabinoid_panel]").expect("write");
    writeln!(txt, "count = {}", rows.len()).expect("write");
    writeln!(txt, "temperature_k = {:.6}", temperature_k).expect("write");
    writeln!(txt, "baseline_firing_rate_hz = {:.6}", coupling.baseline_firing_rate_hz).expect("write");
    writeln!(txt, "baseline_release_probability = {:.6}", coupling.baseline_release_probability)
        .expect("write");

    let mut csv = String::from(
        "name,class,ki_cb1_nM,ki_cb2_nM,intrinsic_efficacy_cb1,experimental_delta_g_kj_mol,qed_floor_total_kj_mol,residual_required_kj_mol,residual_modeled_total_kj_mol,residual_closure_error_kj_mol,explained_fraction_of_abs_delta_g,occupancy_10nM,occupancy_30nM,occupancy_100nM,firing_scale_10nM,firing_scale_30nM,firing_scale_100nM\n",
    );

    for r in &rows {
        csv.push_str(&format!(
            "{},{},{:.6},{:.6},{:.6},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            r.name,
            r.class,
            r.ki_cb1_nanomolar,
            r.ki_cb2_nanomolar,
            r.intrinsic_efficacy_cb1,
            r.experimental_delta_g_kj_mol,
            r.qed_floor_total_kj_mol,
            r.residual_required_kj_mol,
            r.residual_modeled_total_kj_mol,
            r.residual_closure_error_kj_mol,
            r.explained_fraction_of_abs_delta_g,
            r.occupancy_10nm,
            r.occupancy_30nm,
            r.occupancy_100nm,
            r.firing_scale_10nm,
            r.firing_scale_30nm,
            r.firing_scale_100nm
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    let mean_explained = rows
        .iter()
        .map(|r| r.explained_fraction_of_abs_delta_g)
        .sum::<f64>()
        / rows.len().max(1) as f64;
    let mean_abs_residual_error = rows
        .iter()
        .map(|r| r.residual_closure_error_kj_mol.abs())
        .sum::<f64>()
        / rows.len().max(1) as f64;

    let top_occupancy_100 = rows
        .iter()
        .max_by(|a, b| {
            a.occupancy_100nm
                .partial_cmp(&b.occupancy_100nm)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied();
    let strongest_firing_suppression_100 = rows
        .iter()
        .min_by(|a, b| {
            a.firing_scale_100nm
                .partial_cmp(&b.firing_scale_100nm)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied();

    let payload = json!({
        "meta": {
            "lane": "multi_cannabinoid_cb1_panel",
            "count": rows.len(),
            "temperature_k": temperature_k,
            "notes": [
                "Ki priors are literature-scale seeds and may vary by assay/protocol",
                "this lane is comparative and reduced-order, not full receptor MD/QM closure"
            ]
        },
        "coupling": {
            "max_release_inhibition_fraction": coupling.max_release_inhibition_fraction,
            "max_firing_suppression_fraction": coupling.max_firing_suppression_fraction,
            "hill_coefficient": coupling.hill_coefficient,
            "baseline_release_probability": coupling.baseline_release_probability,
            "baseline_firing_rate_hz": coupling.baseline_firing_rate_hz
        },
        "summary": {
            "mean_explained_fraction_of_abs_delta_g": mean_explained,
            "mean_abs_residual_closure_error_kj_mol": mean_abs_residual_error,
            "top_occupancy_100nM": top_occupancy_100.map(|r| json!({
                "name": r.name, "occupancy_100nM": r.occupancy_100nm
            })),
            "strongest_firing_suppression_100nM": strongest_firing_suppression_100.map(|r| json!({
                "name": r.name, "firing_scale_100nM": r.firing_scale_100nm
            }))
        },
        "rows": rows.iter().map(|r| json!({
            "name": r.name,
            "class": r.class,
            "ki_cb1_nM": r.ki_cb1_nanomolar,
            "ki_cb2_nM": r.ki_cb2_nanomolar,
            "intrinsic_efficacy_cb1": r.intrinsic_efficacy_cb1,
            "experimental_delta_g_kj_mol": r.experimental_delta_g_kj_mol,
            "qed_floor_total_kj_mol": r.qed_floor_total_kj_mol,
            "residual_required_kj_mol": r.residual_required_kj_mol,
            "residual_modeled_total_kj_mol": r.residual_modeled_total_kj_mol,
            "residual_closure_error_kj_mol": r.residual_closure_error_kj_mol,
            "explained_fraction_of_abs_delta_g": r.explained_fraction_of_abs_delta_g,
            "occupancy_10nM": r.occupancy_10nm,
            "occupancy_30nM": r.occupancy_30nm,
            "occupancy_100nM": r.occupancy_100nm,
            "firing_scale_10nM": r.firing_scale_10nm,
            "firing_scale_30nM": r.firing_scale_30nm,
            "firing_scale_100nM": r.firing_scale_100nm
        })).collect::<Vec<_>>()
    });
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "cannabinoid_panel: count={}, mean_explained={:.3}, mean_abs_residual_error={:.3}",
        rows.len(),
        mean_explained,
        mean_abs_residual_error
    );
}
