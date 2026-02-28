//! Phage-host matching report lane.

use gutoe_physics::default_phage_host_matching_panel;
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
    let temperature_k = env_f64("GUTOE_PHAGE_MATCH_TEMP_K", 310.15);
    let mut panel = default_phage_host_matching_panel(temperature_k);
    panel.rows.sort_by(|a, b| {
        a.strain_name
            .cmp(&b.strain_name)
            .then_with(|| {
                b.lysis_potential_score
                    .partial_cmp(&a.lysis_potential_score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| a.phage_name.cmp(&b.phage_name))
    });

    let out_dir = std::env::var("GUTOE_PHAGE_MATCH_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/phage_host_matching".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("phage_host_matching_report.txt");
    let csv_path = out.join("phage_host_matching_report.csv");
    let json_path = out.join("phage_host_matching_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[phage_host_matching]").expect("write");
    writeln!(txt, "temperature_k = {:.6}", temperature_k).expect("write");
    writeln!(txt, "pair_count = {}", panel.rows.len()).expect("write");
    writeln!(
        txt,
        "mean_best_lysis_score = {:.9}",
        panel.mean_best_lysis_score
    )
    .expect("write");
    writeln!(
        txt,
        "resistance_independence_probe_abs_delta = {:.12e}",
        panel.resistance_independence_probe_abs_delta
    )
    .expect("write");
    for b in &panel.best_by_strain {
        writeln!(
            txt,
            "best[{}|{}]: phage={} lysis_score={:.6} kd_nM={:.3}",
            b.strain_name,
            b.resistance_marker,
            b.best_phage_name,
            b.best_lysis_score,
            b.best_predicted_kd_nanomolar
        )
        .expect("write");
    }

    let mut csv = String::from(
        "strain,species,resistance_marker,phage,family,receptor_match_score,predicted_kd_nM,attachment_prob,lysis_potential_score,predicted_delta_g_kj_mol,qed_floor_total_kj_mol,residual_modeled_total_kj_mol\n",
    );
    for r in &panel.rows {
        csv.push_str(&format!(
            "{},{},{},{},{},{:.6},{:.6},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            r.strain_name,
            r.strain_species,
            r.resistance_marker,
            r.phage_name,
            r.phage_family,
            r.receptor_match_score,
            r.predicted_kd_nanomolar,
            r.attachment_prob,
            r.lysis_potential_score,
            r.predicted_delta_g_kj_mol,
            r.qed_floor_total_kj_mol,
            r.residual_modeled_total_kj_mol
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    let payload = json!({
        "meta": {
            "lane": "phage_host_matching",
            "temperature_k": temperature_k,
            "notes": [
                "Reduced-order tail-fiber receptor binding lane",
                "Mechanism bypasses beta-lactamase enzyme activity path",
                "Simulation ranking artifact, not clinical guidance"
            ]
        },
        "summary": {
            "pair_count": panel.rows.len(),
            "mean_best_lysis_score": panel.mean_best_lysis_score,
            "resistance_independence_probe_abs_delta": panel.resistance_independence_probe_abs_delta
        },
        "best_by_strain": panel.best_by_strain.iter().map(|b| json!({
            "strain": b.strain_name,
            "resistance_marker": b.resistance_marker,
            "best_phage": b.best_phage_name,
            "best_lysis_score": b.best_lysis_score,
            "best_predicted_kd_nM": b.best_predicted_kd_nanomolar
        })).collect::<Vec<_>>(),
        "rows": panel.rows.iter().map(|r| json!({
            "strain": r.strain_name,
            "species": r.strain_species,
            "resistance_marker": r.resistance_marker,
            "phage": r.phage_name,
            "family": r.phage_family,
            "receptor_match_score": r.receptor_match_score,
            "predicted_kd_nM": r.predicted_kd_nanomolar,
            "attachment_prob": r.attachment_prob,
            "lysis_potential_score": r.lysis_potential_score,
            "predicted_delta_g_kj_mol": r.predicted_delta_g_kj_mol,
            "qed_ionic_floor_kj_mol": r.qed_ionic_floor_kj_mol,
            "qed_hbond_floor_kj_mol": r.qed_hbond_floor_kj_mol,
            "qed_floor_total_kj_mol": r.qed_floor_total_kj_mol,
            "residual_modeled_total_kj_mol": r.residual_modeled_total_kj_mol
        })).collect::<Vec<_>>()
    });
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("serialize"))
        .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "phage_host_matching: pairs={} mean_best_lysis={:.3} probe_abs_delta={:.3e}",
        panel.rows.len(),
        panel.mean_best_lysis_score,
        panel.resistance_independence_probe_abs_delta
    );
}
