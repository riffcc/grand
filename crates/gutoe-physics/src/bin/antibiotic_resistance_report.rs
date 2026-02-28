//! Beta-lactamase inhibitor ranking report lane.

use gutoe_physics::default_antibiotic_resistance_panel;
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

fn scaffold_name(s: gutoe_physics::InhibitorScaffold) -> &'static str {
    match s {
        gutoe_physics::InhibitorScaffold::BetaLactamSuicide => "beta_lactam_suicide",
        gutoe_physics::InhibitorScaffold::Diazabicyclooctane => "diazabicyclooctane",
        gutoe_physics::InhibitorScaffold::CyclicBoronate => "cyclic_boronate",
    }
}

fn enzyme_class_name(c: gutoe_physics::BetaLactamaseClass) -> &'static str {
    match c {
        gutoe_physics::BetaLactamaseClass::SerineClassA => "serine_class_a",
        gutoe_physics::BetaLactamaseClass::MetalloClassB => "metallo_class_b",
    }
}

fn main() {
    let temperature_k = env_f64("GUTOE_BETA_LACTAMASE_TEMP_K", 310.15);
    let mut panel = default_antibiotic_resistance_panel(temperature_k);
    panel.rows.sort_by(|a, b| {
        a.enzyme_name
            .cmp(b.enzyme_name)
            .then_with(|| {
                a.predicted_nanomolar
                    .partial_cmp(&b.predicted_nanomolar)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| a.inhibitor_name.cmp(b.inhibitor_name))
    });

    let out_dir = std::env::var("GUTOE_BETA_LACTAMASE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/antibiotic_resistance".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("antibiotic_resistance_report.txt");
    let csv_path = out.join("antibiotic_resistance_report.csv");
    let json_path = out.join("antibiotic_resistance_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[antibiotic_resistance]").expect("write");
    writeln!(txt, "temperature_k = {:.6}", temperature_k).expect("write");
    writeln!(txt, "pair_count = {}", panel.rows.len()).expect("write");
    writeln!(
        txt,
        "mean_abs_log10_error_pred_vs_anchor = {:.9}",
        panel.mean_abs_log10_error
    )
    .expect("write");
    writeln!(
        txt,
        "ndm_max_predicted_occupancy_at_1uM = {:.9}",
        panel.ndm_max_predicted_occupancy_at_1u_m
    )
    .expect("write");
    for b in &panel.best_by_enzyme {
        writeln!(
            txt,
            "best[{}]: anchor={} ({:.3} nM), predicted={} ({:.3} nM), winner_match={}",
            b.enzyme_name,
            b.by_anchor_inhibitor,
            b.by_anchor_nanomolar,
            b.by_predicted_inhibitor,
            b.by_predicted_nanomolar,
            b.predicted_match_anchor_winner
        )
        .expect("write");
    }

    let mut csv = String::from(
        "enzyme,inhibitor,scaffold,enzyme_class,anchor_nM,predicted_nM,log10_error,evidence_count,imputed_anchor,occupancy_anchor_1uM,occupancy_predicted_1uM,anchor_delta_g_kj_mol,predicted_delta_g_kj_mol,qed_floor_total_kj_mol,residual_modeled_total_kj_mol\n",
    );
    for r in &panel.rows {
        csv.push_str(&format!(
            "{},{},{},{},{:.6},{:.6},{:.6},{},{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            r.enzyme_name,
            r.inhibitor_name,
            scaffold_name(r.scaffold),
            enzyme_class_name(r.enzyme_class),
            r.anchor_nanomolar,
            r.predicted_nanomolar,
            r.log10_error_pred_vs_anchor,
            r.evidence_count,
            r.imputed_anchor,
            r.occupancy_anchor_at_1u_m,
            r.occupancy_predicted_at_1u_m,
            r.anchor_delta_g_kj_mol,
            r.predicted_delta_g_kj_mol,
            r.qed_floor_total_kj_mol,
            r.residual_modeled_total_kj_mol
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    let payload = json!({
        "meta": {
            "lane": "antibiotic_resistance_beta_lactamase_panel",
            "temperature_k": temperature_k,
            "notes": [
                "Reduced-order mechanistic ranking lane; not clinical guidance",
                "Anchors are ChEMBL snapshot priors with assay/filter uncertainty"
            ]
        },
        "summary": {
            "pair_count": panel.rows.len(),
            "mean_abs_log10_error_pred_vs_anchor": panel.mean_abs_log10_error,
            "ndm_max_predicted_occupancy_at_1uM": panel.ndm_max_predicted_occupancy_at_1u_m
        },
        "best_by_enzyme": panel.best_by_enzyme.iter().map(|b| json!({
            "enzyme": b.enzyme_name,
            "anchor_winner": b.by_anchor_inhibitor,
            "anchor_winner_nM": b.by_anchor_nanomolar,
            "predicted_winner": b.by_predicted_inhibitor,
            "predicted_winner_nM": b.by_predicted_nanomolar,
            "winner_match": b.predicted_match_anchor_winner
        })).collect::<Vec<_>>(),
        "rows": panel.rows.iter().map(|r| json!({
            "enzyme": r.enzyme_name,
            "enzyme_chembl_hint": r.enzyme_chembl_hint,
            "enzyme_class": enzyme_class_name(r.enzyme_class),
            "inhibitor": r.inhibitor_name,
            "inhibitor_chembl_id": r.inhibitor_chembl_id,
            "scaffold": scaffold_name(r.scaffold),
            "evidence_count": r.evidence_count,
            "imputed_anchor": r.imputed_anchor,
            "anchor_nM": r.anchor_nanomolar,
            "predicted_nM": r.predicted_nanomolar,
            "log10_error_pred_vs_anchor": r.log10_error_pred_vs_anchor,
            "occupancy_anchor_at_1uM": r.occupancy_anchor_at_1u_m,
            "occupancy_predicted_at_1uM": r.occupancy_predicted_at_1u_m,
            "anchor_delta_g_kj_mol": r.anchor_delta_g_kj_mol,
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
        "antibiotic_resistance_panel: pairs={} mean_abs_log10_err={:.3} ndm_max_occ_1uM={:.3}",
        panel.rows.len(),
        panel.mean_abs_log10_error,
        panel.ndm_max_predicted_occupancy_at_1u_m
    );
}
