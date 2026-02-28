//! File-ingestion phage-host matching report.
//!
//! Inputs:
//! - `GUTOE_PHAGE_STRAINS_FILE` (CSV or JSON)
//! - `GUTOE_PHAGE_LIBRARY_FILE` (CSV or JSON)
//!
//! Optional:
//! - `GUTOE_PHAGE_MATCH_TEMP_K` (default 310.15)
//! - `GUTOE_PHAGE_TOP_N_PER_STRAIN` (default 3)
//! - `GUTOE_PHAGE_MATCH_OUT` (default /tmp/bh_renders/phage_host_matching_ingest)
//! - `GUTOE_PHAGE_WRITE_EXAMPLES=1` to write template input files in output dir
//! - `GUTOE_PHAGE_USE_DEFAULTS=1` to run without input files

use gutoe_physics::{
    default_phage_host_matching_panel, default_phage_specs, default_strain_specs,
    evaluate_phage_host_matching_panel, load_phages_from_path, load_strains_from_path,
    PhageHostPairResult, PhageMatchingCoefficients,
};
use serde_json::json;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process;

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn env_bool(key: &str, default: bool) -> bool {
    std::env::var(key)
        .ok()
        .map(|v| {
            let k = v.trim().to_ascii_lowercase();
            k == "1" || k == "true" || k == "yes" || k == "y"
        })
        .unwrap_or(default)
}

fn receptor_name(r: gutoe_physics::ReceptorKind) -> &'static str {
    match r {
        gutoe_physics::ReceptorKind::LamB => "LamB",
        gutoe_physics::ReceptorKind::OmpK35 => "OmpK35",
        gutoe_physics::ReceptorKind::OmpK36 => "OmpK36",
        gutoe_physics::ReceptorKind::FhuA => "FhuA",
        gutoe_physics::ReceptorKind::LpsCore => "LpsCore",
        gutoe_physics::ReceptorKind::TypeIvPilus => "TypeIvPilus",
    }
}

fn write_example_templates(out: &Path) -> Result<(), String> {
    let strains = default_strain_specs();
    let phages = default_phage_specs();

    let strains_json = out.join("example_strains.json");
    let phages_json = out.join("example_phages.json");
    let strains_csv = out.join("example_strains.csv");
    let phages_csv = out.join("example_phages.csv");

    let strains_payload = json!({
        "strains": strains.iter().map(|s| json!({
            "name": s.name,
            "species": s.species,
            "resistance_marker": s.resistance_marker,
            "receptors": {
                "lamb": s.receptor_profile.lamb,
                "ompk35": s.receptor_profile.ompk35,
                "ompk36": s.receptor_profile.ompk36,
                "fhua": s.receptor_profile.fhua,
                "lps_core": s.receptor_profile.lps_core,
                "type_iv_pilus": s.receptor_profile.type_iv_pilus
            }
        })).collect::<Vec<_>>()
    });
    fs::write(
        &strains_json,
        serde_json::to_string_pretty(&strains_payload).map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;

    let phages_payload = json!({
        "phages": phages.iter().map(|p| json!({
            "name": p.name,
            "family": p.family,
            "primary_receptor": receptor_name(p.primary_receptor),
            "secondary_receptor": p.secondary_receptor.map(receptor_name),
            "secondary_weight": p.secondary_weight,
            "ionic_contact_count": p.ionic_contact_count,
            "hbond_contact_count": p.hbond_contact_count,
            "hydrophobic_area_a2": p.hydrophobic_area_a2,
            "conformational_entropy_penalty": p.conformational_entropy_penalty,
            "host_takeover_efficiency": p.host_takeover_efficiency
        })).collect::<Vec<_>>()
    });
    fs::write(
        &phages_json,
        serde_json::to_string_pretty(&phages_payload).map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;

    let mut s_csv = String::from(
        "name,species,resistance_marker,lamb,ompk35,ompk36,fhua,lps_core,type_iv_pilus\n",
    );
    for s in &strains {
        s_csv.push_str(&format!(
            "{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
            s.name,
            s.species,
            s.resistance_marker,
            s.receptor_profile.lamb,
            s.receptor_profile.ompk35,
            s.receptor_profile.ompk36,
            s.receptor_profile.fhua,
            s.receptor_profile.lps_core,
            s.receptor_profile.type_iv_pilus
        ));
    }
    fs::write(&strains_csv, s_csv).map_err(|e| e.to_string())?;

    let mut p_csv = String::from(
        "name,family,primary_receptor,secondary_receptor,secondary_weight,ionic_contact_count,hbond_contact_count,hydrophobic_area_a2,conformational_entropy_penalty,host_takeover_efficiency\n",
    );
    for p in &phages {
        p_csv.push_str(&format!(
            "{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
            p.name,
            p.family,
            receptor_name(p.primary_receptor),
            p.secondary_receptor.map(receptor_name).unwrap_or(""),
            p.secondary_weight,
            p.ionic_contact_count,
            p.hbond_contact_count,
            p.hydrophobic_area_a2,
            p.conformational_entropy_penalty,
            p.host_takeover_efficiency
        ));
    }
    fs::write(&phages_csv, p_csv).map_err(|e| e.to_string())?;

    Ok(())
}

fn main() {
    let out_dir = std::env::var("GUTOE_PHAGE_MATCH_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/phage_host_matching_ingest".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    if env_bool("GUTOE_PHAGE_WRITE_EXAMPLES", false) {
        if let Err(e) = write_example_templates(&out) {
            eprintln!("failed to write example templates: {e}");
            process::exit(2);
        }
    }

    let temperature_k = env_f64("GUTOE_PHAGE_MATCH_TEMP_K", 310.15);
    let top_n = env_usize("GUTOE_PHAGE_TOP_N_PER_STRAIN", 3).max(1);

    let strains_path = std::env::var("GUTOE_PHAGE_STRAINS_FILE").ok();
    let phages_path = std::env::var("GUTOE_PHAGE_LIBRARY_FILE").ok();
    let use_defaults = env_bool("GUTOE_PHAGE_USE_DEFAULTS", false);

    let (strains, phages, provenance) = match (strains_path, phages_path, use_defaults) {
        (Some(sp), Some(pp), _) => {
            let strains = match load_strains_from_path(&sp) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("failed loading strains file {}: {}", sp, e);
                    process::exit(2);
                }
            };
            let phages = match load_phages_from_path(&pp) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("failed loading phage library file {}: {}", pp, e);
                    process::exit(2);
                }
            };
            (strains, phages, format!("files: strains={} phages={}", sp, pp))
        }
        (_, _, true) => {
            let panel = default_phage_host_matching_panel(temperature_k);
            let strains = panel
                .best_by_strain
                .iter()
                .map(|b| b.strain_name.clone())
                .collect::<Vec<_>>();
            eprintln!(
                "using defaults (no input files provided); strains from default panel include: {}",
                strains.join(", ")
            );
            (
                default_strain_specs(),
                default_phage_specs(),
                "defaults".to_string(),
            )
        }
        _ => {
            eprintln!("missing input files.");
            eprintln!("set both GUTOE_PHAGE_STRAINS_FILE and GUTOE_PHAGE_LIBRARY_FILE");
            eprintln!("or set GUTOE_PHAGE_USE_DEFAULTS=1 for built-in example set.");
            eprintln!(
                "set GUTOE_PHAGE_WRITE_EXAMPLES=1 to write templates to {}",
                out.display()
            );
            process::exit(2);
        }
    };

    if strains.is_empty() || phages.is_empty() {
        eprintln!(
            "empty input after parsing (strains={}, phages={})",
            strains.len(),
            phages.len()
        );
        process::exit(2);
    }

    let mut panel = evaluate_phage_host_matching_panel(
        &strains,
        &phages,
        temperature_k,
        PhageMatchingCoefficients::default(),
    );

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

    let mut top_by_strain = BTreeMap::<String, Vec<PhageHostPairResult>>::new();
    for row in &panel.rows {
        let entry = top_by_strain.entry(row.strain_name.clone()).or_default();
        if entry.len() < top_n {
            entry.push(row.clone());
        }
    }

    let txt_path = out.join("phage_host_matching_ingest_report.txt");
    let csv_path = out.join("phage_host_matching_ingest_report.csv");
    let ranked_csv_path = out.join("phage_host_ranked_candidates.csv");
    let json_path = out.join("phage_host_matching_ingest_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[phage_host_matching_ingest]").expect("write");
    writeln!(txt, "provenance = {}", provenance).expect("write");
    writeln!(txt, "temperature_k = {:.6}", temperature_k).expect("write");
    writeln!(txt, "input_strain_count = {}", strains.len()).expect("write");
    writeln!(txt, "input_phage_count = {}", phages.len()).expect("write");
    writeln!(txt, "pair_count = {}", panel.rows.len()).expect("write");
    writeln!(txt, "top_n_per_strain = {}", top_n).expect("write");
    writeln!(txt, "mean_best_lysis_score = {:.9}", panel.mean_best_lysis_score).expect("write");
    writeln!(
        txt,
        "resistance_independence_probe_abs_delta = {:.12e}",
        panel.resistance_independence_probe_abs_delta
    )
    .expect("write");
    for b in &panel.best_by_strain {
        writeln!(
            txt,
            "best[{}|{}]: phage={} lysis={:.6} kd_nM={:.6}",
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

    let mut ranked_csv =
        String::from("strain,resistance_marker,rank,phage,predicted_kd_nM,lysis_potential_score\n");
    for (strain, rows) in &top_by_strain {
        for (i, r) in rows.iter().enumerate() {
            ranked_csv.push_str(&format!(
                "{},{},{},{},{:.6},{:.9}\n",
                strain,
                r.resistance_marker,
                i + 1,
                r.phage_name,
                r.predicted_kd_nanomolar,
                r.lysis_potential_score
            ));
        }
    }
    fs::write(&ranked_csv_path, ranked_csv).expect("write ranked csv");

    let payload = json!({
        "meta": {
            "lane": "phage_host_matching_ingest",
            "provenance": provenance,
            "temperature_k": temperature_k,
            "top_n_per_strain": top_n,
            "notes": [
                "File-driven phage-host ranking lane",
                "Bypass path: beta-lactamase marker is metadata, not score input",
                "Simulation ranking artifact, not clinical guidance"
            ]
        },
        "summary": {
            "input_strain_count": strains.len(),
            "input_phage_count": phages.len(),
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
        "top_candidates_by_strain": top_by_strain.iter().map(|(strain, rows)| {
            json!({
                "strain": strain,
                "candidates": rows.iter().enumerate().map(|(i, r)| json!({
                    "rank": i + 1,
                    "phage": r.phage_name,
                    "predicted_kd_nM": r.predicted_kd_nanomolar,
                    "lysis_potential_score": r.lysis_potential_score,
                    "resistance_marker": r.resistance_marker
                })).collect::<Vec<_>>()
            })
        }).collect::<Vec<_>>(),
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
            "qed_floor_total_kj_mol": r.qed_floor_total_kj_mol,
            "residual_modeled_total_kj_mol": r.residual_modeled_total_kj_mol
        })).collect::<Vec<_>>()
    });
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("serialize"))
        .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", ranked_csv_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "phage_host_matching_ingest: strains={} phages={} pairs={} mean_best_lysis={:.3}",
        strains.len(),
        phages.len(),
        panel.rows.len(),
        panel.mean_best_lysis_score
    );
}
