//! MS targeted blocker candidate search.
//!
//! Scans reduced-order blocker parameter space and ranks feasible candidates.
//! Output is a computational shortlist for further in-silico refinement.

use gutoe_physics::{
    default_ms_mimicry_input, evaluate_molecular_mimicry, evaluate_targeted_blocker_candidate,
    TargetedBlockerCandidateInput,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

#[derive(Clone, Copy, Debug)]
struct Archetype {
    name: &'static str,
    target_ki_nanomolar: f64,
    off_target_ki_nanomolar: f64,
    max_energy_shift_kj_mol: f64,
}

fn archetypes() -> Vec<Archetype> {
    vec![
        Archetype {
            name: "peptidomimetic_A",
            target_ki_nanomolar: 6.0,
            off_target_ki_nanomolar: 180.0,
            max_energy_shift_kj_mol: 2.2,
        },
        Archetype {
            name: "peptidomimetic_B",
            target_ki_nanomolar: 10.0,
            off_target_ki_nanomolar: 250.0,
            max_energy_shift_kj_mol: 2.4,
        },
        Archetype {
            name: "macrocycle_A",
            target_ki_nanomolar: 3.0,
            off_target_ki_nanomolar: 120.0,
            max_energy_shift_kj_mol: 3.0,
        },
        Archetype {
            name: "macrocycle_B",
            target_ki_nanomolar: 8.0,
            off_target_ki_nanomolar: 320.0,
            max_energy_shift_kj_mol: 2.8,
        },
        Archetype {
            name: "smallmol_A",
            target_ki_nanomolar: 15.0,
            off_target_ki_nanomolar: 220.0,
            max_energy_shift_kj_mol: 1.8,
        },
        Archetype {
            name: "smallmol_B",
            target_ki_nanomolar: 20.0,
            off_target_ki_nanomolar: 400.0,
            max_energy_shift_kj_mol: 1.6,
        },
    ]
}

fn main() {
    let mimicry = evaluate_molecular_mimicry(default_ms_mimicry_input());
    let activation_excess = mimicry.activation_excess_kj_mol;

    let concentrations = [3.0, 10.0, 20.0, 30.0, 50.0];
    let safety_buffers = [0.3, 0.5, 0.8];
    let mut candidates = Vec::new();

    for arch in archetypes() {
        for &c in &concentrations {
            for &buf in &safety_buffers {
                let label = format!("{}__c{}nM__buf{}", arch.name, c as i32, (buf * 10.0) as i32);
                let c_input = TargetedBlockerCandidateInput {
                    label: Box::leak(label.into_boxed_str()),
                    concentration_nanomolar: c,
                    target_ki_nanomolar: arch.target_ki_nanomolar,
                    off_target_ki_nanomolar: arch.off_target_ki_nanomolar,
                    max_energy_shift_kj_mol: arch.max_energy_shift_kj_mol,
                    safety_buffer_kj_mol: buf,
                };
                let s = evaluate_targeted_blocker_candidate(activation_excess, c_input);
                candidates.push(s);
            }
        }
    }

    candidates.sort_by(|a, b| {
        b.candidate_score
            .partial_cmp(&a.candidate_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let feasible: Vec<_> = candidates.iter().copied().filter(|c| c.feasible).collect();
    let top_k = feasible.len().min(12);
    let top = &feasible[..top_k];

    let out_dir = std::env::var("GUTOE_MS_BLOCKER_SEARCH_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ms_targeted_blocker_search".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("ms_targeted_blocker_search.txt");
    let csv_path = out.join("ms_targeted_blocker_search.csv");
    let json_path = out.join("ms_targeted_blocker_search.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ms_targeted_blocker_search]").expect("write");
    writeln!(txt, "activation_excess_kj_mol = {:.9}", activation_excess).expect("write");
    writeln!(txt, "screened_candidates = {}", candidates.len()).expect("write");
    writeln!(txt, "feasible_candidates = {}", feasible.len()).expect("write");
    if let Some(best) = top.first() {
        writeln!(txt, "top_candidate = {}", best.label).expect("write");
        writeln!(txt, "top_candidate_score = {:.9}", best.candidate_score).expect("write");
    }

    let mut csv = String::from(
        "label,score,feasible,concentration_nM,target_ki_nM,off_target_ki_nM,max_energy_shift_kj_mol,target_occupancy,off_target_occupancy,selectivity_ratio,required_shift_kj_mol,achieved_shift_kj_mol,efficacy_margin_kj_mol\n",
    );
    for c in top {
        csv.push_str(&format!(
            "{},{:.9},{},{:.6},{:.6},{:.6},{:.6},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            c.label,
            c.candidate_score,
            c.feasible,
            c.concentration_nanomolar,
            c.target_ki_nanomolar,
            c.off_target_ki_nanomolar,
            c.max_energy_shift_kj_mol,
            c.target_occupancy_fraction,
            c.off_target_occupancy_fraction,
            c.selectivity_ratio,
            c.required_energy_shift_kj_mol,
            c.achieved_energy_shift_kj_mol,
            c.efficacy_margin_kj_mol
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    let payload = json!({
        "meta": {
            "lane": "ms_targeted_blocker_search",
            "activation_excess_kj_mol": activation_excess
        },
        "summary": {
            "screened_candidates": candidates.len(),
            "feasible_candidates": feasible.len(),
            "top_candidate": top.first().map(|c| c.label)
        },
        "top_candidates": top.iter().map(|c| json!({
            "label": c.label,
            "score": c.candidate_score,
            "concentration_nM": c.concentration_nanomolar,
            "target_ki_nM": c.target_ki_nanomolar,
            "off_target_ki_nM": c.off_target_ki_nanomolar,
            "max_energy_shift_kj_mol": c.max_energy_shift_kj_mol,
            "target_occupancy": c.target_occupancy_fraction,
            "off_target_occupancy": c.off_target_occupancy_fraction,
            "selectivity_ratio": c.selectivity_ratio,
            "required_shift_kj_mol": c.required_energy_shift_kj_mol,
            "achieved_shift_kj_mol": c.achieved_energy_shift_kj_mol,
            "efficacy_margin_kj_mol": c.efficacy_margin_kj_mol,
            "feasible": c.feasible
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
    if let Some(best) = top.first() {
        println!(
            "ms_targeted_blocker_search: top={} score={:.3} target_occ={:.3} off_occ={:.3} margin={:.3} kJ/mol",
            best.label,
            best.candidate_score,
            best.target_occupancy_fraction,
            best.off_target_occupancy_fraction,
            best.efficacy_margin_kj_mol
        );
    } else {
        println!("ms_targeted_blocker_search: no feasible candidates found");
    }
}
