//! Finds the closest molecule proxy match to the validated
//! `macrocycle_A__c20nM__buf3` MS application profile.
//!
//! This is a reduced-order in-silico ranking lane, not clinical guidance.

use gutoe_physics::{
    default_ms_mimicry_input, default_natalizumab_proxy, default_ocrelizumab_proxy,
    evaluate_molecular_mimicry, evaluate_targeted_blocker_candidate, evaluate_therapy_effect,
    TargetedBlockerCandidateInput,
};
use serde_json::json;
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

#[derive(Clone, Copy, Debug)]
struct SimParams {
    years: f64,
    base_relapse_rate_per_year: f64,
    lesion_growth_coeff: f64,
    relapse_lesion_impact: f64,
    repair_rate: f64,
    seasonality_amp: f64,
}

#[derive(Clone, Copy, Debug)]
struct CourseSummary {
    annualized_relapse_rate: f64,
    final_lesion_index: f64,
    final_disability_index: f64,
}

#[derive(Clone, Copy, Debug)]
struct ReferenceProfile {
    target_occupancy: f64,
    off_target_occupancy: f64,
    efficacy_margin_kj_mol: f64,
    lesion_reduction_combo: f64,
    relapse_reduction_combo: f64,
    disability_combo: f64,
}

#[derive(Clone, Copy, Debug)]
struct Archetype {
    name: &'static str,
    target_ki_nanomolar: f64,
    off_target_ki_nanomolar: f64,
    max_energy_shift_kj_mol: f64,
}

#[derive(Clone, Copy, Debug)]
struct RankedCandidate {
    label: &'static str,
    archetype: &'static str,
    concentration_nanomolar: f64,
    target_ki_nanomolar: f64,
    off_target_ki_nanomolar: f64,
    max_energy_shift_kj_mol: f64,
    safety_buffer_kj_mol: f64,
    target_occupancy: f64,
    off_target_occupancy: f64,
    efficacy_margin_kj_mol: f64,
    blocker_score: f64,
    lesion_reduction_combo: f64,
    relapse_reduction_combo: f64,
    disability_combo: f64,
    profile_distance: f64,
    profile_match_score: f64,
}

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn archetypes() -> Vec<Archetype> {
    vec![
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
            name: "macrocycle_C",
            target_ki_nanomolar: 4.0,
            off_target_ki_nanomolar: 180.0,
            max_energy_shift_kj_mol: 3.2,
        },
        Archetype {
            name: "macrocycle_D",
            target_ki_nanomolar: 2.6,
            off_target_ki_nanomolar: 105.0,
            max_energy_shift_kj_mol: 2.9,
        },
        Archetype {
            name: "macrocycle_E",
            target_ki_nanomolar: 3.4,
            off_target_ki_nanomolar: 210.0,
            max_energy_shift_kj_mol: 3.3,
        },
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
            name: "spiro_A",
            target_ki_nanomolar: 5.5,
            off_target_ki_nanomolar: 260.0,
            max_energy_shift_kj_mol: 2.9,
        },
        Archetype {
            name: "spiro_B",
            target_ki_nanomolar: 4.2,
            off_target_ki_nanomolar: 240.0,
            max_energy_shift_kj_mol: 3.1,
        },
        Archetype {
            name: "bicyclic_A",
            target_ki_nanomolar: 7.5,
            off_target_ki_nanomolar: 430.0,
            max_energy_shift_kj_mol: 2.6,
        },
        Archetype {
            name: "bicyclic_B",
            target_ki_nanomolar: 4.8,
            off_target_ki_nanomolar: 360.0,
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

fn simulate_course(base_drive_index: f64, p: SimParams) -> CourseSummary {
    let months = (p.years.max(0.25) * 12.0).round() as u32;
    let mut lesion = 1.0_f64;
    let mut cum_relapses = 0.0_f64;

    for m in 0..=months {
        let t = m as f64 / 12.0;
        let seasonal = 1.0 + p.seasonality_amp * (2.0 * PI * t).sin();
        let micro = 1.0 + 0.08 * (2.0 * PI * t * 2.0 + 0.7).cos();
        let drive = (base_drive_index.max(0.0) * seasonal * micro).clamp(0.0, 1.0);

        let monthly_rate = (p.base_relapse_rate_per_year / 12.0) * (0.35 + 2.2 * drive);
        let relapse_prob = (1.0 - (-monthly_rate).exp()).clamp(0.0, 1.0);
        cum_relapses += relapse_prob;

        let growth = p.lesion_growth_coeff * drive + p.relapse_lesion_impact * relapse_prob;
        let repair = p.repair_rate * lesion * (1.0 - 0.45 * drive);
        lesion = (lesion + growth - repair).max(0.0);
    }

    CourseSummary {
        annualized_relapse_rate: cum_relapses / p.years.max(1.0e-9),
        final_lesion_index: lesion,
        final_disability_index: (1.0 - (-lesion / 9.5).exp()).clamp(0.0, 1.0),
    }
}

fn build_reference_profile(
    activation_excess: f64,
    baseline_drive: f64,
    standard_factor: f64,
    sim: SimParams,
    transduction_efficiency: f64,
) -> (ReferenceProfile, CourseSummary) {
    let reference = TargetedBlockerCandidateInput {
        label: "macrocycle_A__c20nM__buf3",
        concentration_nanomolar: 20.0,
        target_ki_nanomolar: 3.0,
        off_target_ki_nanomolar: 120.0,
        max_energy_shift_kj_mol: 3.0,
        safety_buffer_kj_mol: 0.3,
    };
    let reference_score = evaluate_targeted_blocker_candidate(activation_excess, reference);

    let off_target_penalty = 0.15 * reference_score.off_target_occupancy_fraction;
    let effective_shift = reference_score.achieved_energy_shift_kj_mol * transduction_efficiency;
    let activation_after = (activation_excess - effective_shift + off_target_penalty).max(0.0);
    let activation_score_after = (activation_after / (activation_after + 2.0)).clamp(0.0, 1.0);

    let mimicry = evaluate_molecular_mimicry(default_ms_mimicry_input());
    let blocker_drive = mimicry.overlap_score * activation_score_after;
    let combo_drive = blocker_drive * standard_factor;

    let baseline_summary = simulate_course(baseline_drive, sim);
    let combo_summary = simulate_course(combo_drive, sim);

    let lesion_reduction_combo =
        (1.0 - combo_summary.final_lesion_index / baseline_summary.final_lesion_index.max(1.0e-9))
            .clamp(-5.0, 1.0);
    let relapse_reduction_combo = (1.0
        - combo_summary.annualized_relapse_rate / baseline_summary.annualized_relapse_rate.max(1.0e-9))
        .clamp(-5.0, 1.0);

    (
        ReferenceProfile {
            target_occupancy: reference_score.target_occupancy_fraction,
            off_target_occupancy: reference_score.off_target_occupancy_fraction,
            efficacy_margin_kj_mol: reference_score.efficacy_margin_kj_mol,
            lesion_reduction_combo,
            relapse_reduction_combo,
            disability_combo: combo_summary.final_disability_index,
        },
        baseline_summary,
    )
}

fn profile_distance(c: &RankedCandidate, r: ReferenceProfile) -> f64 {
    let z_occ = (c.target_occupancy - r.target_occupancy) / 0.10;
    let z_off = (c.off_target_occupancy - r.off_target_occupancy) / 0.08;
    let z_margin = (c.efficacy_margin_kj_mol - r.efficacy_margin_kj_mol) / 0.60;
    let z_lesion = (c.lesion_reduction_combo - r.lesion_reduction_combo) / 0.10;
    let z_relapse = (c.relapse_reduction_combo - r.relapse_reduction_combo) / 0.08;
    let z_disability = (c.disability_combo - r.disability_combo) / 0.02;

    (z_occ * z_occ
        + z_off * z_off
        + z_margin * z_margin
        + z_lesion * z_lesion
        + z_relapse * z_relapse
        + z_disability * z_disability)
        .sqrt()
}

fn main() {
    let sim = SimParams {
        years: env_f64("GUTOE_MS_SIM_YEARS", 10.0),
        base_relapse_rate_per_year: env_f64("GUTOE_MS_SIM_BASE_RELAPSE_PER_YEAR", 0.60),
        lesion_growth_coeff: env_f64("GUTOE_MS_SIM_LESION_GROWTH", 0.12),
        relapse_lesion_impact: env_f64("GUTOE_MS_SIM_RELAPSE_IMPACT", 0.25),
        repair_rate: env_f64("GUTOE_MS_SIM_REPAIR_RATE", 0.03),
        seasonality_amp: env_f64("GUTOE_MS_SIM_SEASONALITY_AMP", 0.12),
    };

    let mimicry = evaluate_molecular_mimicry(default_ms_mimicry_input());
    let baseline_drive = mimicry.misrecognition_risk_index;
    let standard = evaluate_therapy_effect(
        baseline_drive,
        default_ocrelizumab_proxy(),
        default_natalizumab_proxy(),
    );
    let standard_drive = standard.residual_drive_index;
    let standard_factor = if baseline_drive > 0.0 {
        (standard_drive / baseline_drive).clamp(0.0, 1.0)
    } else {
        1.0
    };

    let transduction_efficiency = env_f64("GUTOE_MS_BLOCKER_SHIFT_EFFICIENCY", 0.30).clamp(0.0, 1.0);
    let (reference, baseline_summary) = build_reference_profile(
        mimicry.activation_excess_kj_mol,
        baseline_drive,
        standard_factor,
        sim,
        transduction_efficiency,
    );

    let concentrations = [6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 35.0];
    let safety_buffers = [0.3, 0.4, 0.5, 0.6];
    let mut ranked = Vec::<RankedCandidate>::new();

    for arch in archetypes() {
        for &c in &concentrations {
            for &buf in &safety_buffers {
                let label = format!("{}__c{}nM__buf{}", arch.name, c as i32, (buf * 10.0) as i32);
                let input = TargetedBlockerCandidateInput {
                    label: Box::leak(label.into_boxed_str()),
                    concentration_nanomolar: c,
                    target_ki_nanomolar: arch.target_ki_nanomolar,
                    off_target_ki_nanomolar: arch.off_target_ki_nanomolar,
                    max_energy_shift_kj_mol: arch.max_energy_shift_kj_mol,
                    safety_buffer_kj_mol: buf,
                };
                let score = evaluate_targeted_blocker_candidate(mimicry.activation_excess_kj_mol, input);
                if !score.feasible {
                    continue;
                }

                let off_target_penalty = 0.15 * score.off_target_occupancy_fraction;
                let effective_shift = score.achieved_energy_shift_kj_mol * transduction_efficiency;
                let activation_after =
                    (mimicry.activation_excess_kj_mol - effective_shift + off_target_penalty).max(0.0);
                let activation_score_after = (activation_after / (activation_after + 2.0)).clamp(0.0, 1.0);
                let blocker_drive = mimicry.overlap_score * activation_score_after;
                let combo_drive = blocker_drive * standard_factor;

                let combo_summary = simulate_course(combo_drive, sim);
                let lesion_reduction_combo = (1.0
                    - combo_summary.final_lesion_index / baseline_summary.final_lesion_index.max(1.0e-9))
                .clamp(-5.0, 1.0);
                let relapse_reduction_combo = (1.0
                    - combo_summary.annualized_relapse_rate
                        / baseline_summary.annualized_relapse_rate.max(1.0e-9))
                .clamp(-5.0, 1.0);

                let mut row = RankedCandidate {
                    label: score.label,
                    archetype: arch.name,
                    concentration_nanomolar: score.concentration_nanomolar,
                    target_ki_nanomolar: score.target_ki_nanomolar,
                    off_target_ki_nanomolar: score.off_target_ki_nanomolar,
                    max_energy_shift_kj_mol: score.max_energy_shift_kj_mol,
                    safety_buffer_kj_mol: buf,
                    target_occupancy: score.target_occupancy_fraction,
                    off_target_occupancy: score.off_target_occupancy_fraction,
                    efficacy_margin_kj_mol: score.efficacy_margin_kj_mol,
                    blocker_score: score.candidate_score,
                    lesion_reduction_combo,
                    relapse_reduction_combo,
                    disability_combo: combo_summary.final_disability_index,
                    profile_distance: 0.0,
                    profile_match_score: 0.0,
                };
                row.profile_distance = profile_distance(&row, reference);
                row.profile_match_score = 1.0 / (1.0 + row.profile_distance);
                ranked.push(row);
            }
        }
    }

    ranked.sort_by(|a, b| {
        a.profile_distance
            .partial_cmp(&b.profile_distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let top_k = ranked.len().min(20);
    let top = &ranked[..top_k];
    let best = top.first().copied();
    let best_alt = top
        .iter()
        .copied()
        .find(|c| c.label != "macrocycle_A__c20nM__buf3");

    let out_dir = std::env::var("GUTOE_MS_MOLECULE_MATCH_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ms_molecule_match_search".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("ms_molecule_match_search.txt");
    let csv_path = out.join("ms_molecule_match_search.csv");
    let json_path = out.join("ms_molecule_match_search.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ms_molecule_match_search]").expect("write");
    writeln!(txt, "reference = macrocycle_A__c20nM__buf3").expect("write");
    writeln!(txt, "screened_feasible_candidates = {}", ranked.len()).expect("write");
    if let Some(b) = best {
        writeln!(txt, "best_match = {}", b.label).expect("write");
        writeln!(txt, "best_distance = {:.9}", b.profile_distance).expect("write");
    }
    if let Some(a) = best_alt {
        writeln!(txt, "best_alternative = {}", a.label).expect("write");
        writeln!(txt, "best_alternative_distance = {:.9}", a.profile_distance).expect("write");
    }

    let mut csv = String::from(
        "label,archetype,profile_distance,profile_match_score,blocker_score,concentration_nM,target_ki_nM,off_target_ki_nM,max_energy_shift_kj_mol,safety_buffer_kj_mol,target_occupancy,off_target_occupancy,efficacy_margin_kj_mol,lesion_reduction_combo,relapse_reduction_combo,disability_combo\n",
    );
    for r in top {
        csv.push_str(&format!(
            "{},{},{:.9},{:.9},{:.9},{:.3},{:.3},{:.3},{:.3},{:.3},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            r.label,
            r.archetype,
            r.profile_distance,
            r.profile_match_score,
            r.blocker_score,
            r.concentration_nanomolar,
            r.target_ki_nanomolar,
            r.off_target_ki_nanomolar,
            r.max_energy_shift_kj_mol,
            r.safety_buffer_kj_mol,
            r.target_occupancy,
            r.off_target_occupancy,
            r.efficacy_margin_kj_mol,
            r.lesion_reduction_combo,
            r.relapse_reduction_combo,
            r.disability_combo
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    let payload = json!({
        "meta": {
            "lane": "ms_molecule_match_search",
            "reference": "macrocycle_A__c20nM__buf3",
            "transduction_efficiency": transduction_efficiency
        },
        "reference_profile": {
            "target_occupancy": reference.target_occupancy,
            "off_target_occupancy": reference.off_target_occupancy,
            "efficacy_margin_kj_mol": reference.efficacy_margin_kj_mol,
            "lesion_reduction_combo": reference.lesion_reduction_combo,
            "relapse_reduction_combo": reference.relapse_reduction_combo,
            "disability_combo": reference.disability_combo
        },
        "summary": {
            "screened_feasible_candidates": ranked.len(),
            "best_match": best.map(|c| c.label),
            "best_alternative": best_alt.map(|c| c.label)
        },
        "top_matches": top.iter().map(|c| json!({
            "label": c.label,
            "archetype": c.archetype,
            "profile_distance": c.profile_distance,
            "profile_match_score": c.profile_match_score,
            "blocker_score": c.blocker_score,
            "concentration_nM": c.concentration_nanomolar,
            "target_ki_nM": c.target_ki_nanomolar,
            "off_target_ki_nM": c.off_target_ki_nanomolar,
            "max_energy_shift_kj_mol": c.max_energy_shift_kj_mol,
            "safety_buffer_kj_mol": c.safety_buffer_kj_mol,
            "target_occupancy": c.target_occupancy,
            "off_target_occupancy": c.off_target_occupancy,
            "efficacy_margin_kj_mol": c.efficacy_margin_kj_mol,
            "lesion_reduction_combo": c.lesion_reduction_combo,
            "relapse_reduction_combo": c.relapse_reduction_combo,
            "disability_combo": c.disability_combo
        })).collect::<Vec<_>>()
    });
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("serialize"))
        .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    if let Some(b) = best {
        println!(
            "ms_molecule_match_search: best={} distance={:.4} match_score={:.4}",
            b.label, b.profile_distance, b.profile_match_score
        );
    }
    if let Some(a) = best_alt {
        println!(
            "ms_molecule_match_search: best_alternative={} distance={:.4} match_score={:.4}",
            a.label, a.profile_distance, a.profile_match_score
        );
    }
}
