//! Direct simulation/application run for candidate:
//! `macrocycle_A__c20nM__buf3`.
//!
//! This is a reduced-order disease-trajectory simulation over expected values.

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
struct CoursePoint {
    month: u32,
    drive_index: f64,
    relapse_probability: f64,
    cumulative_relapses: f64,
    lesion_index: f64,
    disability_index: f64,
}

#[derive(Clone, Copy, Debug)]
struct CourseSummary {
    label: &'static str,
    mean_drive_index: f64,
    annualized_relapse_rate: f64,
    cumulative_relapses: f64,
    final_lesion_index: f64,
    final_disability_index: f64,
}

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn simulate_course(label: &'static str, base_drive_index: f64, p: SimParams) -> (CourseSummary, Vec<CoursePoint>) {
    let months = (p.years.max(0.25) * 12.0).round() as u32;
    let mut lesion = 1.0_f64;
    let mut cum_relapses = 0.0_f64;
    let mut points = Vec::with_capacity(months as usize + 1);
    let mut drive_sum = 0.0_f64;

    for m in 0..=months {
        let t = m as f64 / 12.0;
        let seasonal = 1.0 + p.seasonality_amp * (2.0 * PI * t).sin();
        let micro = 1.0 + 0.08 * (2.0 * PI * t * 2.0 + 0.7).cos();
        let drive = (base_drive_index.max(0.0) * seasonal * micro).clamp(0.0, 1.0);
        drive_sum += drive;

        let monthly_rate = (p.base_relapse_rate_per_year / 12.0) * (0.35 + 2.2 * drive);
        let relapse_prob = (1.0 - (-monthly_rate).exp()).clamp(0.0, 1.0);
        cum_relapses += relapse_prob;

        // Lesion burden: inflammatory growth + relapse spikes - endogenous repair.
        let growth = p.lesion_growth_coeff * drive + p.relapse_lesion_impact * relapse_prob;
        // Repair scales with present lesion burden; this avoids unrealistic
        // full clearance under mild immunomodulation in a reduced-order model.
        let repair = p.repair_rate * lesion * (1.0 - 0.45 * drive);
        lesion = (lesion + growth - repair).max(0.0);

        // Simple bounded disability proxy (0..1).
        let disability = (1.0 - (-lesion / 9.5).exp()).clamp(0.0, 1.0);
        points.push(CoursePoint {
            month: m,
            drive_index: drive,
            relapse_probability: relapse_prob,
            cumulative_relapses: cum_relapses,
            lesion_index: lesion,
            disability_index: disability,
        });
    }

    let mean_drive = drive_sum / (months as f64 + 1.0);
    let annualized_relapse_rate = cum_relapses / p.years.max(1.0e-9);
    let last = points.last().copied().expect("non-empty trajectory");
    (
        CourseSummary {
            label,
            mean_drive_index: mean_drive,
            annualized_relapse_rate,
            cumulative_relapses: cum_relapses,
            final_lesion_index: last.lesion_index,
            final_disability_index: last.disability_index,
        },
        points,
    )
}

fn summary_json(s: CourseSummary) -> serde_json::Value {
    json!({
        "label": s.label,
        "mean_drive_index": s.mean_drive_index,
        "annualized_relapse_rate": s.annualized_relapse_rate,
        "cumulative_relapses": s.cumulative_relapses,
        "final_lesion_index": s.final_lesion_index,
        "final_disability_index": s.final_disability_index
    })
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

    let candidate_label = std::env::var("GUTOE_MS_CANDIDATE_LABEL")
        .unwrap_or_else(|_| "macrocycle_A__c20nM__buf3".to_string());
    let candidate = TargetedBlockerCandidateInput {
        label: Box::leak(candidate_label.into_boxed_str()),
        concentration_nanomolar: env_f64("GUTOE_MS_CANDIDATE_CONC_NM", 20.0),
        target_ki_nanomolar: env_f64("GUTOE_MS_CANDIDATE_TARGET_KI_NM", 3.0),
        off_target_ki_nanomolar: env_f64("GUTOE_MS_CANDIDATE_OFFTARGET_KI_NM", 120.0),
        max_energy_shift_kj_mol: env_f64("GUTOE_MS_CANDIDATE_MAX_SHIFT_KJ_MOL", 3.0),
        safety_buffer_kj_mol: env_f64("GUTOE_MS_CANDIDATE_SAFETY_BUFFER_KJ_MOL", 0.3),
    };
    let blocker = evaluate_targeted_blocker_candidate(mimicry.activation_excess_kj_mol, candidate);

    // Convert achieved shift back into post-blocker risk drive.
    let off_target_penalty = 0.15 * blocker.off_target_occupancy_fraction;
    let transduction_efficiency = env_f64("GUTOE_MS_BLOCKER_SHIFT_EFFICIENCY", 0.30).clamp(0.0, 1.0);
    let effective_shift = blocker.achieved_energy_shift_kj_mol * transduction_efficiency;
    let activation_after = (mimicry.activation_excess_kj_mol - effective_shift + off_target_penalty)
        .max(0.0);
    let activation_score_after = (activation_after / (activation_after + 2.0)).clamp(0.0, 1.0);
    let blocker_only_drive = mimicry.overlap_score * activation_score_after;

    let standard_factor = if baseline_drive > 0.0 {
        (standard_drive / baseline_drive).clamp(0.0, 1.0)
    } else {
        1.0
    };
    let combo_drive = blocker_only_drive * standard_factor;

    let (baseline_summary, baseline_points) = simulate_course("baseline", baseline_drive, sim);
    let (standard_summary, standard_points) = simulate_course("standard_therapy_proxy", standard_drive, sim);
    let (blocker_summary, blocker_points) = simulate_course("macrocycle_candidate_only", blocker_only_drive, sim);
    let (combo_summary, combo_points) = simulate_course("macrocycle_plus_standard", combo_drive, sim);

    let lesion_reduction_blocker =
        (1.0 - blocker_summary.final_lesion_index / baseline_summary.final_lesion_index.max(1.0e-9)).clamp(-5.0, 1.0);
    let lesion_reduction_combo =
        (1.0 - combo_summary.final_lesion_index / baseline_summary.final_lesion_index.max(1.0e-9)).clamp(-5.0, 1.0);
    let relapse_reduction_blocker = (1.0
        - blocker_summary.annualized_relapse_rate / baseline_summary.annualized_relapse_rate.max(1.0e-9))
        .clamp(-5.0, 1.0);
    let relapse_reduction_combo = (1.0
        - combo_summary.annualized_relapse_rate / baseline_summary.annualized_relapse_rate.max(1.0e-9))
        .clamp(-5.0, 1.0);

    let out_dir = std::env::var("GUTOE_MS_MACROCYCLE_SIM_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ms_macrocycle_application".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);
    let txt_path = out.join("ms_macrocycle_application_report.txt");
    let json_path = out.join("ms_macrocycle_application_report.json");
    let summary_csv_path = out.join("ms_macrocycle_application_summary.csv");
    let traj_csv_path = out.join("ms_macrocycle_application_trajectory.csv");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ms_macrocycle_application]").expect("write");
    writeln!(txt, "candidate = {}", candidate.label).expect("write");
    writeln!(txt, "years = {:.3}", sim.years).expect("write");
    writeln!(txt, "baseline_drive = {:.9}", baseline_drive).expect("write");
    writeln!(txt, "standard_drive = {:.9}", standard_drive).expect("write");
    writeln!(txt, "blocker_only_drive = {:.9}", blocker_only_drive).expect("write");
    writeln!(txt, "combo_drive = {:.9}", combo_drive).expect("write");
    writeln!(txt, "blocker_target_occ = {:.9}", blocker.target_occupancy_fraction).expect("write");
    writeln!(txt, "blocker_off_target_occ = {:.9}", blocker.off_target_occupancy_fraction).expect("write");
    writeln!(txt, "blocker_margin_kj_mol = {:.9}", blocker.efficacy_margin_kj_mol).expect("write");
    writeln!(txt, "blocker_feasible = {}", blocker.feasible).expect("write");
    writeln!(txt, "lesion_reduction_blocker = {:.9}", lesion_reduction_blocker).expect("write");
    writeln!(txt, "lesion_reduction_combo = {:.9}", lesion_reduction_combo).expect("write");
    writeln!(txt, "relapse_reduction_blocker = {:.9}", relapse_reduction_blocker).expect("write");
    writeln!(txt, "relapse_reduction_combo = {:.9}", relapse_reduction_combo).expect("write");

    let mut summary_csv = String::from(
        "scenario,mean_drive_index,annualized_relapse_rate,cumulative_relapses,final_lesion_index,final_disability_index\n",
    );
    for s in [baseline_summary, standard_summary, blocker_summary, combo_summary] {
        summary_csv.push_str(&format!(
            "{},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            s.label,
            s.mean_drive_index,
            s.annualized_relapse_rate,
            s.cumulative_relapses,
            s.final_lesion_index,
            s.final_disability_index
        ));
    }
    fs::write(&summary_csv_path, summary_csv).expect("write summary csv");

    let mut traj_csv = String::from(
        "month,baseline_lesion,baseline_disability,baseline_relapses,standard_lesion,standard_disability,standard_relapses,blocker_lesion,blocker_disability,blocker_relapses,combo_lesion,combo_disability,combo_relapses\n",
    );
    let len = baseline_points
        .len()
        .min(standard_points.len())
        .min(blocker_points.len())
        .min(combo_points.len());
    for i in 0..len {
        let b = baseline_points[i];
        let s = standard_points[i];
        let m = blocker_points[i];
        let c = combo_points[i];
        traj_csv.push_str(&format!(
            "{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            b.month,
            b.lesion_index,
            b.disability_index,
            b.cumulative_relapses,
            s.lesion_index,
            s.disability_index,
            s.cumulative_relapses,
            m.lesion_index,
            m.disability_index,
            m.cumulative_relapses,
            c.lesion_index,
            c.disability_index,
            c.cumulative_relapses
        ));
    }
    fs::write(&traj_csv_path, traj_csv).expect("write trajectory csv");

    let payload = json!({
        "meta": {
            "lane": "ms_macrocycle_application_simulation",
            "candidate": candidate.label,
            "note": "reduced-order expected-value simulation, not clinical guidance"
        },
        "candidate": {
            "label": candidate.label,
            "concentration_nM": candidate.concentration_nanomolar,
            "target_ki_nM": candidate.target_ki_nanomolar,
            "off_target_ki_nM": candidate.off_target_ki_nanomolar,
            "max_energy_shift_kj_mol": candidate.max_energy_shift_kj_mol,
            "safety_buffer_kj_mol": candidate.safety_buffer_kj_mol
        },
        "mimicry": {
            "overlap_score": mimicry.overlap_score,
            "activation_excess_kj_mol": mimicry.activation_excess_kj_mol,
            "baseline_drive": baseline_drive
        },
        "blocker_effect": {
            "target_occupancy": blocker.target_occupancy_fraction,
            "off_target_occupancy": blocker.off_target_occupancy_fraction,
            "achieved_shift_kj_mol": blocker.achieved_energy_shift_kj_mol,
            "effective_shift_kj_mol": effective_shift,
            "transduction_efficiency": transduction_efficiency,
            "required_shift_kj_mol": blocker.required_energy_shift_kj_mol,
            "efficacy_margin_kj_mol": blocker.efficacy_margin_kj_mol,
            "feasible": blocker.feasible,
            "drive_after_blocker": blocker_only_drive
        },
        "summary": {
            "baseline": summary_json(baseline_summary),
            "standard_therapy_proxy": summary_json(standard_summary),
            "macrocycle_candidate_only": summary_json(blocker_summary),
            "macrocycle_plus_standard": summary_json(combo_summary),
            "lesion_reduction_blocker_fraction": lesion_reduction_blocker,
            "lesion_reduction_combo_fraction": lesion_reduction_combo,
            "relapse_reduction_blocker_fraction": relapse_reduction_blocker,
            "relapse_reduction_combo_fraction": relapse_reduction_combo
        }
    });
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("wrote {}", summary_csv_path.display());
    println!("wrote {}", traj_csv_path.display());
    println!(
        "ms_macrocycle_application: candidate={} feasible={} lesion_reduction_blocker={:.3} lesion_reduction_combo={:.3} relapse_reduction_blocker={:.3} relapse_reduction_combo={:.3}",
        candidate.label,
        blocker.feasible,
        lesion_reduction_blocker,
        lesion_reduction_combo,
        relapse_reduction_blocker,
        relapse_reduction_combo
    );
}
