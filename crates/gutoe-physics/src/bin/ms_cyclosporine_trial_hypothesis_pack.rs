//! Cyclosporine MS trial-hypothesis pack (simulation-first).
//!
//! Produces a falsifiable protocol scaffold from:
//! - mechanistic MS lane (molecular mimicry + blocker transduction)
//! - PK bridge uncertainty
//! - safety gate probabilities
//!
//! This artifact is for hypothesis generation only, not clinical guidance.

use gutoe_physics::{
    default_cyclosporine_pk_bridge_input, default_cyclosporine_safety_gate_input,
    default_ms_mimicry_input, default_natalizumab_proxy, default_ocrelizumab_proxy,
    evaluate_cyclosporine_safety_gate, evaluate_molecular_mimicry,
    evaluate_targeted_blocker_candidate, evaluate_therapy_effect, simulate_cyclosporine_pk_bridge,
    summarize_cyclosporine_pk_bridge, TargetedBlockerCandidateInput,
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

fn poisson_sample_size_per_arm(
    lambda_control_per_year: f64,
    lambda_treatment_per_year: f64,
    followup_years: f64,
    z_alpha_two_sided: f64,
    z_power: f64,
) -> f64 {
    let lc = lambda_control_per_year.max(1.0e-9);
    let lt = lambda_treatment_per_year.max(1.0e-9);
    let t = followup_years.max(0.25);
    let delta = (lc - lt).abs().max(1.0e-9);
    ((z_alpha_two_sided + z_power).powi(2) * (lc + lt)) / (delta * delta * t)
}

fn main() {
    let cyclosporine = TargetedBlockerCandidateInput {
        label: "cyclosporine__c20nM__buf3",
        concentration_nanomolar: env_f64("GUTOE_MS_CANDIDATE_CONC_NM", 20.0),
        target_ki_nanomolar: env_f64("GUTOE_MS_CANDIDATE_TARGET_KI_NM", 2.64),
        off_target_ki_nanomolar: env_f64("GUTOE_MS_CANDIDATE_OFFTARGET_KI_NM", 200.0),
        max_energy_shift_kj_mol: env_f64("GUTOE_MS_CANDIDATE_MAX_SHIFT_KJ_MOL", 3.0),
        safety_buffer_kj_mol: env_f64("GUTOE_MS_CANDIDATE_SAFETY_BUFFER_KJ_MOL", 0.3),
    };

    let mimicry = evaluate_molecular_mimicry(default_ms_mimicry_input());
    let baseline_drive = mimicry.misrecognition_risk_index;
    let standard = evaluate_therapy_effect(
        baseline_drive,
        default_ocrelizumab_proxy(),
        default_natalizumab_proxy(),
    );
    let standard_drive = standard.residual_drive_index;

    let blocker = evaluate_targeted_blocker_candidate(mimicry.activation_excess_kj_mol, cyclosporine);
    let transduction_efficiency = env_f64("GUTOE_MS_BLOCKER_SHIFT_EFFICIENCY", 0.30).clamp(0.0, 1.0);
    let off_target_penalty = 0.15 * blocker.off_target_occupancy_fraction;
    let effective_shift = blocker.achieved_energy_shift_kj_mol * transduction_efficiency;
    let activation_after =
        (mimicry.activation_excess_kj_mol - effective_shift + off_target_penalty).max(0.0);
    let activation_score_after = (activation_after / (activation_after + 2.0)).clamp(0.0, 1.0);
    let blocker_drive = mimicry.overlap_score * activation_score_after;

    let standard_factor = if baseline_drive > 0.0 {
        (standard_drive / baseline_drive).clamp(0.0, 1.0)
    } else {
        1.0
    };
    let combo_drive = blocker_drive * standard_factor;

    let ten_year = SimParams {
        years: env_f64("GUTOE_MS_SIM_YEARS", 10.0),
        base_relapse_rate_per_year: env_f64("GUTOE_MS_SIM_BASE_RELAPSE_PER_YEAR", 0.60),
        lesion_growth_coeff: env_f64("GUTOE_MS_SIM_LESION_GROWTH", 0.12),
        relapse_lesion_impact: env_f64("GUTOE_MS_SIM_RELAPSE_IMPACT", 0.25),
        repair_rate: env_f64("GUTOE_MS_SIM_REPAIR_RATE", 0.03),
        seasonality_amp: env_f64("GUTOE_MS_SIM_SEASONALITY_AMP", 0.12),
    };
    let two_year = SimParams {
        years: env_f64("GUTOE_MS_TRIAL_HORIZON_YEARS", 2.0),
        ..ten_year
    };

    let baseline_10y = simulate_course(baseline_drive, ten_year);
    let standard_10y = simulate_course(standard_drive, ten_year);
    let combo_10y = simulate_course(combo_drive, ten_year);

    let baseline_2y = simulate_course(baseline_drive, two_year);
    let standard_2y = simulate_course(standard_drive, two_year);
    let combo_2y = simulate_course(combo_drive, two_year);

    let rr_combo_vs_standard_2y =
        combo_2y.annualized_relapse_rate / standard_2y.annualized_relapse_rate.max(1.0e-9);
    let rr_reduction_combo_vs_standard_2y = (1.0 - rr_combo_vs_standard_2y).clamp(-5.0, 1.0);

    let n_per_arm_2y = poisson_sample_size_per_arm(
        standard_2y.annualized_relapse_rate,
        combo_2y.annualized_relapse_rate,
        two_year.years,
        1.96,
        0.84,
    );

    let mut pk_input = default_cyclosporine_pk_bridge_input();
    pk_input.site_target_nanomolar = cyclosporine.concentration_nanomolar;
    pk_input.samples = env_usize("GUTOE_MS_PK_SAMPLES", pk_input.samples);
    pk_input.seed = env_u64("GUTOE_MS_PK_SEED", pk_input.seed);
    pk_input.blood_to_site_gain_median = env_f64("GUTOE_MS_PK_GAIN_MEDIAN", pk_input.blood_to_site_gain_median);
    pk_input.blood_to_site_gain_gsd = env_f64("GUTOE_MS_PK_GAIN_GSD", pk_input.blood_to_site_gain_gsd);
    let pk_ens = simulate_cyclosporine_pk_bridge(pk_input);
    let pk = summarize_cyclosporine_pk_bridge(&pk_ens);

    let mut safety_gate = default_cyclosporine_safety_gate_input();
    safety_gate.windows.target_zone_low_ng_ml =
        env_f64("GUTOE_MS_TARGET_ZONE_LOW_NG_ML", safety_gate.windows.target_zone_low_ng_ml);
    safety_gate.windows.target_zone_high_ng_ml =
        env_f64("GUTOE_MS_TARGET_ZONE_HIGH_NG_ML", safety_gate.windows.target_zone_high_ng_ml);
    safety_gate.windows.renal_caution_ng_ml =
        env_f64("GUTOE_MS_RENAL_CAUTION_NG_ML", safety_gate.windows.renal_caution_ng_ml);
    safety_gate.windows.renal_high_ng_ml =
        env_f64("GUTOE_MS_RENAL_HIGH_NG_ML", safety_gate.windows.renal_high_ng_ml);
    safety_gate.windows.neuro_caution_ng_ml =
        env_f64("GUTOE_MS_NEURO_CAUTION_NG_ML", safety_gate.windows.neuro_caution_ng_ml);
    safety_gate.max_prob_above_renal_caution = env_f64(
        "GUTOE_MS_MAX_P_ABOVE_RENAL_CAUTION",
        safety_gate.max_prob_above_renal_caution,
    );
    safety_gate.max_prob_above_renal_high =
        env_f64("GUTOE_MS_MAX_P_ABOVE_RENAL_HIGH", safety_gate.max_prob_above_renal_high);
    safety_gate.max_prob_above_neuro_caution = env_f64(
        "GUTOE_MS_MAX_P_ABOVE_NEURO_CAUTION",
        safety_gate.max_prob_above_neuro_caution,
    );
    safety_gate.min_prob_in_target_zone =
        env_f64("GUTOE_MS_MIN_P_IN_TARGET_ZONE", safety_gate.min_prob_in_target_zone);
    let safety = evaluate_cyclosporine_safety_gate(&pk_ens, safety_gate);

    let out_dir = std::env::var("GUTOE_MS_CYCLOSPORINE_TRIAL_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ms_cyclosporine_trial_hypothesis".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("ms_cyclosporine_trial_hypothesis.txt");
    let json_path = out.join("ms_cyclosporine_trial_hypothesis.json");
    let csv_path = out.join("ms_cyclosporine_trial_key_metrics.csv");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ms_cyclosporine_trial_hypothesis]").expect("write");
    writeln!(txt, "candidate = {}", cyclosporine.label).expect("write");
    writeln!(txt, "target_ki_nM = {:.6}", cyclosporine.target_ki_nanomolar).expect("write");
    writeln!(txt, "off_target_ki_nM = {:.6}", cyclosporine.off_target_ki_nanomolar).expect("write");
    writeln!(txt, "target_occupancy = {:.9}", blocker.target_occupancy_fraction).expect("write");
    writeln!(txt, "off_target_occupancy = {:.9}", blocker.off_target_occupancy_fraction).expect("write");
    writeln!(txt, "efficacy_margin_kj_mol = {:.9}", blocker.efficacy_margin_kj_mol).expect("write");
    writeln!(txt, "arr_standard_2y = {:.9}", standard_2y.annualized_relapse_rate).expect("write");
    writeln!(txt, "arr_combo_2y = {:.9}", combo_2y.annualized_relapse_rate).expect("write");
    writeln!(txt, "arr_rr_reduction_combo_vs_standard_2y = {:.9}", rr_reduction_combo_vs_standard_2y)
        .expect("write");
    writeln!(txt, "lesion_reduction_combo_vs_standard_10y = {:.9}",
        (1.0 - combo_10y.final_lesion_index / standard_10y.final_lesion_index.max(1.0e-9)).clamp(-5.0, 1.0)
    ).expect("write");
    writeln!(txt, "pk_target_window_ng_mL = [{:.3}, {:.3}]", pk.p25_ng_ml, pk.p75_ng_ml)
        .expect("write");
    writeln!(txt, "safety_gate_pass = {}", safety.overall_pass).expect("write");
    writeln!(txt, "poisson_n_per_arm_2y_80pct = {:.3}", n_per_arm_2y).expect("write");

    let csv = format!(
        "metric,value\narr_standard_2y,{:.9}\narr_combo_2y,{:.9}\nrr_reduction_combo_vs_standard_2y,{:.9}\nlesion_reduction_combo_vs_standard_10y,{:.9}\ndisability_standard_10y,{:.9}\ndisability_combo_10y,{:.9}\npk_p25_ng_mL,{:.9}\npk_p50_ng_mL,{:.9}\npk_p75_ng_mL,{:.9}\np_renal_caution,{:.9}\np_renal_high,{:.9}\np_neuro_caution,{:.9}\npoisson_n_per_arm_2y_80pct,{:.9}\n",
        standard_2y.annualized_relapse_rate,
        combo_2y.annualized_relapse_rate,
        rr_reduction_combo_vs_standard_2y,
        (1.0 - combo_10y.final_lesion_index / standard_10y.final_lesion_index.max(1.0e-9)).clamp(-5.0, 1.0),
        standard_10y.final_disability_index,
        combo_10y.final_disability_index,
        pk.p25_ng_ml,
        pk.p50_ng_ml,
        pk.p75_ng_ml,
        safety.prob_above_renal_caution,
        safety.prob_above_renal_high,
        safety.prob_above_neuro_caution,
        n_per_arm_2y,
    );
    fs::write(&csv_path, csv).expect("write csv");

    let payload = json!({
        "meta": {
            "lane": "ms_cyclosporine_trial_hypothesis_pack",
            "note": "hypothesis pack for prospective evaluation; not clinical guidance"
        },
        "candidate": {
            "label": cyclosporine.label,
            "concentration_nM": cyclosporine.concentration_nanomolar,
            "target_ki_nM": cyclosporine.target_ki_nanomolar,
            "off_target_ki_nM": cyclosporine.off_target_ki_nanomolar,
            "max_energy_shift_kj_mol": cyclosporine.max_energy_shift_kj_mol,
            "safety_buffer_kj_mol": cyclosporine.safety_buffer_kj_mol,
            "target_occupancy": blocker.target_occupancy_fraction,
            "off_target_occupancy": blocker.off_target_occupancy_fraction,
            "efficacy_margin_kj_mol": blocker.efficacy_margin_kj_mol
        },
        "mechanistic_context": {
            "mimicry_gap_kj_mol": mimicry.mimicry_gap_kj_mol,
            "activation_excess_kj_mol": mimicry.activation_excess_kj_mol,
            "baseline_drive": baseline_drive,
            "standard_drive": standard_drive,
            "combo_drive": combo_drive,
            "transduction_efficiency": transduction_efficiency
        },
        "efficacy_projection": {
            "horizon_years_trial": two_year.years,
            "arr_standard": standard_2y.annualized_relapse_rate,
            "arr_combo": combo_2y.annualized_relapse_rate,
            "arr_rr_reduction_combo_vs_standard": rr_reduction_combo_vs_standard_2y,
            "horizon_years_long": ten_year.years,
            "lesion_standard_long": standard_10y.final_lesion_index,
            "lesion_combo_long": combo_10y.final_lesion_index,
            "disability_standard_long": standard_10y.final_disability_index,
            "disability_combo_long": combo_10y.final_disability_index
        },
        "pk_bridge": {
            "p25_ng_mL": pk.p25_ng_ml,
            "p50_ng_mL": pk.p50_ng_ml,
            "p75_ng_mL": pk.p75_ng_ml,
            "p95_ng_mL": pk.p95_ng_ml,
            "recommended_exposure_window_ng_mL": [pk.p25_ng_ml, pk.p75_ng_ml]
        },
        "safety_gate": {
            "overall_pass": safety.overall_pass,
            "prob_in_target_zone": safety.prob_in_target_zone,
            "prob_above_renal_caution": safety.prob_above_renal_caution,
            "prob_above_renal_high": safety.prob_above_renal_high,
            "prob_above_neuro_caution": safety.prob_above_neuro_caution
        },
        "trial_sizing_poisson": {
            "alpha_two_sided": 0.05,
            "power": 0.80,
            "followup_years": two_year.years,
            "n_per_arm": n_per_arm_2y,
            "n_total": 2.0 * n_per_arm_2y
        },
        "stopping_rule_template": {
            "renal_caution_trigger": format!("observed trough > {:.1} ng/mL in repeated sampling", safety_gate.windows.renal_caution_ng_ml),
            "renal_high_trigger": format!("observed trough > {:.1} ng/mL", safety_gate.windows.renal_high_ng_ml),
            "neuro_caution_trigger": format!("observed trough > {:.1} ng/mL", safety_gate.windows.neuro_caution_ng_ml),
            "action": "hold/escalation freeze and safety review"
        }
    });
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("serialize"))
        .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "ms_cyclosporine_trial_hypothesis: arr_reduction_combo_vs_standard_2y={:.3} n_per_arm_2y={:.1} safety_pass={}",
        rr_reduction_combo_vs_standard_2y,
        n_per_arm_2y,
        safety.overall_pass
    );

    let _ = baseline_10y;
    let _ = baseline_2y;
}
