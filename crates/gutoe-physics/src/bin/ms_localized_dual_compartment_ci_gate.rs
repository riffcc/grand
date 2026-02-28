//! CI gate for the MS localized dual-compartment lane.
//!
//! This gate enforces that:
//! - effect-site efficacy stays above configured thresholds, and
//! - systemic exposure safety probabilities stay within configured bounds.
//!
//! It is a simulation integrity gate, not clinical guidance.

use gutoe_physics::{
    default_cyclosporine_pk_bridge_input, default_cyclosporine_safety_gate_input,
    default_ms_mimicry_input, default_natalizumab_proxy, default_ocrelizumab_proxy,
    default_ms_sim_params, evaluate_cyclosporine_safety_gate, evaluate_molecular_mimicry,
    evaluate_targeted_blocker_candidate, evaluate_therapy_effect, ms_boundary_context,
    poisson_n_per_arm_80pct, simulate_cyclosporine_pk_bridge, simulate_ms_course,
    TargetedBlockerCandidateInput,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
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
    let candidate = TargetedBlockerCandidateInput {
        label: "cyclosporine__c20nM__buf3",
        concentration_nanomolar: env_f64("GUTOE_MS_CANDIDATE_CONC_NM", 20.0),
        target_ki_nanomolar: env_f64("GUTOE_MS_CANDIDATE_TARGET_KI_NM", 2.64),
        off_target_ki_nanomolar: env_f64("GUTOE_MS_CANDIDATE_OFFTARGET_KI_NM", 200.0),
        max_energy_shift_kj_mol: env_f64("GUTOE_MS_CANDIDATE_MAX_SHIFT_KJ_MOL", 3.0),
        safety_buffer_kj_mol: env_f64("GUTOE_MS_CANDIDATE_SAFETY_BUFFER_KJ_MOL", 0.3),
    };

    let localization_factor = env_f64("GUTOE_MS_LOCALIZATION_FACTOR", 0.60).max(0.0);
    let site_enrichment_factor = env_f64("GUTOE_MS_SITE_ENRICHMENT_FACTOR", 1.0).max(0.0);
    let transduction_efficiency = env_f64("GUTOE_MS_BLOCKER_SHIFT_EFFICIENCY", 0.30).clamp(0.0, 1.0);
    let tolerance_shift_kj_mol = env_f64("GUTOE_MS_TOLERANCE_SHIFT_KJ_MOL", 0.0).max(0.0);

    let mut sim10 = default_ms_sim_params();
    sim10.years = env_f64("GUTOE_MS_SIM_YEARS", sim10.years);
    sim10.base_relapse_rate_per_year =
        env_f64("GUTOE_MS_SIM_BASE_RELAPSE_PER_YEAR", sim10.base_relapse_rate_per_year);
    sim10.lesion_growth_coeff = env_f64("GUTOE_MS_SIM_LESION_GROWTH", sim10.lesion_growth_coeff);
    sim10.relapse_lesion_impact =
        env_f64("GUTOE_MS_SIM_RELAPSE_IMPACT", sim10.relapse_lesion_impact);
    sim10.repair_rate = env_f64("GUTOE_MS_SIM_REPAIR_RATE", sim10.repair_rate);
    sim10.seasonality_amp = env_f64("GUTOE_MS_SIM_SEASONALITY_AMP", sim10.seasonality_amp);
    let sim2 = gutoe_physics::MsSimParams {
        years: env_f64("GUTOE_MS_TRIAL_HORIZON_YEARS", 2.0),
        ..sim10
    };

    let mimicry = evaluate_molecular_mimicry(default_ms_mimicry_input());
    let standard = evaluate_therapy_effect(
        mimicry.misrecognition_risk_index,
        default_ocrelizumab_proxy(),
        default_natalizumab_proxy(),
    );

    let blocker = evaluate_targeted_blocker_candidate(mimicry.activation_excess_kj_mol, candidate);
    let boundary = ms_boundary_context(
        mimicry,
        standard.residual_drive_index,
        blocker.achieved_energy_shift_kj_mol * site_enrichment_factor,
        blocker.off_target_occupancy_fraction,
        transduction_efficiency,
        tolerance_shift_kj_mol,
        localization_factor,
    );

    let standard_2y = simulate_ms_course(standard.residual_drive_index, sim2);
    let combo_2y = simulate_ms_course(boundary.combo_drive, sim2);
    let standard_10y = simulate_ms_course(standard.residual_drive_index, sim10);
    let combo_10y = simulate_ms_course(boundary.combo_drive, sim10);

    let arr_reduction_2y = (1.0
        - combo_2y.annualized_relapse_rate / standard_2y.annualized_relapse_rate.max(1.0e-9))
    .clamp(-5.0, 1.0);
    let lesion_reduction_10y =
        (1.0 - combo_10y.final_lesion_index / standard_10y.final_lesion_index.max(1.0e-9))
            .clamp(-5.0, 1.0);

    let min_arr_reduction = env_f64("GUTOE_MS_MIN_ARR_REDUCTION", 0.10);
    let min_lesion_reduction = env_f64("GUTOE_MS_MIN_LESION_REDUCTION", 0.30);
    let efficacy_pass = arr_reduction_2y >= min_arr_reduction && lesion_reduction_10y >= min_lesion_reduction;

    let mut pk_input = default_cyclosporine_pk_bridge_input();
    pk_input.site_target_nanomolar = candidate.concentration_nanomolar;
    pk_input.samples = env_usize("GUTOE_MS_PK_SAMPLES", pk_input.samples);
    pk_input.seed = env_u64("GUTOE_MS_PK_SEED", pk_input.seed);
    pk_input.blood_to_site_gain_median =
        env_f64("GUTOE_MS_PK_GAIN_MEDIAN", pk_input.blood_to_site_gain_median) * localization_factor;
    pk_input.blood_to_site_gain_gsd = env_f64("GUTOE_MS_PK_GAIN_GSD", pk_input.blood_to_site_gain_gsd);
    let pk_ensemble = simulate_cyclosporine_pk_bridge(pk_input);
    let pk_summary = gutoe_physics::summarize_cyclosporine_pk_bridge(&pk_ensemble);

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
    let safety = evaluate_cyclosporine_safety_gate(&pk_ensemble, safety_gate);

    let overall_pass = efficacy_pass && safety.overall_pass;
    let n_per_arm_2y_80pct = poisson_n_per_arm_80pct(
        standard_2y.annualized_relapse_rate,
        combo_2y.annualized_relapse_rate,
        sim2.years,
    );

    let out_dir = std::env::var("GUTOE_MS_LOCALIZED_DUAL_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/ms_localized_dual_compartment_ci_gate.json");

    let payload = json!({
        "overall_pass": overall_pass,
        "meta": {
            "lane": "ms_localized_dual_compartment_ci_gate",
            "note": "dual-gate simulation integrity check; not clinical guidance"
        },
        "candidate": {
            "label": candidate.label,
            "concentration_nM": candidate.concentration_nanomolar,
            "target_ki_nM": candidate.target_ki_nanomolar,
            "off_target_ki_nM": candidate.off_target_ki_nanomolar,
            "target_occupancy": blocker.target_occupancy_fraction,
            "off_target_occupancy": blocker.off_target_occupancy_fraction,
            "efficacy_margin_kj_mol": blocker.efficacy_margin_kj_mol
        },
        "controls": {
            "localization_factor": localization_factor,
            "site_enrichment_factor": site_enrichment_factor,
            "transduction_efficiency": transduction_efficiency,
            "tolerance_shift_kj_mol": tolerance_shift_kj_mol
        },
        "thresholds": {
            "min_arr_reduction_2y": min_arr_reduction,
            "min_lesion_reduction_10y": min_lesion_reduction,
            "max_prob_above_renal_caution": safety_gate.max_prob_above_renal_caution,
            "max_prob_above_renal_high": safety_gate.max_prob_above_renal_high,
            "max_prob_above_neuro_caution": safety_gate.max_prob_above_neuro_caution,
            "min_prob_in_target_zone": safety_gate.min_prob_in_target_zone
        },
        "score": {
            "effective_shift_kj_mol": boundary.effective_shift_kj_mol,
            "activation_after_kj_mol": boundary.activation_after_kj_mol,
            "combo_drive": boundary.combo_drive,
            "arr_reduction_2y": arr_reduction_2y,
            "lesion_reduction_10y": lesion_reduction_10y,
            "disability_standard_10y": standard_10y.final_disability_index,
            "disability_combo_10y": combo_10y.final_disability_index,
            "pk_p50_ng_mL": pk_summary.p50_ng_ml,
            "pk_p95_ng_mL": pk_summary.p95_ng_ml,
            "prob_above_renal_caution": safety.prob_above_renal_caution,
            "prob_above_renal_high": safety.prob_above_renal_high,
            "prob_above_neuro_caution": safety.prob_above_neuro_caution,
            "prob_in_target_zone": safety.prob_in_target_zone,
            "n_per_arm_2y_80pct": n_per_arm_2y_80pct
        },
        "gate": {
            "efficacy_pass": efficacy_pass,
            "safety_pass": safety.overall_pass,
            "overall_pass": overall_pass
        }
    });

    let mut json_file = File::create(&json_path).expect("create gate json");
    writeln!(
        json_file,
        "{}",
        serde_json::to_string_pretty(&payload).expect("serialize")
    )
    .expect("write gate json");

    println!(
        "ms_localized_dual_compartment_ci_gate: pass={} efficacy_pass={} safety_pass={} arr_red_2y={:.3} lesion_red_10y={:.3}",
        overall_pass, efficacy_pass, safety.overall_pass, arr_reduction_2y, lesion_reduction_10y
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: efficacy_pass={} safety_pass={} arr_red_2y={:.6} lesion_red_10y={:.6} p_renal_caution={:.6} p_renal_high={:.6} p_neuro={:.6} p_in_zone={:.6}",
            efficacy_pass,
            safety.overall_pass,
            arr_reduction_2y,
            lesion_reduction_10y,
            safety.prob_above_renal_caution,
            safety.prob_above_renal_high,
            safety.prob_above_neuro_caution,
            safety.prob_in_target_zone
        );
        process::exit(2);
    }
}
