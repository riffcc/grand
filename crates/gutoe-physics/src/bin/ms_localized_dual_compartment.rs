//! Dual-compartment localization lane for MS:
//! - site efficacy compartment (immune decision boundary at target tissue)
//! - systemic compartment (blood exposure safety probabilities)
//!
//! This is a reduced-order simulation, not clinical guidance.

use gutoe_physics::{
    default_cyclosporine_pk_bridge_input, default_cyclosporine_safety_gate_input,
    default_ms_mimicry_input, default_natalizumab_proxy, default_ocrelizumab_proxy,
    default_ms_sim_params, evaluate_cyclosporine_safety_gate, evaluate_molecular_mimicry,
    evaluate_targeted_blocker_candidate, evaluate_therapy_effect, ms_boundary_context,
    poisson_n_per_arm_80pct, simulate_cyclosporine_pk_bridge, simulate_ms_course,
    summarize_cyclosporine_pk_bridge, TargetedBlockerCandidateInput, MsSimParams,
};
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

    let mut sim10: MsSimParams = default_ms_sim_params();
    sim10.years = env_f64("GUTOE_MS_SIM_YEARS", sim10.years);
    sim10.base_relapse_rate_per_year = env_f64("GUTOE_MS_SIM_BASE_RELAPSE_PER_YEAR", sim10.base_relapse_rate_per_year);
    sim10.lesion_growth_coeff = env_f64("GUTOE_MS_SIM_LESION_GROWTH", sim10.lesion_growth_coeff);
    sim10.relapse_lesion_impact = env_f64("GUTOE_MS_SIM_RELAPSE_IMPACT", sim10.relapse_lesion_impact);
    sim10.repair_rate = env_f64("GUTOE_MS_SIM_REPAIR_RATE", sim10.repair_rate);
    sim10.seasonality_amp = env_f64("GUTOE_MS_SIM_SEASONALITY_AMP", sim10.seasonality_amp);
    let sim2 = MsSimParams {
        years: env_f64("GUTOE_MS_TRIAL_HORIZON_YEARS", 2.0),
        ..sim10
    };

    let mimicry = evaluate_molecular_mimicry(default_ms_mimicry_input());
    let baseline_drive = mimicry.misrecognition_risk_index;
    let standard = evaluate_therapy_effect(
        baseline_drive,
        default_ocrelizumab_proxy(),
        default_natalizumab_proxy(),
    );
    let standard_drive = standard.residual_drive_index;

    let blocker = evaluate_targeted_blocker_candidate(mimicry.activation_excess_kj_mol, candidate);
    let boundary = ms_boundary_context(
        mimicry,
        standard_drive,
        blocker.achieved_energy_shift_kj_mol * site_enrichment_factor,
        blocker.off_target_occupancy_fraction,
        transduction_efficiency,
        tolerance_shift_kj_mol,
        localization_factor,
    );

    let standard_2y = simulate_ms_course(standard_drive, sim2);
    let combo_2y = simulate_ms_course(boundary.combo_drive, sim2);
    let standard_10y = simulate_ms_course(standard_drive, sim10);
    let combo_10y = simulate_ms_course(boundary.combo_drive, sim10);

    let arr_reduction =
        (1.0 - combo_2y.annualized_relapse_rate / standard_2y.annualized_relapse_rate.max(1.0e-9))
            .clamp(-5.0, 1.0);
    let lesion_reduction_10y =
        (1.0 - combo_10y.final_lesion_index / standard_10y.final_lesion_index.max(1.0e-9))
            .clamp(-5.0, 1.0);

    let min_arr_reduction = env_f64("GUTOE_MS_MIN_ARR_REDUCTION", 0.10);
    let min_lesion_reduction = env_f64("GUTOE_MS_MIN_LESION_REDUCTION", 0.30);
    let efficacy_pass = arr_reduction >= min_arr_reduction && lesion_reduction_10y >= min_lesion_reduction;

    let mut pk_input = default_cyclosporine_pk_bridge_input();
    pk_input.site_target_nanomolar = candidate.concentration_nanomolar;
    pk_input.samples = env_usize("GUTOE_MS_PK_SAMPLES", pk_input.samples);
    pk_input.seed = env_u64("GUTOE_MS_PK_SEED", pk_input.seed);
    pk_input.blood_to_site_gain_median =
        env_f64("GUTOE_MS_PK_GAIN_MEDIAN", pk_input.blood_to_site_gain_median) * localization_factor;
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

    let overall_pass = efficacy_pass && safety.overall_pass;
    let n_per_arm = poisson_n_per_arm_80pct(
        standard_2y.annualized_relapse_rate,
        combo_2y.annualized_relapse_rate,
        sim2.years,
    );

    let out_dir = std::env::var("GUTOE_MS_LOCALIZED_DUAL_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ms_localized_dual_compartment".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("ms_localized_dual_compartment.txt");
    let json_path = out.join("ms_localized_dual_compartment.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ms_localized_dual_compartment]").expect("write");
    writeln!(txt, "candidate = {}", candidate.label).expect("write");
    writeln!(txt, "localization_factor = {:.6}", localization_factor).expect("write");
    writeln!(txt, "site_enrichment_factor = {:.6}", site_enrichment_factor).expect("write");
    writeln!(txt, "transduction_efficiency = {:.6}", transduction_efficiency).expect("write");
    writeln!(txt, "tolerance_shift_kj_mol = {:.6}", tolerance_shift_kj_mol).expect("write");
    writeln!(txt, "arr_reduction_combo_vs_standard_2y = {:.9}", arr_reduction).expect("write");
    writeln!(txt, "lesion_reduction_combo_vs_standard_10y = {:.9}", lesion_reduction_10y).expect("write");
    writeln!(txt, "disability_standard_10y = {:.9}", standard_10y.final_disability_index).expect("write");
    writeln!(txt, "disability_combo_10y = {:.9}", combo_10y.final_disability_index).expect("write");
    writeln!(txt, "pk_p50_ng_mL = {:.9}", pk.p50_ng_ml).expect("write");
    writeln!(txt, "pk_p95_ng_mL = {:.9}", pk.p95_ng_ml).expect("write");
    writeln!(txt, "safety_pass = {}", safety.overall_pass).expect("write");
    writeln!(txt, "efficacy_pass = {}", efficacy_pass).expect("write");
    writeln!(txt, "overall_pass = {}", overall_pass).expect("write");
    writeln!(txt, "poisson_n_per_arm_2y_80pct = {:.3}", n_per_arm).expect("write");

    let payload = json!({
        "meta": {
            "lane": "ms_localized_dual_compartment",
            "note": "dual-gate localization model; not clinical guidance"
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
        "compartment_controls": {
            "localization_factor": localization_factor,
            "site_enrichment_factor": site_enrichment_factor,
            "transduction_efficiency": transduction_efficiency,
            "tolerance_shift_kj_mol": tolerance_shift_kj_mol
        },
        "site_gate": {
            "effective_shift_kj_mol": boundary.effective_shift_kj_mol,
            "activation_after_kj_mol": boundary.activation_after_kj_mol,
            "blocker_drive": boundary.blocker_drive,
            "combo_drive": boundary.combo_drive,
            "arr_reduction_combo_vs_standard_2y": arr_reduction,
            "lesion_reduction_combo_vs_standard_10y": lesion_reduction_10y,
            "efficacy_pass": efficacy_pass
        },
        "systemic_gate": {
            "pk_p25_ng_mL": pk.p25_ng_ml,
            "pk_p50_ng_mL": pk.p50_ng_ml,
            "pk_p75_ng_mL": pk.p75_ng_ml,
            "pk_p95_ng_mL": pk.p95_ng_ml,
            "prob_above_renal_caution": safety.prob_above_renal_caution,
            "prob_above_renal_high": safety.prob_above_renal_high,
            "prob_above_neuro_caution": safety.prob_above_neuro_caution,
            "prob_in_target_zone": safety.prob_in_target_zone,
            "safety_pass": safety.overall_pass
        },
        "integration": {
            "overall_pass": overall_pass,
            "n_per_arm_2y_80pct": n_per_arm
        }
    });
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("serialize"))
        .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "ms_localized_dual_compartment: overall_pass={} efficacy_pass={} safety_pass={} arr_reduction={:.3}",
        overall_pass, efficacy_pass, safety.overall_pass, arr_reduction
    );
}
