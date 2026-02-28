//! Combined boundary-shift sweep for MS:
//! - transduction efficiency
//! - localization factor (systemic attenuation + off-target penalty attenuation)
//! - tolerance shift (decision-boundary term)
//!
//! Uses split gates:
//! - efficacy gate
//! - systemic safety gate
//!
//! This is a simulation ranking lane, not clinical guidance.

use gutoe_physics::{
    default_cyclosporine_pk_bridge_input, default_cyclosporine_safety_gate_input,
    default_ms_mimicry_input, default_natalizumab_proxy, default_ocrelizumab_proxy,
    default_ms_sim_params, evaluate_cyclosporine_safety_gate, evaluate_molecular_mimicry,
    evaluate_targeted_blocker_candidate, evaluate_therapy_effect, ms_boundary_context,
    poisson_n_per_arm_80pct, simulate_cyclosporine_pk_bridge, simulate_ms_course,
    TargetedBlockerCandidateInput,
};
use serde_json::json;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

#[derive(Clone, Copy, Debug)]
struct SafetyByLocalization {
    localization_factor: f64,
    p50_ng_ml: f64,
    p95_ng_ml: f64,
    p_renal_caution: f64,
    p_renal_high: f64,
    p_neuro_caution: f64,
    p_in_target_zone: f64,
    safety_pass: bool,
}

#[derive(Clone, Copy, Debug)]
struct SweepRow {
    transduction_efficiency: f64,
    localization_factor: f64,
    tolerance_shift_kj_mol: f64,
    combo_drive: f64,
    arr_reduction_2y: f64,
    lesion_reduction_10y: f64,
    disability_10y: f64,
    n_per_arm_2y_80pct: f64,
    p_renal_high: f64,
    p_in_target_zone: f64,
    efficacy_pass: bool,
    safety_pass: bool,
    overall_pass: bool,
    objective_score: f64,
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

fn parse_list(key: &str, default: &[f64]) -> Vec<f64> {
    let raw = std::env::var(key).ok();
    let mut vals = if let Some(s) = raw {
        s.split(',')
            .filter_map(|x| x.trim().parse::<f64>().ok())
            .collect::<Vec<_>>()
    } else {
        default.to_vec()
    };
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    vals.dedup_by(|a, b| (*a - *b).abs() < 1.0e-12);
    vals
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

    let efficiencies = parse_list("GUTOE_MS_EFF_LIST", &[0.15, 0.20, 0.25, 0.30, 0.35, 0.40]);
    let localizations = parse_list("GUTOE_MS_LOCALIZATION_LIST", &[1.0, 0.8, 0.6, 0.5]);
    let tolerance_shifts = parse_list("GUTOE_MS_TOL_SHIFT_LIST", &[0.0, 0.2, 0.4, 0.6]);
    let site_enrichment_factor = env_f64("GUTOE_MS_SITE_ENRICHMENT_FACTOR", 1.0).max(0.0);

    let mut sim10 = default_ms_sim_params();
    sim10.years = env_f64("GUTOE_MS_SIM_YEARS", sim10.years);
    sim10.base_relapse_rate_per_year = env_f64("GUTOE_MS_SIM_BASE_RELAPSE_PER_YEAR", sim10.base_relapse_rate_per_year);
    sim10.lesion_growth_coeff = env_f64("GUTOE_MS_SIM_LESION_GROWTH", sim10.lesion_growth_coeff);
    sim10.relapse_lesion_impact = env_f64("GUTOE_MS_SIM_RELAPSE_IMPACT", sim10.relapse_lesion_impact);
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
    let standard_drive = standard.residual_drive_index;

    let blocker = evaluate_targeted_blocker_candidate(mimicry.activation_excess_kj_mol, candidate);

    let standard_2y = simulate_ms_course(standard_drive, sim2);
    let standard_10y = simulate_ms_course(standard_drive, sim10);

    let min_arr_reduction = env_f64("GUTOE_MS_MIN_ARR_REDUCTION", 0.10);
    let min_lesion_reduction = env_f64("GUTOE_MS_MIN_LESION_REDUCTION", 0.30);

    // Precompute safety for each localization factor.
    let mut safety_map = HashMap::<u64, SafetyByLocalization>::new();
    for &lf in &localizations {
        let mut pk_input = default_cyclosporine_pk_bridge_input();
        pk_input.site_target_nanomolar = candidate.concentration_nanomolar;
        pk_input.samples = env_usize("GUTOE_MS_PK_SAMPLES", pk_input.samples);
        pk_input.seed = env_u64("GUTOE_MS_PK_SEED", pk_input.seed);
        pk_input.blood_to_site_gain_median =
            env_f64("GUTOE_MS_PK_GAIN_MEDIAN", pk_input.blood_to_site_gain_median) * lf.max(0.0);
        pk_input.blood_to_site_gain_gsd = env_f64("GUTOE_MS_PK_GAIN_GSD", pk_input.blood_to_site_gain_gsd);
        let pk_ens = simulate_cyclosporine_pk_bridge(pk_input);
        let pk = gutoe_physics::summarize_cyclosporine_pk_bridge(&pk_ens);

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

        safety_map.insert(
            (lf * 1_000_000.0).round() as u64,
            SafetyByLocalization {
                localization_factor: lf,
                p50_ng_ml: pk.p50_ng_ml,
                p95_ng_ml: pk.p95_ng_ml,
                p_renal_caution: safety.prob_above_renal_caution,
                p_renal_high: safety.prob_above_renal_high,
                p_neuro_caution: safety.prob_above_neuro_caution,
                p_in_target_zone: safety.prob_in_target_zone,
                safety_pass: safety.overall_pass,
            },
        );
    }

    let mut rows = Vec::<SweepRow>::new();
    for &eff in &efficiencies {
        for &lf in &localizations {
            for &tol in &tolerance_shifts {
                let boundary = ms_boundary_context(
                    mimicry,
                    standard_drive,
                    blocker.achieved_energy_shift_kj_mol * site_enrichment_factor,
                    blocker.off_target_occupancy_fraction,
                    eff,
                    tol,
                    lf.max(0.0),
                );
                let combo_2y = simulate_ms_course(boundary.combo_drive, sim2);
                let combo_10y = simulate_ms_course(boundary.combo_drive, sim10);

                let arr_reduction_2y = (1.0
                    - combo_2y.annualized_relapse_rate / standard_2y.annualized_relapse_rate.max(1.0e-9))
                    .clamp(-5.0, 1.0);
                let lesion_reduction_10y = (1.0
                    - combo_10y.final_lesion_index / standard_10y.final_lesion_index.max(1.0e-9))
                    .clamp(-5.0, 1.0);

                let sk = safety_map
                    .get(&((lf * 1_000_000.0).round() as u64))
                    .expect("safety map hit");
                let efficacy_pass =
                    arr_reduction_2y >= min_arr_reduction && lesion_reduction_10y >= min_lesion_reduction;
                let safety_pass = sk.safety_pass;
                let overall_pass = efficacy_pass && safety_pass;

                let objective = arr_reduction_2y
                    + 0.8 * lesion_reduction_10y
                    - 0.4 * sk.p_renal_high
                    - 0.2 * combo_10y.final_disability_index
                    + if overall_pass { 0.5 } else { 0.0 };

                rows.push(SweepRow {
                    transduction_efficiency: eff,
                    localization_factor: lf,
                    tolerance_shift_kj_mol: tol,
                    combo_drive: boundary.combo_drive,
                    arr_reduction_2y,
                    lesion_reduction_10y,
                    disability_10y: combo_10y.final_disability_index,
                    n_per_arm_2y_80pct: poisson_n_per_arm_80pct(
                        standard_2y.annualized_relapse_rate,
                        combo_2y.annualized_relapse_rate,
                        sim2.years,
                    ),
                    p_renal_high: sk.p_renal_high,
                    p_in_target_zone: sk.p_in_target_zone,
                    efficacy_pass,
                    safety_pass,
                    overall_pass,
                    objective_score: objective,
                });
            }
        }
    }

    rows.sort_by(|a, b| {
        b.objective_score
            .partial_cmp(&a.objective_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let out_dir = std::env::var("GUTOE_MS_BOUNDARY_SWEEP_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ms_boundary_shift_combined_sweep".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("ms_boundary_shift_combined_sweep.txt");
    let csv_path = out.join("ms_boundary_shift_combined_sweep.csv");
    let json_path = out.join("ms_boundary_shift_combined_sweep.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ms_boundary_shift_combined_sweep]").expect("write");
    writeln!(txt, "candidate = {}", candidate.label).expect("write");
    writeln!(txt, "rows = {}", rows.len()).expect("write");
    for (i, r) in rows.iter().take(10).enumerate() {
        writeln!(
            txt,
            "rank{} eff={:.2} loc={:.2} tol={:.2} arr_red={:.6} lesion_red={:.6} dis10y={:.6} p_renal_high={:.6} pass={}",
            i + 1,
            r.transduction_efficiency,
            r.localization_factor,
            r.tolerance_shift_kj_mol,
            r.arr_reduction_2y,
            r.lesion_reduction_10y,
            r.disability_10y,
            r.p_renal_high,
            r.overall_pass,
        )
        .expect("write");
    }

    let mut csv = String::from(
        "transduction_efficiency,localization_factor,tolerance_shift_kj_mol,combo_drive,arr_reduction_2y,lesion_reduction_10y,disability_10y,n_per_arm_2y_80pct,p_renal_high,p_in_target_zone,efficacy_pass,safety_pass,overall_pass,objective_score\n",
    );
    for r in &rows {
        csv.push_str(&format!(
            "{:.6},{:.6},{:.6},{:.9},{:.9},{:.9},{:.9},{:.3},{:.9},{:.9},{},{},{},{:.9}\n",
            r.transduction_efficiency,
            r.localization_factor,
            r.tolerance_shift_kj_mol,
            r.combo_drive,
            r.arr_reduction_2y,
            r.lesion_reduction_10y,
            r.disability_10y,
            r.n_per_arm_2y_80pct,
            r.p_renal_high,
            r.p_in_target_zone,
            r.efficacy_pass,
            r.safety_pass,
            r.overall_pass,
            r.objective_score,
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    let payload = json!({
        "meta": {
            "lane": "ms_boundary_shift_combined_sweep",
            "note": "combined localization+tolerance+efficiency sweep with split gates"
        },
        "candidate": {
            "label": candidate.label,
            "target_ki_nM": candidate.target_ki_nanomolar,
            "off_target_ki_nM": candidate.off_target_ki_nanomolar,
            "target_occupancy": blocker.target_occupancy_fraction,
            "off_target_occupancy": blocker.off_target_occupancy_fraction,
            "achieved_shift_kj_mol": blocker.achieved_energy_shift_kj_mol
        },
        "safety_by_localization": safety_map.values().map(|s| json!({
            "localization_factor": s.localization_factor,
            "p50_ng_mL": s.p50_ng_ml,
            "p95_ng_mL": s.p95_ng_ml,
            "p_renal_caution": s.p_renal_caution,
            "p_renal_high": s.p_renal_high,
            "p_neuro_caution": s.p_neuro_caution,
            "p_in_target_zone": s.p_in_target_zone,
            "safety_pass": s.safety_pass
        })).collect::<Vec<_>>(),
        "top_rows": rows.iter().take(20).map(|r| json!({
            "transduction_efficiency": r.transduction_efficiency,
            "localization_factor": r.localization_factor,
            "tolerance_shift_kj_mol": r.tolerance_shift_kj_mol,
            "arr_reduction_2y": r.arr_reduction_2y,
            "lesion_reduction_10y": r.lesion_reduction_10y,
            "disability_10y": r.disability_10y,
            "n_per_arm_2y_80pct": r.n_per_arm_2y_80pct,
            "p_renal_high": r.p_renal_high,
            "overall_pass": r.overall_pass,
            "objective_score": r.objective_score
        })).collect::<Vec<_>>()
    });
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("serialize"))
        .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    if let Some(best) = rows.first() {
        println!(
            "ms_boundary_shift_combined_sweep: best eff={:.2} loc={:.2} tol={:.2} arr_red={:.3} lesion_red={:.3} pass={}",
            best.transduction_efficiency,
            best.localization_factor,
            best.tolerance_shift_kj_mol,
            best.arr_reduction_2y,
            best.lesion_reduction_10y,
            best.overall_pass,
        );
    }
}
