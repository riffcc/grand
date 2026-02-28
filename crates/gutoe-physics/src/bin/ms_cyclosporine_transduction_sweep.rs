//! Fragility probe for cyclosporine MS lane:
//! 1) Efficiency sensitivity sweep (Ki -> effective shift -> ARR proxy)
//! 2) Off-target adverse-event noise amplification stress test
//!
//! This is a simulation sensitivity analysis, not clinical guidance.

use gutoe_physics::{
    default_ms_mimicry_input, default_natalizumab_proxy, default_ocrelizumab_proxy,
    evaluate_molecular_mimicry, evaluate_targeted_blocker_candidate, evaluate_therapy_effect,
    TargetedBlockerCandidateInput,
};
use rand::{rngs::StdRng, Rng, SeedableRng};
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
struct SweepRow {
    transduction_efficiency: f64,
    blocker_drive: f64,
    combo_drive: f64,
    arr_standard_2y: f64,
    arr_combo_2y: f64,
    arr_reduction_combo_vs_standard: f64,
    lesion_standard_10y: f64,
    lesion_combo_10y: f64,
    lesion_reduction_combo_vs_standard: f64,
    disability_standard_10y: f64,
    disability_combo_10y: f64,
    n_per_arm_2y_80pct: f64,
}

#[derive(Clone, Copy, Debug)]
struct AENoiseParams {
    samples: usize,
    seed: u64,
    // Base event rate at off-target occupancy = 1, amplification = 1.
    base_event_rate_per_year: f64,
    // Event impacts (sampled normal, clamped >=0)
    disability_event_mean: f64,
    disability_event_sd: f64,
    lesion_event_mean: f64,
    lesion_event_sd: f64,
    // Persistence of disability burden from adverse events.
    disability_half_life_months: f64,
}

#[derive(Clone, Copy, Debug)]
struct NoiseSweepRow {
    amplification: f64,
    eff_ref: f64,
    lesion_reduction_mean_vs_standard: f64,
    lesion_reduction_p05_vs_standard: f64,
    lesion_reduction_p95_vs_standard: f64,
    disability_mean: f64,
    disability_p95: f64,
    prob_disability_better_than_standard: f64,
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

fn standard_normal(rng: &mut StdRng) -> f64 {
    let u1 = rng.gen::<f64>().clamp(1.0e-12, 1.0 - 1.0e-12);
    let u2 = rng.gen::<f64>().clamp(1.0e-12, 1.0 - 1.0e-12);
    (-2.0_f64 * u1.ln()).sqrt() * (2.0_f64 * PI * u2).cos()
}

fn quantile(mut values: Vec<f64>, q: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let qq = q.clamp(0.0, 1.0);
    let idx = ((values.len() - 1) as f64 * qq).round() as usize;
    values[idx.min(values.len() - 1)]
}

fn parse_efficiencies(default: &[f64]) -> Vec<f64> {
    let raw = std::env::var("GUTOE_MS_EFF_LIST").ok();
    let mut vals = if let Some(s) = raw {
        s.split(',')
            .filter_map(|x| x.trim().parse::<f64>().ok())
            .map(|x| x.clamp(0.0, 1.0))
            .collect::<Vec<_>>()
    } else {
        default.to_vec()
    };
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    vals.dedup_by(|a, b| (*a - *b).abs() < 1.0e-12);
    vals
}

fn parse_amplifications(default: &[f64]) -> Vec<f64> {
    let raw = std::env::var("GUTOE_MS_AE_AMP_LIST").ok();
    let mut vals = if let Some(s) = raw {
        s.split(',')
            .filter_map(|x| x.trim().parse::<f64>().ok())
            .map(|x| x.max(0.0))
            .collect::<Vec<_>>()
    } else {
        default.to_vec()
    };
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    vals.dedup_by(|a, b| (*a - *b).abs() < 1.0e-12);
    vals
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

fn poisson_n_per_arm(lambda_control: f64, lambda_treatment: f64, years: f64) -> f64 {
    let lc = lambda_control.max(1.0e-9);
    let lt = lambda_treatment.max(1.0e-9);
    let t = years.max(0.25);
    let delta = (lc - lt).abs().max(1.0e-9);
    let z_alpha: f64 = 1.96;
    let z_power: f64 = 0.84;
    ((z_alpha + z_power).powi(2) * (lc + lt)) / (delta * delta * t)
}

fn simulate_combo_with_ae_noise(
    combo_drive: f64,
    p: SimParams,
    off_target_occupancy: f64,
    amplification: f64,
    ae: AENoiseParams,
    rng: &mut StdRng,
) -> (f64, f64) {
    let months = (p.years.max(0.25) * 12.0).round() as u32;
    let mut lesion = 1.0_f64;
    let mut disability_burden = 0.0_f64;
    let decay = (-std::f64::consts::LN_2 / ae.disability_half_life_months.max(1.0e-6)).exp();

    for m in 0..=months {
        let t = m as f64 / 12.0;
        let seasonal = 1.0 + p.seasonality_amp * (2.0 * PI * t).sin();
        let micro = 1.0 + 0.08 * (2.0 * PI * t * 2.0 + 0.7).cos();
        let drive = (combo_drive.max(0.0) * seasonal * micro).clamp(0.0, 1.0);

        let monthly_rate = (p.base_relapse_rate_per_year / 12.0) * (0.35 + 2.2 * drive);
        let relapse_prob = (1.0 - (-monthly_rate).exp()).clamp(0.0, 1.0);

        let growth = p.lesion_growth_coeff * drive + p.relapse_lesion_impact * relapse_prob;
        let repair = p.repair_rate * lesion * (1.0 - 0.45 * drive);
        lesion = (lesion + growth - repair).max(0.0);

        disability_burden *= decay;

        let event_rate = ae.base_event_rate_per_year.max(0.0)
            * off_target_occupancy.max(0.0)
            * amplification.max(0.0);
        let p_event = (1.0 - (-event_rate / 12.0).exp()).clamp(0.0, 1.0);

        if rng.gen::<f64>() < p_event {
            let z1 = standard_normal(rng);
            let z2 = standard_normal(rng);
            let dis_impact = (ae.disability_event_mean + ae.disability_event_sd * z1).max(0.0);
            let lesion_impact = (ae.lesion_event_mean + ae.lesion_event_sd * z2).max(0.0);
            disability_burden += dis_impact;
            lesion += lesion_impact;
        }
    }

    let base_disability = (1.0 - (-lesion / 9.5).exp()).clamp(0.0, 1.0);
    let noisy_disability = (base_disability + disability_burden).clamp(0.0, 1.0);
    (lesion, noisy_disability)
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

    let efficiencies = parse_efficiencies(&[0.15, 0.20, 0.25, 0.30, 0.35, 0.40]);
    let amplifications = parse_amplifications(&[1.0, 1.25, 1.50]);

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

    let blocker = evaluate_targeted_blocker_candidate(mimicry.activation_excess_kj_mol, candidate);

    let sim_2y = SimParams {
        years: env_f64("GUTOE_MS_TRIAL_HORIZON_YEARS", 2.0),
        base_relapse_rate_per_year: env_f64("GUTOE_MS_SIM_BASE_RELAPSE_PER_YEAR", 0.60),
        lesion_growth_coeff: env_f64("GUTOE_MS_SIM_LESION_GROWTH", 0.12),
        relapse_lesion_impact: env_f64("GUTOE_MS_SIM_RELAPSE_IMPACT", 0.25),
        repair_rate: env_f64("GUTOE_MS_SIM_REPAIR_RATE", 0.03),
        seasonality_amp: env_f64("GUTOE_MS_SIM_SEASONALITY_AMP", 0.12),
    };
    let sim_10y = SimParams {
        years: env_f64("GUTOE_MS_SIM_YEARS", 10.0),
        ..sim_2y
    };

    let standard_2y = simulate_course(standard_drive, sim_2y);
    let standard_10y = simulate_course(standard_drive, sim_10y);

    let mut rows = Vec::<SweepRow>::new();
    for eff in efficiencies.iter().copied() {
        let off_target_penalty = 0.15 * blocker.off_target_occupancy_fraction;
        let effective_shift = blocker.achieved_energy_shift_kj_mol * eff;
        let activation_after =
            (mimicry.activation_excess_kj_mol - effective_shift + off_target_penalty).max(0.0);
        let activation_score_after = (activation_after / (activation_after + 2.0)).clamp(0.0, 1.0);
        let blocker_drive = mimicry.overlap_score * activation_score_after;
        let combo_drive = blocker_drive * standard_factor;

        let combo_2y = simulate_course(combo_drive, sim_2y);
        let combo_10y = simulate_course(combo_drive, sim_10y);

        let arr_reduction_combo_vs_standard = (1.0
            - combo_2y.annualized_relapse_rate / standard_2y.annualized_relapse_rate.max(1.0e-9))
            .clamp(-5.0, 1.0);
        let lesion_reduction_combo_vs_standard = (1.0
            - combo_10y.final_lesion_index / standard_10y.final_lesion_index.max(1.0e-9))
            .clamp(-5.0, 1.0);

        rows.push(SweepRow {
            transduction_efficiency: eff,
            blocker_drive,
            combo_drive,
            arr_standard_2y: standard_2y.annualized_relapse_rate,
            arr_combo_2y: combo_2y.annualized_relapse_rate,
            arr_reduction_combo_vs_standard,
            lesion_standard_10y: standard_10y.final_lesion_index,
            lesion_combo_10y: combo_10y.final_lesion_index,
            lesion_reduction_combo_vs_standard,
            disability_standard_10y: standard_10y.final_disability_index,
            disability_combo_10y: combo_10y.final_disability_index,
            n_per_arm_2y_80pct: poisson_n_per_arm(
                standard_2y.annualized_relapse_rate,
                combo_2y.annualized_relapse_rate,
                sim_2y.years,
            ),
        });
    }

    let row_at_030 = rows
        .iter()
        .min_by(|a, b| {
            (a.transduction_efficiency - 0.30)
                .abs()
                .partial_cmp(&(b.transduction_efficiency - 0.30).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied()
        .expect("non-empty sweep");

    let ae = AENoiseParams {
        samples: env_usize("GUTOE_MS_AE_SAMPLES", 10_000),
        seed: env_u64("GUTOE_MS_AE_SEED", 4242),
        base_event_rate_per_year: env_f64("GUTOE_MS_AE_BASE_EVENT_RATE_PER_YEAR", 0.60),
        disability_event_mean: env_f64("GUTOE_MS_AE_DIS_IMPACT_MEAN", 0.006),
        disability_event_sd: env_f64("GUTOE_MS_AE_DIS_IMPACT_SD", 0.003),
        lesion_event_mean: env_f64("GUTOE_MS_AE_LESION_IMPACT_MEAN", 0.010),
        lesion_event_sd: env_f64("GUTOE_MS_AE_LESION_IMPACT_SD", 0.005),
        disability_half_life_months: env_f64("GUTOE_MS_AE_DIS_HALF_LIFE_MONTHS", 6.0),
    };

    let mut noise_rows = Vec::<NoiseSweepRow>::new();
    for amp in amplifications.iter().copied() {
        let mut lesion_reductions = Vec::with_capacity(ae.samples);
        let mut disabilities = Vec::with_capacity(ae.samples);
        let mut n_better_disability = 0usize;

        let mut rng = StdRng::seed_from_u64(ae.seed ^ ((amp * 10_000.0).round() as u64));
        for _ in 0..ae.samples {
            let (lesion_combo_noisy, disability_combo_noisy) = simulate_combo_with_ae_noise(
                row_at_030.combo_drive,
                sim_10y,
                blocker.off_target_occupancy_fraction,
                amp,
                ae,
                &mut rng,
            );

            let lesion_red = (1.0 - lesion_combo_noisy / standard_10y.final_lesion_index.max(1.0e-9))
                .clamp(-5.0, 1.0);
            lesion_reductions.push(lesion_red);
            disabilities.push(disability_combo_noisy);

            if disability_combo_noisy < standard_10y.final_disability_index {
                n_better_disability += 1;
            }
        }

        let mean_lesion_reduction = lesion_reductions.iter().sum::<f64>() / lesion_reductions.len() as f64;
        let p05_lesion_reduction = quantile(lesion_reductions.clone(), 0.05);
        let p95_lesion_reduction = quantile(lesion_reductions.clone(), 0.95);

        let mean_disability = disabilities.iter().sum::<f64>() / disabilities.len() as f64;
        let p95_disability = quantile(disabilities, 0.95);
        let p_better_disability = n_better_disability as f64 / ae.samples as f64;

        noise_rows.push(NoiseSweepRow {
            amplification: amp,
            eff_ref: row_at_030.transduction_efficiency,
            lesion_reduction_mean_vs_standard: mean_lesion_reduction,
            lesion_reduction_p05_vs_standard: p05_lesion_reduction,
            lesion_reduction_p95_vs_standard: p95_lesion_reduction,
            disability_mean: mean_disability,
            disability_p95: p95_disability,
            prob_disability_better_than_standard: p_better_disability,
        });
    }

    let out_dir = std::env::var("GUTOE_MS_CYCLOSPORINE_EFF_SWEEP_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ms_cyclosporine_transduction_sweep".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("ms_cyclosporine_transduction_sweep.txt");
    let csv_path = out.join("ms_cyclosporine_transduction_sweep.csv");
    let noise_csv_path = out.join("ms_cyclosporine_offtarget_noise_sweep.csv");
    let json_path = out.join("ms_cyclosporine_transduction_sweep.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ms_cyclosporine_transduction_sweep]").expect("write");
    writeln!(txt, "candidate = {}", candidate.label).expect("write");
    writeln!(txt, "efficiencies = {:?}", efficiencies).expect("write");
    writeln!(txt, "arr_reduction_vs_standard_at_eff:").expect("write");
    for r in &rows {
        writeln!(
            txt,
            "  eff={:.2} -> ARR reduction={:.6} lesion_reduction_10y={:.6}",
            r.transduction_efficiency, r.arr_reduction_combo_vs_standard, r.lesion_reduction_combo_vs_standard
        )
        .expect("write");
    }
    writeln!(txt, "offtarget_noise_ref_eff = {:.6}", row_at_030.transduction_efficiency).expect("write");
    for n in &noise_rows {
        writeln!(
            txt,
            "  amp={:.2} -> lesion_reduction_mean={:.6} (p05={:.6},p95={:.6}) disability_mean={:.6} disability_p95={:.6} p(disability<standard)={:.6}",
            n.amplification,
            n.lesion_reduction_mean_vs_standard,
            n.lesion_reduction_p05_vs_standard,
            n.lesion_reduction_p95_vs_standard,
            n.disability_mean,
            n.disability_p95,
            n.prob_disability_better_than_standard
        )
        .expect("write");
    }

    let mut csv = String::from(
        "transduction_efficiency,blocker_drive,combo_drive,arr_standard_2y,arr_combo_2y,arr_reduction_combo_vs_standard,lesion_standard_10y,lesion_combo_10y,lesion_reduction_combo_vs_standard,disability_standard_10y,disability_combo_10y,n_per_arm_2y_80pct\n",
    );
    for r in &rows {
        csv.push_str(&format!(
            "{:.6},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.3}\n",
            r.transduction_efficiency,
            r.blocker_drive,
            r.combo_drive,
            r.arr_standard_2y,
            r.arr_combo_2y,
            r.arr_reduction_combo_vs_standard,
            r.lesion_standard_10y,
            r.lesion_combo_10y,
            r.lesion_reduction_combo_vs_standard,
            r.disability_standard_10y,
            r.disability_combo_10y,
            r.n_per_arm_2y_80pct,
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    let mut noise_csv = String::from(
        "amplification,eff_ref,lesion_reduction_mean_vs_standard,lesion_reduction_p05_vs_standard,lesion_reduction_p95_vs_standard,disability_mean,disability_p95,prob_disability_better_than_standard\n",
    );
    for n in &noise_rows {
        noise_csv.push_str(&format!(
            "{:.6},{:.6},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            n.amplification,
            n.eff_ref,
            n.lesion_reduction_mean_vs_standard,
            n.lesion_reduction_p05_vs_standard,
            n.lesion_reduction_p95_vs_standard,
            n.disability_mean,
            n.disability_p95,
            n.prob_disability_better_than_standard,
        ));
    }
    fs::write(&noise_csv_path, noise_csv).expect("write noise csv");

    let payload = json!({
        "meta": {
            "lane": "ms_cyclosporine_transduction_sweep",
            "note": "efficiency sensitivity + off-target AE noise amplification"
        },
        "candidate": {
            "label": candidate.label,
            "concentration_nM": candidate.concentration_nanomolar,
            "target_ki_nM": candidate.target_ki_nanomolar,
            "off_target_ki_nM": candidate.off_target_ki_nanomolar,
            "achieved_shift_kj_mol": blocker.achieved_energy_shift_kj_mol,
            "required_shift_kj_mol": blocker.required_energy_shift_kj_mol,
            "efficacy_margin_kj_mol": blocker.efficacy_margin_kj_mol,
            "target_occupancy": blocker.target_occupancy_fraction,
            "off_target_occupancy": blocker.off_target_occupancy_fraction
        },
        "efficiency_sweep": rows.iter().map(|r| json!({
            "transduction_efficiency": r.transduction_efficiency,
            "arr_reduction_combo_vs_standard": r.arr_reduction_combo_vs_standard,
            "lesion_reduction_combo_vs_standard": r.lesion_reduction_combo_vs_standard,
            "n_per_arm_2y_80pct": r.n_per_arm_2y_80pct
        })).collect::<Vec<_>>(),
        "offtarget_noise_sweep": {
            "reference_efficiency": row_at_030.transduction_efficiency,
            "ae_params": {
                "samples": ae.samples,
                "seed": ae.seed,
                "base_event_rate_per_year": ae.base_event_rate_per_year,
                "disability_event_mean": ae.disability_event_mean,
                "disability_event_sd": ae.disability_event_sd,
                "lesion_event_mean": ae.lesion_event_mean,
                "lesion_event_sd": ae.lesion_event_sd,
                "disability_half_life_months": ae.disability_half_life_months
            },
            "rows": noise_rows.iter().map(|n| json!({
                "amplification": n.amplification,
                "lesion_reduction_mean_vs_standard": n.lesion_reduction_mean_vs_standard,
                "lesion_reduction_p05_vs_standard": n.lesion_reduction_p05_vs_standard,
                "lesion_reduction_p95_vs_standard": n.lesion_reduction_p95_vs_standard,
                "disability_mean": n.disability_mean,
                "disability_p95": n.disability_p95,
                "prob_disability_better_than_standard": n.prob_disability_better_than_standard
            })).collect::<Vec<_>>()
        }
    });
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("serialize"))
        .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", noise_csv_path.display());
    println!("wrote {}", json_path.display());
    if let Some(r30) = rows
        .iter()
        .min_by(|a, b| {
            (a.transduction_efficiency - 0.30)
                .abs()
                .partial_cmp(&(b.transduction_efficiency - 0.30).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    {
        println!(
            "ms_cyclosporine_transduction_sweep: ARR reduction @eff≈0.30 = {:.3}, lesion_reduction_10y @eff≈0.30 = {:.3}",
            r30.arr_reduction_combo_vs_standard,
            r30.lesion_reduction_combo_vs_standard
        );
    }
}
