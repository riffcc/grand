//! Tolerance-induction dynamics lane for MS decision-boundary shift.
//!
//! Models a time-dependent tolerance shift term (Δ_tol) that changes immune
//! decision boundary without requiring stronger broad suppression.
//!
//! This is reduced-order simulation scaffolding, not clinical guidance.

use gutoe_physics::{
    default_ms_mimicry_input, default_natalizumab_proxy, default_ocrelizumab_proxy,
    evaluate_molecular_mimicry, evaluate_targeted_blocker_candidate, evaluate_therapy_effect,
    ms_boundary_context, standard_factor_from_drives, TargetedBlockerCandidateInput,
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
struct ToleranceParams {
    induction_months: u32,
    induction_rate_kj_per_month: f64,
    maintenance_rate_kj_per_month: f64,
    half_life_months: f64,
}

#[derive(Clone, Copy, Debug)]
struct Point {
    month: u32,
    tolerance_shift_kj_mol: f64,
    combo_drive: f64,
    relapse_probability: f64,
    cumulative_relapses: f64,
    lesion_index: f64,
    disability_index: f64,
}

#[derive(Clone, Copy, Debug)]
struct Summary {
    annualized_relapse_rate: f64,
    cumulative_relapses: f64,
    final_tolerance_shift_kj_mol: f64,
    mean_tolerance_shift_kj_mol: f64,
    final_lesion_index: f64,
    final_disability_index: f64,
}

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn evolve_tolerance(t_prev: f64, month: u32, p: ToleranceParams) -> f64 {
    let decay = (-std::f64::consts::LN_2 / p.half_life_months.max(1.0e-6)).exp();
    let add = if month <= p.induction_months {
        p.induction_rate_kj_per_month
    } else {
        p.maintenance_rate_kj_per_month
    }
    .max(0.0);
    (t_prev.max(0.0) * decay + add).max(0.0)
}

fn simulate_with_tolerance(
    sim: SimParams,
    tol: ToleranceParams,
    baseline_mimicry: gutoe_physics::MolecularMimicryScore,
    standard_factor: f64,
    achieved_shift_kj_mol: f64,
    off_target_occupancy: f64,
    transduction_efficiency: f64,
    off_target_penalty_scale: f64,
    site_enrichment_factor: f64,
) -> (Summary, Vec<Point>) {
    let months = (sim.years.max(0.25) * 12.0).round() as u32;
    let mut lesion = 1.0_f64;
    let mut cumulative_relapses = 0.0_f64;
    let mut tol_shift = 0.0_f64;
    let mut tol_sum = 0.0_f64;
    let mut points = Vec::with_capacity(months as usize + 1);

    for m in 0..=months {
        if m > 0 {
            tol_shift = evolve_tolerance(tol_shift, m, tol);
        }

        let boundary = ms_boundary_context(
            baseline_mimicry,
            baseline_mimicry.misrecognition_risk_index * standard_factor,
            achieved_shift_kj_mol * site_enrichment_factor,
            off_target_occupancy,
            transduction_efficiency,
            tol_shift,
            off_target_penalty_scale,
        );

        let t = m as f64 / 12.0;
        let seasonal = 1.0 + sim.seasonality_amp * (2.0 * PI * t).sin();
        let micro = 1.0 + 0.08 * (2.0 * PI * t * 2.0 + 0.7).cos();
        let drive = (boundary.combo_drive.max(0.0) * seasonal * micro).clamp(0.0, 1.0);

        let monthly_rate = (sim.base_relapse_rate_per_year / 12.0) * (0.35 + 2.2 * drive);
        let relapse_prob = (1.0 - (-monthly_rate).exp()).clamp(0.0, 1.0);
        cumulative_relapses += relapse_prob;

        let growth = sim.lesion_growth_coeff * drive + sim.relapse_lesion_impact * relapse_prob;
        let repair = sim.repair_rate * lesion * (1.0 - 0.45 * drive);
        lesion = (lesion + growth - repair).max(0.0);

        let disability = (1.0 - (-lesion / 9.5).exp()).clamp(0.0, 1.0);
        tol_sum += tol_shift;
        points.push(Point {
            month: m,
            tolerance_shift_kj_mol: tol_shift,
            combo_drive: drive,
            relapse_probability: relapse_prob,
            cumulative_relapses,
            lesion_index: lesion,
            disability_index: disability,
        });
    }

    let summary = Summary {
        annualized_relapse_rate: cumulative_relapses / sim.years.max(1.0e-9),
        cumulative_relapses,
        final_tolerance_shift_kj_mol: tol_shift,
        mean_tolerance_shift_kj_mol: tol_sum / (months as f64 + 1.0),
        final_lesion_index: lesion,
        final_disability_index: (1.0 - (-lesion / 9.5).exp()).clamp(0.0, 1.0),
    };
    (summary, points)
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

    let sim = SimParams {
        years: env_f64("GUTOE_MS_SIM_YEARS", 10.0),
        base_relapse_rate_per_year: env_f64("GUTOE_MS_SIM_BASE_RELAPSE_PER_YEAR", 0.60),
        lesion_growth_coeff: env_f64("GUTOE_MS_SIM_LESION_GROWTH", 0.12),
        relapse_lesion_impact: env_f64("GUTOE_MS_SIM_RELAPSE_IMPACT", 0.25),
        repair_rate: env_f64("GUTOE_MS_SIM_REPAIR_RATE", 0.03),
        seasonality_amp: env_f64("GUTOE_MS_SIM_SEASONALITY_AMP", 0.12),
    };
    let tol = ToleranceParams {
        induction_months: env_f64("GUTOE_MS_TOL_INDUCTION_MONTHS", 6.0).round().max(0.0) as u32,
        induction_rate_kj_per_month: env_f64("GUTOE_MS_TOL_INDUCTION_RATE_KJ_PER_MONTH", 0.08),
        maintenance_rate_kj_per_month: env_f64("GUTOE_MS_TOL_MAINT_RATE_KJ_PER_MONTH", 0.02),
        half_life_months: env_f64("GUTOE_MS_TOL_HALF_LIFE_MONTHS", 12.0),
    };

    let transduction_efficiency = env_f64("GUTOE_MS_BLOCKER_SHIFT_EFFICIENCY", 0.30).clamp(0.0, 1.0);
    let site_enrichment_factor = env_f64("GUTOE_MS_SITE_ENRICHMENT_FACTOR", 1.0).max(0.0);
    let off_target_penalty_scale = env_f64("GUTOE_MS_OFFTARGET_PENALTY_SCALE", 1.0).max(0.0);

    let mimicry = evaluate_molecular_mimicry(default_ms_mimicry_input());
    let standard = evaluate_therapy_effect(
        mimicry.misrecognition_risk_index,
        default_ocrelizumab_proxy(),
        default_natalizumab_proxy(),
    );
    let standard_factor =
        standard_factor_from_drives(mimicry.misrecognition_risk_index, standard.residual_drive_index);

    let blocker = evaluate_targeted_blocker_candidate(mimicry.activation_excess_kj_mol, candidate);

    // Dynamic tolerance run.
    let (tol_summary, tol_points) = simulate_with_tolerance(
        sim,
        tol,
        mimicry,
        standard_factor,
        blocker.achieved_energy_shift_kj_mol,
        blocker.off_target_occupancy_fraction,
        transduction_efficiency,
        off_target_penalty_scale,
        site_enrichment_factor,
    );

    // Control run: same candidate and dynamics, but no tolerance induction.
    let zero_tol = ToleranceParams {
        induction_months: 0,
        induction_rate_kj_per_month: 0.0,
        maintenance_rate_kj_per_month: 0.0,
        half_life_months: tol.half_life_months,
    };
    let (ctrl_summary, ctrl_points) = simulate_with_tolerance(
        sim,
        zero_tol,
        mimicry,
        standard_factor,
        blocker.achieved_energy_shift_kj_mol,
        blocker.off_target_occupancy_fraction,
        transduction_efficiency,
        off_target_penalty_scale,
        site_enrichment_factor,
    );

    let arr_reduction_vs_ctrl =
        (1.0 - tol_summary.annualized_relapse_rate / ctrl_summary.annualized_relapse_rate.max(1.0e-9))
            .clamp(-5.0, 1.0);
    let lesion_reduction_vs_ctrl =
        (1.0 - tol_summary.final_lesion_index / ctrl_summary.final_lesion_index.max(1.0e-9))
            .clamp(-5.0, 1.0);

    let out_dir = std::env::var("GUTOE_MS_TOLERANCE_DYNAMICS_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ms_tolerance_induction_dynamics".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("ms_tolerance_induction_dynamics.txt");
    let json_path = out.join("ms_tolerance_induction_dynamics.json");
    let csv_path = out.join("ms_tolerance_induction_trajectory.csv");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ms_tolerance_induction_dynamics]").expect("write");
    writeln!(txt, "candidate = {}", candidate.label).expect("write");
    writeln!(txt, "final_tolerance_shift_kj_mol = {:.9}", tol_summary.final_tolerance_shift_kj_mol)
        .expect("write");
    writeln!(txt, "mean_tolerance_shift_kj_mol = {:.9}", tol_summary.mean_tolerance_shift_kj_mol)
        .expect("write");
    writeln!(txt, "arr_reduction_vs_no_tolerance = {:.9}", arr_reduction_vs_ctrl).expect("write");
    writeln!(txt, "lesion_reduction_vs_no_tolerance = {:.9}", lesion_reduction_vs_ctrl).expect("write");
    writeln!(txt, "disability_no_tolerance = {:.9}", ctrl_summary.final_disability_index).expect("write");
    writeln!(txt, "disability_with_tolerance = {:.9}", tol_summary.final_disability_index).expect("write");

    let mut csv = String::from(
        "month,tolerance_shift_kj_mol,drive_with_tolerance,drive_no_tolerance,cum_relapses_with_tolerance,cum_relapses_no_tolerance,lesion_with_tolerance,lesion_no_tolerance,disability_with_tolerance,disability_no_tolerance\n",
    );
    let len = tol_points.len().min(ctrl_points.len());
    for i in 0..len {
        let a = tol_points[i];
        let b = ctrl_points[i];
        csv.push_str(&format!(
            "{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            a.month,
            a.tolerance_shift_kj_mol,
            a.combo_drive,
            b.combo_drive,
            a.cumulative_relapses,
            b.cumulative_relapses,
            a.lesion_index,
            b.lesion_index,
            a.disability_index,
            b.disability_index,
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    let payload = json!({
        "meta": {
            "lane": "ms_tolerance_induction_dynamics",
            "note": "dynamic boundary-shift tolerance simulation; not clinical guidance"
        },
        "candidate": {
            "label": candidate.label,
            "target_ki_nM": candidate.target_ki_nanomolar,
            "off_target_ki_nM": candidate.off_target_ki_nanomolar,
            "target_occupancy": blocker.target_occupancy_fraction,
            "off_target_occupancy": blocker.off_target_occupancy_fraction
        },
        "tolerance_model": {
            "induction_months": tol.induction_months,
            "induction_rate_kj_per_month": tol.induction_rate_kj_per_month,
            "maintenance_rate_kj_per_month": tol.maintenance_rate_kj_per_month,
            "half_life_months": tol.half_life_months,
            "final_tolerance_shift_kj_mol": tol_summary.final_tolerance_shift_kj_mol,
            "mean_tolerance_shift_kj_mol": tol_summary.mean_tolerance_shift_kj_mol
        },
        "comparison": {
            "arr_no_tolerance": ctrl_summary.annualized_relapse_rate,
            "arr_with_tolerance": tol_summary.annualized_relapse_rate,
            "arr_reduction_vs_no_tolerance": arr_reduction_vs_ctrl,
            "lesion_no_tolerance": ctrl_summary.final_lesion_index,
            "lesion_with_tolerance": tol_summary.final_lesion_index,
            "lesion_reduction_vs_no_tolerance": lesion_reduction_vs_ctrl,
            "disability_no_tolerance": ctrl_summary.final_disability_index,
            "disability_with_tolerance": tol_summary.final_disability_index
        }
    });
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("serialize"))
        .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "ms_tolerance_induction_dynamics: arr_reduction_vs_no_tolerance={:.3} lesion_reduction_vs_no_tolerance={:.3}",
        arr_reduction_vs_ctrl, lesion_reduction_vs_ctrl
    );
}
