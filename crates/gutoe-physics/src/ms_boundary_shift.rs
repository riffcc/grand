/*!
 * Shared reduced-order utilities for MS boundary-shift simulations.
 *
 * These helpers keep consistency across localization, tolerance, and combined
 * sweep bins. This is simulation scaffolding, not clinical guidance.
 */

use crate::MolecularMimicryScore;
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
pub struct MsSimParams {
    pub years: f64,
    pub base_relapse_rate_per_year: f64,
    pub lesion_growth_coeff: f64,
    pub relapse_lesion_impact: f64,
    pub repair_rate: f64,
    pub seasonality_amp: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct MsCourseSummary {
    pub annualized_relapse_rate: f64,
    pub cumulative_relapses: f64,
    pub final_lesion_index: f64,
    pub final_disability_index: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct MsBoundaryShiftInput {
    pub transduction_efficiency: f64,
    pub achieved_shift_kj_mol: f64,
    pub activation_excess_kj_mol: f64,
    pub overlap_score: f64,
    pub off_target_occupancy: f64,
    /// Scales off-target penalty contribution to activation threshold.
    pub off_target_penalty_scale: f64,
    /// Additional tolerance term that shifts immune decision boundary.
    pub tolerance_shift_kj_mol: f64,
    /// Converts blocker drive to combo drive under background standard therapy.
    pub standard_factor: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct MsBoundaryShiftScore {
    pub effective_shift_kj_mol: f64,
    pub off_target_penalty_kj_mol: f64,
    pub activation_after_kj_mol: f64,
    pub activation_score_after: f64,
    pub blocker_drive: f64,
    pub combo_drive: f64,
}

pub fn default_ms_sim_params() -> MsSimParams {
    MsSimParams {
        years: 10.0,
        base_relapse_rate_per_year: 0.60,
        lesion_growth_coeff: 0.12,
        relapse_lesion_impact: 0.25,
        repair_rate: 0.03,
        seasonality_amp: 0.12,
    }
}

pub fn boundary_shift_score(input: MsBoundaryShiftInput) -> MsBoundaryShiftScore {
    let eff = input.transduction_efficiency.clamp(0.0, 1.0);
    let effective_shift = input.achieved_shift_kj_mol.max(0.0) * eff;
    let off_target_penalty =
        0.15_f64 * input.off_target_penalty_scale.max(0.0) * input.off_target_occupancy.max(0.0);

    let activation_after = (input.activation_excess_kj_mol.max(0.0)
        - effective_shift
        - input.tolerance_shift_kj_mol.max(0.0)
        + off_target_penalty)
        .max(0.0);
    let activation_score_after = (activation_after / (activation_after + 2.0)).clamp(0.0, 1.0);
    let blocker_drive = input.overlap_score.clamp(0.0, 1.0) * activation_score_after;
    let combo_drive = blocker_drive * input.standard_factor.clamp(0.0, 1.0);

    MsBoundaryShiftScore {
        effective_shift_kj_mol: effective_shift,
        off_target_penalty_kj_mol: off_target_penalty,
        activation_after_kj_mol: activation_after,
        activation_score_after,
        blocker_drive,
        combo_drive,
    }
}

pub fn standard_factor_from_drives(baseline_drive: f64, standard_drive: f64) -> f64 {
    if baseline_drive > 0.0 {
        (standard_drive / baseline_drive).clamp(0.0, 1.0)
    } else {
        1.0
    }
}

pub fn simulate_ms_course(base_drive_index: f64, p: MsSimParams) -> MsCourseSummary {
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

    MsCourseSummary {
        annualized_relapse_rate: cum_relapses / p.years.max(1.0e-9),
        cumulative_relapses: cum_relapses,
        final_lesion_index: lesion,
        final_disability_index: (1.0 - (-lesion / 9.5).exp()).clamp(0.0, 1.0),
    }
}

pub fn poisson_n_per_arm_80pct(control_arr: f64, treatment_arr: f64, years: f64) -> f64 {
    let lc = control_arr.max(1.0e-9);
    let lt = treatment_arr.max(1.0e-9);
    let t = years.max(0.25);
    let delta = (lc - lt).abs().max(1.0e-9);
    let z_alpha: f64 = 1.96;
    let z_power: f64 = 0.84;
    ((z_alpha + z_power).powi(2) * (lc + lt)) / (delta * delta * t)
}

pub fn ms_boundary_context(
    mimicry: MolecularMimicryScore,
    standard_drive: f64,
    achieved_shift_kj_mol: f64,
    off_target_occupancy: f64,
    transduction_efficiency: f64,
    tolerance_shift_kj_mol: f64,
    off_target_penalty_scale: f64,
) -> MsBoundaryShiftScore {
    let baseline = mimicry.misrecognition_risk_index;
    let sf = standard_factor_from_drives(baseline, standard_drive);
    boundary_shift_score(MsBoundaryShiftInput {
        transduction_efficiency,
        achieved_shift_kj_mol,
        activation_excess_kj_mol: mimicry.activation_excess_kj_mol,
        overlap_score: mimicry.overlap_score,
        off_target_occupancy,
        off_target_penalty_scale,
        tolerance_shift_kj_mol,
        standard_factor: sf,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tolerance_shift_lowers_combo_drive() {
        let no_tol = boundary_shift_score(MsBoundaryShiftInput {
            transduction_efficiency: 0.30,
            achieved_shift_kj_mol: 2.6,
            activation_excess_kj_mol: 0.91,
            overlap_score: 0.71,
            off_target_occupancy: 0.10,
            off_target_penalty_scale: 1.0,
            tolerance_shift_kj_mol: 0.0,
            standard_factor: 0.19,
        });
        let with_tol = boundary_shift_score(MsBoundaryShiftInput {
            tolerance_shift_kj_mol: 0.30,
            ..MsBoundaryShiftInput {
                transduction_efficiency: 0.30,
                achieved_shift_kj_mol: 2.6,
                activation_excess_kj_mol: 0.91,
                overlap_score: 0.71,
                off_target_occupancy: 0.10,
                off_target_penalty_scale: 1.0,
                tolerance_shift_kj_mol: 0.0,
                standard_factor: 0.19,
            }
        });
        assert!(with_tol.combo_drive <= no_tol.combo_drive + 1.0e-12);
    }

    #[test]
    fn course_summary_is_bounded() {
        let s = simulate_ms_course(0.1, default_ms_sim_params());
        assert!(s.final_lesion_index >= 0.0);
        assert!((0.0..=1.0).contains(&s.final_disability_index));
        assert!(s.annualized_relapse_rate >= 0.0);
    }
}
