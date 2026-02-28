/*!
 * Multiple-sclerosis molecular-mimicry lane (reduced-order).
 *
 * Scope:
 * - TCR interface energetics for a self epitope (MBP-like) and a pathogen-mimic epitope.
 * - Misrecognition risk from near-degenerate binding + tolerance-threshold excess.
 * - Therapeutic suppression proxy for ocrelizumab/natalizumab-like mechanisms.
 * - Targeted blocker feasibility proxy (interface-specific intervention concept).
 *
 * This is a mechanistic simulation scaffold, not clinical guidance.
 */

use crate::cardiovascular_binding::qed_contact_energy_kj_mol;

#[derive(Clone, Copy, Debug)]
pub struct InterfaceElectrostaticProxyInput {
    pub hbond_contact_count: f64,
    pub hbond_charge_product: f64,
    pub hbond_distance_nm: f64,
    pub hbond_dielectric: f64,
    pub polar_contact_count: f64,
    pub polar_charge_product: f64,
    pub polar_distance_nm: f64,
    pub polar_dielectric: f64,
    pub ionic_contact_count: f64,
    pub ionic_charge_product: f64,
    pub ionic_distance_nm: f64,
    pub ionic_dielectric: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct InterfaceResidualProxyInput {
    pub hydrophobic_area_a2: f64,
    pub hydrophobic_coeff_kj_per_a2: f64,
    pub aromatic_contact_count: f64,
    pub aromatic_contact_stabilization_kj: f64,
    pub water_release_count: f64,
    pub water_release_stabilization_kj: f64,
    pub conformational_entropy_penalty_kj: f64,
    pub polar_desolvation_penalty_kj: f64,
    pub strain_penalty_kj: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct InterfaceEnergyScore {
    pub qed_hbond_kj_mol: f64,
    pub qed_polar_kj_mol: f64,
    pub qed_ionic_kj_mol: f64,
    pub qed_total_kj_mol: f64,
    pub residual_total_kj_mol: f64,
    pub total_delta_g_kj_mol: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct MolecularMimicryInput {
    pub self_epitope_electro: InterfaceElectrostaticProxyInput,
    pub self_epitope_residual: InterfaceResidualProxyInput,
    pub mimic_epitope_electro: InterfaceElectrostaticProxyInput,
    pub mimic_epitope_residual: InterfaceResidualProxyInput,
    /// Self-tolerance energetic threshold; more negative self-binding than this can trigger risk.
    pub tolerance_threshold_kj_mol: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct MolecularMimicryScore {
    pub self_binding_kj_mol: f64,
    pub mimic_binding_kj_mol: f64,
    pub mimicry_gap_kj_mol: f64,
    pub activation_excess_kj_mol: f64,
    pub overlap_score: f64,
    pub activation_score: f64,
    pub misrecognition_risk_index: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct TherapyProxyInput {
    pub name: &'static str,
    pub concentration_nanomolar: f64,
    pub ki_nanomolar: f64,
    pub max_drive_reduction_fraction: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct TherapyEffectScore {
    pub name: &'static str,
    pub occupancy_fraction: f64,
    pub effective_drive_reduction_fraction: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct CombinedTherapyScore {
    pub baseline_drive_index: f64,
    pub ocrelizumab: TherapyEffectScore,
    pub natalizumab: TherapyEffectScore,
    pub residual_drive_index: f64,
    pub relative_drive_reduction_fraction: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct TargetedBlockerInput {
    pub concentration_nanomolar: f64,
    pub ki_nanomolar: f64,
    pub max_energy_shift_kj_mol: f64,
    pub safety_buffer_kj_mol: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct TargetedBlockerScore {
    pub occupancy_fraction: f64,
    pub achieved_energy_shift_kj_mol: f64,
    pub required_energy_shift_kj_mol: f64,
    pub required_occupancy_fraction: f64,
    pub feasible_at_given_concentration: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct TargetedBlockerCandidateInput {
    pub label: &'static str,
    pub concentration_nanomolar: f64,
    pub target_ki_nanomolar: f64,
    pub off_target_ki_nanomolar: f64,
    pub max_energy_shift_kj_mol: f64,
    pub safety_buffer_kj_mol: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct TargetedBlockerCandidateScore {
    pub label: &'static str,
    pub concentration_nanomolar: f64,
    pub target_ki_nanomolar: f64,
    pub off_target_ki_nanomolar: f64,
    pub max_energy_shift_kj_mol: f64,
    pub target_occupancy_fraction: f64,
    pub off_target_occupancy_fraction: f64,
    pub selectivity_ratio: f64,
    pub achieved_energy_shift_kj_mol: f64,
    pub required_energy_shift_kj_mol: f64,
    pub efficacy_margin_kj_mol: f64,
    pub feasible: bool,
    pub candidate_score: f64,
}

pub fn default_ms_mimicry_input() -> MolecularMimicryInput {
    MolecularMimicryInput {
        self_epitope_electro: InterfaceElectrostaticProxyInput {
            hbond_contact_count: 3.8,
            hbond_charge_product: 0.18,
            hbond_distance_nm: 0.30,
            hbond_dielectric: 27.0,
            polar_contact_count: 2.8,
            polar_charge_product: 0.12,
            polar_distance_nm: 0.34,
            polar_dielectric: 30.0,
            ionic_contact_count: 0.7,
            ionic_charge_product: 0.55,
            ionic_distance_nm: 0.33,
            ionic_dielectric: 26.0,
        },
        self_epitope_residual: InterfaceResidualProxyInput {
            hydrophobic_area_a2: 185.0,
            hydrophobic_coeff_kj_per_a2: 0.050,
            aromatic_contact_count: 1.8,
            aromatic_contact_stabilization_kj: 1.35,
            water_release_count: 1.8,
            water_release_stabilization_kj: 1.10,
            conformational_entropy_penalty_kj: 3.9,
            polar_desolvation_penalty_kj: 1.6,
            strain_penalty_kj: 0.8,
        },
        mimic_epitope_electro: InterfaceElectrostaticProxyInput {
            hbond_contact_count: 3.9,
            hbond_charge_product: 0.18,
            hbond_distance_nm: 0.30,
            hbond_dielectric: 27.0,
            polar_contact_count: 2.9,
            polar_charge_product: 0.12,
            polar_distance_nm: 0.34,
            polar_dielectric: 30.0,
            ionic_contact_count: 0.72,
            ionic_charge_product: 0.55,
            ionic_distance_nm: 0.33,
            ionic_dielectric: 26.0,
        },
        mimic_epitope_residual: InterfaceResidualProxyInput {
            hydrophobic_area_a2: 187.0,
            hydrophobic_coeff_kj_per_a2: 0.050,
            aromatic_contact_count: 1.85,
            aromatic_contact_stabilization_kj: 1.35,
            water_release_count: 1.85,
            water_release_stabilization_kj: 1.10,
            conformational_entropy_penalty_kj: 3.95,
            polar_desolvation_penalty_kj: 1.6,
            strain_penalty_kj: 0.8,
        },
        tolerance_threshold_kj_mol: -29.0,
    }
}

pub fn default_ocrelizumab_proxy() -> TherapyProxyInput {
    TherapyProxyInput {
        name: "ocrelizumab_like_bcell_depletion",
        concentration_nanomolar: 10.0,
        ki_nanomolar: 1.5,
        max_drive_reduction_fraction: 0.65,
    }
}

pub fn default_natalizumab_proxy() -> TherapyProxyInput {
    TherapyProxyInput {
        name: "natalizumab_like_trafficking_block",
        concentration_nanomolar: 8.0,
        ki_nanomolar: 2.0,
        max_drive_reduction_fraction: 0.70,
    }
}

pub fn default_targeted_blocker_input() -> TargetedBlockerInput {
    TargetedBlockerInput {
        concentration_nanomolar: 30.0,
        ki_nanomolar: 15.0,
        max_energy_shift_kj_mol: 2.5,
        safety_buffer_kj_mol: 0.5,
    }
}

fn occupancy_fraction(concentration_nanomolar: f64, ki_nanomolar: f64) -> f64 {
    let c = concentration_nanomolar.max(0.0);
    let ki = ki_nanomolar.max(1.0e-12);
    (c / (c + ki)).clamp(0.0, 1.0)
}

pub fn evaluate_interface_energy(
    electro: InterfaceElectrostaticProxyInput,
    residual: InterfaceResidualProxyInput,
) -> InterfaceEnergyScore {
    let qed_hbond = electro.hbond_contact_count.max(0.0)
        * qed_contact_energy_kj_mol(
            electro.hbond_charge_product,
            electro.hbond_distance_nm,
            electro.hbond_dielectric,
        );
    let qed_polar = electro.polar_contact_count.max(0.0)
        * qed_contact_energy_kj_mol(
            electro.polar_charge_product,
            electro.polar_distance_nm,
            electro.polar_dielectric,
        );
    let qed_ionic = electro.ionic_contact_count.max(0.0)
        * qed_contact_energy_kj_mol(
            electro.ionic_charge_product,
            electro.ionic_distance_nm,
            electro.ionic_dielectric,
        );
    let qed_total = qed_hbond + qed_polar + qed_ionic;

    let residual_total = -(residual.hydrophobic_area_a2.max(0.0)
        * residual.hydrophobic_coeff_kj_per_a2.max(0.0))
        - (residual.aromatic_contact_count.max(0.0)
            * residual.aromatic_contact_stabilization_kj.max(0.0))
        - (residual.water_release_count.max(0.0)
            * residual.water_release_stabilization_kj.max(0.0))
        + residual.conformational_entropy_penalty_kj.max(0.0)
        + residual.polar_desolvation_penalty_kj.max(0.0)
        + residual.strain_penalty_kj.max(0.0);

    let total = qed_total + residual_total;
    InterfaceEnergyScore {
        qed_hbond_kj_mol: qed_hbond,
        qed_polar_kj_mol: qed_polar,
        qed_ionic_kj_mol: qed_ionic,
        qed_total_kj_mol: qed_total,
        residual_total_kj_mol: residual_total,
        total_delta_g_kj_mol: total,
    }
}

pub fn evaluate_molecular_mimicry(input: MolecularMimicryInput) -> MolecularMimicryScore {
    let self_energy = evaluate_interface_energy(input.self_epitope_electro, input.self_epitope_residual);
    let mimic_energy = evaluate_interface_energy(input.mimic_epitope_electro, input.mimic_epitope_residual);

    let self_binding = self_energy.total_delta_g_kj_mol;
    let mimic_binding = mimic_energy.total_delta_g_kj_mol;
    let gap = (mimic_binding - self_binding).abs();
    let activation_excess = (input.tolerance_threshold_kj_mol - self_binding).max(0.0);

    let overlap_score = (1.0 / (1.0 + gap / 2.0)).clamp(0.0, 1.0);
    let activation_score = (activation_excess / (activation_excess + 2.0)).clamp(0.0, 1.0);
    let risk = (overlap_score * activation_score).clamp(0.0, 1.0);

    MolecularMimicryScore {
        self_binding_kj_mol: self_binding,
        mimic_binding_kj_mol: mimic_binding,
        mimicry_gap_kj_mol: gap,
        activation_excess_kj_mol: activation_excess,
        overlap_score,
        activation_score,
        misrecognition_risk_index: risk,
    }
}

pub fn evaluate_therapy_effect(
    baseline_drive_index: f64,
    ocrelizumab: TherapyProxyInput,
    natalizumab: TherapyProxyInput,
) -> CombinedTherapyScore {
    let occ_o = occupancy_fraction(ocrelizumab.concentration_nanomolar, ocrelizumab.ki_nanomolar);
    let eff_o = (ocrelizumab.max_drive_reduction_fraction.max(0.0).min(1.0) * occ_o).clamp(0.0, 1.0);

    let occ_n = occupancy_fraction(natalizumab.concentration_nanomolar, natalizumab.ki_nanomolar);
    let eff_n = (natalizumab.max_drive_reduction_fraction.max(0.0).min(1.0) * occ_n).clamp(0.0, 1.0);

    let residual = baseline_drive_index.max(0.0) * (1.0 - eff_o) * (1.0 - eff_n);
    let relative_reduction = if baseline_drive_index > 0.0 {
        1.0 - residual / baseline_drive_index
    } else {
        0.0
    };

    CombinedTherapyScore {
        baseline_drive_index,
        ocrelizumab: TherapyEffectScore {
            name: ocrelizumab.name,
            occupancy_fraction: occ_o,
            effective_drive_reduction_fraction: eff_o,
        },
        natalizumab: TherapyEffectScore {
            name: natalizumab.name,
            occupancy_fraction: occ_n,
            effective_drive_reduction_fraction: eff_n,
        },
        residual_drive_index: residual,
        relative_drive_reduction_fraction: relative_reduction.clamp(0.0, 1.0),
    }
}

pub fn evaluate_targeted_blocker(
    activation_excess_kj_mol: f64,
    blocker: TargetedBlockerInput,
) -> TargetedBlockerScore {
    let occ = occupancy_fraction(blocker.concentration_nanomolar, blocker.ki_nanomolar);
    let achieved_shift = blocker.max_energy_shift_kj_mol.max(0.0) * occ;
    let required_shift = (activation_excess_kj_mol.max(0.0) + blocker.safety_buffer_kj_mol.max(0.0)).max(0.0);
    let required_occ = if blocker.max_energy_shift_kj_mol > 0.0 {
        (required_shift / blocker.max_energy_shift_kj_mol).clamp(0.0, 1.5)
    } else {
        1.5
    };
    let feasible = achieved_shift + 1.0e-12 >= required_shift;

    TargetedBlockerScore {
        occupancy_fraction: occ,
        achieved_energy_shift_kj_mol: achieved_shift,
        required_energy_shift_kj_mol: required_shift,
        required_occupancy_fraction: required_occ,
        feasible_at_given_concentration: feasible,
    }
}

pub fn evaluate_targeted_blocker_candidate(
    activation_excess_kj_mol: f64,
    candidate: TargetedBlockerCandidateInput,
) -> TargetedBlockerCandidateScore {
    let target_occ = occupancy_fraction(candidate.concentration_nanomolar, candidate.target_ki_nanomolar);
    let off_occ = occupancy_fraction(candidate.concentration_nanomolar, candidate.off_target_ki_nanomolar);
    let selectivity_ratio = target_occ / off_occ.max(1.0e-9);

    let required = (activation_excess_kj_mol.max(0.0) + candidate.safety_buffer_kj_mol.max(0.0)).max(0.0);
    let achieved = candidate.max_energy_shift_kj_mol.max(0.0) * target_occ;
    let margin = achieved - required;
    let feasible = margin >= 0.0;

    // Heuristic objective:
    // - reward positive energy margin
    // - reward selectivity
    // - penalize off-target occupancy and concentration burden
    let score = (if feasible { 1.0 } else { 0.0 })
        + 0.8 * margin.max(0.0)
        + 0.25 * selectivity_ratio.min(20.0).ln_1p()
        - 1.2 * off_occ
        - 0.002 * candidate.concentration_nanomolar;

    TargetedBlockerCandidateScore {
        label: candidate.label,
        concentration_nanomolar: candidate.concentration_nanomolar,
        target_ki_nanomolar: candidate.target_ki_nanomolar,
        off_target_ki_nanomolar: candidate.off_target_ki_nanomolar,
        max_energy_shift_kj_mol: candidate.max_energy_shift_kj_mol,
        target_occupancy_fraction: target_occ,
        off_target_occupancy_fraction: off_occ,
        selectivity_ratio,
        achieved_energy_shift_kj_mol: achieved,
        required_energy_shift_kj_mol: required,
        efficacy_margin_kj_mol: margin,
        feasible,
        candidate_score: score,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mimicry_default_is_near_tipping_not_far_from_threshold() {
        let s = evaluate_molecular_mimicry(default_ms_mimicry_input());
        assert!(s.activation_excess_kj_mol > 0.0);
        assert!(s.activation_excess_kj_mol < 3.0);
        assert!(s.mimicry_gap_kj_mol < 2.0);
    }

    #[test]
    fn therapies_reduce_drive() {
        let mim = evaluate_molecular_mimicry(default_ms_mimicry_input());
        let baseline = mim.misrecognition_risk_index;
        let t = evaluate_therapy_effect(
            baseline,
            default_ocrelizumab_proxy(),
            default_natalizumab_proxy(),
        );
        assert!(t.residual_drive_index < baseline);
        assert!(t.relative_drive_reduction_fraction > 0.5);
    }

    #[test]
    fn targeted_blocker_can_be_feasible_in_default_lane() {
        let mim = evaluate_molecular_mimicry(default_ms_mimicry_input());
        let b = evaluate_targeted_blocker(mim.activation_excess_kj_mol, default_targeted_blocker_input());
        assert!(b.required_occupancy_fraction <= 1.0);
        assert!(b.feasible_at_given_concentration);
    }

    #[test]
    fn candidate_with_better_selectivity_scores_higher() {
        let mim = evaluate_molecular_mimicry(default_ms_mimicry_input());
        let base = TargetedBlockerCandidateInput {
            label: "base",
            concentration_nanomolar: 20.0,
            target_ki_nanomolar: 8.0,
            off_target_ki_nanomolar: 40.0,
            max_energy_shift_kj_mol: 2.5,
            safety_buffer_kj_mol: 0.5,
        };
        let alt = TargetedBlockerCandidateInput {
            off_target_ki_nanomolar: 200.0,
            ..base
        };
        let s1 = evaluate_targeted_blocker_candidate(mim.activation_excess_kj_mol, base);
        let s2 = evaluate_targeted_blocker_candidate(mim.activation_excess_kj_mol, alt);
        assert!(s2.selectivity_ratio > s1.selectivity_ratio);
        assert!(s2.candidate_score > s1.candidate_score);
    }
}
