/*!
 * Beta-lactamase inhibitor rescue lane (reduced-order, first-principles scaffold).
 *
 * Scope:
 * - Converts potency anchors (nM) to thermodynamic ΔG.
 * - Computes a QED electrostatic floor from α ħ c / (ε r).
 * - Adds structurally-motivated residual terms tied to inhibitor/enzyme mechanisms.
 * - Produces pairwise rankings for TEM-1, KPC, and NDM-1.
 *
 * Honesty:
 * - This is a simulation lane for ranking hypotheses, not clinical guidance.
 * - Potency anchors are assay-level ChEMBL snapshot priors and include uncertainty.
 */

use crate::chemical_thermo::{AVOGADRO, R_GAS_J_MOL_K};
use crate::{ALPHA_LEADING_ORDER, C, HBAR};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InhibitorScaffold {
    BetaLactamSuicide,
    Diazabicyclooctane,
    CyclicBoronate,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BetaLactamaseClass {
    SerineClassA,
    MetalloClassB,
}

#[derive(Clone, Copy, Debug)]
pub struct InhibitorSpec {
    pub name: &'static str,
    pub chembl_id: &'static str,
    pub scaffold: InhibitorScaffold,
    pub anionic_sites: f64,
    pub hbond_sites: f64,
    pub hydrophobic_surface_a2: f64,
    pub serine_trap_strength: f64,
    pub boronate_reversible_strength: f64,
    pub zinc_chelation_strength: f64,
    pub flexibility_penalty: f64,
    pub polar_desolvation_penalty: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct EnzymeSpec {
    pub name: &'static str,
    pub chembl_hint: &'static str,
    pub class: BetaLactamaseClass,
    pub serine_drive: f64,
    pub zinc_drive: f64,
    pub boronate_drive: f64,
    pub hbond_density: f64,
    pub steric_openness: f64,
    pub ionic_distance_nm: f64,
    pub hbond_distance_nm: f64,
    pub active_site_dielectric: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct PotencyAnchor {
    pub inhibitor_name: &'static str,
    pub enzyme_name: &'static str,
    pub anchor_nanomolar: f64,
    pub evidence_count: usize,
    pub imputed: bool,
    pub notes: &'static str,
}

#[derive(Clone, Copy, Debug)]
pub struct ResistanceModelCoefficients {
    pub hbond_charge_product: f64,
    pub ionic_contact_scale: f64,
    pub hbond_contact_scale: f64,
    pub serine_trap_kj: f64,
    pub boronate_serine_kj: f64,
    pub zinc_match_kj: f64,
    pub zinc_mismatch_penalty_kj: f64,
    pub hydrophobic_coeff_kj_per_a2: f64,
    pub entropy_coeff_kj: f64,
    pub desolv_coeff_kj: f64,
    pub offset_serine_class_kj: f64,
    pub offset_metallo_class_kj: f64,
}

impl Default for ResistanceModelCoefficients {
    fn default() -> Self {
        Self {
            hbond_charge_product: 0.20,
            ionic_contact_scale: 0.80,
            hbond_contact_scale: 0.70,
            serine_trap_kj: 8.0,
            boronate_serine_kj: 9.0,
            zinc_match_kj: 6.5,
            zinc_mismatch_penalty_kj: 8.0,
            hydrophobic_coeff_kj_per_a2: 0.010,
            entropy_coeff_kj: 1.0,
            desolv_coeff_kj: 1.0,
            offset_serine_class_kj: -13.0,
            offset_metallo_class_kj: -19.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct PairwiseResistanceResult {
    pub inhibitor_name: &'static str,
    pub inhibitor_chembl_id: &'static str,
    pub scaffold: InhibitorScaffold,
    pub enzyme_name: &'static str,
    pub enzyme_chembl_hint: &'static str,
    pub enzyme_class: BetaLactamaseClass,
    pub evidence_count: usize,
    pub imputed_anchor: bool,
    pub anchor_nanomolar: f64,
    pub anchor_delta_g_kj_mol: f64,
    pub qed_ionic_floor_kj_mol: f64,
    pub qed_hbond_floor_kj_mol: f64,
    pub qed_floor_total_kj_mol: f64,
    pub residual_modeled_total_kj_mol: f64,
    pub predicted_delta_g_kj_mol: f64,
    pub predicted_nanomolar: f64,
    pub log10_error_pred_vs_anchor: f64,
    pub occupancy_anchor_at_1u_m: f64,
    pub occupancy_predicted_at_1u_m: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct EnzymeBest {
    pub enzyme_name: &'static str,
    pub by_anchor_inhibitor: &'static str,
    pub by_anchor_nanomolar: f64,
    pub by_predicted_inhibitor: &'static str,
    pub by_predicted_nanomolar: f64,
    pub predicted_match_anchor_winner: bool,
}

#[derive(Clone, Debug)]
pub struct AntibioticResistancePanel {
    pub rows: Vec<PairwiseResistanceResult>,
    pub best_by_enzyme: Vec<EnzymeBest>,
    pub mean_abs_log10_error: f64,
    pub ndm_max_predicted_occupancy_at_1u_m: f64,
}

pub fn default_beta_lactamase_inhibitors() -> Vec<InhibitorSpec> {
    vec![
        InhibitorSpec {
            name: "clavulanic_acid",
            chembl_id: "CHEMBL777",
            scaffold: InhibitorScaffold::BetaLactamSuicide,
            anionic_sites: 1.0,
            hbond_sites: 2.0,
            hydrophobic_surface_a2: 260.0,
            serine_trap_strength: 0.65,
            boronate_reversible_strength: 0.0,
            zinc_chelation_strength: 0.05,
            flexibility_penalty: 2.8,
            polar_desolvation_penalty: 2.0,
        },
        InhibitorSpec {
            name: "sulbactam",
            chembl_id: "CHEMBL403",
            scaffold: InhibitorScaffold::BetaLactamSuicide,
            anionic_sites: 1.0,
            hbond_sites: 1.8,
            hydrophobic_surface_a2: 230.0,
            serine_trap_strength: 0.20,
            boronate_reversible_strength: 0.0,
            zinc_chelation_strength: 0.05,
            flexibility_penalty: 3.0,
            polar_desolvation_penalty: 2.1,
        },
        InhibitorSpec {
            name: "tazobactam",
            chembl_id: "CHEMBL404",
            scaffold: InhibitorScaffold::BetaLactamSuicide,
            anionic_sites: 1.0,
            hbond_sites: 2.1,
            hydrophobic_surface_a2: 250.0,
            serine_trap_strength: 0.55,
            boronate_reversible_strength: 0.0,
            zinc_chelation_strength: 0.08,
            flexibility_penalty: 2.6,
            polar_desolvation_penalty: 1.9,
        },
        InhibitorSpec {
            name: "avibactam",
            chembl_id: "CHEMBL1689063",
            scaffold: InhibitorScaffold::Diazabicyclooctane,
            anionic_sites: 1.0,
            hbond_sites: 3.1,
            hydrophobic_surface_a2: 240.0,
            serine_trap_strength: 0.90,
            boronate_reversible_strength: 0.0,
            zinc_chelation_strength: 0.10,
            flexibility_penalty: 2.4,
            polar_desolvation_penalty: 1.6,
        },
        InhibitorSpec {
            name: "vaborbactam",
            chembl_id: "CHEMBL3317857",
            scaffold: InhibitorScaffold::CyclicBoronate,
            anionic_sites: 1.0,
            hbond_sites: 2.4,
            hydrophobic_surface_a2: 270.0,
            serine_trap_strength: 0.08,
            boronate_reversible_strength: 0.95,
            zinc_chelation_strength: 0.05,
            flexibility_penalty: 2.1,
            polar_desolvation_penalty: 1.4,
        },
    ]
}

pub fn default_beta_lactamase_enzymes() -> Vec<EnzymeSpec> {
    vec![
        EnzymeSpec {
            name: "TEM-1",
            chembl_hint: "CHEMBL1287599",
            class: BetaLactamaseClass::SerineClassA,
            serine_drive: 1.00,
            zinc_drive: 0.00,
            boronate_drive: 0.20,
            hbond_density: 1.00,
            steric_openness: 1.00,
            ionic_distance_nm: 0.29,
            hbond_distance_nm: 0.30,
            active_site_dielectric: 24.0,
        },
        EnzymeSpec {
            name: "KPC",
            chembl_hint: "CHEMBL6132",
            class: BetaLactamaseClass::SerineClassA,
            serine_drive: 0.95,
            zinc_drive: 0.00,
            boronate_drive: 1.00,
            hbond_density: 0.95,
            steric_openness: 0.88,
            ionic_distance_nm: 0.30,
            hbond_distance_nm: 0.31,
            active_site_dielectric: 24.0,
        },
        EnzymeSpec {
            name: "NDM-1",
            chembl_hint: "CHEMBL4295540",
            class: BetaLactamaseClass::MetalloClassB,
            serine_drive: 0.00,
            zinc_drive: 1.00,
            boronate_drive: 0.05,
            hbond_density: 0.75,
            steric_openness: 0.90,
            ionic_distance_nm: 0.31,
            hbond_distance_nm: 0.32,
            active_site_dielectric: 26.0,
        },
    ]
}

pub fn default_potency_anchors() -> Vec<PotencyAnchor> {
    vec![
        PotencyAnchor {
            inhibitor_name: "clavulanic_acid",
            enzyme_name: "TEM-1",
            anchor_nanomolar: 79.433,
            evidence_count: 29,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 pChEMBL median; assay-description filtered).",
        },
        PotencyAnchor {
            inhibitor_name: "sulbactam",
            enzyme_name: "TEM-1",
            anchor_nanomolar: 1412.538,
            evidence_count: 13,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 pChEMBL median; assay-description filtered).",
        },
        PotencyAnchor {
            inhibitor_name: "tazobactam",
            enzyme_name: "TEM-1",
            anchor_nanomolar: 113.501,
            evidence_count: 30,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 pChEMBL median; assay-description filtered).",
        },
        PotencyAnchor {
            inhibitor_name: "avibactam",
            enzyme_name: "TEM-1",
            anchor_nanomolar: 7.943,
            evidence_count: 9,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 pChEMBL median; assay-description filtered).",
        },
        PotencyAnchor {
            inhibitor_name: "vaborbactam",
            enzyme_name: "TEM-1",
            anchor_nanomolar: 922.571,
            evidence_count: 4,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 pChEMBL median; assay-description filtered).",
        },
        PotencyAnchor {
            inhibitor_name: "clavulanic_acid",
            enzyme_name: "KPC",
            anchor_nanomolar: 11_885.022,
            evidence_count: 6,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 pChEMBL median; assay-description filtered).",
        },
        PotencyAnchor {
            inhibitor_name: "sulbactam",
            enzyme_name: "KPC",
            anchor_nanomolar: 74_989.421,
            evidence_count: 2,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 pChEMBL median; assay-description filtered).",
        },
        PotencyAnchor {
            inhibitor_name: "tazobactam",
            enzyme_name: "KPC",
            anchor_nanomolar: 8_128.305,
            evidence_count: 10,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 pChEMBL median; assay-description filtered).",
        },
        PotencyAnchor {
            inhibitor_name: "avibactam",
            enzyme_name: "KPC",
            anchor_nanomolar: 60.256,
            evidence_count: 7,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 pChEMBL median; assay-description filtered).",
        },
        PotencyAnchor {
            inhibitor_name: "vaborbactam",
            enzyme_name: "KPC",
            anchor_nanomolar: 69.183,
            evidence_count: 7,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 pChEMBL median; assay-description filtered).",
        },
        PotencyAnchor {
            inhibitor_name: "clavulanic_acid",
            enzyme_name: "NDM-1",
            anchor_nanomolar: 100_000.0,
            evidence_count: 1,
            imputed: false,
            notes: "ChEMBL snapshot (IC50-only anchor in filtered records).",
        },
        PotencyAnchor {
            inhibitor_name: "sulbactam",
            enzyme_name: "NDM-1",
            anchor_nanomolar: 200_000.0,
            evidence_count: 0,
            imputed: true,
            notes: "No direct filtered record; conservative imputation from class-B weak-inhibition regime.",
        },
        PotencyAnchor {
            inhibitor_name: "tazobactam",
            enzyme_name: "NDM-1",
            anchor_nanomolar: 173_205.081,
            evidence_count: 2,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 geometric-center anchor in filtered records).",
        },
        PotencyAnchor {
            inhibitor_name: "avibactam",
            enzyme_name: "NDM-1",
            anchor_nanomolar: 63_245.553,
            evidence_count: 2,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 geometric-center anchor in filtered records).",
        },
        PotencyAnchor {
            inhibitor_name: "vaborbactam",
            enzyme_name: "NDM-1",
            anchor_nanomolar: 199_769.728,
            evidence_count: 4,
            imputed: false,
            notes: "ChEMBL snapshot (Ki/IC50 geometric-center anchor in filtered records).",
        },
    ]
}

pub fn delta_g_from_potency_nanomolar(potency_nanomolar: f64, temperature_k: f64) -> f64 {
    let c_molar = potency_nanomolar.max(1.0e-12) * 1.0e-9;
    (R_GAS_J_MOL_K * temperature_k.max(1.0) * c_molar.ln()) / 1000.0
}

pub fn potency_nanomolar_from_delta_g(delta_g_kj_mol: f64, temperature_k: f64) -> f64 {
    let exponent = delta_g_kj_mol * 1000.0 / (R_GAS_J_MOL_K * temperature_k.max(1.0));
    exponent.exp() * 1.0e9
}

fn qed_contact_energy_kj_mol(charge_product: f64, distance_nm: f64, dielectric: f64) -> f64 {
    let q = charge_product.abs();
    let r_m = distance_nm.max(1.0e-6) * 1.0e-9;
    let eps = dielectric.max(1.0);
    let per_molecule_j = -(q * ALPHA_LEADING_ORDER * HBAR * C) / (eps * r_m);
    per_molecule_j * AVOGADRO / 1000.0
}

fn occupancy(concentration_nanomolar: f64, potency_nanomolar: f64) -> f64 {
    let c = concentration_nanomolar.max(0.0);
    let k = potency_nanomolar.max(1.0e-12);
    c / (c + k)
}

fn pair_result(
    inhibitor: InhibitorSpec,
    enzyme: EnzymeSpec,
    anchor: PotencyAnchor,
    temperature_k: f64,
    c: ResistanceModelCoefficients,
) -> PairwiseResistanceResult {
    let ionic_contacts = inhibitor.anionic_sites.max(0.0)
        * (c.ionic_contact_scale * enzyme.serine_drive + 0.35 * enzyme.zinc_drive);
    let hbond_contacts =
        inhibitor.hbond_sites.max(0.0) * enzyme.hbond_density.max(0.0) * c.hbond_contact_scale;

    let ionic_floor = ionic_contacts
        * qed_contact_energy_kj_mol(1.0, enzyme.ionic_distance_nm, enzyme.active_site_dielectric);
    let hbond_floor = hbond_contacts
        * qed_contact_energy_kj_mol(
            c.hbond_charge_product,
            enzyme.hbond_distance_nm,
            enzyme.active_site_dielectric + 2.0,
        );
    let qed_total = ionic_floor + hbond_floor;

    let residual =
        -c.serine_trap_kj * inhibitor.serine_trap_strength.max(0.0) * enzyme.serine_drive
            - c.boronate_serine_kj
                * inhibitor.boronate_reversible_strength.max(0.0)
                * enzyme.boronate_drive
            - c.zinc_match_kj * inhibitor.zinc_chelation_strength.max(0.0) * enzyme.zinc_drive
            + c.zinc_mismatch_penalty_kj
                * enzyme.zinc_drive
                * (1.0 - inhibitor.zinc_chelation_strength.max(0.0)).clamp(0.0, 1.0)
            - c.hydrophobic_coeff_kj_per_a2
                * inhibitor.hydrophobic_surface_a2.max(0.0)
                * enzyme.steric_openness.max(0.0)
            + c.entropy_coeff_kj * inhibitor.flexibility_penalty.max(0.0)
            + c.desolv_coeff_kj * inhibitor.polar_desolvation_penalty.max(0.0)
            + match enzyme.class {
                BetaLactamaseClass::SerineClassA => c.offset_serine_class_kj,
                BetaLactamaseClass::MetalloClassB => c.offset_metallo_class_kj,
            };

    let predicted_delta_g = qed_total + residual;
    let predicted_nanomolar =
        potency_nanomolar_from_delta_g(predicted_delta_g, temperature_k).max(1.0e-9);
    let anchor_delta_g = delta_g_from_potency_nanomolar(anchor.anchor_nanomolar, temperature_k);
    let log10_error =
        (predicted_nanomolar.log10() - anchor.anchor_nanomolar.max(1.0e-12).log10()).abs();

    PairwiseResistanceResult {
        inhibitor_name: inhibitor.name,
        inhibitor_chembl_id: inhibitor.chembl_id,
        scaffold: inhibitor.scaffold,
        enzyme_name: enzyme.name,
        enzyme_chembl_hint: enzyme.chembl_hint,
        enzyme_class: enzyme.class,
        evidence_count: anchor.evidence_count,
        imputed_anchor: anchor.imputed,
        anchor_nanomolar: anchor.anchor_nanomolar,
        anchor_delta_g_kj_mol: anchor_delta_g,
        qed_ionic_floor_kj_mol: ionic_floor,
        qed_hbond_floor_kj_mol: hbond_floor,
        qed_floor_total_kj_mol: qed_total,
        residual_modeled_total_kj_mol: residual,
        predicted_delta_g_kj_mol: predicted_delta_g,
        predicted_nanomolar,
        log10_error_pred_vs_anchor: log10_error,
        occupancy_anchor_at_1u_m: occupancy(1000.0, anchor.anchor_nanomolar),
        occupancy_predicted_at_1u_m: occupancy(1000.0, predicted_nanomolar),
    }
}

fn best_for_enzyme(rows: &[PairwiseResistanceResult], enzyme_name: &str) -> Option<EnzymeBest> {
    let subset = rows.iter().filter(|r| r.enzyme_name == enzyme_name);
    let by_anchor = subset
        .clone()
        .min_by(|a, b| {
            a.anchor_nanomolar
                .partial_cmp(&b.anchor_nanomolar)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied()?;
    let by_pred = subset
        .min_by(|a, b| {
            a.predicted_nanomolar
                .partial_cmp(&b.predicted_nanomolar)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied()?;
    Some(EnzymeBest {
        enzyme_name: by_anchor.enzyme_name,
        by_anchor_inhibitor: by_anchor.inhibitor_name,
        by_anchor_nanomolar: by_anchor.anchor_nanomolar,
        by_predicted_inhibitor: by_pred.inhibitor_name,
        by_predicted_nanomolar: by_pred.predicted_nanomolar,
        predicted_match_anchor_winner: by_anchor.inhibitor_name == by_pred.inhibitor_name,
    })
}

pub fn evaluate_antibiotic_resistance_panel(
    inhibitors: &[InhibitorSpec],
    enzymes: &[EnzymeSpec],
    anchors: &[PotencyAnchor],
    temperature_k: f64,
    coeffs: ResistanceModelCoefficients,
) -> AntibioticResistancePanel {
    let mut rows = Vec::new();
    for inhibitor in inhibitors {
        for enzyme in enzymes {
            if let Some(anchor) = anchors
                .iter()
                .find(|a| a.inhibitor_name == inhibitor.name && a.enzyme_name == enzyme.name)
            {
                rows.push(pair_result(
                    *inhibitor,
                    *enzyme,
                    *anchor,
                    temperature_k,
                    coeffs,
                ));
            }
        }
    }

    let mean_abs_log10_error = rows
        .iter()
        .map(|r| r.log10_error_pred_vs_anchor.abs())
        .sum::<f64>()
        / rows.len().max(1) as f64;

    let ndm_max_predicted_occupancy_at_1u_m = rows
        .iter()
        .filter(|r| r.enzyme_name == "NDM-1")
        .map(|r| r.occupancy_predicted_at_1u_m)
        .fold(0.0_f64, f64::max);

    let best_by_enzyme = enzymes
        .iter()
        .filter_map(|e| best_for_enzyme(&rows, e.name))
        .collect::<Vec<_>>();

    AntibioticResistancePanel {
        rows,
        best_by_enzyme,
        mean_abs_log10_error,
        ndm_max_predicted_occupancy_at_1u_m,
    }
}

pub fn default_antibiotic_resistance_panel(temperature_k: f64) -> AntibioticResistancePanel {
    evaluate_antibiotic_resistance_panel(
        &default_beta_lactamase_inhibitors(),
        &default_beta_lactamase_enzymes(),
        &default_potency_anchors(),
        temperature_k,
        ResistanceModelCoefficients::default(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn panel_builds_full_matrix() {
        let panel = default_antibiotic_resistance_panel(310.15);
        assert_eq!(panel.rows.len(), 15);
    }

    #[test]
    fn tem_and_kpc_anchor_winners_are_avibactam() {
        let panel = default_antibiotic_resistance_panel(310.15);
        let tem = panel
            .best_by_enzyme
            .iter()
            .find(|r| r.enzyme_name == "TEM-1")
            .unwrap();
        let kpc = panel
            .best_by_enzyme
            .iter()
            .find(|r| r.enzyme_name == "KPC")
            .unwrap();
        assert_eq!(tem.by_anchor_inhibitor, "avibactam");
        assert_eq!(kpc.by_anchor_inhibitor, "avibactam");
    }

    #[test]
    fn ndm_predicted_occupancy_stays_low() {
        let panel = default_antibiotic_resistance_panel(310.15);
        assert!(
            panel.ndm_max_predicted_occupancy_at_1u_m < 0.10,
            "NDM occupancy should stay low at 1uM: {}",
            panel.ndm_max_predicted_occupancy_at_1u_m
        );
    }
}
