/*!
 * Cardiovascular molecular-binding transduction lane.
 *
 * This lane provides a first-principles electrostatic floor for a known
 * cardiometabolic drug-target pair (atorvastatin <-> HMG-CoA reductase),
 * alongside the experimental thermodynamic target from Ki.
 *
 * Scope and honesty:
 * - We compute the exact thermodynamic conversion Ki -> ΔG.
 * - We compute a QED electrostatic contact floor from α ħ c / (ε r).
 * - We do not claim full docking/QM closure yet; any remaining stabilization
 *   term is reported explicitly as residual.
 */

use crate::chemical_thermo::{AVOGADRO, R_GAS_J_MOL_K};
use crate::{ALPHA_LEADING_ORDER, C, HBAR};

#[derive(Clone, Copy, Debug)]
pub struct BindingBenchmarkInput {
    pub ki_nanomolar: f64,
    pub temperature_k: f64,
}

impl Default for BindingBenchmarkInput {
    fn default() -> Self {
        Self {
            // Widely reported ballpark for atorvastatin inhibition potency.
            ki_nanomolar: 8.0,
            // Ambient biochemical baseline for standard-state free-energy reporting.
            temperature_k: 298.15,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ElectrostaticProxyInput {
    /// Effective ionic-contact count (dominant charged contact motifs).
    pub ionic_contact_count: f64,
    /// Mean ionic contact distance (nm).
    pub ionic_distance_nm: f64,
    /// Effective protein dielectric for ionic contacts.
    pub ionic_dielectric: f64,
    /// Effective H-bond contact count.
    pub hbond_contact_count: f64,
    /// Effective partial-charge product |q1 q2| / e^2 for H-bond-like contacts.
    pub hbond_charge_product: f64,
    /// Mean H-bond contact distance (nm).
    pub hbond_distance_nm: f64,
    /// Effective dielectric for H-bond-like contacts.
    pub hbond_dielectric: f64,
}

impl Default for ElectrostaticProxyInput {
    fn default() -> Self {
        Self {
            ionic_contact_count: 1.0,
            ionic_distance_nm: 0.30,
            ionic_dielectric: 28.0,
            hbond_contact_count: 5.0,
            hbond_charge_product: 0.20,
            hbond_distance_nm: 0.29,
            hbond_dielectric: 24.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct BindingEnergyDecomposition {
    pub experimental_delta_g_kj_mol: f64,
    pub qed_ionic_floor_kj_mol: f64,
    pub qed_hbond_floor_kj_mol: f64,
    pub qed_floor_total_kj_mol: f64,
    pub residual_required_kj_mol: f64,
    pub explained_fraction_of_abs_delta_g: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct ResidualProxyInput {
    /// Effective non-polar buried area (A^2) contributing to hydrophobic stabilization.
    pub effective_hydrophobic_area_a2: f64,
    /// Hydrophobic free-energy coefficient (kJ/mol/A^2).
    pub hydrophobic_coeff_kj_per_a2: f64,
    /// Effective aromatic/dispersion packing contacts.
    pub aromatic_contact_count: f64,
    /// Stabilization per aromatic contact (kJ/mol).
    pub aromatic_contact_stabilization_kj: f64,
    /// Structured waters displaced upon binding.
    pub released_water_count: f64,
    /// Entropic stabilization per released water (kJ/mol).
    pub water_release_stabilization_kj: f64,
    /// Rotatable bonds constrained on binding.
    pub constrained_rotatable_bonds: f64,
    /// Entropic penalty per constrained rotor (kJ/mol).
    pub conformational_entropy_penalty_per_rotor_kj: f64,
    /// Polar contacts with incomplete compensation (desolvation cost).
    pub polar_desolvated_contact_count: f64,
    /// Penalty per polar desolvated contact (kJ/mol).
    pub polar_desolvation_penalty_kj: f64,
    /// Residual strain/reorganization penalty (kJ/mol).
    pub ligand_strain_penalty_kj: f64,
}

impl Default for ResidualProxyInput {
    fn default() -> Self {
        Self {
            effective_hydrophobic_area_a2: 225.0,
            hydrophobic_coeff_kj_per_a2: 0.046,
            aromatic_contact_count: 2.5,
            aromatic_contact_stabilization_kj: 1.10,
            released_water_count: 3.0,
            water_release_stabilization_kj: 1.05,
            constrained_rotatable_bonds: 6.0,
            conformational_entropy_penalty_per_rotor_kj: 0.78,
            polar_desolvated_contact_count: 3.0,
            polar_desolvation_penalty_kj: 0.45,
            ligand_strain_penalty_kj: 0.60,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ResidualEnergyBreakdown {
    pub hydrophobic_stabilization_kj_mol: f64,
    pub aromatic_packing_stabilization_kj_mol: f64,
    pub water_release_stabilization_kj_mol: f64,
    pub conformational_entropy_penalty_kj_mol: f64,
    pub polar_desolvation_penalty_kj_mol: f64,
    pub ligand_strain_penalty_kj_mol: f64,
    pub modeled_residual_total_kj_mol: f64,
    pub target_residual_kj_mol: f64,
    pub closure_error_kj_mol: f64,
}

/// Exact standard-state free-energy transduction from Ki:
/// ΔG = R T ln(Ki / 1 M), with Ki in mol/L.
pub fn delta_g_from_ki_nanomolar(ki_nanomolar: f64, temperature_k: f64) -> f64 {
    let ki_molar = (ki_nanomolar.max(1.0e-18)) * 1.0e-9;
    let t = temperature_k.max(1.0);
    (R_GAS_J_MOL_K * t * ki_molar.ln()) / 1000.0
}

/// QED electrostatic pair energy:
/// E_pair = -(q1 q2) α ħ c / (ε r), converted to kJ/mol.
pub fn qed_contact_energy_kj_mol(charge_product: f64, distance_nm: f64, dielectric: f64) -> f64 {
    let q = charge_product.abs();
    let r_m = (distance_nm.max(1.0e-6)) * 1.0e-9;
    let eps = dielectric.max(1.0);
    let per_molecule_j = -(q * ALPHA_LEADING_ORDER * HBAR * C) / (eps * r_m);
    per_molecule_j * AVOGADRO / 1000.0
}

pub fn evaluate_atorvastatin_hmgcr_binding(
    benchmark: BindingBenchmarkInput,
    proxy: ElectrostaticProxyInput,
) -> BindingEnergyDecomposition {
    let experimental = delta_g_from_ki_nanomolar(benchmark.ki_nanomolar, benchmark.temperature_k);

    let ionic = proxy.ionic_contact_count.max(0.0)
        * qed_contact_energy_kj_mol(1.0, proxy.ionic_distance_nm, proxy.ionic_dielectric);
    let hbond = proxy.hbond_contact_count.max(0.0)
        * qed_contact_energy_kj_mol(
            proxy.hbond_charge_product,
            proxy.hbond_distance_nm,
            proxy.hbond_dielectric,
        );
    let floor_total = ionic + hbond;

    let residual_required = experimental - floor_total;
    let explained = (floor_total.abs() / experimental.abs().max(1.0e-12)).clamp(0.0, 1.0);

    BindingEnergyDecomposition {
        experimental_delta_g_kj_mol: experimental,
        qed_ionic_floor_kj_mol: ionic,
        qed_hbond_floor_kj_mol: hbond,
        qed_floor_total_kj_mol: floor_total,
        residual_required_kj_mol: residual_required,
        explained_fraction_of_abs_delta_g: explained,
    }
}

pub fn decompose_non_electrostatic_residual(
    target_residual_kj_mol: f64,
    proxy: ResidualProxyInput,
) -> ResidualEnergyBreakdown {
    let hydrophobic =
        -(proxy.effective_hydrophobic_area_a2.max(0.0) * proxy.hydrophobic_coeff_kj_per_a2.max(0.0));
    let aromatic = -(proxy.aromatic_contact_count.max(0.0)
        * proxy.aromatic_contact_stabilization_kj.max(0.0));
    let water = -(proxy.released_water_count.max(0.0) * proxy.water_release_stabilization_kj.max(0.0));
    let entropy = proxy.constrained_rotatable_bonds.max(0.0)
        * proxy
            .conformational_entropy_penalty_per_rotor_kj
            .max(0.0);
    let desolvation =
        proxy.polar_desolvated_contact_count.max(0.0) * proxy.polar_desolvation_penalty_kj.max(0.0);
    let strain = proxy.ligand_strain_penalty_kj.max(0.0);

    let modeled_total = hydrophobic + aromatic + water + entropy + desolvation + strain;
    let closure_error = modeled_total - target_residual_kj_mol;

    ResidualEnergyBreakdown {
        hydrophobic_stabilization_kj_mol: hydrophobic,
        aromatic_packing_stabilization_kj_mol: aromatic,
        water_release_stabilization_kj_mol: water,
        conformational_entropy_penalty_kj_mol: entropy,
        polar_desolvation_penalty_kj_mol: desolvation,
        ligand_strain_penalty_kj_mol: strain,
        modeled_residual_total_kj_mol: modeled_total,
        target_residual_kj_mol,
        closure_error_kj_mol: closure_error,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ki_to_delta_g_matches_known_ballpark() {
        let dg = delta_g_from_ki_nanomolar(8.0, 298.15);
        assert!(dg < -45.0, "expected strong binding, got {dg}");
        assert!(dg > -48.0, "expected known 8 nM ballpark, got {dg}");
    }

    #[test]
    fn qed_contact_energy_scales_with_distance_and_dielectric() {
        let near = qed_contact_energy_kj_mol(1.0, 0.30, 20.0);
        let far = qed_contact_energy_kj_mol(1.0, 0.60, 20.0);
        let screened = qed_contact_energy_kj_mol(1.0, 0.30, 40.0);
        assert!(near < far, "near contact should be more stabilizing");
        assert!(near < screened, "lower dielectric should strengthen attraction");
    }

    #[test]
    fn atorvastatin_decomposition_is_physically_ordered() {
        let score =
            evaluate_atorvastatin_hmgcr_binding(BindingBenchmarkInput::default(), ElectrostaticProxyInput::default());
        assert!(score.experimental_delta_g_kj_mol < 0.0);
        assert!(score.qed_floor_total_kj_mol < 0.0);
        assert!(score.explained_fraction_of_abs_delta_g > 0.5);
        assert!(score.explained_fraction_of_abs_delta_g < 1.0);
        // Residual remains negative: extra stabilization is still needed beyond the floor.
        assert!(score.residual_required_kj_mol < 0.0);
    }

    #[test]
    fn residual_decomposition_has_expected_term_signs() {
        let residual = decompose_non_electrostatic_residual(-9.706, ResidualProxyInput::default());
        assert!(residual.hydrophobic_stabilization_kj_mol < 0.0);
        assert!(residual.aromatic_packing_stabilization_kj_mol < 0.0);
        assert!(residual.water_release_stabilization_kj_mol < 0.0);
        assert!(residual.conformational_entropy_penalty_kj_mol > 0.0);
        assert!(residual.polar_desolvation_penalty_kj_mol > 0.0);
        assert!(residual.ligand_strain_penalty_kj_mol > 0.0);
    }

    #[test]
    fn residual_default_proxies_close_most_of_required_gap() {
        let score = evaluate_atorvastatin_hmgcr_binding(
            BindingBenchmarkInput::default(),
            ElectrostaticProxyInput::default(),
        );
        let residual =
            decompose_non_electrostatic_residual(score.residual_required_kj_mol, ResidualProxyInput::default());
        assert!(
            residual.closure_error_kj_mol.abs() < 2.0,
            "residual proxy should be close-order accurate: {residual:?}"
        );
    }
}
