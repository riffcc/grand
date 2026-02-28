/*!
 * THC-CB1 neuronal response lane.
 *
 * Scope:
 * - Binding thermodynamics for THC <-> CB1 from Ki.
 * - QED electrostatic floor contribution.
 * - Explicit non-electrostatic residual decomposition.
 * - Neuronal response transduction via CB1 occupancy.
 *
 * This is a reduced biophysical lane, not a full receptor MD/QM solver.
 */

use crate::cardiovascular_binding::{delta_g_from_ki_nanomolar, qed_contact_energy_kj_mol};

#[derive(Clone, Copy, Debug)]
pub struct ThcCb1BindingInput {
    pub ki_nanomolar: f64,
    pub temperature_k: f64,
}

pub type Cb1LigandBindingInput = ThcCb1BindingInput;

impl Default for ThcCb1BindingInput {
    fn default() -> Self {
        Self {
            // Tunable literature-scale potency baseline for THC at CB1.
            ki_nanomolar: 40.0,
            // Neuronal physiological temperature.
            temperature_k: 310.15,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ThcElectrostaticProxyInput {
    pub hbond_contact_count: f64,
    pub hbond_charge_product: f64,
    pub hbond_distance_nm: f64,
    pub hbond_dielectric: f64,
    pub polar_dipole_contact_count: f64,
    pub polar_dipole_charge_product: f64,
    pub polar_dipole_distance_nm: f64,
    pub polar_dipole_dielectric: f64,
}

pub type Cb1ElectrostaticProxyInput = ThcElectrostaticProxyInput;

impl Default for ThcElectrostaticProxyInput {
    fn default() -> Self {
        Self {
            hbond_contact_count: 1.4,
            hbond_charge_product: 0.16,
            hbond_distance_nm: 0.30,
            hbond_dielectric: 28.0,
            polar_dipole_contact_count: 1.2,
            polar_dipole_charge_product: 0.10,
            polar_dipole_distance_nm: 0.34,
            polar_dipole_dielectric: 32.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ThcCb1BindingScore {
    pub experimental_delta_g_kj_mol: f64,
    pub qed_hbond_floor_kj_mol: f64,
    pub qed_polar_floor_kj_mol: f64,
    pub qed_floor_total_kj_mol: f64,
    pub residual_required_kj_mol: f64,
    pub explained_fraction_of_abs_delta_g: f64,
}

pub type Cb1BindingScore = ThcCb1BindingScore;

#[derive(Clone, Copy, Debug)]
pub struct ThcResidualProxyInput {
    pub effective_hydrophobic_area_a2: f64,
    pub hydrophobic_coeff_kj_per_a2: f64,
    pub aromatic_contact_count: f64,
    pub aromatic_contact_stabilization_kj: f64,
    pub released_water_count: f64,
    pub water_release_stabilization_kj: f64,
    pub constrained_rotatable_bonds: f64,
    pub conformational_entropy_penalty_per_rotor_kj: f64,
    pub polar_desolvated_contact_count: f64,
    pub polar_desolvation_penalty_kj: f64,
    pub ligand_strain_penalty_kj: f64,
}

pub type Cb1ResidualProxyInput = ThcResidualProxyInput;

impl Default for ThcResidualProxyInput {
    fn default() -> Self {
        Self {
            effective_hydrophobic_area_a2: 700.0,
            hydrophobic_coeff_kj_per_a2: 0.052,
            aromatic_contact_count: 3.0,
            aromatic_contact_stabilization_kj: 1.55,
            released_water_count: 2.8,
            water_release_stabilization_kj: 1.20,
            constrained_rotatable_bonds: 5.5,
            conformational_entropy_penalty_per_rotor_kj: 0.60,
            polar_desolvated_contact_count: 1.4,
            polar_desolvation_penalty_kj: 0.55,
            ligand_strain_penalty_kj: 0.75,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ThcResidualBreakdown {
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

pub type Cb1ResidualBreakdown = ThcResidualBreakdown;

#[derive(Clone, Copy, Debug)]
pub struct NeuronCouplingInput {
    pub intrinsic_efficacy: f64,
    pub max_release_inhibition_fraction: f64,
    pub max_firing_suppression_fraction: f64,
    pub hill_coefficient: f64,
    pub baseline_release_probability: f64,
    pub baseline_firing_rate_hz: f64,
}

impl Default for NeuronCouplingInput {
    fn default() -> Self {
        Self {
            intrinsic_efficacy: 0.55,
            max_release_inhibition_fraction: 0.75,
            max_firing_suppression_fraction: 0.45,
            hill_coefficient: 1.0,
            baseline_release_probability: 0.35,
            baseline_firing_rate_hz: 8.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct NeuronResponsePoint {
    pub concentration_nanomolar: f64,
    pub occupancy_fraction: f64,
    pub effective_activation_fraction: f64,
    pub release_probability: f64,
    pub firing_rate_hz: f64,
    pub relative_firing_scale: f64,
}

pub fn thc_cb1_occupancy_fraction(
    concentration_nanomolar: f64,
    ki_nanomolar: f64,
    hill_coefficient: f64,
) -> f64 {
    let c = concentration_nanomolar.max(0.0);
    let ki = ki_nanomolar.max(1.0e-12);
    let n = hill_coefficient.clamp(0.2, 3.0);
    if c <= 0.0 {
        return 0.0;
    }
    let c_n = c.powf(n);
    let ki_n = ki.powf(n);
    (c_n / (c_n + ki_n)).clamp(0.0, 1.0)
}

pub fn evaluate_thc_cb1_binding(
    binding: ThcCb1BindingInput,
    electro: ThcElectrostaticProxyInput,
) -> ThcCb1BindingScore {
    let experimental = delta_g_from_ki_nanomolar(binding.ki_nanomolar, binding.temperature_k);

    let hbond = electro.hbond_contact_count.max(0.0)
        * qed_contact_energy_kj_mol(
            electro.hbond_charge_product,
            electro.hbond_distance_nm,
            electro.hbond_dielectric,
        );
    let polar = electro.polar_dipole_contact_count.max(0.0)
        * qed_contact_energy_kj_mol(
            electro.polar_dipole_charge_product,
            electro.polar_dipole_distance_nm,
            electro.polar_dipole_dielectric,
        );
    let floor_total = hbond + polar;
    let residual_required = experimental - floor_total;
    let explained = (floor_total.abs() / experimental.abs().max(1.0e-12)).clamp(0.0, 1.0);

    ThcCb1BindingScore {
        experimental_delta_g_kj_mol: experimental,
        qed_hbond_floor_kj_mol: hbond,
        qed_polar_floor_kj_mol: polar,
        qed_floor_total_kj_mol: floor_total,
        residual_required_kj_mol: residual_required,
        explained_fraction_of_abs_delta_g: explained,
    }
}

pub fn evaluate_cb1_ligand_binding(
    binding: Cb1LigandBindingInput,
    electro: Cb1ElectrostaticProxyInput,
) -> Cb1BindingScore {
    evaluate_thc_cb1_binding(binding, electro)
}

pub fn decompose_thc_cb1_non_electrostatic_residual(
    target_residual_kj_mol: f64,
    proxy: ThcResidualProxyInput,
) -> ThcResidualBreakdown {
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

    ThcResidualBreakdown {
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

pub fn decompose_cb1_non_electrostatic_residual(
    target_residual_kj_mol: f64,
    proxy: Cb1ResidualProxyInput,
) -> Cb1ResidualBreakdown {
    decompose_thc_cb1_non_electrostatic_residual(target_residual_kj_mol, proxy)
}

pub fn simulate_thc_cb1_neuron_response(
    binding: ThcCb1BindingInput,
    coupling: NeuronCouplingInput,
    concentrations_nanomolar: &[f64],
) -> Vec<NeuronResponsePoint> {
    concentrations_nanomolar
        .iter()
        .map(|&c_nm| {
            let occupancy =
                thc_cb1_occupancy_fraction(c_nm, binding.ki_nanomolar, coupling.hill_coefficient);
            let activation = (occupancy * coupling.intrinsic_efficacy).clamp(0.0, 1.0);

            let release_scale =
                (1.0 - coupling.max_release_inhibition_fraction.clamp(0.0, 1.0) * activation)
                    .clamp(0.0, 1.0);
            let firing_scale =
                (1.0 - coupling.max_firing_suppression_fraction.clamp(0.0, 1.0) * activation)
                    .clamp(0.0, 1.0);

            NeuronResponsePoint {
                concentration_nanomolar: c_nm.max(0.0),
                occupancy_fraction: occupancy,
                effective_activation_fraction: activation,
                release_probability: coupling.baseline_release_probability.max(0.0) * release_scale,
                firing_rate_hz: coupling.baseline_firing_rate_hz.max(0.0) * firing_scale,
                relative_firing_scale: firing_scale,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thc_ki_to_delta_g_is_reasonable() {
        let binding = ThcCb1BindingInput::default();
        let score = evaluate_thc_cb1_binding(binding, ThcElectrostaticProxyInput::default());
        assert!(score.experimental_delta_g_kj_mol < -40.0);
        assert!(score.experimental_delta_g_kj_mol > -48.0);
    }

    #[test]
    fn occupancy_is_monotone_in_concentration() {
        let binding = ThcCb1BindingInput::default();
        let coupling = NeuronCouplingInput::default();
        let xs = [0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0];
        let pts = simulate_thc_cb1_neuron_response(binding, coupling, &xs);
        for i in 1..pts.len() {
            assert!(pts[i].occupancy_fraction >= pts[i - 1].occupancy_fraction);
            assert!(pts[i].release_probability <= pts[i - 1].release_probability + 1.0e-12);
            assert!(pts[i].firing_rate_hz <= pts[i - 1].firing_rate_hz + 1.0e-12);
        }
    }

    #[test]
    fn thc_residual_proxy_closes_majority_of_gap() {
        let score = evaluate_thc_cb1_binding(
            ThcCb1BindingInput::default(),
            ThcElectrostaticProxyInput::default(),
        );
        let residual = decompose_thc_cb1_non_electrostatic_residual(
            score.residual_required_kj_mol,
            ThcResidualProxyInput::default(),
        );
        assert!(residual.modeled_residual_total_kj_mol < -30.0);
        assert!(residual.closure_error_kj_mol.abs() < 3.0);
    }
}
