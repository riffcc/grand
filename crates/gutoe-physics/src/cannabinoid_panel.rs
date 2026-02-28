/*!
 * Multi-cannabinoid CB1 panel built on the reduced biophysical lane.
 *
 * This is a comparative inference layer over:
 * - Ki -> ΔG thermodynamics
 * - QED electrostatic floor
 * - residual decomposition
 * - occupancy-to-neuron suppression mapping
 */

use crate::{
    decompose_cb1_non_electrostatic_residual, evaluate_cb1_ligand_binding,
    simulate_thc_cb1_neuron_response, Cb1ElectrostaticProxyInput, Cb1LigandBindingInput,
    Cb1ResidualProxyInput, NeuronCouplingInput,
};

#[derive(Clone, Copy, Debug)]
pub struct CannabinoidSpec {
    pub name: &'static str,
    pub class: &'static str,
    pub ki_cb1_nanomolar: f64,
    pub ki_cb2_nanomolar: f64,
    pub intrinsic_efficacy_cb1: f64,
    pub electro: Cb1ElectrostaticProxyInput,
    pub residual: Cb1ResidualProxyInput,
}

#[derive(Clone, Copy, Debug)]
pub struct CannabinoidPanelResult {
    pub name: &'static str,
    pub class: &'static str,
    pub ki_cb1_nanomolar: f64,
    pub ki_cb2_nanomolar: f64,
    pub intrinsic_efficacy_cb1: f64,
    pub experimental_delta_g_kj_mol: f64,
    pub qed_floor_total_kj_mol: f64,
    pub residual_required_kj_mol: f64,
    pub residual_modeled_total_kj_mol: f64,
    pub residual_closure_error_kj_mol: f64,
    pub explained_fraction_of_abs_delta_g: f64,
    pub occupancy_10nm: f64,
    pub occupancy_30nm: f64,
    pub occupancy_100nm: f64,
    pub firing_scale_10nm: f64,
    pub firing_scale_30nm: f64,
    pub firing_scale_100nm: f64,
}

fn nearest_point(
    points: &[crate::NeuronResponsePoint],
    target_nm: f64,
) -> Option<crate::NeuronResponsePoint> {
    points.iter().copied().min_by(|a, b| {
        let da = (a.concentration_nanomolar - target_nm).abs();
        let db = (b.concentration_nanomolar - target_nm).abs();
        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
    })
}

pub fn default_cannabinoid_specs() -> Vec<CannabinoidSpec> {
    vec![
        CannabinoidSpec {
            name: "delta9_thc",
            class: "phytocannabinoid",
            ki_cb1_nanomolar: 40.0,
            ki_cb2_nanomolar: 36.0,
            intrinsic_efficacy_cb1: 0.55,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 1.4,
                hbond_charge_product: 0.16,
                hbond_distance_nm: 0.30,
                hbond_dielectric: 28.0,
                polar_dipole_contact_count: 1.2,
                polar_dipole_charge_product: 0.10,
                polar_dipole_distance_nm: 0.34,
                polar_dipole_dielectric: 32.0,
            },
            residual: Cb1ResidualProxyInput {
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
            },
        },
        CannabinoidSpec {
            name: "delta8_thc",
            class: "phytocannabinoid",
            ki_cb1_nanomolar: 44.0,
            ki_cb2_nanomolar: 44.0,
            intrinsic_efficacy_cb1: 0.50,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 1.3,
                hbond_charge_product: 0.15,
                hbond_distance_nm: 0.30,
                hbond_dielectric: 28.0,
                polar_dipole_contact_count: 1.1,
                polar_dipole_charge_product: 0.10,
                polar_dipole_distance_nm: 0.34,
                polar_dipole_dielectric: 32.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 690.0,
                hydrophobic_coeff_kj_per_a2: 0.052,
                aromatic_contact_count: 3.0,
                aromatic_contact_stabilization_kj: 1.50,
                released_water_count: 2.7,
                water_release_stabilization_kj: 1.20,
                constrained_rotatable_bonds: 5.4,
                conformational_entropy_penalty_per_rotor_kj: 0.60,
                polar_desolvated_contact_count: 1.3,
                polar_desolvation_penalty_kj: 0.55,
                ligand_strain_penalty_kj: 0.75,
            },
        },
        CannabinoidSpec {
            name: "11_oh_thc",
            class: "metabolite",
            ki_cb1_nanomolar: 20.0,
            ki_cb2_nanomolar: 25.0,
            intrinsic_efficacy_cb1: 0.70,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 2.0,
                hbond_charge_product: 0.18,
                hbond_distance_nm: 0.29,
                hbond_dielectric: 27.0,
                polar_dipole_contact_count: 1.6,
                polar_dipole_charge_product: 0.11,
                polar_dipole_distance_nm: 0.33,
                polar_dipole_dielectric: 31.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 660.0,
                hydrophobic_coeff_kj_per_a2: 0.051,
                aromatic_contact_count: 3.0,
                aromatic_contact_stabilization_kj: 1.50,
                released_water_count: 3.0,
                water_release_stabilization_kj: 1.20,
                constrained_rotatable_bonds: 5.8,
                conformational_entropy_penalty_per_rotor_kj: 0.60,
                polar_desolvated_contact_count: 1.6,
                polar_desolvation_penalty_kj: 0.55,
                ligand_strain_penalty_kj: 0.80,
            },
        },
        CannabinoidSpec {
            name: "cbd",
            class: "phytocannabinoid",
            ki_cb1_nanomolar: 4350.0,
            ki_cb2_nanomolar: 2860.0,
            intrinsic_efficacy_cb1: 0.05,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 1.7,
                hbond_charge_product: 0.17,
                hbond_distance_nm: 0.30,
                hbond_dielectric: 29.0,
                polar_dipole_contact_count: 1.6,
                polar_dipole_charge_product: 0.12,
                polar_dipole_distance_nm: 0.34,
                polar_dipole_dielectric: 32.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 720.0,
                hydrophobic_coeff_kj_per_a2: 0.050,
                aromatic_contact_count: 2.8,
                aromatic_contact_stabilization_kj: 1.40,
                released_water_count: 2.9,
                water_release_stabilization_kj: 1.10,
                constrained_rotatable_bonds: 6.8,
                conformational_entropy_penalty_per_rotor_kj: 0.62,
                polar_desolvated_contact_count: 1.9,
                polar_desolvation_penalty_kj: 0.55,
                ligand_strain_penalty_kj: 0.70,
            },
        },
        CannabinoidSpec {
            name: "cbn",
            class: "phytocannabinoid",
            ki_cb1_nanomolar: 310.0,
            ki_cb2_nanomolar: 126.0,
            intrinsic_efficacy_cb1: 0.30,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 1.2,
                hbond_charge_product: 0.15,
                hbond_distance_nm: 0.31,
                hbond_dielectric: 29.0,
                polar_dipole_contact_count: 1.1,
                polar_dipole_charge_product: 0.10,
                polar_dipole_distance_nm: 0.35,
                polar_dipole_dielectric: 32.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 690.0,
                hydrophobic_coeff_kj_per_a2: 0.051,
                aromatic_contact_count: 3.1,
                aromatic_contact_stabilization_kj: 1.50,
                released_water_count: 2.6,
                water_release_stabilization_kj: 1.15,
                constrained_rotatable_bonds: 5.2,
                conformational_entropy_penalty_per_rotor_kj: 0.60,
                polar_desolvated_contact_count: 1.2,
                polar_desolvation_penalty_kj: 0.54,
                ligand_strain_penalty_kj: 0.70,
            },
        },
        CannabinoidSpec {
            name: "cbg",
            class: "phytocannabinoid",
            ki_cb1_nanomolar: 381.0,
            ki_cb2_nanomolar: 260.0,
            intrinsic_efficacy_cb1: 0.25,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 1.5,
                hbond_charge_product: 0.16,
                hbond_distance_nm: 0.30,
                hbond_dielectric: 29.0,
                polar_dipole_contact_count: 1.3,
                polar_dipole_charge_product: 0.11,
                polar_dipole_distance_nm: 0.34,
                polar_dipole_dielectric: 32.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 700.0,
                hydrophobic_coeff_kj_per_a2: 0.051,
                aromatic_contact_count: 2.8,
                aromatic_contact_stabilization_kj: 1.45,
                released_water_count: 2.7,
                water_release_stabilization_kj: 1.15,
                constrained_rotatable_bonds: 6.0,
                conformational_entropy_penalty_per_rotor_kj: 0.60,
                polar_desolvated_contact_count: 1.5,
                polar_desolvation_penalty_kj: 0.55,
                ligand_strain_penalty_kj: 0.70,
            },
        },
        CannabinoidSpec {
            name: "cbc",
            class: "phytocannabinoid",
            ki_cb1_nanomolar: 1200.0,
            ki_cb2_nanomolar: 700.0,
            intrinsic_efficacy_cb1: 0.20,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 1.3,
                hbond_charge_product: 0.15,
                hbond_distance_nm: 0.31,
                hbond_dielectric: 29.0,
                polar_dipole_contact_count: 1.2,
                polar_dipole_charge_product: 0.10,
                polar_dipole_distance_nm: 0.35,
                polar_dipole_dielectric: 32.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 680.0,
                hydrophobic_coeff_kj_per_a2: 0.050,
                aromatic_contact_count: 2.6,
                aromatic_contact_stabilization_kj: 1.40,
                released_water_count: 2.5,
                water_release_stabilization_kj: 1.10,
                constrained_rotatable_bonds: 5.8,
                conformational_entropy_penalty_per_rotor_kj: 0.60,
                polar_desolvated_contact_count: 1.3,
                polar_desolvation_penalty_kj: 0.54,
                ligand_strain_penalty_kj: 0.70,
            },
        },
        CannabinoidSpec {
            name: "thcv",
            class: "phytocannabinoid",
            ki_cb1_nanomolar: 75.0,
            ki_cb2_nanomolar: 63.0,
            intrinsic_efficacy_cb1: 0.15,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 1.3,
                hbond_charge_product: 0.16,
                hbond_distance_nm: 0.30,
                hbond_dielectric: 28.0,
                polar_dipole_contact_count: 1.1,
                polar_dipole_charge_product: 0.10,
                polar_dipole_distance_nm: 0.34,
                polar_dipole_dielectric: 32.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 640.0,
                hydrophobic_coeff_kj_per_a2: 0.052,
                aromatic_contact_count: 2.7,
                aromatic_contact_stabilization_kj: 1.45,
                released_water_count: 2.4,
                water_release_stabilization_kj: 1.15,
                constrained_rotatable_bonds: 4.8,
                conformational_entropy_penalty_per_rotor_kj: 0.60,
                polar_desolvated_contact_count: 1.2,
                polar_desolvation_penalty_kj: 0.54,
                ligand_strain_penalty_kj: 0.70,
            },
        },
        CannabinoidSpec {
            name: "cbdv",
            class: "phytocannabinoid",
            ki_cb1_nanomolar: 950.0,
            ki_cb2_nanomolar: 1400.0,
            intrinsic_efficacy_cb1: 0.08,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 1.6,
                hbond_charge_product: 0.17,
                hbond_distance_nm: 0.30,
                hbond_dielectric: 29.0,
                polar_dipole_contact_count: 1.5,
                polar_dipole_charge_product: 0.11,
                polar_dipole_distance_nm: 0.34,
                polar_dipole_dielectric: 32.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 690.0,
                hydrophobic_coeff_kj_per_a2: 0.050,
                aromatic_contact_count: 2.6,
                aromatic_contact_stabilization_kj: 1.35,
                released_water_count: 2.7,
                water_release_stabilization_kj: 1.10,
                constrained_rotatable_bonds: 6.3,
                conformational_entropy_penalty_per_rotor_kj: 0.61,
                polar_desolvated_contact_count: 1.7,
                polar_desolvation_penalty_kj: 0.55,
                ligand_strain_penalty_kj: 0.70,
            },
        },
        CannabinoidSpec {
            name: "thca",
            class: "acidic_phytocannabinoid",
            ki_cb1_nanomolar: 630.0,
            ki_cb2_nanomolar: 1400.0,
            intrinsic_efficacy_cb1: 0.10,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 2.2,
                hbond_charge_product: 0.19,
                hbond_distance_nm: 0.29,
                hbond_dielectric: 27.0,
                polar_dipole_contact_count: 2.0,
                polar_dipole_charge_product: 0.14,
                polar_dipole_distance_nm: 0.33,
                polar_dipole_dielectric: 31.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 620.0,
                hydrophobic_coeff_kj_per_a2: 0.049,
                aromatic_contact_count: 2.7,
                aromatic_contact_stabilization_kj: 1.35,
                released_water_count: 3.4,
                water_release_stabilization_kj: 1.15,
                constrained_rotatable_bonds: 6.9,
                conformational_entropy_penalty_per_rotor_kj: 0.62,
                polar_desolvated_contact_count: 2.4,
                polar_desolvation_penalty_kj: 0.60,
                ligand_strain_penalty_kj: 0.80,
            },
        },
        CannabinoidSpec {
            name: "cbda",
            class: "acidic_phytocannabinoid",
            ki_cb1_nanomolar: 5000.0,
            ki_cb2_nanomolar: 10000.0,
            intrinsic_efficacy_cb1: 0.05,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 2.4,
                hbond_charge_product: 0.20,
                hbond_distance_nm: 0.29,
                hbond_dielectric: 27.0,
                polar_dipole_contact_count: 2.2,
                polar_dipole_charge_product: 0.14,
                polar_dipole_distance_nm: 0.33,
                polar_dipole_dielectric: 31.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 610.0,
                hydrophobic_coeff_kj_per_a2: 0.048,
                aromatic_contact_count: 2.5,
                aromatic_contact_stabilization_kj: 1.30,
                released_water_count: 3.6,
                water_release_stabilization_kj: 1.15,
                constrained_rotatable_bonds: 7.1,
                conformational_entropy_penalty_per_rotor_kj: 0.62,
                polar_desolvated_contact_count: 2.6,
                polar_desolvation_penalty_kj: 0.60,
                ligand_strain_penalty_kj: 0.80,
            },
        },
        CannabinoidSpec {
            name: "anandamide",
            class: "endocannabinoid",
            ki_cb1_nanomolar: 61.0,
            ki_cb2_nanomolar: 193.0,
            intrinsic_efficacy_cb1: 0.80,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 1.8,
                hbond_charge_product: 0.17,
                hbond_distance_nm: 0.30,
                hbond_dielectric: 28.0,
                polar_dipole_contact_count: 1.8,
                polar_dipole_charge_product: 0.12,
                polar_dipole_distance_nm: 0.34,
                polar_dipole_dielectric: 32.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 540.0,
                hydrophobic_coeff_kj_per_a2: 0.052,
                aromatic_contact_count: 0.4,
                aromatic_contact_stabilization_kj: 1.10,
                released_water_count: 2.6,
                water_release_stabilization_kj: 1.20,
                constrained_rotatable_bonds: 10.0,
                conformational_entropy_penalty_per_rotor_kj: 0.60,
                polar_desolvated_contact_count: 1.8,
                polar_desolvation_penalty_kj: 0.57,
                ligand_strain_penalty_kj: 0.85,
            },
        },
        CannabinoidSpec {
            name: "2_ag",
            class: "endocannabinoid",
            ki_cb1_nanomolar: 450.0,
            ki_cb2_nanomolar: 1400.0,
            intrinsic_efficacy_cb1: 0.95,
            electro: Cb1ElectrostaticProxyInput {
                hbond_contact_count: 1.9,
                hbond_charge_product: 0.17,
                hbond_distance_nm: 0.30,
                hbond_dielectric: 28.0,
                polar_dipole_contact_count: 1.9,
                polar_dipole_charge_product: 0.12,
                polar_dipole_distance_nm: 0.34,
                polar_dipole_dielectric: 32.0,
            },
            residual: Cb1ResidualProxyInput {
                effective_hydrophobic_area_a2: 510.0,
                hydrophobic_coeff_kj_per_a2: 0.051,
                aromatic_contact_count: 0.2,
                aromatic_contact_stabilization_kj: 1.00,
                released_water_count: 2.4,
                water_release_stabilization_kj: 1.20,
                constrained_rotatable_bonds: 9.0,
                conformational_entropy_penalty_per_rotor_kj: 0.60,
                polar_desolvated_contact_count: 1.7,
                polar_desolvation_penalty_kj: 0.57,
                ligand_strain_penalty_kj: 0.85,
            },
        },
    ]
}

pub fn evaluate_cannabinoid_panel(
    specs: &[CannabinoidSpec],
    temperature_k: f64,
    coupling_template: NeuronCouplingInput,
) -> Vec<CannabinoidPanelResult> {
    let concs = [0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0];
    specs
        .iter()
        .map(|spec| {
            let binding = Cb1LigandBindingInput {
                ki_nanomolar: spec.ki_cb1_nanomolar,
                temperature_k,
            };
            let score = evaluate_cb1_ligand_binding(binding, spec.electro);
            let residual =
                decompose_cb1_non_electrostatic_residual(score.residual_required_kj_mol, spec.residual);

            let mut coupling = coupling_template;
            coupling.intrinsic_efficacy = spec.intrinsic_efficacy_cb1;
            let points = simulate_thc_cb1_neuron_response(binding, coupling, &concs);

            let p10 = nearest_point(&points, 10.0).expect("10 nM point");
            let p30 = nearest_point(&points, 30.0).expect("30 nM point");
            let p100 = nearest_point(&points, 100.0).expect("100 nM point");

            CannabinoidPanelResult {
                name: spec.name,
                class: spec.class,
                ki_cb1_nanomolar: spec.ki_cb1_nanomolar,
                ki_cb2_nanomolar: spec.ki_cb2_nanomolar,
                intrinsic_efficacy_cb1: spec.intrinsic_efficacy_cb1,
                experimental_delta_g_kj_mol: score.experimental_delta_g_kj_mol,
                qed_floor_total_kj_mol: score.qed_floor_total_kj_mol,
                residual_required_kj_mol: score.residual_required_kj_mol,
                residual_modeled_total_kj_mol: residual.modeled_residual_total_kj_mol,
                residual_closure_error_kj_mol: residual.closure_error_kj_mol,
                explained_fraction_of_abs_delta_g: score.explained_fraction_of_abs_delta_g,
                occupancy_10nm: p10.occupancy_fraction,
                occupancy_30nm: p30.occupancy_fraction,
                occupancy_100nm: p100.occupancy_fraction,
                firing_scale_10nm: p10.relative_firing_scale,
                firing_scale_30nm: p30.relative_firing_scale,
                firing_scale_100nm: p100.relative_firing_scale,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn panel_runs_and_returns_expected_size() {
        let specs = default_cannabinoid_specs();
        let out = evaluate_cannabinoid_panel(&specs, 310.15, NeuronCouplingInput::default());
        assert_eq!(out.len(), specs.len());
    }

    #[test]
    fn potency_order_reflects_ki_for_thc_vs_cbd() {
        let specs = default_cannabinoid_specs();
        let out = evaluate_cannabinoid_panel(&specs, 310.15, NeuronCouplingInput::default());
        let thc = out.iter().find(|r| r.name == "delta9_thc").unwrap();
        let cbd = out.iter().find(|r| r.name == "cbd").unwrap();
        assert!(thc.occupancy_100nm > cbd.occupancy_100nm);
    }
}
