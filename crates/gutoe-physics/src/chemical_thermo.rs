/*!
 * Chemical thermodynamics proxy lane for periodic-table scale modeling.
 *
 * This lane uses physically motivated, explicit assumptions to estimate
 * bulk thermodynamic properties per element:
 * - cohesive-energy proxy
 * - Debye-temperature proxy
 * - latent heats
 * - melting/boiling temperatures
 * - phase at ambient conditions
 *
 * It is not a quantum-chemistry solver; the formulas are closed-form
 * transduction rules intended for broad trend modeling and extrapolation.
 */

use std::f64::consts::PI;

pub const AVOGADRO: f64 = 6.022_140_76e23;
pub const R_GAS_J_MOL_K: f64 = 8.314_462_618;
pub const BOHR_RADIUS_PM: f64 = 52.917_721;
pub const EV_TO_KJ_MOL: f64 = 96.485_332_123;
pub const ENTROPY_FUSION_J_MOL_K: f64 = 10.0;
pub const ENTROPY_VAPORIZATION_J_MOL_K: f64 = 85.0;
pub const P_REF_PA: f64 = 101_325.0;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ChemicalFamily {
    Alkali,
    AlkalineEarth,
    Transition,
    PostTransition,
    Metalloid,
    Nonmetal,
    Halogen,
    NobleGas,
    Lanthanide,
    Actinide,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MatterState {
    Solid,
    Liquid,
    Gas,
}

#[derive(Clone, Copy, Debug)]
pub struct PhaseGibbsKjMol {
    pub solid: f64,
    pub liquid: f64,
    pub gas: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct ElementThermoPrediction {
    pub z: u16,
    pub a: u16,
    pub family: ChemicalFamily,
    pub period: u8,
    pub molar_mass_g_mol: f64,
    pub atomic_radius_pm: f64,
    pub molar_volume_cm3_mol: f64,
    pub density_g_cm3: f64,
    pub cohesive_energy_ev_per_atom: f64,
    pub debye_temperature_k: f64,
    pub latent_fusion_kj_mol: f64,
    pub latent_vaporization_kj_mol: f64,
    pub melting_temperature_k: f64,
    pub boiling_temperature_k: f64,
    pub vapor_pressure_pa_298k: f64,
    pub cp_solid_j_mol_k: f64,
    pub cp_liquid_j_mol_k: f64,
    pub cp_gas_j_mol_k: f64,
    pub bulk_modulus_gpa: f64,
    pub thermal_expansion_1_per_k: f64,
    pub ambient_state_298k: MatterState,
}

pub fn family_of_z(z: u16) -> ChemicalFamily {
    match z {
        1 | 3 | 11 | 19 | 37 | 55 | 87 => ChemicalFamily::Alkali,
        4 | 12 | 20 | 38 | 56 | 88 => ChemicalFamily::AlkalineEarth,
        2 | 10 | 18 | 36 | 54 | 86 | 118 => ChemicalFamily::NobleGas,
        9 | 17 | 35 | 53 | 85 | 117 => ChemicalFamily::Halogen,
        57..=71 => ChemicalFamily::Lanthanide,
        89..=103 => ChemicalFamily::Actinide,
        5 | 14 | 32 | 33 | 51 | 52 | 84 => ChemicalFamily::Metalloid,
        6 | 7 | 8 | 15 | 16 | 34 => ChemicalFamily::Nonmetal,
        21..=30 | 39..=48 | 72..=80 | 104..=112 => ChemicalFamily::Transition,
        _ => ChemicalFamily::PostTransition,
    }
}

pub fn period_of_z(z: u16) -> u8 {
    match z {
        0 => 1,
        1..=2 => 1,
        3..=10 => 2,
        11..=18 => 3,
        19..=36 => 4,
        37..=54 => 5,
        55..=86 => 6,
        _ => 7,
    }
}

fn noble_core_electrons(z: u16) -> f64 {
    if z > 86 {
        86.0
    } else if z > 54 {
        54.0
    } else if z > 36 {
        36.0
    } else if z > 18 {
        18.0
    } else if z > 10 {
        10.0
    } else if z > 2 {
        2.0
    } else {
        0.0
    }
}

fn valence_proxy(family: ChemicalFamily) -> f64 {
    match family {
        ChemicalFamily::Alkali => 1.0,
        ChemicalFamily::AlkalineEarth => 2.0,
        ChemicalFamily::Transition => 2.5,
        ChemicalFamily::PostTransition => 3.0,
        ChemicalFamily::Metalloid => 4.0,
        ChemicalFamily::Nonmetal => 2.0,
        ChemicalFamily::Halogen => 1.0,
        ChemicalFamily::NobleGas => 0.5,
        ChemicalFamily::Lanthanide => 3.0,
        ChemicalFamily::Actinide => 3.0,
    }
}

fn radius_family_factor(family: ChemicalFamily) -> f64 {
    match family {
        ChemicalFamily::Alkali => 1.85,
        ChemicalFamily::AlkalineEarth => 1.55,
        ChemicalFamily::Transition => 1.20,
        ChemicalFamily::PostTransition => 1.30,
        ChemicalFamily::Metalloid => 1.10,
        ChemicalFamily::Nonmetal => 0.92,
        ChemicalFamily::Halogen => 0.86,
        ChemicalFamily::NobleGas => 1.02,
        ChemicalFamily::Lanthanide => 1.42,
        ChemicalFamily::Actinide => 1.46,
    }
}

fn cohesive_multiplier(family: ChemicalFamily) -> f64 {
    match family {
        ChemicalFamily::Alkali => 0.35,
        ChemicalFamily::AlkalineEarth => 0.45,
        ChemicalFamily::Transition => 0.80,
        ChemicalFamily::PostTransition => 0.55,
        ChemicalFamily::Metalloid => 0.65,
        ChemicalFamily::Nonmetal => 0.50,
        ChemicalFamily::Halogen => 0.30,
        ChemicalFamily::NobleGas => 0.08,
        ChemicalFamily::Lanthanide => 0.75,
        ChemicalFamily::Actinide => 0.75,
    }
}

fn latent_fusion_fraction(family: ChemicalFamily) -> f64 {
    match family {
        ChemicalFamily::Alkali => 0.020,
        ChemicalFamily::AlkalineEarth => 0.025,
        ChemicalFamily::Transition => 0.035,
        ChemicalFamily::PostTransition => 0.030,
        ChemicalFamily::Metalloid => 0.030,
        ChemicalFamily::Nonmetal => 0.028,
        ChemicalFamily::Halogen => 0.015,
        ChemicalFamily::NobleGas => 0.005,
        ChemicalFamily::Lanthanide => 0.032,
        ChemicalFamily::Actinide => 0.035,
    }
}

fn latent_vapor_fraction(family: ChemicalFamily) -> f64 {
    match family {
        ChemicalFamily::Alkali => 0.55,
        ChemicalFamily::AlkalineEarth => 0.65,
        ChemicalFamily::Transition => 0.82,
        ChemicalFamily::PostTransition => 0.70,
        ChemicalFamily::Metalloid => 0.60,
        ChemicalFamily::Nonmetal => 0.45,
        ChemicalFamily::Halogen => 0.35,
        ChemicalFamily::NobleGas => 0.12,
        ChemicalFamily::Lanthanide => 0.78,
        ChemicalFamily::Actinide => 0.78,
    }
}

/// Clausius-Clapeyron transduction using the element's boiling point as anchor
/// (P_sat(T_b) = 1 atm).
pub fn vapor_pressure_clausius_pa(
    latent_vaporization_kj_mol: f64,
    boiling_temperature_k: f64,
    t_k: f64,
) -> f64 {
    let t = t_k.max(1.0);
    let t_b = boiling_temperature_k.max(2.0);
    let delta_h = (latent_vaporization_kj_mol * 1000.0).max(1.0);
    let exponent = (-delta_h / R_GAS_J_MOL_K * (1.0 / t - 1.0 / t_b)).clamp(-80.0, 80.0);
    (P_REF_PA * exponent.exp()).clamp(1.0e-9, 1.0e12)
}

/// Reduced Gibbs offsets for solid/liquid/gas at (T, P). Solid is used as the
/// zero reference and liquid/gas are computed from latent-heat transductions.
pub fn phase_gibbs_offsets_kj_mol(
    latent_fusion_kj_mol: f64,
    latent_vaporization_kj_mol: f64,
    melting_temperature_k: f64,
    boiling_temperature_k: f64,
    t_k: f64,
    pressure_pa: f64,
) -> PhaseGibbsKjMol {
    let t = t_k.max(1.0);
    let t_m = melting_temperature_k.max(1.0);
    let t_b = boiling_temperature_k.max(t_m + 1.0e-6);
    let delta_h_f = latent_fusion_kj_mol.max(1.0e-6);
    let delta_h_v = latent_vaporization_kj_mol.max(1.0e-6);
    let p = pressure_pa.max(1.0e-12);

    let g_solid = 0.0;
    let g_liquid = delta_h_f * (1.0 - t / t_m);
    let g_gas_ref = g_liquid + delta_h_v * (1.0 - t / t_b);
    let pressure_term = (R_GAS_J_MOL_K * t / 1000.0) * (p / P_REF_PA).ln();
    let g_gas = g_gas_ref + pressure_term;

    PhaseGibbsKjMol {
        solid: g_solid,
        liquid: g_liquid,
        gas: g_gas,
    }
}

/// Select phase by minimum reduced Gibbs free energy at (T, P).
pub fn phase_from_gibbs(
    latent_fusion_kj_mol: f64,
    latent_vaporization_kj_mol: f64,
    melting_temperature_k: f64,
    boiling_temperature_k: f64,
    t_k: f64,
    pressure_pa: f64,
) -> MatterState {
    let g = phase_gibbs_offsets_kj_mol(
        latent_fusion_kj_mol,
        latent_vaporization_kj_mol,
        melting_temperature_k,
        boiling_temperature_k,
        t_k,
        pressure_pa,
    );
    let mut state = MatterState::Solid;
    let mut best = g.solid;
    if g.liquid < best {
        best = g.liquid;
        state = MatterState::Liquid;
    }
    if g.gas < best {
        state = MatterState::Gas;
    }
    state
}

pub fn predict_element_thermo(z: u16, a: u16) -> ElementThermoPrediction {
    let family = family_of_z(z);
    let period = period_of_z(z);
    let period_f = period as f64;
    let valence = valence_proxy(family);

    let z_eff = ((z as f64) - noble_core_electrons(z) + 0.5 * valence).max(1.0);
    let radius_pm = (BOHR_RADIUS_PM * (period_f * period_f) / z_eff * radius_family_factor(family))
        .clamp(30.0, 320.0);

    // Spherical atom proxy with close-packing correction.
    let r_cm = radius_pm * 1.0e-10;
    let atom_vol_cm3 = (4.0 / 3.0) * PI * r_cm.powi(3) / 0.74;
    let molar_volume_cm3_mol = atom_vol_cm3 * AVOGADRO;

    let molar_mass_g_mol = a as f64;
    let density_g_cm3 = (molar_mass_g_mol / molar_volume_cm3_mol).clamp(0.0005, 40.0);

    // Hydrogenic binding-energy style cohesive proxy.
    let cohesive_energy_ev_per_atom = (13.605_693 * valence.powi(2) / period_f.powi(2)
        * cohesive_multiplier(family))
    .clamp(0.03, 12.0);

    let latent_fusion_kj_mol =
        cohesive_energy_ev_per_atom * EV_TO_KJ_MOL * latent_fusion_fraction(family);
    let latent_vaporization_kj_mol =
        cohesive_energy_ev_per_atom * EV_TO_KJ_MOL * latent_vapor_fraction(family);

    let melting_temperature_k =
        (latent_fusion_kj_mol * 1000.0 / ENTROPY_FUSION_J_MOL_K).clamp(2.0, 8000.0);
    let boiling_temperature_k =
        (latent_vaporization_kj_mol * 1000.0 / ENTROPY_VAPORIZATION_J_MOL_K).clamp(4.0, 12000.0);

    let debye_temperature_k =
        (120.0 * cohesive_energy_ev_per_atom.sqrt() * (density_g_cm3 / 5.0).powf(0.25))
            .clamp(20.0, 2200.0);

    // Debye-like saturation toward Dulong-Petit near room temperature.
    let cp_solid_j_mol_k =
        (3.0 * R_GAS_J_MOL_K * (1.0 - (-298.15 / (0.35 * debye_temperature_k)).exp()))
            .clamp(1.5 * R_GAS_J_MOL_K, 3.0 * R_GAS_J_MOL_K);
    let cp_liquid_j_mol_k = 3.5 * R_GAS_J_MOL_K;
    let cp_gas_j_mol_k = match family {
        ChemicalFamily::Nonmetal | ChemicalFamily::Halogen => 3.5 * R_GAS_J_MOL_K,
        _ => 2.5 * R_GAS_J_MOL_K,
    };

    let bulk_modulus_gpa =
        (20.0 * cohesive_energy_ev_per_atom * density_g_cm3.powf(0.7)).clamp(0.1, 500.0);
    let thermal_expansion_1_per_k =
        (2.2e-5 * (300.0 / debye_temperature_k).powf(0.7) * (30.0 / bulk_modulus_gpa).powf(0.3))
            .clamp(1.0e-6, 2.5e-4);

    let vapor_pressure_pa_298k =
        vapor_pressure_clausius_pa(latent_vaporization_kj_mol, boiling_temperature_k, 298.15);
    let ambient_state_298k = phase_from_gibbs(
        latent_fusion_kj_mol,
        latent_vaporization_kj_mol,
        melting_temperature_k,
        boiling_temperature_k,
        298.15,
        P_REF_PA,
    );

    ElementThermoPrediction {
        z,
        a,
        family,
        period,
        molar_mass_g_mol,
        atomic_radius_pm: radius_pm,
        molar_volume_cm3_mol,
        density_g_cm3,
        cohesive_energy_ev_per_atom,
        debye_temperature_k,
        latent_fusion_kj_mol,
        latent_vaporization_kj_mol,
        melting_temperature_k,
        boiling_temperature_k,
        vapor_pressure_pa_298k,
        cp_solid_j_mol_k,
        cp_liquid_j_mol_k,
        cp_gas_j_mol_k,
        bulk_modulus_gpa,
        thermal_expansion_1_per_k,
        ambient_state_298k,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boiling_is_above_melting_for_all_families() {
        for z in 1..=118 {
            let p = predict_element_thermo(z, (2.5 * z as f64).round() as u16);
            assert!(p.boiling_temperature_k > p.melting_temperature_k);
        }
    }

    #[test]
    fn noble_gases_are_low_temperature_volatiles() {
        for &z in &[2_u16, 10, 18, 36, 54, 86, 118] {
            let p = predict_element_thermo(z, (2.5 * z as f64).round() as u16);
            assert!(
                p.boiling_temperature_k < 400.0,
                "noble gas Z={z} expected low boiling point, got {} K",
                p.boiling_temperature_k
            );
        }
    }

    #[test]
    fn transition_family_is_more_cohesive_than_alkali_proxy() {
        // Fe-like vs Na-like proxy comparison.
        let fe_like = predict_element_thermo(26, 56);
        let na_like = predict_element_thermo(11, 23);
        assert!(fe_like.cohesive_energy_ev_per_atom > na_like.cohesive_energy_ev_per_atom);
    }

    #[test]
    fn clapeyron_anchor_hits_one_atm_at_boiling_point() {
        for z in 1..=118 {
            let p = predict_element_thermo(z, (2.5 * z as f64).round() as u16);
            let p_sat = vapor_pressure_clausius_pa(
                p.latent_vaporization_kj_mol,
                p.boiling_temperature_k,
                p.boiling_temperature_k,
            );
            assert!(((p_sat / P_REF_PA) - 1.0).abs() < 1.0e-10);
        }
    }

    #[test]
    fn gibbs_phase_matches_threshold_rule_at_reference_pressure() {
        for z in 1..=118 {
            let p = predict_element_thermo(z, (2.5 * z as f64).round() as u16);
            let threshold_state = if 298.15 < p.melting_temperature_k {
                MatterState::Solid
            } else if 298.15 < p.boiling_temperature_k {
                MatterState::Liquid
            } else {
                MatterState::Gas
            };
            let gibbs_state = phase_from_gibbs(
                p.latent_fusion_kj_mol,
                p.latent_vaporization_kj_mol,
                p.melting_temperature_k,
                p.boiling_temperature_k,
                298.15,
                P_REF_PA,
            );
            assert_eq!(gibbs_state, threshold_state);
        }
    }

    #[test]
    fn high_pressure_can_condense_noble_gas_proxy() {
        let he_like = predict_element_thermo(2, 4);
        let ambient = phase_from_gibbs(
            he_like.latent_fusion_kj_mol,
            he_like.latent_vaporization_kj_mol,
            he_like.melting_temperature_k,
            he_like.boiling_temperature_k,
            298.15,
            P_REF_PA,
        );
        assert_eq!(ambient, MatterState::Gas);

        let extreme_pressure = phase_from_gibbs(
            he_like.latent_fusion_kj_mol,
            he_like.latent_vaporization_kj_mol,
            he_like.melting_temperature_k,
            he_like.boiling_temperature_k,
            298.15,
            1.0e9,
        );
        assert!(extreme_pressure != MatterState::Gas);
    }
}
