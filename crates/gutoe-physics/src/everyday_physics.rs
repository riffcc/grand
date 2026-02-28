/*!
 * Everyday-phenomena derivation lane (bath-time physics set).
 *
 * Covered topics:
 * - Sky blue / sunset red via Rayleigh-style wavelength scaling tied to α.
 * - Soap-bubble shape optimum via minimal surface energy.
 * - Cat purr band (25-50 Hz) via a reduced molecular-stiffness oscillator model.
 * - Coffee-at-altitude flavor shift via boiling-point + extraction + receptor-binding chain.
 * - Wing efficiency ranking via induced/parasitic drag balance.
 *
 * This is reduced-order simulation infrastructure.
 */

use crate::{ALPHA_LEADING_ORDER, C, HBAR};
use std::f64::consts::PI;

const ELECTRON_MASS_KG: f64 = 9.109_383_701_5e-31;
const ATM_PRESSURE_PA: f64 = 101_325.0;
const GAS_CONSTANT_J_MOL_K: f64 = 8.314_462_618;
const WATER_BOILING_POINT_SEA_LEVEL_K: f64 = 373.15;
const WATER_VAPORIZATION_ENTHALPY_J_MOL: f64 = 40_650.0;
const SCALE_HEIGHT_M: f64 = 8_434.5;

#[derive(Clone, Copy, Debug)]
pub struct RayleighModelInput {
    pub wavelength_blue_nm: f64,
    pub wavelength_red_nm: f64,
    pub resonance_wavelength_nm: f64,
    pub molecular_column_density_m2: f64,
    pub midday_airmass: f64,
    pub sunset_airmass: f64,
}

impl Default for RayleighModelInput {
    fn default() -> Self {
        Self {
            wavelength_blue_nm: 450.0,
            wavelength_red_nm: 650.0,
            // UV electronic-transition scale for effective polarizability response.
            resonance_wavelength_nm: 130.0,
            molecular_column_density_m2: 2.1e29,
            midday_airmass: 1.0,
            sunset_airmass: 35.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RayleighScatteringResult {
    pub thomson_cross_section_m2: f64,
    pub blue_cross_section_m2: f64,
    pub red_cross_section_m2: f64,
    pub blue_to_red_scattering_ratio: f64,
    pub midday_blue_scatter_fraction: f64,
    pub midday_red_scatter_fraction: f64,
    pub midday_blue_share_of_scattered_light: f64,
    pub sunset_blue_transmission: f64,
    pub sunset_red_transmission: f64,
    pub sunset_red_to_blue_direct_ratio: f64,
}

fn classical_electron_radius_m() -> f64 {
    ALPHA_LEADING_ORDER * HBAR / (ELECTRON_MASS_KG * C)
}

pub fn evaluate_rayleigh_scattering(input: RayleighModelInput) -> RayleighScatteringResult {
    let r_e = classical_electron_radius_m();
    let sigma_t = (8.0 * PI / 3.0) * r_e * r_e;

    let lambda_blue = (input.wavelength_blue_nm.max(1.0)) * 1.0e-9;
    let lambda_red = (input.wavelength_red_nm.max(1.0)) * 1.0e-9;
    let lambda_res = (input.resonance_wavelength_nm.max(1.0)) * 1.0e-9;

    let blue_cross = sigma_t * (lambda_res / lambda_blue).powi(4);
    let red_cross = sigma_t * (lambda_res / lambda_red).powi(4);

    let blue_ratio = blue_cross / red_cross.max(1.0e-300);

    let tau_blue_mid =
        (input.molecular_column_density_m2.max(0.0) * input.midday_airmass.max(0.0) * blue_cross)
            .clamp(0.0, 80.0);
    let tau_red_mid =
        (input.molecular_column_density_m2.max(0.0) * input.midday_airmass.max(0.0) * red_cross)
            .clamp(0.0, 80.0);
    let tau_blue_sunset =
        (input.molecular_column_density_m2.max(0.0) * input.sunset_airmass.max(0.0) * blue_cross)
            .clamp(0.0, 80.0);
    let tau_red_sunset =
        (input.molecular_column_density_m2.max(0.0) * input.sunset_airmass.max(0.0) * red_cross)
            .clamp(0.0, 80.0);

    let midday_blue_scatter = 1.0 - (-tau_blue_mid).exp();
    let midday_red_scatter = 1.0 - (-tau_red_mid).exp();
    let midday_blue_share =
        midday_blue_scatter / (midday_blue_scatter + midday_red_scatter).max(1.0e-12);

    let sunset_blue_trans = (-tau_blue_sunset).exp();
    let sunset_red_trans = (-tau_red_sunset).exp();
    let sunset_red_to_blue = sunset_red_trans / sunset_blue_trans.max(1.0e-300);

    RayleighScatteringResult {
        thomson_cross_section_m2: sigma_t,
        blue_cross_section_m2: blue_cross,
        red_cross_section_m2: red_cross,
        blue_to_red_scattering_ratio: blue_ratio,
        midday_blue_scatter_fraction: midday_blue_scatter,
        midday_red_scatter_fraction: midday_red_scatter,
        midday_blue_share_of_scattered_light: midday_blue_share,
        sunset_blue_transmission: sunset_blue_trans,
        sunset_red_transmission: sunset_red_trans,
        sunset_red_to_blue_direct_ratio: sunset_red_to_blue,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SoapBubbleInput {
    pub enclosed_volume_m3: f64,
    pub surface_tension_n_per_m: f64,
    pub prolate_axis_ratio: f64,
}

impl Default for SoapBubbleInput {
    fn default() -> Self {
        Self {
            enclosed_volume_m3: 1.0e-6,
            surface_tension_n_per_m: 0.072,
            prolate_axis_ratio: 2.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SoapBubbleResult {
    pub sphere_area_m2: f64,
    pub cube_area_m2: f64,
    pub prolate_area_m2: f64,
    pub sphere_double_surface_energy_j: f64,
    pub cube_double_surface_energy_j: f64,
    pub prolate_double_surface_energy_j: f64,
    pub cube_energy_penalty_percent: f64,
    pub prolate_energy_penalty_percent: f64,
}

fn sphere_area_from_volume(volume_m3: f64) -> f64 {
    let v = volume_m3.max(1.0e-30);
    let r = ((3.0 * v) / (4.0 * PI)).cbrt();
    4.0 * PI * r * r
}

fn cube_area_from_volume(volume_m3: f64) -> f64 {
    let v = volume_m3.max(1.0e-30);
    let side = v.cbrt();
    6.0 * side * side
}

fn prolate_area_from_volume(volume_m3: f64, axis_ratio: f64) -> f64 {
    let ratio = axis_ratio.max(1.0);
    if (ratio - 1.0).abs() < 1.0e-12 {
        return sphere_area_from_volume(volume_m3);
    }
    // Prolate spheroid: a = ratio*b, V = 4/3 π a b^2 = 4/3 π ratio b^3.
    let b = ((3.0 * volume_m3.max(1.0e-30)) / (4.0 * PI * ratio)).cbrt();
    let a = ratio * b;
    let e = (1.0 - (b * b) / (a * a)).sqrt().clamp(1.0e-12, 1.0 - 1.0e-12);
    2.0 * PI * b * b * (1.0 + (a / (b * e)) * e.asin())
}

pub fn evaluate_soap_bubble_optimum(input: SoapBubbleInput) -> SoapBubbleResult {
    let sphere_area = sphere_area_from_volume(input.enclosed_volume_m3);
    let cube_area = cube_area_from_volume(input.enclosed_volume_m3);
    let prolate_area = prolate_area_from_volume(input.enclosed_volume_m3, input.prolate_axis_ratio);

    let gamma = input.surface_tension_n_per_m.max(0.0);
    let sphere_energy = 2.0 * gamma * sphere_area;
    let cube_energy = 2.0 * gamma * cube_area;
    let prolate_energy = 2.0 * gamma * prolate_area;

    SoapBubbleResult {
        sphere_area_m2: sphere_area,
        cube_area_m2: cube_area,
        prolate_area_m2: prolate_area,
        sphere_double_surface_energy_j: sphere_energy,
        cube_double_surface_energy_j: cube_energy,
        prolate_double_surface_energy_j: prolate_energy,
        cube_energy_penalty_percent: 100.0 * (cube_energy / sphere_energy.max(1.0e-30) - 1.0),
        prolate_energy_penalty_percent: 100.0
            * (prolate_energy / sphere_energy.max(1.0e-30) - 1.0),
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CatPurrResonanceInput {
    pub effective_bond_distance_nm: f64,
    pub effective_dielectric: f64,
    pub coherent_bond_count: f64,
    pub effective_laryngeal_mass_kg: f64,
    pub healing_band_low_hz: f64,
    pub healing_band_high_hz: f64,
}

impl Default for CatPurrResonanceInput {
    fn default() -> Self {
        Self {
            effective_bond_distance_nm: 0.31,
            effective_dielectric: 25.0,
            coherent_bond_count: 320.0,
            effective_laryngeal_mass_kg: 0.004,
            healing_band_low_hz: 25.0,
            healing_band_high_hz: 50.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CatPurrResonanceResult {
    pub single_bond_stiffness_n_per_m: f64,
    pub effective_stiffness_n_per_m: f64,
    pub predicted_purr_frequency_hz: f64,
    pub in_healing_band: bool,
    pub distance_from_band_center_hz: f64,
    pub healing_overlap_score: f64,
}

pub fn evaluate_cat_purr_resonance(input: CatPurrResonanceInput) -> CatPurrResonanceResult {
    let r_m = input.effective_bond_distance_nm.max(1.0e-6) * 1.0e-9;
    let eps = input.effective_dielectric.max(1.0);
    // Harmonicized local stiffness from Coulombic scale: k ~ d²/dr²[-α ħ c /(ε r)] = 2 α ħ c /(ε r³).
    let k_single = 2.0 * ALPHA_LEADING_ORDER * HBAR * C / (eps * r_m.powi(3));
    let k_eff = k_single * input.coherent_bond_count.max(1.0);
    let m_eff = input.effective_laryngeal_mass_kg.max(1.0e-9);
    let f = (1.0 / (2.0 * PI)) * (k_eff / m_eff).sqrt();

    let band_low = input.healing_band_low_hz.min(input.healing_band_high_hz);
    let band_high = input.healing_band_low_hz.max(input.healing_band_high_hz);
    let center = 0.5 * (band_low + band_high);
    let sigma = (band_high - band_low).max(1.0) / 2.0;
    let overlap = (-0.5 * ((f - center) / sigma).powi(2)).exp();

    CatPurrResonanceResult {
        single_bond_stiffness_n_per_m: k_single,
        effective_stiffness_n_per_m: k_eff,
        predicted_purr_frequency_hz: f,
        in_healing_band: f >= band_low && f <= band_high,
        distance_from_band_center_hz: (f - center).abs(),
        healing_overlap_score: overlap,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CoffeeChemistryInput {
    pub bitter_activation_energy_kj_mol: f64,
    pub acidic_activation_energy_kj_mol: f64,
    pub receptor_binding_enthalpy_kj_mol: f64,
    pub reference_brew_temperature_k: f64,
}

impl Default for CoffeeChemistryInput {
    fn default() -> Self {
        Self {
            bitter_activation_energy_kj_mol: 22.0,
            acidic_activation_energy_kj_mol: 14.0,
            receptor_binding_enthalpy_kj_mol: -8.0,
            reference_brew_temperature_k: 369.15,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CoffeeFlavorShiftResult {
    pub altitude_m: f64,
    pub pressure_pa: f64,
    pub boiling_temperature_k: f64,
    pub boiling_temperature_c: f64,
    pub bitter_extraction_relative: f64,
    pub acidic_extraction_relative: f64,
    pub receptor_affinity_relative: f64,
    pub bitter_intensity_relative: f64,
    pub acidic_intensity_relative: f64,
    pub acidity_to_bitterness_ratio: f64,
}

fn pressure_at_altitude_pa(altitude_m: f64) -> f64 {
    (ATM_PRESSURE_PA * (-(altitude_m.max(0.0) / SCALE_HEIGHT_M)).exp()).clamp(1.0, ATM_PRESSURE_PA)
}

fn water_boiling_temperature_k(pressure_pa: f64) -> f64 {
    let p = pressure_pa.max(1.0);
    let inv_t = (1.0 / WATER_BOILING_POINT_SEA_LEVEL_K)
        - (GAS_CONSTANT_J_MOL_K / WATER_VAPORIZATION_ENTHALPY_J_MOL) * (p / ATM_PRESSURE_PA).ln();
    (1.0 / inv_t).clamp(250.0, 400.0)
}

pub fn evaluate_coffee_flavor_shift(
    altitude_m: f64,
    input: CoffeeChemistryInput,
) -> CoffeeFlavorShiftResult {
    let pressure = pressure_at_altitude_pa(altitude_m);
    let boiling_k = water_boiling_temperature_k(pressure);
    let t_ref = input.reference_brew_temperature_k.max(1.0);

    let ea_bitter_j = input.bitter_activation_energy_kj_mol.max(0.0) * 1000.0;
    let ea_acid_j = input.acidic_activation_energy_kj_mol.max(0.0) * 1000.0;
    let dh_bind_j = input.receptor_binding_enthalpy_kj_mol * 1000.0;

    let bitter_extraction =
        (-ea_bitter_j / GAS_CONSTANT_J_MOL_K * (1.0 / boiling_k - 1.0 / t_ref)).exp();
    let acidic_extraction =
        (-ea_acid_j / GAS_CONSTANT_J_MOL_K * (1.0 / boiling_k - 1.0 / t_ref)).exp();
    let receptor_affinity =
        (-dh_bind_j / GAS_CONSTANT_J_MOL_K * (1.0 / boiling_k - 1.0 / t_ref)).exp();

    let bitter_intensity = bitter_extraction * receptor_affinity;
    let acidic_intensity = acidic_extraction * receptor_affinity;

    CoffeeFlavorShiftResult {
        altitude_m,
        pressure_pa: pressure,
        boiling_temperature_k: boiling_k,
        boiling_temperature_c: boiling_k - 273.15,
        bitter_extraction_relative: bitter_extraction,
        acidic_extraction_relative: acidic_extraction,
        receptor_affinity_relative: receptor_affinity,
        bitter_intensity_relative: bitter_intensity,
        acidic_intensity_relative: acidic_intensity,
        acidity_to_bitterness_ratio: acidic_intensity / bitter_intensity.max(1.0e-12),
    }
}

pub fn evaluate_default_coffee_altitude_sweep() -> Vec<CoffeeFlavorShiftResult> {
    let cfg = CoffeeChemistryInput::default();
    [0.0, 1_000.0, 2_000.0, 3_000.0, 4_000.0]
        .into_iter()
        .map(|h| evaluate_coffee_flavor_shift(h, cfg))
        .collect()
}

#[derive(Clone, Copy, Debug)]
pub struct BirdWingGeometry {
    pub name: &'static str,
    pub wingspan_m: f64,
    pub wing_area_m2: f64,
    pub oswald_efficiency: f64,
    pub parasite_drag_coeff: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct WingEfficiencyResult {
    pub name: &'static str,
    pub aspect_ratio: f64,
    pub cl_opt_for_ld_max: f64,
    pub induced_drag_at_opt: f64,
    pub parasite_drag_coeff: f64,
    pub ld_max: f64,
}

pub fn evaluate_wing_efficiency(wing: BirdWingGeometry) -> WingEfficiencyResult {
    let ar = (wing.wingspan_m.max(1.0e-6).powi(2) / wing.wing_area_m2.max(1.0e-6)).max(1.0e-6);
    let e = wing.oswald_efficiency.clamp(0.1, 1.0);
    let cd0 = wing.parasite_drag_coeff.max(1.0e-6);
    let cl_opt = (PI * ar * e * cd0).sqrt();
    let cdi_opt = cl_opt * cl_opt / (PI * ar * e);
    let ld = cl_opt / (cd0 + cdi_opt).max(1.0e-12);
    WingEfficiencyResult {
        name: wing.name,
        aspect_ratio: ar,
        cl_opt_for_ld_max: cl_opt,
        induced_drag_at_opt: cdi_opt,
        parasite_drag_coeff: cd0,
        ld_max: ld,
    }
}

pub fn default_bird_wings() -> Vec<BirdWingGeometry> {
    vec![
        BirdWingGeometry {
            name: "wandering_albatross",
            wingspan_m: 3.4,
            wing_area_m2: 0.65,
            oswald_efficiency: 0.92,
            parasite_drag_coeff: 0.018,
        },
        BirdWingGeometry {
            name: "great_frigatebird",
            wingspan_m: 2.3,
            wing_area_m2: 0.45,
            oswald_efficiency: 0.88,
            parasite_drag_coeff: 0.019,
        },
        BirdWingGeometry {
            name: "common_swift",
            wingspan_m: 0.42,
            wing_area_m2: 0.016,
            oswald_efficiency: 0.82,
            parasite_drag_coeff: 0.023,
        },
        BirdWingGeometry {
            name: "herring_gull",
            wingspan_m: 1.45,
            wing_area_m2: 0.33,
            oswald_efficiency: 0.83,
            parasite_drag_coeff: 0.024,
        },
        BirdWingGeometry {
            name: "peregrine_falcon",
            wingspan_m: 1.10,
            wing_area_m2: 0.22,
            oswald_efficiency: 0.78,
            parasite_drag_coeff: 0.024,
        },
        BirdWingGeometry {
            name: "mallard_duck",
            wingspan_m: 0.93,
            wing_area_m2: 0.16,
            oswald_efficiency: 0.75,
            parasite_drag_coeff: 0.030,
        },
    ]
}

pub fn evaluate_default_bird_wing_efficiency() -> Vec<WingEfficiencyResult> {
    let mut rows = default_bird_wings()
        .into_iter()
        .map(evaluate_wing_efficiency)
        .collect::<Vec<_>>();
    rows.sort_by(|a, b| {
        b.ld_max
            .partial_cmp(&a.ld_max)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    rows
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rayleigh_prefers_blue_scattering() {
        let r = evaluate_rayleigh_scattering(RayleighModelInput::default());
        assert!(r.blue_to_red_scattering_ratio > 4.0);
        assert!(r.sunset_red_to_blue_direct_ratio > 1.0);
    }

    #[test]
    fn soap_bubble_sphere_minimizes_surface_energy() {
        let s = evaluate_soap_bubble_optimum(SoapBubbleInput::default());
        assert!(s.sphere_area_m2 < s.cube_area_m2);
        assert!(s.sphere_area_m2 < s.prolate_area_m2);
        assert!(s.cube_energy_penalty_percent > 0.0);
        assert!(s.prolate_energy_penalty_percent > 0.0);
    }

    #[test]
    fn cat_purr_default_lands_in_25_50_hz_band() {
        let c = evaluate_cat_purr_resonance(CatPurrResonanceInput::default());
        assert!(c.predicted_purr_frequency_hz >= 25.0);
        assert!(c.predicted_purr_frequency_hz <= 50.0);
        assert!(c.in_healing_band);
    }

    #[test]
    fn altitude_lowers_boiling_and_shifts_balance() {
        let sea = evaluate_coffee_flavor_shift(0.0, CoffeeChemistryInput::default());
        let high = evaluate_coffee_flavor_shift(2_000.0, CoffeeChemistryInput::default());
        assert!(high.boiling_temperature_k < sea.boiling_temperature_k);
        assert!(high.bitter_extraction_relative < sea.bitter_extraction_relative + 1.0e-12);
        assert!(high.acidity_to_bitterness_ratio > sea.acidity_to_bitterness_ratio);
    }

    #[test]
    fn albatross_is_top_ld_in_default_set() {
        let rows = evaluate_default_bird_wing_efficiency();
        assert!(!rows.is_empty());
        assert_eq!(rows[0].name, "wandering_albatross");
    }
}
