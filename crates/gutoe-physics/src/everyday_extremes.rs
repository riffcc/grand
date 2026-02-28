/*!
 * Everyday extremes lane:
 * - Ice slipperiness
 * - Popcorn popping threshold
 * - Raindrop shape optimum
 * - Mpemba-effect regime map
 *
 * Reduced-order physics transduction intended for trend-level behavior.
 */

use std::f64::consts::PI;

const R_GAS_J_MOL_K: f64 = 8.314_462_618;
const WATER_CP_J_KG_K: f64 = 4180.0;
const WATER_LATENT_J_KG: f64 = 333_550.0;
const WATER_DENSITY_KG_M3: f64 = 997.0;
const WATER_SURFACE_TENSION_N_M: f64 = 0.072;
const AIR_DENSITY_KG_M3: f64 = 1.225;
const GRAVITY_M_S2: f64 = 9.80665;
const ATM_PRESSURE_PA: f64 = 101_325.0;

#[derive(Clone, Copy, Debug)]
pub struct IceSlipperinessInput {
    pub ice_temperature_c: f64,
    pub contact_pressure_mpa: f64,
    pub sliding_speed_m_s: f64,
}

impl Default for IceSlipperinessInput {
    fn default() -> Self {
        Self {
            ice_temperature_c: -5.0,
            contact_pressure_mpa: 2.0,
            sliding_speed_m_s: 4.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct IceSlipperinessResult {
    pub effective_surface_temperature_c: f64,
    pub quasi_liquid_layer_thickness_nm: f64,
    pub friction_coefficient: f64,
    pub pressure_melting_shift_c: f64,
    pub frictional_heating_shift_c: f64,
}

pub fn evaluate_ice_slipperiness(input: IceSlipperinessInput) -> IceSlipperinessResult {
    let pressure_shift_c = 0.0074 * input.contact_pressure_mpa.max(0.0);
    let friction_shift_c = 0.9 * (input.sliding_speed_m_s.max(0.0) / (input.sliding_speed_m_s + 5.0));
    let eff_temp_c = input.ice_temperature_c + pressure_shift_c + friction_shift_c;

    // Quasi-liquid layer proxy: thin at deep subzero, thicker near melting.
    let qll_nm = 0.25 + 8.0 * (eff_temp_c / 6.0).exp();
    let qll_nm = qll_nm.clamp(0.2, 20.0);

    // Friction model: lubrication lowers μ; very thick films add viscous shear.
    let mu_dry = 0.35;
    let lubrication_term = mu_dry / (1.0 + 0.7 * qll_nm);
    let viscous_term = 0.002 * qll_nm * input.sliding_speed_m_s.max(0.0);
    let friction_coefficient = (lubrication_term + viscous_term).clamp(0.01, 0.6);

    IceSlipperinessResult {
        effective_surface_temperature_c: eff_temp_c,
        quasi_liquid_layer_thickness_nm: qll_nm,
        friction_coefficient,
        pressure_melting_shift_c: pressure_shift_c,
        frictional_heating_shift_c: friction_shift_c,
    }
}

pub fn default_ice_temperature_sweep() -> Vec<(f64, IceSlipperinessResult)> {
    [-30.0, -20.0, -10.0, -5.0, -2.0, -1.0]
        .into_iter()
        .map(|t| {
            let out = evaluate_ice_slipperiness(IceSlipperinessInput {
                ice_temperature_c: t,
                ..IceSlipperinessInput::default()
            });
            (t, out)
        })
        .collect()
}

#[derive(Clone, Copy, Debug)]
pub struct PopcornInput {
    pub kernel_radius_mm: f64,
    pub hull_thickness_mm: f64,
    pub hull_strength_mpa: f64,
    pub kernel_density_kg_m3: f64,
    pub moisture_mass_fraction: f64,
    pub vaporized_moisture_fraction: f64,
    pub void_fraction: f64,
    pub initial_kernel_temperature_c: f64,
    pub heating_rate_c_per_s: f64,
    pub thermal_lag_tau_s: f64,
    pub vapor_kinetics_preexp_s_inv: f64,
    pub vapor_activation_energy_kj_mol: f64,
    pub hull_softening_onset_c: f64,
    pub hull_softening_per_c: f64,
    pub hull_min_strength_fraction: f64,
    pub damage_timescale_s: f64,
    pub damage_exponent: f64,
    pub max_temperature_c: f64,
    pub time_step_s: f64,
}

impl Default for PopcornInput {
    fn default() -> Self {
        Self {
            kernel_radius_mm: 3.0,
            hull_thickness_mm: 0.20,
            hull_strength_mpa: 8.0,
            kernel_density_kg_m3: 1_250.0,
            moisture_mass_fraction: 0.14,
            vaporized_moisture_fraction: 0.40,
            void_fraction: 0.15,
            initial_kernel_temperature_c: 25.0,
            heating_rate_c_per_s: 1.3,
            thermal_lag_tau_s: 22.0,
            vapor_kinetics_preexp_s_inv: 900_000.0,
            vapor_activation_energy_kj_mol: 53.0,
            hull_softening_onset_c: 125.0,
            hull_softening_per_c: 0.0018,
            hull_min_strength_fraction: 0.72,
            damage_timescale_s: 14.0,
            damage_exponent: 1.9,
            max_temperature_c: 240.0,
            time_step_s: 0.2,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct PopcornResult {
    pub ready_temperature_c: f64,
    pub burst_temperature_c: f64,
    pub ready_time_s: f64,
    pub burst_time_s: f64,
    pub internal_pressure_ready_mpa: f64,
    pub internal_pressure_burst_mpa: f64,
    pub rupture_threshold_ready_mpa: f64,
    pub rupture_threshold_burst_mpa: f64,
    pub pressure_margin_ready_mpa: f64,
    pub pressure_margin_burst_mpa: f64,
    pub hysteresis_delta_c: f64,
    pub estimated_expansion_ratio: f64,
    pub pops: bool,
}

fn saturation_pressure_water_pa(temp_k: f64) -> f64 {
    // Clausius-Clapeyron anchored at 100 C and 1 atm.
    let t = temp_k.max(250.0);
    let d_h_vap = 40_650.0;
    let exponent =
        (-d_h_vap / R_GAS_J_MOL_K * (1.0 / t - 1.0 / 373.15)).clamp(-60.0, 60.0);
    ATM_PRESSURE_PA * exponent.exp()
}

fn popcorn_internal_pressure_mpa(input: PopcornInput, temp_c: f64, vapor_fraction: f64) -> f64 {
    let r_m = input.kernel_radius_mm.max(0.5) * 1.0e-3;
    let v_kernel = (4.0 / 3.0) * PI * r_m.powi(3);
    let m_kernel = input.kernel_density_kg_m3.max(100.0) * v_kernel;
    let m_water = m_kernel * input.moisture_mass_fraction.clamp(0.01, 0.4);
    let n_water_total = m_water / 0.018_015_28;
    let n_vap = n_water_total * vapor_fraction.clamp(0.001, input.vaporized_moisture_fraction.clamp(0.05, 0.95));
    let v_void = (v_kernel * input.void_fraction.clamp(0.05, 0.4)).max(1.0e-12);
    let t_k = (temp_c + 273.15).max(273.15);

    let p_ideal = n_vap * R_GAS_J_MOL_K * t_k / v_void;
    let p_sat = saturation_pressure_water_pa(t_k);
    // Effective internal pressure cannot exceed ideal vapor-limited estimate,
    // and tends to be bounded by saturation + trapped gas compression.
    let p_eff = p_ideal.min(2.0 * p_sat);
    p_eff * 1.0e-6
}

fn popcorn_rupture_threshold_mpa(input: PopcornInput) -> f64 {
    // Thin-shell hoop stress threshold: sigma = P r / (2 t) => P = 2 t sigma / r.
    let r_m = input.kernel_radius_mm.max(0.5) * 1.0e-3;
    let t_m = input.hull_thickness_mm.max(0.05) * 1.0e-3;
    let sigma_pa = input.hull_strength_mpa.max(1.0) * 1.0e6;
    (2.0 * t_m * sigma_pa / r_m) * 1.0e-6
}

fn popcorn_strength_factor(input: PopcornInput, kernel_temperature_c: f64) -> f64 {
    let excess = (kernel_temperature_c - input.hull_softening_onset_c).max(0.0);
    let factor = 1.0 - input.hull_softening_per_c.max(0.0) * excess;
    factor.clamp(input.hull_min_strength_fraction.clamp(0.1, 1.0), 1.0)
}

pub fn evaluate_popcorn_popping(input: PopcornInput) -> PopcornResult {
    let dt = input.time_step_s.clamp(0.02, 2.0);
    let base_thresh = popcorn_rupture_threshold_mpa(input);

    let mut t_s = 0.0;
    let mut t_ext_c = input.initial_kernel_temperature_c;
    let mut t_core_c = input.initial_kernel_temperature_c;
    let mut vapor_fraction = 0.01f64;
    let mut damage = 0.0f64;

    let mut ready_temperature_c = input.max_temperature_c;
    let mut ready_time_s = f64::NAN;
    let mut ready_pressure_mpa = 0.0;
    let mut ready_thresh_mpa = 0.0;
    let mut ready_found = false;

    let mut burst_temperature_c = input.max_temperature_c;
    let mut burst_time_s = f64::NAN;
    let mut burst_pressure_mpa = 0.0;
    let mut burst_thresh_mpa = 0.0;
    let mut popped = false;

    let max_time_s =
        ((input.max_temperature_c - input.initial_kernel_temperature_c).max(0.0)
            / input.heating_rate_c_per_s.max(1.0e-6))
            * 1.5;

    while t_s <= max_time_s && t_ext_c <= input.max_temperature_c + 1.0 {
        t_ext_c += input.heating_rate_c_per_s.max(0.05) * dt;
        let tau = input.thermal_lag_tau_s.max(0.1);
        t_core_c += (t_ext_c - t_core_c) * (dt / tau);

        let t_k = (t_core_c + 273.15).max(250.0);
        let k = input.vapor_kinetics_preexp_s_inv.max(0.0)
            * (-(input.vapor_activation_energy_kj_mol.max(1.0) * 1000.0)
                / (R_GAS_J_MOL_K * t_k))
                .exp();
        let dx = k * (input.vaporized_moisture_fraction.clamp(0.05, 0.95) - vapor_fraction).max(0.0) * dt;
        vapor_fraction = (vapor_fraction + dx).clamp(0.001, input.vaporized_moisture_fraction.clamp(0.05, 0.95));

        let p_int = popcorn_internal_pressure_mpa(input, t_core_c, vapor_fraction);
        let strength_factor = popcorn_strength_factor(input, t_core_c);
        let p_thresh = base_thresh * strength_factor;

        if !ready_found && p_int >= p_thresh {
            ready_found = true;
            ready_temperature_c = t_core_c;
            ready_time_s = t_s;
            ready_pressure_mpa = p_int;
            ready_thresh_mpa = p_thresh;
        }

        let overstress = (p_int / p_thresh.max(1.0e-9) - 1.0).max(0.0);
        let damage_rate = overstress.powf(input.damage_exponent.clamp(1.0, 4.0))
            / input.damage_timescale_s.max(0.2);
        damage += damage_rate * dt;

        if ready_found && damage >= 1.0 && p_int >= p_thresh {
            popped = true;
            burst_temperature_c = t_core_c;
            burst_time_s = t_s;
            burst_pressure_mpa = p_int;
            burst_thresh_mpa = p_thresh;
            break;
        }

        t_s += dt;
    }

    if !ready_found {
        let p_int = popcorn_internal_pressure_mpa(input, t_core_c, vapor_fraction);
        let p_thresh = base_thresh * popcorn_strength_factor(input, t_core_c);
        ready_temperature_c = t_core_c;
        ready_time_s = t_s;
        ready_pressure_mpa = p_int;
        ready_thresh_mpa = p_thresh;
    }

    if !popped {
        let p_int = popcorn_internal_pressure_mpa(input, t_core_c, vapor_fraction);
        let p_thresh = base_thresh * popcorn_strength_factor(input, t_core_c);
        burst_temperature_c = t_core_c;
        burst_time_s = t_s;
        burst_pressure_mpa = p_int;
        burst_thresh_mpa = p_thresh;
    }

    let margin_ready = ready_pressure_mpa - ready_thresh_mpa;
    let margin_burst = burst_pressure_mpa - burst_thresh_mpa;
    let expansion_ratio = (2.8 * (burst_pressure_mpa / 0.101_325).powf(0.9)).clamp(5.0, 60.0);

    PopcornResult {
        ready_temperature_c,
        burst_temperature_c,
        ready_time_s,
        burst_time_s,
        internal_pressure_ready_mpa: ready_pressure_mpa,
        internal_pressure_burst_mpa: burst_pressure_mpa,
        rupture_threshold_ready_mpa: ready_thresh_mpa,
        rupture_threshold_burst_mpa: burst_thresh_mpa,
        pressure_margin_ready_mpa: margin_ready,
        pressure_margin_burst_mpa: margin_burst,
        hysteresis_delta_c: (burst_temperature_c - ready_temperature_c).max(0.0),
        estimated_expansion_ratio: expansion_ratio,
        pops: popped,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RaindropPoint {
    pub diameter_mm: f64,
    pub aspect_ratio: f64,
    pub bond_number: f64,
    pub weber_number: f64,
    pub terminal_velocity_m_s: f64,
    pub drag_coefficient: f64,
    pub transport_score: f64,
    pub stable: bool,
}

#[derive(Clone, Debug)]
pub struct RaindropSweepResult {
    pub points: Vec<RaindropPoint>,
    pub optimal: RaindropPoint,
}

fn raindrop_terminal_velocity_m_s(diameter_mm: f64) -> f64 {
    // Atlas-style empirical fit (trend-level).
    let d = diameter_mm.max(0.05);
    9.65 - 10.3 * (-0.6 * d).exp()
}

pub fn evaluate_raindrop_shape_sweep() -> RaindropSweepResult {
    let mut points = Vec::new();
    for d_mm in (2..=70).map(|i| i as f64 * 0.1) {
        let d_m = d_mm * 1.0e-3;
        let bo = WATER_DENSITY_KG_M3 * GRAVITY_M_S2 * d_m * d_m / WATER_SURFACE_TENSION_N_M;
        let aspect = 1.0 / (1.0 + 0.12 * bo.powf(0.9));
        let aspect = aspect.clamp(0.45, 1.0);
        let vt = raindrop_terminal_velocity_m_s(d_mm);
        let we = AIR_DENSITY_KG_M3 * vt * vt * d_m / WATER_SURFACE_TENSION_N_M;
        let cd = (0.47 + 0.2 * (1.0 - aspect)).clamp(0.47, 0.8);
        let risk = (-((we - 10.0).max(0.0).powi(2) / 8.0)).exp();
        let score = (d_m.powi(3) * vt / cd.max(1.0e-6)) * risk;
        let stable = we < 12.0 && d_mm < 6.0;
        points.push(RaindropPoint {
            diameter_mm: d_mm,
            aspect_ratio: aspect,
            bond_number: bo,
            weber_number: we,
            terminal_velocity_m_s: vt,
            drag_coefficient: cd,
            transport_score: score,
            stable,
        });
    }

    let optimal = points
        .iter()
        .filter(|p| p.stable)
        .max_by(|a, b| {
            a.transport_score
                .partial_cmp(&b.transport_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied()
        .unwrap_or(points[0]);

    RaindropSweepResult { points, optimal }
}

#[derive(Clone, Copy, Debug)]
pub struct MpembaInput {
    pub initial_hot_c: f64,
    pub initial_cold_c: f64,
    pub ambient_c: f64,
    pub initial_mass_kg: f64,
    pub area_m2: f64,
    pub h_base_w_m2_k: f64,
    pub evap_fraction_hot: f64,
    pub evap_fraction_cold: f64,
    pub convection_boost_hot: f64,
    pub convection_boost_cold: f64,
    pub supercool_hot_c: f64,
    pub supercool_cold_c: f64,
    pub freezing_flux_boost_hot: f64,
    pub freezing_flux_boost_cold: f64,
}

impl Default for MpembaInput {
    fn default() -> Self {
        Self {
            initial_hot_c: 80.0,
            initial_cold_c: 30.0,
            ambient_c: -18.0,
            initial_mass_kg: 0.25,
            area_m2: 0.015,
            h_base_w_m2_k: 22.0,
            evap_fraction_hot: 0.09,
            evap_fraction_cold: 0.02,
            convection_boost_hot: 0.35,
            convection_boost_cold: 0.08,
            supercool_hot_c: 0.8,
            supercool_cold_c: 5.0,
            freezing_flux_boost_hot: 0.15,
            freezing_flux_boost_cold: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct MpembaResult {
    pub hot_total_freeze_time_s: f64,
    pub cold_total_freeze_time_s: f64,
    pub hot_faster: bool,
    pub time_advantage_minutes: f64,
    pub hot_mass_after_evap_kg: f64,
    pub cold_mass_after_evap_kg: f64,
}

fn cooling_time_to_target_s(
    mass_kg: f64,
    c_p: f64,
    h_w_m2_k: f64,
    area_m2: f64,
    t_start_c: f64,
    t_target_c: f64,
    t_env_c: f64,
) -> f64 {
    let delta_start = (t_start_c - t_env_c).max(1.0e-6);
    let delta_target = (t_target_c - t_env_c).max(1.0e-6);
    let tau = mass_kg.max(1.0e-9) * c_p.max(1.0) / (h_w_m2_k.max(1.0e-9) * area_m2.max(1.0e-9));
    tau * (delta_start / delta_target).ln().max(0.0)
}

fn freezing_time_s(
    mass_kg: f64,
    latent_j_kg: f64,
    h_w_m2_k: f64,
    area_m2: f64,
    t_env_c: f64,
    flux_boost: f64,
) -> f64 {
    let qdot = h_w_m2_k.max(1.0e-9)
        * area_m2.max(1.0e-9)
        * (0.0 - t_env_c).max(1.0)
        * (1.0 + flux_boost.max(0.0));
    mass_kg.max(1.0e-9) * latent_j_kg.max(1.0) / qdot
}

pub fn evaluate_mpemba(input: MpembaInput) -> MpembaResult {
    let m_hot = input.initial_mass_kg.max(1.0e-6) * (1.0 - input.evap_fraction_hot.clamp(0.0, 0.5));
    let m_cold =
        input.initial_mass_kg.max(1.0e-6) * (1.0 - input.evap_fraction_cold.clamp(0.0, 0.5));
    let h_hot = input.h_base_w_m2_k.max(1.0e-3) * (1.0 + input.convection_boost_hot.max(0.0));
    let h_cold = input.h_base_w_m2_k.max(1.0e-3) * (1.0 + input.convection_boost_cold.max(0.0));

    let hot_target_c = -input.supercool_hot_c.max(0.0);
    let cold_target_c = -input.supercool_cold_c.max(0.0);

    let hot_cool = cooling_time_to_target_s(
        m_hot,
        WATER_CP_J_KG_K,
        h_hot,
        input.area_m2,
        input.initial_hot_c,
        hot_target_c,
        input.ambient_c,
    );
    let cold_cool = cooling_time_to_target_s(
        m_cold,
        WATER_CP_J_KG_K,
        h_cold,
        input.area_m2,
        input.initial_cold_c,
        cold_target_c,
        input.ambient_c,
    );

    let hot_freeze = freezing_time_s(
        m_hot,
        WATER_LATENT_J_KG,
        h_hot,
        input.area_m2,
        input.ambient_c,
        input.freezing_flux_boost_hot,
    );
    let cold_freeze = freezing_time_s(
        m_cold,
        WATER_LATENT_J_KG,
        h_cold,
        input.area_m2,
        input.ambient_c,
        input.freezing_flux_boost_cold,
    );

    let t_hot = hot_cool + hot_freeze;
    let t_cold = cold_cool + cold_freeze;
    let hot_faster = t_hot < t_cold;

    MpembaResult {
        hot_total_freeze_time_s: t_hot,
        cold_total_freeze_time_s: t_cold,
        hot_faster,
        time_advantage_minutes: (t_cold - t_hot) / 60.0,
        hot_mass_after_evap_kg: m_hot,
        cold_mass_after_evap_kg: m_cold,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct MpembaSweepSummary {
    pub sample_count: usize,
    pub hot_faster_count: usize,
    pub hot_faster_fraction: f64,
}

pub fn evaluate_mpemba_small_sweep() -> MpembaSweepSummary {
    let mut total = 0usize;
    let mut faster = 0usize;
    for evap_hot in [0.02, 0.05, 0.08, 0.12] {
        for supercool_cold in [1.0, 3.0, 5.0, 7.0] {
            let cfg = MpembaInput {
                evap_fraction_hot: evap_hot,
                supercool_cold_c: supercool_cold,
                ..MpembaInput::default()
            };
            let out = evaluate_mpemba(cfg);
            total += 1;
            if out.hot_faster {
                faster += 1;
            }
        }
    }
    MpembaSweepSummary {
        sample_count: total,
        hot_faster_count: faster,
        hot_faster_fraction: faster as f64 / (total as f64).max(1.0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ice_is_slipperier_near_melting() {
        let cold = evaluate_ice_slipperiness(IceSlipperinessInput {
            ice_temperature_c: -20.0,
            ..IceSlipperinessInput::default()
        });
        let warm = evaluate_ice_slipperiness(IceSlipperinessInput {
            ice_temperature_c: -2.0,
            ..IceSlipperinessInput::default()
        });
        assert!(warm.friction_coefficient < cold.friction_coefficient);
        assert!(warm.quasi_liquid_layer_thickness_nm > cold.quasi_liquid_layer_thickness_nm);
    }

    #[test]
    fn popcorn_default_has_ready_then_burst_in_observed_band() {
        let out = evaluate_popcorn_popping(PopcornInput::default());
        assert!(out.pops);
        assert!(out.ready_temperature_c >= 145.0);
        assert!(out.ready_temperature_c <= 175.0);
        assert!(out.burst_temperature_c >= 170.0);
        assert!(out.burst_temperature_c <= 195.0);
        assert!(out.burst_temperature_c > out.ready_temperature_c);
        assert!(out.hysteresis_delta_c >= 5.0);
        assert!(out.estimated_expansion_ratio > 10.0);
    }

    #[test]
    fn raindrop_optimum_is_oblate_and_stable() {
        let sweep = evaluate_raindrop_shape_sweep();
        let opt = sweep.optimal;
        assert!(opt.diameter_mm > 2.5);
        assert!(opt.diameter_mm < 6.0);
        assert!(opt.aspect_ratio < 1.0);
        assert!(opt.stable);
    }

    #[test]
    fn mpemba_needs_conditions_not_universal() {
        let default_case = evaluate_mpemba(MpembaInput::default());
        assert!(default_case.hot_faster);

        let control = evaluate_mpemba(MpembaInput {
            evap_fraction_hot: 0.0,
            evap_fraction_cold: 0.0,
            convection_boost_hot: 0.0,
            convection_boost_cold: 0.0,
            supercool_hot_c: 0.0,
            supercool_cold_c: 0.0,
            freezing_flux_boost_hot: 0.0,
            freezing_flux_boost_cold: 0.0,
            ..MpembaInput::default()
        });
        assert!(!control.hot_faster);
    }
}
