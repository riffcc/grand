//! GRAND-346: dark-sector phenomenology harness (particle + geometric branches).

use gutoe_physics::constants::{
    DARK_FRACTION_GEOMETRIC_STRUCTURAL, DARK_FRACTION_TOTAL_STATE_SPLIT,
    DARK_GEOMETRIC_AMPLIFICATION, DARK_TO_VISIBLE_COUNT_RATIO, DARK_TO_VISIBLE_GEOMETRIC_RATIO,
};
use gutoe_physics::dark_sector::{
    circular_velocity, curvature_factor_from_einstein_cosmology, dark_density, enclosed_mass_constant_density,
    lensing_deflection, DarkSectorBranch,
};
use std::fs::{self, File};
use std::io::Write;

const OMEGA_BARYON_OBS: f64 = 0.0493;
const OMEGA_DM_OBS: f64 = 0.264;
const OMEGA_MATTER_OBS: f64 = OMEGA_BARYON_OBS + OMEGA_DM_OBS;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let rho_visible = env_f64("GUTOE_DM_RHO_VISIBLE", 1.0e-21);
    let radius = env_f64("GUTOE_DM_RADIUS_M", 3.0e20);
    let impact = env_f64("GUTOE_DM_IMPACT_M", 3.0e20);

    let kappa = curvature_factor_from_einstein_cosmology(rho_visible, radius);
    let rho_dark_particle = dark_density(DarkSectorBranch::Particle, rho_visible, 1.0);
    let rho_dark_geometric = dark_density(DarkSectorBranch::Geometric, rho_visible, kappa);
    let rho_dark_unified = dark_density(DarkSectorBranch::Unified, rho_visible, kappa);

    let m_visible = enclosed_mass_constant_density(rho_visible, radius);
    let m_total_particle = enclosed_mass_constant_density(rho_visible + rho_dark_particle, radius);
    let m_total_geometric = enclosed_mass_constant_density(rho_visible + rho_dark_geometric, radius);
    let m_total_unified = enclosed_mass_constant_density(rho_visible + rho_dark_unified, radius);

    let v_visible = circular_velocity(m_visible, radius).unwrap_or(f64::NAN);
    let v_particle = circular_velocity(m_total_particle, radius).unwrap_or(f64::NAN);
    let v_geometric = circular_velocity(m_total_geometric, radius).unwrap_or(f64::NAN);
    let v_unified = circular_velocity(m_total_unified, radius).unwrap_or(f64::NAN);

    let alpha_visible = lensing_deflection(m_visible, impact).unwrap_or(f64::NAN);
    let alpha_particle = lensing_deflection(m_total_particle, impact).unwrap_or(f64::NAN);
    let alpha_geometric = lensing_deflection(m_total_geometric, impact).unwrap_or(f64::NAN);
    let alpha_unified = lensing_deflection(m_total_unified, impact).unwrap_or(f64::NAN);

    let omega_dm_particle = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_COUNT_RATIO;
    let omega_m_particle = OMEGA_BARYON_OBS + omega_dm_particle;
    let dm_fraction_particle = omega_dm_particle / omega_m_particle;
    let omega_dm_geometric = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m_geometric = OMEGA_BARYON_OBS + omega_dm_geometric;
    let dm_fraction_geometric = omega_dm_geometric / omega_m_geometric;
    let dm_fraction_geometric_with_curvature = rho_dark_geometric / (rho_visible + rho_dark_geometric);
    let dm_fraction_unified_local = rho_dark_unified / (rho_visible + rho_dark_unified);
    let dm_fraction_obs = OMEGA_DM_OBS / OMEGA_MATTER_OBS;

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let txt_path = format!("{out_dir}/dark_matter_report.txt");
    let json_path = format!("{out_dir}/dark_matter_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "GRAND-346 dark-sector harness").expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[structural_split]").expect("write");
    writeln!(txt, "dark_to_visible_ratio = {:.12}", DARK_TO_VISIBLE_COUNT_RATIO).expect("write");
    writeln!(txt, "dark_fraction_total_split = {:.12}", DARK_FRACTION_TOTAL_STATE_SPLIT).expect("write");
    writeln!(txt, "dark_geometric_amplification = {:.12}", DARK_GEOMETRIC_AMPLIFICATION).expect("write");
    writeln!(txt, "dark_to_visible_geometric_ratio = {:.12}", DARK_TO_VISIBLE_GEOMETRIC_RATIO).expect("write");
    writeln!(txt, "dark_fraction_geometric_structural = {:.12}", DARK_FRACTION_GEOMETRIC_STRUCTURAL)
        .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[inputs]").expect("write");
    writeln!(txt, "rho_visible = {:.6e} kg/m^3", rho_visible).expect("write");
    writeln!(txt, "radius = {:.6e} m", radius).expect("write");
    writeln!(txt, "impact = {:.6e} m", impact).expect("write");
    writeln!(txt, "kappa_einstein_cosmology = {:.9}", kappa).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[branch_densities]").expect("write");
    writeln!(txt, "rho_dark_particle = {:.6e}", rho_dark_particle).expect("write");
    writeln!(txt, "rho_dark_geometric = {:.6e}", rho_dark_geometric).expect("write");
    writeln!(txt, "rho_dark_unified_local = {:.6e}", rho_dark_unified).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[rotation_proxy]").expect("write");
    writeln!(txt, "v_visible = {:.6e} m/s", v_visible).expect("write");
    writeln!(txt, "v_particle = {:.6e} m/s", v_particle).expect("write");
    writeln!(txt, "v_geometric = {:.6e} m/s", v_geometric).expect("write");
    writeln!(txt, "v_unified_local = {:.6e} m/s", v_unified).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[lensing_proxy]").expect("write");
    writeln!(txt, "alpha_visible = {:.6e} rad", alpha_visible).expect("write");
    writeln!(txt, "alpha_particle = {:.6e} rad", alpha_particle).expect("write");
    writeln!(txt, "alpha_geometric = {:.6e} rad", alpha_geometric).expect("write");
    writeln!(txt, "alpha_unified_local = {:.6e} rad", alpha_unified).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[cmb_matter_fraction_check]").expect("write");
    writeln!(txt, "omega_baryon_obs = {:.9}", OMEGA_BARYON_OBS).expect("write");
    writeln!(txt, "omega_dm_obs = {:.9}", OMEGA_DM_OBS).expect("write");
    writeln!(txt, "omega_dm_particle_from_ratio = {:.9}", omega_dm_particle).expect("write");
    writeln!(txt, "omega_dm_geometric_from_ratio = {:.9}", omega_dm_geometric).expect("write");
    writeln!(txt, "dm_fraction_obs = {:.9}", dm_fraction_obs).expect("write");
    writeln!(txt, "dm_fraction_particle = {:.9}", dm_fraction_particle).expect("write");
    writeln!(txt, "dm_fraction_geometric = {:.9}", dm_fraction_geometric).expect("write");
    writeln!(
        txt,
        "dm_fraction_geometric_with_curvature = {:.9}",
        dm_fraction_geometric_with_curvature
    )
    .expect("write");
    writeln!(txt, "dm_fraction_unified_local = {:.9}", dm_fraction_unified_local).expect("write");
    writeln!(
        txt,
        "dm_fraction_particle_delta = {:.9}",
        dm_fraction_particle - dm_fraction_obs
    )
    .expect("write");
    writeln!(
        txt,
        "dm_fraction_geometric_delta = {:.9}",
        dm_fraction_geometric - dm_fraction_obs
    )
    .expect("write");
    writeln!(
        txt,
        "dm_fraction_geometric_with_curvature_delta = {:.9}",
        dm_fraction_geometric_with_curvature - dm_fraction_obs
    )
    .expect("write");
    writeln!(
        txt,
        "dm_fraction_unified_local_delta = {:.9}",
        dm_fraction_unified_local - dm_fraction_obs
    )
    .expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"structural_split\": {{\"dark_to_visible_ratio\": {:.12}, \"dark_fraction_total_split\": {:.12}}},",
        DARK_TO_VISIBLE_COUNT_RATIO, DARK_FRACTION_TOTAL_STATE_SPLIT
    )
    .expect("write");
    writeln!(
        json,
        "  \"structural_geometric\": {{\"amplification\": {:.12}, \"dark_to_visible_ratio\": {:.12}, \"dark_fraction\": {:.12}}},",
        DARK_GEOMETRIC_AMPLIFICATION,
        DARK_TO_VISIBLE_GEOMETRIC_RATIO,
        DARK_FRACTION_GEOMETRIC_STRUCTURAL
    )
    .expect("write");
    writeln!(
        json,
        "  \"inputs\": {{\"rho_visible\": {:.12e}, \"radius_m\": {:.12e}, \"impact_m\": {:.12e}, \"kappa_einstein_cosmology\": {:.12}}},",
        rho_visible, radius, impact, kappa
    )
    .expect("write");
    writeln!(
        json,
        "  \"densities\": {{\"rho_dark_particle\": {:.12e}, \"rho_dark_geometric\": {:.12e}}},",
        rho_dark_particle, rho_dark_geometric
    )
    .expect("write");
    writeln!(
        json,
        "  \"unified_local_density\": {{\"rho_dark_unified\": {:.12e}}},",
        rho_dark_unified
    )
    .expect("write");
    writeln!(
        json,
        "  \"rotation_proxy\": {{\"v_visible_m_s\": {:.12e}, \"v_particle_m_s\": {:.12e}, \"v_geometric_m_s\": {:.12e}}},",
        v_visible, v_particle, v_geometric
    )
    .expect("write");
    writeln!(
        json,
        "  \"rotation_unified_local\": {{\"v_unified_m_s\": {:.12e}}},",
        v_unified
    )
    .expect("write");
    writeln!(
        json,
        "  \"lensing_proxy\": {{\"alpha_visible_rad\": {:.12e}, \"alpha_particle_rad\": {:.12e}, \"alpha_geometric_rad\": {:.12e}}},",
        alpha_visible, alpha_particle, alpha_geometric
    )
    .expect("write");
    writeln!(
        json,
        "  \"lensing_unified_local\": {{\"alpha_unified_rad\": {:.12e}}},",
        alpha_unified
    )
    .expect("write");
    writeln!(
        json,
        "  \"cmb_check\": {{\"omega_baryon_obs\": {:.12}, \"omega_dm_obs\": {:.12}, \"omega_dm_particle\": {:.12}, \"omega_dm_geometric\": {:.12}, \"dm_fraction_obs\": {:.12}, \"dm_fraction_particle\": {:.12}, \"dm_fraction_geometric\": {:.12}, \"dm_fraction_geometric_with_curvature\": {:.12}, \"dm_fraction_unified_local\": {:.12}, \"dm_fraction_particle_delta\": {:.12}, \"dm_fraction_geometric_delta\": {:.12}, \"dm_fraction_geometric_with_curvature_delta\": {:.12}, \"dm_fraction_unified_local_delta\": {:.12}}}\n}}",
        OMEGA_BARYON_OBS,
        OMEGA_DM_OBS,
        omega_dm_particle,
        omega_dm_geometric,
        dm_fraction_obs,
        dm_fraction_particle,
        dm_fraction_geometric,
        dm_fraction_geometric_with_curvature,
        dm_fraction_unified_local,
        dm_fraction_particle - dm_fraction_obs,
        dm_fraction_geometric - dm_fraction_obs,
        dm_fraction_geometric_with_curvature - dm_fraction_obs,
        dm_fraction_unified_local - dm_fraction_obs
    )
    .expect("write");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
}
