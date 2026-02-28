//! Neutrino oscillation consistency gate for absolute-mass transduction.
//!
//! This gate checks whether the current absolute neutrino masses reproduce
//! oscillation-scale mass-squared splittings.

use gutoe_em::{
    neutrino_absolute_masses_from_texture, neutrino_hierarchy_prediction,
};
use std::fs::{self, File};
use std::io::Write;
use std::process;

const SOLAR_DM21_TARGET_EV2: f64 = 7.53e-5;
const ATMOSPHERIC_DM32_TARGET_EV2: f64 = 2.453e-3;
const REL_TOL: f64 = 0.05;

fn rel_err(observed: f64, target: f64) -> f64 {
    if target.abs() < 1.0e-30 {
        0.0
    } else {
        (observed - target) / target
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_NEUTRINO_OSC_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/neutrino_oscillation_ci_gate.json");

    let hierarchy = neutrino_hierarchy_prediction();
    let abs = neutrino_absolute_masses_from_texture();
    let m1_ev = abs.m1_ev;
    let m2_ev = abs.m2_ev;
    let m3_ev = abs.m3_ev;
    let sum_ev = abs.sum_ev;

    let dm21_ev2 = abs.dm21_ev2;
    let dm32_ev2 = abs.dm32_ev2;
    let dm31_ev2 = abs.dm31_ev2;

    let dm21_rel_err = rel_err(dm21_ev2, SOLAR_DM21_TARGET_EV2);
    let dm32_rel_err = rel_err(dm32_ev2.abs(), ATMOSPHERIC_DM32_TARGET_EV2);

    let solar_ok = dm21_rel_err.abs() <= REL_TOL;
    let atmospheric_ok = dm32_rel_err.abs() <= REL_TOL;
    let ordering_ok = dm21_ev2 > 0.0 && dm32_ev2 > 0.0;
    let hierarchy_ok = hierarchy == "normal";

    let overall_pass = solar_ok && atmospheric_ok && ordering_ok && hierarchy_ok;

    let mut json = File::create(&json_path).expect("create neutrino oscillation gate json");
    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"windows\": {{\"solar_dm21_target_ev2\": {:.12e}, \"atmospheric_dm32_target_ev2\": {:.12e}, \"relative_tolerance\": {:.6}}},\n  \"observed\": {{\"hierarchy\": \"{}\", \"hierarchy_exponent\": {:.12}, \"m1_ev\": {:.12e}, \"m2_ev\": {:.12e}, \"m3_ev\": {:.12e}, \"sum_ev\": {:.12e}, \"dm21_ev2\": {:.12e}, \"dm32_ev2\": {:.12e}, \"dm31_ev2\": {:.12e}, \"dm21_rel_err\": {:.12e}, \"dm32_rel_err\": {:.12e}}},\n  \"checks\": {{\"hierarchy_ok\": {}, \"ordering_ok\": {}, \"solar_ok\": {}, \"atmospheric_ok\": {}}}\n}}",
        if overall_pass { "true" } else { "false" },
        SOLAR_DM21_TARGET_EV2,
        ATMOSPHERIC_DM32_TARGET_EV2,
        REL_TOL,
        hierarchy,
        abs.hierarchy_exponent,
        m1_ev,
        m2_ev,
        m3_ev,
        sum_ev,
        dm21_ev2,
        dm32_ev2,
        dm31_ev2,
        dm21_rel_err,
        dm32_rel_err,
        if hierarchy_ok { "true" } else { "false" },
        if ordering_ok { "true" } else { "false" },
        if solar_ok { "true" } else { "false" },
        if atmospheric_ok { "true" } else { "false" }
    )
    .expect("write neutrino oscillation gate json");

    println!(
        "neutrino_oscillation_ci_gate: pass={} hierarchy={} dm21={:.3e} dm32={:.3e}",
        overall_pass, hierarchy, dm21_ev2, dm32_ev2
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: hierarchy_ok={} ordering_ok={} solar_ok={} atmospheric_ok={} dm21_rel_err={:.3e} dm32_rel_err={:.3e}",
            hierarchy_ok, ordering_ok, solar_ok, atmospheric_ok, dm21_rel_err, dm32_rel_err
        );
        process::exit(2);
    }
}
