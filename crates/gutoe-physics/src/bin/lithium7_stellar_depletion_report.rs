//! Pop-II stellar lithium depletion report.
//!
//! This lane compares:
//! - required post-BBN Li-7 survival from the BBN closure transduction
//! - simulated pre-main-sequence Li-7 survival from a Pop-II envelope burn model
//!   using the stellar reaction-rate engine (`li7_burn` channel).

use gutoe_physics::{
    evaluate_lithium7_stellar_depletion_default, ENVELOPE_HYDROGEN_FRACTION,
    INTEGRATION_STEP_YEARS, LITHIUM_BURN_ONSET_TEMPERATURE_K, PRE_MAIN_SEQUENCE_YEARS,
    SOLAR_CORE_TEMPERATURE_K, SOLAR_METALLICITY_Z,
};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let out_dir = std::env::var("GUTOE_LI7_STELLAR_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/lithium7_stellar_depletion".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let report = evaluate_lithium7_stellar_depletion_default();
    let best = report.best_match;

    let txt_path = out.join("lithium7_stellar_depletion_report.txt");
    let json_path = out.join("lithium7_stellar_depletion_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[target_from_bbn]").expect("write");
    writeln!(txt, "eta10 = {:.12}", report.eta10).expect("write");
    writeln!(
        txt,
        "required_survival_factor = {:.12}",
        report.required_survival_factor
    )
    .expect("write");
    writeln!(
        txt,
        "required_depletion_percent = {:.12}",
        report.required_depletion_percent
    )
    .expect("write");
    writeln!(
        txt,
        "lithium_burn_onset_temperature_k = {:.3}",
        LITHIUM_BURN_ONSET_TEMPERATURE_K
    )
    .expect("write");
    writeln!(
        txt,
        "convective_exposure_factor_structural = {:.12}",
        report.convective_exposure_factor
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[popii_cases]").expect("write");
    for r in &report.cases {
        writeln!(
            txt,
            "{}: mass={:.3} z={:.6} core_k={:.3e} bce_k={:.3e} t9={:.6} rate_per_year={:.6e} survival={:.12} depletion_pct={:.6} closure_delta={:.12}",
            r.input.label,
            r.input.mass_solar,
            r.input.metallicity_z,
            r.core_temperature_k,
            r.convective_base_temperature_k,
            r.t9,
            r.li7_burn_rate_per_year,
            r.survival_factor,
            r.depletion_percent,
            r.closure_delta
        )
        .expect("write");
    }
    writeln!(txt).expect("write");
    writeln!(txt, "[best_match]").expect("write");
    writeln!(txt, "label = {}", best.input.label).expect("write");
    writeln!(txt, "survival_factor = {:.12}", best.survival_factor).expect("write");
    writeln!(txt, "depletion_percent = {:.12}", best.depletion_percent).expect("write");
    writeln!(txt, "closure_delta = {:.12}", best.closure_delta).expect("write");
    writeln!(
        txt,
        "agreement_with_required = {:.12}",
        report.agreement_with_required
    )
    .expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"target_from_bbn\": {{\"eta10\": {:.12}, \"required_survival_factor\": {:.12}, \"required_depletion_percent\": {:.12}}},\n  \"model\": {{\"solar_core_temperature_k\": {:.3}, \"solar_metallicity_z\": {:.6}, \"lithium_burn_onset_temperature_k\": {:.3}, \"pre_main_sequence_years\": {:.3e}, \"integration_step_years\": {:.3e}, \"envelope_hydrogen_fraction\": {:.6}, \"convective_exposure_factor_structural\": {:.12}}},\n  \"cases\": [",
        report.eta10,
        report.required_survival_factor,
        report.required_depletion_percent,
        SOLAR_CORE_TEMPERATURE_K,
        SOLAR_METALLICITY_Z,
        LITHIUM_BURN_ONSET_TEMPERATURE_K,
        PRE_MAIN_SEQUENCE_YEARS,
        INTEGRATION_STEP_YEARS,
        ENVELOPE_HYDROGEN_FRACTION,
        report.convective_exposure_factor,
    )
    .expect("write");

    for (i, r) in report.cases.iter().enumerate() {
        if i > 0 {
            writeln!(json, ",").expect("write");
        }
        write!(
            json,
            "\n    {{\"label\":\"{}\",\"mass_solar\":{:.6},\"metallicity_z\":{:.8},\"core_temperature_k\":{:.6e},\"convective_base_temperature_k\":{:.6e},\"t9\":{:.6},\"li7_burn_rate_per_year\":{:.6e},\"survival_factor\":{:.12},\"depletion_percent\":{:.12},\"closure_delta\":{:.12}}}",
            r.input.label,
            r.input.mass_solar,
            r.input.metallicity_z,
            r.core_temperature_k,
            r.convective_base_temperature_k,
            r.t9,
            r.li7_burn_rate_per_year,
            r.survival_factor,
            r.depletion_percent,
            r.closure_delta
        )
        .expect("write");
    }

    writeln!(
        json,
        "\n  ],\n  \"best_match\": {{\"label\":\"{}\", \"survival_factor\": {:.12}, \"depletion_percent\": {:.12}, \"closure_delta\": {:.12}, \"agreement_with_required\": {:.12}}}\n}}",
        best.input.label,
        best.survival_factor,
        best.depletion_percent,
        best.closure_delta,
        report.agreement_with_required
    )
    .expect("write");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "Li-7 stellar depletion: required={:.6}, best_case={} survival={:.6} (Δ={:.6})",
        report.required_survival_factor,
        best.input.label,
        best.survival_factor,
        best.closure_delta
    );
}
