//! GRAND-353 hadron transduction CI gate.

use gutoe_physics::{
    evaluate_hadron_transduction_gate, HadronReferenceAnchors, HadronTransductionWindows,
    HadronUncertaintyAssumptions,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_u64(key: &str, default: u64) -> u64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(default)
}

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn main() {
    let mut anchors = HadronReferenceAnchors::default();
    anchors.electron_mass_mev =
        env_f64("GUTOE_HADRON_ELECTRON_MASS_MEV", anchors.electron_mass_mev);
    anchors.proton_mass_mev_obs = env_f64(
        "GUTOE_HADRON_PROTON_MASS_OBS_MEV",
        anchors.proton_mass_mev_obs,
    );
    anchors.neutron_mass_mev_obs = env_f64(
        "GUTOE_HADRON_NEUTRON_MASS_OBS_MEV",
        anchors.neutron_mass_mev_obs,
    );
    anchors.pion_mass_mev_obs =
        env_f64("GUTOE_HADRON_PION_MASS_OBS_MEV", anchors.pion_mass_mev_obs);
    anchors.kaon_mass_mev_obs =
        env_f64("GUTOE_HADRON_KAON_MASS_OBS_MEV", anchors.kaon_mass_mev_obs);
    anchors.alpha_s_mz = env_f64("GUTOE_HADRON_ALPHA_S_MZ", anchors.alpha_s_mz);
    anchors.q_ref_gev = env_f64("GUTOE_HADRON_QREF_GEV", anchors.q_ref_gev);
    anchors.mb_gev = env_f64("GUTOE_HADRON_MB_GEV", anchors.mb_gev);
    anchors.mc_gev = env_f64("GUTOE_HADRON_MC_GEV", anchors.mc_gev);

    let mut assumptions = HadronUncertaintyAssumptions::default();
    assumptions.samples = env_usize("GUTOE_HADRON_SAMPLES", assumptions.samples);
    assumptions.seed = env_u64("GUTOE_HADRON_SEED", assumptions.seed);
    assumptions.electron_mass_sigma_mev = env_f64(
        "GUTOE_HADRON_ELECTRON_SIGMA_MEV",
        assumptions.electron_mass_sigma_mev,
    );
    assumptions.alpha_s_mz_sigma =
        env_f64("GUTOE_HADRON_ALPHA_S_SIGMA", assumptions.alpha_s_mz_sigma);
    assumptions.q_ref_gev_sigma =
        env_f64("GUTOE_HADRON_QREF_SIGMA_GEV", assumptions.q_ref_gev_sigma);
    assumptions.mb_gev_sigma = env_f64("GUTOE_HADRON_MB_SIGMA_GEV", assumptions.mb_gev_sigma);
    assumptions.mc_gev_sigma = env_f64("GUTOE_HADRON_MC_SIGMA_GEV", assumptions.mc_gev_sigma);
    assumptions.transduction_rel_sigma = env_f64(
        "GUTOE_HADRON_TRANSDUCTION_REL_SIGMA",
        assumptions.transduction_rel_sigma,
    );

    let mut windows = HadronTransductionWindows::default();
    windows.proton_rel_error_abs_max = env_f64(
        "GUTOE_HADRON_PROTON_REL_ERROR_ABS_MAX",
        windows.proton_rel_error_abs_max,
    );
    windows.neutron_rel_error_abs_max = env_f64(
        "GUTOE_HADRON_NEUTRON_REL_ERROR_ABS_MAX",
        windows.neutron_rel_error_abs_max,
    );
    windows.pion_rel_error_abs_max = env_f64(
        "GUTOE_HADRON_PION_REL_ERROR_ABS_MAX",
        windows.pion_rel_error_abs_max,
    );
    windows.kaon_rel_error_abs_max = env_f64(
        "GUTOE_HADRON_KAON_REL_ERROR_ABS_MAX",
        windows.kaon_rel_error_abs_max,
    );
    windows.min_valid_fraction = env_f64(
        "GUTOE_HADRON_MIN_VALID_FRACTION",
        windows.min_valid_fraction,
    );
    windows.pion_rel_span95_max = env_f64(
        "GUTOE_HADRON_PION_REL_SPAN95_MAX",
        windows.pion_rel_span95_max,
    );
    windows.kaon_rel_span95_max = env_f64(
        "GUTOE_HADRON_KAON_REL_SPAN95_MAX",
        windows.kaon_rel_span95_max,
    );

    let gate = evaluate_hadron_transduction_gate(anchors, assumptions, windows);
    let overall_pass = gate.passes_all();

    let out_dir =
        std::env::var("GUTOE_HADRON_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/hadron_transduction_ci_gate.json");

    let payload = json!({
        "overall_pass": overall_pass,
        "windows": {
            "proton_rel_error_abs_max": windows.proton_rel_error_abs_max,
            "neutron_rel_error_abs_max": windows.neutron_rel_error_abs_max,
            "pion_rel_error_abs_max": windows.pion_rel_error_abs_max,
            "kaon_rel_error_abs_max": windows.kaon_rel_error_abs_max,
            "min_valid_fraction": windows.min_valid_fraction,
            "pion_rel_span95_max": windows.pion_rel_span95_max,
            "kaon_rel_span95_max": windows.kaon_rel_span95_max
        },
        "gate": {
            "proton_rel_error_ok": gate.proton_rel_error_ok,
            "neutron_rel_error_ok": gate.neutron_rel_error_ok,
            "pion_rel_error_ok": gate.pion_rel_error_ok,
            "kaon_rel_error_ok": gate.kaon_rel_error_ok,
            "valid_fraction_ok": gate.valid_fraction_ok,
            "pion_span95_ok": gate.pion_span95_ok,
            "kaon_span95_ok": gate.kaon_span95_ok,
            "proton_obs_in_p95": gate.proton_obs_in_p95,
            "neutron_obs_in_p95": gate.neutron_obs_in_p95,
            "pion_obs_in_p95": gate.pion_obs_in_p95,
            "kaon_obs_in_p95": gate.kaon_obs_in_p95
        },
        "summary": {
            "valid_fraction": gate.score.uncertainty.valid_fraction,
            "proton_rel_error": gate.score.residuals.proton_rel_error,
            "neutron_rel_error": gate.score.residuals.neutron_rel_error,
            "pion_rel_error": gate.score.residuals.pion_rel_error,
            "kaon_rel_error": gate.score.residuals.kaon_rel_error,
            "pion_rel_span95": gate.score.uncertainty.pion_mev.rel_span95(),
            "kaon_rel_span95": gate.score.uncertainty.kaon_mev.rel_span95()
        },
        "prediction_mev": {
            "qcd_scale_nf3_mev": gate.score.central.qcd_scale_nf3_mev,
            "qcd_scale_effective_mev": gate.score.central.qcd_scale_effective_mev,
            "proton_mev": gate.score.central.proton_mev,
            "neutron_mev": gate.score.central.neutron_mev,
            "pion_mev": gate.score.central.pion_mev,
            "kaon_mev": gate.score.central.kaon_mev,
            "neutron_proton_split_mev": gate.score.central.neutron_proton_split_mev
        },
        "uncertainty": {
            "requested_samples": gate.score.uncertainty.requested_samples,
            "valid_samples": gate.score.uncertainty.valid_samples,
            "valid_fraction": gate.score.uncertainty.valid_fraction,
            "pion_mev_p05": gate.score.uncertainty.pion_mev.p05,
            "pion_mev_p50": gate.score.uncertainty.pion_mev.p50,
            "pion_mev_p95": gate.score.uncertainty.pion_mev.p95,
            "kaon_mev_p05": gate.score.uncertainty.kaon_mev.p05,
            "kaon_mev_p50": gate.score.uncertainty.kaon_mev.p50,
            "kaon_mev_p95": gate.score.uncertainty.kaon_mev.p95
        }
    });

    let mut json_file = File::create(&json_path).expect("create gate json");
    writeln!(
        json_file,
        "{}",
        serde_json::to_string_pretty(&payload).expect("serialize")
    )
    .expect("write gate json");

    println!("wrote {json_path}");
    println!("hadron_transduction_ci_gate overall_pass={overall_pass}");

    if !overall_pass {
        process::exit(1);
    }
}
