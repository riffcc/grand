//! GRAND-353 hadron transduction report.
//!
//! Produces:
//! - explicit structural transduction map
//! - central p/n/pi/K predictions in MeV
//! - uncertainty summaries from anchor propagation

use gutoe_physics::{
    evaluate_hadron_transduction, HadronReferenceAnchors, HadronUncertaintyAssumptions,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

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

    let score = evaluate_hadron_transduction(anchors, assumptions);

    let out_dir = std::env::var("GUTOE_HADRON_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/hadron_transduction".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("hadron_transduction_report.txt");
    let json_path = out.join("hadron_transduction_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[anchors]").expect("write");
    writeln!(
        txt,
        "electron_mass_mev = {:.12}",
        score.anchors.electron_mass_mev
    )
    .expect("write");
    writeln!(txt, "alpha_s_mz = {:.12}", score.anchors.alpha_s_mz).expect("write");
    writeln!(txt, "q_ref_gev = {:.12}", score.anchors.q_ref_gev).expect("write");
    writeln!(txt, "mb_gev = {:.12}", score.anchors.mb_gev).expect("write");
    writeln!(txt, "mc_gev = {:.12}", score.anchors.mc_gev).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[structural_transduction_map]").expect("write");
    writeln!(
        txt,
        "mp_me_structural_ratio = {:.12}",
        score.structural.mp_me_structural_ratio
    )
    .expect("write");
    writeln!(txt, "pion_proxy = {:.12e}", score.structural.pion_proxy).expect("write");
    writeln!(
        txt,
        "pion_transduction_factor = {:.12}",
        score.structural.pion_transduction_factor
    )
    .expect("write");
    writeln!(
        txt,
        "qcd_visibility_damping_factor = {:.12}",
        score.structural.qcd_visibility_damping_factor
    )
    .expect("write");
    writeln!(
        txt,
        "delta_np_from_pion_factor = {:.12e}",
        score.structural.delta_np_from_pion_factor
    )
    .expect("write");
    writeln!(
        txt,
        "corrected_dark_to_visible_ratio = {:.12}",
        score.structural.corrected_dark_to_visible_ratio
    )
    .expect("write");
    writeln!(
        txt,
        "kaon_to_pion_factor = {:.12}",
        score.structural.kaon_to_pion_factor
    )
    .expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[central_prediction_mev]").expect("write");
    writeln!(
        txt,
        "qcd_scale_nf3_mev = {:.12}",
        score.central.qcd_scale_nf3_mev
    )
    .expect("write");
    writeln!(
        txt,
        "qcd_scale_effective_mev = {:.12}",
        score.central.qcd_scale_effective_mev
    )
    .expect("write");
    writeln!(txt, "proton_mev = {:.12}", score.central.proton_mev).expect("write");
    writeln!(txt, "neutron_mev = {:.12}", score.central.neutron_mev).expect("write");
    writeln!(txt, "pion_mev = {:.12}", score.central.pion_mev).expect("write");
    writeln!(txt, "kaon_mev = {:.12}", score.central.kaon_mev).expect("write");
    writeln!(
        txt,
        "neutron_proton_split_mev = {:.12}",
        score.central.neutron_proton_split_mev
    )
    .expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[central_residual]").expect("write");
    writeln!(
        txt,
        "proton_rel_error = {:.12e}",
        score.residuals.proton_rel_error
    )
    .expect("write");
    writeln!(
        txt,
        "neutron_rel_error = {:.12e}",
        score.residuals.neutron_rel_error
    )
    .expect("write");
    writeln!(
        txt,
        "pion_rel_error = {:.12e}",
        score.residuals.pion_rel_error
    )
    .expect("write");
    writeln!(
        txt,
        "kaon_rel_error = {:.12e}",
        score.residuals.kaon_rel_error
    )
    .expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[uncertainty]").expect("write");
    writeln!(
        txt,
        "requested_samples = {}",
        score.uncertainty.requested_samples
    )
    .expect("write");
    writeln!(txt, "valid_samples = {}", score.uncertainty.valid_samples).expect("write");
    writeln!(
        txt,
        "valid_fraction = {:.12}",
        score.uncertainty.valid_fraction
    )
    .expect("write");
    writeln!(
        txt,
        "pion_mev_p05_p50_p95 = {:.9},{:.9},{:.9}",
        score.uncertainty.pion_mev.p05,
        score.uncertainty.pion_mev.p50,
        score.uncertainty.pion_mev.p95
    )
    .expect("write");
    writeln!(
        txt,
        "kaon_mev_p05_p50_p95 = {:.9},{:.9},{:.9}",
        score.uncertainty.kaon_mev.p05,
        score.uncertainty.kaon_mev.p50,
        score.uncertainty.kaon_mev.p95
    )
    .expect("write");

    let payload = json!({
        "meta": {
            "lane": "hadron_transduction_grand_353",
            "note": "structural transduction + quantitative uncertainty, reduced-order scaffold"
        },
        "anchors": {
            "electron_mass_mev": score.anchors.electron_mass_mev,
            "proton_mass_mev_obs": score.anchors.proton_mass_mev_obs,
            "neutron_mass_mev_obs": score.anchors.neutron_mass_mev_obs,
            "pion_mass_mev_obs": score.anchors.pion_mass_mev_obs,
            "kaon_mass_mev_obs": score.anchors.kaon_mass_mev_obs,
            "alpha_s_mz": score.anchors.alpha_s_mz,
            "q_ref_gev": score.anchors.q_ref_gev,
            "mb_gev": score.anchors.mb_gev,
            "mc_gev": score.anchors.mc_gev
        },
        "structural_transduction_map": {
            "mp_me_structural_ratio": score.structural.mp_me_structural_ratio,
            "pion_proxy": score.structural.pion_proxy,
            "pion_transduction_factor": score.structural.pion_transduction_factor,
            "qcd_visibility_damping_factor": score.structural.qcd_visibility_damping_factor,
            "delta_np_from_pion_factor": score.structural.delta_np_from_pion_factor,
            "corrected_dark_to_visible_ratio": score.structural.corrected_dark_to_visible_ratio,
            "kaon_to_pion_factor": score.structural.kaon_to_pion_factor
        },
        "prediction_mev": {
            "qcd_scale_nf3_mev": score.central.qcd_scale_nf3_mev,
            "qcd_scale_effective_mev": score.central.qcd_scale_effective_mev,
            "proton_mev": score.central.proton_mev,
            "neutron_mev": score.central.neutron_mev,
            "pion_mev": score.central.pion_mev,
            "kaon_mev": score.central.kaon_mev,
            "neutron_proton_split_mev": score.central.neutron_proton_split_mev
        },
        "residuals": {
            "proton_rel_error": score.residuals.proton_rel_error,
            "neutron_rel_error": score.residuals.neutron_rel_error,
            "pion_rel_error": score.residuals.pion_rel_error,
            "kaon_rel_error": score.residuals.kaon_rel_error
        },
        "uncertainty": {
            "requested_samples": score.uncertainty.requested_samples,
            "valid_samples": score.uncertainty.valid_samples,
            "valid_fraction": score.uncertainty.valid_fraction,
            "qcd_scale_nf3_mev": {
                "mean": score.uncertainty.qcd_scale_nf3_mev.mean,
                "std": score.uncertainty.qcd_scale_nf3_mev.std,
                "p05": score.uncertainty.qcd_scale_nf3_mev.p05,
                "p50": score.uncertainty.qcd_scale_nf3_mev.p50,
                "p95": score.uncertainty.qcd_scale_nf3_mev.p95,
                "min": score.uncertainty.qcd_scale_nf3_mev.min,
                "max": score.uncertainty.qcd_scale_nf3_mev.max
            },
            "proton_mev": {
                "mean": score.uncertainty.proton_mev.mean,
                "std": score.uncertainty.proton_mev.std,
                "p05": score.uncertainty.proton_mev.p05,
                "p50": score.uncertainty.proton_mev.p50,
                "p95": score.uncertainty.proton_mev.p95,
                "min": score.uncertainty.proton_mev.min,
                "max": score.uncertainty.proton_mev.max
            },
            "neutron_mev": {
                "mean": score.uncertainty.neutron_mev.mean,
                "std": score.uncertainty.neutron_mev.std,
                "p05": score.uncertainty.neutron_mev.p05,
                "p50": score.uncertainty.neutron_mev.p50,
                "p95": score.uncertainty.neutron_mev.p95,
                "min": score.uncertainty.neutron_mev.min,
                "max": score.uncertainty.neutron_mev.max
            },
            "pion_mev": {
                "mean": score.uncertainty.pion_mev.mean,
                "std": score.uncertainty.pion_mev.std,
                "p05": score.uncertainty.pion_mev.p05,
                "p50": score.uncertainty.pion_mev.p50,
                "p95": score.uncertainty.pion_mev.p95,
                "min": score.uncertainty.pion_mev.min,
                "max": score.uncertainty.pion_mev.max
            },
            "kaon_mev": {
                "mean": score.uncertainty.kaon_mev.mean,
                "std": score.uncertainty.kaon_mev.std,
                "p05": score.uncertainty.kaon_mev.p05,
                "p50": score.uncertainty.kaon_mev.p50,
                "p95": score.uncertainty.kaon_mev.p95,
                "min": score.uncertainty.kaon_mev.min,
                "max": score.uncertainty.kaon_mev.max
            },
            "neutron_proton_split_mev": {
                "mean": score.uncertainty.neutron_proton_split_mev.mean,
                "std": score.uncertainty.neutron_proton_split_mev.std,
                "p05": score.uncertainty.neutron_proton_split_mev.p05,
                "p50": score.uncertainty.neutron_proton_split_mev.p50,
                "p95": score.uncertainty.neutron_proton_split_mev.p95,
                "min": score.uncertainty.neutron_proton_split_mev.min,
                "max": score.uncertainty.neutron_proton_split_mev.max
            }
        }
    });

    let mut json_file = File::create(&json_path).expect("create json");
    writeln!(
        json_file,
        "{}",
        serde_json::to_string_pretty(&payload).expect("serialize")
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "hadron_transduction: p={:.6} n={:.6} pi={:.6} K={:.6}",
        score.central.proton_mev,
        score.central.neutron_mev,
        score.central.pion_mev,
        score.central.kaon_mev
    );
}
