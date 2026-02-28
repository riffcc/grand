//! Blind prediction register (single-candidate lock).
//!
//! Purpose:
//! - Freeze one high-value, low-knob prediction before external measurement.
//! - Keep a machine-readable artifact that can be audited for post-hoc drift.
//!
//! Candidate selected:
//! - Neutrino hierarchy + mass character (normal ordering, Dirac-like).

use gutoe_em::{
    neutrino_absolute_masses_from_texture, neutrino_dirac_majorana_prediction,
    neutrino_hierarchy_prediction, neutrino_majorana_symmetry_residual, neutrino_texture_eigenvalues,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

const BLIND_PREDICTION_ID: &str = "BLIND-NEUTRINO-001";
const FREEZE_DATE_UTC: &str = "2026-02-28";

fn main() {
    let hierarchy = neutrino_hierarchy_prediction();
    let mass_character = neutrino_dirac_majorana_prediction();
    let majorana_residual = neutrino_majorana_symmetry_residual();
    let tex = neutrino_texture_eigenvalues();
    let dm21 = tex[1] - tex[0];
    let dm31 = tex[2] - tex[0];
    let abs = neutrino_absolute_masses_from_texture();
    let alpha = abs.alpha_physical;
    let electron_ev = abs.electron_mass_anchor_ev;
    let scale_ev = abs.mass_scale_ev;
    let m1_ev = abs.m1_ev;
    let m2_ev = abs.m2_ev;
    let m3_ev = abs.m3_ev;
    let sum_ev = abs.sum_ev;

    let out_dir = std::env::var("GUTOE_BLIND_PRED_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/blind_prediction_register".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("blind_prediction_register.txt");
    let json_path = out.join("blind_prediction_register.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[blind_prediction_register]").expect("write");
    writeln!(txt, "id = {}", BLIND_PREDICTION_ID).expect("write");
    writeln!(txt, "freeze_date_utc = {}", FREEZE_DATE_UTC).expect("write");
    writeln!(txt, "candidate = neutrino_hierarchy_mass_character").expect("write");
    writeln!(txt, "hierarchy_prediction = {}", hierarchy).expect("write");
    writeln!(txt, "mass_character_prediction = {}", mass_character).expect("write");
    writeln!(txt, "majorana_symmetry_residual = {:.12e}", majorana_residual).expect("write");
    writeln!(txt, "texture_m1 = {:.12e}", tex[0]).expect("write");
    writeln!(txt, "texture_m2 = {:.12e}", tex[1]).expect("write");
    writeln!(txt, "texture_m3 = {:.12e}", tex[2]).expect("write");
    writeln!(txt, "texture_dm21 = {:.12e}", dm21).expect("write");
    writeln!(txt, "texture_dm31 = {:.12e}", dm31).expect("write");
    writeln!(txt, "absolute_scale_ev = {:.12e}", scale_ev).expect("write");
    writeln!(txt, "absolute_hierarchy_exponent = {:.12}", abs.hierarchy_exponent).expect("write");
    writeln!(txt, "absolute_m1_ev = {:.12e}", m1_ev).expect("write");
    writeln!(txt, "absolute_m2_ev = {:.12e}", m2_ev).expect("write");
    writeln!(txt, "absolute_m3_ev = {:.12e}", m3_ev).expect("write");
    writeln!(txt, "absolute_sum_ev = {:.12e}", sum_ev).expect("write");
    writeln!(txt, "alpha_physical = {:.12}", alpha).expect("write");
    writeln!(txt, "electron_anchor_ev = {:.12e}", electron_ev).expect("write");
    writeln!(txt, "falsification_criterion_1 = measured_hierarchy != normal").expect("write");
    writeln!(
        txt,
        "falsification_criterion_2 = robust_0nu_beta_beta_detection (Dirac lane falsified)"
    )
    .expect("write");

    let payload = json!({
        "meta": {
            "id": BLIND_PREDICTION_ID,
            "freeze_date_utc": FREEZE_DATE_UTC,
            "candidate": "neutrino_hierarchy_mass_character",
            "notes": [
                "Chosen as low-knob binary candidate: normal vs inverted, dirac vs majorana-like",
                "No seeded Ki or pharmacology anchors in this lane",
                "Report is a freeze artifact; amend only via a new numbered finding"
            ]
        },
        "prediction": {
            "hierarchy_prediction": hierarchy,
            "mass_character_prediction": mass_character,
            "majorana_symmetry_residual": majorana_residual
        },
        "texture_lane": {
            "m1": tex[0],
            "m2": tex[1],
            "m3": tex[2],
            "delta_m21": dm21,
            "delta_m31": dm31
        },
        "absolute_lane": {
            "alpha": alpha,
            "electron_mass_anchor_ev": electron_ev,
            "mass_scale_ev": scale_ev,
            "hierarchy_exponent": abs.hierarchy_exponent,
            "m1_ev": m1_ev,
            "m2_ev": m2_ev,
            "m3_ev": m3_ev,
            "sum_ev": sum_ev
        },
        "falsification": {
            "criterion_1": "hierarchy_measurement_must_be_normal",
            "criterion_2": "majorana_signal_detection_falsifies_dirac_prediction"
        }
    });
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "blind_prediction_register: id={} hierarchy={} type={} sum_ev={:.6e}",
        BLIND_PREDICTION_ID, hierarchy, mass_character, sum_ev
    );
}
