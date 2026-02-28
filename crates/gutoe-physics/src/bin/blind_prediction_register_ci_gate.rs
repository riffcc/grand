//! CI freeze gate for BLIND-NEUTRINO-001.
//!
//! Ensures the registered blind candidate cannot silently drift.

use gutoe_em::{
    neutrino_absolute_masses_from_texture, neutrino_dirac_majorana_prediction,
    neutrino_hierarchy_prediction, neutrino_majorana_symmetry_residual, neutrino_texture_eigenvalues,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::process;

const EXPECTED_HIERARCHY: &str = "normal";
const EXPECTED_MASS_CHARACTER: &str = "dirac";

fn main() {
    let hierarchy = neutrino_hierarchy_prediction();
    let mass_character = neutrino_dirac_majorana_prediction();
    let majorana_residual = neutrino_majorana_symmetry_residual();
    let tex = neutrino_texture_eigenvalues();
    let dm21 = tex[1] - tex[0];
    let dm31 = tex[2] - tex[0];
    let abs = neutrino_absolute_masses_from_texture();
    let m1_ev = abs.m1_ev;
    let m2_ev = abs.m2_ev;
    let m3_ev = abs.m3_ev;
    let sum_ev = abs.sum_ev;

    let hierarchy_ok = hierarchy == EXPECTED_HIERARCHY;
    let mass_character_ok = mass_character == EXPECTED_MASS_CHARACTER;
    let overall_pass = hierarchy_ok && mass_character_ok;

    let out_dir = std::env::var("GUTOE_BLIND_PRED_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/blind_prediction_register_ci_gate.json");

    let payload = json!({
        "overall_pass": overall_pass,
        "expected": {
            "hierarchy": EXPECTED_HIERARCHY,
            "mass_character": EXPECTED_MASS_CHARACTER,
            "lock_scope": "hierarchy+mass_character_only"
        },
        "observed": {
            "hierarchy": hierarchy,
            "mass_character": mass_character,
            "majorana_symmetry_residual": majorana_residual,
            "texture_dm21": dm21,
            "texture_dm31": dm31,
            "hierarchy_exponent": abs.hierarchy_exponent,
            "m1_ev": m1_ev,
            "m2_ev": m2_ev,
            "m3_ev": m3_ev,
            "sum_ev": sum_ev
        },
        "checks": {
            "hierarchy_ok": hierarchy_ok,
            "mass_character_ok": mass_character_ok
        },
        "advisory": {
            "majorana_symmetry_residual": majorana_residual,
            "texture_dm21": dm21,
            "texture_dm31": dm31,
            "absolute_m1_ev": m1_ev,
            "absolute_m2_ev": m2_ev,
            "absolute_m3_ev": m3_ev,
            "absolute_sum_ev": sum_ev
        }
    });

    let mut f = File::create(&json_path).expect("create gate json");
    writeln!(
        f,
        "{}",
        serde_json::to_string_pretty(&payload).expect("serialize")
    )
    .expect("write gate");

    println!(
        "blind_prediction_register_gate: pass={} hierarchy={} type={} sum_ev={:.6e}",
        overall_pass, hierarchy, mass_character, sum_ev
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: hierarchy_ok={} mass_character_ok={}",
            hierarchy_ok, mass_character_ok
        );
        process::exit(2);
    }
}
