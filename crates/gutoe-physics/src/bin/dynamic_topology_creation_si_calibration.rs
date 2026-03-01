//! SI calibration for dynamic topology creation threshold.
//!
//! Converts structural gate units `(3/16)|R||T|` into Joules using a derived
//! worldsheet calibration:
//!   E_sheet = sigma_rear * (2*pi*R) * (c*T)
//! where sigma_rear comes from derived EW gradient tension with structural
//! constants.

use gutoe_physics::constants::{C, HBAR, HIGGS_QUARTIC_STRUCTURAL};
use serde_json::json;
use std::f64::consts::PI;
use std::fs;
use std::path::PathBuf;

const FC_VOID: f64 = 3.0 / 16.0;
const REAR_FACE_FACTOR: f64 = 1.0 / 10.0;
const V_EWSB_GEV: f64 = 245.3;
const GEV_TO_J: f64 = 1.602_176_634e-10;
const HBARC_GEV_M: f64 = 0.197_326_980_4e-15;

fn wall_tension_front_j_m2(delta_theta: f64, thickness_m: f64) -> f64 {
    let l_nat = thickness_m / HBARC_GEV_M;
    let sigma_gev3 = FC_VOID * V_EWSB_GEV.powi(2) * delta_theta.powi(2) / (2.0 * l_nat);
    let gev3_to_j_m2 = GEV_TO_J / HBARC_GEV_M.powi(2);
    sigma_gev3 * gev3_to_j_m2
}

fn main() {
    let out_dir = std::env::var("GUTOE_DYNAMIC_TOPOLOGY_SI_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/dynamic_topology_creation_si_calibration".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Structural thickness scale: l = 1/v (in natural units) converted to meters.
    let thickness_m = HBARC_GEV_M / V_EWSB_GEV;
    let delta_theta = PI;

    let sigma_front = wall_tension_front_j_m2(delta_theta, thickness_m);
    let sigma_rear = REAR_FACE_FACTOR * sigma_front;

    // Calibration coefficient kappa in J / (m*s):
    // E_gate_SI = kappa * ((3/16)|R||T|)
    // and E_sheet = sigma_rear*(2*pi*R)*(c*T)
    // => kappa = E_sheet / ((3/16)RT) = 2*pi*c*sigma_rear/(3/16).
    let kappa_j_per_m_s = 2.0 * PI * C * sigma_rear / FC_VOID;

    let cases = [
        (1.0_f64, 1.0_f64),
        (2.0, 10.0),
        (5.0, 60.0),
        (10.0, 60.0),
        (20.0, 50.0),
        (40.0, 120.0),
    ];

    let mut rows = Vec::new();
    for (r, t) in cases {
        let gate_units = FC_VOID * r.abs() * t.abs();
        let e_gate_si = kappa_j_per_m_s * gate_units;
        let e_sheet_direct = sigma_rear * (2.0 * PI * r.abs()) * (C * t.abs());
        rows.push(json!({
            "R_m": r,
            "T_s": t,
            "gate_units": gate_units,
            "energy_j_calibrated": e_gate_si,
            "energy_j_worldsheet_direct": e_sheet_direct,
            "match_ratio": if e_sheet_direct > 0.0 { e_gate_si / e_sheet_direct } else { f64::NAN }
        }));
    }

    let payload = json!({
        "model": {
            "gate_units": "(3/16)|R||T|",
            "calibration": "worldsheet_energy = sigma_rear*(2*pi*R)*(c*T)",
            "thickness_m": thickness_m,
            "delta_theta": delta_theta,
            "sigma_front_j_m2": sigma_front,
            "sigma_rear_j_m2": sigma_rear,
            "kappa_j_per_m_s": kappa_j_per_m_s,
            "lambda_structural": HIGGS_QUARTIC_STRUCTURAL,
            "hbar_si": HBAR
        },
        "cases": rows
    });

    let txt_path = out.join("dynamic_topology_creation_si_calibration.txt");
    let json_path = out.join("dynamic_topology_creation_si_calibration.json");

    let mut txt = String::new();
    txt.push_str("[dynamic_topology_creation_si_calibration]\n");
    txt.push_str(&format!("thickness_m = {:.12e}\n", thickness_m));
    txt.push_str(&format!("sigma_front_j_m2 = {:.12e}\n", sigma_front));
    txt.push_str(&format!("sigma_rear_j_m2 = {:.12e}\n", sigma_rear));
    txt.push_str(&format!("kappa_j_per_m_s = {:.12e}\n", kappa_j_per_m_s));
    txt.push_str("\n[cases]\n");
    for c in payload["cases"].as_array().expect("array") {
        txt.push_str(&format!(
            "R={:.3e} m, T={:.3e} s -> units={:.3e}, E={:.3e} J\n",
            c["R_m"].as_f64().unwrap_or(f64::NAN),
            c["T_s"].as_f64().unwrap_or(f64::NAN),
            c["gate_units"].as_f64().unwrap_or(f64::NAN),
            c["energy_j_calibrated"].as_f64().unwrap_or(f64::NAN)
        ));
    }

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
