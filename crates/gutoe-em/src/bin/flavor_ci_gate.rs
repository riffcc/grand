//! CKM/PMNS falsification CI gate.
//!
//! Fails fast if either direct or texture-derived observables drift outside
//! fixed PDG envelope bands.

use gutoe_em::{
    ckm_from_clifford, ckm_from_textures, pmns_from_clifford, pmns_from_textures,
    within_envelope, MixingEnvelope, MixingObservables, CKM_PDG_ENVELOPE, PMNS_PDG_ENVELOPE,
};
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn status_line(label: &str, obs: MixingObservables, env: MixingEnvelope) -> Result<String, String> {
    within_envelope(obs, env)?;
    Ok(format!(
        "{label}: pass (θ12={:.3}°, θ23={:.3}°, θ13={:.3}°, δ={:.3}°, J={:.3e})",
        obs.theta12_deg, obs.theta23_deg, obs.theta13_deg, obs.delta_deg, obs.jarlskog
    ))
}

fn main() {
    let ckm_direct = ckm_from_clifford();
    let ckm_texture = ckm_from_textures();
    let pmns_direct = pmns_from_clifford();
    let pmns_texture = pmns_from_textures();

    let checks = [
        (
            "ckm_direct",
            status_line("ckm_direct", ckm_direct, CKM_PDG_ENVELOPE),
            ckm_direct,
        ),
        (
            "ckm_texture",
            status_line("ckm_texture", ckm_texture, CKM_PDG_ENVELOPE),
            ckm_texture,
        ),
        (
            "pmns_direct",
            status_line("pmns_direct", pmns_direct, PMNS_PDG_ENVELOPE),
            pmns_direct,
        ),
        (
            "pmns_texture",
            status_line("pmns_texture", pmns_texture, PMNS_PDG_ENVELOPE),
            pmns_texture,
        ),
    ];

    let out_dir = std::env::var("GUTOE_FLAVOR_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/flavor_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    let mut overall_pass = true;
    let mut rows: Vec<String> = Vec::new();

    for (name, result, obs) in checks {
        match result {
            Ok(msg) => {
                println!("{msg}");
                rows.push(format!(
                    "    \"{name}\": {{ \"pass\": true, \"theta12_deg\": {:.12}, \"theta23_deg\": {:.12}, \"theta13_deg\": {:.12}, \"delta_deg\": {:.12}, \"jarlskog\": {:.12e} }}",
                    obs.theta12_deg, obs.theta23_deg, obs.theta13_deg, obs.delta_deg, obs.jarlskog
                ));
            }
            Err(err) => {
                overall_pass = false;
                eprintln!("{name}: FAIL - {err}");
                rows.push(format!(
                    "    \"{name}\": {{ \"pass\": false, \"error\": {:?}, \"theta12_deg\": {:.12}, \"theta23_deg\": {:.12}, \"theta13_deg\": {:.12}, \"delta_deg\": {:.12}, \"jarlskog\": {:.12e} }}",
                    err, obs.theta12_deg, obs.theta23_deg, obs.theta13_deg, obs.delta_deg, obs.jarlskog
                ));
            }
        }
    }

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n{}\n}}",
        if overall_pass { "true" } else { "false" },
        rows.join(",\n")
    )
    .expect("write gate json");

    println!("wrote {json_path}");

    if !overall_pass {
        process::exit(2);
    }
}
