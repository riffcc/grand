//! CKM/PMNS falsification CI gate.
//!
//! Fails fast if either direct or texture-derived observables drift outside
//! fixed PDG envelope bands.

use gutoe_em::{
    ckm_from_clifford, ckm_from_textures, cp_violation_witness, pmns_from_clifford,
    pmns_from_clifford_theta23_alpha2, pmns_from_textures, residuals, within_envelope,
    MixingEnvelope, MixingObservables, CKM_CP_J_MIN, CKM_PDG_ENVELOPE, CP_PHASE_TOL_DEG,
    PMNS_CP_J_MIN, PMNS_PDG_ENVELOPE, PMNS_TARGET, PMNS_THETA23_ALPHA2_COEFF_STRUCTURAL,
};
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::process;
use std::time::{SystemTime, UNIX_EPOCH};

fn status_line(
    label: &str,
    obs: MixingObservables,
    env: MixingEnvelope,
    cp_j_min: f64,
) -> Result<String, String> {
    within_envelope(obs, env)?;
    cp_violation_witness(obs, cp_j_min, CP_PHASE_TOL_DEG)?;
    Ok(format!(
        "{label}: pass (θ12={:.3}°, θ23={:.3}°, θ13={:.3}°, δ={:.3}°, J={:.3e})",
        obs.theta12_deg, obs.theta23_deg, obs.theta13_deg, obs.delta_deg, obs.jarlskog
    ))
}

fn unix_timestamp() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0)
}

fn append_gate_trend_row(
    out_dir: &str,
    ts: i64,
    overall_pass: bool,
    th23_improvement_ok: bool,
    rd: f64,
    rc: f64,
    ckm_direct_delta_deg: f64,
    ckm_texture_delta_deg: f64,
    pmns_direct_delta_deg: f64,
    pmns_texture_delta_deg: f64,
) {
    let header = "timestamp_unix,overall_pass,pmns_theta23_improvement_pass,pmns_theta23_direct_abs_residual_deg,pmns_theta23_corr_abs_residual_deg,ckm_direct_delta_deg,ckm_texture_delta_deg,ckm_texture_delta_drift_deg,pmns_direct_delta_deg,pmns_texture_delta_deg,pmns_texture_delta_drift_deg";
    let trend_path = format!("{out_dir}/flavor_ci_gate_trend.csv");
    let trend_exists = Path::new(&trend_path).exists();
    let mut needs_header = !trend_exists;

    if trend_exists {
        let existing_header = fs::read_to_string(&trend_path)
            .ok()
            .and_then(|s| s.lines().next().map(|line| line.trim().to_string()))
            .unwrap_or_default();
        if existing_header != header {
            let legacy_path = format!("{out_dir}/flavor_ci_gate_trend.legacy_{ts}.csv");
            let _ = fs::rename(&trend_path, &legacy_path);
            needs_header = true;
        }
    }

    let mut trend = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&trend_path)
        .expect("open flavor_ci_gate_trend.csv");

    if needs_header {
        writeln!(trend, "{header}").expect("write trend header");
    }

    writeln!(
        trend,
        "{},{},{},{:.12},{:.12},{:.12},{:.12},{:.12},{:.12},{:.12},{:.12}",
        ts,
        if overall_pass { "true" } else { "false" },
        if th23_improvement_ok { "true" } else { "false" },
        rd,
        rc,
        ckm_direct_delta_deg,
        ckm_texture_delta_deg,
        ckm_texture_delta_deg - ckm_direct_delta_deg,
        pmns_direct_delta_deg,
        pmns_texture_delta_deg,
        pmns_texture_delta_deg - pmns_direct_delta_deg,
    )
    .expect("write trend row");
}

fn main() {
    let ckm_direct = ckm_from_clifford();
    let ckm_texture = ckm_from_textures();
    let pmns_direct = pmns_from_clifford();
    let pmns_corr_c = std::env::var("GUTOE_PMNS_TH23_ALPHA2_C")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(PMNS_THETA23_ALPHA2_COEFF_STRUCTURAL);
    let pmns_direct_corr = pmns_from_clifford_theta23_alpha2(pmns_corr_c);
    let pmns_texture = pmns_from_textures();

    let checks = [
        (
            "ckm_direct",
            status_line("ckm_direct", ckm_direct, CKM_PDG_ENVELOPE, CKM_CP_J_MIN),
            ckm_direct,
        ),
        (
            "ckm_texture",
            status_line("ckm_texture", ckm_texture, CKM_PDG_ENVELOPE, CKM_CP_J_MIN),
            ckm_texture,
        ),
        (
            "pmns_direct",
            status_line("pmns_direct", pmns_direct, PMNS_PDG_ENVELOPE, PMNS_CP_J_MIN),
            pmns_direct,
        ),
        (
            "pmns_direct_theta23_alpha2_corrected",
            status_line(
                "pmns_direct_theta23_alpha2_corrected",
                pmns_direct_corr,
                PMNS_PDG_ENVELOPE,
                PMNS_CP_J_MIN,
            ),
            pmns_direct_corr,
        ),
        (
            "pmns_texture",
            status_line(
                "pmns_texture",
                pmns_texture,
                PMNS_PDG_ENVELOPE,
                PMNS_CP_J_MIN,
            ),
            pmns_texture,
        ),
    ];

    let out_dir =
        std::env::var("GUTOE_FLAVOR_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
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

    // Hard improvement gate: corrected theta23 residual must improve 10x vs direct.
    let rd = residuals(pmns_direct, PMNS_TARGET).d_theta23_deg.abs();
    let rc = residuals(pmns_direct_corr, PMNS_TARGET).d_theta23_deg.abs();
    let th23_improvement_ok = rc <= rd / 10.0;
    if !th23_improvement_ok {
        overall_pass = false;
        eprintln!(
            "pmns_theta23_improvement: FAIL - corrected residual {:.9} not <= direct/10 {:.9}",
            rc,
            rd / 10.0
        );
    } else {
        println!(
            "pmns_theta23_improvement: pass (direct={:.6}°, corrected={:.6}°, c={:.6})",
            rd, rc, pmns_corr_c
        );
    }

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"pmns_theta23_improvement\": {{\"pass\": {}, \"c_alpha2\": {:.12}, \"direct_abs_residual_deg\": {:.12}, \"corrected_abs_residual_deg\": {:.12}}},\n{}\n}}",
        if overall_pass { "true" } else { "false" },
        if th23_improvement_ok { "true" } else { "false" },
        pmns_corr_c,
        rd,
        rc,
        rows.join(",\n")
    )
    .expect("write gate json");

    let ts = unix_timestamp();
    let json_snapshot_path = format!("{out_dir}/flavor_ci_gate.{ts}.json");
    let _ = fs::copy(&json_path, &json_snapshot_path);
    append_gate_trend_row(
        &out_dir,
        ts,
        overall_pass,
        th23_improvement_ok,
        rd,
        rc,
        ckm_direct.delta_deg,
        ckm_texture.delta_deg,
        pmns_direct.delta_deg,
        pmns_texture.delta_deg,
    );

    println!("wrote {json_path}");
    println!("snapshotted {json_snapshot_path}");
    println!("appended {out_dir}/flavor_ci_gate_trend.csv");

    if !overall_pass {
        process::exit(2);
    }
}
