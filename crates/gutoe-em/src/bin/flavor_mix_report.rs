/*!
 * CKM/PMNS observables report from Cl(1,3) algebraic definitions.
 */

use gutoe_em::{
    ckm_from_clifford, ckm_from_textures, pmns_from_clifford, pmns_from_clifford_theta23_alpha2,
    pmns_from_textures, residuals, PMNS_THETA23_ALPHA2_COEFF_STRUCTURAL, CKM_TARGET, PMNS_TARGET,
};
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

fn write_block(
    f: &mut File,
    name: &str,
    theta12: f64,
    theta23: f64,
    theta13: f64,
    delta: f64,
    j: f64,
    d12: f64,
    d23: f64,
    d13: f64,
    ddelta: f64,
    dj: f64,
) {
    writeln!(f, "[{name}]").expect("write section");
    writeln!(f, "theta12_deg = {theta12:.9}").expect("write theta12");
    writeln!(f, "theta23_deg = {theta23:.9}").expect("write theta23");
    writeln!(f, "theta13_deg = {theta13:.9}").expect("write theta13");
    writeln!(f, "delta_deg = {delta:.9}").expect("write delta");
    writeln!(f, "jarlskog = {j:.12e}").expect("write J");
    writeln!(f, "delta_theta12_deg = {d12:.9}").expect("write d12");
    writeln!(f, "delta_theta23_deg = {d23:.9}").expect("write d23");
    writeln!(f, "delta_theta13_deg = {d13:.9}").expect("write d13");
    writeln!(f, "delta_delta_deg = {ddelta:.9}").expect("write ddelta");
    writeln!(f, "delta_jarlskog = {dj:.12e}").expect("write dJ");
    writeln!(f).expect("newline");
}

fn unix_timestamp() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0)
}

fn append_trend_row(
    out_dir: &str,
    ts: i64,
    ckm_direct_delta_deg: f64,
    ckm_texture_delta_deg: f64,
    pmns_direct_delta_deg: f64,
    pmns_texture_delta_deg: f64,
    pmns_theta23_direct_abs: f64,
    pmns_theta23_corr_abs: f64,
    pmns_theta23_improvement_factor: f64,
    pmns_theta23_improves_10x: bool,
) {
    let header = "timestamp_unix,ckm_direct_delta_deg,ckm_texture_delta_deg,ckm_texture_delta_drift_deg,pmns_direct_delta_deg,pmns_texture_delta_deg,pmns_texture_delta_drift_deg,pmns_theta23_direct_abs_residual_deg,pmns_theta23_corr_abs_residual_deg,pmns_theta23_improvement_factor,pmns_theta23_improves_10x";
    let trend_path = format!("{out_dir}/flavor_mix_trend.csv");
    let trend_exists = Path::new(&trend_path).exists();
    let mut needs_header = !trend_exists;

    if trend_exists {
        let existing_header = fs::read_to_string(&trend_path)
            .ok()
            .and_then(|s| s.lines().next().map(|line| line.trim().to_string()))
            .unwrap_or_default();
        if existing_header != header {
            let legacy_path = format!("{out_dir}/flavor_mix_trend.legacy_{ts}.csv");
            let _ = fs::rename(&trend_path, &legacy_path);
            needs_header = true;
        }
    }

    let mut trend = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&trend_path)
        .expect("open flavor_mix_trend.csv");

    if needs_header {
        writeln!(trend, "{header}").expect("write trend header");
    }

    let ckm_texture_delta_drift = ckm_texture_delta_deg - ckm_direct_delta_deg;
    let pmns_texture_delta_drift = pmns_texture_delta_deg - pmns_direct_delta_deg;
    writeln!(
        trend,
        "{ts},{ckm_direct_delta_deg:.12},{ckm_texture_delta_deg:.12},{ckm_texture_delta_drift:.12},{pmns_direct_delta_deg:.12},{pmns_texture_delta_deg:.12},{pmns_texture_delta_drift:.12},{pmns_theta23_direct_abs:.12},{pmns_theta23_corr_abs:.12},{pmns_theta23_improvement_factor:.12},{}",
        if pmns_theta23_improves_10x { "true" } else { "false" }
    )
    .expect("write trend row");
}

fn main() {
    let pmns_th23_corr_c = std::env::var("GUTOE_PMNS_TH23_ALPHA2_C")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(PMNS_THETA23_ALPHA2_COEFF_STRUCTURAL);

    let ckm = ckm_from_clifford();
    let pmns = pmns_from_clifford();
    let pmns_corr = pmns_from_clifford_theta23_alpha2(pmns_th23_corr_c);
    let ckm_tex = ckm_from_textures();
    let pmns_tex = pmns_from_textures();
    let ckm_r = residuals(ckm, CKM_TARGET);
    let pmns_r = residuals(pmns, PMNS_TARGET);
    let pmns_corr_r = residuals(pmns_corr, PMNS_TARGET);
    let pmns_theta23_direct_abs = pmns_r.d_theta23_deg.abs();
    let pmns_theta23_corr_abs = pmns_corr_r.d_theta23_deg.abs();
    let pmns_theta23_improvement_factor = if pmns_theta23_corr_abs > 0.0 {
        pmns_theta23_direct_abs / pmns_theta23_corr_abs
    } else {
        f64::INFINITY
    };
    let pmns_theta23_improves_10x = pmns_theta23_corr_abs <= pmns_theta23_direct_abs / 10.0;
    let ckm_tex_r = residuals(ckm_tex, CKM_TARGET);
    let pmns_tex_r = residuals(pmns_tex, PMNS_TARGET);

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let txt_path = format!("{out_dir}/flavor_mix_report.txt");
    let json_path = format!("{out_dir}/flavor_mix_report.json");

    let mut txt = File::create(&txt_path).expect("create report txt");
    write_block(
        &mut txt,
        "CKM (direct algebraic)",
        ckm.theta12_deg,
        ckm.theta23_deg,
        ckm.theta13_deg,
        ckm.delta_deg,
        ckm.jarlskog,
        ckm_r.d_theta12_deg,
        ckm_r.d_theta23_deg,
        ckm_r.d_theta13_deg,
        ckm_r.d_delta_deg,
        ckm_r.d_jarlskog,
    );
    write_block(
        &mut txt,
        "PMNS (direct algebraic)",
        pmns.theta12_deg,
        pmns.theta23_deg,
        pmns.theta13_deg,
        pmns.delta_deg,
        pmns.jarlskog,
        pmns_r.d_theta12_deg,
        pmns_r.d_theta23_deg,
        pmns_r.d_theta13_deg,
        pmns_r.d_delta_deg,
        pmns_r.d_jarlskog,
    );
    write_block(
        &mut txt,
        "PMNS (direct algebraic, theta23 alpha^2 corrected)",
        pmns_corr.theta12_deg,
        pmns_corr.theta23_deg,
        pmns_corr.theta13_deg,
        pmns_corr.delta_deg,
        pmns_corr.jarlskog,
        pmns_corr_r.d_theta12_deg,
        pmns_corr_r.d_theta23_deg,
        pmns_corr_r.d_theta13_deg,
        pmns_corr_r.d_delta_deg,
        pmns_corr_r.d_jarlskog,
    );
    writeln!(txt, "[PMNS theta23 alpha^2 improvement]").expect("write section");
    writeln!(txt, "c_alpha2 = {pmns_th23_corr_c:.12}").expect("write c");
    writeln!(txt, "direct_abs_residual_deg = {pmns_theta23_direct_abs:.9}").expect("write rd");
    writeln!(txt, "corrected_abs_residual_deg = {pmns_theta23_corr_abs:.9}").expect("write rc");
    writeln!(
        txt,
        "improvement_factor = {pmns_theta23_improvement_factor:.6}"
    )
    .expect("write factor");
    writeln!(txt, "improves_10x = {pmns_theta23_improves_10x}").expect("write bool");
    writeln!(txt).expect("newline");
    write_block(
        &mut txt,
        "CKM (texture diagonalization)",
        ckm_tex.theta12_deg,
        ckm_tex.theta23_deg,
        ckm_tex.theta13_deg,
        ckm_tex.delta_deg,
        ckm_tex.jarlskog,
        ckm_tex_r.d_theta12_deg,
        ckm_tex_r.d_theta23_deg,
        ckm_tex_r.d_theta13_deg,
        ckm_tex_r.d_delta_deg,
        ckm_tex_r.d_jarlskog,
    );
    write_block(
        &mut txt,
        "PMNS (texture diagonalization)",
        pmns_tex.theta12_deg,
        pmns_tex.theta23_deg,
        pmns_tex.theta13_deg,
        pmns_tex.delta_deg,
        pmns_tex.jarlskog,
        pmns_tex_r.d_theta12_deg,
        pmns_tex_r.d_theta23_deg,
        pmns_tex_r.d_theta13_deg,
        pmns_tex_r.d_delta_deg,
        pmns_tex_r.d_jarlskog,
    );

    let mut json = File::create(&json_path).expect("create report json");
    writeln!(
        json,
        "{{\n  \"ckm\": {{\n    \"direct\": {{\n      \"theta12_deg\": {:.12},\n      \"theta23_deg\": {:.12},\n      \"theta13_deg\": {:.12},\n      \"delta_deg\": {:.12},\n      \"jarlskog\": {:.12e},\n      \"delta_theta12_deg\": {:.12},\n      \"delta_theta23_deg\": {:.12},\n      \"delta_theta13_deg\": {:.12},\n      \"delta_delta_deg\": {:.12},\n      \"delta_jarlskog\": {:.12e}\n    }},\n    \"texture\": {{\n      \"theta12_deg\": {:.12},\n      \"theta23_deg\": {:.12},\n      \"theta13_deg\": {:.12},\n      \"delta_deg\": {:.12},\n      \"jarlskog\": {:.12e},\n      \"delta_theta12_deg\": {:.12},\n      \"delta_theta23_deg\": {:.12},\n      \"delta_theta13_deg\": {:.12},\n      \"delta_delta_deg\": {:.12},\n      \"delta_jarlskog\": {:.12e}\n    }}\n  }},\n  \"pmns\": {{\n    \"direct\": {{\n      \"theta12_deg\": {:.12},\n      \"theta23_deg\": {:.12},\n      \"theta13_deg\": {:.12},\n      \"delta_deg\": {:.12},\n      \"jarlskog\": {:.12e},\n      \"delta_theta12_deg\": {:.12},\n      \"delta_theta23_deg\": {:.12},\n      \"delta_theta13_deg\": {:.12},\n      \"delta_delta_deg\": {:.12},\n      \"delta_jarlskog\": {:.12e}\n    }},\n    \"direct_theta23_alpha2_corrected\": {{\n      \"c_alpha2\": {:.12},\n      \"theta12_deg\": {:.12},\n      \"theta23_deg\": {:.12},\n      \"theta13_deg\": {:.12},\n      \"delta_deg\": {:.12},\n      \"jarlskog\": {:.12e},\n      \"delta_theta12_deg\": {:.12},\n      \"delta_theta23_deg\": {:.12},\n      \"delta_theta13_deg\": {:.12},\n      \"delta_delta_deg\": {:.12},\n      \"delta_jarlskog\": {:.12e}\n    }},\n    \"theta23_alpha2_improvement\": {{\n      \"c_alpha2\": {:.12},\n      \"direct_abs_residual_deg\": {:.12},\n      \"corrected_abs_residual_deg\": {:.12},\n      \"improvement_factor\": {:.12},\n      \"improves_10x\": {}\n    }},\n    \"texture\": {{\n      \"theta12_deg\": {:.12},\n      \"theta23_deg\": {:.12},\n      \"theta13_deg\": {:.12},\n      \"delta_deg\": {:.12},\n      \"jarlskog\": {:.12e},\n      \"delta_theta12_deg\": {:.12},\n      \"delta_theta23_deg\": {:.12},\n      \"delta_theta13_deg\": {:.12},\n      \"delta_delta_deg\": {:.12},\n      \"delta_jarlskog\": {:.12e}\n    }}\n  }}\n}}",
        ckm.theta12_deg,
        ckm.theta23_deg,
        ckm.theta13_deg,
        ckm.delta_deg,
        ckm.jarlskog,
        ckm_r.d_theta12_deg,
        ckm_r.d_theta23_deg,
        ckm_r.d_theta13_deg,
        ckm_r.d_delta_deg,
        ckm_r.d_jarlskog,
        ckm_tex.theta12_deg,
        ckm_tex.theta23_deg,
        ckm_tex.theta13_deg,
        ckm_tex.delta_deg,
        ckm_tex.jarlskog,
        ckm_tex_r.d_theta12_deg,
        ckm_tex_r.d_theta23_deg,
        ckm_tex_r.d_theta13_deg,
        ckm_tex_r.d_delta_deg,
        ckm_tex_r.d_jarlskog,
        pmns.theta12_deg,
        pmns.theta23_deg,
        pmns.theta13_deg,
        pmns.delta_deg,
        pmns.jarlskog,
        pmns_r.d_theta12_deg,
        pmns_r.d_theta23_deg,
        pmns_r.d_theta13_deg,
        pmns_r.d_delta_deg,
        pmns_r.d_jarlskog,
        pmns_th23_corr_c,
        pmns_corr.theta12_deg,
        pmns_corr.theta23_deg,
        pmns_corr.theta13_deg,
        pmns_corr.delta_deg,
        pmns_corr.jarlskog,
        pmns_corr_r.d_theta12_deg,
        pmns_corr_r.d_theta23_deg,
        pmns_corr_r.d_theta13_deg,
        pmns_corr_r.d_delta_deg,
        pmns_corr_r.d_jarlskog,
        pmns_th23_corr_c,
        pmns_theta23_direct_abs,
        pmns_theta23_corr_abs,
        pmns_theta23_improvement_factor,
        pmns_theta23_improves_10x,
        pmns_tex.theta12_deg,
        pmns_tex.theta23_deg,
        pmns_tex.theta13_deg,
        pmns_tex.delta_deg,
        pmns_tex.jarlskog,
        pmns_tex_r.d_theta12_deg,
        pmns_tex_r.d_theta23_deg,
        pmns_tex_r.d_theta13_deg,
        pmns_tex_r.d_delta_deg,
        pmns_tex_r.d_jarlskog
    )
    .expect("write report json");

    let ts = unix_timestamp();
    let txt_snapshot_path = format!("{out_dir}/flavor_mix_report.{ts}.txt");
    let json_snapshot_path = format!("{out_dir}/flavor_mix_report.{ts}.json");
    let _ = fs::copy(&txt_path, &txt_snapshot_path);
    let _ = fs::copy(&json_path, &json_snapshot_path);

    append_trend_row(
        out_dir,
        ts,
        ckm.delta_deg,
        ckm_tex.delta_deg,
        pmns.delta_deg,
        pmns_tex.delta_deg,
        pmns_theta23_direct_abs,
        pmns_theta23_corr_abs,
        pmns_theta23_improvement_factor,
        pmns_theta23_improves_10x,
    );

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!("snapshotted {txt_snapshot_path}");
    println!("snapshotted {json_snapshot_path}");
    println!("appended {out_dir}/flavor_mix_trend.csv");
    println!(
        "CKM  θ12={:.3}° θ23={:.3}° θ13={:.3}° δ={:.3}° J={:.3e}",
        ckm.theta12_deg, ckm.theta23_deg, ckm.theta13_deg, ckm.delta_deg, ckm.jarlskog
    );
    println!(
        "PMNS θ12={:.3}° θ23={:.3}° θ13={:.3}° δ={:.3}° J={:.3e}",
        pmns.theta12_deg, pmns.theta23_deg, pmns.theta13_deg, pmns.delta_deg, pmns.jarlskog
    );
    println!(
        "PMNS(corr c={:.3}) θ12={:.3}° θ23={:.3}° θ13={:.3}° δ={:.3}° J={:.3e}",
        pmns_th23_corr_c,
        pmns_corr.theta12_deg,
        pmns_corr.theta23_deg,
        pmns_corr.theta13_deg,
        pmns_corr.delta_deg,
        pmns_corr.jarlskog
    );
    println!(
        "CKM(texture)  θ12={:.3}° θ23={:.3}° θ13={:.3}° δ={:.3}° J={:.3e}",
        ckm_tex.theta12_deg,
        ckm_tex.theta23_deg,
        ckm_tex.theta13_deg,
        ckm_tex.delta_deg,
        ckm_tex.jarlskog
    );
    println!(
        "PMNS(texture) θ12={:.3}° θ23={:.3}° θ13={:.3}° δ={:.3}° J={:.3e}",
        pmns_tex.theta12_deg,
        pmns_tex.theta23_deg,
        pmns_tex.theta13_deg,
        pmns_tex.delta_deg,
        pmns_tex.jarlskog
    );
}
