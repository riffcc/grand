use gutoe_physics::bbn::{evaluate_bbn_gate, BbnWindows};
use gutoe_physics::cmb_reionization::derive_tau_reio;
use gutoe_physics::constants::{lambda_cosmological_full_candidate, C, DARK_TO_VISIBLE_GEOMETRIC_RATIO};
use gutoe_physics::dark_matter_falsification::OMEGA_BARYON_OBS;
use gutoe_physics::microphysics::MicrophysicsAssumptions;
use gutoe_physics::{evaluate_inflation_gate, InflationWindows};
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

const SIGMA8_TARGET_PLANCK: f64 = 0.811;
const OMEGA_CDM_H2_TARGET_PLANCK: f64 = 0.1200;
const ETA10_TO_OMEGA_B_H2: f64 = 273.9;

#[derive(Clone, Copy, Debug)]
struct Cosmo {
    h: f64,
    omega_b_h2: f64,
    omega_cdm_h2: f64,
    omega_k: f64,
    n_s: f64,
    a_s: f64,
    tau_reio: f64,
}

fn h0_from_lambda_and_omega_lambda(lambda: f64, omega_lambda: f64) -> f64 {
    let meter_per_mpc = 3.085_677_581_491_367e22;
    let h0_s_inv = C * (lambda / (3.0 * omega_lambda)).sqrt();
    h0_s_inv * meter_per_mpc / 1_000.0
}

fn derived_cosmo() -> Result<(Cosmo, f64), String> {
    let inflation = evaluate_inflation_gate(InflationWindows::default());
    let bbn = evaluate_bbn_gate(BbnWindows::default());

    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_cdm0;
    let omega_r0 = 9.0e-5;
    let omega_k0 = 0.0;
    let omega_lambda0 = 1.0 - omega_m0 - omega_r0 - omega_k0;

    let h0 = h0_from_lambda_and_omega_lambda(lambda_cosmological_full_candidate(), omega_lambda0);
    let h = h0 / 100.0;
    let omega_b_h2 = omega_b0 * h * h;
    let omega_cdm_h2 = omega_cdm0 * h * h;

    let micro = MicrophysicsAssumptions {
        h0_km_s_mpc: h0,
        omega_b0,
        omega_m0,
        omega_r0,
        omega_k0,
        omega_lambda0,
        eta10: bbn.eta10,
    };
    let tau_reio = derive_tau_reio(micro, bbn.eta10)
        .map_err(|e| format!("derive_tau_reio failed: {e}"))?
        .tau_reio;

    Ok((
        Cosmo {
            h,
            omega_b_h2,
            omega_cdm_h2,
            omega_k: omega_k0,
            n_s: inflation.n_s,
            a_s: inflation.a_s,
            tau_reio,
        },
        bbn.eta10,
    ))
}

fn run_class_sigma8(class_bin: &str, out_dir: &Path, tag: &str, c: Cosmo) -> Result<f64, String> {
    let run_dir = out_dir.join(tag);
    fs::create_dir_all(&run_dir).map_err(|e| format!("mkdir {:?}: {e}", run_dir))?;
    let ini = run_dir.join("in.ini");
    let root = run_dir.join("g_");
    let ini_txt = format!(
        "h = {h}\nomega_b = {ob}\nomega_cdm = {oc}\nOmega_k = {ok}\nA_s = {as_}\nn_s = {ns}\ntau_reio = {tau}\noutput = mPk\nP_k_max_h/Mpc = 50\nz_pk = 0\nroot = {root}\n",
        h = c.h,
        ob = c.omega_b_h2,
        oc = c.omega_cdm_h2,
        ok = c.omega_k,
        as_ = c.a_s,
        ns = c.n_s,
        tau = c.tau_reio,
        root = root.to_string_lossy()
    );
    fs::write(&ini, ini_txt).map_err(|e| format!("write ini {:?}: {e}", ini))?;
    let status = Command::new(class_bin)
        .arg(&ini)
        .status()
        .map_err(|e| format!("run class '{}': {e}", class_bin))?;
    if !status.success() {
        return Err(format!("CLASS failed status {status}"));
    }
    let mut pk_files: Vec<PathBuf> = fs::read_dir(&run_dir)
        .map_err(|e| format!("read run dir: {e}"))?
        .filter_map(|e| e.ok().map(|x| x.path()))
        .filter(|p| {
            p.extension()
                .and_then(|x| x.to_str())
                .is_some_and(|ext| ext.eq_ignore_ascii_case("dat"))
                && p.file_name()
                    .and_then(|x| x.to_str())
                    .is_some_and(|n| n.to_ascii_lowercase().contains("pk"))
        })
        .collect();
    pk_files.sort();
    let pk = pk_files
        .first()
        .ok_or_else(|| "no pk dat file".to_string())?
        .clone();
    sigma8_from_pk(&pk, 8.0)
}

fn sigma8_from_pk(pk_path: &Path, r_hinv_mpc: f64) -> Result<f64, String> {
    let f = File::open(pk_path).map_err(|e| format!("open pk {:?}: {e}", pk_path))?;
    let mut pairs: Vec<(f64, f64)> = Vec::new();
    for line in BufReader::new(f).lines() {
        let s = line.map_err(|e| format!("read pk line: {e}"))?;
        let s = s.trim();
        if s.is_empty() || s.starts_with('#') {
            continue;
        }
        let sp: Vec<&str> = s.split_whitespace().collect();
        if sp.len() < 2 {
            continue;
        }
        let k: f64 = match sp[0].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let p: f64 = match sp[1].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if k > 0.0 && p > 0.0 {
            pairs.push((k, p));
        }
    }
    if pairs.len() < 16 {
        return Err("insufficient pk rows".to_string());
    }
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut acc = 0.0;
    for w in pairs.windows(2) {
        let (k0, p0) = w[0];
        let (k1, p1) = w[1];
        let f = |k: f64, p: f64| {
            let x = k * r_hinv_mpc;
            let w = if x.abs() < 1e-8 {
                1.0
            } else {
                3.0 * (x.sin() - x * x.cos()) / (x * x * x)
            };
            p * w * w * k * k
        };
        acc += 0.5 * (k1 - k0) * (f(k0, p0) + f(k1, p1));
    }
    Ok((acc / (2.0 * PI * PI)).max(0.0).sqrt())
}

fn elasticity(
    class_bin: &str,
    out: &Path,
    base: Cosmo,
    base_sigma8: f64,
    name: &str,
    frac_step: f64,
) -> Result<(f64, f64, f64, f64, f64), String> {
    let mut plus = base;
    let mut minus = base;
    match name {
        "h" => {
            plus.h *= 1.0 + frac_step;
            minus.h *= 1.0 - frac_step;
        }
        "n_s" => {
            plus.n_s *= 1.0 + frac_step;
            minus.n_s *= 1.0 - frac_step;
        }
        "omega_cdm_h2" => {
            plus.omega_cdm_h2 *= 1.0 + frac_step;
            minus.omega_cdm_h2 *= 1.0 - frac_step;
        }
        _ => return Err(format!("unknown parameter '{name}'")),
    }

    let s_plus = run_class_sigma8(class_bin, out, &format!("{name}_p_{frac_step:.4}"), plus)?;
    let s_minus = run_class_sigma8(class_bin, out, &format!("{name}_m_{frac_step:.4}"), minus)?;

    let dlns = (s_plus.ln() - s_minus.ln()) / ((1.0 + frac_step).ln() - (1.0 - frac_step).ln());
    let ds = (s_plus - s_minus) / (2.0 * frac_step);
    let dsdx = ds / match name {
        "h" => base.h,
        "n_s" => base.n_s,
        "omega_cdm_h2" => base.omega_cdm_h2,
        _ => 1.0,
    };
    let unit_slope = dsdx;
    let elasticity = dlns;
    let recon = elasticity * (base_sigma8 / match name {
        "h" => base.h,
        "n_s" => base.n_s,
        "omega_cdm_h2" => base.omega_cdm_h2,
        _ => 1.0,
    });
    Ok((s_plus, s_minus, unit_slope, elasticity, recon))
}

fn main() {
    let class_bin = std::env::var("GUTOE_CLASS_BIN")
        .unwrap_or_else(|_| "/tmp/class_public/class".to_string());
    let out_dir =
        std::env::var("GUTOE_SIGMA8_OUT").unwrap_or_else(|_| "/tmp/bh_renders/sigma8_torsion_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let (base, eta10) = match derived_cosmo() {
        Ok(v) => v,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    };
    let sigma0 = match run_class_sigma8(&class_bin, &out, "base", base) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    };

    let step_small = 0.005;
    let step_large = 0.01;
    let e_h_s = elasticity(&class_bin, &out, base, sigma0, "h", step_small).unwrap_or((f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN));
    let e_h_l = elasticity(&class_bin, &out, base, sigma0, "h", step_large).unwrap_or((f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN));
    let e_ns_s = elasticity(&class_bin, &out, base, sigma0, "n_s", step_small).unwrap_or((f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN));
    let e_ns_l = elasticity(&class_bin, &out, base, sigma0, "n_s", step_large).unwrap_or((f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN));
    let e_oc_s =
        elasticity(&class_bin, &out, base, sigma0, "omega_cdm_h2", step_small).unwrap_or((f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN));
    let e_oc_l =
        elasticity(&class_bin, &out, base, sigma0, "omega_cdm_h2", step_large).unwrap_or((f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN));

    let ratio_now = base.omega_cdm_h2 / base.omega_b_h2;
    let delta_for_omega = 60.0 - 11.0 * (OMEGA_CDM_H2_TARGET_PLANCK / base.omega_b_h2);

    let mut c_delta_omega = base;
    c_delta_omega.omega_cdm_h2 = OMEGA_CDM_H2_TARGET_PLANCK;
    let sigma_delta_omega =
        run_class_sigma8(&class_bin, &out, "delta_target_omega", c_delta_omega).unwrap_or(f64::NAN);

    let delta_for_sigma = 60.0 - 11.0 * ((base.omega_cdm_h2 / base.omega_b_h2) * (SIGMA8_TARGET_PLANCK / sigma0));
    let mut c_delta_sigma = base;
    c_delta_sigma.omega_cdm_h2 = base.omega_b_h2 * ((60.0 - delta_for_sigma) / 11.0);
    let sigma_delta_sigma =
        run_class_sigma8(&class_bin, &out, "delta_target_sigma_linearized", c_delta_sigma).unwrap_or(f64::NAN);

    let omega_b_h2_bbn = eta10 / ETA10_TO_OMEGA_B_H2;
    let omega_b0_bbn = omega_b_h2_bbn / (base.h * base.h);
    let omega_cdm_h2_bbn_anchor = omega_b_h2_bbn * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let mut c_bbn_anchor = base;
    c_bbn_anchor.omega_b_h2 = omega_b_h2_bbn;
    c_bbn_anchor.omega_cdm_h2 = omega_cdm_h2_bbn_anchor;
    let sigma_bbn_anchor =
        run_class_sigma8(&class_bin, &out, "bbn_anchor", c_bbn_anchor).unwrap_or(f64::NAN);

    let report = out.join("sigma8_torsion_probe_report.json");
    let mut f = File::create(&report).expect("create report");
    writeln!(f, "{{").expect("write");
    writeln!(
        f,
        "  \"base\": {{\"h\": {:.12}, \"omega_b_h2\": {:.12}, \"omega_cdm_h2\": {:.12}, \"ratio_omega_cdm_to_omega_b\": {:.12}, \"n_s\": {:.12}, \"A_s\": {:.12e}, \"tau_reio\": {:.12}, \"sigma8\": {:.12}}},",
        base.h, base.omega_b_h2, base.omega_cdm_h2, ratio_now, base.n_s, base.a_s, base.tau_reio, sigma0
    )
    .expect("write");
    writeln!(f, "  \"elasticity\": {{").expect("write");
    writeln!(
        f,
        "    \"h\": {{\"step_0p5pct\": {{\"d_sigma8_dh\": {:.12}, \"dln_sigma8_dln_h\": {:.12}}}, \"step_1pct\": {{\"d_sigma8_dh\": {:.12}, \"dln_sigma8_dln_h\": {:.12}}}}},",
        e_h_s.2, e_h_s.3, e_h_l.2, e_h_l.3
    )
    .expect("write");
    writeln!(
        f,
        "    \"n_s\": {{\"step_0p5pct\": {{\"d_sigma8_dns\": {:.12}, \"dln_sigma8_dln_ns\": {:.12}}}, \"step_1pct\": {{\"d_sigma8_dns\": {:.12}, \"dln_sigma8_dln_ns\": {:.12}}}}},",
        e_ns_s.2, e_ns_s.3, e_ns_l.2, e_ns_l.3
    )
    .expect("write");
    writeln!(
        f,
        "    \"omega_cdm_h2\": {{\"step_0p5pct\": {{\"d_sigma8_doc\": {:.12}, \"dln_sigma8_dln_oc\": {:.12}}}, \"step_1pct\": {{\"d_sigma8_doc\": {:.12}, \"dln_sigma8_dln_oc\": {:.12}}}}}",
        e_oc_s.2, e_oc_s.3, e_oc_l.2, e_oc_l.3
    )
    .expect("write");
    writeln!(f, "  }},").expect("write");
    writeln!(
        f,
        "  \"delta_scan\": {{\"delta_for_omega_cdm_h2_target\": {:.12}, \"sigma8_at_omega_cdm_target\": {:.12}, \"delta_for_sigma8_target_linearized\": {:.12}, \"sigma8_at_delta_sigma_linearized\": {:.12}}},",
        delta_for_omega, sigma_delta_omega, delta_for_sigma, sigma_delta_sigma
    )
    .expect("write");
    writeln!(
        f,
        "  \"bbn_anchor\": {{\"eta10\": {:.12}, \"eta10_to_omega_b_h2_factor\": {:.3}, \"omega_b_h2_from_eta10\": {:.12}, \"omega_b0_from_eta10\": {:.12}, \"omega_cdm_h2_with_ratio_60_over_11\": {:.12}, \"sigma8\": {:.12}}}",
        eta10, ETA10_TO_OMEGA_B_H2, omega_b_h2_bbn, omega_b0_bbn, omega_cdm_h2_bbn_anchor, sigma_bbn_anchor
    )
    .expect("write");
    writeln!(f, "}}").expect("write");

    println!("wrote {}", report.display());
    println!(
        "sigma8={:.6} | dlnσ8/dln(h)~{:.4} | dlnσ8/dln(ns)~{:.4} | dlnσ8/dln(ωcdm)~{:.4}",
        sigma0, e_h_l.3, e_ns_l.3, e_oc_l.3
    );
    println!(
        "delta_for_omega_target={:.4}, delta_for_sigma_target_linearized={:.4}, sigma8_bbn_anchor={:.6}",
        delta_for_omega, delta_for_sigma, sigma_bbn_anchor
    );
}
