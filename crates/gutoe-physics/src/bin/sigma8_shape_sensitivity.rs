use gutoe_physics::cmb_reionization::derive_tau_reio;
use gutoe_physics::constants::{lambda_cosmological_full_candidate, C, DARK_TO_VISIBLE_GEOMETRIC_RATIO};
use gutoe_physics::dark_matter_falsification::OMEGA_BARYON_OBS;
use gutoe_physics::microphysics::MicrophysicsAssumptions;
use gutoe_physics::{evaluate_bbn_gate, evaluate_inflation_gate, BbnWindows, InflationWindows};
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Clone, Copy, Debug)]
struct Cosmo {
    h: f64,
    omega_b: f64,
    omega_cdm: f64,
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

fn derived_cosmo() -> Result<Cosmo, String> {
    let inflation = evaluate_inflation_gate(InflationWindows::default());
    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_cdm0;
    let omega_r0 = 9.0e-5;
    let omega_k0 = 0.0;
    let omega_lambda0 = 1.0 - omega_m0 - omega_r0 - omega_k0;
    let h0 = h0_from_lambda_and_omega_lambda(lambda_cosmological_full_candidate(), omega_lambda0);
    let h = h0 / 100.0;
    let omega_b = omega_b0 * h * h;
    let omega_cdm = omega_cdm0 * h * h;

    let bbn = evaluate_bbn_gate(BbnWindows::default());
    let micro = MicrophysicsAssumptions {
        h0_km_s_mpc: h0,
        omega_b0,
        omega_m0,
        omega_r0: 9.0e-5,
        omega_k0,
        omega_lambda0,
        eta10: bbn.eta10,
    };
    let tau = derive_tau_reio(micro, bbn.eta10)
        .map_err(|e| format!("derive_tau_reio failed: {e}"))?
        .tau_reio;

    Ok(Cosmo {
        h,
        omega_b,
        omega_cdm,
        omega_k: omega_k0,
        n_s: inflation.n_s,
        a_s: inflation.a_s,
        tau_reio: tau,
    })
}

fn run_class_sigma8(class_bin: &str, out_dir: &Path, tag: &str, c: Cosmo) -> Result<f64, String> {
    let run_dir = out_dir.join(tag);
    fs::create_dir_all(&run_dir).map_err(|e| format!("mkdir {:?}: {e}", run_dir))?;
    let ini = run_dir.join("in.ini");
    let root = run_dir.join("g_");
    let ini_txt = format!(
        "h = {h}\nomega_b = {ob}\nomega_cdm = {oc}\nOmega_k = {ok}\nA_s = {as_}\nn_s = {ns}\ntau_reio = {tau}\noutput = mPk\nP_k_max_h/Mpc = 50\nz_pk = 0\nroot = {root}\n",
        h = c.h,
        ob = c.omega_b,
        oc = c.omega_cdm,
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

fn main() {
    let class_bin = std::env::var("GUTOE_CLASS_BIN")
        .unwrap_or_else(|_| "/tmp/class_public/class".to_string());
    let out_dir = std::env::var("GUTOE_SIGMA8_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/sigma8_shape_sensitivity".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let base = match derived_cosmo() {
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

    let d_oc = base.omega_cdm * 0.01;
    let d_h = base.h * 0.005;
    let d_ns = 0.002;

    let mut c = base;
    c.omega_cdm = base.omega_cdm + d_oc;
    let s_oc_p = run_class_sigma8(&class_bin, &out, "oc_plus", c).unwrap_or(f64::NAN);
    c.omega_cdm = base.omega_cdm - d_oc;
    let s_oc_m = run_class_sigma8(&class_bin, &out, "oc_minus", c).unwrap_or(f64::NAN);
    let ds_doc = (s_oc_p - s_oc_m) / (2.0 * d_oc);

    c = base;
    c.h = base.h + d_h;
    let s_h_p = run_class_sigma8(&class_bin, &out, "h_plus", c).unwrap_or(f64::NAN);
    c.h = base.h - d_h;
    let s_h_m = run_class_sigma8(&class_bin, &out, "h_minus", c).unwrap_or(f64::NAN);
    let ds_dh = (s_h_p - s_h_m) / (2.0 * d_h);

    c = base;
    c.n_s = base.n_s + d_ns;
    let s_ns_p = run_class_sigma8(&class_bin, &out, "ns_plus", c).unwrap_or(f64::NAN);
    c.n_s = base.n_s - d_ns;
    let s_ns_m = run_class_sigma8(&class_bin, &out, "ns_minus", c).unwrap_or(f64::NAN);
    let ds_dns = (s_ns_p - s_ns_m) / (2.0 * d_ns);

    let report = out.join("sigma8_shape_sensitivity_report.json");
    let mut f = File::create(&report).expect("create report");
    writeln!(
        f,
        "{{\n  \"base\": {{\"h\": {:.12}, \"omega_b\": {:.12}, \"omega_cdm\": {:.12}, \"n_s\": {:.12}, \"A_s\": {:.12e}, \"tau_reio\": {:.12}, \"sigma8\": {:.12}}},\n  \"steps\": {{\"d_omega_cdm\": {:.12e}, \"d_h\": {:.12e}, \"d_n_s\": {:.12e}}},\n  \"finite_difference\": {{\"d_sigma8_d_omega_cdm\": {:.12}, \"d_sigma8_d_h\": {:.12}, \"d_sigma8_d_n_s\": {:.12}}},\n  \"samples\": {{\"sigma8_oc_plus\": {:.12}, \"sigma8_oc_minus\": {:.12}, \"sigma8_h_plus\": {:.12}, \"sigma8_h_minus\": {:.12}, \"sigma8_ns_plus\": {:.12}, \"sigma8_ns_minus\": {:.12}}}\n}}",
        base.h,
        base.omega_b,
        base.omega_cdm,
        base.n_s,
        base.a_s,
        base.tau_reio,
        sigma0,
        d_oc,
        d_h,
        d_ns,
        ds_doc,
        ds_dh,
        ds_dns,
        s_oc_p,
        s_oc_m,
        s_h_p,
        s_h_m,
        s_ns_p,
        s_ns_m
    )
    .expect("write");
    println!("wrote {}", report.display());
    println!(
        "sigma8={:.6} | dσ8/d(ω_cdm)={:.3}, dσ8/dh={:.3}, dσ8/dn_s={:.3}",
        sigma0, ds_doc, ds_dh, ds_dns
    );
}
