//! GRAND-355: CLASS full-shape likelihood contour scan.
//!
//! Produces fast contour CSVs for the dominant CMB-shape levers:
//! - A_s vs tau_reio
//! - omega_b h^2 vs omega_cdm h^2

use gutoe_physics::cmb_class::{
    compare_class_to_planck, read_class_tt_camb, read_planck_tt_csv, run_class,
    run_classy_fallback, write_class_ini, ClassRunInputs,
};
use gutoe_physics::constants::{
    lambda_cosmological_full_candidate, C, DARK_TO_VISIBLE_GEOMETRIC_RATIO,
};
use gutoe_physics::dark_matter_falsification::OMEGA_BARYON_OBS;
use gutoe_physics::{evaluate_inflation_gate, InflationWindows};
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

fn h0_from_lambda_and_omega_lambda(lambda: f64, omega_lambda: f64) -> f64 {
    let meter_per_mpc = 3.085_677_581_491_367e22;
    let h0_s_inv = C * (lambda / (3.0 * omega_lambda)).sqrt();
    h0_s_inv * meter_per_mpc / 1_000.0
}

fn derived_class_inputs() -> ClassRunInputs {
    let inflation = evaluate_inflation_gate(InflationWindows::default());
    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_cdm0;
    let omega_r0 = 9.0e-5;
    let omega_k0 = 0.0;
    let omega_lambda0 = 1.0 - omega_m0 - omega_r0 - omega_k0;
    let h0 = h0_from_lambda_and_omega_lambda(lambda_cosmological_full_candidate(), omega_lambda0);
    ClassRunInputs {
        h: h0 / 100.0,
        omega_b: omega_b0 * (h0 / 100.0).powi(2),
        omega_cdm: omega_cdm0 * (h0 / 100.0).powi(2),
        omega_k: omega_k0,
        omega_lambda: omega_lambda0,
        n_s: inflation.n_s,
        a_s: inflation.a_s,
        // Explicit assumption until tau derivation closes (GRAND-355 / GRAND-343).
        tau_reio: 0.054,
    }
}

fn parse_env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|x| x.parse::<f64>().ok())
        .unwrap_or(default)
}

fn parse_env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|x| x.parse::<usize>().ok())
        .unwrap_or(default)
}

fn linspace(n: usize, lo: f64, hi: f64) -> Vec<f64> {
    if n <= 1 {
        return vec![0.5 * (lo + hi)];
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        out.push(lo * (1.0 - t) + hi * t);
    }
    out
}

fn find_class_tt_output(run_dir: &Path) -> Result<PathBuf, String> {
    let mut candidates: Vec<PathBuf> = fs::read_dir(run_dir)
        .map_err(|e| format!("read CLASS run dir {:?}: {e}", run_dir))?
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|p| {
            p.extension()
                .and_then(|x| x.to_str())
                .is_some_and(|ext| ext.eq_ignore_ascii_case("dat"))
        })
        .filter(|p| {
            p.file_name()
                .and_then(|x| x.to_str())
                .is_some_and(|name| name.to_ascii_lowercase().contains("cl"))
        })
        .collect();
    if candidates.is_empty() {
        return Err(format!(
            "no CLASS .dat C_l files found in {:?} (expected *cl*.dat)",
            run_dir
        ));
    }
    candidates.sort_by_key(|p| {
        let name = p
            .file_name()
            .and_then(|x| x.to_str())
            .map(|x| x.to_ascii_lowercase())
            .unwrap_or_default();
        if name.contains("lensedcls") {
            0
        } else if name.ends_with("cl.dat") {
            1
        } else {
            2
        }
    });
    Ok(candidates[0].clone())
}

fn run_one_fit(
    base_run_dir: &Path,
    class_bin: &str,
    classy_python: &str,
    planck_path: &Path,
    inputs: ClassRunInputs,
    tag: &str,
) -> Result<(f64, f64), String> {
    let run_dir = base_run_dir.join(tag);
    let _ = fs::create_dir_all(&run_dir);
    let ini_path = run_dir.join("run.ini");
    let root = run_dir.join("g_");
    write_class_ini(&ini_path, &root.to_string_lossy(), 2_500, inputs)?;

    let class_tt_path = match run_class(class_bin, &ini_path) {
        Ok(_) => find_class_tt_output(&run_dir)?,
        Err(class_err) => {
            let fallback_path = run_dir.join("g_classy_cl.dat");
            run_classy_fallback(classy_python, &fallback_path, 2_500, inputs).map_err(|e| {
                format!(
                    "CLASS failed ({class_err}); classy fallback also failed ({e}) for tag={tag}"
                )
            })?;
            fallback_path
        }
    };

    let class_tt = read_class_tt_camb(&class_tt_path, 2, 2_500)?;
    let planck_tt_all = read_planck_tt_csv(planck_path)?;
    let planck_tt: Vec<_> = planck_tt_all
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();
    if planck_tt.is_empty() {
        return Err(format!(
            "no Planck points in supported ell range [2, 2500] for {:?}",
            planck_path
        ));
    }
    let fit = compare_class_to_planck(&class_tt, &planck_tt)?;
    Ok((fit.chi2, fit.reduced_chi2))
}

fn main() {
    let out_dir = std::env::var("GUTOE_CMB_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let planck_path = std::env::var("GUTOE_PLANCK_TT").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt".to_string()
    });
    let class_bin = std::env::var("GUTOE_CLASS_BIN").unwrap_or_else(|_| "class".to_string());
    let classy_python =
        std::env::var("GUTOE_CLASSY_PYTHON").unwrap_or_else(|_| "python3".to_string());

    let as_lo = parse_env_f64("GUTOE_AS_FACTOR_MIN", 0.94);
    let as_hi = parse_env_f64("GUTOE_AS_FACTOR_MAX", 1.02);
    let tau_lo = parse_env_f64("GUTOE_TAU_MIN", 0.045);
    let tau_hi = parse_env_f64("GUTOE_TAU_MAX", 0.070);
    let ob_lo = parse_env_f64("GUTOE_OB_FACTOR_MIN", 0.97);
    let ob_hi = parse_env_f64("GUTOE_OB_FACTOR_MAX", 1.03);
    let oc_lo = parse_env_f64("GUTOE_OC_FACTOR_MIN", 0.96);
    let oc_hi = parse_env_f64("GUTOE_OC_FACTOR_MAX", 1.04);
    let n_as_tau = parse_env_usize("GUTOE_AS_TAU_STEPS", 9);
    let n_ob_oc = parse_env_usize("GUTOE_OB_OC_STEPS", 9);

    let base = derived_class_inputs();
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    let run_root = std::env::temp_dir().join(format!("gutoe_cmb_scan_{stamp}"));
    let _ = fs::create_dir_all(&run_root);

    let as_vals = linspace(n_as_tau, as_lo, as_hi);
    let tau_vals = linspace(n_as_tau, tau_lo, tau_hi);
    let ob_vals = linspace(n_ob_oc, ob_lo, ob_hi);
    let oc_vals = linspace(n_ob_oc, oc_lo, oc_hi);

    let as_tau_csv = format!("{out_dir}/cmb_likelihood_as_tau.csv");
    let ob_oc_csv = format!("{out_dir}/cmb_likelihood_ob_oc.csv");
    let summary_json = format!("{out_dir}/cmb_likelihood_scan.json");
    let mut as_tau_f = File::create(&as_tau_csv).expect("create as_tau csv");
    let mut ob_oc_f = File::create(&ob_oc_csv).expect("create ob_oc csv");
    writeln!(as_tau_f, "as_factor,tau_reio,chi2,reduced_chi2").expect("write");
    writeln!(ob_oc_f, "ob_factor,oc_factor,chi2,reduced_chi2").expect("write");

    let mut best_as_tau = (f64::INFINITY, base.a_s, base.tau_reio, f64::INFINITY);
    for af in &as_vals {
        for tau in &tau_vals {
            let mut i = base;
            i.a_s = base.a_s * af;
            i.tau_reio = *tau;
            let tag = format!("as_{:.5}_tau_{:.5}", af, tau).replace('.', "p");
            let (chi2, red) = run_one_fit(
                &run_root,
                &class_bin,
                &classy_python,
                Path::new(&planck_path),
                i,
                &tag,
            )
            .expect("run as/tau point");
            if chi2 < best_as_tau.0 {
                best_as_tau = (chi2, i.a_s, i.tau_reio, red);
            }
            writeln!(as_tau_f, "{:.8},{:.8},{:.10},{:.10}", af, tau, chi2, red).expect("write");
        }
    }

    let mut best_ob_oc = (f64::INFINITY, base.omega_b, base.omega_cdm, f64::INFINITY);
    for obf in &ob_vals {
        for ocf in &oc_vals {
            let mut i = base;
            i.omega_b = base.omega_b * obf;
            i.omega_cdm = base.omega_cdm * ocf;
            let tag = format!("ob_{:.5}_oc_{:.5}", obf, ocf).replace('.', "p");
            let (chi2, red) = run_one_fit(
                &run_root,
                &class_bin,
                &classy_python,
                Path::new(&planck_path),
                i,
                &tag,
            )
            .expect("run ob/oc point");
            if chi2 < best_ob_oc.0 {
                best_ob_oc = (chi2, i.omega_b, i.omega_cdm, red);
            }
            writeln!(ob_oc_f, "{:.8},{:.8},{:.10},{:.10}", obf, ocf, chi2, red).expect("write");
        }
    }

    let mut j = File::create(&summary_json).expect("create summary json");
    writeln!(
        j,
        "{{\n  \"inputs\": {{\"class_bin\":\"{}\", \"classy_python\":\"{}\", \"planck_tt\":\"{}\"}},\n  \"base\": {{\"h\": {:.12}, \"omega_b\": {:.12}, \"omega_cdm\": {:.12}, \"n_s\": {:.12}, \"a_s\": {:.12e}, \"tau_reio\": {:.12}}},\n  \"scan\": {{\"as_factor_range\":[{:.6},{:.6}], \"tau_range\":[{:.6},{:.6}], \"ob_factor_range\":[{:.6},{:.6}], \"oc_factor_range\":[{:.6},{:.6}], \"as_tau_steps\": {}, \"ob_oc_steps\": {}}},\n  \"best_as_tau\": {{\"chi2\": {:.12}, \"reduced_chi2\": {:.12}, \"a_s\": {:.12e}, \"tau_reio\": {:.12}}},\n  \"best_ob_oc\": {{\"chi2\": {:.12}, \"reduced_chi2\": {:.12}, \"omega_b\": {:.12}, \"omega_cdm\": {:.12}}},\n  \"artifacts\": {{\"as_tau_csv\":\"{}\", \"ob_oc_csv\":\"{}\"}}\n}}",
        class_bin,
        classy_python,
        planck_path,
        base.h,
        base.omega_b,
        base.omega_cdm,
        base.n_s,
        base.a_s,
        base.tau_reio,
        as_lo,
        as_hi,
        tau_lo,
        tau_hi,
        ob_lo,
        ob_hi,
        oc_lo,
        oc_hi,
        n_as_tau,
        n_ob_oc,
        best_as_tau.0,
        best_as_tau.3,
        best_as_tau.1,
        best_as_tau.2,
        best_ob_oc.0,
        best_ob_oc.3,
        best_ob_oc.1,
        best_ob_oc.2,
        as_tau_csv,
        ob_oc_csv
    )
    .expect("write summary json");

    println!("wrote {}", as_tau_csv);
    println!("wrote {}", ob_oc_csv);
    println!("wrote {}", summary_json);
    println!(
        "best A_s/tau: chi2={:.3} red={:.3} A_s={:.6e} tau={:.5}",
        best_as_tau.0, best_as_tau.3, best_as_tau.1, best_as_tau.2
    );
    println!(
        "best ob/oc: chi2={:.3} red={:.3} omega_b={:.8} omega_cdm={:.8}",
        best_ob_oc.0, best_ob_oc.3, best_ob_oc.1, best_ob_oc.2
    );
}
