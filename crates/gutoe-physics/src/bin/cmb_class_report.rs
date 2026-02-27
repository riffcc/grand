//! GRAND-355: Full-shape CMB TT check via CLASS + Planck binned data.

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

fn derived_class_inputs(tau_reio: f64) -> ClassRunInputs {
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
        tau_reio,
    }
}

fn parse_env_f64(name: &str) -> Option<f64> {
    let raw = std::env::var(name).ok()?;
    raw.parse::<f64>().ok()
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

fn main() {
    let out_dir = std::env::var("GUTOE_CMB_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);

    let planck_path = std::env::var("GUTOE_PLANCK_TT").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt".to_string()
    });
    let class_bin = std::env::var("GUTOE_CLASS_BIN").unwrap_or_else(|_| "class".to_string());
    let classy_python =
        std::env::var("GUTOE_CLASSY_PYTHON").unwrap_or_else(|_| "python3".to_string());
    let tau_from_env = parse_env_f64("GUTOE_TAU_REIO");
    let tau_reio = tau_from_env.unwrap_or(0.054);
    let tau_assumption = tau_from_env.is_none();

    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    let run_dir = std::env::temp_dir().join(format!("gutoe_class_run_{stamp}"));
    let _ = fs::create_dir_all(&run_dir);
    let ini_path = run_dir.join("gutoe_class.ini");
    let root = run_dir.join("gutoe_");
    let root_str = root.to_string_lossy().to_string();

    let inputs = derived_class_inputs(tau_reio);
    if let Err(e) = write_class_ini(&ini_path, &root_str, 2_500, inputs) {
        eprintln!("failed to write CLASS ini: {e}");
        std::process::exit(2);
    }
    let (backend_used, class_tt_path) = match run_class(&class_bin, &ini_path) {
        Ok(_) => match find_class_tt_output(&run_dir) {
            Ok(p) => ("class", p),
            Err(e) => {
                eprintln!("failed to locate CLASS TT output: {e}");
                std::process::exit(2);
            }
        },
        Err(class_err) => {
            let fallback_path = run_dir.join("gutoe_classy_cl.dat");
            match run_classy_fallback(&classy_python, &fallback_path, 2_500, inputs) {
                Ok(_) => {
                    eprintln!(
                        "CLASS binary failed ({}); used classy fallback via {}",
                        class_err, classy_python
                    );
                    ("classy", fallback_path)
                }
                Err(classy_err) => {
                    eprintln!("failed to execute CLASS: {class_err}");
                    eprintln!("failed to execute classy fallback: {classy_err}");
                    eprintln!(
                        "hint: install CLASS and/or python classy.\n  class bin: '{}'\n  python: '{}'",
                        class_bin, classy_python
                    );
                    std::process::exit(2);
                }
            }
        }
    };
    let class_tt = match read_class_tt_camb(&class_tt_path, 2, 2_500) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("failed to parse CLASS TT output {:?}: {e}", class_tt_path);
            std::process::exit(2);
        }
    };
    let planck_tt = match read_planck_tt_csv(Path::new(&planck_path)) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("failed to parse Planck TT dataset {}: {e}", planck_path);
            eprintln!("expected CSV columns: ell,d_ell_tt_uk2,sigma_uk2");
            std::process::exit(2);
        }
    };
    let fit = match compare_class_to_planck(&class_tt, &planck_tt) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("failed to compare CLASS vs Planck: {e}");
            std::process::exit(2);
        }
    };

    let txt_path = format!("{out_dir}/cmb_class_report.txt");
    let json_path = format!("{out_dir}/cmb_class_report.json");
    let residual_csv_path = format!("{out_dir}/cmb_class_residuals.csv");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[inputs]").expect("write");
    writeln!(txt, "class_bin = {}", class_bin).expect("write");
    writeln!(txt, "classy_python = {}", classy_python).expect("write");
    writeln!(txt, "backend_used = {}", backend_used).expect("write");
    writeln!(txt, "planck_tt_csv = {}", planck_path).expect("write");
    writeln!(txt, "h = {:.12}", inputs.h).expect("write");
    writeln!(txt, "omega_b = {:.12}", inputs.omega_b).expect("write");
    writeln!(txt, "omega_cdm = {:.12}", inputs.omega_cdm).expect("write");
    writeln!(txt, "omega_k = {:.12}", inputs.omega_k).expect("write");
    writeln!(txt, "omega_lambda = {:.12}", inputs.omega_lambda).expect("write");
    writeln!(txt, "n_s = {:.12}", inputs.n_s).expect("write");
    writeln!(txt, "A_s = {:.12e}", inputs.a_s).expect("write");
    writeln!(txt, "tau_reio = {:.12}", inputs.tau_reio).expect("write");
    writeln!(txt, "tau_reio_is_assumption = {}", tau_assumption).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[fit]").expect("write");
    writeln!(txt, "n_points = {}", fit.n_points).expect("write");
    writeln!(txt, "chi2 = {:.9}", fit.chi2).expect("write");
    writeln!(txt, "reduced_chi2 = {:.9}", fit.reduced_chi2).expect("write");
    writeln!(txt, "mean_abs_pull = {:.9}", fit.mean_abs_pull).expect("write");
    writeln!(txt, "max_abs_pull = {:.9}", fit.max_abs_pull).expect("write");
    writeln!(txt, "rms_residual_uk2 = {:.9}", fit.rms_residual_uk2).expect("write");

    let mut residual_csv = File::create(&residual_csv_path).expect("create residual csv");
    writeln!(
        residual_csv,
        "ell,observed_uk2,predicted_uk2,sigma_uk2,pull"
    )
    .expect("write");
    for r in &fit.residuals {
        writeln!(
            residual_csv,
            "{},{:.12},{:.12},{:.12},{:.12}",
            r.ell, r.observed_uk2, r.predicted_uk2, r.sigma_uk2, r.pull
        )
        .expect("write");
    }

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"inputs\": {{\"class_bin\":\"{}\", \"planck_tt_csv\":\"{}\", \"h\": {:.12}, \"omega_b\": {:.12}, \"omega_cdm\": {:.12}, \"omega_k\": {:.12}, \"omega_lambda\": {:.12}, \"n_s\": {:.12}, \"a_s\": {:.12e}, \"tau_reio\": {:.12}, \"tau_reio_is_assumption\": {}}},\n  \"fit\": {{\"n_points\": {}, \"chi2\": {:.12}, \"reduced_chi2\": {:.12}, \"mean_abs_pull\": {:.12}, \"max_abs_pull\": {:.12}, \"rms_residual_uk2\": {:.12}}},\n  \"artifacts\": {{\"class_tt_path\":\"{}\", \"residual_csv\":\"{}\"}}\n}}",
        class_bin,
        planck_path,
        inputs.h,
        inputs.omega_b,
        inputs.omega_cdm,
        inputs.omega_k,
        inputs.omega_lambda,
        inputs.n_s,
        inputs.a_s,
        inputs.tau_reio,
        tau_assumption,
        fit.n_points,
        fit.chi2,
        fit.reduced_chi2,
        fit.mean_abs_pull,
        fit.max_abs_pull,
        fit.rms_residual_uk2,
        class_tt_path.to_string_lossy(),
        residual_csv_path
    )
    .expect("write json");

    println!("wrote {}", txt_path);
    println!("wrote {}", json_path);
    println!("wrote {}", residual_csv_path);
    println!(
        "CMB TT full-shape fit ({backend_used}): n={} chi2={:.3} red_chi2={:.3} mean|pull|={:.3} max|pull|={:.3}",
        fit.n_points,
        fit.chi2,
        fit.reduced_chi2,
        fit.mean_abs_pull,
        fit.max_abs_pull
    );
    if tau_assumption {
        println!(
            "NOTE: tau_reio not yet derived in-framework; using explicit assumption tau_reio={:.6}",
            inputs.tau_reio
        );
    }
}
