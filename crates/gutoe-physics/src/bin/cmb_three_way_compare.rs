//! GRAND-355: Three-way CMB consistency check.
//!
//! Compares:
//! 1) CLASS prediction vs Planck binned TT
//! 2) CLASS prediction vs Planck full (unbinned) TT
//! 3) Planck rebinned-from-full vs Planck published binned TT

use gutoe_physics::cmb_class::{
    compare_class_to_planck, read_class_tt_camb, read_planck_tt_csv, run_class,
    run_classy_fallback, write_class_ini, ClassRunInputs, PlanckTtPoint,
};
use gutoe_physics::cmb_reionization::derive_tau_reio;
use gutoe_physics::constants::{
    lambda_cosmological_full_candidate, C, DARK_TO_VISIBLE_GEOMETRIC_RATIO,
};
use gutoe_physics::dark_matter_falsification::OMEGA_BARYON_OBS;
use gutoe_physics::microphysics::MicrophysicsAssumptions;
use gutoe_physics::{evaluate_bbn_gate, evaluate_inflation_gate, BbnWindows, InflationWindows};
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

fn derived_tau_from_structure(inputs: ClassRunInputs) -> Result<(f64, f64), String> {
    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_cdm0;
    let bbn = evaluate_bbn_gate(BbnWindows::default());
    let micro = MicrophysicsAssumptions {
        h0_km_s_mpc: inputs.h * 100.0,
        omega_b0,
        omega_m0,
        omega_r0: 9.0e-5,
        omega_k0: inputs.omega_k,
        omega_lambda0: inputs.omega_lambda,
        eta10: bbn.eta10,
    };
    let reion = derive_tau_reio(micro, bbn.eta10)?;
    Ok((reion.tau_reio, reion.z_reion_structural))
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

fn as_chi2(
    points_a: &[PlanckTtPoint],
    points_b: &[PlanckTtPoint],
) -> Result<(usize, f64, f64), String> {
    let mut n = 0usize;
    let mut chi2 = 0.0_f64;
    let mut ia = 0usize;
    let mut ib = 0usize;
    while ia < points_a.len() && ib < points_b.len() {
        let a = points_a[ia];
        let b = points_b[ib];
        if a.ell == b.ell {
            let sigma2 = a.sigma_uk2 * a.sigma_uk2 + b.sigma_uk2 * b.sigma_uk2;
            if sigma2 > 0.0 {
                let pull = (a.d_ell_tt_uk2 - b.d_ell_tt_uk2) / sigma2.sqrt();
                chi2 += pull * pull;
                n += 1;
            }
            ia += 1;
            ib += 1;
        } else if a.ell < b.ell {
            ia += 1;
        } else {
            ib += 1;
        }
    }
    if n == 0 {
        return Err("no overlapping multipoles for chi2 comparison".to_string());
    }
    let red = chi2 / ((n as i64 - 1).max(1) as f64);
    Ok((n, chi2, red))
}

fn rebin_full_to_binned_centers(
    full: &[PlanckTtPoint],
    binned: &[PlanckTtPoint],
) -> Result<Vec<PlanckTtPoint>, String> {
    if binned.len() < 3 {
        return Err("need at least 3 binned points for stable bin windows".to_string());
    }
    let mut out = Vec::with_capacity(binned.len());
    for i in 0..binned.len() {
        let c = binned[i].ell as f64;
        let lo = if i == 0 {
            c - 0.5 * (binned[i + 1].ell as f64 - c)
        } else {
            0.5 * (binned[i - 1].ell as f64 + c)
        };
        let hi = if i + 1 == binned.len() {
            c + 0.5 * (c - binned[i - 1].ell as f64)
        } else {
            0.5 * (c + binned[i + 1].ell as f64)
        };

        let mut sum_w = 0.0_f64;
        let mut sum_wx = 0.0_f64;
        for p in full {
            let e = p.ell as f64;
            if e >= lo && e < hi && p.sigma_uk2 > 0.0 {
                let w = 1.0 / (p.sigma_uk2 * p.sigma_uk2);
                sum_w += w;
                sum_wx += w * p.d_ell_tt_uk2;
            }
        }
        if sum_w <= 0.0 {
            // If no full points fall in this window, skip this bin.
            continue;
        }
        out.push(PlanckTtPoint {
            ell: binned[i].ell,
            d_ell_tt_uk2: sum_wx / sum_w,
            sigma_uk2: (1.0 / sum_w).sqrt(),
        });
    }
    if out.is_empty() {
        return Err("rebinning produced no points".to_string());
    }
    Ok(out)
}

fn main() {
    let out_dir = std::env::var("GUTOE_CMB_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);

    let planck_binned_path = std::env::var("GUTOE_PLANCK_TT_BINNED").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt".to_string()
    });
    let planck_full_path = std::env::var("GUTOE_PLANCK_TT_FULL").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-full_R3.01.txt".to_string()
    });
    let class_bin = std::env::var("GUTOE_CLASS_BIN").unwrap_or_else(|_| "class".to_string());
    let classy_python =
        std::env::var("GUTOE_CLASSY_PYTHON").unwrap_or_else(|_| "python3".to_string());
    let tau_env = std::env::var("GUTOE_TAU_REIO")
        .ok()
        .and_then(|x| x.parse::<f64>().ok());

    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    let run_dir = std::env::temp_dir().join(format!("gutoe_cmb_threeway_{stamp}"));
    let _ = fs::create_dir_all(&run_dir);
    let ini_path = run_dir.join("run.ini");
    let root = run_dir.join("g_");

    let mut inputs = derived_class_inputs(0.054);
    let (tau_reio, z_reion_structural_opt) = if let Some(tau) = tau_env {
        (tau, None)
    } else {
        match derived_tau_from_structure(inputs) {
            Ok((tau, z_reion)) => (tau, Some(z_reion)),
            Err(e) => {
                eprintln!("failed to derive tau_reio structurally: {e}");
                std::process::exit(2);
            }
        }
    };
    inputs.tau_reio = tau_reio;
    write_class_ini(&ini_path, &root.to_string_lossy(), 2_500, inputs).expect("write ini");
    let (backend_used, class_tt_path) = match run_class(&class_bin, &ini_path) {
        Ok(_) => (
            "class",
            find_class_tt_output(&run_dir).expect("find class output"),
        ),
        Err(class_err) => {
            let fallback_path = run_dir.join("g_classy_cl.dat");
            run_classy_fallback(&classy_python, &fallback_path, 2_500, inputs).unwrap_or_else(
                |e| panic!("CLASS failed ({class_err}); classy fallback failed ({e})"),
            );
            ("classy", fallback_path)
        }
    };

    let class_tt = read_class_tt_camb(&class_tt_path, 2, 2_500).expect("parse class");
    let planck_binned =
        read_planck_tt_csv(Path::new(&planck_binned_path)).expect("parse planck binned");
    let planck_full = read_planck_tt_csv(Path::new(&planck_full_path)).expect("parse planck full");
    let planck_full_cut: Vec<_> = planck_full
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();

    let fit_pred_binned = compare_class_to_planck(&class_tt, &planck_binned).expect("fit binned");
    let fit_pred_full = compare_class_to_planck(&class_tt, &planck_full_cut).expect("fit full");

    let rebinned_from_full =
        rebin_full_to_binned_centers(&planck_full_cut, &planck_binned).expect("rebin full");
    let (n_bin_vs_rebinned, chi2_bin_vs_rebinned, red_bin_vs_rebinned) =
        as_chi2(&planck_binned, &rebinned_from_full).expect("binned vs rebinned");

    let txt_path = format!("{out_dir}/cmb_three_way_report.txt");
    let json_path = format!("{out_dir}/cmb_three_way_report.json");
    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[inputs]").expect("write");
    writeln!(txt, "backend = {}", backend_used).expect("write");
    writeln!(txt, "class_bin = {}", class_bin).expect("write");
    writeln!(txt, "classy_python = {}", classy_python).expect("write");
    writeln!(txt, "tau_reio = {:.6}", tau_reio).expect("write");
    writeln!(txt, "planck_binned = {}", planck_binned_path).expect("write");
    writeln!(txt, "planck_full = {}", planck_full_path).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[class_vs_planck_binned]").expect("write");
    writeln!(txt, "n = {}", fit_pred_binned.n_points).expect("write");
    writeln!(txt, "chi2 = {:.9}", fit_pred_binned.chi2).expect("write");
    writeln!(txt, "reduced_chi2 = {:.9}", fit_pred_binned.reduced_chi2).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[class_vs_planck_full]").expect("write");
    writeln!(txt, "n = {}", fit_pred_full.n_points).expect("write");
    writeln!(txt, "chi2 = {:.9}", fit_pred_full.chi2).expect("write");
    writeln!(txt, "reduced_chi2 = {:.9}", fit_pred_full.reduced_chi2).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[planck_binned_vs_planck_rebinned_from_full]").expect("write");
    writeln!(txt, "n = {}", n_bin_vs_rebinned).expect("write");
    writeln!(txt, "chi2 = {:.9}", chi2_bin_vs_rebinned).expect("write");
    writeln!(txt, "reduced_chi2 = {:.9}", red_bin_vs_rebinned).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"inputs\": {{\"backend\":\"{}\", \"class_bin\":\"{}\", \"classy_python\":\"{}\", \"tau_reio\": {:.6}, \"z_reion_structural\": {}, \"planck_binned\":\"{}\", \"planck_full\":\"{}\"}},\n  \"class_vs_planck_binned\": {{\"n\": {}, \"chi2\": {:.12}, \"reduced_chi2\": {:.12}}},\n  \"class_vs_planck_full\": {{\"n\": {}, \"chi2\": {:.12}, \"reduced_chi2\": {:.12}}},\n  \"planck_binned_vs_planck_rebinned_from_full\": {{\"n\": {}, \"chi2\": {:.12}, \"reduced_chi2\": {:.12}}}\n}}",
        backend_used,
        class_bin,
        classy_python,
        tau_reio,
        z_reion_structural_opt
            .map(|v| format!("{v:.12}"))
            .unwrap_or_else(|| "null".to_string()),
        planck_binned_path,
        planck_full_path,
        fit_pred_binned.n_points,
        fit_pred_binned.chi2,
        fit_pred_binned.reduced_chi2,
        fit_pred_full.n_points,
        fit_pred_full.chi2,
        fit_pred_full.reduced_chi2,
        n_bin_vs_rebinned,
        chi2_bin_vs_rebinned,
        red_bin_vs_rebinned
    )
    .expect("write json");

    println!("wrote {}", txt_path);
    println!("wrote {}", json_path);
    println!(
        "3-way: pred-vs-binned chi2={:.1} (n={}), pred-vs-full chi2={:.1} (n={}), binned-vs-rebinned chi2={:.1} (n={})",
        fit_pred_binned.chi2,
        fit_pred_binned.n_points,
        fit_pred_full.chi2,
        fit_pred_full.n_points,
        chi2_bin_vs_rebinned,
        n_bin_vs_rebinned
    );
    if let Some(zr) = z_reion_structural_opt {
        println!(
            "tau_reio derived structurally from reionization timing: tau={:.6}, z_reion={:.3}",
            tau_reio, zr
        );
    }
}
