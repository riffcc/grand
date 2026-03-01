//! Combined TT+TE+EE A_s–tau profile scan (diagonal-chi2 lane).
//! Purpose: minimum kill-test harness in current infrastructure.

use gutoe_physics::cmb_class::{
    compare_class_to_planck, read_class_dl_camb_column, read_planck_dl_csv, run_class,
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

#[derive(Debug, Clone, Copy)]
struct ChannelScore {
    chi2: f64,
    red: f64,
    n: usize,
}

#[derive(Debug, Clone, Copy)]
struct CombinedScore {
    tt: ChannelScore,
    te: ChannelScore,
    ee: ChannelScore,
    chi2_total: f64,
    red_total: f64,
}

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

fn find_class_cl_output(run_dir: &Path) -> Result<PathBuf, String> {
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
    tt_full: &Path,
    te_full: &Path,
    ee_full: &Path,
    inputs: ClassRunInputs,
    tag: &str,
) -> Result<CombinedScore, String> {
    let run_dir = base_run_dir.join(tag);
    let _ = fs::create_dir_all(&run_dir);
    let ini_path = run_dir.join("run.ini");
    let root = run_dir.join("g_");
    write_class_ini(&ini_path, &root.to_string_lossy(), 2_500, inputs)?;

    let class_cl_path = match run_class(class_bin, &ini_path) {
        Ok(_) => find_class_cl_output(&run_dir)?,
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

    let tt = read_class_dl_camb_column(&class_cl_path, 2, 2_500, 2)?;
    let ee = read_class_dl_camb_column(&class_cl_path, 2, 2_500, 3)?;
    let te = read_class_dl_camb_column(&class_cl_path, 2, 2_500, 5)?;

    let tt_obs: Vec<_> = read_planck_dl_csv(tt_full)?
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();
    let te_obs: Vec<_> = read_planck_dl_csv(te_full)?
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();
    let ee_obs: Vec<_> = read_planck_dl_csv(ee_full)?
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();

    let fit_tt = compare_class_to_planck(&tt, &tt_obs)?;
    let fit_te = compare_class_to_planck(&te, &te_obs)?;
    let fit_ee = compare_class_to_planck(&ee, &ee_obs)?;

    let tt_score = ChannelScore {
        chi2: fit_tt.chi2,
        red: fit_tt.reduced_chi2,
        n: fit_tt.n_points,
    };
    let te_score = ChannelScore {
        chi2: fit_te.chi2,
        red: fit_te.reduced_chi2,
        n: fit_te.n_points,
    };
    let ee_score = ChannelScore {
        chi2: fit_ee.chi2,
        red: fit_ee.reduced_chi2,
        n: fit_ee.n_points,
    };

    let chi2_total = tt_score.chi2 + te_score.chi2 + ee_score.chi2;
    let n_total = tt_score.n + te_score.n + ee_score.n;
    let ndof_total = (n_total as i64 - 1).max(1) as f64;
    let red_total = chi2_total / ndof_total;

    Ok(CombinedScore {
        tt: tt_score,
        te: te_score,
        ee: ee_score,
        chi2_total,
        red_total,
    })
}

fn main() {
    let out_dir = std::env::var("GUTOE_CMB_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);

    let class_bin = std::env::var("GUTOE_CLASS_BIN").unwrap_or_else(|_| "class".to_string());
    let classy_python =
        std::env::var("GUTOE_CLASSY_PYTHON").unwrap_or_else(|_| "python3".to_string());

    let tt_full = std::env::var("GUTOE_PLANCK_TT_FULL").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-full_R3.01.txt".to_string()
    });
    let te_full = std::env::var("GUTOE_PLANCK_TE_FULL").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TE-full_R3.01.txt".to_string()
    });
    let ee_full = std::env::var("GUTOE_PLANCK_EE_FULL").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-full_R3.01.txt".to_string()
    });

    let as_lo = parse_env_f64("GUTOE_AS_FACTOR_MIN", 0.90);
    let as_hi = parse_env_f64("GUTOE_AS_FACTOR_MAX", 1.10);
    let tau_lo = parse_env_f64("GUTOE_TAU_MIN", 0.040);
    let tau_hi = parse_env_f64("GUTOE_TAU_MAX", 0.080);
    let n_as_tau = parse_env_usize("GUTOE_AS_TAU_STEPS", 17);

    let base = derived_class_inputs();

    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    let run_root = std::env::temp_dir().join(format!("gutoe_cmb_ttee_scan_{stamp}"));
    let _ = fs::create_dir_all(&run_root);

    let as_vals = linspace(n_as_tau, as_lo, as_hi);
    let tau_vals = linspace(n_as_tau, tau_lo, tau_hi);

    let csv_path = format!("{out_dir}/cmb_ttee_as_tau_profile.csv");
    let json_path = format!("{out_dir}/cmb_ttee_as_tau_profile.json");
    let mut out = File::create(&csv_path).expect("create profile csv");
    writeln!(
        out,
        "as_factor,tau_reio,chi2_total,reduced_total,chi2_tt,chi2_te,chi2_ee,red_tt,red_te,red_ee"
    )
    .expect("write");

    let mut best = (
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        ChannelScore { chi2: 0.0, red: 0.0, n: 0 },
        ChannelScore { chi2: 0.0, red: 0.0, n: 0 },
        ChannelScore { chi2: 0.0, red: 0.0, n: 0 },
    );

    for af in &as_vals {
        for tau in &tau_vals {
            let mut i = base;
            i.a_s = base.a_s * af;
            i.tau_reio = *tau;
            let tag = format!("as_{:.5}_tau_{:.5}", af, tau).replace('.', "p");
            let s = run_one_fit(
                &run_root,
                &class_bin,
                &classy_python,
                Path::new(&tt_full),
                Path::new(&te_full),
                Path::new(&ee_full),
                i,
                &tag,
            )
            .expect("run ttee point");

            if s.chi2_total < best.0 {
                best = (s.chi2_total, s.red_total, i.a_s, i.tau_reio, s.tt, s.te, s.ee);
            }

            writeln!(
                out,
                "{:.8},{:.8},{:.10},{:.10},{:.10},{:.10},{:.10},{:.10},{:.10},{:.10}",
                af,
                tau,
                s.chi2_total,
                s.red_total,
                s.tt.chi2,
                s.te.chi2,
                s.ee.chi2,
                s.tt.red,
                s.te.red,
                s.ee.red
            )
            .expect("write");
        }
    }

    let mut j = File::create(&json_path).expect("create profile json");
    writeln!(
        j,
        "{{\n  \"inputs\": {{\"class_bin\":\"{}\", \"classy_python\":\"{}\", \"tt_full\":\"{}\", \"te_full\":\"{}\", \"ee_full\":\"{}\"}},\n  \"base\": {{\"h\": {:.12}, \"omega_b\": {:.12}, \"omega_cdm\": {:.12}, \"n_s\": {:.12}, \"a_s\": {:.12e}, \"tau_reio\": {:.12}}},\n  \"scan\": {{\"as_factor_range\":[{:.6},{:.6}], \"tau_range\":[{:.6},{:.6}], \"as_tau_steps\": {}}},\n  \"best\": {{\"chi2_total\": {:.12}, \"reduced_total\": {:.12}, \"a_s\": {:.12e}, \"tau_reio\": {:.12}, \"chi2_tt\": {:.12}, \"chi2_te\": {:.12}, \"chi2_ee\": {:.12}}},\n  \"artifacts\": {{\"profile_csv\":\"{}\"}}\n}}",
        class_bin,
        classy_python,
        tt_full,
        te_full,
        ee_full,
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
        n_as_tau,
        best.0,
        best.1,
        best.2,
        best.3,
        best.4.chi2,
        best.5.chi2,
        best.6.chi2,
        csv_path
    )
    .expect("write summary json");

    println!("wrote {}", csv_path);
    println!("wrote {}", json_path);
    println!(
        "best combined: chi2={:.3} red={:.6} A_s={:.6e} tau={:.5}",
        best.0, best.1, best.2, best.3
    );
}
