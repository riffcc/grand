//! GRAND-344 differential-envelope scale sweep:
//! test whether one ℓ_D correction moves TT/TE/EE in the same direction.

use gutoe_physics::bbn::eta10_from_baryogenesis;
use gutoe_physics::cmb_class::{
    compare_class_to_planck, read_class_dl_camb_column, read_planck_dl_csv, run_class,
    run_classy_fallback, write_class_ini, ClassRunInputs, ClassTtPoint, PlanckTtPoint,
};
use gutoe_physics::cmb_differential::{
    apply_differential_envelope, default_transition_scale, estimate_class_ell_diff,
    DifferentialEnvelope,
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

#[derive(Debug, Clone, Copy)]
struct ChannelFit {
    binned_red: f64,
    full_red: f64,
}

#[derive(Debug, Clone, Copy)]
struct SweepRow {
    ell_struct: f64,
    tt_full_red: f64,
    te_full_red: f64,
    ee_full_red: f64,
    tt_binned_red: f64,
    te_binned_red: f64,
    ee_binned_red: f64,
    ee_band_pull: f64,
}

fn h0_from_lambda_and_omega_lambda(lambda: f64, omega_lambda: f64) -> f64 {
    let meter_per_mpc = 3.085_677_581_491_367e22;
    let h0_s_inv = C * (lambda / (3.0 * omega_lambda)).sqrt();
    h0_s_inv * meter_per_mpc / 1_000.0
}

fn class_inputs_with_tau(tau_reio: f64) -> ClassRunInputs {
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

fn derived_tau_and_micro() -> Result<(f64, f64, MicrophysicsAssumptions), String> {
    let base = class_inputs_with_tau(0.054);
    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_cdm0;
    let bbn = evaluate_bbn_gate(BbnWindows::default());
    let micro = MicrophysicsAssumptions {
        h0_km_s_mpc: base.h * 100.0,
        omega_b0,
        omega_m0,
        omega_r0: 9.0e-5,
        omega_k0: base.omega_k,
        omega_lambda0: base.omega_lambda,
        eta10: bbn.eta10,
    };
    let reion = derive_tau_reio(micro, eta10_from_baryogenesis())?;
    Ok((reion.tau_reio, reion.z_reion_structural, micro))
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

fn parse_env_u32(name: &str, default: u32) -> u32 {
    std::env::var(name)
        .ok()
        .and_then(|x| x.parse::<u32>().ok())
        .unwrap_or(default)
}

fn parse_env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|x| x.parse::<f64>().ok())
        .unwrap_or(default)
}

fn fit_channel(
    pred: &[ClassTtPoint],
    binned: &[PlanckTtPoint],
    full: &[PlanckTtPoint],
) -> Result<ChannelFit, String> {
    let fb = compare_class_to_planck(pred, binned)?;
    let ff = compare_class_to_planck(pred, full)?;
    Ok(ChannelFit {
        binned_red: fb.reduced_chi2,
        full_red: ff.reduced_chi2,
    })
}

fn mean_abs_pull_band(
    pred: &[ClassTtPoint],
    obs: &[PlanckTtPoint],
    lo: u32,
    hi: u32,
) -> Result<f64, String> {
    let fit = compare_class_to_planck(pred, obs)?;
    let mut acc = 0.0;
    let mut n = 0usize;
    for r in fit.residuals {
        if r.ell >= lo && r.ell <= hi {
            acc += r.pull.abs();
            n += 1;
        }
    }
    if n == 0 {
        return Err(format!("no points in band [{lo},{hi}]"));
    }
    Ok(acc / n as f64)
}

fn main() {
    let out_dir = std::env::var("GUTOE_CMB_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);

    let class_bin = std::env::var("GUTOE_CLASS_BIN").unwrap_or_else(|_| "class".to_string());
    let classy_python =
        std::env::var("GUTOE_CLASSY_PYTHON").unwrap_or_else(|_| "python3".to_string());

    let tt_binned_path = std::env::var("GUTOE_PLANCK_TT_BINNED").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt".to_string()
    });
    let tt_full_path = std::env::var("GUTOE_PLANCK_TT_FULL").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-full_R3.01.txt".to_string()
    });
    let te_binned_path = std::env::var("GUTOE_PLANCK_TE_BINNED").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TE-binned_R3.02.txt".to_string()
    });
    let te_full_path = std::env::var("GUTOE_PLANCK_TE_FULL").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TE-full_R3.01.txt".to_string()
    });
    let ee_binned_path = std::env::var("GUTOE_PLANCK_EE_BINNED").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-binned_R3.02.txt".to_string()
    });
    let ee_full_path = std::env::var("GUTOE_PLANCK_EE_FULL").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-full_R3.01.txt".to_string()
    });

    let ell_fit_min = parse_env_u32("GUTOE_DIFF_ELLFIT_MIN", 1200);
    let ell_fit_max = parse_env_u32("GUTOE_DIFF_ELLFIT_MAX", 2200);
    let band_lo = parse_env_u32("GUTOE_EE_BAND_LO", 1200);
    let band_hi = parse_env_u32("GUTOE_EE_BAND_HI", 1600);

    let ell_min = parse_env_f64("GUTOE_DIFF_SWEEP_MIN", 1100.0);
    let ell_max = parse_env_f64("GUTOE_DIFF_SWEEP_MAX", 1800.0);
    let ell_step = parse_env_f64("GUTOE_DIFF_SWEEP_STEP", 10.0).max(1.0);

    let (tau, z_reion, _micro) = match derived_tau_and_micro() {
        Ok(v) => v,
        Err(e) => {
            eprintln!("failed to derive tau/micro assumptions: {e}");
            std::process::exit(2);
        }
    };

    let inputs = class_inputs_with_tau(tau);
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    let run_dir = std::env::temp_dir().join(format!("gutoe_cmb_diff_sweep_{stamp}"));
    let _ = fs::create_dir_all(&run_dir);
    let ini = run_dir.join("run.ini");
    let root = run_dir.join("g_");

    if let Err(e) = write_class_ini(&ini, &root.to_string_lossy(), 2_500, inputs) {
        eprintln!("failed to write CLASS ini: {e}");
        std::process::exit(2);
    }

    let class_cl = match run_class(&class_bin, &ini) {
        Ok(_) => match find_class_cl_output(&run_dir) {
            Ok(v) => v,
            Err(e) => {
                eprintln!("failed to locate CLASS output: {e}");
                std::process::exit(2);
            }
        },
        Err(class_err) => {
            let fb = run_dir.join("g_classy_cl.dat");
            if let Err(e) = run_classy_fallback(&classy_python, &fb, 2_500, inputs) {
                eprintln!("CLASS failed ({class_err}); classy fallback failed ({e})");
                std::process::exit(2);
            }
            fb
        }
    };

    let tt = read_class_dl_camb_column(&class_cl, 2, 2_500, 2).expect("parse TT");
    let ee = read_class_dl_camb_column(&class_cl, 2, 2_500, 3).expect("parse EE");
    let te = read_class_dl_camb_column(&class_cl, 2, 2_500, 5).expect("parse TE");

    let tt_binned = read_planck_dl_csv(Path::new(&tt_binned_path)).expect("parse TT binned");
    let te_binned = read_planck_dl_csv(Path::new(&te_binned_path)).expect("parse TE binned");
    let ee_binned = read_planck_dl_csv(Path::new(&ee_binned_path)).expect("parse EE binned");

    let tt_full: Vec<_> = read_planck_dl_csv(Path::new(&tt_full_path))
        .expect("parse TT full")
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();
    let te_full: Vec<_> = read_planck_dl_csv(Path::new(&te_full_path))
        .expect("parse TE full")
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();
    let ee_full: Vec<_> = read_planck_dl_csv(Path::new(&ee_full_path))
        .expect("parse EE full")
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();

    let base_tt = fit_channel(&tt, &tt_binned, &tt_full).expect("fit TT base");
    let base_te = fit_channel(&te, &te_binned, &te_full).expect("fit TE base");
    let base_ee = fit_channel(&ee, &ee_binned, &ee_full).expect("fit EE base");

    let ell_class = estimate_class_ell_diff(&tt, ell_fit_min, ell_fit_max).expect("fit class ell_diff");
    let ell_transition = default_transition_scale();

    let mut rows = Vec::<SweepRow>::new();

    // Baseline row for convenience
    rows.push(SweepRow {
        ell_struct: f64::INFINITY,
        tt_full_red: base_tt.full_red,
        te_full_red: base_te.full_red,
        ee_full_red: base_ee.full_red,
        tt_binned_red: base_tt.binned_red,
        te_binned_red: base_te.binned_red,
        ee_binned_red: base_ee.binned_red,
        ee_band_pull: mean_abs_pull_band(&ee, &ee_binned, band_lo, band_hi)
            .expect("base band pull"),
    });

    let mut ell = ell_min;
    while ell <= ell_max + 1e-9 {
        let env = DifferentialEnvelope {
            ell_diff_struct: ell,
            ell_diff_class: ell_class,
            ell_transition,
        };
        let tt_corr = apply_differential_envelope(&tt, env);
        let te_corr = apply_differential_envelope(&te, env);
        let ee_corr = apply_differential_envelope(&ee, env);

        let fit_tt = fit_channel(&tt_corr, &tt_binned, &tt_full).expect("fit TT corr");
        let fit_te = fit_channel(&te_corr, &te_binned, &te_full).expect("fit TE corr");
        let fit_ee = fit_channel(&ee_corr, &ee_binned, &ee_full).expect("fit EE corr");
        let ee_band = mean_abs_pull_band(&ee_corr, &ee_binned, band_lo, band_hi)
            .expect("corr band pull");

        rows.push(SweepRow {
            ell_struct: ell,
            tt_full_red: fit_tt.full_red,
            te_full_red: fit_te.full_red,
            ee_full_red: fit_ee.full_red,
            tt_binned_red: fit_tt.binned_red,
            te_binned_red: fit_te.binned_red,
            ee_binned_red: fit_ee.binned_red,
            ee_band_pull: ee_band,
        });

        ell += ell_step;
    }

    let base = rows[0];
    let mut best_all: Option<SweepRow> = None;
    let mut best_sum = f64::INFINITY;

    for row in rows.iter().skip(1) {
        let all_improve = row.tt_full_red <= base.tt_full_red
            && row.te_full_red <= base.te_full_red
            && row.ee_full_red <= base.ee_full_red;
        if all_improve {
            let s = row.tt_full_red + row.te_full_red + row.ee_full_red;
            if s < best_sum {
                best_sum = s;
                best_all = Some(*row);
            }
        }
    }

    let csv_path = format!("{out_dir}/cmb_differential_sweep.csv");
    let mut csv = File::create(&csv_path).expect("create csv");
    writeln!(csv, "ell_struct,tt_full_red,te_full_red,ee_full_red,tt_binned_red,te_binned_red,ee_binned_red,ee_band_pull").expect("write csv header");
    for r in &rows {
        if r.ell_struct.is_infinite() {
            writeln!(
                csv,
                "baseline,{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}",
                r.tt_full_red,
                r.te_full_red,
                r.ee_full_red,
                r.tt_binned_red,
                r.te_binned_red,
                r.ee_binned_red,
                r.ee_band_pull,
            )
            .expect("write csv row");
        } else {
            writeln!(
                csv,
                "{:.3},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}",
                r.ell_struct,
                r.tt_full_red,
                r.te_full_red,
                r.ee_full_red,
                r.tt_binned_red,
                r.te_binned_red,
                r.ee_binned_red,
                r.ee_band_pull,
            )
            .expect("write csv row");
        }
    }

    let txt_path = format!("{out_dir}/cmb_differential_sweep.txt");
    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[inputs]").expect("write");
    writeln!(txt, "tau_reio = {:.9}", tau).expect("write");
    writeln!(txt, "z_reion_structural = {:.9}", z_reion).expect("write");
    writeln!(txt, "ell_diff_class = {:.9}", ell_class).expect("write");
    writeln!(txt, "ell_transition = {:.9}", ell_transition).expect("write");
    writeln!(txt, "sweep = [{:.3}, {:.3}] step {:.3}", ell_min, ell_max, ell_step).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[baseline]").expect("write");
    writeln!(txt, "tt_full_red = {:.9}", base.tt_full_red).expect("write");
    writeln!(txt, "te_full_red = {:.9}", base.te_full_red).expect("write");
    writeln!(txt, "ee_full_red = {:.9}", base.ee_full_red).expect("write");
    writeln!(txt, "ee_band_pull = {:.9}", base.ee_band_pull).expect("write");
    writeln!(txt).expect("write");

    if let Some(best) = best_all {
        writeln!(txt, "[best_all_channels]").expect("write");
        writeln!(txt, "ell_struct = {:.9}", best.ell_struct).expect("write");
        writeln!(txt, "tt_full_red = {:.9}", best.tt_full_red).expect("write");
        writeln!(txt, "te_full_red = {:.9}", best.te_full_red).expect("write");
        writeln!(txt, "ee_full_red = {:.9}", best.ee_full_red).expect("write");
        writeln!(txt, "sum_full_red = {:.9}", best_sum).expect("write");
        writeln!(txt, "ee_band_pull = {:.9}", best.ee_band_pull).expect("write");
    } else {
        writeln!(txt, "[best_all_channels]").expect("write");
        writeln!(txt, "none = true").expect("write");
    }

    println!("wrote {}", csv_path);
    println!("wrote {}", txt_path);
    if let Some(best) = best_all {
        println!(
            "Best single-scale all-channel improvement: ell_struct={:.1}; TT {:.3}->{:.3}, TE {:.3}->{:.3}, EE {:.3}->{:.3}",
            best.ell_struct,
            base.tt_full_red,
            best.tt_full_red,
            base.te_full_red,
            best.te_full_red,
            base.ee_full_red,
            best.ee_full_red,
        );
    } else {
        println!(
            "No single ell_struct in sweep range improved TT/TE/EE simultaneously (baseline TT {:.3}, TE {:.3}, EE {:.3}).",
            base.tt_full_red, base.te_full_red, base.ee_full_red,
        );
    }
}
