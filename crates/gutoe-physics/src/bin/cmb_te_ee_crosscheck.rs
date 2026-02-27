//! GRAND-344/355 next step:
//! Cross-channel validation (TT/TE/EE) with structurally derived tau_reio,
//! using zero additional channel-specific tuning.

use gutoe_physics::bbn::eta10_from_baryogenesis;
use gutoe_physics::cmb_class::{
    compare_class_to_planck, read_class_dl_camb_column, read_planck_dl_csv, run_class,
    run_classy_fallback, write_class_ini, ClassRunInputs,
};
use gutoe_physics::cmb_reionization::derive_tau_reio;
use gutoe_physics::constants::{
    lambda_cosmological_full_candidate, C, DARK_TO_VISIBLE_GEOMETRIC_RATIO,
};
use gutoe_physics::dark_matter_falsification::OMEGA_BARYON_OBS;
use gutoe_physics::microphysics::MicrophysicsAssumptions;
use gutoe_physics::{evaluate_inflation_gate, InflationWindows};
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Debug, Clone, Copy)]
struct ChannelFit {
    binned_chi2: f64,
    binned_red: f64,
    full_chi2: f64,
    full_red: f64,
}

#[derive(Debug, Clone, Copy)]
struct SpectraFitSet {
    tt: ChannelFit,
    te: ChannelFit,
    ee: ChannelFit,
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

fn derived_tau_from_structure(inputs: ClassRunInputs) -> Result<(f64, f64), String> {
    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_cdm0;
    let eta10 = eta10_from_baryogenesis();
    let micro = MicrophysicsAssumptions {
        h0_km_s_mpc: inputs.h * 100.0,
        omega_b0,
        omega_m0,
        omega_r0: 9.0e-5,
        omega_k0: inputs.omega_k,
        omega_lambda0: inputs.omega_lambda,
        eta10,
    };
    let reion = derive_tau_reio(micro, eta10)?;
    Ok((reion.tau_reio, reion.z_reion_structural))
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

fn run_fit_set(
    class_bin: &str,
    classy_python: &str,
    inputs: ClassRunInputs,
    tt_binned: &Path,
    tt_full: &Path,
    te_binned: &Path,
    te_full: &Path,
    ee_binned: &Path,
    ee_full: &Path,
    tag: &str,
) -> Result<SpectraFitSet, String> {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let run_dir = std::env::temp_dir().join(format!("gutoe_teee_{tag}_{stamp}"));
    let _ = fs::create_dir_all(&run_dir);
    let ini = run_dir.join("run.ini");
    let root = run_dir.join("g_");

    write_class_ini(&ini, &root.to_string_lossy(), 2_500, inputs)?;
    let class_cl_path = match run_class(class_bin, &ini) {
        Ok(_) => find_class_cl_output(&run_dir)?,
        Err(class_err) => {
            let fallback_path = run_dir.join("g_classy_cl.dat");
            run_classy_fallback(classy_python, &fallback_path, 2_500, inputs).map_err(|e| {
                format!(
                    "CLASS failed ({class_err}); classy fallback failed ({e}) for tag={tag}"
                )
            })?;
            fallback_path
        }
    };

    let tt = read_class_dl_camb_column(&class_cl_path, 2, 2_500, 2)?;
    let ee = read_class_dl_camb_column(&class_cl_path, 2, 2_500, 3)?;
    let te = read_class_dl_camb_column(&class_cl_path, 2, 2_500, 5)?;

    let tt_b = read_planck_dl_csv(tt_binned)?;
    let tt_f: Vec<_> = read_planck_dl_csv(tt_full)?
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();
    let te_b = read_planck_dl_csv(te_binned)?;
    let te_f: Vec<_> = read_planck_dl_csv(te_full)?
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();
    let ee_b = read_planck_dl_csv(ee_binned)?;
    let ee_f: Vec<_> = read_planck_dl_csv(ee_full)?
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();

    let fit_tt_b = compare_class_to_planck(&tt, &tt_b)?;
    let fit_tt_f = compare_class_to_planck(&tt, &tt_f)?;
    let fit_te_b = compare_class_to_planck(&te, &te_b)?;
    let fit_te_f = compare_class_to_planck(&te, &te_f)?;
    let fit_ee_b = compare_class_to_planck(&ee, &ee_b)?;
    let fit_ee_f = compare_class_to_planck(&ee, &ee_f)?;

    Ok(SpectraFitSet {
        tt: ChannelFit {
            binned_chi2: fit_tt_b.chi2,
            binned_red: fit_tt_b.reduced_chi2,
            full_chi2: fit_tt_f.chi2,
            full_red: fit_tt_f.reduced_chi2,
        },
        te: ChannelFit {
            binned_chi2: fit_te_b.chi2,
            binned_red: fit_te_b.reduced_chi2,
            full_chi2: fit_te_f.chi2,
            full_red: fit_te_f.reduced_chi2,
        },
        ee: ChannelFit {
            binned_chi2: fit_ee_b.chi2,
            binned_red: fit_ee_b.reduced_chi2,
            full_chi2: fit_ee_f.chi2,
            full_red: fit_ee_f.reduced_chi2,
        },
    })
}

fn main() {
    let out_dir = std::env::var("GUTOE_CMB_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);

    let class_bin = std::env::var("GUTOE_CLASS_BIN").unwrap_or_else(|_| "class".to_string());
    let classy_python =
        std::env::var("GUTOE_CLASSY_PYTHON").unwrap_or_else(|_| "python3".to_string());

    let tt_binned = std::env::var("GUTOE_PLANCK_TT_BINNED").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt".to_string()
    });
    let tt_full = std::env::var("GUTOE_PLANCK_TT_FULL").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-full_R3.01.txt".to_string()
    });
    let te_binned = std::env::var("GUTOE_PLANCK_TE_BINNED").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TE-binned_R3.02.txt".to_string()
    });
    let te_full = std::env::var("GUTOE_PLANCK_TE_FULL").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TE-full_R3.01.txt".to_string()
    });
    let ee_binned = std::env::var("GUTOE_PLANCK_EE_BINNED").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-binned_R3.02.txt".to_string()
    });
    let ee_full = std::env::var("GUTOE_PLANCK_EE_FULL").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-full_R3.01.txt".to_string()
    });

    let tau_assumed = std::env::var("GUTOE_TAU_ASSUMED")
        .ok()
        .and_then(|x| x.parse::<f64>().ok())
        .unwrap_or(0.054);

    let baseline_inputs = class_inputs_with_tau(tau_assumed);
    let (tau_derived, z_reion_structural) =
        derived_tau_from_structure(baseline_inputs).expect("derive tau_reio structurally");
    let derived_inputs = class_inputs_with_tau(tau_derived);

    let fit_assumed = run_fit_set(
        &class_bin,
        &classy_python,
        baseline_inputs,
        Path::new(&tt_binned),
        Path::new(&tt_full),
        Path::new(&te_binned),
        Path::new(&te_full),
        Path::new(&ee_binned),
        Path::new(&ee_full),
        "assumed",
    )
    .expect("assumed fit set");

    let fit_derived = run_fit_set(
        &class_bin,
        &classy_python,
        derived_inputs,
        Path::new(&tt_binned),
        Path::new(&tt_full),
        Path::new(&te_binned),
        Path::new(&te_full),
        Path::new(&ee_binned),
        Path::new(&ee_full),
        "derived",
    )
    .expect("derived fit set");

    let txt_path = format!("{out_dir}/cmb_te_ee_crosscheck.txt");
    let json_path = format!("{out_dir}/cmb_te_ee_crosscheck.json");
    let mut txt = File::create(&txt_path).expect("create txt");

    writeln!(txt, "[tau]").expect("write");
    writeln!(txt, "tau_assumed = {:.9}", tau_assumed).expect("write");
    writeln!(txt, "tau_derived = {:.9}", tau_derived).expect("write");
    writeln!(txt, "z_reion_structural = {:.9}", z_reion_structural).expect("write");
    writeln!(txt).expect("write");

    let print_channel = |f: &mut File, name: &str, a: ChannelFit, d: ChannelFit| {
        writeln!(f, "[{}]", name).expect("write");
        writeln!(f, "assumed_binned_chi2 = {:.9}", a.binned_chi2).expect("write");
        writeln!(f, "assumed_binned_red = {:.9}", a.binned_red).expect("write");
        writeln!(f, "derived_binned_chi2 = {:.9}", d.binned_chi2).expect("write");
        writeln!(f, "derived_binned_red = {:.9}", d.binned_red).expect("write");
        writeln!(f, "delta_binned_chi2 = {:.9}", d.binned_chi2 - a.binned_chi2).expect("write");
        writeln!(f, "assumed_full_chi2 = {:.9}", a.full_chi2).expect("write");
        writeln!(f, "derived_full_chi2 = {:.9}", d.full_chi2).expect("write");
        writeln!(f, "delta_full_chi2 = {:.9}", d.full_chi2 - a.full_chi2).expect("write");
        writeln!(f, "assumed_full_red = {:.9}", a.full_red).expect("write");
        writeln!(f, "derived_full_red = {:.9}", d.full_red).expect("write");
        writeln!(f).expect("write");
    };

    print_channel(&mut txt, "TT", fit_assumed.tt, fit_derived.tt);
    print_channel(&mut txt, "TE", fit_assumed.te, fit_derived.te);
    print_channel(&mut txt, "EE", fit_assumed.ee, fit_derived.ee);

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"tau\": {{\"assumed\": {:.12}, \"derived\": {:.12}, \"z_reion_structural\": {:.12}}},\n  \"tt\": {{\"assumed_binned_chi2\": {:.12}, \"assumed_binned_red\": {:.12}, \"derived_binned_chi2\": {:.12}, \"derived_binned_red\": {:.12}, \"delta_binned_chi2\": {:.12}, \"assumed_full_chi2\": {:.12}, \"derived_full_chi2\": {:.12}, \"delta_full_chi2\": {:.12}, \"assumed_full_red\": {:.12}, \"derived_full_red\": {:.12}}},\n  \"te\": {{\"assumed_binned_chi2\": {:.12}, \"assumed_binned_red\": {:.12}, \"derived_binned_chi2\": {:.12}, \"derived_binned_red\": {:.12}, \"delta_binned_chi2\": {:.12}, \"assumed_full_chi2\": {:.12}, \"derived_full_chi2\": {:.12}, \"delta_full_chi2\": {:.12}, \"assumed_full_red\": {:.12}, \"derived_full_red\": {:.12}}},\n  \"ee\": {{\"assumed_binned_chi2\": {:.12}, \"assumed_binned_red\": {:.12}, \"derived_binned_chi2\": {:.12}, \"derived_binned_red\": {:.12}, \"delta_binned_chi2\": {:.12}, \"assumed_full_chi2\": {:.12}, \"derived_full_chi2\": {:.12}, \"delta_full_chi2\": {:.12}, \"assumed_full_red\": {:.12}, \"derived_full_red\": {:.12}}}\n}}",
        tau_assumed,
        tau_derived,
        z_reion_structural,
        fit_assumed.tt.binned_chi2,
        fit_assumed.tt.binned_red,
        fit_derived.tt.binned_chi2,
        fit_derived.tt.binned_red,
        fit_derived.tt.binned_chi2 - fit_assumed.tt.binned_chi2,
        fit_assumed.tt.full_chi2,
        fit_derived.tt.full_chi2,
        fit_derived.tt.full_chi2 - fit_assumed.tt.full_chi2,
        fit_assumed.tt.full_red,
        fit_derived.tt.full_red,
        fit_assumed.te.binned_chi2,
        fit_assumed.te.binned_red,
        fit_derived.te.binned_chi2,
        fit_derived.te.binned_red,
        fit_derived.te.binned_chi2 - fit_assumed.te.binned_chi2,
        fit_assumed.te.full_chi2,
        fit_derived.te.full_chi2,
        fit_derived.te.full_chi2 - fit_assumed.te.full_chi2,
        fit_assumed.te.full_red,
        fit_derived.te.full_red,
        fit_assumed.ee.binned_chi2,
        fit_assumed.ee.binned_red,
        fit_derived.ee.binned_chi2,
        fit_derived.ee.binned_red,
        fit_derived.ee.binned_chi2 - fit_assumed.ee.binned_chi2,
        fit_assumed.ee.full_chi2,
        fit_derived.ee.full_chi2,
        fit_derived.ee.full_chi2 - fit_assumed.ee.full_chi2,
        fit_assumed.ee.full_red,
        fit_derived.ee.full_red,
    )
    .expect("write json");

    println!("wrote {}", txt_path);
    println!("wrote {}", json_path);
    println!(
        "tau derived: {:.6} (z_reion={:.3}) | TT full red {:.3} -> {:.3} | TE full red {:.3} -> {:.3} | EE full red {:.3} -> {:.3}",
        tau_derived,
        z_reion_structural,
        fit_assumed.tt.full_red,
        fit_derived.tt.full_red,
        fit_assumed.te.full_red,
        fit_derived.te.full_red,
        fit_assumed.ee.full_red,
        fit_derived.ee.full_red,
    );
}
