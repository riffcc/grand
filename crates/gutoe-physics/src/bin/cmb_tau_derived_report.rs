//! GRAND-343 / GRAND-355:
//! Derive tau_reio from structural reionization and compare TT fit impact
//! against the explicit tau assumption lane.

use gutoe_physics::bbn::eta10_from_baryogenesis;
use gutoe_physics::cmb_class::{
    compare_class_to_planck, read_class_tt_camb, read_planck_tt_csv, run_class,
    run_classy_fallback, write_class_ini, ClassRunInputs,
};
use gutoe_physics::cmb_reionization::{derive_tau_reio, ReionizationDerived};
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

fn micro_assumptions() -> MicrophysicsAssumptions {
    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_cdm0;
    let omega_r0 = 9.0e-5;
    let omega_k0 = 0.0;
    let omega_lambda0 = 1.0 - omega_m0 - omega_r0 - omega_k0;
    let h0 = h0_from_lambda_and_omega_lambda(lambda_cosmological_full_candidate(), omega_lambda0);
    MicrophysicsAssumptions {
        h0_km_s_mpc: h0,
        omega_b0,
        omega_m0,
        omega_r0,
        omega_k0,
        omega_lambda0,
        eta10: eta10_from_baryogenesis(),
    }
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

fn run_fit(
    class_bin: &str,
    classy_python: &str,
    inputs: ClassRunInputs,
    planck_binned_path: &Path,
    planck_full_path: &Path,
    label: &str,
) -> (f64, f64, f64, f64) {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let run_dir = std::env::temp_dir().join(format!("gutoe_tau_cmp_{label}_{stamp}"));
    let _ = fs::create_dir_all(&run_dir);
    let ini_path = run_dir.join("run.ini");
    let root = run_dir.join("g_");

    write_class_ini(&ini_path, &root.to_string_lossy(), 2_500, inputs).expect("write ini");
    let class_tt_path = match run_class(class_bin, &ini_path) {
        Ok(_) => find_class_tt_output(&run_dir).expect("find class output"),
        Err(class_err) => {
            let fallback_path = run_dir.join("g_classy_cl.dat");
            run_classy_fallback(classy_python, &fallback_path, 2_500, inputs).unwrap_or_else(
                |e| panic!("CLASS failed ({class_err}); classy fallback failed ({e})"),
            );
            fallback_path
        }
    };

    let class_tt = read_class_tt_camb(&class_tt_path, 2, 2_500).expect("parse class");
    let planck_binned = read_planck_tt_csv(planck_binned_path).expect("parse binned");
    let planck_full = read_planck_tt_csv(planck_full_path).expect("parse full");
    let planck_full_cut: Vec<_> = planck_full
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();

    let fb = compare_class_to_planck(&class_tt, &planck_binned).expect("fit binned");
    let ff = compare_class_to_planck(&class_tt, &planck_full_cut).expect("fit full");
    (fb.chi2, fb.reduced_chi2, ff.chi2, ff.reduced_chi2)
}

fn main() {
    let out_dir = std::env::var("GUTOE_CMB_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);

    let class_bin = std::env::var("GUTOE_CLASS_BIN").unwrap_or_else(|_| "class".to_string());
    let classy_python =
        std::env::var("GUTOE_CLASSY_PYTHON").unwrap_or_else(|_| "python3".to_string());
    let planck_binned = std::env::var("GUTOE_PLANCK_TT_BINNED").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt".to_string()
    });
    let planck_full = std::env::var("GUTOE_PLANCK_TT_FULL").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-full_R3.01.txt".to_string()
    });
    let tau_assumed = std::env::var("GUTOE_TAU_ASSUMED")
        .ok()
        .and_then(|x| x.parse::<f64>().ok())
        .unwrap_or(0.054);

    let reion: ReionizationDerived =
        derive_tau_reio(micro_assumptions(), eta10_from_baryogenesis()).expect("derive tau");

    let assumed_inputs = class_inputs_with_tau(tau_assumed);
    let derived_inputs = class_inputs_with_tau(reion.tau_reio);

    let (chi2b_assumed, redb_assumed, chi2f_assumed, redf_assumed) = run_fit(
        &class_bin,
        &classy_python,
        assumed_inputs,
        Path::new(&planck_binned),
        Path::new(&planck_full),
        "assumed",
    );
    let (chi2b_derived, redb_derived, chi2f_derived, redf_derived) = run_fit(
        &class_bin,
        &classy_python,
        derived_inputs,
        Path::new(&planck_binned),
        Path::new(&planck_full),
        "derived",
    );

    let txt_path = format!("{out_dir}/cmb_tau_derived_report.txt");
    let json_path = format!("{out_dir}/cmb_tau_derived_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[reionization]").expect("write");
    writeln!(txt, "z_reion_structural = {:.9}", reion.z_reion_structural).expect("write");
    writeln!(txt, "tau_reio_derived = {:.9}", reion.tau_reio).expect("write");
    writeln!(txt, "suppression_e2tau = {:.9}", reion.suppression_e2tau).expect("write");
    writeln!(txt, "tau_reio_assumed = {:.9}", tau_assumed).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[fit_binned]").expect("write");
    writeln!(txt, "assumed_chi2 = {:.9}", chi2b_assumed).expect("write");
    writeln!(txt, "assumed_reduced_chi2 = {:.9}", redb_assumed).expect("write");
    writeln!(txt, "derived_chi2 = {:.9}", chi2b_derived).expect("write");
    writeln!(txt, "derived_reduced_chi2 = {:.9}", redb_derived).expect("write");
    writeln!(txt, "delta_chi2 = {:.9}", chi2b_derived - chi2b_assumed).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[fit_full]").expect("write");
    writeln!(txt, "assumed_chi2 = {:.9}", chi2f_assumed).expect("write");
    writeln!(txt, "assumed_reduced_chi2 = {:.9}", redf_assumed).expect("write");
    writeln!(txt, "derived_chi2 = {:.9}", chi2f_derived).expect("write");
    writeln!(txt, "derived_reduced_chi2 = {:.9}", redf_derived).expect("write");
    writeln!(txt, "delta_chi2 = {:.9}", chi2f_derived - chi2f_assumed).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"reionization\": {{\"z_reion_structural\": {:.12}, \"tau_reio_derived\": {:.12}, \"suppression_e2tau\": {:.12}, \"tau_reio_assumed\": {:.12}}},\n  \"fit_binned\": {{\"assumed_chi2\": {:.12}, \"assumed_reduced_chi2\": {:.12}, \"derived_chi2\": {:.12}, \"derived_reduced_chi2\": {:.12}, \"delta_chi2\": {:.12}}},\n  \"fit_full\": {{\"assumed_chi2\": {:.12}, \"assumed_reduced_chi2\": {:.12}, \"derived_chi2\": {:.12}, \"derived_reduced_chi2\": {:.12}, \"delta_chi2\": {:.12}}}\n}}",
        reion.z_reion_structural,
        reion.tau_reio,
        reion.suppression_e2tau,
        tau_assumed,
        chi2b_assumed,
        redb_assumed,
        chi2b_derived,
        redb_derived,
        chi2b_derived - chi2b_assumed,
        chi2f_assumed,
        redf_assumed,
        chi2f_derived,
        redf_derived,
        chi2f_derived - chi2f_assumed,
    )
    .expect("write json");

    println!("wrote {}", txt_path);
    println!("wrote {}", json_path);
    println!(
        "Derived tau: z_reion={:.3}, tau={:.6}, e^(-2tau)={:.6}",
        reion.z_reion_structural, reion.tau_reio, reion.suppression_e2tau
    );
    println!(
        "Binned chi2: assumed {:.1} -> derived {:.1} (Δ={:+.1})",
        chi2b_assumed,
        chi2b_derived,
        chi2b_derived - chi2b_assumed
    );
    println!(
        "Full chi2: assumed {:.1} -> derived {:.1} (Δ={:+.1})",
        chi2f_assumed,
        chi2f_derived,
        chi2f_derived - chi2f_assumed
    );
}
