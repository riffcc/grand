//! GRAND-344 / GRAND-355:
//! Derive CMB damping envelope from in-framework microphysics (no TT fitting),
//! apply it to the baseline CLASS TT spectrum, and report fit deltas.

use gutoe_physics::cmb_class::{
    compare_class_to_planck, read_class_tt_camb, read_planck_tt_csv, run_class,
    run_classy_fallback, write_class_ini, ClassRunInputs, ClassTtPoint,
};
use gutoe_physics::cmb_damping::{apply_microphysics_damping, derive_silk_damping_envelope};
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

fn derived_microphysics_assumptions(class_i: ClassRunInputs) -> MicrophysicsAssumptions {
    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_cdm0;
    let omega_r0 = 9.0e-5;
    let omega_k0 = class_i.omega_k;
    let omega_lambda0 = class_i.omega_lambda;
    let bbn = evaluate_bbn_gate(BbnWindows::default());
    MicrophysicsAssumptions {
        h0_km_s_mpc: class_i.h * 100.0,
        omega_b0,
        omega_m0,
        omega_r0,
        omega_k0,
        omega_lambda0,
        eta10: bbn.eta10,
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

fn parse_env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|x| x.parse::<f64>().ok())
        .unwrap_or(default)
}

fn write_overlay_csv(path: &str, baseline: &[ClassTtPoint], damped: &[ClassTtPoint]) -> Result<(), String> {
    let mut f = File::create(path).map_err(|e| format!("create overlay csv: {e}"))?;
    writeln!(f, "ell,baseline_uk2,damped_uk2,damping_factor")
        .map_err(|e| format!("write overlay csv header: {e}"))?;

    for (b, d) in baseline.iter().zip(damped.iter()) {
        let factor = if b.d_ell_tt_uk2 > 0.0 {
            d.d_ell_tt_uk2 / b.d_ell_tt_uk2
        } else {
            0.0
        };
        writeln!(
            f,
            "{},{:.12},{:.12},{:.12}",
            b.ell, b.d_ell_tt_uk2, d.d_ell_tt_uk2, factor
        )
        .map_err(|e| format!("write overlay csv row: {e}"))?;
    }
    Ok(())
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

    let tau_reio = parse_env_f64("GUTOE_TAU_REIO", 0.054);

    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    let run_dir = std::env::temp_dir().join(format!("gutoe_cmb_damping_{stamp}"));
    let _ = fs::create_dir_all(&run_dir);
    let ini_path = run_dir.join("run.ini");
    let root = run_dir.join("g_");

    let class_inputs = derived_class_inputs(tau_reio);
    write_class_ini(&ini_path, &root.to_string_lossy(), 2_500, class_inputs).expect("write ini");

    let (backend_used, class_tt_path) = match run_class(&class_bin, &ini_path) {
        Ok(_) => (
            "class",
            find_class_tt_output(&run_dir).expect("find class output"),
        ),
        Err(class_err) => {
            let fallback_path = run_dir.join("g_classy_cl.dat");
            run_classy_fallback(&classy_python, &fallback_path, 2_500, class_inputs).unwrap_or_else(
                |e| panic!("CLASS failed ({class_err}); classy fallback failed ({e})"),
            );
            ("classy", fallback_path)
        }
    };

    let baseline_tt = read_class_tt_camb(&class_tt_path, 2, 2_500).expect("parse class output");
    let micro_a = derived_microphysics_assumptions(class_inputs);
    let damping = derive_silk_damping_envelope(micro_a).expect("derive damping envelope");
    let damped_tt = apply_microphysics_damping(&baseline_tt, damping);

    let planck_binned = read_planck_tt_csv(Path::new(&planck_binned_path)).expect("parse binned");
    let planck_full = read_planck_tt_csv(Path::new(&planck_full_path)).expect("parse full");
    let planck_full_cut: Vec<_> = planck_full
        .into_iter()
        .filter(|p| p.ell >= 2 && p.ell <= 2_500)
        .collect();

    let fit_base_binned = compare_class_to_planck(&baseline_tt, &planck_binned).expect("fit base binned");
    let fit_damp_binned = compare_class_to_planck(&damped_tt, &planck_binned).expect("fit damp binned");

    let fit_base_full = compare_class_to_planck(&baseline_tt, &planck_full_cut).expect("fit base full");
    let fit_damp_full = compare_class_to_planck(&damped_tt, &planck_full_cut).expect("fit damp full");

    let txt_path = format!("{out_dir}/cmb_derived_damping_report.txt");
    let json_path = format!("{out_dir}/cmb_derived_damping_report.json");
    let overlay_csv = format!("{out_dir}/cmb_derived_damping_overlay.csv");

    write_overlay_csv(&overlay_csv, &baseline_tt, &damped_tt).expect("write overlay csv");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[inputs]").expect("write");
    writeln!(txt, "backend_used = {}", backend_used).expect("write");
    writeln!(txt, "class_bin = {}", class_bin).expect("write");
    writeln!(txt, "classy_python = {}", classy_python).expect("write");
    writeln!(txt, "planck_binned = {}", planck_binned_path).expect("write");
    writeln!(txt, "planck_full = {}", planck_full_path).expect("write");
    writeln!(txt, "tau_reio = {:.12}", class_inputs.tau_reio).expect("write");
    writeln!(txt, "tau_reio_is_assumption = true").expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[derived_damping]").expect("write");
    writeln!(txt, "z_star = {:.6}", damping.z_star).expect("write");
    writeln!(txt, "sigma_z = {:.6}", damping.sigma_z).expect("write");
    writeln!(txt, "diffusion_length_mpc = {:.9}", damping.diffusion_length_mpc).expect("write");
    writeln!(txt, "visibility_width_mpc = {:.9}", damping.visibility_width_mpc).expect("write");
    writeln!(txt, "k_diff_mpc_inv = {:.9}", damping.k_diff_mpc_inv).expect("write");
    writeln!(txt, "ell_diff = {:.6}", damping.ell_diff).expect("write");
    writeln!(txt, "ell_vis = {:.6}", damping.ell_vis).expect("write");
    writeln!(txt, "d_m_star_mpc = {:.6}", damping.d_m_star_mpc).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[fit_binned]").expect("write");
    writeln!(txt, "baseline_chi2 = {:.9}", fit_base_binned.chi2).expect("write");
    writeln!(txt, "baseline_reduced_chi2 = {:.9}", fit_base_binned.reduced_chi2).expect("write");
    writeln!(txt, "damped_chi2 = {:.9}", fit_damp_binned.chi2).expect("write");
    writeln!(txt, "damped_reduced_chi2 = {:.9}", fit_damp_binned.reduced_chi2).expect("write");
    writeln!(txt, "delta_chi2 = {:.9}", fit_damp_binned.chi2 - fit_base_binned.chi2).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[fit_full]").expect("write");
    writeln!(txt, "baseline_chi2 = {:.9}", fit_base_full.chi2).expect("write");
    writeln!(txt, "baseline_reduced_chi2 = {:.9}", fit_base_full.reduced_chi2).expect("write");
    writeln!(txt, "damped_chi2 = {:.9}", fit_damp_full.chi2).expect("write");
    writeln!(txt, "damped_reduced_chi2 = {:.9}", fit_damp_full.reduced_chi2).expect("write");
    writeln!(txt, "delta_chi2 = {:.9}", fit_damp_full.chi2 - fit_base_full.chi2).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"inputs\": {{\"backend\": \"{}\", \"class_bin\": \"{}\", \"classy_python\": \"{}\", \"planck_binned\": \"{}\", \"planck_full\": \"{}\", \"tau_reio\": {:.12}, \"tau_reio_is_assumption\": true}},\n  \"derived_damping\": {{\"z_star\": {:.12}, \"sigma_z\": {:.12}, \"diffusion_length_mpc\": {:.12}, \"visibility_width_mpc\": {:.12}, \"k_diff_mpc_inv\": {:.12}, \"ell_diff\": {:.12}, \"ell_vis\": {:.12}, \"d_m_star_mpc\": {:.12}}},\n  \"fit_binned\": {{\"baseline_chi2\": {:.12}, \"baseline_reduced_chi2\": {:.12}, \"damped_chi2\": {:.12}, \"damped_reduced_chi2\": {:.12}, \"delta_chi2\": {:.12}}},\n  \"fit_full\": {{\"baseline_chi2\": {:.12}, \"baseline_reduced_chi2\": {:.12}, \"damped_chi2\": {:.12}, \"damped_reduced_chi2\": {:.12}, \"delta_chi2\": {:.12}}},\n  \"artifacts\": {{\"overlay_csv\": \"{}\"}}\n}}",
        backend_used,
        class_bin,
        classy_python,
        planck_binned_path,
        planck_full_path,
        class_inputs.tau_reio,
        damping.z_star,
        damping.sigma_z,
        damping.diffusion_length_mpc,
        damping.visibility_width_mpc,
        damping.k_diff_mpc_inv,
        damping.ell_diff,
        damping.ell_vis,
        damping.d_m_star_mpc,
        fit_base_binned.chi2,
        fit_base_binned.reduced_chi2,
        fit_damp_binned.chi2,
        fit_damp_binned.reduced_chi2,
        fit_damp_binned.chi2 - fit_base_binned.chi2,
        fit_base_full.chi2,
        fit_base_full.reduced_chi2,
        fit_damp_full.chi2,
        fit_damp_full.reduced_chi2,
        fit_damp_full.chi2 - fit_base_full.chi2,
        overlay_csv,
    )
    .expect("write json");

    println!("wrote {}", txt_path);
    println!("wrote {}", json_path);
    println!("wrote {}", overlay_csv);
    println!(
        "Derived damping scales: z*= {:.2}, ell_diff= {:.1}, ell_vis= {:.1}",
        damping.z_star, damping.ell_diff, damping.ell_vis
    );
    println!(
        "Binned chi2: baseline {:.1} -> damped {:.1} (Δ={:+.1})",
        fit_base_binned.chi2,
        fit_damp_binned.chi2,
        fit_damp_binned.chi2 - fit_base_binned.chi2
    );
    println!(
        "Full chi2: baseline {:.1} -> damped {:.1} (Δ={:+.1})",
        fit_base_full.chi2,
        fit_damp_full.chi2,
        fit_damp_full.chi2 - fit_base_full.chi2
    );
}
