use gutoe_physics::bbn::{eta10_from_baryogenesis, evaluate_bbn_gate, BbnWindows};
use gutoe_physics::cmb_class::{
    compare_class_to_planck, read_class_dl_camb_column, read_planck_dl_csv, run_class,
    run_classy_fallback, write_class_ini, ClassRunInputs, ClassTtPoint, PlanckTtPoint,
};
use gutoe_physics::cmb_reionization::derive_tau_reio;
use gutoe_physics::constants::{
    lambda_cosmological_full_candidate, ALPHA_LEADING_ORDER, C, LAMBDA_QG, PLANCK_MASS,
};
use gutoe_physics::dark_matter_falsification::OMEGA_BARYON_OBS;
use gutoe_physics::dynamics_map::StandardModelDynamicsMap;
use gutoe_physics::inflation::{
    evaluate_inflation_gate, inflation_hubble_ratio_structural, scalar_amplitude, InflationWindows,
};
use gutoe_physics::microphysics::MicrophysicsAssumptions;
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

const DELTA_STRUCT: f64 = 2.5;
const C_INF_STRUCT: f64 = 1.0 + 1.0 / 66.0;
const ELECTRON_MASS_MEV_OBS: f64 = 0.510_998_950;
const KG_TO_MEV: f64 = 5.609_588_603e29;

#[derive(Debug, Clone, Copy)]
struct ChannelFit {
    binned_chi2: f64,
    binned_red: f64,
    full_chi2: f64,
    full_red: f64,
}

fn h0_from_lambda_and_omega_lambda(lambda: f64, omega_lambda: f64) -> f64 {
    let meter_per_mpc = 3.085_677_581_491_367e22;
    let h0_s_inv = C * (lambda / (3.0 * omega_lambda)).sqrt();
    h0_s_inv * meter_per_mpc / 1_000.0
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
        return Err("no CLASS cl dat files found".to_string());
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

fn find_pk_output(run_dir: &Path) -> Result<PathBuf, String> {
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
                .is_some_and(|name| name.to_ascii_lowercase().contains("pk"))
        })
        .collect();
    if candidates.is_empty() {
        return Err("no CLASS pk dat files found".to_string());
    }
    candidates.sort();
    Ok(candidates[0].clone())
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

fn fit_channel(pred: &[ClassTtPoint], binned: &[PlanckTtPoint], full: &[PlanckTtPoint]) -> Result<ChannelFit, String> {
    let fb = compare_class_to_planck(pred, binned)?;
    let ff = compare_class_to_planck(pred, full)?;
    Ok(ChannelFit {
        binned_chi2: fb.chi2,
        binned_red: fb.reduced_chi2,
        full_chi2: ff.chi2,
        full_red: ff.reduced_chi2,
    })
}

fn env_enabled(var: &str, default: bool) -> bool {
    match std::env::var(var) {
        Ok(v) => {
            let s = v.trim().to_ascii_lowercase();
            !matches!(s.as_str(), "0" | "false" | "off" | "no")
        }
        Err(_) => default,
    }
}

fn electron_mass_pred_flagged_mev() -> f64 {
    let planck_mev = PLANCK_MASS * KG_TO_MEV;
    let ratio = (60.0 - DELTA_STRUCT) / 11.0; // corrected 115/22
    planck_mev
        * ALPHA_LEADING_ORDER.powi(13)
        * ratio.powi(3)
        * C_INF_STRUCT
        * LAMBDA_QG.powi(-3)
}

fn build_full_inputs() -> Result<(ClassRunInputs, f64, f64, f64, f64, String, f64, f64, f64), String>
{
    let infl = evaluate_inflation_gate(InflationWindows::default());
    let bbn = evaluate_bbn_gate(BbnWindows::default());

    let ratio = (60.0 - DELTA_STRUCT) / 11.0;
    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = omega_b0 * ratio;
    let omega_m0 = omega_b0 + omega_cdm0;
    let omega_r0 = 9.0e-5;
    let omega_k0 = 0.0;
    let omega_lambda0 = 1.0 - omega_m0 - omega_r0 - omega_k0;
    let h0 = h0_from_lambda_and_omega_lambda(lambda_cosmological_full_candidate(), omega_lambda0);
    let h = h0 / 100.0;

    let n = infl.n_efolds;
    let h_base = inflation_hubble_ratio_structural();
    let ratio_scale = ratio / (60.0 / 11.0);
    let h_corr = h_base * ratio_scale * C_INF_STRUCT;
    let mut a_s_corr = scalar_amplitude(n, h_corr);
    let hook = std::env::var("GUTOE_ELECTRON_SCALE_HOOK").unwrap_or_else(|_| "none".to_string());
    let mut hook_scale = 1.0_f64;
    if hook.eq_ignore_ascii_case("flagged_as") {
        let me_pred = electron_mass_pred_flagged_mev();
        // Experimental hook: project electron absolute-scale candidate into A_s.
        hook_scale = me_pred / ELECTRON_MASS_MEV_OBS;
        a_s_corr *= hook_scale;
    }

    let micro = MicrophysicsAssumptions {
        h0_km_s_mpc: h0,
        omega_b0,
        omega_m0,
        omega_r0,
        omega_k0,
        omega_lambda0,
        eta10: bbn.eta10,
    };
    let reion = derive_tau_reio(micro, eta10_from_baryogenesis())?;
    let sm = StandardModelDynamicsMap::from_clifford_z3();
    let sin2_mz_coupled = sm.sin2_theta_w_at_mz();
    let sin2_mz_legacy = (3.0 / 13.0) + ALPHA_LEADING_ORDER.powi(2) * (sm.clifford_dim as f64 / 2.0);
    let enable_ew_coupling = env_enabled("GUTOE_EW_CMB_COUPLING", true);
    // Structural coupling from the coupled EW bridge into reionization optics.
    let ew_cmb_coupling = if enable_ew_coupling {
        sin2_mz_coupled / sin2_mz_legacy
    } else {
        1.0
    };
    let tau_reio_coupled = reion.tau_reio * ew_cmb_coupling;

    Ok((
        ClassRunInputs {
            h,
            omega_b: omega_b0 * h * h,
            omega_cdm: omega_cdm0 * h * h,
            omega_k: omega_k0,
            omega_lambda: omega_lambda0,
            n_s: infl.n_s,
            a_s: a_s_corr,
            tau_reio: tau_reio_coupled,
        },
        ratio,
        a_s_corr,
        tau_reio_coupled,
        reion.z_reion_structural,
        hook,
        hook_scale,
        ew_cmb_coupling,
        sin2_mz_coupled,
    ))
}

fn main() {
    let out_dir = std::env::var("GUTOE_CMB_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/cmb_full_derived".to_string());
    let _ = fs::create_dir_all(&out_dir);

    let class_bin = std::env::var("GUTOE_CLASS_BIN")
        .unwrap_or_else(|_| "/tmp/class_public/class".to_string());
    let classy_python =
        std::env::var("GUTOE_CLASSY_PYTHON").unwrap_or_else(|_| "python3".to_string());

    let tt_binned = read_planck_dl_csv(Path::new(
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt",
    ))
    .expect("tt binned");
    let te_binned = read_planck_dl_csv(Path::new(
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TE-binned_R3.02.txt",
    ))
    .expect("te binned");
    let ee_binned = read_planck_dl_csv(Path::new(
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-binned_R3.02.txt",
    ))
    .expect("ee binned");

    let tt_full: Vec<_> = read_planck_dl_csv(Path::new(
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-full_R3.01.txt",
    ))
    .expect("tt full")
    .into_iter()
    .filter(|p| p.ell >= 2 && p.ell <= 2_500)
    .collect();
    let te_full: Vec<_> = read_planck_dl_csv(Path::new(
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-TE-full_R3.01.txt",
    ))
    .expect("te full")
    .into_iter()
    .filter(|p| p.ell >= 2 && p.ell <= 2_500)
    .collect();
    let ee_full: Vec<_> = read_planck_dl_csv(Path::new(
        "crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-full_R3.01.txt",
    ))
    .expect("ee full")
    .into_iter()
    .filter(|p| p.ell >= 2 && p.ell <= 2_500)
    .collect();

    let (inputs, ratio, a_s_corr, tau, z_reion, hook, hook_scale, ew_cmb_coupling, sin2_mz_coupled) =
        build_full_inputs().expect("build inputs");

    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    let run_dir = std::env::temp_dir().join(format!("gutoe_cmb_full_{stamp}"));
    let _ = fs::create_dir_all(&run_dir);
    let ini = run_dir.join("run.ini");
    let root = run_dir.join("g_");
    write_class_ini(&ini, &root.to_string_lossy(), 2_500, inputs).expect("write ini");

    let mut ini_txt = fs::read_to_string(&ini).expect("read ini");
    if !ini_txt.contains("output =") {
        ini_txt.push_str("\noutput = tCl,mPk\n");
    } else {
        ini_txt = ini_txt.replace("output = tCl", "output = tCl,mPk");
    }
    if !ini_txt.contains("P_k_max_h/Mpc") {
        ini_txt.push_str("\nP_k_max_h/Mpc = 50\nz_pk = 0\n");
    }
    fs::write(&ini, ini_txt).expect("rewrite ini");

    let class_cl = match run_class(&class_bin, &ini) {
        Ok(_) => find_class_cl_output(&run_dir).expect("find cl"),
        Err(class_err) => {
            let fb = run_dir.join("g_classy_cl.dat");
            run_classy_fallback(&classy_python, &fb, 2_500, inputs).unwrap_or_else(|e| {
                panic!("CLASS failed ({class_err}); classy fallback failed ({e})")
            });
            fb
        }
    };
    let pk_path = find_pk_output(&run_dir).expect("find pk");

    let tt = read_class_dl_camb_column(&class_cl, 2, 2_500, 2).expect("parse tt");
    let te = read_class_dl_camb_column(&class_cl, 2, 2_500, 5).expect("parse te");
    let ee = read_class_dl_camb_column(&class_cl, 2, 2_500, 3).expect("parse ee");

    let base_tt = fit_channel(&tt, &tt_binned, &tt_full).expect("fit tt");
    let base_te = fit_channel(&te, &te_binned, &te_full).expect("fit te");
    let base_ee = fit_channel(&ee, &ee_binned, &ee_full).expect("fit ee");

    let sigma8 = sigma8_from_pk(&pk_path, 8.0).expect("sigma8");

    let json_path = format!("{out_dir}/cmb_full_derived_report.json");
    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"inputs\": {{\"delta\": {:.12}, \"c_inf\": {:.12}, \"ratio_corrected\": {:.12}, \"h\": {:.12}, \"omega_b\": {:.12}, \"omega_cdm\": {:.12}, \"n_s\": {:.12}, \"A_s\": {:.12e}, \"tau_reio\": {:.12}, \"z_reion_structural\": {:.12}, \"electron_scale_hook\": \"{}\", \"electron_hook_scale\": {:.12}, \"ew_cmb_coupling\": {:.12}, \"sin2_theta_w_mz_coupled\": {:.12}}},\n  \"tt\": {{\"full_red\": {:.12}}},\n  \"te\": {{\"full_red\": {:.12}}},\n  \"ee\": {{\"full_red\": {:.12}}},\n  \"sigma8\": {{\"value\": {:.12}}}\n}}",
        DELTA_STRUCT,
        C_INF_STRUCT,
        ratio,
        inputs.h,
        inputs.omega_b,
        inputs.omega_cdm,
        inputs.n_s,
        a_s_corr,
        tau,
        z_reion,
        hook,
        hook_scale,
        ew_cmb_coupling,
        sin2_mz_coupled,
        base_tt.full_red,
        base_te.full_red,
        base_ee.full_red,
        sigma8
    )
    .expect("write json");

    println!("wrote {}", json_path);
    println!(
        "full-derived (envelope-free default): TT {:.3}, TE {:.3}, EE {:.3}, sigma8={:.6}",
        base_tt.full_red,
        base_te.full_red,
        base_ee.full_red,
        sigma8
    );
}
