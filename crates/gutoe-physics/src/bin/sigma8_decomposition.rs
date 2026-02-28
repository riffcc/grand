//! End-to-end sigma8/S8 pipeline:
//! derived parameters -> CLASS CMB fit -> P(k,z=0) -> sigma8 integral

use gutoe_physics::cmb_class::{
    compare_class_to_planck, read_class_tt_camb, read_planck_tt_csv, run_class, write_class_ini,
    ClassRunInputs,
};
use gutoe_physics::bbn::eta10_from_baryogenesis;
use gutoe_physics::cmb_reionization::derive_tau_reio;
use gutoe_physics::constants::{lambda_cosmological_full_candidate, C};
use gutoe_physics::dark_matter_falsification::OMEGA_BARYON_OBS;
use gutoe_physics::inflation::{inflation_hubble_ratio_structural, scalar_amplitude};
use gutoe_physics::microphysics::MicrophysicsAssumptions;
use gutoe_physics::{evaluate_bbn_gate, evaluate_inflation_gate, BbnWindows, InflationWindows};
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

const SIGMA8_TARGET_PLANCK: f64 = 0.811;
const AS_PLANCK_REF: f64 = 2.10e-9;
const DELTA_STRUCT: f64 = 2.5;
const C_INF_STRUCT: f64 = 1.0 + 1.0 / 66.0;

fn h0_from_lambda_and_omega_lambda(lambda: f64, omega_lambda: f64) -> f64 {
    let meter_per_mpc = 3.085_677_581_491_367e22;
    let h0_s_inv = C * (lambda / (3.0 * omega_lambda)).sqrt();
    h0_s_inv * meter_per_mpc / 1_000.0
}

fn derived_inputs(tau_reio: f64) -> ClassRunInputs {
    let inflation = evaluate_inflation_gate(InflationWindows::default());
    let ratio = (60.0 - DELTA_STRUCT) / 11.0;
    let n = inflation.n_efolds;
    let h_base = inflation_hubble_ratio_structural();
    let h_corr = h_base * (ratio / (60.0 / 11.0)) * C_INF_STRUCT;
    let a_s_corr = scalar_amplitude(n, h_corr);
    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = OMEGA_BARYON_OBS * ratio;
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
        a_s: a_s_corr,
        tau_reio,
    }
}

fn derive_tau(inputs: ClassRunInputs) -> Result<f64, String> {
    let ratio = (60.0 - DELTA_STRUCT) / 11.0;
    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = OMEGA_BARYON_OBS * ratio;
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
    derive_tau_reio(micro, eta10_from_baryogenesis()).map(|r| r.tau_reio)
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

fn run_class_for(
    class_bin: &str,
    out_dir: &Path,
    tag: &str,
    inputs: ClassRunInputs,
) -> Result<(f64, f64, PathBuf, PathBuf), String> {
    let run_dir = out_dir.join(tag);
    if run_dir.exists() {
        fs::remove_dir_all(&run_dir)
            .map_err(|e| format!("clear CLASS run dir {:?}: {e}", run_dir))?;
    }
    fs::create_dir_all(&run_dir).map_err(|e| format!("create CLASS run dir {:?}: {e}", run_dir))?;
    let ini = run_dir.join("in.ini");
    let root = run_dir.join("g_");
    write_class_ini(&ini, &root.to_string_lossy(), 2_500, inputs)?;
    // Ensure matter power output is enabled in the same run.
    let mut ini_txt = fs::read_to_string(&ini).map_err(|e| format!("read ini: {e}"))?;
    if !ini_txt.contains("output =") {
        ini_txt.push_str("\noutput = tCl,mPk\n");
    } else {
        ini_txt = ini_txt.replace("output = tCl", "output = tCl,mPk");
    }
    if !ini_txt.contains("P_k_max_h/Mpc") {
        ini_txt.push_str("\nP_k_max_h/Mpc = 50\nz_pk = 0\n");
    }
    fs::write(&ini, ini_txt).map_err(|e| format!("write ini: {e}"))?;
    run_class(class_bin, &ini)?;
    let cl = find_class_tt_output(&run_dir)?;
    let pk = find_pk_output(&run_dir)?;

    let planck_path = Path::new("crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt");
    let class_tt = read_class_tt_camb(&cl, 2, 2_500)?;
    let planck = read_planck_tt_csv(planck_path)?;
    let fit = compare_class_to_planck(&class_tt, &planck)?;
    Ok((fit.chi2, fit.reduced_chi2, cl, pk))
}

fn sigma8_from_pk(pk_path: &Path, r_hinv_mpc: f64) -> Result<f64, String> {
    let f = File::open(pk_path).map_err(|e| format!("open pk {:?}: {e}", pk_path))?;
    let mut k = Vec::new();
    let mut p = Vec::new();
    for (idx, line) in BufReader::new(f).lines().enumerate() {
        let line = line.map_err(|e| format!("read pk line {}: {e}", idx + 1))?;
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = s.split_whitespace().collect();
        if fields.len() < 2 {
            continue;
        }
        let kv: f64 = match fields[0].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let pv: f64 = match fields[1].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if kv > 0.0 && pv > 0.0 {
            k.push(kv);
            p.push(pv);
        }
    }
    if k.len() < 16 {
        return Err("pk file too short".to_string());
    }
    let mut pairs: Vec<(f64, f64)> = k.into_iter().zip(p).collect();
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut acc = 0.0;
    for w in pairs.windows(2) {
        let (k0, p0) = w[0];
        let (k1, p1) = w[1];
        let f = |kk: f64, pp: f64| {
            let x = kk * r_hinv_mpc;
            let w = if x.abs() < 1e-8 {
                1.0
            } else {
                3.0 * (x.sin() - x * x.cos()) / (x * x * x)
            };
            pp * w * w * kk * kk
        };
        acc += 0.5 * (k1 - k0) * (f(k0, p0) + f(k1, p1));
    }
    let sigma2 = acc / (2.0 * PI * PI);
    Ok(sigma2.max(0.0).sqrt())
}

fn main() {
    let class_bin = std::env::var("GUTOE_CLASS_BIN")
        .unwrap_or_else(|_| "/tmp/class_public/class".to_string());
    let out_dir = std::env::var("GUTOE_SIGMA8_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/sigma8_decomposition".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let base_inputs = derived_inputs(0.054);
    let tau = derive_tau(base_inputs).unwrap_or(0.054);
    let derived = derived_inputs(tau);
    let om0 = (derived.omega_b + derived.omega_cdm) / (derived.h * derived.h);

    let (chi2, red, cl, pk) = match run_class_for(&class_bin, &out, "derived", derived) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("sigma8_decomposition failed in derived run: {e}");
            std::process::exit(2);
        }
    };
    let sigma8 = sigma8_from_pk(&pk, 8.0).unwrap_or(f64::NAN);
    let s8 = sigma8 * (om0 / 0.3).sqrt();

    let mut planck_inputs = derived;
    planck_inputs.a_s = AS_PLANCK_REF;
    let (_chi2_p, _red_p, _cl_p, pk_p) = match run_class_for(&class_bin, &out, "planck_as_ref", planck_inputs) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("sigma8_decomposition failed in planck_as_ref run: {e}");
            std::process::exit(2);
        }
    };
    let sigma8_planck_as = sigma8_from_pk(&pk_p, 8.0).unwrap_or(f64::NAN);

    let as_for_target = derived.a_s * (SIGMA8_TARGET_PLANCK / sigma8).powi(2);
    let mut target_inputs = derived;
    target_inputs.a_s = as_for_target;
    let (_chi2_t, _red_t, _cl_t, pk_t) = match run_class_for(&class_bin, &out, "as_target_sigma8", target_inputs) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("sigma8_decomposition failed in as_target_sigma8 run: {e}");
            std::process::exit(2);
        }
    };
    let sigma8_target_check = sigma8_from_pk(&pk_t, 8.0).unwrap_or(f64::NAN);

    let normalization_piece = sigma8 - sigma8_planck_as;
    let shape_piece = sigma8_planck_as - SIGMA8_TARGET_PLANCK;

    let report = out.join("sigma8_decomposition_report.json");
    let mut f = File::create(&report).expect("create report");
    writeln!(
        f,
        "{{\n  \"inputs\": {{\"class_bin\": \"{}\", \"h\": {:.12}, \"omega_b\": {:.12}, \"omega_cdm\": {:.12}, \"omega_m0\": {:.12}, \"n_s\": {:.12}, \"A_s_derived\": {:.12e}, \"tau_reio\": {:.12}}},\n  \"derived\": {{\"chi2\": {:.12}, \"reduced_chi2\": {:.12}, \"sigma8\": {:.12}, \"S8\": {:.12}, \"cl_path\": \"{}\", \"pk_path\": \"{}\"}},\n  \"decomposition\": {{\"sigma8_target_planck\": {:.12}, \"A_s_planck_ref\": {:.12e}, \"sigma8_at_A_s_planck_ref\": {:.12}, \"A_s_for_sigma8_target_fixed_shape\": {:.12e}, \"sigma8_target_check\": {:.12}, \"normalization_piece_sigma8\": {:.12}, \"shape_piece_sigma8\": {:.12}}}\n}}",
        class_bin,
        derived.h,
        derived.omega_b,
        derived.omega_cdm,
        om0,
        derived.n_s,
        derived.a_s,
        derived.tau_reio,
        chi2,
        red,
        sigma8,
        s8,
        cl.display(),
        pk.display(),
        SIGMA8_TARGET_PLANCK,
        AS_PLANCK_REF,
        sigma8_planck_as,
        as_for_target,
        sigma8_target_check,
        normalization_piece,
        shape_piece
    )
    .expect("write report");

    println!("wrote {}", report.display());
    println!(
        "sigma8={:.6}, S8={:.6}, sigma8_planckAs={:.6}, As_target={:.6e}",
        sigma8, s8, sigma8_planck_as, as_for_target
    );
}
