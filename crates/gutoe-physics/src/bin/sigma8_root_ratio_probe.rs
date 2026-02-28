use gutoe_physics::bbn::{evaluate_bbn_gate, BbnWindows};
use gutoe_physics::cmb_reionization::derive_tau_reio;
use gutoe_physics::constants::{lambda_cosmological_full_candidate, C};
use gutoe_physics::dark_matter_falsification::OMEGA_BARYON_OBS;
use gutoe_physics::inflation::{
    evaluate_inflation_gate, inflation_hubble_ratio_structural, scalar_amplitude, AS_OBSERVED,
    InflationWindows,
};
use gutoe_physics::microphysics::MicrophysicsAssumptions;
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

const OMEGA_CDM_H2_TARGET: f64 = 0.1200;
const SIGMA8_TARGET: f64 = 0.811;
const RATIO_BASE: f64 = 60.0 / 11.0;

#[derive(Clone, Copy)]
struct Inputs {
    h: f64,
    omega_b_h2: f64,
    omega_cdm_h2: f64,
    omega_k: f64,
    n_s: f64,
    a_s: f64,
    tau_reio: f64,
}

fn h0_from_lambda_and_omega_lambda(lambda: f64, omega_lambda: f64) -> f64 {
    let meter_per_mpc = 3.085_677_581_491_367e22;
    let h0_s_inv = C * (lambda / (3.0 * omega_lambda)).sqrt();
    h0_s_inv * meter_per_mpc / 1_000.0
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

fn run_class_sigma8(class_bin: &str, out_dir: &Path, tag: &str, i: Inputs) -> Result<f64, String> {
    let run_dir = out_dir.join(tag);
    fs::create_dir_all(&run_dir).map_err(|e| format!("mkdir {:?}: {e}", run_dir))?;
    let ini = run_dir.join("in.ini");
    let root = run_dir.join("g_");
    let ini_txt = format!(
        "h = {h}\nomega_b = {ob}\nomega_cdm = {oc}\nOmega_k = {ok}\nA_s = {as_}\nn_s = {ns}\ntau_reio = {tau}\noutput = mPk\nP_k_max_h/Mpc = 50\nz_pk = 0\nroot = {root}\n",
        h = i.h,
        ob = i.omega_b_h2,
        oc = i.omega_cdm_h2,
        ok = i.omega_k,
        as_ = i.a_s,
        ns = i.n_s,
        tau = i.tau_reio,
        root = root.to_string_lossy()
    );
    fs::write(&ini, ini_txt).map_err(|e| format!("write ini {:?}: {e}", ini))?;
    let status = Command::new(class_bin)
        .arg(&ini)
        .status()
        .map_err(|e| format!("run class '{}': {e}", class_bin))?;
    if !status.success() {
        return Err(format!("CLASS failed status {status}"));
    }
    let mut pk_files: Vec<PathBuf> = fs::read_dir(&run_dir)
        .map_err(|e| format!("read run dir: {e}"))?
        .filter_map(|e| e.ok().map(|x| x.path()))
        .filter(|p| {
            p.extension()
                .and_then(|x| x.to_str())
                .is_some_and(|ext| ext.eq_ignore_ascii_case("dat"))
                && p.file_name()
                    .and_then(|x| x.to_str())
                    .is_some_and(|n| n.to_ascii_lowercase().contains("pk"))
        })
        .collect();
    pk_files.sort();
    let pk = pk_files
        .first()
        .ok_or_else(|| "no pk dat file".to_string())?
        .clone();
    sigma8_from_pk(&pk, 8.0)
}

fn build_inputs(delta: f64, c_inf: f64) -> Result<(Inputs, f64, f64), String> {
    let infl = evaluate_inflation_gate(InflationWindows::default());
    let bbn = evaluate_bbn_gate(BbnWindows::default());

    let ratio = (60.0 - delta) / 11.0;
    let ratio_scale = ratio / RATIO_BASE;

    let n = infl.n_efolds;
    let h_base = inflation_hubble_ratio_structural();
    let h_new = h_base * ratio_scale * c_inf;
    let a_s_new = scalar_amplitude(n, h_new);

    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_cdm0 = omega_b0 * ratio;
    let omega_m0 = omega_b0 + omega_cdm0;
    let omega_r0 = 9.0e-5;
    let omega_k0 = 0.0;
    let omega_lambda0 = 1.0 - omega_m0 - omega_r0 - omega_k0;
    let h0 = h0_from_lambda_and_omega_lambda(lambda_cosmological_full_candidate(), omega_lambda0);
    let h = h0 / 100.0;

    let micro = MicrophysicsAssumptions {
        h0_km_s_mpc: h0,
        omega_b0,
        omega_m0,
        omega_r0,
        omega_k0,
        omega_lambda0,
        eta10: bbn.eta10,
    };
    let tau_reio = derive_tau_reio(micro, bbn.eta10)
        .map_err(|e| format!("derive_tau_reio failed: {e}"))?
        .tau_reio;

    let out = Inputs {
        h,
        omega_b_h2: omega_b0 * h * h,
        omega_cdm_h2: omega_cdm0 * h * h,
        omega_k: omega_k0,
        n_s: infl.n_s,
        a_s: a_s_new,
        tau_reio,
    };
    Ok((out, a_s_new, ratio))
}

fn main() {
    let class_bin = std::env::var("GUTOE_CLASS_BIN")
        .unwrap_or_else(|_| "/tmp/class_public/class".to_string());
    let out_dir =
        std::env::var("GUTOE_SIGMA8_OUT").unwrap_or_else(|_| "/tmp/bh_renders/sigma8_root_ratio_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let (base_inputs, base_as, _base_ratio) = match build_inputs(0.0, 1.0) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    };
    let sigma8_base = run_class_sigma8(&class_bin, &out, "base", base_inputs).unwrap_or(f64::NAN);

    // Delta that matches A_s with c_inf fixed to 1.
    let delta_for_as_only = 60.0 * (1.0 - (AS_OBSERVED / base_as).sqrt());
    let (as_inputs, as_val, _) = build_inputs(delta_for_as_only, 1.0).expect("as-only build");
    let sigma8_as_only = run_class_sigma8(&class_bin, &out, "delta_for_as_only", as_inputs).unwrap_or(f64::NAN);

    // Delta that matches omega_cdm h^2 target with c_inf fixed to 1 (linearized on ratio).
    let delta_for_oc_only = 60.0 - 11.0 * (OMEGA_CDM_H2_TARGET / base_inputs.omega_b_h2);
    let (oc_inputs, oc_as, _) = build_inputs(delta_for_oc_only, 1.0).expect("oc-only build");
    let sigma8_oc_only = run_class_sigma8(&class_bin, &out, "delta_for_oc_only", oc_inputs).unwrap_or(f64::NAN);
    let (d25_inputs, d25_as, _) = build_inputs(2.5, 1.0).expect("d25 build");
    let sigma8_d25_no_cinf =
        run_class_sigma8(&class_bin, &out, "delta_5_over_2_no_cinf", d25_inputs).unwrap_or(f64::NAN);
    let fixed_cinf_66 = 1.015151515152_f64;
    let (d25_c66_inputs, d25_c66_as, _) = build_inputs(2.5, fixed_cinf_66).expect("d25 c66 build");
    let sigma8_d25_c66 = run_class_sigma8(
        &class_bin,
        &out,
        "delta_5_over_2_cinf_1p015151515152",
        d25_c66_inputs,
    )
    .unwrap_or(f64::NAN);

    // Shared root scan with c_inf auto-solved so A_s is exactly restored at each delta.
    let mut best_sigma = (f64::INFINITY, 0.0, 1.0, f64::NAN, f64::NAN, f64::NAN);
    for step in 180..=320 {
        let delta = step as f64 * 0.01;
        let (tmp_i, tmp_as, _) = match build_inputs(delta, 1.0) {
            Ok(v) => v,
            Err(_) => continue,
        };
        if !(tmp_as > 0.0) {
            continue;
        }
        let c_inf = (AS_OBSERVED / tmp_as).sqrt();
        let (i, as_fixed, _) = match build_inputs(delta, c_inf) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let sigma = match run_class_sigma8(
            &class_bin,
            &out,
            &format!("scan_d{:04}", (delta * 100.0).round() as i32),
            i,
        ) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let err = (sigma - SIGMA8_TARGET).abs();
        if err < best_sigma.0 {
            best_sigma = (err, delta, c_inf, sigma, as_fixed, i.omega_cdm_h2);
        }
        let _ = tmp_i; // keep structure explicit, avoid accidental simplification.
    }

    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let candidates: [(&str, f64); 7] = [
        ("5/2", 2.5),
        ("8/3", 8.0 / 3.0),
        ("1+phi", 1.0 + phi),
        ("3-1/e", 3.0 - 1.0 / std::f64::consts::E),
        ("3", 3.0),
        ("delta_omega_match", delta_for_oc_only),
        ("delta_prev_best_coarse", 2.6),
    ];
    let mut cand_rows: Vec<(&str, f64, f64, f64, f64)> = Vec::new();
    for (name, delta) in candidates {
        let (_, as_tmp, _) = match build_inputs(delta, 1.0) {
            Ok(v) => v,
            Err(_) => {
                cand_rows.push((name, delta, f64::NAN, f64::NAN, f64::NAN));
                continue;
            }
        };
        let c_inf = (AS_OBSERVED / as_tmp).sqrt();
        let (i, as_fixed, _) = match build_inputs(delta, c_inf) {
            Ok(v) => v,
            Err(_) => {
                cand_rows.push((name, delta, c_inf, f64::NAN, f64::NAN));
                continue;
            }
        };
        let sigma = run_class_sigma8(
            &class_bin,
            &out,
            &format!("cand_{}", name.replace('/', "_")),
            i,
        )
        .unwrap_or(f64::NAN);
        cand_rows.push((name, delta, c_inf, sigma, as_fixed));
    }

    let report = out.join("sigma8_root_ratio_probe_report.json");
    let mut f = File::create(&report).expect("create report");
    writeln!(f, "{{").expect("write");
    writeln!(
        f,
        "  \"base\": {{\"A_s\": {:.12e}, \"omega_b_h2\": {:.12}, \"omega_cdm_h2\": {:.12}, \"sigma8\": {:.12}}},",
        base_as, base_inputs.omega_b_h2, base_inputs.omega_cdm_h2, sigma8_base
    )
    .expect("write");
    writeln!(
        f,
        "  \"single_delta_without_c_inf\": {{\"delta_for_A_s\": {:.12}, \"A_s_at_delta_for_A_s\": {:.12e}, \"omega_cdm_h2_at_delta_for_A_s\": {:.12}, \"sigma8_at_delta_for_A_s\": {:.12}, \"delta_for_omega_cdm_h2\": {:.12}, \"A_s_at_delta_for_omega_cdm_h2\": {:.12e}, \"omega_cdm_h2_at_delta_for_omega_cdm_h2\": {:.12}, \"sigma8_at_delta_for_omega_cdm_h2\": {:.12}}},",
        delta_for_as_only,
        as_val,
        as_inputs.omega_cdm_h2,
        sigma8_as_only,
        delta_for_oc_only,
        oc_as,
        oc_inputs.omega_cdm_h2,
        sigma8_oc_only
    )
    .expect("write");
    writeln!(
        f,
        "  \"delta_5_over_2_without_c_inf\": {{\"delta\": 2.500000000000, \"A_s\": {:.12e}, \"omega_cdm_h2\": {:.12}, \"sigma8\": {:.12}, \"abs_sigma8_error\": {:.12}}},",
        d25_as,
        d25_inputs.omega_cdm_h2,
        sigma8_d25_no_cinf,
        (sigma8_d25_no_cinf - SIGMA8_TARGET).abs()
    )
    .expect("write");
    writeln!(
        f,
        "  \"delta_5_over_2_with_fixed_c_inf_1p015151515152\": {{\"delta\": 2.500000000000, \"c_inf\": {:.12}, \"A_s\": {:.12e}, \"omega_cdm_h2\": {:.12}, \"sigma8\": {:.12}, \"abs_sigma8_error\": {:.12}}},",
        fixed_cinf_66,
        d25_c66_as,
        d25_c66_inputs.omega_cdm_h2,
        sigma8_d25_c66,
        (sigma8_d25_c66 - SIGMA8_TARGET).abs()
    )
    .expect("write");
    writeln!(
        f,
        "  \"shared_delta_with_c_inf_autosolved\": {{\"sigma8_target\": {:.12}, \"best_delta\": {:.12}, \"best_c_inf\": {:.12}, \"sigma8\": {:.12}, \"A_s\": {:.12e}, \"omega_cdm_h2\": {:.12}, \"abs_sigma8_error\": {:.12}}}",
        SIGMA8_TARGET, best_sigma.1, best_sigma.2, best_sigma.3, best_sigma.4, best_sigma.5, best_sigma.0
    )
    .expect("write");
    writeln!(f, ",").expect("write");
    writeln!(f, "  \"candidates\": [").expect("write");
    for (idx, row) in cand_rows.iter().enumerate() {
        let comma = if idx + 1 == cand_rows.len() { "" } else { "," };
        writeln!(
            f,
            "    {{\"name\": \"{}\", \"delta\": {:.12}, \"c_inf\": {:.12}, \"sigma8\": {:.12}, \"A_s\": {:.12e}, \"abs_sigma8_error\": {:.12}}}{}",
            row.0,
            row.1,
            row.2,
            row.3,
            row.4,
            (row.3 - SIGMA8_TARGET).abs(),
            comma
        )
        .expect("write");
    }
    writeln!(f, "  ]").expect("write");
    writeln!(f, "}}").expect("write");

    println!("wrote {}", report.display());
    println!(
        "base sigma8={:.6}; delta(A_s-only)={:.4} -> sigma8={:.6}; delta(ocdm-only)={:.4} -> sigma8={:.6}",
        sigma8_base, delta_for_as_only, sigma8_as_only, delta_for_oc_only, sigma8_oc_only
    );
    println!(
        "best shared delta + c_inf: delta={:.4}, c_inf={:.6}, sigma8={:.6}, A_s={:.6e}, omega_cdm_h2={:.6}",
        best_sigma.1, best_sigma.2, best_sigma.3, best_sigma.4, best_sigma.5
    );
    for row in cand_rows {
        println!(
            "candidate {}: delta={:.6}, c_inf={:.6}, sigma8={:.6}, |err|={:.6}",
            row.0,
            row.1,
            row.2,
            row.3,
            (row.3 - SIGMA8_TARGET).abs()
        );
    }
}
