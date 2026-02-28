use serde_json::{json, Value};
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::process::Command;

fn run_with_env(cmd: &str, args: &[&str], env: &[(&str, &str)]) -> Result<(), String> {
    let mut command = Command::new(cmd);
    command.args(args);
    for (k, v) in env {
        command.env(k, v);
    }
    let status = command
        .status()
        .map_err(|e| format!("spawn {cmd}: {e}"))?;
    if !status.success() {
        return Err(format!("command failed: {} {} (status={status})", cmd, args.join(" ")));
    }
    Ok(())
}

fn run(cmd: &str, args: &[&str]) -> Result<(), String> {
    run_with_env(cmd, args, &[])
}

fn run_or_exit(label: &str, cmd: &str, args: &[&str], env: &[(&str, &str)]) {
    if let Err(e) = run_with_env(cmd, args, env) {
        eprintln!("global_gate: {label} failed: {e}");
        std::process::exit(2);
    }
}

fn read_json(path: &str) -> Result<Value, String> {
    let s = fs::read_to_string(path).map_err(|e| format!("read {path}: {e}"))?;
    serde_json::from_str(&s).map_err(|e| format!("parse {path}: {e}"))
}

fn v_f64(v: &Value, path: &[&str]) -> Result<f64, String> {
    let mut cur = v;
    for k in path {
        cur = cur
            .get(*k)
            .ok_or_else(|| format!("missing key {} in {}", k, path.join(".")))?;
    }
    cur.as_f64()
        .ok_or_else(|| format!("non-f64 at {}", path.join(".")))
}

fn v_bool(v: &Value, path: &[&str]) -> Result<bool, String> {
    let mut cur = v;
    for k in path {
        cur = cur
            .get(*k)
            .ok_or_else(|| format!("missing key {} in {}", k, path.join(".")))?;
    }
    cur.as_bool()
        .ok_or_else(|| format!("non-bool at {}", path.join(".")))
}

fn main() {
    const PMNS_GAIN_ENV: &str = "GUTOE_LEPTOGENESIS_PMNS_GAIN";
    const PMNS_GAIN_BASELINE: &str = "0";
    const PMNS_GAIN_STRUCTURAL: &str = "1";
    const PMNS_BARYO_DELTA_MIN: f64 = 1.0e-15;
    const PMNS_CMB_DELTA_MIN: f64 = 1.0e-12;

    let out_dir = std::env::var("GUTOE_GLOBAL_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/global_gate".to_string());
    let _ = fs::create_dir_all(&out_dir);

    // Run fresh artifacts for each critical lane.
    if let Err(e) = run("cargo", &["run", "-q", "-p", "gutoe-em", "--bin", "flavor_ci_gate"]) {
        eprintln!("global_gate: flavor_ci_gate failed: {e}");
        std::process::exit(2);
    }
    if let Err(e) = run(
        "cargo",
        &["run", "-q", "-p", "gutoe-physics", "--bin", "proton_mass_report"],
    ) {
        eprintln!("global_gate: proton_mass_report failed: {e}");
        std::process::exit(2);
    }
    if let Err(e) = run(
        "cargo",
        &["run", "-q", "-p", "gutoe-physics", "--bin", "alpha_web_ci_report"],
    ) {
        eprintln!("global_gate: alpha_web_ci_report failed: {e}");
        std::process::exit(2);
    }
    if let Err(e) = run(
        "cargo",
        &[
            "run",
            "-q",
            "-p",
            "gutoe-physics",
            "--bin",
            "chiral_symmetry_breaking_report",
        ],
    ) {
        eprintln!("global_gate: chiral_symmetry_breaking_report failed: {e}");
        std::process::exit(2);
    }

    run_or_exit(
        "baryogenesis_report[gain=0]",
        "cargo",
        &["run", "-q", "-p", "gutoe-physics", "--bin", "baryogenesis_report"],
        &[(PMNS_GAIN_ENV, PMNS_GAIN_BASELINE)],
    );
    let baryo0 = read_json("/tmp/bh_renders/baryogenesis_report.json").expect("baryogenesis gain0 json");

    run_or_exit(
        "baryogenesis_report[gain=1]",
        "cargo",
        &["run", "-q", "-p", "gutoe-physics", "--bin", "baryogenesis_report"],
        &[(PMNS_GAIN_ENV, PMNS_GAIN_STRUCTURAL)],
    );
    let baryo1 = read_json("/tmp/bh_renders/baryogenesis_report.json").expect("baryogenesis gain1 json");

    run_or_exit(
        "cmb_full_derived_report[gain=0]",
        "cargo",
        &["run", "-q", "-p", "gutoe-physics", "--bin", "cmb_full_derived_report"],
        &[(PMNS_GAIN_ENV, PMNS_GAIN_BASELINE)],
    );
    let cmb0 =
        read_json("/tmp/bh_renders/cmb_full_derived/cmb_full_derived_report.json").expect("cmb gain0 json");

    run_or_exit(
        "cmb_full_derived_report[gain=1]",
        "cargo",
        &["run", "-q", "-p", "gutoe-physics", "--bin", "cmb_full_derived_report"],
        &[(PMNS_GAIN_ENV, PMNS_GAIN_STRUCTURAL)],
    );
    let cmb =
        read_json("/tmp/bh_renders/cmb_full_derived/cmb_full_derived_report.json").expect("cmb gain1 json");

    run_or_exit(
        "sigma8_decomposition[gain=1]",
        "cargo",
        &["run", "-q", "-p", "gutoe-physics", "--bin", "sigma8_decomposition"],
        &[(PMNS_GAIN_ENV, PMNS_GAIN_STRUCTURAL)],
    );

    let flavor = read_json("/tmp/bh_renders/flavor_ci_gate.json").expect("flavor json");
    let proton = read_json("/tmp/bh_renders/proton_mass_report/proton_mass_report.json").expect("proton json");
    let alpha = read_json("/tmp/bh_renders/alpha_web_ci_report/alpha_web_ci_report.json").expect("alpha json");
    let chiral = read_json("/tmp/bh_renders/chiral_symmetry_breaking/chiral_symmetry_breaking_report.json")
        .expect("chiral json");
    let sigma = read_json("/tmp/bh_renders/sigma8_decomposition/sigma8_decomposition_report.json")
        .expect("sigma json");

    let pmns_corr_res = v_f64(&flavor, &["pmns_theta23_improvement", "corrected_abs_residual_deg"]).unwrap();
    let pmns_corr_pass = v_bool(&flavor, &["pmns_theta23_improvement", "pass"]).unwrap() && pmns_corr_res <= 0.01;

    let eta0 = v_f64(&baryo0, &["eta_predicted"]).expect("baryogenesis eta gain0");
    let eta1 = v_f64(&baryo1, &["eta_predicted"]).expect("baryogenesis eta gain1");
    let pmns_baryo_eta_delta = eta1 - eta0;
    let pmns_baryo_coupling_pass = pmns_baryo_eta_delta > PMNS_BARYO_DELTA_MIN;

    let tt0 = v_f64(&cmb0, &["tt", "full_red"]).expect("tt gain0");
    let te0 = v_f64(&cmb0, &["te", "full_red"]).expect("te gain0");
    let ee0 = v_f64(&cmb0, &["ee", "full_red"]).expect("ee gain0");
    let sigma8_cmb0 = v_f64(&cmb0, &["sigma8", "value"]).expect("sigma8 gain0");

    let tt = v_f64(&cmb, &["tt", "full_red"]).expect("tt gain1");
    let te = v_f64(&cmb, &["te", "full_red"]).expect("te gain1");
    let ee = v_f64(&cmb, &["ee", "full_red"]).expect("ee gain1");
    let sigma8_cmb = v_f64(&cmb, &["sigma8", "value"]).expect("sigma8 gain1");

    let pmns_cmb_tt_delta = tt - tt0;
    let pmns_cmb_te_delta = te - te0;
    let pmns_cmb_ee_delta = ee - ee0;
    let pmns_cmb_sigma8_delta = sigma8_cmb - sigma8_cmb0;
    let pmns_cmb_delta_max = pmns_cmb_tt_delta
        .abs()
        .max(pmns_cmb_te_delta.abs())
        .max(pmns_cmb_ee_delta.abs())
        .max(pmns_cmb_sigma8_delta.abs());
    let pmns_cmb_response_pass = pmns_cmb_delta_max > PMNS_CMB_DELTA_MIN;
    let pmns_propagation_pass = pmns_baryo_coupling_pass && pmns_cmb_response_pass;

    let proton_rel_err = v_f64(&proton, &["route_a_electron_anchor", "proton_rel_error"]).unwrap();
    let proton_pass = proton_rel_err.abs() <= 1.0e-3;

    let alpha_pass = v_bool(&alpha, &["ci_gate", "passes_all"]).unwrap();
    let chiral_pass = v_bool(&chiral, &["gate", "passes_all"]).unwrap();

    let cmb_pass = tt <= 1.30 && te <= 1.20 && ee <= 1.10;

    let sigma8_dec = v_f64(&sigma, &["derived", "sigma8"]).unwrap();
    let sigma8_match = (sigma8_dec - sigma8_cmb).abs();
    let sigma8_target = 0.811;
    let sigma8_pass = (sigma8_cmb - sigma8_target).abs() <= 0.005 && sigma8_match <= 1.0e-9;

    let overall_pass = pmns_corr_pass
        && pmns_propagation_pass
        && proton_pass
        && alpha_pass
        && chiral_pass
        && cmb_pass
        && sigma8_pass;

    let txt_path = Path::new(&out_dir).join("global_gate_report.txt");
    let json_path = Path::new(&out_dir).join("global_gate_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[global_gate]").ok();
    writeln!(txt, "pmns_theta23_corrected_residual_deg = {:.12}", pmns_corr_res).ok();
    writeln!(txt, "pmns_theta23_pass = {}", pmns_corr_pass).ok();
    writeln!(txt, "pmns_gain0_eta_predicted = {:.12e}", eta0).ok();
    writeln!(txt, "pmns_gain1_eta_predicted = {:.12e}", eta1).ok();
    writeln!(txt, "pmns_baryogenesis_eta_delta = {:.12e}", pmns_baryo_eta_delta).ok();
    writeln!(txt, "pmns_baryogenesis_coupling_pass = {}", pmns_baryo_coupling_pass).ok();
    writeln!(txt, "pmns_cmb_tt_delta = {:.12e}", pmns_cmb_tt_delta).ok();
    writeln!(txt, "pmns_cmb_te_delta = {:.12e}", pmns_cmb_te_delta).ok();
    writeln!(txt, "pmns_cmb_ee_delta = {:.12e}", pmns_cmb_ee_delta).ok();
    writeln!(txt, "pmns_cmb_sigma8_delta = {:.12e}", pmns_cmb_sigma8_delta).ok();
    writeln!(txt, "pmns_cmb_delta_max_abs = {:.12e}", pmns_cmb_delta_max).ok();
    writeln!(txt, "pmns_cmb_response_pass = {}", pmns_cmb_response_pass).ok();
    writeln!(txt, "pmns_propagation_pass = {}", pmns_propagation_pass).ok();
    writeln!(txt, "proton_rel_error = {:.12e}", proton_rel_err).ok();
    writeln!(txt, "proton_pass = {}", proton_pass).ok();
    writeln!(txt, "alpha_web_ci_pass = {}", alpha_pass).ok();
    writeln!(txt, "chiral_symmetry_breaking_pass = {}", chiral_pass).ok();
    writeln!(txt, "cmb_tt_red = {:.12}", tt).ok();
    writeln!(txt, "cmb_te_red = {:.12}", te).ok();
    writeln!(txt, "cmb_ee_red = {:.12}", ee).ok();
    writeln!(txt, "cmb_pass = {}", cmb_pass).ok();
    writeln!(txt, "sigma8_cmb = {:.12}", sigma8_cmb).ok();
    writeln!(txt, "sigma8_decomp = {:.12}", sigma8_dec).ok();
    writeln!(txt, "sigma8_match_abs = {:.3e}", sigma8_match).ok();
    writeln!(txt, "sigma8_pass = {}", sigma8_pass).ok();
    writeln!(txt, "overall_pass = {}", overall_pass).ok();

    let report = json!({
        "pmns": {
            "theta23_corrected_residual_deg": pmns_corr_res,
            "pass": pmns_corr_pass
        },
        "pmns_propagation": {
            "gain0": {
                "eta_predicted": eta0,
                "tt_full_red": tt0,
                "te_full_red": te0,
                "ee_full_red": ee0,
                "sigma8": sigma8_cmb0
            },
            "gain1": {
                "eta_predicted": eta1,
                "tt_full_red": tt,
                "te_full_red": te,
                "ee_full_red": ee,
                "sigma8": sigma8_cmb
            },
            "delta": {
                "eta_predicted": pmns_baryo_eta_delta,
                "tt_full_red": pmns_cmb_tt_delta,
                "te_full_red": pmns_cmb_te_delta,
                "ee_full_red": pmns_cmb_ee_delta,
                "sigma8": pmns_cmb_sigma8_delta,
                "max_abs": pmns_cmb_delta_max
            },
            "baryogenesis_coupling_pass": pmns_baryo_coupling_pass,
            "cmb_response_pass": pmns_cmb_response_pass,
            "pass": pmns_propagation_pass
        },
        "proton": {
            "rel_error": proton_rel_err,
            "pass": proton_pass
        },
        "alpha_web_ci": {
            "pass": alpha_pass
        },
        "chiral_symmetry_breaking": {
            "pass": chiral_pass
        },
        "cmb": {
            "tt_full_red": tt,
            "te_full_red": te,
            "ee_full_red": ee,
            "pass": cmb_pass
        },
        "sigma8": {
            "from_cmb": sigma8_cmb,
            "from_decomposition": sigma8_dec,
            "target": sigma8_target,
            "match_abs": sigma8_match,
            "pass": sigma8_pass
        },
        "overall_pass": overall_pass
    });

    let mut json_file = File::create(&json_path).expect("create json");
    writeln!(
        json_file,
        "{}",
        serde_json::to_string_pretty(&report).expect("serialize json")
    )
    .ok();

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("global_gate overall_pass={}", overall_pass);

    if std::env::var("GUTOE_GLOBAL_GATE_STRICT").ok().as_deref() == Some("1") && !overall_pass {
        std::process::exit(3);
    }
}
