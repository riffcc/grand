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
    let status = command.status().map_err(|e| format!("spawn {cmd}: {e}"))?;
    if !status.success() {
        return Err(format!(
            "command failed: {} {} (status={status})",
            cmd,
            args.join(" ")
        ));
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

fn v_str(v: &Value, path: &[&str]) -> Result<String, String> {
    let mut cur = v;
    for k in path {
        cur = cur
            .get(*k)
            .ok_or_else(|| format!("missing key {} in {}", k, path.join(".")))?;
    }
    cur.as_str()
        .map(str::to_string)
        .ok_or_else(|| format!("non-string at {}", path.join(".")))
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
    if let Err(e) = run(
        "cargo",
        &["run", "-q", "-p", "gutoe-em", "--bin", "flavor_ci_gate"],
    ) {
        eprintln!("global_gate: flavor_ci_gate failed: {e}");
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
            "proton_mass_report",
        ],
    ) {
        eprintln!("global_gate: proton_mass_report failed: {e}");
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
            "alpha_web_ci_report",
        ],
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
    if let Err(e) = run(
        "cargo",
        &["run", "-q", "-p", "gutoe-em", "--bin", "neutrino_ci_gate"],
    ) {
        eprintln!("global_gate: neutrino_ci_gate failed: {e}");
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
            "lithium7_stellar_ci_gate",
        ],
    ) {
        eprintln!("global_gate: lithium7_stellar_ci_gate failed: {e}");
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
            "parameter_degeneracy_ci_gate",
        ],
    ) {
        eprintln!("global_gate: parameter_degeneracy_ci_gate failed: {e}");
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
            "abiogenesis_ci_gate",
        ],
    ) {
        eprintln!("global_gate: abiogenesis_ci_gate failed: {e}");
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
            "entropy_progression_ci_gate",
        ],
    ) {
        eprintln!("global_gate: entropy_progression_ci_gate failed: {e}");
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
            "cardiovascular_binding_ci_gate",
        ],
    ) {
        eprintln!("global_gate: cardiovascular_binding_ci_gate failed: {e}");
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
            "ms_localized_dual_compartment_ci_gate",
        ],
    ) {
        eprintln!("global_gate: ms_localized_dual_compartment_ci_gate failed: {e}");
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
            "antibiotic_resistance_ci_gate",
        ],
    ) {
        eprintln!("global_gate: antibiotic_resistance_ci_gate failed: {e}");
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
            "phage_host_matching_ci_gate",
        ],
    ) {
        eprintln!("global_gate: phage_host_matching_ci_gate failed: {e}");
        std::process::exit(2);
    }

    run_or_exit(
        "baryogenesis_report[gain=0]",
        "cargo",
        &[
            "run",
            "-q",
            "-p",
            "gutoe-physics",
            "--bin",
            "baryogenesis_report",
        ],
        &[(PMNS_GAIN_ENV, PMNS_GAIN_BASELINE)],
    );
    let baryo0 =
        read_json("/tmp/bh_renders/baryogenesis_report.json").expect("baryogenesis gain0 json");

    run_or_exit(
        "baryogenesis_report[gain=1]",
        "cargo",
        &[
            "run",
            "-q",
            "-p",
            "gutoe-physics",
            "--bin",
            "baryogenesis_report",
        ],
        &[(PMNS_GAIN_ENV, PMNS_GAIN_STRUCTURAL)],
    );
    let baryo1 =
        read_json("/tmp/bh_renders/baryogenesis_report.json").expect("baryogenesis gain1 json");

    run_or_exit(
        "cmb_full_derived_report[gain=0]",
        "cargo",
        &[
            "run",
            "-q",
            "-p",
            "gutoe-physics",
            "--bin",
            "cmb_full_derived_report",
        ],
        &[(PMNS_GAIN_ENV, PMNS_GAIN_BASELINE)],
    );
    let cmb0 = read_json("/tmp/bh_renders/cmb_full_derived/cmb_full_derived_report.json")
        .expect("cmb gain0 json");

    run_or_exit(
        "cmb_full_derived_report[gain=1]",
        "cargo",
        &[
            "run",
            "-q",
            "-p",
            "gutoe-physics",
            "--bin",
            "cmb_full_derived_report",
        ],
        &[(PMNS_GAIN_ENV, PMNS_GAIN_STRUCTURAL)],
    );
    let cmb = read_json("/tmp/bh_renders/cmb_full_derived/cmb_full_derived_report.json")
        .expect("cmb gain1 json");

    run_or_exit(
        "sigma8_decomposition[gain=1]",
        "cargo",
        &[
            "run",
            "-q",
            "-p",
            "gutoe-physics",
            "--bin",
            "sigma8_decomposition",
        ],
        &[(PMNS_GAIN_ENV, PMNS_GAIN_STRUCTURAL)],
    );

    let flavor = read_json("/tmp/bh_renders/flavor_ci_gate.json").expect("flavor json");
    let proton = read_json("/tmp/bh_renders/proton_mass_report/proton_mass_report.json")
        .expect("proton json");
    let alpha = read_json("/tmp/bh_renders/alpha_web_ci_report/alpha_web_ci_report.json")
        .expect("alpha json");
    let chiral =
        read_json("/tmp/bh_renders/chiral_symmetry_breaking/chiral_symmetry_breaking_report.json")
            .expect("chiral json");
    let neutrino = read_json("/tmp/bh_renders/neutrino_ci_gate.json").expect("neutrino json");
    let li7_stellar =
        read_json("/tmp/bh_renders/lithium7_stellar_ci_gate.json").expect("lithium7 stellar json");
    let degeneracy = read_json("/tmp/bh_renders/parameter_degeneracy_ci_gate.json")
        .expect("parameter degeneracy json");
    let abiogenesis =
        read_json("/tmp/bh_renders/abiogenesis_ci_gate.json").expect("abiogenesis json");
    let entropy_progression = read_json("/tmp/bh_renders/entropy_progression_ci_gate.json")
        .expect("entropy progression json");
    let ms_localized = read_json("/tmp/bh_renders/ms_localized_dual_compartment_ci_gate.json")
        .expect("ms localized dual-compartment gate json");
    let antibiotic_resistance = read_json("/tmp/bh_renders/antibiotic_resistance_ci_gate.json")
        .expect("antibiotic resistance gate json");
    let phage_host_matching = read_json("/tmp/bh_renders/phage_host_matching_ci_gate.json")
        .expect("phage host matching gate json");
    let sigma = read_json("/tmp/bh_renders/sigma8_decomposition/sigma8_decomposition_report.json")
        .expect("sigma json");

    let pmns_corr_res = v_f64(
        &flavor,
        &["pmns_theta23_improvement", "corrected_abs_residual_deg"],
    )
    .unwrap();
    let pmns_corr_pass =
        v_bool(&flavor, &["pmns_theta23_improvement", "pass"]).unwrap() && pmns_corr_res <= 0.01;

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
    let neutrino_pass = v_bool(&neutrino, &["overall_pass"]).unwrap();
    let neutrino_hierarchy = v_str(&neutrino, &["texture_lane", "hierarchy_prediction"]).unwrap();
    let neutrino_mass_character =
        v_str(&neutrino, &["texture_lane", "mass_character_prediction"]).unwrap();
    let neutrino_m3_ev = v_f64(&neutrino, &["absolute_lane", "m3_ev"]).unwrap();
    let neutrino_sum_ev = v_f64(&neutrino, &["absolute_lane", "sum_ev"]).unwrap();
    let neutrino_majorana_resid =
        v_f64(&neutrino, &["texture_lane", "majorana_symmetry_residual"]).unwrap();
    let li7_stellar_pass = v_bool(&li7_stellar, &["overall_pass"]).unwrap();
    let li7_stellar_delta_abs = v_f64(&li7_stellar, &["best_match", "closure_delta_abs"]).unwrap();
    let li7_stellar_delta_abs_max = v_f64(&li7_stellar, &["closure_delta_abs_max"]).unwrap();
    let degeneracy_pass = v_bool(&degeneracy, &["overall_pass"]).unwrap();
    let degeneracy_free = v_f64(&degeneracy, &["counts", "free_parameters"]).unwrap();
    let degeneracy_rank_tunable =
        v_f64(&degeneracy, &["linear_algebra", "tunable_only", "rank"]).unwrap();
    let degeneracy_transfer_coupling_max = v_f64(
        &degeneracy,
        &[
            "hidden_reencoding_checks",
            "tunable_to_transfer_max_abs_sensitivity",
        ],
    )
    .unwrap();
    let degeneracy_verdict = v_str(&degeneracy, &["verdict"]).unwrap();
    let abiogenesis_pass = v_bool(&abiogenesis, &["overall_pass"]).unwrap();
    let abiogenesis_n_times_p = v_f64(&abiogenesis, &["closure", "n_times_p"]).unwrap();
    let abiogenesis_lower_3sigma =
        v_f64(&abiogenesis, &["inevitability", "n_times_p_lower_3sigma"]).unwrap();
    let abiogenesis_margin = v_f64(&abiogenesis, &["inevitability", "robust_margin"]).unwrap();
    let entropy_progression_pass = v_bool(&entropy_progression, &["overall_pass"]).unwrap();
    let entropy_progression_final_per_area = v_f64(
        &entropy_progression,
        &["summary", "final_total_per_area_w_m2_k"],
    )
    .unwrap();
    let entropy_progression_final_universe = v_f64(
        &entropy_progression,
        &["summary", "final_total_universe_w_k"],
    )
    .unwrap();
    let entropy_progression_maxima =
        v_f64(&entropy_progression, &["summary", "local_maxima_count"]).unwrap();
    let entropy_progression_minima =
        v_f64(&entropy_progression, &["summary", "local_minima_count"]).unwrap();
    let ms_localized_pass = v_bool(&ms_localized, &["overall_pass"]).unwrap();
    let ms_localized_efficacy_pass = v_bool(&ms_localized, &["gate", "efficacy_pass"]).unwrap();
    let ms_localized_safety_pass = v_bool(&ms_localized, &["gate", "safety_pass"]).unwrap();
    let ms_localized_arr_reduction_2y =
        v_f64(&ms_localized, &["score", "arr_reduction_2y"]).unwrap();
    let ms_localized_lesion_reduction_10y =
        v_f64(&ms_localized, &["score", "lesion_reduction_10y"]).unwrap();
    let ms_localized_prob_above_renal_high =
        v_f64(&ms_localized, &["score", "prob_above_renal_high"]).unwrap();
    let ms_localized_prob_in_target_zone =
        v_f64(&ms_localized, &["score", "prob_in_target_zone"]).unwrap();
    let ms_localized_localization_factor =
        v_f64(&ms_localized, &["controls", "localization_factor"]).unwrap();
    let ms_localized_transduction_efficiency =
        v_f64(&ms_localized, &["controls", "transduction_efficiency"]).unwrap();
    let antibiotic_resistance_pass = v_bool(&antibiotic_resistance, &["overall_pass"]).unwrap();
    let antibiotic_resistance_pair_count =
        v_f64(&antibiotic_resistance, &["summary", "pair_count"]).unwrap();
    let antibiotic_resistance_mean_abs_log10_error = v_f64(
        &antibiotic_resistance,
        &["summary", "mean_abs_log10_error_pred_vs_anchor"],
    )
    .unwrap();
    let antibiotic_resistance_ndm_occ_1um = v_f64(
        &antibiotic_resistance,
        &["summary", "ndm_max_predicted_occupancy_at_1uM"],
    )
    .unwrap();
    let antibiotic_resistance_tem_winner =
        v_str(&antibiotic_resistance, &["summary", "tem_predicted_winner"]).unwrap();
    let antibiotic_resistance_kpc_winner =
        v_str(&antibiotic_resistance, &["summary", "kpc_predicted_winner"]).unwrap();
    let phage_host_matching_pass = v_bool(&phage_host_matching, &["overall_pass"]).unwrap();
    let phage_host_matching_pair_count =
        v_f64(&phage_host_matching, &["summary", "pair_count"]).unwrap();
    let phage_host_matching_mean_best_lysis =
        v_f64(&phage_host_matching, &["summary", "mean_best_lysis_score"]).unwrap();
    let phage_host_matching_probe_delta = v_f64(
        &phage_host_matching,
        &["summary", "resistance_independence_probe_abs_delta"],
    )
    .unwrap();
    let phage_host_matching_ndm_best =
        v_str(&phage_host_matching, &["summary", "ndm_best_phage"]).unwrap();

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
        && neutrino_pass
        && li7_stellar_pass
        && degeneracy_pass
        && abiogenesis_pass
        && entropy_progression_pass
        && ms_localized_pass
        && antibiotic_resistance_pass
        && phage_host_matching_pass
        && cmb_pass
        && sigma8_pass;

    let txt_path = Path::new(&out_dir).join("global_gate_report.txt");
    let json_path = Path::new(&out_dir).join("global_gate_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[global_gate]").ok();
    writeln!(
        txt,
        "pmns_theta23_corrected_residual_deg = {:.12}",
        pmns_corr_res
    )
    .ok();
    writeln!(txt, "pmns_theta23_pass = {}", pmns_corr_pass).ok();
    writeln!(txt, "pmns_gain0_eta_predicted = {:.12e}", eta0).ok();
    writeln!(txt, "pmns_gain1_eta_predicted = {:.12e}", eta1).ok();
    writeln!(
        txt,
        "pmns_baryogenesis_eta_delta = {:.12e}",
        pmns_baryo_eta_delta
    )
    .ok();
    writeln!(
        txt,
        "pmns_baryogenesis_coupling_pass = {}",
        pmns_baryo_coupling_pass
    )
    .ok();
    writeln!(txt, "pmns_cmb_tt_delta = {:.12e}", pmns_cmb_tt_delta).ok();
    writeln!(txt, "pmns_cmb_te_delta = {:.12e}", pmns_cmb_te_delta).ok();
    writeln!(txt, "pmns_cmb_ee_delta = {:.12e}", pmns_cmb_ee_delta).ok();
    writeln!(
        txt,
        "pmns_cmb_sigma8_delta = {:.12e}",
        pmns_cmb_sigma8_delta
    )
    .ok();
    writeln!(txt, "pmns_cmb_delta_max_abs = {:.12e}", pmns_cmb_delta_max).ok();
    writeln!(txt, "pmns_cmb_response_pass = {}", pmns_cmb_response_pass).ok();
    writeln!(txt, "pmns_propagation_pass = {}", pmns_propagation_pass).ok();
    writeln!(txt, "proton_rel_error = {:.12e}", proton_rel_err).ok();
    writeln!(txt, "proton_pass = {}", proton_pass).ok();
    writeln!(txt, "alpha_web_ci_pass = {}", alpha_pass).ok();
    writeln!(txt, "chiral_symmetry_breaking_pass = {}", chiral_pass).ok();
    writeln!(txt, "neutrino_hierarchy = {}", neutrino_hierarchy).ok();
    writeln!(txt, "neutrino_mass_character = {}", neutrino_mass_character).ok();
    writeln!(
        txt,
        "neutrino_majorana_symmetry_residual = {:.12e}",
        neutrino_majorana_resid
    )
    .ok();
    writeln!(txt, "neutrino_m3_ev = {:.12e}", neutrino_m3_ev).ok();
    writeln!(txt, "neutrino_sum_ev = {:.12e}", neutrino_sum_ev).ok();
    writeln!(txt, "neutrino_pass = {}", neutrino_pass).ok();
    writeln!(
        txt,
        "li7_stellar_closure_delta_abs = {:.12}",
        li7_stellar_delta_abs
    )
    .ok();
    writeln!(
        txt,
        "li7_stellar_closure_delta_abs_max = {:.12}",
        li7_stellar_delta_abs_max
    )
    .ok();
    writeln!(txt, "li7_stellar_pass = {}", li7_stellar_pass).ok();
    writeln!(txt, "degeneracy_verdict = {}", degeneracy_verdict).ok();
    writeln!(txt, "degeneracy_free_parameters = {:.0}", degeneracy_free).ok();
    writeln!(
        txt,
        "degeneracy_rank_tunable = {:.0}",
        degeneracy_rank_tunable
    )
    .ok();
    writeln!(
        txt,
        "degeneracy_transfer_coupling_max = {:.12e}",
        degeneracy_transfer_coupling_max
    )
    .ok();
    writeln!(txt, "degeneracy_pass = {}", degeneracy_pass).ok();
    writeln!(txt, "abiogenesis_n_times_p = {:.12}", abiogenesis_n_times_p).ok();
    writeln!(
        txt,
        "abiogenesis_n_times_p_lower_3sigma = {:.12}",
        abiogenesis_lower_3sigma
    )
    .ok();
    writeln!(
        txt,
        "abiogenesis_robust_margin = {:.12}",
        abiogenesis_margin
    )
    .ok();
    writeln!(txt, "abiogenesis_pass = {}", abiogenesis_pass).ok();
    writeln!(
        txt,
        "entropy_progression_final_per_area_w_m2_k = {:.12e}",
        entropy_progression_final_per_area
    )
    .ok();
    writeln!(
        txt,
        "entropy_progression_final_universe_w_k = {:.12e}",
        entropy_progression_final_universe
    )
    .ok();
    writeln!(
        txt,
        "entropy_progression_local_maxima = {:.0}",
        entropy_progression_maxima
    )
    .ok();
    writeln!(
        txt,
        "entropy_progression_local_minima = {:.0}",
        entropy_progression_minima
    )
    .ok();
    writeln!(
        txt,
        "entropy_progression_pass = {}",
        entropy_progression_pass
    )
    .ok();
    writeln!(
        txt,
        "ms_localized_localization_factor = {:.6}",
        ms_localized_localization_factor
    )
    .ok();
    writeln!(
        txt,
        "ms_localized_transduction_efficiency = {:.6}",
        ms_localized_transduction_efficiency
    )
    .ok();
    writeln!(
        txt,
        "ms_localized_arr_reduction_2y = {:.12}",
        ms_localized_arr_reduction_2y
    )
    .ok();
    writeln!(
        txt,
        "ms_localized_lesion_reduction_10y = {:.12}",
        ms_localized_lesion_reduction_10y
    )
    .ok();
    writeln!(
        txt,
        "ms_localized_prob_above_renal_high = {:.12e}",
        ms_localized_prob_above_renal_high
    )
    .ok();
    writeln!(
        txt,
        "ms_localized_prob_in_target_zone = {:.12}",
        ms_localized_prob_in_target_zone
    )
    .ok();
    writeln!(
        txt,
        "ms_localized_efficacy_pass = {}",
        ms_localized_efficacy_pass
    )
    .ok();
    writeln!(txt, "ms_localized_safety_pass = {}", ms_localized_safety_pass).ok();
    writeln!(txt, "ms_localized_pass = {}", ms_localized_pass).ok();
    writeln!(
        txt,
        "antibiotic_resistance_pair_count = {:.0}",
        antibiotic_resistance_pair_count
    )
    .ok();
    writeln!(
        txt,
        "antibiotic_resistance_mean_abs_log10_error = {:.12}",
        antibiotic_resistance_mean_abs_log10_error
    )
    .ok();
    writeln!(
        txt,
        "antibiotic_resistance_ndm_max_occ_1uM = {:.12}",
        antibiotic_resistance_ndm_occ_1um
    )
    .ok();
    writeln!(
        txt,
        "antibiotic_resistance_tem_predicted_winner = {}",
        antibiotic_resistance_tem_winner
    )
    .ok();
    writeln!(
        txt,
        "antibiotic_resistance_kpc_predicted_winner = {}",
        antibiotic_resistance_kpc_winner
    )
    .ok();
    writeln!(
        txt,
        "antibiotic_resistance_pass = {}",
        antibiotic_resistance_pass
    )
    .ok();
    writeln!(
        txt,
        "phage_host_matching_pair_count = {:.0}",
        phage_host_matching_pair_count
    )
    .ok();
    writeln!(
        txt,
        "phage_host_matching_mean_best_lysis_score = {:.12}",
        phage_host_matching_mean_best_lysis
    )
    .ok();
    writeln!(
        txt,
        "phage_host_matching_probe_abs_delta = {:.12e}",
        phage_host_matching_probe_delta
    )
    .ok();
    writeln!(
        txt,
        "phage_host_matching_ndm_best_phage = {}",
        phage_host_matching_ndm_best
    )
    .ok();
    writeln!(txt, "phage_host_matching_pass = {}", phage_host_matching_pass).ok();
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
        "neutrino": {
            "hierarchy_prediction": neutrino_hierarchy,
            "mass_character_prediction": neutrino_mass_character,
            "majorana_symmetry_residual": neutrino_majorana_resid,
            "m3_ev": neutrino_m3_ev,
            "sum_ev": neutrino_sum_ev,
            "pass": neutrino_pass
        },
        "lithium7_stellar_depletion": {
            "closure_delta_abs": li7_stellar_delta_abs,
            "closure_delta_abs_max": li7_stellar_delta_abs_max,
            "pass": li7_stellar_pass
        },
        "parameter_degeneracy_audit": {
            "verdict": degeneracy_verdict,
            "free_parameters": degeneracy_free,
            "rank_tunable": degeneracy_rank_tunable,
            "transfer_coupling_max": degeneracy_transfer_coupling_max,
            "pass": degeneracy_pass
        },
        "abiogenesis": {
            "n_times_p": abiogenesis_n_times_p,
            "n_times_p_lower_3sigma": abiogenesis_lower_3sigma,
            "robust_margin": abiogenesis_margin,
            "pass": abiogenesis_pass
        },
        "entropy_progression": {
            "final_per_area_w_m2_k": entropy_progression_final_per_area,
            "final_universe_w_k": entropy_progression_final_universe,
            "local_maxima_count": entropy_progression_maxima,
            "local_minima_count": entropy_progression_minima,
            "pass": entropy_progression_pass
        },
        "ms_localized_dual_compartment": {
            "localization_factor": ms_localized_localization_factor,
            "transduction_efficiency": ms_localized_transduction_efficiency,
            "arr_reduction_2y": ms_localized_arr_reduction_2y,
            "lesion_reduction_10y": ms_localized_lesion_reduction_10y,
            "prob_above_renal_high": ms_localized_prob_above_renal_high,
            "prob_in_target_zone": ms_localized_prob_in_target_zone,
            "efficacy_pass": ms_localized_efficacy_pass,
            "safety_pass": ms_localized_safety_pass,
            "pass": ms_localized_pass
        },
        "antibiotic_resistance": {
            "pair_count": antibiotic_resistance_pair_count,
            "mean_abs_log10_error_pred_vs_anchor": antibiotic_resistance_mean_abs_log10_error,
            "ndm_max_predicted_occupancy_at_1uM": antibiotic_resistance_ndm_occ_1um,
            "tem_predicted_winner": antibiotic_resistance_tem_winner,
            "kpc_predicted_winner": antibiotic_resistance_kpc_winner,
            "pass": antibiotic_resistance_pass
        },
        "phage_host_matching": {
            "pair_count": phage_host_matching_pair_count,
            "mean_best_lysis_score": phage_host_matching_mean_best_lysis,
            "resistance_independence_probe_abs_delta": phage_host_matching_probe_delta,
            "ndm_best_phage": phage_host_matching_ndm_best,
            "pass": phage_host_matching_pass
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
