//! Quantitative baryogenesis report for GRAND-348.

use gutoe_physics::{
    baryogenesis_structural_prefactor, evaluate_baryogenesis_gate, nonequilibrium_survival_factor,
    BaryogenesisWindows,
};
use std::fs::{self, File};
use std::io::Write;

fn main() {
    let windows = BaryogenesisWindows::default();
    let score = evaluate_baryogenesis_gate(windows);
    let pref = baryogenesis_structural_prefactor();
    let f_neq = nonequilibrium_survival_factor();

    let out_dir =
        std::env::var("GUTOE_BARYO_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/baryogenesis_report.txt");
    let json_path = format!("{out_dir}/baryogenesis_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[baryogenesis_structural]").expect("write");
    writeln!(
        txt,
        "jarlskog_ckm_direct = {:.12e}",
        score.jarlskog_ckm_direct
    )
    .expect("write");
    writeln!(
        txt,
        "jarlskog_ckm_texture = {:.12e}",
        score.jarlskog_ckm_texture
    )
    .expect("write");
    writeln!(
        txt,
        "pmns_theta23_alpha2_c = {:.12e}",
        score.pmns_theta23_alpha2_c
    )
    .expect("write");
    writeln!(
        txt,
        "leptogenesis_pmns_gain = {:.12e}",
        score.leptogenesis_pmns_gain
    )
    .expect("write");
    writeln!(
        txt,
        "leptogenesis_pmns_scalar = {:.12e}",
        score.leptogenesis_pmns_scalar
    )
    .expect("write");
    writeln!(
        txt,
        "leptogenesis_multiplier = {:.12e}",
        score.leptogenesis_multiplier
    )
    .expect("write");
    writeln!(txt, "prefactor = {:.12e}", pref).expect("write");
    writeln!(txt, "nonequilibrium_survival = {:.12e}", f_neq).expect("write");
    writeln!(txt, "eta_predicted = {:.12e}", score.eta_predicted).expect("write");
    writeln!(txt, "eta_observed = {:.12e}", score.eta_observed).expect("write");
    writeln!(txt, "eta_rel_error = {:.12}", score.eta_rel_error).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[sakharov_checks]").expect("write");
    writeln!(txt, "cp_violation_ok = {}", score.cp_violation_ok).expect("write");
    writeln!(
        txt,
        "baryon_violation_channel_ok = {}",
        score.baryon_violation_channel_ok
    )
    .expect("write");
    writeln!(txt, "nonequilibrium_ok = {}", score.nonequilibrium_ok).expect("write");
    writeln!(txt, "eta_window_ok = {}", score.eta_window_ok).expect("write");
    writeln!(txt, "sakharov_ok = {}", score.sakharov_ok()).expect("write");
    writeln!(txt, "passes_all = {}", score.passes_all()).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"jarlskog_ckm_direct\": {:.12e},\n  \"jarlskog_ckm_texture\": {:.12e},\n  \"pmns_theta23_alpha2_c\": {:.12e},\n  \"leptogenesis_pmns_gain\": {:.12e},\n  \"leptogenesis_pmns_scalar\": {:.12e},\n  \"leptogenesis_multiplier\": {:.12e},\n  \"prefactor\": {:.12e},\n  \"nonequilibrium_survival\": {:.12e},\n  \"eta_predicted\": {:.12e},\n  \"eta_observed\": {:.12e},\n  \"eta_rel_error\": {:.12},\n  \"windows\": {{\"eta_rel_error_max\": {:.12}}},\n  \"cp_violation_ok\": {},\n  \"baryon_violation_channel_ok\": {},\n  \"nonequilibrium_ok\": {},\n  \"eta_window_ok\": {},\n  \"sakharov_ok\": {},\n  \"passes_all\": {}\n}}",
        score.jarlskog_ckm_direct,
        score.jarlskog_ckm_texture,
        score.pmns_theta23_alpha2_c,
        score.leptogenesis_pmns_gain,
        score.leptogenesis_pmns_scalar,
        score.leptogenesis_multiplier,
        pref,
        f_neq,
        score.eta_predicted,
        score.eta_observed,
        score.eta_rel_error,
        windows.eta_rel_error_max,
        score.cp_violation_ok,
        score.baryon_violation_channel_ok,
        score.nonequilibrium_ok,
        score.eta_window_ok,
        score.sakharov_ok(),
        score.passes_all()
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "η_B(pred)={:.6e}, η_B(obs)={:.6e}, rel_err={:.4}, passes_all={}",
        score.eta_predicted,
        score.eta_observed,
        score.eta_rel_error,
        score.passes_all()
    );
}
