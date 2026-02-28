//! Quantitative Big Bang Nucleosynthesis report for GRAND-349.

use gutoe_physics::{
    evaluate_bbn_gate, lithium7_be7_component_dark_coupled, lithium7_be7_component_unsuppressed,
    lithium7_direct_component, lithium7_reaction_network_source, BbnWindows, DEUTERIUM_ETA_EXP,
    ETA10_REF, HELIUM3_ETA_EXP, LITHIUM7_BE7_DARK_SUPPRESSION, LITHIUM7_CHANNEL_COUPLED_FACTOR,
    LITHIUM7_REACTION_NETWORK_GAIN, LITHIUM7_TENSION_AMPLIFICATION, LITHIUM7_VISIBLE_FRACTION,
};
use std::fs::{self, File};
use std::io::Write;

fn main() {
    let windows = BbnWindows::default();
    let score = evaluate_bbn_gate(windows);
    let li7_source = lithium7_reaction_network_source(score.eta10);
    let li7_direct = lithium7_direct_component(score.eta10);
    let li7_be7_raw = lithium7_be7_component_unsuppressed(score.eta10);
    let li7_be7_dark = lithium7_be7_component_dark_coupled(score.eta10);

    let out_dir = std::env::var("GUTOE_BBN_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/bbn_report.txt");
    let json_path = format!("{out_dir}/bbn_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[bbn_structural]").expect("write");
    writeln!(txt, "eta10 = {:.12}", score.eta10).expect("write");
    writeln!(txt, "eta10_ref = {:.12}", ETA10_REF).expect("write");
    writeln!(txt, "deuterium_eta_exp = {:.12}", DEUTERIUM_ETA_EXP).expect("write");
    writeln!(txt, "helium3_eta_exp = {:.12}", HELIUM3_ETA_EXP).expect("write");
    writeln!(
        txt,
        "lithium7_tension_amplification = {:.12}",
        LITHIUM7_TENSION_AMPLIFICATION
    )
    .expect("write");
    writeln!(txt, "lithium7_visible_fraction = {:.12}", LITHIUM7_VISIBLE_FRACTION).expect("write");
    writeln!(
        txt,
        "lithium7_reaction_network_gain = {:.12}",
        LITHIUM7_REACTION_NETWORK_GAIN
    )
    .expect("write");
    writeln!(txt, "Y_p_pred = {:.12e}", score.yp_pred).expect("write");
    writeln!(txt, "D/H_pred = {:.12e}", score.dh_pred).expect("write");
    writeln!(txt, "3He/H_pred = {:.12e}", score.he3h_pred).expect("write");
    writeln!(txt, "7Li/H_pred = {:.12e}", score.li7h_pred).expect("write");
    writeln!(txt, "Li7_source = {:.12e}", li7_source).expect("write");
    writeln!(txt, "Li7_direct_component = {:.12e}", li7_direct).expect("write");
    writeln!(txt, "Li7_be7_component_raw = {:.12e}", li7_be7_raw).expect("write");
    writeln!(txt, "Li7_be7_component_dark = {:.12e}", li7_be7_dark).expect("write");
    writeln!(txt, "Li7_be7_dark_suppression = {:.12}", LITHIUM7_BE7_DARK_SUPPRESSION).expect("write");
    writeln!(txt, "Li7_channel_coupled_factor = {:.12}", LITHIUM7_CHANNEL_COUPLED_FACTOR).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[residuals]").expect("write");
    writeln!(txt, "Y_p_delta = {:.12e}", score.yp_delta).expect("write");
    writeln!(txt, "D/H_rel_error = {:.12}", score.dh_rel_error).expect("write");
    writeln!(txt, "3He/H_rel_error = {:.12}", score.he3_rel_error).expect("write");
    writeln!(txt, "Li7_tension_ratio = {:.12}", score.li_tension_ratio).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[gate]").expect("write");
    writeln!(txt, "yp_ok = {}", score.yp_ok).expect("write");
    writeln!(txt, "dh_ok = {}", score.dh_ok).expect("write");
    writeln!(txt, "he3_ok = {}", score.he3_ok).expect("write");
    writeln!(txt, "li_tension_ok = {}", score.li_tension_ok).expect("write");
    writeln!(txt, "passes_primary = {}", score.passes_primary()).expect("write");
    writeln!(txt, "passes_all = {}", score.passes_all()).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"eta10\": {:.12},\n  \"eta10_ref\": {:.12},\n  \"deuterium_eta_exp\": {:.12},\n  \"helium3_eta_exp\": {:.12},\n  \"lithium7_tension_amplification\": {:.12},\n  \"lithium7_visible_fraction\": {:.12},\n  \"lithium7_reaction_network_gain\": {:.12},\n  \"lithium7_be7_dark_suppression\": {:.12},\n  \"lithium7_channel_coupled_factor\": {:.12},\n  \"yp_pred\": {:.12e},\n  \"dh_pred\": {:.12e},\n  \"he3h_pred\": {:.12e},\n  \"li7h_pred\": {:.12e},\n  \"li7_source\": {:.12e},\n  \"li7_direct_component\": {:.12e},\n  \"li7_be7_component_raw\": {:.12e},\n  \"li7_be7_component_dark\": {:.12e},\n  \"yp_delta\": {:.12e},\n  \"dh_rel_error\": {:.12},\n  \"he3_rel_error\": {:.12},\n  \"li_tension_ratio\": {:.12},\n  \"windows\": {{\"yp_abs_max\": {:.12}, \"dh_rel_max\": {:.12}, \"he3_rel_max\": {:.12}, \"li_tension_ratio_min\": {:.12}, \"li_tension_ratio_max\": {:.12}}},\n  \"yp_ok\": {},\n  \"dh_ok\": {},\n  \"he3_ok\": {},\n  \"li_tension_ok\": {},\n  \"passes_primary\": {},\n  \"passes_all\": {}\n}}",
        score.eta10,
        ETA10_REF,
        DEUTERIUM_ETA_EXP,
        HELIUM3_ETA_EXP,
        LITHIUM7_TENSION_AMPLIFICATION,
        LITHIUM7_VISIBLE_FRACTION,
        LITHIUM7_REACTION_NETWORK_GAIN,
        LITHIUM7_BE7_DARK_SUPPRESSION,
        LITHIUM7_CHANNEL_COUPLED_FACTOR,
        score.yp_pred,
        score.dh_pred,
        score.he3h_pred,
        score.li7h_pred,
        li7_source,
        li7_direct,
        li7_be7_raw,
        li7_be7_dark,
        score.yp_delta,
        score.dh_rel_error,
        score.he3_rel_error,
        score.li_tension_ratio,
        windows.yp_abs_max,
        windows.dh_rel_max,
        windows.he3_rel_max,
        windows.li_tension_ratio_min,
        windows.li_tension_ratio_max,
        score.yp_ok,
        score.dh_ok,
        score.he3_ok,
        score.li_tension_ok,
        score.passes_primary(),
        score.passes_all()
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "BBN: η10={:.4}, Yp={:.6}, D/H={:.4e}, 3He/H={:.4e}, 7Li/H={:.4e}, primary_pass={}, all_pass={}",
        score.eta10,
        score.yp_pred,
        score.dh_pred,
        score.he3h_pred,
        score.li7h_pred,
        score.passes_primary(),
        score.passes_all()
    );
}
