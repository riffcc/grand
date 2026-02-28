//! Lithium-7 report: observed-anchored diagnostic vs channel-derived predictive lane.

use gutoe_physics::bbn::{
    evaluate_bbn_gate, lithium7_be7_component_dark_coupled, lithium7_be7_component_unsuppressed,
    lithium7_direct_component, lithium7_reaction_network_source, lithium7_tension_ratio_channel_coupled,
    primordial_lithium7_ratio_channel_coupled, primordial_lithium7_ratio_observed_anchored,
    BbnWindows, LI7H_OBSERVED, LITHIUM7_BE7_CHANNEL_FRACTION, LITHIUM7_BE7_DARK_SUPPRESSION,
    LITHIUM7_CHANNEL_COUPLED_FACTOR, LITHIUM7_DIRECT_CHANNEL_FRACTION, LITHIUM7_REACTION_NETWORK_GAIN,
    LITHIUM7_TENSION_AMPLIFICATION, LITHIUM7_VISIBLE_FRACTION, LITHIUM7_VOID_CORRECTION,
};
use gutoe_physics::constants::DARK_FRACTION_TOTAL_STATE_SPLIT;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let out_dir =
        std::env::var("GUTOE_LITHIUM_OUT").unwrap_or_else(|_| "/tmp/bh_renders/lithium7_report".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let score = evaluate_bbn_gate(BbnWindows::default());
    let li7h_predictive = score.li7h_pred;
    let ratio_predictive = score.li_tension_ratio;
    let li7h_observed_anchor = primordial_lithium7_ratio_observed_anchored(score.eta10);
    let ratio_observed_anchor = li7h_observed_anchor / LI7H_OBSERVED;

    let li7h_channel = primordial_lithium7_ratio_channel_coupled(score.eta10);
    let ratio_channel = lithium7_tension_ratio_channel_coupled(score.eta10);
    let li7_source = lithium7_reaction_network_source(score.eta10);
    let li7_direct = lithium7_direct_component(score.eta10);
    let li7_be7_raw = lithium7_be7_component_unsuppressed(score.eta10);
    let li7_be7_dark = lithium7_be7_component_dark_coupled(score.eta10);
    let closure_factor =
        LITHIUM7_TENSION_AMPLIFICATION * DARK_FRACTION_TOTAL_STATE_SPLIT * LITHIUM7_VOID_CORRECTION;

    let txt_path = out.join("lithium7_report.txt");
    let json_path = out.join("lithium7_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "eta10 = {:.12}", score.eta10).expect("write");
    writeln!(txt, "li7h_observed = {:.12e}", LI7H_OBSERVED).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[observed_anchored_diagnostic]").expect("write");
    writeln!(
        txt,
        "amplification = {:.12} (expected tension lane)",
        LITHIUM7_TENSION_AMPLIFICATION
    )
    .expect("write");
    writeln!(txt, "li7h_pred = {:.12e}", li7h_observed_anchor).expect("write");
    writeln!(txt, "li7_tension_ratio = {:.12}", ratio_observed_anchor).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[channel_derived_predictive]").expect("write");
    writeln!(txt, "li7_visible_fraction_11_over_16 = {:.12}", LITHIUM7_VISIBLE_FRACTION).expect("write");
    writeln!(txt, "li7_reaction_network_gain_33_over_16 = {:.12}", LITHIUM7_REACTION_NETWORK_GAIN)
        .expect("write");
    writeln!(txt, "li7_source = {:.12e}", li7_source).expect("write");
    writeln!(txt, "li7_direct_component = {:.12e}", li7_direct).expect("write");
    writeln!(txt, "li7_be7_component_raw = {:.12e}", li7_be7_raw).expect("write");
    writeln!(txt, "li7_be7_component_dark = {:.12e}", li7_be7_dark).expect("write");
    writeln!(txt, "occupancy_factor_5_over_16 = {:.12}", DARK_FRACTION_TOTAL_STATE_SPLIT).expect("write");
    writeln!(txt, "void_factor_66_over_67 = {:.12}", LITHIUM7_VOID_CORRECTION).expect("write");
    writeln!(txt, "closure_factor = {:.12}", closure_factor).expect("write");
    writeln!(txt, "li7h_pred = {:.12e}", li7h_predictive).expect("write");
    writeln!(txt, "li7_tension_ratio = {:.12}", ratio_predictive).expect("write");
    writeln!(txt, "unity_window_0p8_1p4 = {}", (0.8..=1.4).contains(&ratio_predictive))
        .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[channel_coupled]").expect("write");
    writeln!(txt, "direct_fraction_1_over_16 = {:.12}", LITHIUM7_DIRECT_CHANNEL_FRACTION)
        .expect("write");
    writeln!(txt, "be7_fraction_15_over_16 = {:.12}", LITHIUM7_BE7_CHANNEL_FRACTION).expect("write");
    writeln!(txt, "be7_dark_suppression = {:.12}", LITHIUM7_BE7_DARK_SUPPRESSION).expect("write");
    writeln!(txt, "channel_coupled_factor = {:.12}", LITHIUM7_CHANNEL_COUPLED_FACTOR).expect("write");
    writeln!(txt, "li7h_pred = {:.12e}", li7h_channel).expect("write");
    writeln!(txt, "li7_tension_ratio = {:.12}", ratio_channel).expect("write");
    writeln!(txt, "unity_window_0p8_1p4 = {}", (0.8..=1.4).contains(&ratio_channel)).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"eta10\": {:.12},\n  \"li7_observed\": {:.12e},\n  \"observed_anchored_diagnostic\": {{\"amplification\": {:.12}, \"li7h_pred\": {:.12e}, \"tension_ratio\": {:.12}}},\n  \"channel_derived_predictive\": {{\"visible_fraction\": {:.12}, \"reaction_network_gain\": {:.12}, \"source\": {:.12e}, \"direct_component\": {:.12e}, \"be7_component_raw\": {:.12e}, \"be7_component_dark\": {:.12e}, \"occupancy_factor\": {:.12}, \"void_factor\": {:.12}, \"closure_factor\": {:.12}, \"li7h_pred\": {:.12e}, \"tension_ratio\": {:.12}, \"unity_window_0p8_1p4\": {}}},\n  \"channel_coupled\": {{\"direct_fraction\": {:.12}, \"be7_fraction\": {:.12}, \"be7_dark_suppression\": {:.12}, \"channel_factor\": {:.12}, \"li7h_pred\": {:.12e}, \"tension_ratio\": {:.12}, \"unity_window_0p8_1p4\": {}}}\n}}",
        score.eta10,
        LI7H_OBSERVED,
        LITHIUM7_TENSION_AMPLIFICATION,
        li7h_observed_anchor,
        ratio_observed_anchor,
        LITHIUM7_VISIBLE_FRACTION,
        LITHIUM7_REACTION_NETWORK_GAIN,
        li7_source,
        li7_direct,
        li7_be7_raw,
        li7_be7_dark,
        DARK_FRACTION_TOTAL_STATE_SPLIT,
        LITHIUM7_VOID_CORRECTION,
        closure_factor,
        li7h_predictive,
        ratio_predictive,
        (0.8..=1.4).contains(&ratio_predictive),
        LITHIUM7_DIRECT_CHANNEL_FRACTION,
        LITHIUM7_BE7_CHANNEL_FRACTION,
        LITHIUM7_BE7_DARK_SUPPRESSION,
        LITHIUM7_CHANNEL_COUPLED_FACTOR,
        li7h_channel,
        ratio_channel,
        (0.8..=1.4).contains(&ratio_channel)
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "Li-7 ratio observed-anchor={:.6}, predictive={:.6}, channel_coupled={:.6} (eta10={:.4})",
        ratio_observed_anchor, ratio_predictive, ratio_channel, score.eta10
    );
}
