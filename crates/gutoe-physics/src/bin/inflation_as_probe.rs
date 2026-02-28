use gutoe_physics::constants::{
    ALPHA_LEADING_ORDER, DARK_TO_VISIBLE_GEOMETRIC_RATIO, LAMBDA_QG,
};
use gutoe_physics::inflation::{
    evaluate_inflation_gate, inflation_hubble_ratio_structural, scalar_amplitude, slow_roll_epsilon,
    AS_OBSERVED, InflationWindows,
};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let out_dir =
        std::env::var("GUTOE_AS_OUT").unwrap_or_else(|_| "/tmp/bh_renders/inflation_as_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let s = evaluate_inflation_gate(InflationWindows::default());
    let n = s.n_efolds;
    let eps = slow_roll_epsilon(n);
    let h = inflation_hubble_ratio_structural();
    let a_s = scalar_amplitude(n, h);

    let geometric_budget = DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let survival = 1.0 - LAMBDA_QG;
    let signature_split = 3.0 / 6.0;
    let micro_dilution = 1.0 / 486.0_f64.sqrt();
    let alpha_sq = ALPHA_LEADING_ORDER.powi(2);

    let ratio = a_s / AS_OBSERVED;
    let needed_h_factor = (AS_OBSERVED / a_s).sqrt();
    let needed_h_percent = (needed_h_factor - 1.0) * 100.0;

    // If only one multiplicative factor in h changes:
    // h -> f*h  =>  A_s -> f^2*A_s
    let needed_single_factor = needed_h_factor;
    // If correction is shared equally across k independent multiplicative factors in h:
    let k2 = needed_h_factor.sqrt();
    let k3 = needed_h_factor.powf(1.0 / 3.0);
    let k4 = needed_h_factor.powf(0.25);
    let k5 = needed_h_factor.powf(0.2);

    let report = out.join("inflation_as_probe_report.json");
    let mut f = File::create(&report).expect("create report");
    writeln!(f, "{{").expect("write");
    writeln!(
        f,
        "  \"base\": {{\"N\": {:.12}, \"epsilon\": {:.12e}, \"h_over_mpl\": {:.12e}, \"A_s\": {:.12e}, \"A_s_observed_ref\": {:.12e}, \"A_s_ratio\": {:.12}, \"A_s_excess_percent\": {:.9}}},",
        n,
        eps,
        h,
        a_s,
        AS_OBSERVED,
        ratio,
        (ratio - 1.0) * 100.0
    )
    .expect("write");
    writeln!(
        f,
        "  \"h_chain\": {{\"alpha_sq\": {:.12e}, \"geometric_budget_60_over_11\": {:.12}, \"survival_11_over_12\": {:.12}, \"signature_split_half\": {:.12}, \"micro_dilution_1_over_sqrt486\": {:.12}}},",
        alpha_sq, geometric_budget, survival, signature_split, micro_dilution
    )
    .expect("write");
    writeln!(
        f,
        "  \"needed_h_correction\": {{\"multiplicative_factor\": {:.12}, \"percent\": {:.9}, \"if_single_factor_changes\": {:.12}, \"if_shared_2_factors_each\": {:.12}, \"if_shared_3_factors_each\": {:.12}, \"if_shared_4_factors_each\": {:.12}, \"if_shared_5_factors_each\": {:.12}}}",
        needed_h_factor,
        needed_h_percent,
        needed_single_factor,
        k2,
        k3,
        k4,
        k5
    )
    .expect("write");
    writeln!(f, "}}").expect("write");

    println!("wrote {}", report.display());
    println!(
        "A_s={:.6e} (ref {:.6e}) | ratio={:.6} | need h-factor {:.6} ({:.3}%)",
        a_s, AS_OBSERVED, ratio, needed_h_factor, needed_h_percent
    );
}
