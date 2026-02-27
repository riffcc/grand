//! Structural inflation report for GRAND-347.

use gutoe_physics::{evaluate_inflation_gate, InflationWindows};
use std::fs::{self, File};
use std::io::Write;

fn main() {
    let windows = InflationWindows::default();
    let score = evaluate_inflation_gate(windows);

    let out_dir =
        std::env::var("GUTOE_INFLATION_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/inflation_report.txt");
    let json_path = format!("{out_dir}/inflation_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[inflation_structural]").expect("write");
    writeln!(txt, "n_efolds = {:.12}", score.n_efolds).expect("write");
    writeln!(txt, "epsilon = {:.12e}", score.epsilon).expect("write");
    writeln!(txt, "eta = {:.12e}", score.eta).expect("write");
    writeln!(txt, "n_s = {:.12}", score.n_s).expect("write");
    writeln!(txt, "r = {:.12}", score.r).expect("write");
    writeln!(txt, "a_s = {:.12e}", score.a_s).expect("write");
    writeln!(txt, "h_over_mpl = {:.12e}", score.h_over_mpl).expect("write");
    writeln!(txt, "n_reheat = {:.12}", score.n_reheat).expect("write");
    writeln!(txt, "w_reheat = {:.12}", score.w_reheat).expect("write");
    writeln!(txt, "t_reheat_gev = {:.12e}", score.t_reheat_gev).expect("write");
    writeln!(txt, "cmb_proxy_chi2 = {:.12}", score.cmb_proxy_chi2).expect("write");
    writeln!(txt, "cmb_proxy_loglike = {:.12}", score.cmb_proxy_loglike).expect("write");
    writeln!(txt, "n_end = {:.12}", score.n_end).expect("write");
    writeln!(txt, "expansion_factor = {:.12e}", score.expansion_factor).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[windows]").expect("write");
    writeln!(
        txt,
        "n_efolds_range = [{:.3}, {:.3}]",
        windows.n_efolds_min, windows.n_efolds_max
    )
    .expect("write");
    writeln!(txt, "n_s_center = {:.6}", windows.ns_center).expect("write");
    writeln!(txt, "n_s_tol = {:.6}", windows.ns_tol).expect("write");
    writeln!(txt, "r_max = {:.6}", windows.r_max).expect("write");
    writeln!(txt, "a_s_center = {:.12e}", windows.as_center).expect("write");
    writeln!(txt, "a_s_tol = {:.12e}", windows.as_tol).expect("write");
    writeln!(txt, "t_reheat_min_gev = {:.12e}", windows.t_reheat_min_gev).expect("write");
    writeln!(txt, "chi2_max = {:.6}", windows.chi2_max).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[gate]").expect("write");
    writeln!(txt, "n_efolds_ok = {}", score.n_efolds_ok).expect("write");
    writeln!(txt, "n_s_ok = {}", score.n_s_ok).expect("write");
    writeln!(txt, "r_ok = {}", score.r_ok).expect("write");
    writeln!(txt, "a_s_ok = {}", score.a_s_ok).expect("write");
    writeln!(txt, "reheating_ok = {}", score.reheating_ok).expect("write");
    writeln!(txt, "cmb_like_ok = {}", score.cmb_like_ok).expect("write");
    writeln!(txt, "graceful_exit_ok = {}", score.graceful_exit_ok).expect("write");
    writeln!(txt, "passes_all = {}", score.passes_all()).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"n_efolds\": {:.12},\n  \"epsilon\": {:.12e},\n  \"eta\": {:.12e},\n  \"n_s\": {:.12},\n  \"r\": {:.12},\n  \"a_s\": {:.12e},\n  \"h_over_mpl\": {:.12e},\n  \"n_reheat\": {:.12},\n  \"w_reheat\": {:.12},\n  \"t_reheat_gev\": {:.12e},\n  \"cmb_proxy_chi2\": {:.12},\n  \"cmb_proxy_loglike\": {:.12},\n  \"n_end\": {:.12},\n  \"expansion_factor\": {:.12e},\n  \"windows\": {{\"n_efolds_min\": {:.6}, \"n_efolds_max\": {:.6}, \"n_s_center\": {:.6}, \"n_s_tol\": {:.6}, \"r_max\": {:.6}, \"a_s_center\": {:.12e}, \"a_s_tol\": {:.12e}, \"t_reheat_min_gev\": {:.12e}, \"chi2_max\": {:.6}}},\n  \"n_efolds_ok\": {},\n  \"n_s_ok\": {},\n  \"r_ok\": {},\n  \"a_s_ok\": {},\n  \"reheating_ok\": {},\n  \"cmb_like_ok\": {},\n  \"graceful_exit_ok\": {},\n  \"passes_all\": {}\n}}",
        score.n_efolds,
        score.epsilon,
        score.eta,
        score.n_s,
        score.r,
        score.a_s,
        score.h_over_mpl,
        score.n_reheat,
        score.w_reheat,
        score.t_reheat_gev,
        score.cmb_proxy_chi2,
        score.cmb_proxy_loglike,
        score.n_end,
        score.expansion_factor,
        windows.n_efolds_min,
        windows.n_efolds_max,
        windows.ns_center,
        windows.ns_tol,
        windows.r_max,
        windows.as_center,
        windows.as_tol,
        windows.t_reheat_min_gev,
        windows.chi2_max,
        score.n_efolds_ok,
        score.n_s_ok,
        score.r_ok,
        score.a_s_ok,
        score.reheating_ok,
        score.cmb_like_ok,
        score.graceful_exit_ok,
        score.passes_all()
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "inflation: N={:.3}, n_s={:.6}, r={:.6}, A_s={:.3e}, chi2={:.3}, passes_all={}",
        score.n_efolds,
        score.n_s,
        score.r,
        score.a_s,
        score.cmb_proxy_chi2,
        score.passes_all()
    );
}
