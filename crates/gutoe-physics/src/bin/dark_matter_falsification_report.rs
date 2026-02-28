//! GRAND-346: dataset-backed dark-matter falsification report.

use gutoe_physics::dark_matter_falsification::{
    evaluate_dark_matter_gate, DarkMatterBranchScorecard, DarkMatterFalsificationWindows,
};
use gutoe_physics::dark_sector::DarkSectorBranch;
use std::fs::{self, File};
use std::io::Write;

fn branch_name(branch: DarkSectorBranch) -> &'static str {
    match branch {
        DarkSectorBranch::Particle => "particle",
        DarkSectorBranch::Geometric => "geometric",
        DarkSectorBranch::Unified => "unified",
    }
}

fn write_branch_txt(out: &mut File, s: DarkMatterBranchScorecard) {
    let m = s.metrics;
    writeln!(out, "[{}]", branch_name(s.branch)).expect("write");
    writeln!(out, "n_points = {}", m.n_points).expect("write");
    writeln!(out, "rotation_rmse_kms = {:.6}", m.rotation_rmse_kms).expect("write");
    writeln!(out, "rotation_mape = {:.9}", m.rotation_mape).expect("write");
    writeln!(out, "rotation_chi2_ndof = {:.6}", m.rotation_chi2_ndof).expect("write");
    writeln!(
        out,
        "lensing_proxy_rmse_rad = {:.12e}",
        m.lensing_proxy_rmse_rad
    )
    .expect("write");
    writeln!(out, "lensing_proxy_mape = {:.9}", m.lensing_proxy_mape).expect("write");
    writeln!(
        out,
        "predicted_dm_fraction = {:.9}",
        m.predicted_dm_fraction
    )
    .expect("write");
    writeln!(out, "observed_dm_fraction = {:.9}", m.observed_dm_fraction).expect("write");
    writeln!(out, "dm_fraction_delta = {:.9}", m.dm_fraction_delta).expect("write");
    writeln!(out, "rotation_ok = {}", s.rotation_ok).expect("write");
    writeln!(out, "lensing_ok = {}", s.lensing_ok).expect("write");
    writeln!(out, "cmb_fraction_ok = {}", s.cmb_fraction_ok).expect("write");
    writeln!(out, "passes_all = {}", s.passes_all()).expect("write");
    writeln!(out).expect("write");
}

fn main() {
    let windows = DarkMatterFalsificationWindows::default();
    let scorecards = evaluate_dark_matter_gate(windows);

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let txt_path = format!("{out_dir}/dark_matter_falsification_report.txt");
    let json_path = format!("{out_dir}/dark_matter_falsification_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "GRAND-346 dark-matter falsification gate").expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[windows]").expect("write");
    writeln!(txt, "rotation_mape_max = {:.6}", windows.rotation_mape_max).expect("write");
    writeln!(
        txt,
        "lensing_proxy_mape_max = {:.6}",
        windows.lensing_proxy_mape_max
    )
    .expect("write");
    writeln!(
        txt,
        "dm_fraction_delta_abs_max = {:.6}",
        windows.dm_fraction_delta_abs_max
    )
    .expect("write");
    writeln!(txt).expect("write");
    for s in &scorecards {
        write_branch_txt(&mut txt, *s);
    }

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"windows\": {{\"rotation_mape_max\": {:.12}, \"lensing_proxy_mape_max\": {:.12}, \"dm_fraction_delta_abs_max\": {:.12}}},",
        windows.rotation_mape_max, windows.lensing_proxy_mape_max, windows.dm_fraction_delta_abs_max
    )
    .expect("write");
    let write_branch = |json: &mut File, s: DarkMatterBranchScorecard, is_last: bool| {
        let m = s.metrics;
        writeln!(
            json,
            "  \"{}\": {{\"n_points\": {}, \"rotation_rmse_kms\": {:.12}, \"rotation_mape\": {:.12}, \"rotation_chi2_ndof\": {:.12}, \"lensing_proxy_rmse_rad\": {:.12e}, \"lensing_proxy_mape\": {:.12}, \"predicted_dm_fraction\": {:.12}, \"observed_dm_fraction\": {:.12}, \"dm_fraction_delta\": {:.12}, \"rotation_ok\": {}, \"lensing_ok\": {}, \"cmb_fraction_ok\": {}, \"passes_all\": {}}}{}",
            branch_name(s.branch),
            m.n_points,
            m.rotation_rmse_kms,
            m.rotation_mape,
            m.rotation_chi2_ndof,
            m.lensing_proxy_rmse_rad,
            m.lensing_proxy_mape,
            m.predicted_dm_fraction,
            m.observed_dm_fraction,
            m.dm_fraction_delta,
            s.rotation_ok,
            s.lensing_ok,
            s.cmb_fraction_ok,
            s.passes_all(),
            if is_last { "" } else { "," },
        )
        .expect("write");
    };

    for (i, s) in scorecards.iter().enumerate() {
        write_branch(&mut json, *s, i + 1 == scorecards.len());
    }
    writeln!(
        json,
        ",\n  \"summary\": {{\"at_least_one_branch_passes_all\": {}}}\n}}",
        scorecards.iter().any(|s| s.passes_all())
    )
    .expect("write");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
}
