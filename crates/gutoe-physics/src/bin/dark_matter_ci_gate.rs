//! Dark-matter falsification CI gate.
//!
//! Hard-gates the unified dark-sector branch against the SPARC + CMB windows
//! used in GRAND-346. If the unified branch fails any window, exit non-zero.

use gutoe_physics::{
    evaluate_dark_matter_gate, DarkMatterBranchScorecard, DarkMatterFalsificationWindows,
    DarkSectorBranch,
};
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn branch_label(branch: DarkSectorBranch) -> &'static str {
    match branch {
        DarkSectorBranch::Particle => "particle",
        DarkSectorBranch::Geometric => "geometric",
        DarkSectorBranch::Unified => "unified",
    }
}

fn parse_branch(s: &str) -> Option<DarkSectorBranch> {
    match s.trim().to_ascii_lowercase().as_str() {
        "particle" => Some(DarkSectorBranch::Particle),
        "geometric" => Some(DarkSectorBranch::Geometric),
        "unified" => Some(DarkSectorBranch::Unified),
        _ => None,
    }
}

fn print_scorecard(s: &DarkMatterBranchScorecard) {
    println!(
        "{}: pass={} (rot_mape={:.6}, lens_mape={:.6}, dm_delta={:+.6})",
        branch_label(s.branch),
        s.passes_all(),
        s.metrics.rotation_mape,
        s.metrics.lensing_proxy_mape,
        s.metrics.dm_fraction_delta
    );
}

fn write_json(
    file: &mut File,
    target: DarkSectorBranch,
    windows: DarkMatterFalsificationWindows,
    scorecards: &[DarkMatterBranchScorecard],
    overall_pass: bool,
) {
    writeln!(file, "{{").expect("write json");
    writeln!(file, "  \"target_branch\": \"{}\",", branch_label(target)).expect("write json");
    writeln!(
        file,
        "  \"windows\": {{\"rotation_mape_max\": {:.12}, \"lensing_proxy_mape_max\": {:.12}, \"dm_fraction_delta_abs_max\": {:.12}}},",
        windows.rotation_mape_max, windows.lensing_proxy_mape_max, windows.dm_fraction_delta_abs_max
    )
    .expect("write json");
    writeln!(file, "  \"overall_pass\": {},", overall_pass).expect("write json");
    writeln!(file, "  \"branches\": {{").expect("write json");

    for (i, s) in scorecards.iter().enumerate() {
        let comma = if i + 1 == scorecards.len() { "" } else { "," };
        writeln!(
            file,
            "    \"{}\": {{\"pass\": {}, \"rotation_ok\": {}, \"lensing_ok\": {}, \"cmb_fraction_ok\": {}, \"rotation_mape\": {:.12}, \"lensing_proxy_mape\": {:.12}, \"dm_fraction_delta\": {:.12}, \"predicted_dm_fraction\": {:.12}, \"observed_dm_fraction\": {:.12}}}{}",
            branch_label(s.branch),
            s.passes_all(),
            s.rotation_ok,
            s.lensing_ok,
            s.cmb_fraction_ok,
            s.metrics.rotation_mape,
            s.metrics.lensing_proxy_mape,
            s.metrics.dm_fraction_delta,
            s.metrics.predicted_dm_fraction,
            s.metrics.observed_dm_fraction,
            comma
        )
        .expect("write json");
    }

    writeln!(file, "  }}").expect("write json");
    writeln!(file, "}}").expect("write json");
}

fn main() {
    let windows = DarkMatterFalsificationWindows::default();
    let target_branch = std::env::var("GUTOE_DARK_GATE_BRANCH")
        .ok()
        .and_then(|s| parse_branch(&s))
        .unwrap_or(DarkSectorBranch::Unified);

    let scorecards = evaluate_dark_matter_gate(windows);
    for s in &scorecards {
        print_scorecard(s);
    }

    let Some(target) = scorecards.iter().find(|s| s.branch == target_branch) else {
        eprintln!(
            "target branch '{}' missing in scorecards",
            branch_label(target_branch)
        );
        process::exit(2);
    };

    let out_dir =
        std::env::var("GUTOE_DARK_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/dark_matter_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    let overall_pass = target.passes_all();
    write_json(&mut json, target_branch, windows, &scorecards, overall_pass);
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "dark-matter gate FAIL for '{}': rotation_ok={}, lensing_ok={}, cmb_fraction_ok={}",
            branch_label(target_branch),
            target.rotation_ok,
            target.lensing_ok,
            target.cmb_fraction_ok
        );
        process::exit(2);
    }
}
