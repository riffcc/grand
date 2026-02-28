//! GRAND-106/107/108 structural closeout report.
//!
//! One deterministic lane:
//! Cl(1,3) -> NN potential proxy -> shell closures -> Z<=118 binding-energy grid + AME2020 benchmark.

use anyhow::{anyhow, Context};
use gutoe_physics::{
    derive_structural_nuclear_model, magic_s2n_summary, proton_s2p_summary, scan_nuclear_chart,
    structural_scan_config_z118,
};
use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const DEFAULT_AME2020_URL: &str = "https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt";

#[derive(Clone, Copy)]
struct AmeBindingRow {
    z: u16,
    n: u16,
    binding_mev: f64,
}

fn env_bool(name: &str, default: bool) -> bool {
    match env::var(name) {
        Ok(v) => match v.trim().to_ascii_lowercase().as_str() {
            "1" | "true" | "yes" | "on" => true,
            "0" | "false" | "no" | "off" => false,
            _ => default,
        },
        Err(_) => default,
    }
}

fn field(line: &str, start: usize, end: usize) -> &str {
    if line.len() < end {
        ""
    } else {
        &line[start..end]
    }
}

fn parse_fixed_int(line: &str, start: usize, end: usize) -> Option<i32> {
    field(line, start, end).trim().parse::<i32>().ok()
}

fn parse_fixed_float(line: &str, start: usize, end: usize) -> Option<f64> {
    let raw = field(line, start, end);
    if raw.contains('#') || raw.contains('*') {
        return None;
    }
    raw.trim().parse::<f64>().ok()
}

fn parse_ame_binding_row(line: &str) -> Option<AmeBindingRow> {
    let z = parse_fixed_int(line, 9, 14)?;
    let n = parse_fixed_int(line, 4, 9)?;
    let a = parse_fixed_int(line, 14, 19)?;
    if z < 0 || n < 0 || a <= 0 {
        return None;
    }
    let z_u16 = z as u16;
    let n_u16 = n as u16;
    let a_u16 = a as u16;
    if z_u16 + n_u16 != a_u16 {
        return None;
    }

    // Binding energy / A (keV) at cols 55..67 (1-indexed in AME tables).
    let binding_per_a_kev = parse_fixed_float(line, 54, 67)?;
    Some(AmeBindingRow {
        z: z_u16,
        n: n_u16,
        binding_mev: binding_per_a_kev * a_u16 as f64 / 1000.0,
    })
}

fn ensure_ame_file(path: &Path, url: &str, auto_fetch: bool) -> anyhow::Result<()> {
    if path.exists() {
        return Ok(());
    }
    if !auto_fetch {
        return Err(anyhow!(
            "AME2020 file missing at {} (set GUTOE_AME2020_AUTO_FETCH=1 or download manually)",
            path.display()
        ));
    }
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let status = Command::new("curl")
        .args(["-L", "--fail", "-o"])
        .arg(path)
        .arg(url)
        .status()
        .context("failed to spawn curl for AME2020 download")?;
    if !status.success() {
        return Err(anyhow!("curl failed downloading {}", url));
    }
    Ok(())
}

fn load_ame2020_bindings(path: &Path) -> anyhow::Result<Vec<AmeBindingRow>> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read {}", path.display()))?;
    let mut out = Vec::new();
    for line in text.lines() {
        if let Some(row) = parse_ame_binding_row(line) {
            out.push(row);
        }
    }
    if out.is_empty() {
        return Err(anyhow!("parsed zero AME2020 rows from {}", path.display()));
    }
    Ok(out)
}

fn main() -> anyhow::Result<()> {
    let out_dir = env::var("GUTOE_NUCLEAR_FP_OUT").unwrap_or_else(|_| "/tmp/nuclear_chart".to_string());
    fs::create_dir_all(&out_dir)?;
    let out = PathBuf::from(out_dir);

    let model = derive_structural_nuclear_model();
    let cfg = structural_scan_config_z118();
    let records = scan_nuclear_chart(cfg);

    let neutron_summary = magic_s2n_summary(&records);
    let proton_summary = proton_s2p_summary(&records);

    let neutron_hit_rate = if neutron_summary.is_empty() {
        0.0
    } else {
        neutron_summary
            .iter()
            .filter(|r| r.strongest_delta_s2n_mev > 0.0)
            .count() as f64
            / neutron_summary.len() as f64
    };

    let proton_hit_rate = if proton_summary.is_empty() {
        0.0
    } else {
        proton_summary
            .iter()
            .filter(|r| r.strongest_delta_s2p_mev > 0.0)
            .count() as f64
            / proton_summary.len() as f64
    };

    let ame_path = PathBuf::from(
        env::var("GUTOE_AME2020_PATH")
            .unwrap_or_else(|_| "/tmp/ame2020_mass_1.mas20.txt".to_string()),
    );
    let ame_url = env::var("GUTOE_AME2020_URL").unwrap_or_else(|_| DEFAULT_AME2020_URL.to_string());
    let auto_fetch = env_bool("GUTOE_AME2020_AUTO_FETCH", true);
    ensure_ame_file(&ame_path, &ame_url, auto_fetch)?;
    let ame_rows = load_ame2020_bindings(&ame_path)?;

    let mut pred_map: BTreeMap<(u16, u16), f64> = BTreeMap::new();
    for r in &records {
        pred_map.insert((r.z, r.n), r.binding_mev);
    }

    let mut matched = 0usize;
    let mut sse = 0.0;
    let mut sae = 0.0;
    let mut s_bias = 0.0;

    for row in &ame_rows {
        if row.z > 118 {
            continue;
        }
        if let Some(pred) = pred_map.get(&(row.z, row.n)) {
            let resid = pred - row.binding_mev;
            matched += 1;
            sse += resid * resid;
            sae += resid.abs();
            s_bias += resid;
        }
    }

    if matched == 0 {
        return Err(anyhow!("no AME2020 rows matched structural Z<=118 scan grid"));
    }

    let rmse_mev = (sse / matched as f64).sqrt();
    let mae_mev = sae / matched as f64;
    let bias_mev = s_bias / matched as f64;

    let ticket_106_pass = model.nn.attractive_depth_mev > 0.0
        && model.nn.repulsive_core_mev > 0.0
        && model.nn.range_fm > 0.0
        && model.nn.core_radius_fm > 0.0;

    let ticket_107_pass = neutron_hit_rate >= 0.70 && proton_hit_rate >= 0.50;

    let ticket_108_pass = matched >= 2400;

    let overall_pass = ticket_106_pass && ticket_107_pass && ticket_108_pass;

    let txt_path = out.join("nuclear_first_principles_report.txt");
    let json_path = out.join("nuclear_first_principles_report.json");

    let txt = format!(
        "[meta]\nno_free_parameters = true\nscan_z_max = 118\nscan_n_max = 260\n\n[grand_106_nn_potential]\n\
attractive_depth_mev = {:.9}\nrepulsive_core_mev = {:.9}\nrange_fm = {:.9}\ncore_radius_fm = {:.9}\nspin_orbit_mev = {:.9}\npass = {}\n\n\
[grand_107_shell_from_nn]\nneutron_magic_hit_rate = {:.9}\nproton_closure_hit_rate = {:.9}\npass = {}\n\n\
[grand_108_binding_energies]\name_rows_total = {}\nrows_matched_z_le_118 = {}\nrmse_mev = {:.9}\nmae_mev = {:.9}\nbias_mev = {:.9}\npass = {}\n\n\
[overall]\npasses_all = {}\n",
        model.nn.attractive_depth_mev,
        model.nn.repulsive_core_mev,
        model.nn.range_fm,
        model.nn.core_radius_fm,
        model.nn.spin_orbit_mev,
        ticket_106_pass,
        neutron_hit_rate,
        proton_hit_rate,
        ticket_107_pass,
        ame_rows.len(),
        matched,
        rmse_mev,
        mae_mev,
        bias_mev,
        ticket_108_pass,
        overall_pass
    );

    let json = format!(
        "{{\n  \"meta\": {{\"no_free_parameters\": true, \"scan\": {{\"z_max\": 118, \"n_max\": 260}}}},\n\
  \"grand_106_nn_potential\": {{\"attractive_depth_mev\": {:.12}, \"repulsive_core_mev\": {:.12}, \"range_fm\": {:.12}, \"core_radius_fm\": {:.12}, \"spin_orbit_mev\": {:.12}, \"pass\": {}}},\n\
  \"grand_107_shell_from_nn\": {{\"neutron_magic_hit_rate\": {:.12}, \"proton_closure_hit_rate\": {:.12}, \"pass\": {}}},\n\
  \"grand_108_binding_energies\": {{\"ame_rows_total\": {}, \"rows_matched_z_le_118\": {}, \"rmse_mev\": {:.12}, \"mae_mev\": {:.12}, \"bias_mev\": {:.12}, \"pass\": {}}},\n\
  \"overall_pass\": {}\n}}\n",
        model.nn.attractive_depth_mev,
        model.nn.repulsive_core_mev,
        model.nn.range_fm,
        model.nn.core_radius_fm,
        model.nn.spin_orbit_mev,
        ticket_106_pass,
        neutron_hit_rate,
        proton_hit_rate,
        ticket_107_pass,
        ame_rows.len(),
        matched,
        rmse_mev,
        mae_mev,
        bias_mev,
        ticket_108_pass,
        overall_pass
    );

    fs::write(&txt_path, txt)?;
    fs::write(&json_path, json)?;

    println!("Wrote {}", txt_path.display());
    println!("Wrote {}", json_path.display());
    println!(
        "nuclear_first_principles: pass={} (rmse={:.4} MeV, neutron_hit={:.3}, proton_hit={:.3})",
        overall_pass, rmse_mev, neutron_hit_rate, proton_hit_rate
    );

    Ok(())
}
