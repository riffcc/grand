use anyhow::{anyhow, Context};
use gutoe_physics::{scan_nuclear_chart, ScanConfig, ShellParams};
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
    a: u16,
    binding_mev: f64,
    binding_unc_mev: f64,
}

#[derive(Clone, Copy)]
struct ResidualRow {
    z: u16,
    n: u16,
    a: u16,
    pred_binding_mev: f64,
    obs_binding_mev: f64,
    obs_unc_mev: f64,
    residual_mev: f64,
    abs_residual_mev: f64,
}

fn env_f64(name: &str, default: f64) -> f64 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_u16(name: &str, default: u16) -> u16 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<u16>().ok())
        .unwrap_or(default)
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
    // Header lines in AME2020 are not fixed-width numeric rows.
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

    // Binding energy / A is in keV (cols 55..67, 1-indexed), with uncertainty in keV (cols 69..78).
    let binding_per_a_kev = parse_fixed_float(line, 54, 67)?;
    let binding_per_a_unc_kev = parse_fixed_float(line, 68, 78).unwrap_or(0.0);
    Some(AmeBindingRow {
        z: z_u16,
        n: n_u16,
        a: a_u16,
        binding_mev: binding_per_a_kev * a_u16 as f64 / 1000.0,
        binding_unc_mev: binding_per_a_unc_kev * a_u16 as f64 / 1000.0,
    })
}

fn load_ame2020_bindings(path: &Path) -> anyhow::Result<Vec<AmeBindingRow>> {
    let text =
        fs::read_to_string(path).with_context(|| format!("failed to read {}", path.display()))?;
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

fn main() -> anyhow::Result<()> {
    let out_dir =
        env::var("GUTOE_MASS_PERIODIC_OUT").unwrap_or_else(|_| "/tmp/nuclear_chart".to_string());
    fs::create_dir_all(&out_dir)?;
    let out = PathBuf::from(out_dir);

    let ame_path = PathBuf::from(
        env::var("GUTOE_AME2020_PATH")
            .unwrap_or_else(|_| "/tmp/ame2020_mass_1.mas20.txt".to_string()),
    );
    let ame_url = env::var("GUTOE_AME2020_URL").unwrap_or_else(|_| DEFAULT_AME2020_URL.to_string());
    let auto_fetch = env_bool("GUTOE_AME2020_AUTO_FETCH", true);
    ensure_ame_file(&ame_path, &ame_url, auto_fetch)?;
    let ame_rows = load_ame2020_bindings(&ame_path)?;

    let default_shell = ShellParams::default();
    let cfg = ScanConfig {
        z_min: env_u16("GUTOE_NUCLEAR_Z_MIN", 1),
        z_max: env_u16("GUTOE_NUCLEAR_Z_MAX", 140),
        n_min: env_u16("GUTOE_NUCLEAR_N_MIN", 1),
        n_max: env_u16("GUTOE_NUCLEAR_N_MAX", 260),
        shell: ShellParams {
            amplitude_z: env_f64("GUTOE_NUCLEAR_AMP_Z", default_shell.amplitude_z),
            amplitude_n: env_f64("GUTOE_NUCLEAR_AMP_N", default_shell.amplitude_n),
            shell_amp: env_f64("GUTOE_NUCLEAR_SHELL_AMP", default_shell.shell_amp),
            shell_scale_exp: env_f64(
                "GUTOE_NUCLEAR_SHELL_SCALE_EXP",
                default_shell.shell_scale_exp,
            ),
            use_strutinsky: env_bool("GUTOE_NUCLEAR_USE_STRUTINSKY", default_shell.use_strutinsky),
            strutinsky_gamma: env_f64(
                "GUTOE_NUCLEAR_STRUTINSKY_GAMMA",
                default_shell.strutinsky_gamma,
            ),
            strutinsky_spacing_mev: env_f64(
                "GUTOE_NUCLEAR_STRUTINSKY_SPACING_MEV",
                default_shell.strutinsky_spacing_mev,
            ),
            strutinsky_spin_orbit_mev: env_f64(
                "GUTOE_NUCLEAR_STRUTINSKY_SPIN_ORBIT_MEV",
                default_shell.strutinsky_spin_orbit_mev,
            ),
            strutinsky_coulomb_shift_mev: env_f64(
                "GUTOE_NUCLEAR_STRUTINSKY_COULOMB_SHIFT_MEV",
                default_shell.strutinsky_coulomb_shift_mev,
            ),
            strutinsky_ws_depth_mev: env_f64(
                "GUTOE_NUCLEAR_STRUTINSKY_WS_DEPTH_MEV",
                default_shell.strutinsky_ws_depth_mev,
            ),
            strutinsky_ws_r0_fm: env_f64(
                "GUTOE_NUCLEAR_STRUTINSKY_WS_R0_FM",
                default_shell.strutinsky_ws_r0_fm,
            ),
            strutinsky_ws_diffuseness_fm: env_f64(
                "GUTOE_NUCLEAR_STRUTINSKY_WS_DIFFUSENESS_FM",
                default_shell.strutinsky_ws_diffuseness_fm,
            ),
            strutinsky_ws_a_ref: env_f64(
                "GUTOE_NUCLEAR_STRUTINSKY_WS_A_REF",
                default_shell.strutinsky_ws_a_ref,
            ),
            strutinsky_ws_ref_nosc: env_f64(
                "GUTOE_NUCLEAR_STRUTINSKY_WS_REF_NOSC",
                default_shell.strutinsky_ws_ref_nosc,
            ),
            strutinsky_ws_coulomb_z_ref: env_f64(
                "GUTOE_NUCLEAR_STRUTINSKY_WS_COULOMB_Z_REF",
                default_shell.strutinsky_ws_coulomb_z_ref,
            ),
            strutinsky_mix: env_f64("GUTOE_NUCLEAR_STRUTINSKY_MIX", default_shell.strutinsky_mix),
            sigma_z: env_f64("GUTOE_NUCLEAR_SIGMA_Z", default_shell.sigma_z),
            sigma_n: env_f64("GUTOE_NUCLEAR_SIGMA_N", default_shell.sigma_n),
            proton_magic_weight_coeff: env_f64(
                "GUTOE_NUCLEAR_PROTON_MAGIC_WEIGHT_COEFF",
                default_shell.proton_magic_weight_coeff,
            ),
            proton_magic_weight_cap: env_f64(
                "GUTOE_NUCLEAR_PROTON_MAGIC_WEIGHT_CAP",
                default_shell.proton_magic_weight_cap,
            ),
            neutron_magic_weight_coeff: env_f64(
                "GUTOE_NUCLEAR_NEUTRON_MAGIC_WEIGHT_COEFF",
                default_shell.neutron_magic_weight_coeff,
            ),
            neutron_magic_weight_cap: env_f64(
                "GUTOE_NUCLEAR_NEUTRON_MAGIC_WEIGHT_CAP",
                default_shell.neutron_magic_weight_cap,
            ),
            closure_index_attenuation: env_f64(
                "GUTOE_NUCLEAR_CLOSURE_INDEX_ATTENUATION",
                default_shell.closure_index_attenuation,
            ),
            superheavy_proton_amplitude: env_f64(
                "GUTOE_NUCLEAR_SUPERHEAVY_PROTON_AMP",
                default_shell.superheavy_proton_amplitude,
            ),
            superheavy_proton_sigma: env_f64(
                "GUTOE_NUCLEAR_SUPERHEAVY_PROTON_SIGMA",
                default_shell.superheavy_proton_sigma,
            ),
            superheavy_proton_gate_n_sigma: env_f64(
                "GUTOE_NUCLEAR_SUPERHEAVY_PROTON_GATE_N_SIGMA",
                default_shell.superheavy_proton_gate_n_sigma,
            ),
            heavy_target_z: env_f64("GUTOE_NUCLEAR_HEAVY_TARGET_Z", default_shell.heavy_target_z),
            heavy_target_n: env_f64("GUTOE_NUCLEAR_HEAVY_TARGET_N", default_shell.heavy_target_n),
            heavy_sigma_z: env_f64("GUTOE_NUCLEAR_HEAVY_SIGMA_Z", default_shell.heavy_sigma_z),
            heavy_sigma_n: env_f64("GUTOE_NUCLEAR_HEAVY_SIGMA_N", default_shell.heavy_sigma_n),
            heavy_amplitude: env_f64("GUTOE_NUCLEAR_HEAVY_AMP", default_shell.heavy_amplitude),
            heavy_gate_z_min: env_u16(
                "GUTOE_NUCLEAR_HEAVY_GATE_Z_MIN",
                default_shell.heavy_gate_z_min,
            ),
            heavy_gate_n_min: env_u16(
                "GUTOE_NUCLEAR_HEAVY_GATE_N_MIN",
                default_shell.heavy_gate_n_min,
            ),
            z50_isovector_valley_amplitude: env_f64(
                "GUTOE_NUCLEAR_Z50_ISOVECTOR_VALLEY_AMP",
                default_shell.z50_isovector_valley_amplitude,
            ),
            z50_isovector_beta_coeff: env_f64(
                "GUTOE_NUCLEAR_Z50_ISOVECTOR_BETA_COEFF",
                default_shell.z50_isovector_beta_coeff,
            ),
        },
        ..ScanConfig::default()
    };

    let records = scan_nuclear_chart(cfg);
    let pred_by_zn: BTreeMap<(u16, u16), f64> = records
        .iter()
        .map(|r| ((r.z, r.n), r.binding_mev))
        .collect();

    let mut residuals = Vec::new();
    for row in ame_rows.iter().copied() {
        if let Some(pred) = pred_by_zn.get(&(row.z, row.n)).copied() {
            let residual = pred - row.binding_mev;
            residuals.push(ResidualRow {
                z: row.z,
                n: row.n,
                a: row.a,
                pred_binding_mev: pred,
                obs_binding_mev: row.binding_mev,
                obs_unc_mev: row.binding_unc_mev,
                residual_mev: residual,
                abs_residual_mev: residual.abs(),
            });
        }
    }

    if residuals.is_empty() {
        return Err(anyhow!(
            "no AME2020 rows matched scan range (z={}..{}, n={}..{})",
            cfg.z_min,
            cfg.z_max,
            cfg.n_min,
            cfg.n_max
        ));
    }

    let n = residuals.len() as f64;
    let rmse = (residuals
        .iter()
        .map(|r| r.residual_mev * r.residual_mev)
        .sum::<f64>()
        / n)
        .sqrt();
    let mae = residuals.iter().map(|r| r.abs_residual_mev).sum::<f64>() / n;
    let bias = residuals.iter().map(|r| r.residual_mev).sum::<f64>() / n;

    let mut weighted_num = 0.0;
    let mut weighted_den = 0.0;
    let mut chi2_num = 0.0;
    let mut chi2_count = 0usize;
    for r in &residuals {
        if r.obs_unc_mev > 0.0 {
            let w = 1.0 / (r.obs_unc_mev * r.obs_unc_mev);
            weighted_num += w * r.residual_mev * r.residual_mev;
            weighted_den += w;
            chi2_num += (r.residual_mev / r.obs_unc_mev).powi(2);
            chi2_count += 1;
        }
    }
    let wrmse = if weighted_den > 0.0 {
        (weighted_num / weighted_den).sqrt()
    } else {
        f64::NAN
    };
    let reduced_chi2 = if chi2_count > 0 {
        chi2_num / chi2_count as f64
    } else {
        f64::NAN
    };

    let mut all_csv = String::from(
        "Z,N,A,pred_binding_mev,obs_binding_mev,obs_unc_mev,residual_mev,abs_residual_mev\n",
    );
    for r in &residuals {
        all_csv.push_str(&format!(
            "{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
            r.z,
            r.n,
            r.a,
            r.pred_binding_mev,
            r.obs_binding_mev,
            r.obs_unc_mev,
            r.residual_mev,
            r.abs_residual_mev
        ));
    }
    fs::write(out.join("ame2020_residuals.csv"), all_csv)?;

    let mut top = residuals.clone();
    top.sort_by(|a, b| b.abs_residual_mev.total_cmp(&a.abs_residual_mev));
    top.truncate(50);
    let mut top_csv = String::from(
        "Z,N,A,pred_binding_mev,obs_binding_mev,obs_unc_mev,residual_mev,abs_residual_mev\n",
    );
    for r in &top {
        top_csv.push_str(&format!(
            "{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
            r.z,
            r.n,
            r.a,
            r.pred_binding_mev,
            r.obs_binding_mev,
            r.obs_unc_mev,
            r.residual_mev,
            r.abs_residual_mev
        ));
    }
    fs::write(out.join("ame2020_residuals_top50.csv"), top_csv)?;

    let summary = format!(
        "{{\n  \"ame2020_source\": \"{}\",\n  \"ame2020_rows_parsed\": {},\n  \"rows_matched\": {},\n  \"scan_range\": {{\"z_min\": {}, \"z_max\": {}, \"n_min\": {}, \"n_max\": {}}},\n  \"rmse_mev\": {:.6},\n  \"mae_mev\": {:.6},\n  \"bias_mev\": {:.6},\n  \"weighted_rmse_mev\": {:.6},\n  \"reduced_chi2\": {:.6}\n}}\n",
        ame_path.display(),
        ame_rows.len(),
        residuals.len(),
        cfg.z_min,
        cfg.z_max,
        cfg.n_min,
        cfg.n_max,
        rmse,
        mae,
        bias,
        wrmse,
        reduced_chi2
    );
    fs::write(out.join("ame2020_benchmark.json"), summary)?;

    println!("Wrote {}", out.join("ame2020_benchmark.json").display());
    println!("Wrote {}", out.join("ame2020_residuals.csv").display());
    println!(
        "Wrote {}",
        out.join("ame2020_residuals_top50.csv").display()
    );
    println!(
        "AME2020 benchmark: matched={} rmse={:.4} MeV mae={:.4} MeV bias={:.4} MeV wrmse={:.4} MeV chi2ν={:.3}",
        residuals.len(),
        rmse,
        mae,
        bias,
        wrmse,
        reduced_chi2
    );
    Ok(())
}
