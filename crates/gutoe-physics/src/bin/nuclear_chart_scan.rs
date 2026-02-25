use gutoe_physics::{
    magic_s2n_discontinuities, magic_s2n_summary, proton_s2p_discontinuities, proton_s2p_summary,
    rank_island_candidates, scan_nuclear_chart, write_magic_discontinuities_csv, write_magic_summary_csv,
    write_proton_discontinuities_csv, write_proton_summary_csv, write_records_csv, NucleusRecord, ScanConfig,
};
use std::env;
use std::fs;
use std::path::PathBuf;

fn env_u16(name: &str, default: u16) -> u16 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<u16>().ok())
        .unwrap_or(default)
}

fn env_usize(name: &str, default: usize) -> usize {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
}

fn env_f64(name: &str, default: f64) -> f64 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn write_top_islands_csv(path: PathBuf, rows: &[NucleusRecord]) -> std::io::Result<()> {
    let mut csv = String::from(
        "rank,Z,N,A,binding_per_nucleon_mev,s2n_mev,s2p_mev,fissility,fission_barrier_mev,sf_log10_half_life_s,stability_score\n",
    );
    for (idx, r) in rows.iter().enumerate() {
        csv.push_str(&format!(
            "{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
            idx + 1,
            r.z,
            r.n,
            r.a,
            r.binding_per_nucleon_mev,
            r.s2n_mev.unwrap_or(f64::NAN),
            r.s2p_mev.unwrap_or(f64::NAN),
            r.fissility,
            r.fission_barrier_mev,
            r.sf_log10_half_life_s,
            r.stability_score
        ));
    }
    fs::write(path, csv)
}

fn main() -> anyhow::Result<()> {
    let output_dir = env::var("GUTOE_NUCLEAR_OUT").unwrap_or_else(|_| "/tmp/nuclear_chart".to_string());
    let default_shell = gutoe_physics::ShellParams::default();
    let cfg = ScanConfig {
        z_min: env_u16("GUTOE_NUCLEAR_Z_MIN", 1),
        z_max: env_u16("GUTOE_NUCLEAR_Z_MAX", 140),
        n_min: env_u16("GUTOE_NUCLEAR_N_MIN", 1),
        n_max: env_u16("GUTOE_NUCLEAR_N_MAX", 260),
        shell: gutoe_physics::ShellParams {
            amplitude_z: env_f64("GUTOE_NUCLEAR_AMP_Z", default_shell.amplitude_z),
            amplitude_n: env_f64("GUTOE_NUCLEAR_AMP_N", default_shell.amplitude_n),
            shell_amp: env_f64("GUTOE_NUCLEAR_SHELL_AMP", default_shell.shell_amp),
            shell_scale_exp: env_f64("GUTOE_NUCLEAR_SHELL_SCALE_EXP", default_shell.shell_scale_exp),
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
            heavy_gate_z_min: env_u16("GUTOE_NUCLEAR_HEAVY_GATE_Z_MIN", default_shell.heavy_gate_z_min),
            heavy_gate_n_min: env_u16("GUTOE_NUCLEAR_HEAVY_GATE_N_MIN", default_shell.heavy_gate_n_min),
        },
        ..ScanConfig::default()
    };
    let top_k = env_usize("GUTOE_NUCLEAR_TOP_K", 40);
    let min_island_z = env_u16("GUTOE_NUCLEAR_MIN_ISLAND_Z", 104);

    fs::create_dir_all(&output_dir)?;
    let out = PathBuf::from(output_dir);

    println!(
        "GUTOE nuclear chart scan: Z={}..{}, N={}..{}",
        cfg.z_min, cfg.z_max, cfg.n_min, cfg.n_max
    );

    let records = scan_nuclear_chart(cfg);
    let magic = magic_s2n_discontinuities(&records, top_k);
    let magic_summary = magic_s2n_summary(&records);
    let proton = proton_s2p_discontinuities(&records, top_k);
    let proton_summary = proton_s2p_summary(&records);
    let islands = rank_island_candidates(&records, min_island_z, top_k);

    write_records_csv(out.join("nuclides.csv"), &records)?;
    write_magic_discontinuities_csv(out.join("magic_s2n_discontinuities.csv"), &magic)?;
    write_magic_summary_csv(out.join("magic_s2n_summary.csv"), &magic_summary)?;
    write_proton_discontinuities_csv(out.join("proton_s2p_discontinuities.csv"), &proton)?;
    write_proton_summary_csv(out.join("proton_s2p_summary.csv"), &proton_summary)?;
    write_top_islands_csv(out.join("top_islands.csv"), &islands)?;

    let slice_z_min = env_u16("GUTOE_NUCLEAR_SLICE_Z_MIN", 104);
    let slice_z_max = env_u16("GUTOE_NUCLEAR_SLICE_Z_MAX", 126);
    let mut slice_csv = String::from(
        "Z,N,A,binding_per_nucleon_mev,s2n_mev,s2p_mev,shell_bonus_mev,shell_bonus_baseline_mev,shell_bonus_heavy_mev,shell_bonus_superheavy_proton_mev,fissility,fission_barrier_mev,sf_log10_half_life_s,stability_score\n",
    );
    for r in records.iter().filter(|r| r.z >= slice_z_min && r.z <= slice_z_max) {
        slice_csv.push_str(&format!(
            "{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
            r.z,
            r.n,
            r.a,
            r.binding_per_nucleon_mev,
            r.s2n_mev.unwrap_or(f64::NAN),
            r.s2p_mev.unwrap_or(f64::NAN),
            r.shell_bonus_mev,
            r.shell_bonus_baseline_mev,
            r.shell_bonus_heavy_mev,
            r.shell_bonus_superheavy_proton_mev,
            r.fissility,
            r.fission_barrier_mev,
            r.sf_log10_half_life_s,
            r.stability_score
        ));
    }
    fs::write(out.join("zslice_superheavy.csv"), slice_csv)?;

    let valley: Vec<&NucleusRecord> = records.iter().filter(|r| r.beta_optimal_for_a).collect();
    let mut valley_csv = String::from("A,Z,N,binding_per_nucleon_mev\n");
    for r in &valley {
        valley_csv.push_str(&format!(
            "{},{},{},{:.6}\n",
            r.a, r.z, r.n, r.binding_per_nucleon_mev
        ));
    }
    fs::write(out.join("valley_of_stability.csv"), valley_csv)?;

    let summary = format!(
        "rows={}\nvalley_rows={}\nmagic_rows={}\nisland_rows={}\noutput_dir={}\n",
        records.len(),
        valley.len(),
        magic.len(),
        islands.len(),
        out.display()
    );
    fs::write(out.join("summary.txt"), summary)?;

    println!("Wrote:");
    println!("  {}", out.join("nuclides.csv").display());
    println!("  {}", out.join("valley_of_stability.csv").display());
    println!("  {}", out.join("magic_s2n_discontinuities.csv").display());
    println!("  {}", out.join("magic_s2n_summary.csv").display());
    println!("  {}", out.join("proton_s2p_discontinuities.csv").display());
    println!("  {}", out.join("proton_s2p_summary.csv").display());
    println!("  {}", out.join("top_islands.csv").display());
    println!("  {}", out.join("zslice_superheavy.csv").display());
    println!("  {}", out.join("summary.txt").display());
    Ok(())
}
