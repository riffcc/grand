use gutoe_physics::{
    magic_s2n_discontinuities, magic_s2n_summary, rank_island_candidates, scan_nuclear_chart,
    write_magic_discontinuities_csv, write_magic_summary_csv, write_records_csv, NucleusRecord, ScanConfig,
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
    let cfg = ScanConfig {
        z_min: env_u16("GUTOE_NUCLEAR_Z_MIN", 1),
        z_max: env_u16("GUTOE_NUCLEAR_Z_MAX", 140),
        n_min: env_u16("GUTOE_NUCLEAR_N_MIN", 1),
        n_max: env_u16("GUTOE_NUCLEAR_N_MAX", 260),
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
    let islands = rank_island_candidates(&records, min_island_z, top_k);

    write_records_csv(out.join("nuclides.csv"), &records)?;
    write_magic_discontinuities_csv(out.join("magic_s2n_discontinuities.csv"), &magic)?;
    write_magic_summary_csv(out.join("magic_s2n_summary.csv"), &magic_summary)?;
    write_top_islands_csv(out.join("top_islands.csv"), &islands)?;

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
    println!("  {}", out.join("top_islands.csv").display());
    println!("  {}", out.join("summary.txt").display());
    Ok(())
}
