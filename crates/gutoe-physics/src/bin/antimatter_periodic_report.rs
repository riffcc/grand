use gutoe_physics::{
    rank_island_candidates_with_config, scan_nuclear_chart, IslandRankingConfig, NucleusRecord,
    ScanConfig,
};
use std::collections::BTreeMap;
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

fn is_local_peak(best_by_z: &BTreeMap<u16, NucleusRecord>, z: u16) -> bool {
    let Some(cur) = best_by_z.get(&z) else {
        return false;
    };
    let Some(prev) = z.checked_sub(1).and_then(|zz| best_by_z.get(&zz)) else {
        return false;
    };
    let Some(next) = z.checked_add(1).and_then(|zz| best_by_z.get(&zz)) else {
        return false;
    };
    cur.stability_score >= prev.stability_score && cur.stability_score >= next.stability_score
}

fn first_sustained_supercritical(
    best_by_z: &BTreeMap<u16, NucleusRecord>,
    z_min: u16,
    window: usize,
) -> Option<u16> {
    let mut zs: Vec<u16> = best_by_z.keys().copied().filter(|&z| z >= z_min).collect();
    zs.sort_unstable();
    for i in 0..zs.len() {
        if i + window > zs.len() {
            break;
        }
        let ok = zs[i..i + window]
            .iter()
            .all(|z| best_by_z.get(z).map(|r| r.fissility > 1.0).unwrap_or(false));
        if ok {
            return Some(zs[i]);
        }
    }
    None
}

fn main() -> anyhow::Result<()> {
    let out_dir = env::var("GUTOE_ANTIMATTER_OUT")
        .unwrap_or_else(|_| "/tmp/nuclear_chart_antimatter".to_string());
    fs::create_dir_all(&out_dir)?;
    let out = PathBuf::from(out_dir);

    let cfg = ScanConfig {
        z_min: env_u16("GUTOE_NUCLEAR_Z_MIN", 1),
        z_max: env_u16("GUTOE_NUCLEAR_Z_MAX", 300),
        n_min: env_u16("GUTOE_NUCLEAR_N_MIN", 1),
        n_max: env_u16("GUTOE_NUCLEAR_N_MAX", 520),
        ..ScanConfig::default()
    };
    let top_k = env_usize("GUTOE_NUCLEAR_TOP_K", 40);
    let _min_island_z = env_u16("GUTOE_NUCLEAR_MIN_ISLAND_Z", 104);

    // Matter surface.
    let matter = scan_nuclear_chart(cfg);

    // Antimatter surface under strict CPT map:
    // charge signs flip, but nuclear energies depend on even powers of charge.
    // Current lane therefore maps one-to-one with identical energies.
    let antimatter = scan_nuclear_chart(cfg);

    let matter_by_za: BTreeMap<(u16, u16), NucleusRecord> =
        matter.iter().map(|r| ((r.z, r.a), *r)).collect();
    let anti_by_za: BTreeMap<(u16, u16), NucleusRecord> =
        antimatter.iter().map(|r| ((r.z, r.a), *r)).collect();

    let mut max_abs_binding_delta_mev = 0.0_f64;
    let mut max_abs_stability_delta = 0.0_f64;
    let mut compared_rows = 0usize;
    let mut delta_csv = String::from(
        "Z,A,N,binding_delta_mev,stability_delta,fissility_delta,fission_barrier_delta_mev\n",
    );
    for (za, mr) in &matter_by_za {
        if let Some(ar) = anti_by_za.get(za) {
            compared_rows += 1;
            let d_b = ar.binding_mev - mr.binding_mev;
            let d_s = ar.stability_score - mr.stability_score;
            let d_fiss = ar.fissility - mr.fissility;
            let d_bar = ar.fission_barrier_mev - mr.fission_barrier_mev;
            max_abs_binding_delta_mev = max_abs_binding_delta_mev.max(d_b.abs());
            max_abs_stability_delta = max_abs_stability_delta.max(d_s.abs());
            delta_csv.push_str(&format!(
                "{},{},{},{:.12e},{:.12e},{:.12e},{:.12e}\n",
                mr.z, mr.a, mr.n, d_b, d_s, d_fiss, d_bar
            ));
        }
    }
    fs::write(out.join("matter_vs_antimatter_delta.csv"), delta_csv)?;

    let matter_islands =
        rank_island_candidates_with_config(&matter, IslandRankingConfig::default(), top_k);
    let antimatter_islands =
        rank_island_candidates_with_config(&antimatter, IslandRankingConfig::default(), top_k);

    let mut best_by_z_matter: BTreeMap<u16, NucleusRecord> = BTreeMap::new();
    for r in &matter {
        best_by_z_matter
            .entry(r.z)
            .and_modify(|cur| {
                if r.stability_score > cur.stability_score {
                    *cur = *r;
                }
            })
            .or_insert(*r);
    }

    let mut local_peaks: Vec<(u16, u16, f64)> = Vec::new();
    for z in cfg.z_min.max(82)..=cfg.z_max {
        if is_local_peak(&best_by_z_matter, z) {
            if let Some(r) = best_by_z_matter.get(&z) {
                local_peaks.push((r.z, r.n, r.stability_score));
            }
        }
    }

    let first_20z_supercritical = first_sustained_supercritical(&best_by_z_matter, 90, 20);
    let highest_z_best_subcritical = best_by_z_matter
        .iter()
        .filter(|(_, r)| r.fissility < 1.0)
        .map(|(z, _)| *z)
        .max();
    let highest_z_best_barrier_gt_1 = best_by_z_matter
        .iter()
        .filter(|(_, r)| r.fission_barrier_mev > 1.0)
        .map(|(z, _)| *z)
        .max();

    let strongest_114_184 = best_by_z_matter
        .get(&114)
        .copied()
        .or_else(|| matter.iter().find(|r| r.z == 114 && r.n == 184).copied());

    let mut txt = String::new();
    txt.push_str("[antimatter_periodic_report]\n");
    txt.push_str(&format!(
        "scan: Z={}..{}, N={}..{}\n",
        cfg.z_min, cfg.z_max, cfg.n_min, cfg.n_max
    ));
    txt.push_str(&format!("compared_rows = {}\n", compared_rows));
    txt.push_str(&format!(
        "max_abs_binding_delta_mev = {:.12e}\n",
        max_abs_binding_delta_mev
    ));
    txt.push_str(&format!(
        "max_abs_stability_delta = {:.12e}\n",
        max_abs_stability_delta
    ));
    txt.push_str(&format!(
        "matter_top_island = Z={},N={},score={:.6}\n",
        matter_islands.first().map(|r| r.z).unwrap_or(0),
        matter_islands.first().map(|r| r.n).unwrap_or(0),
        matter_islands
            .first()
            .map(|r| r.stability_score)
            .unwrap_or(0.0)
    ));
    txt.push_str(&format!(
        "antimatter_top_island = Z={},N={},score={:.6}\n",
        antimatter_islands.first().map(|r| r.z).unwrap_or(0),
        antimatter_islands.first().map(|r| r.n).unwrap_or(0),
        antimatter_islands
            .first()
            .map(|r| r.stability_score)
            .unwrap_or(0.0)
    ));
    txt.push_str(&format!(
        "first_20z_window_best_fissility_gt1 = {}\n",
        first_20z_supercritical
            .map(|z| z.to_string())
            .unwrap_or_else(|| "none".to_string())
    ));
    txt.push_str(&format!(
        "highest_z_best_subcritical = {}\n",
        highest_z_best_subcritical
            .map(|z| z.to_string())
            .unwrap_or_else(|| "none".to_string())
    ));
    txt.push_str(&format!(
        "highest_z_best_barrier_gt_1mev = {}\n",
        highest_z_best_barrier_gt_1
            .map(|z| z.to_string())
            .unwrap_or_else(|| "none".to_string())
    ));
    if let Some(r) = strongest_114_184 {
        txt.push_str(&format!(
            "z114_anchor = N={}, score={:.6}, fissility={:.6}, barrier_mev={:.6}\n",
            r.n, r.stability_score, r.fissility, r.fission_barrier_mev
        ));
    }
    txt.push_str("local_peaks_ge82 = ");
    if local_peaks.is_empty() {
        txt.push_str("[]\n");
    } else {
        let formatted = local_peaks
            .iter()
            .map(|(z, n, s)| format!("(Z={},N={},score={:.6})", z, n, s))
            .collect::<Vec<_>>()
            .join(", ");
        txt.push_str(&formatted);
        txt.push('\n');
    }

    fs::write(out.join("antimatter_periodic_report.txt"), &txt)?;

    let json = format!(
        concat!(
            "{{\n",
            "  \"scan\": {{\"z_min\": {}, \"z_max\": {}, \"n_min\": {}, \"n_max\": {}}},\n",
            "  \"comparisons\": {{\"rows\": {}, \"max_abs_binding_delta_mev\": {:.12e}, \"max_abs_stability_delta\": {:.12e}}},\n",
            "  \"frontier\": {{\"first_20z_window_best_fissility_gt1\": {}, \"highest_z_best_subcritical\": {}, \"highest_z_best_barrier_gt_1mev\": {}}},\n",
            "  \"matter_top_island\": {{\"z\": {}, \"n\": {}, \"score\": {:.12}}},\n",
            "  \"antimatter_top_island\": {{\"z\": {}, \"n\": {}, \"score\": {:.12}}}\n",
            "}}\n"
        ),
        cfg.z_min,
        cfg.z_max,
        cfg.n_min,
        cfg.n_max,
        compared_rows,
        max_abs_binding_delta_mev,
        max_abs_stability_delta,
        first_20z_supercritical
            .map(|z| z.to_string())
            .unwrap_or_else(|| "null".to_string()),
        highest_z_best_subcritical
            .map(|z| z.to_string())
            .unwrap_or_else(|| "null".to_string()),
        highest_z_best_barrier_gt_1
            .map(|z| z.to_string())
            .unwrap_or_else(|| "null".to_string()),
        matter_islands.first().map(|r| r.z).unwrap_or(0),
        matter_islands.first().map(|r| r.n).unwrap_or(0),
        matter_islands
            .first()
            .map(|r| r.stability_score)
            .unwrap_or(0.0),
        antimatter_islands.first().map(|r| r.z).unwrap_or(0),
        antimatter_islands.first().map(|r| r.n).unwrap_or(0),
        antimatter_islands
            .first()
            .map(|r| r.stability_score)
            .unwrap_or(0.0)
    );
    fs::write(out.join("antimatter_periodic_report.json"), json)?;

    println!("Wrote {}", out.join("antimatter_periodic_report.txt").display());
    println!(
        "Wrote {}",
        out.join("antimatter_periodic_report.json").display()
    );
    println!(
        "Wrote {}",
        out.join("matter_vs_antimatter_delta.csv").display()
    );
    println!(
        "antimatter_periodic: compared_rows={} max_abs_binding_delta={:.3e} max_abs_stability_delta={:.3e}",
        compared_rows, max_abs_binding_delta_mev, max_abs_stability_delta
    );
    Ok(())
}
