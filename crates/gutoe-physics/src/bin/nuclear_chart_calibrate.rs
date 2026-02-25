use gutoe_physics::{
    closest_to_target_island, magic_s2n_summary, rank_island_candidates_with_config, scan_nuclear_chart,
    shell_gate_metrics, write_magic_discontinuities_csv, write_magic_summary_csv, write_records_csv,
    IslandRankingConfig, ScanConfig, ShellParams,
};
use std::env;
use std::fs;
use std::path::PathBuf;

fn parse_list(name: &str, default: &[f64]) -> Vec<f64> {
    match env::var(name) {
        Ok(v) => {
            let vals: Vec<f64> = v
                .split(',')
                .map(|x| x.trim())
                .filter(|x| !x.is_empty())
                .filter_map(|x| x.parse::<f64>().ok())
                .collect();
            if vals.is_empty() {
                default.to_vec()
            } else {
                vals
            }
        }
        Err(_) => default.to_vec(),
    }
}

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

fn main() -> anyhow::Result<()> {
    let output_dir = env::var("GUTOE_NUCLEAR_CAL_OUT").unwrap_or_else(|_| "/tmp/nuclear_chart_cal".to_string());
    fs::create_dir_all(&output_dir)?;
    let out = PathBuf::from(output_dir);

    let z_min = env_u16("GUTOE_NUCLEAR_Z_MIN", 1);
    let z_max = env_u16("GUTOE_NUCLEAR_Z_MAX", 140);
    let n_min = env_u16("GUTOE_NUCLEAR_N_MIN", 1);
    let n_max = env_u16("GUTOE_NUCLEAR_N_MAX", 260);
    let target_z = env_u16("GUTOE_NUCLEAR_TARGET_Z", 114);
    let target_n = env_u16("GUTOE_NUCLEAR_TARGET_N", 184);
    let top_k = env_usize("GUTOE_NUCLEAR_TOP_K", 30);

    let amp_z_grid = parse_list("GUTOE_NUCLEAR_AMP_Z_GRID", &[1.8, 2.2, 2.8, 3.4]);
    let amp_n_grid = parse_list("GUTOE_NUCLEAR_AMP_N_GRID", &[2.2, 2.8, 3.4, 4.2]);
    let shell_amp_grid = parse_list("GUTOE_NUCLEAR_SHELL_AMP_GRID", &[12.0]);
    let shell_scale_exp_grid = parse_list("GUTOE_NUCLEAR_SHELL_SCALE_EXP_GRID", &[0.33]);
    let shell_sigma_grid = parse_list("GUTOE_NUCLEAR_SHELL_SIGMA_GRID", &[]);
    let sigma_z_grid = parse_list("GUTOE_NUCLEAR_SIGMA_Z_GRID", &[3.0, 4.0, 5.0, 6.0]);
    let sigma_n_grid = parse_list("GUTOE_NUCLEAR_SIGMA_N_GRID", &[4.0, 5.0, 6.0, 7.0]);
    let sigma_pairs: Vec<(f64, f64)> = if shell_sigma_grid.is_empty() {
        sigma_z_grid
            .iter()
            .flat_map(|&sz| sigma_n_grid.iter().map(move |&sn| (sz, sn)))
            .collect()
    } else {
        shell_sigma_grid.iter().map(|&s| (s, s)).collect()
    };
    let superheavy_proton_amp_grid = parse_list("GUTOE_NUCLEAR_SUPERHEAVY_PROTON_AMP_GRID", &[2.0]);
    let superheavy_proton_sigma_grid = parse_list("GUTOE_NUCLEAR_SUPERHEAVY_PROTON_SIGMA_GRID", &[5.0]);
    let heavy_amp_grid = parse_list("GUTOE_NUCLEAR_HEAVY_AMP_GRID", &[0.0, 1.2, 2.4, 3.6]);
    let heavy_sigma_z_grid = parse_list("GUTOE_NUCLEAR_HEAVY_SIGMA_Z_GRID", &[7.0, 9.0, 12.0]);
    let heavy_sigma_n_grid = parse_list("GUTOE_NUCLEAR_HEAVY_SIGMA_N_GRID", &[10.0, 14.0, 18.0]);
    let heavy_target_z_grid = parse_list(
        "GUTOE_NUCLEAR_HEAVY_TARGET_Z_GRID",
        &[target_z.saturating_sub(6) as f64, target_z as f64, target_z.saturating_add(6) as f64],
    );
    let heavy_target_n_grid = parse_list(
        "GUTOE_NUCLEAR_HEAVY_TARGET_N_GRID",
        &[target_n.saturating_sub(8) as f64, target_n as f64, target_n.saturating_add(8) as f64],
    );

    #[derive(Clone, Copy)]
    struct Row {
        score: f64,
        amp_z: f64,
        amp_n: f64,
        shell_amp: f64,
        shell_scale_exp: f64,
        sigma_z: f64,
        sigma_n: f64,
        superheavy_proton_amp: f64,
        superheavy_proton_sigma: f64,
        heavy_amp: f64,
        heavy_sigma_z: f64,
        heavy_sigma_n: f64,
        heavy_target_z: f64,
        heavy_target_n: f64,
        top_delta_s2n: f64,
        avg_top5_delta_s2n: f64,
        n184_delta: f64,
        proton_avg_delta_s2p: f64,
        proton_min_delta_s2p: f64,
        closest_z: u16,
        closest_n: u16,
        closest_score: f64,
        top_candidate_z: u16,
        top_candidate_n: u16,
        top_candidate_score: f64,
        top_candidate_barrier: f64,
        top_candidate_sf_log10: f64,
    }

    let mut leaderboard = String::from(
        "rank,score,amp_z,amp_n,shell_amp,shell_scale_exp,sigma_z,sigma_n,superheavy_proton_amp,superheavy_proton_sigma,heavy_amp,heavy_sigma_z,heavy_sigma_n,heavy_target_z,heavy_target_n,top_delta_s2n,avg_top5_delta_s2n,n184_delta,proton_avg_delta_s2p,proton_min_delta_s2p,closest_z,closest_n,closest_score,top_candidate_z,top_candidate_n,top_candidate_score,top_candidate_barrier,top_candidate_sf_log10\n",
    );
    let mut rows: Vec<Row> = Vec::new();
    let mut best_score = f64::NEG_INFINITY;
    let mut best_records = Vec::new();
    let mut best_ranked = Vec::new();
    let mut best_metrics = None;
    let mut best_params = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    for &amp_z in &amp_z_grid {
        for &amp_n in &amp_n_grid {
            for &shell_amp in &shell_amp_grid {
                for &shell_scale_exp in &shell_scale_exp_grid {
                    for &(sigma_z, sigma_n) in &sigma_pairs {
                            for &superheavy_proton_amp in &superheavy_proton_amp_grid {
                                for &superheavy_proton_sigma in &superheavy_proton_sigma_grid {
                                    for &heavy_amp in &heavy_amp_grid {
                                        for &heavy_sigma_z in &heavy_sigma_z_grid {
                                            for &heavy_sigma_n in &heavy_sigma_n_grid {
                                                for &heavy_target_z_f in &heavy_target_z_grid {
                                                    for &heavy_target_n_f in &heavy_target_n_grid {
                                                        let cfg = ScanConfig {
                                                            z_min,
                                                            z_max,
                                                            n_min,
                                                            n_max,
                                                            shell: ShellParams {
                                                                amplitude_z: amp_z,
                                                                amplitude_n: amp_n,
                                                                shell_amp,
                                                                shell_scale_exp,
                                                                sigma_z,
                                                                sigma_n,
                                                                superheavy_proton_amplitude: superheavy_proton_amp,
                                                                superheavy_proton_sigma,
                                                                heavy_target_z: heavy_target_z_f,
                                                                heavy_target_n: heavy_target_n_f,
                                                                heavy_sigma_z,
                                                                heavy_sigma_n,
                                                                heavy_amplitude: heavy_amp,
                                                                ..ShellParams::default()
                                                            },
                                                            ..ScanConfig::default()
                                                        };
                                                let records = scan_nuclear_chart(cfg);
                                                let metrics = shell_gate_metrics(&records);
                                                let ranking_cfg = IslandRankingConfig {
                                                    min_z: 104,
                                                    target_z,
                                                    target_n,
                                                    ..IslandRankingConfig::default()
                                                };
                                                let ranked =
                                                    rank_island_candidates_with_config(&records, ranking_cfg, top_k);
                                                let closest =
                                                    closest_to_target_island(&records, target_z, target_n);
                                                let top = ranked.first().copied();

                                                let top_score =
                                                    top.map(|r| r.stability_score).unwrap_or(-1.0e9);
                                                let closest_score =
                                                    closest.map(|r| r.stability_score).unwrap_or(-1.0e9);
                                                let objective = 0.42 * metrics.strongest_n184_delta_s2n_mev
                                                    + 0.14 * metrics.avg_top5_delta_s2n_mev
                                                    + 0.20 * metrics.avg_superheavy_proton_delta_s2p_mev
                                                    + 0.10 * metrics.min_superheavy_proton_delta_s2p_mev
                                                    + 0.08 * top_score
                                                    + 0.06 * closest_score;

                                                let (top_candidate_z, top_candidate_n) =
                                                    top.map(|r| (r.z, r.n)).unwrap_or((0, 0));
                                                let top_candidate_barrier =
                                                    top.map(|r| r.fission_barrier_mev).unwrap_or(0.0);
                                                let top_candidate_sf_log10 =
                                                    top.map(|r| r.sf_log10_half_life_s).unwrap_or(-1.0e9);
                                                let (closest_z, closest_n) =
                                                    closest.map(|r| (r.z, r.n)).unwrap_or((0, 0));
                                                        rows.push(Row {
                                                            score: objective,
                                                            amp_z,
                                                            amp_n,
                                                            shell_amp,
                                                            shell_scale_exp,
                                                            sigma_z,
                                                            sigma_n,
                                                            superheavy_proton_amp,
                                                            superheavy_proton_sigma,
                                                            heavy_amp,
                                                            heavy_sigma_z,
                                                            heavy_sigma_n,
                                                            heavy_target_z: heavy_target_z_f,
                                                            heavy_target_n: heavy_target_n_f,
                                                            top_delta_s2n: metrics.top_delta_s2n_mev,
                                                            avg_top5_delta_s2n: metrics.avg_top5_delta_s2n_mev,
                                                            n184_delta: metrics.strongest_n184_delta_s2n_mev,
                                                            proton_avg_delta_s2p: metrics.avg_superheavy_proton_delta_s2p_mev,
                                                            proton_min_delta_s2p: metrics.min_superheavy_proton_delta_s2p_mev,
                                                            closest_z,
                                                            closest_n,
                                                            closest_score,
                                                            top_candidate_z,
                                                            top_candidate_n,
                                                            top_candidate_score: top_score,
                                                            top_candidate_barrier,
                                                            top_candidate_sf_log10,
                                                        });

                                                        if objective > best_score {
                                                            best_score = objective;
                                                            best_records = records;
                                                            best_ranked = ranked;
                                                            best_metrics = Some(metrics);
                                                            best_params = (
                                                                amp_z,
                                                                amp_n,
                                                                shell_amp,
                                                                shell_scale_exp,
                                                                sigma_z,
                                                                sigma_n,
                                                                superheavy_proton_amp,
                                                                superheavy_proton_sigma,
                                                                heavy_amp,
                                                                heavy_sigma_z,
                                                                heavy_sigma_n,
                                                                heavy_target_z_f,
                                                                heavy_target_n_f,
                                                            );
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                    }
                }
            }
        }
    }

    rows.sort_by(|a, b| b.score.total_cmp(&a.score));

    for (idx, row) in rows.iter().enumerate() {
        leaderboard.push_str(&format!(
            "{},{:.6},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.6},{:.6},{:.6},{:.6},{:.6},{},{},{:.6},{},{},{:.6},{:.6},{:.6}\n",
            idx + 1,
            row.score,
            row.amp_z,
            row.amp_n,
            row.shell_amp,
            row.shell_scale_exp,
            row.sigma_z,
            row.sigma_n,
            row.superheavy_proton_amp,
            row.superheavy_proton_sigma,
            row.heavy_amp,
            row.heavy_sigma_z,
            row.heavy_sigma_n,
            row.heavy_target_z,
            row.heavy_target_n,
            row.top_delta_s2n,
            row.avg_top5_delta_s2n,
            row.n184_delta,
            row.proton_avg_delta_s2p,
            row.proton_min_delta_s2p,
            row.closest_z,
            row.closest_n,
            row.closest_score,
            row.top_candidate_z,
            row.top_candidate_n,
            row.top_candidate_score,
            row.top_candidate_barrier,
            row.top_candidate_sf_log10
        ));
    }

    let best_magic = gutoe_physics::magic_s2n_discontinuities(&best_records, top_k);
    let best_magic_summary = magic_s2n_summary(&best_records);
    let (
        amp_z,
        amp_n,
        shell_amp,
        shell_scale_exp,
        sigma_z,
        sigma_n,
        superheavy_proton_amp,
        superheavy_proton_sigma,
        heavy_amp,
        heavy_sigma_z,
        heavy_sigma_n,
        heavy_target_z,
        heavy_target_n,
    ) =
        best_params;
    let metrics = best_metrics.unwrap_or(gutoe_physics::ShellGateMetrics {
        top_delta_s2n_mev: 0.0,
        avg_top5_delta_s2n_mev: 0.0,
        strongest_n184_delta_s2n_mev: 0.0,
        strongest_superheavy_proton_delta_s2p_mev: 0.0,
        avg_superheavy_proton_delta_s2p_mev: 0.0,
        min_superheavy_proton_delta_s2p_mev: 0.0,
    });
    let best_desc = format!(
        "best_score={:.6}\nbest_shell=amp_z:{:.3},amp_n:{:.3},shell_amp:{:.3},shell_scale_exp:{:.3},sigma_z:{:.3},sigma_n:{:.3},superheavy_proton_amp:{:.3},superheavy_proton_sigma:{:.3},heavy_amp:{:.3},heavy_sigma_z:{:.3},heavy_sigma_n:{:.3},heavy_target_z:{:.3},heavy_target_n:{:.3}\ntop_delta_s2n={:.6}\navg_top5_delta_s2n={:.6}\nn184_delta={:.6}\nproton_avg_delta_s2p={:.6}\nproton_min_delta_s2p={:.6}\n",
        best_score,
        amp_z,
        amp_n,
        shell_amp,
        shell_scale_exp,
        sigma_z,
        sigma_n,
        superheavy_proton_amp,
        superheavy_proton_sigma,
        heavy_amp,
        heavy_sigma_z,
        heavy_sigma_n,
        heavy_target_z,
        heavy_target_n,
        metrics.top_delta_s2n_mev,
        metrics.avg_top5_delta_s2n_mev,
        metrics.strongest_n184_delta_s2n_mev,
        metrics.avg_superheavy_proton_delta_s2p_mev,
        metrics.min_superheavy_proton_delta_s2p_mev
    );
    let mut top_csv = String::from(
        "rank,Z,N,A,binding_per_nucleon_mev,s2n_mev,s2p_mev,fissility,fission_barrier_mev,sf_log10_half_life_s,stability_score\n",
    );
    for (idx, r) in best_ranked.iter().enumerate() {
        top_csv.push_str(&format!(
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
    fs::write(out.join("best_top_islands.csv"), top_csv)?;

    write_records_csv(out.join("best_nuclides.csv"), &best_records)?;
    write_magic_discontinuities_csv(out.join("best_magic_discontinuities.csv"), &best_magic)?;
    write_magic_summary_csv(out.join("best_magic_summary.csv"), &best_magic_summary)?;
    fs::write(out.join("leaderboard.csv"), leaderboard)?;
    fs::write(
        out.join("summary.txt"),
        format!(
            "target=({}, {})\n{}\nobjective_top={:.6}\nrows={}\n",
            target_z, target_n, best_desc, best_score, best_records.len()
        ),
    )?;

    println!("Calibration output:");
    println!("  {}", out.join("summary.txt").display());
    println!("  {}", out.join("leaderboard.csv").display());
    println!("  {}", out.join("best_top_islands.csv").display());
    println!("  {}", out.join("best_nuclides.csv").display());
    println!("  {}", out.join("best_magic_discontinuities.csv").display());
    println!("  {}", out.join("best_magic_summary.csv").display());
    Ok(())
}
