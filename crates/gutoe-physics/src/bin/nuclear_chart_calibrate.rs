use gutoe_physics::{
    closest_to_target_island, magic_s2n_summary, rank_island_candidates_with_config,
    scan_nuclear_chart, shell_gate_metrics, write_magic_discontinuities_csv,
    write_magic_summary_csv, write_records_csv, IslandRankingConfig, MagicSummaryRow, ScanConfig,
    ShellParams, MONITORED_SUPERHEAVY_PROTON_CLOSURES,
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

fn parse_bool_list(name: &str, default: &[bool]) -> Vec<bool> {
    match env::var(name) {
        Ok(v) => {
            let vals: Vec<bool> = v
                .split(',')
                .map(|x| x.trim().to_ascii_lowercase())
                .filter(|x| !x.is_empty())
                .filter_map(|x| match x.as_str() {
                    "1" | "true" | "yes" | "on" => Some(true),
                    "0" | "false" | "no" | "off" => Some(false),
                    _ => None,
                })
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

fn reference_shell_gap_bounds_mev(magic_n: u16) -> Option<(f64, f64)> {
    match magic_n {
        50 => Some((5.0, 8.0)),
        82 => Some((4.0, 6.5)),
        126 => Some((3.0, 5.5)),
        _ => None,
    }
}

fn strongest_ratio_for_magic(summary: &[MagicSummaryRow], magic_n: u16) -> f64 {
    let Some((ref_min, ref_max)) = reference_shell_gap_bounds_mev(magic_n) else {
        return 0.0;
    };
    let ref_mid = 0.5 * (ref_min + ref_max);
    summary
        .iter()
        .find(|row| row.magic_n == magic_n)
        .map(|row| row.strongest_delta_s2n_mev / ref_mid)
        .unwrap_or(0.0)
}

fn band_penalty(value: f64, low: f64, high: f64) -> f64 {
    if value < low {
        (low - value).powi(2)
    } else if value > high {
        (value - high).powi(2)
    } else {
        0.0
    }
}

fn main() -> anyhow::Result<()> {
    let output_dir =
        env::var("GUTOE_NUCLEAR_CAL_OUT").unwrap_or_else(|_| "/tmp/nuclear_chart_cal".to_string());
    fs::create_dir_all(&output_dir)?;
    let out = PathBuf::from(output_dir);

    let z_min = env_u16("GUTOE_NUCLEAR_Z_MIN", 1);
    let z_max = env_u16("GUTOE_NUCLEAR_Z_MAX", 140);
    let n_min = env_u16("GUTOE_NUCLEAR_N_MIN", 1);
    let n_max = env_u16("GUTOE_NUCLEAR_N_MAX", 260);
    let target_z = env_u16("GUTOE_NUCLEAR_TARGET_Z", 114);
    let target_n = env_u16("GUTOE_NUCLEAR_TARGET_N", 184);
    let top_k = env_usize("GUTOE_NUCLEAR_TOP_K", 30);
    let ratio_band_low = env::var("GUTOE_NUCLEAR_RATIO_BAND_LOW")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(0.80);
    let ratio_band_high = env::var("GUTOE_NUCLEAR_RATIO_BAND_HIGH")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(1.20);
    let ratio_penalty_weight = env::var("GUTOE_NUCLEAR_RATIO_BAND_WEIGHT")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(12.0);
    let default_shell = ShellParams::default();

    let amp_z_grid = parse_list("GUTOE_NUCLEAR_AMP_Z_GRID", &[1.8, 2.2, 2.8, 3.4]);
    let amp_n_grid = parse_list("GUTOE_NUCLEAR_AMP_N_GRID", &[2.2, 2.8, 3.4, 4.2]);
    let shell_amp_grid = parse_list("GUTOE_NUCLEAR_SHELL_AMP_GRID", &[12.0]);
    let shell_scale_exp_grid = parse_list("GUTOE_NUCLEAR_SHELL_SCALE_EXP_GRID", &[0.33]);
    let shell_sigma_grid = parse_list("GUTOE_NUCLEAR_SHELL_SIGMA_GRID", &[]);
    let use_strutinsky_grid = parse_bool_list(
        "GUTOE_NUCLEAR_USE_STRUTINSKY_GRID",
        &[default_shell.use_strutinsky],
    );
    let strutinsky_mix_grid = parse_list(
        "GUTOE_NUCLEAR_STRUTINSKY_MIX_GRID",
        &[default_shell.strutinsky_mix],
    );
    let strutinsky_gamma_grid = parse_list(
        "GUTOE_NUCLEAR_STRUTINSKY_GAMMA_GRID",
        &[default_shell.strutinsky_gamma],
    );
    let strutinsky_spacing_grid = parse_list(
        "GUTOE_NUCLEAR_STRUTINSKY_SPACING_GRID",
        &[default_shell.strutinsky_spacing_mev],
    );
    let strutinsky_spin_orbit_grid = parse_list(
        "GUTOE_NUCLEAR_STRUTINSKY_SPIN_ORBIT_GRID",
        &[default_shell.strutinsky_spin_orbit_mev],
    );
    let strutinsky_coulomb_shift_grid = parse_list(
        "GUTOE_NUCLEAR_STRUTINSKY_COULOMB_SHIFT_GRID",
        &[default_shell.strutinsky_coulomb_shift_mev],
    );
    let strutinsky_ws_depth_grid = parse_list(
        "GUTOE_NUCLEAR_STRUTINSKY_WS_DEPTH_GRID",
        &[default_shell.strutinsky_ws_depth_mev],
    );
    let strutinsky_ws_diffuseness_grid = parse_list(
        "GUTOE_NUCLEAR_STRUTINSKY_WS_DIFFUSENESS_GRID",
        &[default_shell.strutinsky_ws_diffuseness_fm],
    );
    let strutinsky_ws_a_ref_grid = parse_list(
        "GUTOE_NUCLEAR_STRUTINSKY_WS_A_REF_GRID",
        &[default_shell.strutinsky_ws_a_ref],
    );
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
    let superheavy_proton_sigma_grid =
        parse_list("GUTOE_NUCLEAR_SUPERHEAVY_PROTON_SIGMA_GRID", &[5.0]);
    let heavy_amp_grid = parse_list("GUTOE_NUCLEAR_HEAVY_AMP_GRID", &[0.0, 1.2, 2.4, 3.6]);
    let heavy_sigma_z_grid = parse_list("GUTOE_NUCLEAR_HEAVY_SIGMA_Z_GRID", &[7.0, 9.0, 12.0]);
    let heavy_sigma_n_grid = parse_list("GUTOE_NUCLEAR_HEAVY_SIGMA_N_GRID", &[10.0, 14.0, 18.0]);
    let heavy_target_z_grid = parse_list(
        "GUTOE_NUCLEAR_HEAVY_TARGET_Z_GRID",
        &[
            target_z.saturating_sub(6) as f64,
            target_z as f64,
            target_z.saturating_add(6) as f64,
        ],
    );
    let heavy_target_n_grid = parse_list(
        "GUTOE_NUCLEAR_HEAVY_TARGET_N_GRID",
        &[
            target_n.saturating_sub(8) as f64,
            target_n as f64,
            target_n.saturating_add(8) as f64,
        ],
    );

    #[derive(Clone, Copy)]
    struct Row {
        score: f64,
        amp_z: f64,
        amp_n: f64,
        shell_amp: f64,
        shell_scale_exp: f64,
        use_strutinsky: bool,
        strutinsky_mix: f64,
        strutinsky_gamma: f64,
        strutinsky_spacing_mev: f64,
        strutinsky_spin_orbit_mev: f64,
        strutinsky_coulomb_shift_mev: f64,
        strutinsky_ws_depth_mev: f64,
        strutinsky_ws_diffuseness_fm: f64,
        strutinsky_ws_a_ref: f64,
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
        proton_monitored_avg_delta_s2p: f64,
        proton_monitored_min_delta_s2p: f64,
        n50_ratio: f64,
        n82_ratio: f64,
        n126_ratio: f64,
        ratio_penalty: f64,
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
        "rank,score,amp_z,amp_n,shell_amp,shell_scale_exp,use_strutinsky,strutinsky_mix,strutinsky_gamma,strutinsky_spacing_mev,strutinsky_spin_orbit_mev,strutinsky_coulomb_shift_mev,strutinsky_ws_depth_mev,strutinsky_ws_diffuseness_fm,strutinsky_ws_a_ref,sigma_z,sigma_n,superheavy_proton_amp,superheavy_proton_sigma,heavy_amp,heavy_sigma_z,heavy_sigma_n,heavy_target_z,heavy_target_n,top_delta_s2n,avg_top5_delta_s2n,n184_delta,proton_avg_delta_s2p,proton_min_delta_s2p,proton_monitored_avg_delta_s2p,proton_monitored_min_delta_s2p,n50_ratio,n82_ratio,n126_ratio,ratio_penalty,closest_z,closest_n,closest_score,top_candidate_z,top_candidate_n,top_candidate_score,top_candidate_barrier,top_candidate_sf_log10\n",
    );
    let mut rows: Vec<Row> = Vec::new();
    let mut best_score = f64::NEG_INFINITY;
    let mut best_records = Vec::new();
    let mut best_ranked = Vec::new();
    let mut best_metrics = None;
    let mut best_shell = ShellParams::default();

    for &amp_z in &amp_z_grid {
        for &amp_n in &amp_n_grid {
            for &shell_amp in &shell_amp_grid {
                for &shell_scale_exp in &shell_scale_exp_grid {
                    for &use_strutinsky in &use_strutinsky_grid {
                        for &strutinsky_mix in &strutinsky_mix_grid {
                            for &strutinsky_gamma in &strutinsky_gamma_grid {
                                for &strutinsky_spacing_mev in &strutinsky_spacing_grid {
                                    for &strutinsky_spin_orbit_mev in &strutinsky_spin_orbit_grid {
                                        for &strutinsky_coulomb_shift_mev in
                                            &strutinsky_coulomb_shift_grid
                                        {
                                            for &strutinsky_ws_depth_mev in
                                                &strutinsky_ws_depth_grid
                                            {
                                                for &strutinsky_ws_diffuseness_fm in
                                                    &strutinsky_ws_diffuseness_grid
                                                {
                                                    for &strutinsky_ws_a_ref in
                                                        &strutinsky_ws_a_ref_grid
                                                    {
                                                        for &(sigma_z, sigma_n) in &sigma_pairs {
                                                            for &superheavy_proton_amp in
                                                                &superheavy_proton_amp_grid
                                                            {
                                                                for &superheavy_proton_sigma in
                                                                    &superheavy_proton_sigma_grid
                                                                {
                                                                    for &heavy_amp in
                                                                        &heavy_amp_grid
                                                                    {
                                                                        for &heavy_sigma_z in
                                                                            &heavy_sigma_z_grid
                                                                        {
                                                                            for &heavy_sigma_n in
                                                                                &heavy_sigma_n_grid
                                                                            {
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
                                                                                    use_strutinsky,
                                                                                    strutinsky_mix,
                                                                                    strutinsky_gamma,
                                                                                    strutinsky_spacing_mev,
                                                                                    strutinsky_spin_orbit_mev,
                                                                                    strutinsky_coulomb_shift_mev,
                                                                                    strutinsky_ws_depth_mev,
                                                                                    strutinsky_ws_diffuseness_fm,
                                                                                    strutinsky_ws_a_ref,
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
                                                                            let summary = magic_s2n_summary(&records);
                                                                            let n50_ratio = strongest_ratio_for_magic(&summary, 50);
                                                                            let n82_ratio = strongest_ratio_for_magic(&summary, 82);
                                                                            let n126_ratio = strongest_ratio_for_magic(&summary, 126);
                                                                            let ratio_penalty = band_penalty(
                                                                                n50_ratio,
                                                                                ratio_band_low,
                                                                                ratio_band_high,
                                                                            ) + band_penalty(
                                                                                n82_ratio,
                                                                                ratio_band_low,
                                                                                ratio_band_high,
                                                                            ) + band_penalty(
                                                                                n126_ratio,
                                                                                ratio_band_low,
                                                                                ratio_band_high,
                                                                            );
                                                                            let ranking_cfg = IslandRankingConfig {
                                                                                min_z: 104,
                                                                                target_z,
                                                                                target_n,
                                                                                ..IslandRankingConfig::default()
                                                                            };
                                                                            let ranked = rank_island_candidates_with_config(
                                                                                &records,
                                                                                ranking_cfg,
                                                                                top_k,
                                                                            );
                                                                            let closest = closest_to_target_island(
                                                                                &records, target_z, target_n,
                                                                            );
                                                                            let top = ranked.first().copied();

                                                                            let top_score = top
                                                                                .map(|r| r.stability_score)
                                                                                .unwrap_or(-1.0e9);
                                                                            let closest_score = closest
                                                                                .map(|r| r.stability_score)
                                                                                .unwrap_or(-1.0e9);
                                                                            let objective_base = 0.40 * metrics.strongest_n184_delta_s2n_mev
                                                                                + 0.14 * metrics.avg_top5_delta_s2n_mev
                                                                                + 0.16 * metrics.avg_monitored_proton_delta_s2p_mev
                                                                                + 0.12 * metrics.min_monitored_proton_delta_s2p_mev
                                                                                + 0.08 * metrics.avg_superheavy_proton_delta_s2p_mev
                                                                                + 0.04 * metrics.min_superheavy_proton_delta_s2p_mev
                                                                                + 0.04 * top_score
                                                                                + 0.02 * closest_score;
                                                                            let objective = objective_base
                                                                                - ratio_penalty_weight * ratio_penalty;

                                                                            let (top_candidate_z, top_candidate_n) =
                                                                                top.map(|r| (r.z, r.n)).unwrap_or((0, 0));
                                                                            let top_candidate_barrier = top
                                                                                .map(|r| r.fission_barrier_mev)
                                                                                .unwrap_or(0.0);
                                                                            let top_candidate_sf_log10 = top
                                                                                .map(|r| r.sf_log10_half_life_s)
                                                                                .unwrap_or(-1.0e9);
                                                                            let (closest_z, closest_n) =
                                                                                closest.map(|r| (r.z, r.n)).unwrap_or((0, 0));

                                                                            rows.push(Row {
                                                                                score: objective,
                                                                                amp_z,
                                                                                amp_n,
                                                                                shell_amp,
                                                                                shell_scale_exp,
                                                                                use_strutinsky,
                                                                                strutinsky_mix,
                                                                                strutinsky_gamma,
                                                                                strutinsky_spacing_mev,
                                                                                strutinsky_spin_orbit_mev,
                                                                                strutinsky_coulomb_shift_mev,
                                                                                strutinsky_ws_depth_mev,
                                                                                strutinsky_ws_diffuseness_fm,
                                                                                strutinsky_ws_a_ref,
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
                                                                                proton_monitored_avg_delta_s2p: metrics.avg_monitored_proton_delta_s2p_mev,
                                                                                proton_monitored_min_delta_s2p: metrics.min_monitored_proton_delta_s2p_mev,
                                                                                n50_ratio,
                                                                                n82_ratio,
                                                                                n126_ratio,
                                                                                ratio_penalty,
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
                                                                                best_shell = cfg.shell;
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
            "{},{:.6},{:.3},{:.3},{:.3},{:.3},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{},{},{:.6},{},{},{:.6},{:.6},{:.6}\n",
            idx + 1,
            row.score,
            row.amp_z,
            row.amp_n,
            row.shell_amp,
            row.shell_scale_exp,
            row.use_strutinsky,
            row.strutinsky_mix,
            row.strutinsky_gamma,
            row.strutinsky_spacing_mev,
            row.strutinsky_spin_orbit_mev,
            row.strutinsky_coulomb_shift_mev,
            row.strutinsky_ws_depth_mev,
            row.strutinsky_ws_diffuseness_fm,
            row.strutinsky_ws_a_ref,
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
            row.proton_monitored_avg_delta_s2p,
            row.proton_monitored_min_delta_s2p,
            row.n50_ratio,
            row.n82_ratio,
            row.n126_ratio,
            row.ratio_penalty,
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
    let metrics = best_metrics.unwrap_or(gutoe_physics::ShellGateMetrics {
        top_delta_s2n_mev: 0.0,
        avg_top5_delta_s2n_mev: 0.0,
        strongest_n184_delta_s2n_mev: 0.0,
        strongest_superheavy_proton_delta_s2p_mev: 0.0,
        avg_superheavy_proton_delta_s2p_mev: 0.0,
        min_superheavy_proton_delta_s2p_mev: 0.0,
        strongest_monitored_proton_delta_s2p_mev: 0.0,
        avg_monitored_proton_delta_s2p_mev: 0.0,
        min_monitored_proton_delta_s2p_mev: 0.0,
    });
    let best_desc = format!(
        "best_score={:.6}\nbest_shell=amp_z:{:.3},amp_n:{:.3},shell_amp:{:.3},shell_scale_exp:{:.3},use_strutinsky:{},strutinsky_mix:{:.3},strutinsky_gamma:{:.3},strutinsky_spacing:{:.3},strutinsky_spin_orbit:{:.3},strutinsky_coulomb_shift:{:.3},strutinsky_ws_depth:{:.3},strutinsky_ws_diffuseness:{:.3},strutinsky_ws_a_ref:{:.3},sigma_z:{:.3},sigma_n:{:.3},superheavy_proton_amp:{:.3},superheavy_proton_sigma:{:.3},heavy_amp:{:.3},heavy_sigma_z:{:.3},heavy_sigma_n:{:.3},heavy_target_z:{:.3},heavy_target_n:{:.3}\nmonitored_proton_closures={:?}\ntop_delta_s2n={:.6}\navg_top5_delta_s2n={:.6}\nn184_delta={:.6}\nproton_avg_delta_s2p={:.6}\nproton_min_delta_s2p={:.6}\nproton_monitored_avg_delta_s2p={:.6}\nproton_monitored_min_delta_s2p={:.6}\n",
        best_score,
        best_shell.amplitude_z,
        best_shell.amplitude_n,
        best_shell.shell_amp,
        best_shell.shell_scale_exp,
        best_shell.use_strutinsky,
        best_shell.strutinsky_mix,
        best_shell.strutinsky_gamma,
        best_shell.strutinsky_spacing_mev,
        best_shell.strutinsky_spin_orbit_mev,
        best_shell.strutinsky_coulomb_shift_mev,
        best_shell.strutinsky_ws_depth_mev,
        best_shell.strutinsky_ws_diffuseness_fm,
        best_shell.strutinsky_ws_a_ref,
        best_shell.sigma_z,
        best_shell.sigma_n,
        best_shell.superheavy_proton_amplitude,
        best_shell.superheavy_proton_sigma,
        best_shell.heavy_amplitude,
        best_shell.heavy_sigma_z,
        best_shell.heavy_sigma_n,
        best_shell.heavy_target_z,
        best_shell.heavy_target_n,
        MONITORED_SUPERHEAVY_PROTON_CLOSURES,
        metrics.top_delta_s2n_mev,
        metrics.avg_top5_delta_s2n_mev,
        metrics.strongest_n184_delta_s2n_mev,
        metrics.avg_superheavy_proton_delta_s2p_mev,
        metrics.min_superheavy_proton_delta_s2p_mev,
        metrics.avg_monitored_proton_delta_s2p_mev,
        metrics.min_monitored_proton_delta_s2p_mev
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
            target_z,
            target_n,
            best_desc,
            best_score,
            best_records.len()
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
