use gutoe_physics::{
    closest_to_target_island, derived_superheavy_proton_candidates, magic_s2n_summary, proton_s2p_summary,
    rank_island_candidates_with_config, scan_nuclear_chart, IslandRankingConfig, NucleusRecord, ScanConfig,
    StandardModelDynamicsMap,
};
use std::collections::BTreeMap;
use std::env;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

const ELECTRON_MASS_MEV_OBS: f64 = 0.510_998_950;
const PROTON_MASS_MEV_OBS: f64 = 938.272_088_16;
const NEUTRON_MASS_MEV_OBS: f64 = 939.565_420_52;
const BETA_MASS_COEFF_Z_MEV: f64 = (PROTON_MASS_MEV_OBS + ELECTRON_MASS_MEV_OBS) - NEUTRON_MASS_MEV_OBS;

fn reference_shell_gap_bounds_mev(magic_n: u16) -> Option<(f64, f64)> {
    // Broad experimental windows (MeV) used only for attenuation diagnostics.
    // This is for calibration visibility, not parameter fitting.
    match magic_n {
        8 => Some((10.0, 14.0)),
        20 => Some((5.0, 8.0)),
        28 => Some((4.0, 6.0)),
        50 => Some((5.0, 8.0)),
        82 => Some((4.0, 6.5)),
        126 => Some((3.0, 5.5)),
        _ => None,
    }
}

fn triangular(n: u32) -> u32 {
    n * (n + 1) / 2
}

fn observed_stable_isotope_count(z: u16) -> Option<u16> {
    // Stable-isotope count reference (strictly stable nuclei; no long-lived radioisotopes).
    match z {
        1 => Some(2),
        2 => Some(2),
        3 => Some(2),
        4 => Some(1),
        5 => Some(2),
        6 => Some(2),
        7 => Some(2),
        8 => Some(3),
        9 => Some(1),
        10 => Some(3),
        11 => Some(1),
        12 => Some(3),
        13 => Some(1),
        14 => Some(3),
        15 => Some(1),
        16 => Some(4),
        17 => Some(2),
        18 => Some(3),
        19 => Some(2),
        20 => Some(6),
        21 => Some(1),
        22 => Some(5),
        23 => Some(1),
        24 => Some(4),
        25 => Some(1),
        26 => Some(4),
        27 => Some(1),
        28 => Some(5),
        29 => Some(2),
        30 => Some(5),
        31 => Some(2),
        32 => Some(5),
        33 => Some(1),
        34 => Some(6),
        35 => Some(2),
        36 => Some(6),
        37 => Some(2),
        38 => Some(4),
        39 => Some(1),
        40 => Some(5),
        41 => Some(1),
        42 => Some(7),
        43 => Some(0),
        44 => Some(7),
        45 => Some(1),
        46 => Some(6),
        47 => Some(2),
        48 => Some(8),
        49 => Some(2),
        50 => Some(10),
        51 => Some(2),
        52 => Some(8),
        53 => Some(1),
        54 => Some(9),
        55 => Some(1),
        56 => Some(7),
        57 => Some(1),
        58 => Some(4),
        59 => Some(1),
        60 => Some(7),
        61 => Some(0),
        62 => Some(7),
        63 => Some(2),
        64 => Some(7),
        65 => Some(1),
        66 => Some(7),
        67 => Some(1),
        68 => Some(6),
        69 => Some(1),
        70 => Some(7),
        71 => Some(1),
        72 => Some(6),
        73 => Some(1),
        74 => Some(5),
        75 => Some(0),
        76 => Some(7),
        77 => Some(2),
        78 => Some(6),
        79 => Some(1),
        80 => Some(7),
        81 => Some(2),
        82 => Some(4),
        83 => Some(0),
        84 => Some(0),
        85 => Some(0),
        86 => Some(0),
        87 => Some(0),
        88 => Some(0),
        89 => Some(0),
        90 => Some(0),
        91 => Some(0),
        92 => Some(0),
        93 => Some(0),
        94 => Some(0),
        _ => None,
    }
}

fn observed_stable_mass_numbers_for_z(z: u16) -> Option<&'static [u16]> {
    match z {
        // Tin has the most stable isotopes; strict reference set.
        50 => Some(&[112, 114, 115, 116, 117, 118, 119, 120, 122, 124]),
        _ => None,
    }
}

fn classify_long_lived(r: &NucleusRecord, beta_stable_local_min: bool) -> (bool, bool, bool, bool, bool, bool) {
    let fail_beta_optimal = !beta_stable_local_min;
    let fail_fissility = r.fissility > 1.0;
    let fail_s2n = if r.n <= 2 {
        false
    } else {
        !r.s2n_mev.map(|v| v > 0.0).unwrap_or(false)
    };
    let fail_s2p = if r.z <= 2 {
        false
    } else {
        !r.s2p_mev.map(|v| v > 0.0).unwrap_or(false)
    };
    let fail_sf = r.z > 82 && r.sf_log10_half_life_s < 20.0;
    let predicted = !(fail_beta_optimal || fail_fissility || fail_s2n || fail_s2p || fail_sf);
    (
        predicted,
        fail_beta_optimal,
        fail_fissility,
        fail_s2n,
        fail_s2p,
        fail_sf,
    )
}

fn build_beta_local_min_map(records: &[NucleusRecord]) -> BTreeMap<(u16, u16), bool> {
    let mut mass_proxy_by_az: BTreeMap<(u16, u16), f64> = BTreeMap::new();
    for r in records {
        // Atomic mass at fixed A differs by Z * ((m_p + m_e) - m_n) - B(Z,N).
        // Local minima of this proxy correspond to beta-stable isobars.
        let mass_proxy = BETA_MASS_COEFF_Z_MEV * r.z as f64 - r.binding_mev;
        mass_proxy_by_az.insert((r.a, r.z), mass_proxy);
    }

    let mut out: BTreeMap<(u16, u16), bool> = BTreeMap::new();
    for r in records {
        let Some(&m0) = mass_proxy_by_az.get(&(r.a, r.z)) else {
            out.insert((r.z, r.n), false);
            continue;
        };
        let left_ok = if r.z > 1 {
            mass_proxy_by_az
                .get(&(r.a, r.z - 1))
                .map(|&ml| m0 <= ml + 1e-9)
                .unwrap_or(true)
        } else {
            true
        };
        let right_ok = mass_proxy_by_az
            .get(&(r.a, r.z + 1))
            .map(|&mr| m0 <= mr + 1e-9)
            .unwrap_or(true);
        out.insert((r.z, r.n), left_ok && right_ok);
    }
    out
}

fn now_unix_seconds() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

fn main() -> anyhow::Result<()> {
    let out_dir = env::var("GUTOE_MASS_PERIODIC_OUT").unwrap_or_else(|_| "/tmp/nuclear_chart".to_string());
    fs::create_dir_all(&out_dir)?;
    let out = PathBuf::from(out_dir);

    let sm = StandardModelDynamicsMap::from_clifford_z3();
    let alpha_inv_struct = (1.0 / sm.alpha_leading_order).round() as u32;
    let mp_me_struct = sm.total_gauge_generators * triangular(sm.clifford_dim + 1);
    let proton_pred_from_e = ELECTRON_MASS_MEV_OBS * mp_me_struct as f64;
    let electron_pred_from_p = PROTON_MASS_MEV_OBS / mp_me_struct as f64;
    let neutron_minus_proton_struct = sm.lambda_qg * sm.total_gauge_generators as f64;
    let neutron_pred = proton_pred_from_e + neutron_minus_proton_struct;

    let cfg = ScanConfig::default();
    let records = scan_nuclear_chart(cfg);
    let ranked = rank_island_candidates_with_config(
        &records,
        IslandRankingConfig {
            target_z: 114,
            target_n: 184,
            ..IslandRankingConfig::default()
        },
        40,
    );
    let beta_local_min = build_beta_local_min_map(&records);
    let stable_like: Vec<_> = records
        .iter()
        .filter(|r| {
            let beta_ok = beta_local_min.get(&(r.z, r.n)).copied().unwrap_or(false);
            classify_long_lived(r, beta_ok).0
        })
        .collect();
    let valley: Vec<_> = records.iter().filter(|r| r.beta_optimal_for_a).collect();
    let closest_target = closest_to_target_island(&records, 114, 184);
    let top = ranked.first().copied();

    let mut isotopes_per_z: BTreeMap<u16, usize> = BTreeMap::new();
    for r in &stable_like {
        *isotopes_per_z.entry(r.z).or_insert(0) += 1;
    }

    // Tin diagnostics (magic-proton showcase): compare exact stable A-set.
    let mut tin_predicted_a: Vec<u16> = stable_like
        .iter()
        .filter(|r| r.z == 50)
        .map(|r| r.a)
        .collect();
    tin_predicted_a.sort_unstable();
    let tin_observed_a: Vec<u16> = observed_stable_mass_numbers_for_z(50)
        .unwrap_or(&[])
        .iter()
        .copied()
        .collect();
    let tin_missing: Vec<u16> = tin_observed_a
        .iter()
        .copied()
        .filter(|a| !tin_predicted_a.contains(a))
        .collect();
    let tin_extra: Vec<u16> = tin_predicted_a
        .iter()
        .copied()
        .filter(|a| !tin_observed_a.contains(a))
        .collect();
    let mut tin_csv = String::from(
        "A,N,predicted_long_lived,observed_stable,fail_beta_optimal,fail_fissility,fail_s2n,fail_s2p,fail_sf,stability_score,s2n_mev,s2p_mev,fissility,sf_log10_half_life_s\n",
    );
    for r in records.iter().filter(|r| r.z == 50 && (100..=130).contains(&r.a)) {
        let beta_ok = beta_local_min.get(&(r.z, r.n)).copied().unwrap_or(false);
        let (pred, fail_beta, fail_fiss, fail_s2n, fail_s2p, fail_sf) = classify_long_lived(r, beta_ok);
        let observed = tin_observed_a.contains(&r.a);
        tin_csv.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
            r.a,
            r.n,
            pred,
            observed,
            fail_beta,
            fail_fiss,
            fail_s2n,
            fail_s2p,
            fail_sf,
            r.stability_score,
            r.s2n_mev.unwrap_or(f64::NAN),
            r.s2p_mev.unwrap_or(f64::NAN),
            r.fissility,
            r.sf_log10_half_life_s
        ));
    }
    fs::write(out.join("tin_isotope_diagnostics.csv"), tin_csv)?;

    let elements_with_stable_like = isotopes_per_z.len();
    let max_z_with_stable_like = isotopes_per_z.keys().max().copied().unwrap_or(0);
    let avg_isotopes_per_element = if elements_with_stable_like > 0 {
        isotopes_per_z.values().copied().sum::<usize>() as f64 / elements_with_stable_like as f64
    } else {
        0.0
    };

    let derived_closure_candidates = derived_superheavy_proton_candidates();
    let neutron_magic = magic_s2n_summary(&records);
    let proton_magic = proton_s2p_summary(&records);
    let neutron_hit_rate = if neutron_magic.is_empty() {
        0.0
    } else {
        neutron_magic
            .iter()
            .filter(|row| row.strongest_delta_s2n_mev > 1.0)
            .count() as f64
            / neutron_magic.len() as f64
    };
    let proton_hit_rate = if proton_magic.is_empty() {
        0.0
    } else {
        proton_magic
            .iter()
            .filter(|row| row.strongest_delta_s2p_mev > 1.0)
            .count() as f64
            / proton_magic.len() as f64
    };

    let mut shell_gap_csv = String::from(
        "magic_n,strongest_delta_s2n_mev,mean_delta_s2n_mev,ref_min_mev,ref_max_mev,ref_mid_mev,strongest_over_ref_mid,mean_over_ref_mid\n",
    );
    let mut heavy_gap_ratios: Vec<f64> = Vec::new();
    let mut n50_ratio = 0.0;
    let mut n82_ratio = 0.0;
    let mut n126_ratio = 0.0;
    for row in &neutron_magic {
        if let Some((ref_min, ref_max)) = reference_shell_gap_bounds_mev(row.magic_n) {
            let ref_mid = 0.5 * (ref_min + ref_max);
            let strongest_ratio = if ref_mid > 0.0 {
                row.strongest_delta_s2n_mev / ref_mid
            } else {
                0.0
            };
            let mean_ratio = if ref_mid > 0.0 {
                row.mean_delta_s2n_mev / ref_mid
            } else {
                0.0
            };
            shell_gap_csv.push_str(&format!(
                "{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
                row.magic_n,
                row.strongest_delta_s2n_mev,
                row.mean_delta_s2n_mev,
                ref_min,
                ref_max,
                ref_mid,
                strongest_ratio,
                mean_ratio
            ));
            if row.magic_n >= 50 && row.magic_n <= 126 {
                heavy_gap_ratios.push(strongest_ratio);
            }
            if row.magic_n == 50 {
                n50_ratio = strongest_ratio;
            } else if row.magic_n == 82 {
                n82_ratio = strongest_ratio;
            } else if row.magic_n == 126 {
                n126_ratio = strongest_ratio;
            }
        } else {
            shell_gap_csv.push_str(&format!(
                "{},{:.6},{:.6},,,,,\n",
                row.magic_n, row.strongest_delta_s2n_mev, row.mean_delta_s2n_mev
            ));
        }
    }
    fs::write(out.join("shell_gap_attenuation.csv"), shell_gap_csv)?;
    let heavy_gap_mean_ratio = if heavy_gap_ratios.is_empty() {
        0.0
    } else {
        heavy_gap_ratios.iter().sum::<f64>() / heavy_gap_ratios.len() as f64
    };
    let heavy_gap_min_ratio = heavy_gap_ratios
        .iter()
        .copied()
        .min_by(|a, b| a.total_cmp(b))
        .unwrap_or(0.0);

    let mut stable_presence_correct = 0usize;
    let mut stable_presence_total = 0usize;
    let mut ref_count_abs_error_sum = 0.0;
    let mut ref_count_samples = 0usize;
    let mut scoreboard_csv = String::from(
        "Z,predicted_stable_like_isotopes,predicted_has_stable,observed_has_stable,observed_stable_isotopes_ref,abs_drift_isotope_count\n",
    );
    for z in cfg.z_min..=cfg.z_max {
        let pred_count = isotopes_per_z.get(&z).copied().unwrap_or(0);
        let pred_has = pred_count > 0;
        let (obs_ref_s, obs_has, drift_s) = match observed_stable_isotope_count(z) {
            Some(obs_ref) => {
                let obs_has = obs_ref > 0;
                stable_presence_total += 1;
                if pred_has == obs_has {
                    stable_presence_correct += 1;
                }
                let drift = (pred_count as f64 - obs_ref as f64).abs();
                ref_count_abs_error_sum += drift;
                ref_count_samples += 1;
                (obs_ref.to_string(), obs_has, format!("{drift:.3}"))
            }
            None => (String::new(), false, String::new()),
        };
        scoreboard_csv.push_str(&format!(
            "{},{},{},{},{},{}\n",
            z, pred_count, pred_has, obs_has, obs_ref_s, drift_s
        ));
    }
    fs::write(out.join("periodic_table_scoreboard.csv"), scoreboard_csv)?;

    let stable_presence_accuracy = if stable_presence_total > 0 {
        stable_presence_correct as f64 / stable_presence_total as f64
    } else {
        0.0
    };
    let ref_count_mae = if ref_count_samples > 0 {
        ref_count_abs_error_sum / ref_count_samples as f64
    } else {
        0.0
    };

    let json = format!(
        concat!(
            "{{\n",
            "  \"mass_predictions\": {{\n",
            "    \"alpha_inv_struct\": {},\n",
            "    \"mp_me_struct\": {},\n",
            "    \"electron_mass_mev_pred_from_proton_anchor\": {:.9},\n",
            "    \"electron_mass_mev_obs\": {:.9},\n",
            "    \"electron_rel_error\": {:.6},\n",
            "    \"proton_mass_mev_pred_from_electron_anchor\": {:.9},\n",
            "    \"proton_mass_mev_obs\": {:.9},\n",
            "    \"proton_rel_error\": {:.6},\n",
            "    \"neutron_minus_proton_struct_mev\": {:.9},\n",
            "    \"neutron_mass_mev_pred\": {:.9},\n",
            "    \"neutron_mass_mev_obs\": {:.9},\n",
            "    \"neutron_rel_error\": {:.6}\n",
            "  }},\n",
            "  \"periodic_stats\": {{\n",
            "    \"rows\": {},\n",
            "    \"stable_like_rows\": {},\n",
            "    \"valley_rows\": {},\n",
            "    \"elements_with_stable_like\": {},\n",
            "    \"max_z_with_stable_like\": {},\n",
            "    \"avg_isotopes_per_element\": {:.3},\n",
            "    \"stable_presence_accuracy_z_le_94\": {:.6},\n",
            "    \"stable_isotope_count_mae_z_le_94\": {:.6},\n",
            "    \"top_island\": {{\"z\": {}, \"n\": {}, \"score\": {:.6}}},\n",
            "    \"closest_to_114_184\": {{\"z\": {}, \"n\": {}, \"score\": {:.6}}}\n",
            "  }},\n",
            "  \"closure_stats\": {{\n",
            "    \"neutron_magic_hit_rate\": {:.6},\n",
            "    \"proton_closure_hit_rate\": {:.6}\n",
            "  }},\n",
            "  \"shell_gap_attenuation\": {{\n",
            "    \"heavy_magic_mean_ratio\": {:.6},\n",
            "    \"heavy_magic_min_ratio\": {:.6},\n",
            "    \"n50_ratio\": {:.6},\n",
            "    \"n82_ratio\": {:.6},\n",
            "    \"n126_ratio\": {:.6}\n",
            "  }},\n",
            "  \"derived_superheavy_proton_candidates\": [{}]\n",
            "  ,\"tin_diagnostics\": {{\n",
            "    \"observed_stable_a\": [{}],\n",
            "    \"predicted_stable_like_a\": [{}],\n",
            "    \"missing_from_prediction\": [{}],\n",
            "    \"extra_in_prediction\": [{}],\n",
            "    \"observed_count\": {},\n",
            "    \"predicted_count\": {}\n",
            "  }}\n",
            "}}\n"
        ),
        alpha_inv_struct,
        mp_me_struct,
        electron_pred_from_p,
        ELECTRON_MASS_MEV_OBS,
        ((electron_pred_from_p - ELECTRON_MASS_MEV_OBS) / ELECTRON_MASS_MEV_OBS).abs(),
        proton_pred_from_e,
        PROTON_MASS_MEV_OBS,
        ((proton_pred_from_e - PROTON_MASS_MEV_OBS) / PROTON_MASS_MEV_OBS).abs(),
        neutron_minus_proton_struct,
        neutron_pred,
        NEUTRON_MASS_MEV_OBS,
        ((neutron_pred - NEUTRON_MASS_MEV_OBS) / NEUTRON_MASS_MEV_OBS).abs(),
        records.len(),
        stable_like.len(),
        valley.len(),
        elements_with_stable_like,
        max_z_with_stable_like,
        avg_isotopes_per_element,
        stable_presence_accuracy,
        ref_count_mae,
        top.map(|r| r.z).unwrap_or(0),
        top.map(|r| r.n).unwrap_or(0),
        top.map(|r| r.stability_score).unwrap_or(0.0),
        closest_target.map(|r| r.z).unwrap_or(0),
        closest_target.map(|r| r.n).unwrap_or(0),
        closest_target.map(|r| r.stability_score).unwrap_or(0.0),
        neutron_hit_rate,
        proton_hit_rate,
        heavy_gap_mean_ratio,
        heavy_gap_min_ratio,
        n50_ratio,
        n82_ratio,
        n126_ratio,
        derived_closure_candidates
            .iter()
            .map(|z| z.to_string())
            .collect::<Vec<_>>()
            .join(", "),
        tin_observed_a
            .iter()
            .map(|a| a.to_string())
            .collect::<Vec<_>>()
            .join(", "),
        tin_predicted_a
            .iter()
            .map(|a| a.to_string())
            .collect::<Vec<_>>()
            .join(", "),
        tin_missing
            .iter()
            .map(|a| a.to_string())
            .collect::<Vec<_>>()
            .join(", "),
        tin_extra
            .iter()
            .map(|a| a.to_string())
            .collect::<Vec<_>>()
            .join(", "),
        tin_observed_a.len(),
        tin_predicted_a.len()
    );
    fs::write(out.join("mass_periodic_report.json"), json)?;

    let trend_path = out.join("periodic_table_trend.csv");
    let trend_exists = trend_path.exists();
    let mut trend = OpenOptions::new().create(true).append(true).open(&trend_path)?;
    if !trend_exists {
        writeln!(
            trend,
            "timestamp_unix,rows,stable_like_rows,elements_with_stable_like,max_z_with_stable_like,stable_presence_accuracy,stable_isotope_count_mae,neutron_magic_hit_rate,proton_closure_hit_rate,top_island_z,top_island_n,top_island_score,closest_114_184_score,mp_me_struct,electron_rel_error,proton_rel_error,neutron_rel_error"
        )?;
    }
    writeln!(
        trend,
        "{},{},{},{},{},{:.6},{:.6},{:.6},{:.6},{},{},{:.6},{:.6},{},{:.6},{:.6},{:.6}",
        now_unix_seconds(),
        records.len(),
        stable_like.len(),
        elements_with_stable_like,
        max_z_with_stable_like,
        stable_presence_accuracy,
        ref_count_mae,
        neutron_hit_rate,
        proton_hit_rate,
        top.map(|r| r.z).unwrap_or(0),
        top.map(|r| r.n).unwrap_or(0),
        top.map(|r| r.stability_score).unwrap_or(0.0),
        closest_target.map(|r| r.stability_score).unwrap_or(0.0),
        mp_me_struct,
        ((electron_pred_from_p - ELECTRON_MASS_MEV_OBS) / ELECTRON_MASS_MEV_OBS).abs(),
        ((proton_pred_from_e - PROTON_MASS_MEV_OBS) / PROTON_MASS_MEV_OBS).abs(),
        ((neutron_pred - NEUTRON_MASS_MEV_OBS) / NEUTRON_MASS_MEV_OBS).abs(),
    )?;

    println!("Wrote {}", out.join("mass_periodic_report.json").display());
    println!("Wrote {}", out.join("periodic_table_scoreboard.csv").display());
    println!("Wrote {}", out.join("shell_gap_attenuation.csv").display());
    println!("Wrote {}", out.join("tin_isotope_diagnostics.csv").display());
    println!("Appended {}", trend_path.display());
    Ok(())
}
