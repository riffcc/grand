use gutoe_physics::{
    closest_to_target_island, derived_superheavy_proton_candidates, magic_s2n_summary, proton_s2p_summary,
    rank_island_candidates_with_config, scan_nuclear_chart, IslandRankingConfig, ScanConfig,
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

fn triangular(n: u32) -> u32 {
    n * (n + 1) / 2
}

fn observed_has_stable_isotope(z: u16) -> bool {
    z <= 82 && z != 43 && z != 61
}

fn observed_stable_isotope_count_ref(z: u16) -> Option<u16> {
    // High-confidence anchor set for drift metrics (Z <= 20).
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
        _ => None,
    }
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
    let stable_like: Vec<_> = records
        .iter()
        .filter(|r| r.s2n_mev.unwrap_or(-1.0) > 0.0 && r.s2p_mev.unwrap_or(-1.0) > 0.0)
        .collect();
    let valley: Vec<_> = records.iter().filter(|r| r.beta_optimal_for_a).collect();
    let closest_target = closest_to_target_island(&records, 114, 184);
    let top = ranked.first().copied();

    let mut isotopes_per_z: BTreeMap<u16, usize> = BTreeMap::new();
    for r in &stable_like {
        *isotopes_per_z.entry(r.z).or_insert(0) += 1;
    }

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
        let obs_has = observed_has_stable_isotope(z);
        if z <= 94 {
            stable_presence_total += 1;
            if pred_has == obs_has {
                stable_presence_correct += 1;
            }
        }
        let (obs_ref_s, drift_s) = match observed_stable_isotope_count_ref(z) {
            Some(obs_ref) => {
                let drift = (pred_count as f64 - obs_ref as f64).abs();
                ref_count_abs_error_sum += drift;
                ref_count_samples += 1;
                (obs_ref.to_string(), format!("{drift:.3}"))
            }
            None => (String::new(), String::new()),
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
            "    \"stable_isotope_count_mae_anchor_z_le_20\": {:.6},\n",
            "    \"top_island\": {{\"z\": {}, \"n\": {}, \"score\": {:.6}}},\n",
            "    \"closest_to_114_184\": {{\"z\": {}, \"n\": {}, \"score\": {:.6}}}\n",
            "  }},\n",
            "  \"closure_stats\": {{\n",
            "    \"neutron_magic_hit_rate\": {:.6},\n",
            "    \"proton_closure_hit_rate\": {:.6}\n",
            "  }},\n",
            "  \"derived_superheavy_proton_candidates\": [{}]\n",
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
        derived_closure_candidates
            .iter()
            .map(|z| z.to_string())
            .collect::<Vec<_>>()
            .join(", ")
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
    println!("Appended {}", trend_path.display());
    Ok(())
}
