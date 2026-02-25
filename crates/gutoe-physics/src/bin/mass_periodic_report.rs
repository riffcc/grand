use gutoe_physics::{
    closest_to_target_island, derived_superheavy_proton_candidates, rank_island_candidates_with_config,
    scan_nuclear_chart, IslandRankingConfig, ScanConfig, StandardModelDynamicsMap,
};
use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::PathBuf;

const ELECTRON_MASS_MEV_OBS: f64 = 0.510_998_950;
const PROTON_MASS_MEV_OBS: f64 = 938.272_088_16;
const NEUTRON_MASS_MEV_OBS: f64 = 939.565_420_52;

fn triangular(n: u32) -> u32 {
    n * (n + 1) / 2
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

    // Provisional no-fit structural n-p split from λ_QG × gauge-generator count.
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

    let top = ranked.first().copied();
    let closure_candidates = derived_superheavy_proton_candidates();

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
            "    \"top_island\": {{\"z\": {}, \"n\": {}, \"score\": {:.6}}},\n",
            "    \"closest_to_114_184\": {{\"z\": {}, \"n\": {}, \"score\": {:.6}}}\n",
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
        top.map(|r| r.z).unwrap_or(0),
        top.map(|r| r.n).unwrap_or(0),
        top.map(|r| r.stability_score).unwrap_or(0.0),
        closest_target.map(|r| r.z).unwrap_or(0),
        closest_target.map(|r| r.n).unwrap_or(0),
        closest_target.map(|r| r.stability_score).unwrap_or(0.0),
        closure_candidates
            .iter()
            .map(|z| z.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    );

    let path = out.join("mass_periodic_report.json");
    fs::write(&path, json)?;
    println!("Wrote {}", path.display());
    Ok(())
}
