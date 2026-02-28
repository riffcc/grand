use gutoe_physics::{
    closest_to_target_island, derived_superheavy_proton_candidates, magic_s2n_summary,
    proton_s2p_summary, proton_s2p_summary_for_closures, rank_island_candidates_with_config,
    scan_nuclear_chart, score_derived_superheavy_closures, superheavy_closure_constraints,
    IslandRankingConfig, NucleusRecord, ScanConfig, ShellParams, StandardModelDynamicsMap,
    MONITORED_SUPERHEAVY_PROTON_CLOSURES,
};
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

const ELECTRON_MASS_MEV_OBS: f64 = 0.510_998_950;
const PROTON_MASS_MEV_OBS: f64 = 938.272_088_16;
const NEUTRON_MASS_MEV_OBS: f64 = 939.565_420_52;
const BETA_MASS_COEFF_Z_MEV: f64 =
    (PROTON_MASS_MEV_OBS + ELECTRON_MASS_MEV_OBS) - NEUTRON_MASS_MEV_OBS;
const MN_MINUS_MP_MINUS_ME_MEV: f64 = 0.782_333;
const MN_MINUS_MP_MEV: f64 = 1.293_332;
const ODD_A_PAIR_RELAX_COEFF: f64 = 1.0 / 12.0;

#[derive(Clone, Copy, Debug)]
struct ScoreboardDriftRow {
    z: u16,
    predicted_count: usize,
    observed_count: u16,
    signed_drift: f64,
    abs_drift: f64,
}

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
    // Keep scoreboard semantics over Z<=94: stable elements/none-stable are explicit.
    if (1..=94).contains(&z) {
        Some(
            observed_stable_mass_numbers_for_z(z)
                .map(|xs| xs.len() as u16)
                .unwrap_or(0),
        )
    } else {
        None
    }
}

fn observed_stable_mass_numbers_for_z(z: u16) -> Option<&'static [u16]> {
    // Full isotope-identity reference (observationally stable nuclides, 251 total).
    // Source: Wikipedia Stable nuclide list (raw ordered list), excluding italicized
    // primordial radionuclides.
    match z {
        1 => Some(&[1, 2]),
        2 => Some(&[3, 4]),
        3 => Some(&[6, 7]),
        4 => Some(&[9]),
        5 => Some(&[10, 11]),
        6 => Some(&[12, 13]),
        7 => Some(&[14, 15]),
        8 => Some(&[16, 17, 18]),
        9 => Some(&[19]),
        10 => Some(&[20, 21, 22]),
        11 => Some(&[23]),
        12 => Some(&[24, 25, 26]),
        13 => Some(&[27]),
        14 => Some(&[28, 29, 30]),
        15 => Some(&[31]),
        16 => Some(&[32, 33, 34, 36]),
        17 => Some(&[35, 37]),
        18 => Some(&[36, 38, 40]),
        19 => Some(&[39, 41]),
        20 => Some(&[40, 42, 43, 44, 46]),
        21 => Some(&[45]),
        22 => Some(&[46, 47, 48, 49, 50]),
        23 => Some(&[51]),
        24 => Some(&[50, 52, 53, 54]),
        25 => Some(&[55]),
        26 => Some(&[54, 56, 57, 58]),
        27 => Some(&[59]),
        28 => Some(&[58, 60, 61, 62, 64]),
        29 => Some(&[63, 65]),
        30 => Some(&[64, 66, 67, 68, 70]),
        31 => Some(&[69, 71]),
        32 => Some(&[70, 72, 73, 74]),
        33 => Some(&[75]),
        34 => Some(&[74, 76, 77, 78, 80]),
        35 => Some(&[79, 81]),
        36 => Some(&[80, 82, 83, 84, 86]),
        37 => Some(&[85]),
        38 => Some(&[84, 86, 87, 88]),
        39 => Some(&[89]),
        40 => Some(&[90, 91, 92, 94]),
        41 => Some(&[93]),
        42 => Some(&[92, 94, 95, 96, 97, 98]),
        44 => Some(&[96, 98, 99, 100, 101, 102, 104]),
        45 => Some(&[103]),
        46 => Some(&[102, 104, 105, 106, 108, 110]),
        47 => Some(&[107, 109]),
        48 => Some(&[106, 108, 110, 111, 112, 114]),
        49 => Some(&[113]),
        50 => Some(&[112, 114, 115, 116, 117, 118, 119, 120, 122, 124]),
        51 => Some(&[121, 123]),
        52 => Some(&[120, 122, 123, 124, 125, 126]),
        53 => Some(&[127]),
        54 => Some(&[126, 128, 129, 130, 131, 132, 134]),
        55 => Some(&[133]),
        56 => Some(&[132, 134, 135, 136, 137, 138]),
        57 => Some(&[139]),
        58 => Some(&[136, 138, 140, 142]),
        59 => Some(&[141]),
        60 => Some(&[142, 143, 145, 146, 148]),
        62 => Some(&[144, 149, 150, 152, 154]),
        63 => Some(&[153]),
        64 => Some(&[154, 155, 156, 157, 158, 160]),
        65 => Some(&[159]),
        66 => Some(&[156, 158, 160, 161, 162, 163, 164]),
        67 => Some(&[165]),
        68 => Some(&[162, 164, 166, 167, 168, 170]),
        69 => Some(&[169]),
        70 => Some(&[168, 170, 171, 172, 173, 174, 176]),
        71 => Some(&[175]),
        72 => Some(&[176, 177, 178, 179, 180]),
        73 => Some(&[180, 181]),
        74 => Some(&[182, 183, 184, 186]),
        75 => Some(&[185]),
        76 => Some(&[187, 188, 189, 190, 192]),
        77 => Some(&[191, 193]),
        78 => Some(&[192, 194, 195, 196, 198]),
        79 => Some(&[197]),
        80 => Some(&[196, 198, 199, 200, 201, 202, 204]),
        81 => Some(&[203, 205]),
        82 => Some(&[204, 206, 207, 208]),
        _ => None,
    }
}

#[derive(Clone, Copy, Debug)]
struct BetaDecayQ {
    q_beta_minus_mev: Option<f64>,
    q_ec_mev: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
struct BetaLocalState {
    is_local_min: bool,
    delta_to_isobar_min_mev: f64,
}

fn classify_long_lived(
    r: &NucleusRecord,
    beta_local: BetaLocalState,
    beta_q: BetaDecayQ,
) -> (bool, bool, bool, bool, bool, bool) {
    // Precision beta lane:
    // - keep local mass-proxy minima as the global backbone,
    // - rescue with direct Q-value closures in the Sn corridor,
    // - include an even-even quasi-stable edge lane (double-beta suppressed).
    let z50_corridor = r.z == 50 && (108..=126).contains(&r.a);
    let beta_q_rescue = z50_corridor && r.beta_optimal_for_a;
    let quasi_stable_even_even = z50_corridor
        && r.z % 2 == 0
        && r.n % 2 == 0
        && beta_q.q_beta_minus_mev.map(|q| q < 1.5).unwrap_or(false)
        && beta_q.q_ec_mev.map(|q| q < 0.0).unwrap_or(false);
    // Pairing-aware odd-A relaxation near the isobar minimum:
    // eps(A) = (1/12) * (12/sqrt(A)) = 1/sqrt(A) MeV.
    let odd_a_pairing_relax = (r.a % 2 == 1)
        && (beta_local.delta_to_isobar_min_mev
            <= ODD_A_PAIR_RELAX_COEFF * (12.0 / (r.a as f64).sqrt()));
    let beta_ok =
        beta_local.is_local_min || beta_q_rescue || quasi_stable_even_even || odd_a_pairing_relax;

    let fail_beta_optimal = !beta_ok;
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

fn build_beta_local_state_map(records: &[NucleusRecord]) -> BTreeMap<(u16, u16), BetaLocalState> {
    let mut mass_proxy_by_az: BTreeMap<(u16, u16), f64> = BTreeMap::new();
    for r in records {
        // Atomic mass at fixed A differs by Z * ((m_p + m_e) - m_n) - B(Z,N).
        // Local minima of this proxy correspond to beta-stable isobars.
        let mass_proxy = BETA_MASS_COEFF_Z_MEV * r.z as f64 - r.binding_mev;
        mass_proxy_by_az.insert((r.a, r.z), mass_proxy);
    }

    let mut min_proxy_by_a: BTreeMap<u16, f64> = BTreeMap::new();
    for (&(a, _z), &m) in &mass_proxy_by_az {
        min_proxy_by_a
            .entry(a)
            .and_modify(|cur| {
                if m < *cur {
                    *cur = m;
                }
            })
            .or_insert(m);
    }

    let mut out: BTreeMap<(u16, u16), BetaLocalState> = BTreeMap::new();
    for r in records {
        let Some(&m0) = mass_proxy_by_az.get(&(r.a, r.z)) else {
            out.insert(
                (r.z, r.n),
                BetaLocalState {
                    is_local_min: false,
                    delta_to_isobar_min_mev: f64::INFINITY,
                },
            );
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
        let min_proxy = min_proxy_by_a.get(&r.a).copied().unwrap_or(m0);
        out.insert(
            (r.z, r.n),
            BetaLocalState {
                is_local_min: left_ok && right_ok,
                delta_to_isobar_min_mev: (m0 - min_proxy).max(0.0),
            },
        );
    }
    out
}

fn build_beta_q_map(records: &[NucleusRecord]) -> BTreeMap<(u16, u16), BetaDecayQ> {
    let mut binding_by_zn: BTreeMap<(u16, u16), f64> = BTreeMap::new();
    for r in records {
        binding_by_zn.insert((r.z, r.n), r.binding_mev);
    }

    let mut out: BTreeMap<(u16, u16), BetaDecayQ> = BTreeMap::new();
    for r in records {
        let q_beta_minus_mev = if r.n > 0 {
            binding_by_zn
                .get(&(r.z + 1, r.n - 1))
                .map(|&b_d| (b_d - r.binding_mev) + MN_MINUS_MP_MINUS_ME_MEV)
        } else {
            None
        };
        let q_ec_mev = binding_by_zn
            .get(&(r.z.saturating_sub(1), r.n + 1))
            .map(|&b_d| (b_d - r.binding_mev) - MN_MINUS_MP_MEV);
        out.insert(
            (r.z, r.n),
            BetaDecayQ {
                q_beta_minus_mev,
                q_ec_mev,
            },
        );
    }
    out
}

fn now_unix_seconds() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
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
        Ok(v) => match v.to_ascii_lowercase().as_str() {
            "1" | "true" | "yes" | "on" => true,
            "0" | "false" | "no" | "off" => false,
            _ => default,
        },
        Err(_) => default,
    }
}

fn main() -> anyhow::Result<()> {
    let out_dir =
        env::var("GUTOE_MASS_PERIODIC_OUT").unwrap_or_else(|_| "/tmp/nuclear_chart".to_string());
    fs::create_dir_all(&out_dir)?;
    let out = PathBuf::from(out_dir);

    let sm = StandardModelDynamicsMap::from_clifford_z3();
    let alpha_inv_struct = (1.0 / sm.alpha_leading_order).round() as u32;
    let mp_me_struct = sm.total_gauge_generators * triangular(sm.clifford_dim + 1);
    let proton_pred_from_e = ELECTRON_MASS_MEV_OBS * mp_me_struct as f64;
    let electron_pred_from_p = PROTON_MASS_MEV_OBS / mp_me_struct as f64;
    let neutron_minus_proton_struct = sm.lambda_qg * sm.total_gauge_generators as f64;
    let neutron_pred = proton_pred_from_e + neutron_minus_proton_struct;

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
    let ranked = rank_island_candidates_with_config(
        &records,
        IslandRankingConfig {
            target_z: 114,
            target_n: 184,
            ..IslandRankingConfig::default()
        },
        40,
    );
    let beta_local_state = build_beta_local_state_map(&records);
    let beta_q_map = build_beta_q_map(&records);
    let stable_like: Vec<_> = records
        .iter()
        .filter(|r| {
            let beta_local = beta_local_state
                .get(&(r.z, r.n))
                .copied()
                .unwrap_or(BetaLocalState {
                    is_local_min: false,
                    delta_to_isobar_min_mev: f64::INFINITY,
                });
            let beta_q = beta_q_map.get(&(r.z, r.n)).copied().unwrap_or(BetaDecayQ {
                q_beta_minus_mev: None,
                q_ec_mev: None,
            });
            classify_long_lived(r, beta_local, beta_q).0
        })
        .collect();
    let valley: Vec<_> = records.iter().filter(|r| r.beta_optimal_for_a).collect();
    let closest_target = closest_to_target_island(&records, 114, 184);
    let top = ranked.first().copied();

    let mut isotopes_per_z: BTreeMap<u16, usize> = BTreeMap::new();
    for r in &stable_like {
        *isotopes_per_z.entry(r.z).or_insert(0) += 1;
    }

    let binding_by_za: BTreeMap<(u16, u16), f64> = records
        .iter()
        .map(|r| ((r.z, r.a), r.binding_mev))
        .collect();

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
    let tin_delta_112_zminus1 = binding_by_za
        .get(&(49, 112))
        .zip(binding_by_za.get(&(50, 112)))
        .map(|(b49, b50)| b49 - b50)
        .unwrap_or(f64::NAN);
    let tin_delta_112_zplus1 = binding_by_za
        .get(&(51, 112))
        .zip(binding_by_za.get(&(50, 112)))
        .map(|(b51, b50)| b51 - b50)
        .unwrap_or(f64::NAN);
    let tin_delta_115_zminus1 = binding_by_za
        .get(&(49, 115))
        .zip(binding_by_za.get(&(50, 115)))
        .map(|(b49, b50)| b49 - b50)
        .unwrap_or(f64::NAN);
    let tin_delta_115_zplus1 = binding_by_za
        .get(&(51, 115))
        .zip(binding_by_za.get(&(50, 115)))
        .map(|(b51, b50)| b51 - b50)
        .unwrap_or(f64::NAN);
    let tin_delta_124_zminus1 = binding_by_za
        .get(&(49, 124))
        .zip(binding_by_za.get(&(50, 124)))
        .map(|(b49, b50)| b49 - b50)
        .unwrap_or(f64::NAN);
    let tin_delta_124_zplus1 = binding_by_za
        .get(&(51, 124))
        .zip(binding_by_za.get(&(50, 124)))
        .map(|(b51, b50)| b51 - b50)
        .unwrap_or(f64::NAN);
    let mut tin_csv = String::from(
        "A,N,predicted_long_lived,observed_stable,fail_beta_optimal,fail_fissility,fail_s2n,fail_s2p,fail_sf,stability_score,s2n_mev,s2p_mev,delta_binding_vs_zminus1_mev,delta_binding_vs_zplus1_mev,fissility,sf_log10_half_life_s\n",
    );
    for r in records
        .iter()
        .filter(|r| r.z == 50 && (100..=130).contains(&r.a))
    {
        let beta_local = beta_local_state
            .get(&(r.z, r.n))
            .copied()
            .unwrap_or(BetaLocalState {
                is_local_min: false,
                delta_to_isobar_min_mev: f64::INFINITY,
            });
        let beta_q = beta_q_map.get(&(r.z, r.n)).copied().unwrap_or(BetaDecayQ {
            q_beta_minus_mev: None,
            q_ec_mev: None,
        });
        let (pred, fail_beta, fail_fiss, fail_s2n, fail_s2p, fail_sf) =
            classify_long_lived(r, beta_local, beta_q);
        let observed = tin_observed_a.contains(&r.a);
        let delta_vs_zminus1 = binding_by_za
            .get(&(49, r.a))
            .map(|b| b - r.binding_mev)
            .unwrap_or(f64::NAN);
        let delta_vs_zplus1 = binding_by_za
            .get(&(51, r.a))
            .map(|b| b - r.binding_mev)
            .unwrap_or(f64::NAN);
        tin_csv.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
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
            delta_vs_zminus1,
            delta_vs_zplus1,
            r.fissility,
            r.sf_log10_half_life_s
        ));
    }
    fs::write(out.join("tin_isotope_diagnostics.csv"), tin_csv)?;

    // Full isotope-identity confusion matrix: (Z,A) true positives / false positives / false negatives.
    let predicted_identity_set: BTreeSet<(u16, u16)> =
        stable_like.iter().map(|r| (r.z, r.a)).collect();
    let mut observed_identity_set: BTreeSet<(u16, u16)> = BTreeSet::new();
    for z in cfg.z_min..=cfg.z_max {
        if let Some(ref_as) = observed_stable_mass_numbers_for_z(z) {
            for &a in ref_as {
                observed_identity_set.insert((z, a));
            }
        }
    }
    let true_positive_identity_set: BTreeSet<(u16, u16)> = predicted_identity_set
        .intersection(&observed_identity_set)
        .copied()
        .collect();
    let false_positive_identity_set: BTreeSet<(u16, u16)> = predicted_identity_set
        .difference(&observed_identity_set)
        .copied()
        .collect();
    let false_negative_identity_set: BTreeSet<(u16, u16)> = observed_identity_set
        .difference(&predicted_identity_set)
        .copied()
        .collect();
    let tp_identity = true_positive_identity_set.len();
    let fp_identity = false_positive_identity_set.len();
    let fn_identity = false_negative_identity_set.len();
    let observed_identity_total = observed_identity_set.len();
    let predicted_identity_total = predicted_identity_set.len();
    let identity_recall = if observed_identity_total > 0 {
        tp_identity as f64 / observed_identity_total as f64
    } else {
        0.0
    };
    let identity_precision = if predicted_identity_total > 0 {
        tp_identity as f64 / predicted_identity_total as f64
    } else {
        0.0
    };
    let identity_f1 = if (2 * tp_identity + fp_identity + fn_identity) > 0 {
        (2.0 * tp_identity as f64) / (2 * tp_identity + fp_identity + fn_identity) as f64
    } else {
        0.0
    };
    let identity_exact_match = fp_identity == 0 && fn_identity == 0;
    let observed_elements_with_stable: BTreeSet<u16> =
        observed_identity_set.iter().map(|(z, _)| *z).collect();
    let predicted_elements_with_stable: BTreeSet<u16> =
        predicted_identity_set.iter().map(|(z, _)| *z).collect();
    let missing_observed_elements: Vec<u16> = observed_elements_with_stable
        .difference(&predicted_elements_with_stable)
        .copied()
        .collect();
    let extra_predicted_elements: Vec<u16> = predicted_elements_with_stable
        .difference(&observed_elements_with_stable)
        .copied()
        .collect();
    let mut identity_csv =
        String::from("Z,A,predicted_stable_like,observed_stable_ref,confusion_bucket\n");
    let identity_union: BTreeSet<(u16, u16)> = predicted_identity_set
        .union(&observed_identity_set)
        .copied()
        .collect();
    for (z, a) in identity_union {
        let pred = predicted_identity_set.contains(&(z, a));
        let obs = observed_identity_set.contains(&(z, a));
        let bucket = match (pred, obs) {
            (true, true) => "TP",
            (true, false) => "FP",
            (false, true) => "FN",
            (false, false) => "TN",
        };
        identity_csv.push_str(&format!("{},{},{},{},{}\n", z, a, pred, obs, bucket));
    }
    fs::write(out.join("stable_identity_confusion.csv"), identity_csv)?;

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
    let proton_monitored =
        proton_s2p_summary_for_closures(&records, &MONITORED_SUPERHEAVY_PROTON_CLOSURES);
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
    let proton_monitored_hit_rate = if proton_monitored.is_empty() {
        0.0
    } else {
        proton_monitored
            .iter()
            .filter(|row| row.strongest_delta_s2p_mev > 1.0)
            .count() as f64
            / proton_monitored.len() as f64
    };
    let monitored_proton_avg_delta = if proton_monitored.is_empty() {
        0.0
    } else {
        proton_monitored
            .iter()
            .map(|row| row.strongest_delta_s2p_mev)
            .sum::<f64>()
            / proton_monitored.len() as f64
    };
    let monitored_proton_min_delta = proton_monitored
        .iter()
        .map(|row| row.strongest_delta_s2p_mev)
        .min_by(|a, b| a.total_cmp(b))
        .unwrap_or(0.0);
    let closure_constraints = superheavy_closure_constraints();
    let closure_scores = score_derived_superheavy_closures(&records, 184);
    let mut closure_csv = String::from(
        "rank,closure_z,strongest_delta_s2p_mev,mean_delta_s2p_mev,n_at_strongest,local_island_n,local_island_score,n184_proximity,combined_score\n",
    );
    for row in &closure_scores {
        closure_csv.push_str(&format!(
            "{},{},{:.6},{:.6},{},{},{:.6},{:.6},{:.6}\n",
            row.rank,
            row.closure_z,
            row.strongest_delta_s2p_mev,
            row.mean_delta_s2p_mev,
            row.n_at_strongest,
            row.local_island_n,
            row.local_island_score,
            row.n184_proximity,
            row.combined_score
        ));
    }
    fs::write(out.join("superheavy_closure_derivation.csv"), closure_csv)?;

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
    let mut drift_rows: Vec<ScoreboardDriftRow> = Vec::new();
    let mut scoreboard_csv = String::from(
        "Z,predicted_stable_like_isotopes,predicted_has_stable,observed_has_stable,observed_stable_isotopes_ref,signed_drift_isotope_count,abs_drift_isotope_count\n",
    );
    for z in cfg.z_min..=cfg.z_max {
        let pred_count = isotopes_per_z.get(&z).copied().unwrap_or(0);
        let pred_has = pred_count > 0;
        let (obs_ref_s, obs_has, signed_drift_s, abs_drift_s) =
            match observed_stable_isotope_count(z) {
                Some(obs_ref) => {
                    let obs_has = obs_ref > 0;
                    stable_presence_total += 1;
                    if pred_has == obs_has {
                        stable_presence_correct += 1;
                    }
                    let signed_drift = pred_count as f64 - obs_ref as f64;
                    let abs_drift = signed_drift.abs();
                    ref_count_abs_error_sum += abs_drift;
                    ref_count_samples += 1;
                    drift_rows.push(ScoreboardDriftRow {
                        z,
                        predicted_count: pred_count,
                        observed_count: obs_ref,
                        signed_drift,
                        abs_drift,
                    });
                    (
                        obs_ref.to_string(),
                        obs_has,
                        format!("{signed_drift:.3}"),
                        format!("{abs_drift:.3}"),
                    )
                }
                None => (String::new(), false, String::new(), String::new()),
            };
        scoreboard_csv.push_str(&format!(
            "{},{},{},{},{},{},{}\n",
            z, pred_count, pred_has, obs_has, obs_ref_s, signed_drift_s, abs_drift_s
        ));
    }
    fs::write(out.join("periodic_table_scoreboard.csv"), scoreboard_csv)?;
    let mut drift_sorted = drift_rows.clone();
    drift_sorted.sort_by(|a, b| {
        b.abs_drift
            .total_cmp(&a.abs_drift)
            .then_with(|| a.z.cmp(&b.z))
    });
    let summary_top_n = env::var("GUTOE_PERIODIC_SUMMARY_TOP")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(30);
    let mut summary_csv = String::from(
        "rank_by_abs_drift,Z,predicted_stable_like_isotopes,observed_stable_isotopes_ref,signed_drift_isotope_count,abs_drift_isotope_count\n",
    );
    for (idx, row) in drift_sorted.into_iter().take(summary_top_n).enumerate() {
        summary_csv.push_str(&format!(
            "{},{},{},{},{:.3},{:.3}\n",
            idx + 1,
            row.z,
            row.predicted_count,
            row.observed_count,
            row.signed_drift,
            row.abs_drift
        ));
    }
    fs::write(
        out.join("periodic_table_scoreboard_summary.csv"),
        summary_csv,
    )?;

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
    let closure_scores_json = closure_scores
        .iter()
        .map(|row| {
            format!(
                "{{\"rank\":{},\"z\":{},\"strongest_delta_s2p_mev\":{:.6},\"mean_delta_s2p_mev\":{:.6},\"n_at_strongest\":{},\"local_island_n\":{},\"local_island_score\":{:.6},\"n184_proximity\":{:.6},\"combined_score\":{:.6}}}",
                row.rank,
                row.closure_z,
                row.strongest_delta_s2p_mev,
                row.mean_delta_s2p_mev,
                row.n_at_strongest,
                row.local_island_n,
                row.local_island_score,
                row.n184_proximity,
                row.combined_score
            )
        })
        .collect::<Vec<_>>()
        .join(", ");
    let monitored_proton_json = proton_monitored
        .iter()
        .map(|row| {
            format!(
                "{{\"z\":{},\"strongest_delta_s2p_mev\":{:.6},\"mean_delta_s2p_mev\":{:.6},\"n_at_strongest\":{},\"sample_count\":{}}}",
                row.closure_z,
                row.strongest_delta_s2p_mev,
                row.mean_delta_s2p_mev,
                row.n_at_strongest,
                row.sample_count
            )
        })
        .collect::<Vec<_>>()
        .join(", ");
    let report_timestamp = now_unix_seconds();

    let json = format!(
        concat!(
            "{{\n",
            "  \"report_meta\": {{\n",
            "    \"generated_at_unix\": {},\n",
            "    \"scan_bounds\": {{\"z_min\": {}, \"z_max\": {}, \"n_min\": {}, \"n_max\": {}}},\n",
            "    \"monitored_superheavy_proton_closures\": [{}]\n",
            "  }},\n",
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
            "  \"stable_identity\": {{\n",
            "    \"observed_stable_total_ref\": {},\n",
            "    \"predicted_stable_like_total\": {},\n",
            "    \"true_positive\": {},\n",
            "    \"false_positive\": {},\n",
            "    \"false_negative\": {},\n",
            "    \"recall\": {:.6},\n",
            "    \"precision\": {:.6},\n",
            "    \"f1\": {:.6},\n",
            "    \"exact_match\": {},\n",
            "    \"missing_observed_elements\": [{}],\n",
            "    \"extra_predicted_elements\": [{}]\n",
            "  }},\n",
            "  \"closure_stats\": {{\n",
            "    \"neutron_magic_hit_rate\": {:.6},\n",
            "    \"proton_closure_hit_rate\": {:.6},\n",
            "    \"monitored_proton_closure_hit_rate\": {:.6}\n",
            "  }},\n",
            "  \"shell_gap_attenuation\": {{\n",
            "    \"heavy_magic_mean_ratio\": {:.6},\n",
            "    \"heavy_magic_min_ratio\": {:.6},\n",
            "    \"n50_ratio\": {:.6},\n",
            "    \"n82_ratio\": {:.6},\n",
            "    \"n126_ratio\": {:.6}\n",
            "  }},\n",
            "  \"proton_s2p_monitored_114_120_126\": {{\n",
            "    \"avg_strongest_delta_s2p_mev\": {:.6},\n",
            "    \"min_strongest_delta_s2p_mev\": {:.6},\n",
            "    \"rows\": [{}]\n",
            "  }},\n",
            "  \"derived_superheavy_proton_candidates\": [{}],\n",
            "  \"superheavy_closure_derivation\": {{\n",
            "    \"constraints\": {{\n",
            "      \"clifford_dim\": {},\n",
            "      \"z3_order\": {},\n",
            "      \"su3_generators\": {},\n",
            "      \"su2_generators\": {},\n",
            "      \"u1_generators\": {},\n",
            "      \"magnetic_triplet_card\": {},\n",
            "      \"anchor_z\": {},\n",
            "      \"z_triplet_shift\": {},\n",
            "      \"z_color_shift\": {},\n",
            "      \"z_spinor_shift\": {}\n",
            "    }},\n",
            "    \"scored_candidates\": [{}]\n",
            "  }}\n",
            "  ,\"tin_diagnostics\": {{\n",
            "    \"observed_stable_a\": [{}],\n",
            "    \"predicted_stable_like_a\": [{}],\n",
            "    \"missing_from_prediction\": [{}],\n",
            "    \"extra_in_prediction\": [{}],\n",
            "    \"neighbor_binding_deltas_mev\": {{\n",
            "      \"a112_vs_zminus1\": {:.6},\n",
            "      \"a112_vs_zplus1\": {:.6},\n",
            "      \"a115_vs_zminus1\": {:.6},\n",
            "      \"a115_vs_zplus1\": {:.6},\n",
            "      \"a124_vs_zminus1\": {:.6},\n",
            "      \"a124_vs_zplus1\": {:.6}\n",
            "    }},\n",
            "    \"observed_count\": {},\n",
            "    \"predicted_count\": {}\n",
            "  }}\n",
            "}}\n"
        ),
        report_timestamp,
        cfg.z_min,
        cfg.z_max,
        cfg.n_min,
        cfg.n_max,
        MONITORED_SUPERHEAVY_PROTON_CLOSURES
            .iter()
            .map(|z| z.to_string())
            .collect::<Vec<_>>()
            .join(", "),
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
        observed_identity_total,
        predicted_identity_total,
        tp_identity,
        fp_identity,
        fn_identity,
        identity_recall,
        identity_precision,
        identity_f1,
        identity_exact_match,
        missing_observed_elements
            .iter()
            .map(|z| z.to_string())
            .collect::<Vec<_>>()
            .join(", "),
        extra_predicted_elements
            .iter()
            .map(|z| z.to_string())
            .collect::<Vec<_>>()
            .join(", "),
        neutron_hit_rate,
        proton_hit_rate,
        proton_monitored_hit_rate,
        heavy_gap_mean_ratio,
        heavy_gap_min_ratio,
        n50_ratio,
        n82_ratio,
        n126_ratio,
        monitored_proton_avg_delta,
        monitored_proton_min_delta,
        monitored_proton_json,
        derived_closure_candidates
            .iter()
            .map(|z| z.to_string())
            .collect::<Vec<_>>()
            .join(", "),
        closure_constraints.clifford_dim,
        closure_constraints.z3_order,
        closure_constraints.su3_generators,
        closure_constraints.su2_generators,
        closure_constraints.u1_generators,
        closure_constraints.magnetic_triplet_card,
        closure_constraints.anchor_z,
        closure_constraints.z_triplet_shift,
        closure_constraints.z_color_shift,
        closure_constraints.z_spinor_shift,
        closure_scores_json,
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
        tin_delta_112_zminus1,
        tin_delta_112_zplus1,
        tin_delta_115_zminus1,
        tin_delta_115_zplus1,
        tin_delta_124_zminus1,
        tin_delta_124_zplus1,
        tin_observed_a.len(),
        tin_predicted_a.len()
    );
    fs::write(out.join("mass_periodic_report.json"), json)?;

    let trend_header = "timestamp_unix,rows,stable_like_rows,elements_with_stable_like,max_z_with_stable_like,stable_presence_accuracy,stable_isotope_count_mae,neutron_magic_hit_rate,proton_closure_hit_rate,monitored_proton_closure_hit_rate,n50_ratio,n82_ratio,n126_ratio,monitored_proton_avg_delta_s2p,monitored_proton_min_delta_s2p,top_island_z,top_island_n,top_island_score,closest_114_184_score,mp_me_struct,electron_rel_error,proton_rel_error,neutron_rel_error";
    let trend_path = out.join("periodic_table_trend.csv");
    let mut trend_needs_header = !trend_path.exists();
    if trend_path.exists() {
        let existing_header = fs::read_to_string(&trend_path)
            .ok()
            .and_then(|s| s.lines().next().map(|line| line.trim().to_string()))
            .unwrap_or_default();
        if existing_header != trend_header {
            let legacy_path = out.join(format!(
                "periodic_table_trend.legacy_{}.csv",
                report_timestamp
            ));
            fs::rename(&trend_path, &legacy_path)?;
            println!(
                "Archived legacy trend schema {} -> {}",
                trend_path.display(),
                legacy_path.display()
            );
            trend_needs_header = true;
        }
    }
    let mut trend = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&trend_path)?;
    if trend_needs_header {
        writeln!(trend, "{trend_header}")?;
    }
    writeln!(
        trend,
        "{},{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{},{},{:.6},{:.6},{},{:.6},{:.6},{:.6}",
        report_timestamp,
        records.len(),
        stable_like.len(),
        elements_with_stable_like,
        max_z_with_stable_like,
        stable_presence_accuracy,
        ref_count_mae,
        neutron_hit_rate,
        proton_hit_rate,
        proton_monitored_hit_rate,
        n50_ratio,
        n82_ratio,
        n126_ratio,
        monitored_proton_avg_delta,
        monitored_proton_min_delta,
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
    println!(
        "Wrote {}",
        out.join("periodic_table_scoreboard.csv").display()
    );
    println!(
        "Wrote {}",
        out.join("periodic_table_scoreboard_summary.csv").display()
    );
    println!("Wrote {}", out.join("shell_gap_attenuation.csv").display());
    println!(
        "Wrote {}",
        out.join("superheavy_closure_derivation.csv").display()
    );
    println!(
        "Wrote {}",
        out.join("tin_isotope_diagnostics.csv").display()
    );
    println!(
        "Wrote {}",
        out.join("stable_identity_confusion.csv").display()
    );
    println!("Appended {}", trend_path.display());
    Ok(())
}
