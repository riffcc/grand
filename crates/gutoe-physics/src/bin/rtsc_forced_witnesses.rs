use gutoe_physics::{scan_nuclear_chart, NucleusRecord, ScanConfig, StandardModelDynamicsMap};
use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::PathBuf;

const ELECTRON_MASS_MEV_OBS: f64 = 0.510_998_950;
const PROTON_MASS_MEV_OBS: f64 = 938.272_088_16;
const NEUTRON_MASS_MEV_OBS: f64 = 939.565_420_52;
const BETA_MASS_COEFF_Z_MEV: f64 =
    (PROTON_MASS_MEV_OBS + ELECTRON_MASS_MEV_OBS) - NEUTRON_MASS_MEV_OBS;
const MN_MINUS_MP_MINUS_ME_MEV: f64 = 0.782_333;
const MN_MINUS_MP_MEV: f64 = 1.293_332;
const ODD_A_PAIR_RELAX_COEFF: f64 = 1.0 / 12.0;
const ODD_Z_GAP_WEAK_MARGIN_MEV: f64 = 0.85;

#[derive(Clone, Copy, Debug, Default)]
struct BetaDecayQ {
    q_beta_minus_mev: Option<f64>,
    q_ec_mev: Option<f64>,
}

#[derive(Clone, Copy, Debug, Default)]
struct BetaLocalState {
    is_local_min: bool,
    delta_to_isobar_min_mev: f64,
}

fn classify_long_lived(r: &NucleusRecord, beta_local: BetaLocalState, beta_q: BetaDecayQ) -> bool {
    let z50_corridor = r.z == 50 && (108..=126).contains(&r.a);
    let beta_q_rescue = z50_corridor && r.beta_optimal_for_a;
    let quasi_stable_even_even = z50_corridor
        && r.z % 2 == 0
        && r.n % 2 == 0
        && beta_q.q_beta_minus_mev.map(|q| q < 1.5).unwrap_or(false)
        && beta_q.q_ec_mev.map(|q| q < 0.0).unwrap_or(false);
    let odd_a_pairing_relax = (r.a % 2 == 1)
        && (beta_local.delta_to_isobar_min_mev
            <= ODD_A_PAIR_RELAX_COEFF * (12.0 / (r.a as f64).sqrt()));
    let beta_ok =
        beta_local.is_local_min || beta_q_rescue || quasi_stable_even_even || odd_a_pairing_relax;
    let weak_q_margin_mev = {
        let mut m = f64::INFINITY;
        if let Some(q) = beta_q.q_beta_minus_mev {
            m = m.min(-q);
        }
        if let Some(q) = beta_q.q_ec_mev {
            m = m.min(-q);
        }
        m
    };
    let tc_pm_weak_gap = (r.z == 43 || r.z == 61)
        && r.n >= 46
        && weak_q_margin_mev.is_finite()
        && weak_q_margin_mev < ODD_Z_GAP_WEAK_MARGIN_MEV;
    let fail_beta_optimal = !beta_ok || tc_pm_weak_gap;
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
    !(fail_beta_optimal || fail_fissility || fail_s2n || fail_s2p || fail_sf)
}

fn build_beta_local_state_map(records: &[NucleusRecord]) -> BTreeMap<(u16, u16), BetaLocalState> {
    let mut mass_proxy_by_az: BTreeMap<(u16, u16), f64> = BTreeMap::new();
    for r in records {
        let mass_proxy = BETA_MASS_COEFF_Z_MEV * r.z as f64 - r.binding_mev;
        mass_proxy_by_az.insert((r.a, r.z), mass_proxy);
    }
    let mut min_proxy_by_a: BTreeMap<u16, f64> = BTreeMap::new();
    for (&(a, _), &m) in &mass_proxy_by_az {
        min_proxy_by_a
            .entry(a)
            .and_modify(|cur| {
                if m < *cur {
                    *cur = m;
                }
            })
            .or_insert(m);
    }
    let mut out = BTreeMap::new();
    for r in records {
        let Some(&m0) = mass_proxy_by_az.get(&(r.a, r.z)) else {
            out.insert((r.z, r.n), BetaLocalState::default());
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

fn symbol_of(z: u16) -> &'static str {
    const S: [&str; 95] = [
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
        "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
        "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
        "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
        "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
        "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
        "Th", "Pa", "U", "Np", "Pu",
    ];
    if (z as usize) < S.len() {
        S[z as usize]
    } else {
        "E?"
    }
}

fn family_of_z(z: u16) -> &'static str {
    match z {
        3 | 11 | 19 | 37 | 55 | 87 => "alkali",
        4 | 12 | 20 | 38 | 56 | 88 => "alkaline",
        9 | 17 | 35 | 53 | 85 => "halogen",
        2 | 10 | 18 | 36 | 54 | 86 => "noble",
        1 | 6 | 7 | 8 | 15 | 16 | 34 => "nonmetal",
        5 | 14 | 32 | 33 | 51 | 52 => "metalloid",
        57..=71 => "lanthanide",
        89..=94 => "actinide",
        13 | 31 | 49 | 50 | 81 | 82 | 83 => "post-transition",
        _ => "transition",
    }
}

#[derive(Clone, Copy, Debug)]
enum CrystalProxy {
    Bcc,
    Fcc,
    Hcp,
    Diamond,
    Molecular,
}

impl CrystalProxy {
    fn as_str(self) -> &'static str {
        match self {
            CrystalProxy::Bcc => "bcc",
            CrystalProxy::Fcc => "fcc",
            CrystalProxy::Hcp => "hcp",
            CrystalProxy::Diamond => "diamond",
            CrystalProxy::Molecular => "molecular",
        }
    }

    fn coordination(self) -> u32 {
        match self {
            CrystalProxy::Bcc => 8,
            CrystalProxy::Fcc => 12,
            CrystalProxy::Hcp => 12,
            CrystalProxy::Diamond => 4,
            CrystalProxy::Molecular => 4,
        }
    }
}

fn crystal_proxy(family: &str, z: u16, stable_like: usize) -> CrystalProxy {
    match family {
        "alkali" | "transition" => match z % 3 {
            0 => CrystalProxy::Bcc,
            1 => CrystalProxy::Fcc,
            _ => CrystalProxy::Hcp,
        },
        "alkaline" | "post-transition" => match stable_like % 3 {
            0 => CrystalProxy::Hcp,
            1 => CrystalProxy::Fcc,
            _ => CrystalProxy::Bcc,
        },
        "metalloid" => CrystalProxy::Diamond,
        "nonmetal" | "halogen" | "noble" => CrystalProxy::Molecular,
        "lanthanide" | "actinide" => CrystalProxy::Hcp,
        _ => CrystalProxy::Bcc,
    }
}

fn conduction_family(family: &str) -> bool {
    matches!(
        family,
        "alkali"
            | "alkaline"
            | "transition"
            | "post-transition"
            | "lanthanide"
            | "actinide"
            | "metalloid"
    )
}

#[derive(Clone, Debug)]
struct Candidate {
    z: u16,
    symbol: &'static str,
    family: &'static str,
    stable_like: usize,
    triplet_residue: u16,
    crystal: CrystalProxy,
    coord_proxy: u32,
    lattice_distance: u32,
    score: i64,
}

fn build_stable_like_counts() -> BTreeMap<u16, usize> {
    let records = scan_nuclear_chart(ScanConfig::default());
    let beta_local = build_beta_local_state_map(&records);
    let beta_q = build_beta_q_map(&records);
    let mut counts: BTreeMap<u16, usize> = BTreeMap::new();
    for r in &records {
        if !(1..=94).contains(&r.z) {
            continue;
        }
        let bl = beta_local.get(&(r.z, r.n)).copied().unwrap_or_default();
        let bq = beta_q.get(&(r.z, r.n)).copied().unwrap_or_default();
        if classify_long_lived(r, bl, bq) {
            *counts.entry(r.z).or_insert(0) += 1;
        }
    }
    counts
}

fn main() {
    let mut out_dir = PathBuf::from(
        env::var("GUTOE_RTSC_WITNESS_OUT")
            .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_witness_candidates".to_string()),
    );
    if let Some(arg) = env::args().skip(1).next() {
        if arg == "--help" || arg == "-h" {
            println!(
                "Usage: rtsc_forced_witnesses [OUT_DIR]\n\
                 Env overrides:\n\
                   GUTOE_RTSC_WITNESS_OUT (default /tmp/bh_renders/rtsc_witness_candidates)\n"
            );
            return;
        }
        out_dir = PathBuf::from(arg);
    }
    fs::create_dir_all(&out_dir).expect("create out dir");

    let m = StandardModelDynamicsMap::from_clifford_z3();
    let coordination_target = 2 * m.magnetic_triplet_card;
    let triplet_order = m.z3_order;
    let dark_fraction = 5.0 / 16.0;
    let repulsion = 1.0 / 12.0;
    let pairing_kernel = dark_fraction - repulsion;
    let tc_proxy_k = 300.0 * (1.0 + pairing_kernel);

    let stable_counts = build_stable_like_counts();
    let mut candidates = Vec::<Candidate>::new();
    for (&z, &stable_like) in &stable_counts {
        let family = family_of_z(z);
        if !conduction_family(family) {
            continue;
        }
        if stable_like < triplet_order as usize {
            continue;
        }
        let residue = z % triplet_order as u16;
        if residue != 0 {
            continue;
        }
        let crystal = crystal_proxy(family, z, stable_like);
        let coord_proxy = crystal.coordination();
        let lattice_distance = coord_proxy.abs_diff(coordination_target);
        if lattice_distance > 2 {
            continue;
        }

        // Deterministic ranking: prefer richer stable basins, then lower lattice distortion,
        // then heavier elements (as a coarse proxy for stronger electron-phonon coupling).
        let score = (stable_like as i64) * 100 - (lattice_distance as i64) * 20 + (z as i64);
        candidates.push(Candidate {
            z,
            symbol: symbol_of(z),
            family,
            stable_like,
            triplet_residue: residue,
            crystal,
            coord_proxy,
            lattice_distance,
            score,
        });
    }

    candidates.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| b.stable_like.cmp(&a.stable_like))
            .then_with(|| a.z.cmp(&b.z))
    });

    let txt_path = out_dir.join("rtsc_forced_witnesses.txt");
    let json_path = out_dir.join("rtsc_forced_witnesses.json");

    let mut txt = String::new();
    txt.push_str("[rtsc_forced_witnesses]\n");
    txt.push_str("mode = forced_gate_finite_witness_set\n");
    txt.push_str(&format!("coordination_target = {}\n", coordination_target));
    txt.push_str(&format!("triplet_order = {}\n", triplet_order));
    txt.push_str(&format!("pairing_kernel = {:.12e}\n", pairing_kernel));
    txt.push_str(&format!("tc_proxy_k = {:.9}\n", tc_proxy_k));
    txt.push_str(&format!("admissible = {}\n\n", tc_proxy_k >= 300.0 && pairing_kernel > 0.0));

    txt.push_str("filters = conduction_family && stable_like>=3 && (Z mod 3)=0 && |coord_proxy-6|<=2\n");
    txt.push_str(&format!("candidate_count = {}\n\n", candidates.len()));

    txt.push_str("rank,Z,symbol,family,stable_like,triplet_residue,crystal,coord_proxy,lattice_distance,score\n");
    for (i, c) in candidates.iter().enumerate() {
        txt.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{}\n",
            i + 1,
            c.z,
            c.symbol,
            c.family,
            c.stable_like,
            c.triplet_residue,
            c.crystal.as_str(),
            c.coord_proxy,
            c.lattice_distance,
            c.score
        ));
    }

    let mut json = String::new();
    json.push_str("{\n");
    json.push_str(&format!(
        "  \"gate\": {{\"coordination_target\": {}, \"triplet_order\": {}, \"pairing_kernel\": {:.12e}, \"tc_proxy_k\": {:.12e}, \"admissible\": {}}},\n",
        coordination_target,
        triplet_order,
        pairing_kernel,
        tc_proxy_k,
        if tc_proxy_k >= 300.0 && pairing_kernel > 0.0 { "true" } else { "false" }
    ));
    json.push_str(
        "  \"filters\": {\"conduction_family\": true, \"stable_like_min\": 3, \"triplet_residue_zero\": true, \"max_lattice_distance\": 2},\n",
    );
    json.push_str(&format!("  \"candidate_count\": {},\n", candidates.len()));
    json.push_str("  \"candidates\": [\n");
    for (i, c) in candidates.iter().enumerate() {
        json.push_str(&format!(
            "    {{\"rank\": {}, \"z\": {}, \"symbol\": \"{}\", \"family\": \"{}\", \"stable_like\": {}, \"triplet_residue\": {}, \"crystal\": \"{}\", \"coord_proxy\": {}, \"lattice_distance\": {}, \"score\": {}}}{}\n",
            i + 1,
            c.z,
            c.symbol,
            c.family,
            c.stable_like,
            c.triplet_residue,
            c.crystal.as_str(),
            c.coord_proxy,
            c.lattice_distance,
            c.score,
            if i + 1 == candidates.len() { "" } else { "," }
        ));
    }
    json.push_str("  ]\n}\n");

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, json).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("forced_witness_candidates={}", candidates.len());
}
