use gutoe_physics::{gamow_factor, ALPHA_LEADING_ORDER};
use std::env;
use std::fs;

const NUCLEON_MASS_MEV: f64 = 938.272_088_16;

#[derive(Clone, Copy, Debug)]
struct ReactionModel {
    id: &'static str,
    z1: u16,
    a1: u16,
    z2: u16,
    a2: u16,
    q_mev: f64,
    branch_weight: f64,
    neutron_fraction: f64,
}

#[derive(Clone, Copy, Debug)]
enum PairMode {
    MatterMatter,
    MatterAntimatter,
    AntimatterAntimatter,
}

impl PairMode {
    fn id(self) -> &'static str {
        match self {
            PairMode::MatterMatter => "matter+matter",
            PairMode::MatterAntimatter => "matter+antimatter",
            PairMode::AntimatterAntimatter => "antimatter+antimatter",
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct ReactionScore {
    t_kev: f64,
    gamow: f64,
    fusion_raw: f64,
    fusion_clean: f64,
    annihilation_raw: f64,
    total_raw: f64,
    total_clean: f64,
}

fn reduced_mass_mev(a1: u16, a2: u16) -> f64 {
    let a1f = a1 as f64;
    let a2f = a2 as f64;
    (a1f * a2f / (a1f + a2f)) * 931.494
}

fn gamow_with_charges(z1: u16, z2: u16, a1: u16, a2: u16, e_cm_kev: f64) -> f64 {
    let e_cm_mev = e_cm_kev / 1000.0;
    if e_cm_mev <= 0.0 {
        return 0.0;
    }
    let alpha_eff = ALPHA_LEADING_ORDER * (z1 as f64) * (z2 as f64);
    let m_reduced = reduced_mass_mev(a1, a2);
    gamow_factor(alpha_eff, m_reduced, e_cm_mev).unwrap_or(0.0)
}

fn dd_gamma_branch_suppression(q_gamma_mev: f64) -> f64 {
    (ALPHA_LEADING_ORDER * (q_gamma_mev / NUCLEON_MASS_MEV).powi(3)).clamp(0.0, 1.0)
}

fn model_set() -> Vec<ReactionModel> {
    let q_dd_gamma = 23.846_530;
    let dd_gamma = dd_gamma_branch_suppression(q_dd_gamma);
    vec![
        ReactionModel {
            id: "D+T->n+He4",
            z1: 1,
            a1: 2,
            z2: 1,
            a2: 3,
            q_mev: 17.589_300,
            branch_weight: 1.0,
            neutron_fraction: 0.80,
        },
        ReactionModel {
            id: "D+He3->p+He4",
            z1: 1,
            a1: 2,
            z2: 2,
            a2: 3,
            q_mev: 18.353_055,
            branch_weight: 1.0,
            neutron_fraction: 0.0,
        },
        ReactionModel {
            id: "D+D->p+T",
            z1: 1,
            a1: 2,
            z2: 1,
            a2: 2,
            q_mev: 4.033,
            branch_weight: 0.5,
            neutron_fraction: 0.45,
        },
        ReactionModel {
            id: "D+D->n+He3",
            z1: 1,
            a1: 2,
            z2: 1,
            a2: 2,
            q_mev: 3.269,
            branch_weight: 0.5,
            neutron_fraction: 0.90,
        },
        ReactionModel {
            id: "D+D->gamma+He4",
            z1: 1,
            a1: 2,
            z2: 1,
            a2: 2,
            q_mev: q_dd_gamma,
            branch_weight: dd_gamma,
            neutron_fraction: 0.0,
        },
        ReactionModel {
            id: "p+B11->3alpha",
            z1: 1,
            a1: 1,
            z2: 5,
            a2: 11,
            q_mev: 8.680,
            branch_weight: 1.0,
            neutron_fraction: 0.0,
        },
        ReactionModel {
            id: "B10+B10->O16+He4",
            z1: 5,
            a1: 10,
            z2: 5,
            a2: 10,
            q_mev: 26.413_308,
            branch_weight: 1.0,
            neutron_fraction: 0.0,
        },
        ReactionModel {
            id: "D+Li6->He4+He4",
            z1: 1,
            a1: 2,
            z2: 3,
            a2: 6,
            q_mev: 22.372_771,
            branch_weight: 1.0,
            neutron_fraction: 0.0,
        },
    ]
}

fn temperature_grid_kev() -> Vec<f64> {
    vec![10.0, 30.0, 100.0, 300.0, 1000.0]
}

fn annihilation_energy_mev(a1: u16, a2: u16) -> f64 {
    // First-order structural proxy: convert total rest mass of both nuclei to energy.
    // E ~= (A1 + A2) * m_N, ignoring O(1%) binding corrections.
    (a1 as f64 + a2 as f64) * NUCLEON_MASS_MEV
}

fn score_reaction(
    r: ReactionModel,
    t_kev: f64,
    mode: PairMode,
    neutron_penalty_kappa: f64,
) -> ReactionScore {
    let gamow = match mode {
        PairMode::MatterMatter | PairMode::AntimatterAntimatter => {
            gamow_with_charges(r.z1, r.z2, r.a1, r.a2, t_kev)
        }
        PairMode::MatterAntimatter => {
            // Opposite charges attract; no Coulomb tunneling suppression.
            1.0
        }
    };

    let fusion_raw = r.branch_weight * r.q_mev * gamow;
    let fusion_clean = fusion_raw * (-neutron_penalty_kappa.max(0.0) * r.neutron_fraction).exp();
    let annihilation_raw = match mode {
        PairMode::MatterAntimatter => annihilation_energy_mev(r.a1, r.a2),
        _ => 0.0,
    };

    ReactionScore {
        t_kev,
        gamow,
        fusion_raw,
        fusion_clean,
        annihilation_raw,
        total_raw: fusion_raw + annihilation_raw,
        total_clean: fusion_clean + annihilation_raw,
    }
}

fn main() {
    let out_dir = env::var("GUTOE_FUSION_MA_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/fusion_matter_antimatter".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");

    let reactions = model_set();
    let temps = temperature_grid_kev();
    let modes = [
        PairMode::MatterMatter,
        PairMode::MatterAntimatter,
        PairMode::AntimatterAntimatter,
    ];
    let neutron_penalty_kappa = env::var("GUTOE_FUSION_NEUTRON_PENALTY_KAPPA")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(5.0)
        .max(0.0);

    let mut csv = String::from(
        "mode,reaction_id,t_kev,z1,a1,z2,a2,q_fusion_mev,annihilation_mev,branch_weight,neutron_fraction,gamow,fusion_raw,fusion_clean,total_raw,total_clean\n",
    );
    let mut txt = String::new();
    txt.push_str("[fusion_matter_antimatter_report]\n");
    txt.push_str("fusion_raw   = branch_weight * Q_fusion * Gamow\n");
    txt.push_str("fusion_clean = fusion_raw * exp(-kappa * neutron_fraction)\n");
    txt.push_str("total_raw    = fusion_raw + annihilation_raw\n");
    txt.push_str("total_clean  = fusion_clean + annihilation_raw\n");
    txt.push_str(&format!(
        "neutron_penalty_kappa = {:.3}\n\n",
        neutron_penalty_kappa
    ));

    for mode in modes {
        txt.push_str(&format!("[mode: {}]\n", mode.id()));
        for &t in &temps {
            let mut scored: Vec<(ReactionModel, ReactionScore)> = reactions
                .iter()
                .copied()
                .map(|r| (r, score_reaction(r, t, mode, neutron_penalty_kappa)))
                .collect();
            for (r, s) in &scored {
                csv.push_str(&format!(
                    "{},{},{:.1},{},{},{},{},{:.6},{:.6},{:.12e},{:.3},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}\n",
                    mode.id(),
                    r.id,
                    s.t_kev,
                    r.z1,
                    r.a1,
                    r.z2,
                    r.a2,
                    r.q_mev,
                    s.annihilation_raw,
                    r.branch_weight,
                    r.neutron_fraction,
                    s.gamow,
                    s.fusion_raw,
                    s.fusion_clean,
                    s.total_raw,
                    s.total_clean
                ));
            }

            scored.sort_by(|a, b| b.1.fusion_raw.total_cmp(&a.1.fusion_raw));
            let (top_fusion_r, top_fusion_s) = scored[0];
            scored.sort_by(|a, b| b.1.total_raw.total_cmp(&a.1.total_raw));
            let (top_total_r, top_total_s) = scored[0];
            txt.push_str(&format!(
                "T_keV={:.1} | top_fusion={} (fusion_raw={:.6e}, gamow={:.3e}) | top_total={} (total_raw={:.6e}, annihilation={:.3e})\n",
                t,
                top_fusion_r.id,
                top_fusion_s.fusion_raw,
                top_fusion_s.gamow,
                top_total_r.id,
                top_total_s.total_raw,
                top_total_s.annihilation_raw
            ));
        }
        txt.push('\n');
    }

    let t_ref = 100.0;
    txt.push_str("[reference_T_100keV]\n");
    for mode in modes {
        let mut scored: Vec<(ReactionModel, ReactionScore)> = reactions
            .iter()
            .copied()
            .map(|r| (r, score_reaction(r, t_ref, mode, neutron_penalty_kappa)))
            .collect();
        scored.sort_by(|a, b| b.1.fusion_raw.total_cmp(&a.1.fusion_raw));
        let (best_fusion_r, best_fusion_s) = scored[0];
        scored.sort_by(|a, b| b.1.total_raw.total_cmp(&a.1.total_raw));
        let (best_total_r, best_total_s) = scored[0];
        txt.push_str(&format!(
            "mode={} | best_fusion={} (fusion_raw={:.6e}) | best_total={} (total_raw={:.6e}, annihilation={:.6e})\n",
            mode.id(),
            best_fusion_r.id,
            best_fusion_s.fusion_raw,
            best_total_r.id,
            best_total_s.total_raw,
            best_total_s.annihilation_raw
        ));
    }

    let csv_path = format!("{out_dir}/fusion_matter_antimatter.csv");
    let txt_path = format!("{out_dir}/fusion_matter_antimatter.txt");
    fs::write(&csv_path, csv).expect("write csv");
    fs::write(&txt_path, txt).expect("write txt");
    println!("wrote {}", csv_path);
    println!("wrote {}", txt_path);
}

