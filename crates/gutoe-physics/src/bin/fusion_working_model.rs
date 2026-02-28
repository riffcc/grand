use gutoe_physics::{gamow_factor, ALPHA_LEADING_ORDER};

const NUCLEON_MASS_MEV: f64 = 938.272_088_16;

#[derive(Clone, Copy, Debug)]
struct ReactionModel {
    id: &'static str,
    z1: u16,
    a1: u16,
    z2: u16,
    a2: u16,
    q_mev: f64,
    // Structural branch fraction/proxy weight in [0,1].
    branch_weight: f64,
    // Fraction of released channel energy carried by neutron-bearing products.
    neutron_fraction: f64,
    notes: &'static str,
}

#[derive(Clone, Debug)]
struct ReactionScore {
    t_kev: f64,
    gamow: f64,
    barrier_weighted_q: f64,
    raw_score: f64,
    clean_score: f64,
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
    // Leading radiative suppression proxy from earlier derivation lane:
    // B_gamma ~ alpha * (E_gamma / m_N)^3
    (ALPHA_LEADING_ORDER * (q_gamma_mev / NUCLEON_MASS_MEV).powi(3)).clamp(0.0, 1.0)
}

fn model_set() -> Vec<ReactionModel> {
    // Q values (MeV) from AME-based scan lane already run in repo context.
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
            notes: "baseline high-reactivity branch",
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
            notes: "aneutronic primary winner",
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
            notes: "dominant DD strong branch",
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
            notes: "dominant DD strong branch",
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
            notes: "radiative branch (strongly suppressed)",
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
            notes: "classic aneutronic candidate",
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
            notes: "high-Q but extreme Coulomb barrier",
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
            notes: "matrix high-Q light-fuel channel",
        },
    ]
}

fn temperature_grid_kev() -> Vec<f64> {
    vec![
        5.0, 10.0, 20.0, 30.0, 50.0, 80.0, 100.0, 150.0, 200.0, 300.0, 500.0, 1000.0, 2000.0,
    ]
}

fn score_reaction(r: ReactionModel, t_kev: f64, neutron_penalty_kappa: f64) -> ReactionScore {
    // Proxy score: branch_weight * Q * Gamow.
    // This is a structural ignition proxy, not an absolute cross-section.
    let g = gamow_with_charges(r.z1, r.z2, r.a1, r.a2, t_kev);
    let bq = r.branch_weight * r.q_mev;
    let raw = bq * g;
    let cleanliness = (-neutron_penalty_kappa.max(0.0) * r.neutron_fraction).exp();
    ReactionScore {
        t_kev,
        gamow: g,
        barrier_weighted_q: bq,
        raw_score: raw,
        clean_score: raw * cleanliness,
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_FUSION_MODEL_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/fusion_working_model".to_string());
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    let reactions = model_set();
    let temps = temperature_grid_kev();
    let neutron_penalty_kappa = std::env::var("GUTOE_FUSION_NEUTRON_PENALTY_KAPPA")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(5.0)
        .max(0.0);

    let mut csv = String::from(
        "reaction_id,t_kev,z1,a1,z2,a2,q_mev,branch_weight,neutron_fraction,gamow,barrier_weighted_q,raw_score,clean_score,notes\n",
    );

    let mut best_overall_raw: Option<(ReactionModel, ReactionScore)> = None;
    let mut best_overall_clean: Option<(ReactionModel, ReactionScore)> = None;
    let mut best_per_t_raw: Vec<(f64, ReactionModel, ReactionScore)> = Vec::new();
    let mut best_per_t_clean: Vec<(f64, ReactionModel, ReactionScore)> = Vec::new();

    for &t in &temps {
        let scores: Vec<(ReactionModel, ReactionScore)> = reactions
            .iter()
            .copied()
            .map(|r| {
                let s = score_reaction(r, t, neutron_penalty_kappa);
                (r, s)
            })
            .collect();

        for (r, s) in &scores {
            csv.push_str(&format!(
                "{},{:.1},{},{},{},{},{:.6},{:.12e},{:.6},{:.12e},{:.12e},{:.12e},{:.12e},{}\n",
                r.id,
                s.t_kev,
                r.z1,
                r.a1,
                r.z2,
                r.a2,
                r.q_mev,
                r.branch_weight,
                r.neutron_fraction,
                s.gamow,
                s.barrier_weighted_q,
                s.raw_score,
                s.clean_score,
                r.notes
            ));
        }

        let mut raw_sorted = scores.clone();
        raw_sorted.sort_by(|a, b| b.1.raw_score.total_cmp(&a.1.raw_score));
        if let Some((r, s)) = raw_sorted.first().cloned() {
            best_per_t_raw.push((t, r, s.clone()));
            if best_overall_raw
                .as_ref()
                .map(|(_, cur)| s.raw_score > cur.raw_score)
                .unwrap_or(true)
            {
                best_overall_raw = Some((r, s));
            }
        }

        let mut clean_sorted = scores;
        clean_sorted.sort_by(|a, b| b.1.clean_score.total_cmp(&a.1.clean_score));
        if let Some((r, s)) = clean_sorted.first().cloned() {
            best_per_t_clean.push((t, r, s.clone()));
            if best_overall_clean
                .as_ref()
                .map(|(_, cur)| s.clean_score > cur.clean_score)
                .unwrap_or(true)
            {
                best_overall_clean = Some((r, s));
            }
        }
    }

    let csv_path = format!("{out_dir}/fusion_working_model.csv");
    std::fs::write(&csv_path, csv).expect("write csv");

    let mut txt = String::new();
    txt.push_str("[fusion_working_model]\n");
    txt.push_str("raw_score   = branch_weight * Q * Gamow(Z1,Z2,mu,E_cm=T)\n");
    txt.push_str("clean_score = raw_score * exp(-neutron_penalty_kappa * neutron_fraction)\n");
    txt.push_str("(relative ignition/power proxy; not absolute sigma)\n\n");
    txt.push_str(&format!(
        "neutron_penalty_kappa = {:.3}\n\n",
        neutron_penalty_kappa
    ));

    // Branch suppression checkpoint for DD radiative lane.
    if let Some(ddg) = reactions.iter().find(|r| r.id == "D+D->gamma+He4") {
        txt.push_str(&format!(
            "dd_gamma_branch_weight = {:.12e} (inverse {:.3e})\n\n",
            ddg.branch_weight,
            1.0 / ddg.branch_weight.max(1e-30)
        ));
    }

    txt.push_str("[best_per_temperature_raw]\n");
    for (t, r, s) in &best_per_t_raw {
        txt.push_str(&format!(
            "T_keV={:.1}: {} | raw_score={:.6e} | Gamow={:.6e} | Q={:.3} | branch={:.3e}\n",
            t, r.id, s.raw_score, s.gamow, r.q_mev, r.branch_weight
        ));
    }
    txt.push_str("\n[best_per_temperature_clean]\n");
    for (t, r, s) in &best_per_t_clean {
        txt.push_str(&format!(
            "T_keV={:.1}: {} | clean_score={:.6e} | raw={:.6e} | neutron_fraction={:.3}\n",
            t, r.id, s.clean_score, s.raw_score, r.neutron_fraction
        ));
    }

    if let Some((r, s)) = &best_overall_raw {
        txt.push_str("\n[overall_best_raw]\n");
        txt.push_str(&format!(
            "reaction={}\nT_keV={:.1}\nraw_score={:.12e}\nGamow={:.12e}\nQ_mev={:.6}\nbranch_weight={:.12e}\nneutron_fraction={:.3}\nnotes={}\n",
            r.id,
            s.t_kev,
            s.raw_score,
            s.gamow,
            r.q_mev,
            r.branch_weight,
            r.neutron_fraction,
            r.notes
        ));
    }
    if let Some((r, s)) = &best_overall_clean {
        txt.push_str("\n[overall_best_clean]\n");
        txt.push_str(&format!(
            "reaction={}\nT_keV={:.1}\nclean_score={:.12e}\nraw_score={:.12e}\nneutron_fraction={:.3}\nnotes={}\n",
            r.id, s.t_kev, s.clean_score, s.raw_score, r.neutron_fraction, r.notes
        ));
    }

    // Explicit comparison rows requested by current direction.
    let t_ref = 100.0;
    let d_he3 = score_reaction(
        reactions
            .iter()
            .copied()
            .find(|r| r.id == "D+He3->p+He4")
            .expect("D+He3 exists"),
        t_ref,
        neutron_penalty_kappa,
    );
    let b10b10 = score_reaction(
        reactions
            .iter()
            .copied()
            .find(|r| r.id == "B10+B10->O16+He4")
            .expect("B10+B10 exists"),
        t_ref,
        neutron_penalty_kappa,
    );
    let ratio = b10b10.raw_score / d_he3.raw_score.max(1e-300);
    txt.push_str("\n[barrier_reality_check]\n");
    txt.push_str(&format!(
        "T_keV={:.1}\nraw_score_ratio_B10B10_over_DHe3={:.12e}\n",
        t_ref, ratio
    ));

    let txt_path = format!("{out_dir}/fusion_working_model.txt");
    std::fs::write(&txt_path, txt).expect("write txt");

    println!("wrote {}", csv_path);
    println!("wrote {}", txt_path);
}
