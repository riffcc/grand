//! Yukawa acid test report.
//!
//! Goal: quantify what is already constrained by the current Cl(1,3) lanes
//! and where absolute Yukawa closure is still missing physics.

use gutoe_em::alpha::{
    lepton_masses_from_electron_structural_alpha, ALPHA_INVERSE_PHYSICAL, MP_ME_CLIFFORD,
};
use gutoe_em::flavor::neutrino_absolute_masses_from_texture;
use gutoe_em::weak::{electron_mass_from_proton_anchor, electroweak_vev_from_fermi};
use serde::Serialize;
use std::fs::{self, File};
use std::io::Write;

const SQRT2: f64 = std::f64::consts::SQRT_2;
const G_F: f64 = 1.166_378_7e-5; // GeV^-2

// Charged-lepton references (MeV).
const ME_OBS: f64 = 0.510_998_95;
const MMU_OBS: f64 = 105.658_375_5;
const MTAU_OBS: f64 = 1776.93;

// Quark reference masses (MeV), coarse PDG-like central values.
const MU_OBS: f64 = 2.16;
const MD_OBS: f64 = 4.67;
const MS_OBS: f64 = 93.0;
const MC_OBS: f64 = 1270.0;
const MB_OBS: f64 = 4180.0;
const MT_OBS: f64 = 172_760.0;

#[derive(Debug, Clone, Copy, Serialize)]
struct MassRow {
    name: &'static str,
    mass_pred_mev: f64,
    mass_ref_mev: f64,
    rel_err: f64,
    yukawa_pred: f64,
    yukawa_ref: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct RatioRow {
    name: &'static str,
    pred: f64,
    ref_from_masses: f64,
    rel_err_vs_ref: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct RatioConsistency {
    name: &'static str,
    lhs: f64,
    rhs: f64,
    rel_mismatch: f64,
}

#[derive(Debug, Clone, Serialize)]
struct Report {
    alpha_inv_physical: f64,
    mp_me_structural: f64,
    vev_gev: f64,
    charged_leptons: Vec<MassRow>,
    quark_ratios: Vec<RatioRow>,
    ratio_consistency_cycles: Vec<RatioConsistency>,
    quarks_one_anchor_d_lstsq: Vec<MassRow>,
    neutrino_lane: NeutrinoLane,
    summary: Summary,
}

#[derive(Debug, Clone, Serialize)]
struct NeutrinoLane {
    m1_ev: f64,
    m2_ev: f64,
    m3_ev: f64,
    sum_ev: f64,
    y1: f64,
    y2: f64,
    y3: f64,
    splitting_ratio_32_over_21: f64,
}

#[derive(Debug, Clone, Serialize)]
struct Summary {
    lepton_max_rel_err_pct: f64,
    quark_ratio_max_rel_err_pct: f64,
    quark_cycle_max_mismatch_pct: f64,
    quark_lstsq_max_rel_err_pct: f64,
    dof_note: String,
}

fn rel_err(pred: f64, reference: f64) -> f64 {
    if reference == 0.0 {
        0.0
    } else {
        (pred - reference).abs() / reference.abs()
    }
}

fn yukawa_from_mev(mass_mev: f64, vev_gev: f64) -> f64 {
    let mass_gev = mass_mev / 1000.0;
    SQRT2 * mass_gev / vev_gev
}

fn yukawa_from_ev(mass_ev: f64, vev_gev: f64) -> f64 {
    let mass_gev = mass_ev * 1.0e-9;
    SQRT2 * mass_gev / vev_gev
}

fn build_quark_ratio_rows() -> (Vec<RatioRow>, [f64; 7]) {
    // Shared structural lanes from existing report/theorem chain.
    let lambda_inv2 = 19.0; // (16 + 3)
    let c_inf = 67.0 / 66.0;

    let mu_over_md = 8.0 / 17.0;
    let mc_over_ms = (13.0 / 21.0) * lambda_inv2 * c_inf;
    let mt_over_mb = (13.0 / 6.0) * lambda_inv2 * c_inf;
    let mc_over_mu = (8.0 / 5.0) * lambda_inv2 * lambda_inv2 * c_inf;
    let mt_over_mc = 8.0 * 17.0;
    let ms_over_md = lambda_inv2;
    let mb_over_ms = (8.0 / 3.0) * lambda_inv2 * c_inf;

    let rows = vec![
        RatioRow {
            name: "m_u/m_d",
            pred: mu_over_md,
            ref_from_masses: MU_OBS / MD_OBS,
            rel_err_vs_ref: rel_err(mu_over_md, MU_OBS / MD_OBS),
        },
        RatioRow {
            name: "m_c/m_s",
            pred: mc_over_ms,
            ref_from_masses: MC_OBS / MS_OBS,
            rel_err_vs_ref: rel_err(mc_over_ms, MC_OBS / MS_OBS),
        },
        RatioRow {
            name: "m_t/m_b",
            pred: mt_over_mb,
            ref_from_masses: MT_OBS / MB_OBS,
            rel_err_vs_ref: rel_err(mt_over_mb, MT_OBS / MB_OBS),
        },
        RatioRow {
            name: "m_c/m_u",
            pred: mc_over_mu,
            ref_from_masses: MC_OBS / MU_OBS,
            rel_err_vs_ref: rel_err(mc_over_mu, MC_OBS / MU_OBS),
        },
        RatioRow {
            name: "m_t/m_c",
            pred: mt_over_mc,
            ref_from_masses: MT_OBS / MC_OBS,
            rel_err_vs_ref: rel_err(mt_over_mc, MT_OBS / MC_OBS),
        },
        RatioRow {
            name: "m_s/m_d",
            pred: ms_over_md,
            ref_from_masses: MS_OBS / MD_OBS,
            rel_err_vs_ref: rel_err(ms_over_md, MS_OBS / MD_OBS),
        },
        RatioRow {
            name: "m_b/m_s",
            pred: mb_over_ms,
            ref_from_masses: MB_OBS / MS_OBS,
            rel_err_vs_ref: rel_err(mb_over_ms, MB_OBS / MS_OBS),
        },
    ];

    (
        rows,
        [
            mu_over_md,
            mc_over_ms,
            mt_over_mb,
            mc_over_mu,
            mt_over_mc,
            ms_over_md,
            mb_over_ms,
        ],
    )
}

fn ratio_consistency(r: [f64; 7]) -> Vec<RatioConsistency> {
    let mu_md = r[0];
    let mc_ms = r[1];
    let mt_mb = r[2];
    let mc_mu = r[3];
    let mt_mc = r[4];
    let ms_md = r[5];
    let mb_ms = r[6];

    // Cycles that should be 1 in a fully self-consistent ratio graph.
    let lhs1 = mc_mu * mu_md;
    let rhs1 = mc_ms * ms_md;

    let lhs2 = mt_mc * mc_ms;
    let rhs2 = mt_mb * mb_ms;

    let lhs3 = mt_mc * mc_mu * mu_md;
    let rhs3 = mt_mb * mb_ms * ms_md;

    vec![
        RatioConsistency {
            name: "(mc/mu)*(mu/md) ?= (mc/ms)*(ms/md)",
            lhs: lhs1,
            rhs: rhs1,
            rel_mismatch: rel_err(lhs1, rhs1),
        },
        RatioConsistency {
            name: "(mt/mc)*(mc/ms) ?= (mt/mb)*(mb/ms)",
            lhs: lhs2,
            rhs: rhs2,
            rel_mismatch: rel_err(lhs2, rhs2),
        },
        RatioConsistency {
            name: "(mt/mc)*(mc/mu)*(mu/md) ?= (mt/mb)*(mb/ms)*(ms/md)",
            lhs: lhs3,
            rhs: rhs3,
            rel_mismatch: rel_err(lhs3, rhs3),
        },
    ]
}

fn solve_linear_system_5(mut a: [[f64; 5]; 5], mut b: [f64; 5]) -> Option<[f64; 5]> {
    for col in 0..5 {
        // Pivot.
        let mut pivot = col;
        let mut pivot_abs = a[col][col].abs();
        for row in (col + 1)..5 {
            let v = a[row][col].abs();
            if v > pivot_abs {
                pivot_abs = v;
                pivot = row;
            }
        }
        if pivot_abs < 1.0e-14 {
            return None;
        }
        if pivot != col {
            a.swap(col, pivot);
            b.swap(col, pivot);
        }

        // Normalize pivot row.
        let diag = a[col][col];
        for j in col..5 {
            a[col][j] /= diag;
        }
        b[col] /= diag;

        // Eliminate.
        for row in 0..5 {
            if row == col {
                continue;
            }
            let factor = a[row][col];
            if factor == 0.0 {
                continue;
            }
            for j in col..5 {
                a[row][j] -= factor * a[col][j];
            }
            b[row] -= factor * b[col];
        }
    }
    Some(b)
}

fn quark_masses_one_anchor_d_lstsq(
    r: [f64; 7],
    d_anchor_mev: f64,
    vev_gev: f64,
) -> Vec<MassRow> {
    // Unknown z = [x_u, x_s, x_c, x_b, x_t] with x = ln(mass_mev).
    // Equations:
    //  1) x_u             = x_d + ln(mu/md)
    //  2) x_c - x_s       = ln(mc/ms)
    //  3) x_t - x_b       = ln(mt/mb)
    //  4) x_c - x_u       = ln(mc/mu)
    //  5) x_t - x_c       = ln(mt/mc)
    //  6) x_s             = x_d + ln(ms/md)
    //  7) x_b - x_s       = ln(mb/ms)
    let x_d = d_anchor_mev.ln();
    let ln = |v: f64| v.ln();

    let rows = [
        ([1.0, 0.0, 0.0, 0.0, 0.0], x_d + ln(r[0])),
        ([0.0, -1.0, 1.0, 0.0, 0.0], ln(r[1])),
        ([0.0, 0.0, 0.0, -1.0, 1.0], ln(r[2])),
        ([-1.0, 0.0, 1.0, 0.0, 0.0], ln(r[3])),
        ([0.0, 0.0, -1.0, 0.0, 1.0], ln(r[4])),
        ([0.0, 1.0, 0.0, 0.0, 0.0], x_d + ln(r[5])),
        ([0.0, -1.0, 0.0, 1.0, 0.0], ln(r[6])),
    ];

    // Normal equations AtA z = Atb.
    let mut ata = [[0.0_f64; 5]; 5];
    let mut atb = [0.0_f64; 5];
    for (a, b) in rows {
        for i in 0..5 {
            atb[i] += a[i] * b;
            for j in 0..5 {
                ata[i][j] += a[i] * a[j];
            }
        }
    }

    let z = solve_linear_system_5(ata, atb).unwrap_or([0.0; 5]);
    let masses = [z[0].exp(), d_anchor_mev, z[1].exp(), z[2].exp(), z[3].exp(), z[4].exp()];

    let refs = [MU_OBS, MD_OBS, MS_OBS, MC_OBS, MB_OBS, MT_OBS];
    let names = ["u", "d", "s", "c", "b", "t"];

    names
        .iter()
        .enumerate()
        .map(|(i, &name)| MassRow {
            name,
            mass_pred_mev: masses[i],
            mass_ref_mev: refs[i],
            rel_err: rel_err(masses[i], refs[i]),
            yukawa_pred: yukawa_from_mev(masses[i], vev_gev),
            yukawa_ref: yukawa_from_mev(refs[i], vev_gev),
        })
        .collect()
}

fn main() {
    let out_dir =
        std::env::var("GUTOE_YUKAWA_ACID_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/yukawa_acid_test_report.txt");
    let json_path = format!("{out_dir}/yukawa_acid_test_report.json");

    let vev_gev = electroweak_vev_from_fermi(G_F);

    // Charged lepton lane: one absolute anchor from mp/me structural lane.
    let m_e_pred = electron_mass_from_proton_anchor();
    let [me_s, mmu_s, mtau_s] = lepton_masses_from_electron_structural_alpha(m_e_pred);
    let charged_leptons = vec![
        MassRow {
            name: "e",
            mass_pred_mev: me_s,
            mass_ref_mev: ME_OBS,
            rel_err: rel_err(me_s, ME_OBS),
            yukawa_pred: yukawa_from_mev(me_s, vev_gev),
            yukawa_ref: yukawa_from_mev(ME_OBS, vev_gev),
        },
        MassRow {
            name: "mu",
            mass_pred_mev: mmu_s,
            mass_ref_mev: MMU_OBS,
            rel_err: rel_err(mmu_s, MMU_OBS),
            yukawa_pred: yukawa_from_mev(mmu_s, vev_gev),
            yukawa_ref: yukawa_from_mev(MMU_OBS, vev_gev),
        },
        MassRow {
            name: "tau",
            mass_pred_mev: mtau_s,
            mass_ref_mev: MTAU_OBS,
            rel_err: rel_err(mtau_s, MTAU_OBS),
            yukawa_pred: yukawa_from_mev(mtau_s, vev_gev),
            yukawa_ref: yukawa_from_mev(MTAU_OBS, vev_gev),
        },
    ];

    let (quark_ratios, rvals) = build_quark_ratio_rows();
    let ratio_consistency_cycles = ratio_consistency(rvals);
    let quarks_one_anchor_d_lstsq = quark_masses_one_anchor_d_lstsq(rvals, MD_OBS, vev_gev);

    let nu_abs = neutrino_absolute_masses_from_texture();
    let neutrino_lane = NeutrinoLane {
        m1_ev: nu_abs.m1_ev,
        m2_ev: nu_abs.m2_ev,
        m3_ev: nu_abs.m3_ev,
        sum_ev: nu_abs.sum_ev,
        y1: yukawa_from_ev(nu_abs.m1_ev, vev_gev),
        y2: yukawa_from_ev(nu_abs.m2_ev, vev_gev),
        y3: yukawa_from_ev(nu_abs.m3_ev, vev_gev),
        splitting_ratio_32_over_21: nu_abs.splitting_ratio_32_over_21,
    };

    let lepton_max_rel_err = charged_leptons
        .iter()
        .map(|r| r.rel_err)
        .fold(0.0_f64, f64::max);
    let quark_ratio_max_rel_err = quark_ratios
        .iter()
        .map(|r| r.rel_err_vs_ref)
        .fold(0.0_f64, f64::max);
    let quark_cycle_max_mismatch = ratio_consistency_cycles
        .iter()
        .map(|r| r.rel_mismatch)
        .fold(0.0_f64, f64::max);
    let quark_lstsq_max_rel_err = quarks_one_anchor_d_lstsq
        .iter()
        .map(|r| r.rel_err)
        .fold(0.0_f64, f64::max);

    let summary = Summary {
        lepton_max_rel_err_pct: lepton_max_rel_err * 100.0,
        quark_ratio_max_rel_err_pct: quark_ratio_max_rel_err * 100.0,
        quark_cycle_max_mismatch_pct: quark_cycle_max_mismatch * 100.0,
        quark_lstsq_max_rel_err_pct: quark_lstsq_max_rel_err * 100.0,
        dof_note: String::from(
            "Current lane: strong ratio structure + one-anchor absolute closure; quark ratio cycles are not yet fully integrable (missing higher-order physics).",
        ),
    };

    let report = Report {
        alpha_inv_physical: ALPHA_INVERSE_PHYSICAL,
        mp_me_structural: MP_ME_CLIFFORD as f64,
        vev_gev,
        charged_leptons,
        quark_ratios,
        ratio_consistency_cycles,
        quarks_one_anchor_d_lstsq,
        neutrino_lane,
        summary,
    };

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[yukawa_acid_test]").expect("write");
    writeln!(txt, "alpha_inv_physical = {:.12}", report.alpha_inv_physical).expect("write");
    writeln!(txt, "mp_me_structural = {:.6}", report.mp_me_structural).expect("write");
    writeln!(txt, "vev_gev = {:.9}", report.vev_gev).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[charged_leptons_structural_alpha]").expect("write");
    for row in &report.charged_leptons {
        writeln!(
            txt,
            "{}: m_pred={:.9} MeV  m_ref={:.9} MeV  rel_err={:.6}%  y_pred={:.12e}  y_ref={:.12e}",
            row.name,
            row.mass_pred_mev,
            row.mass_ref_mev,
            row.rel_err * 100.0,
            row.yukawa_pred,
            row.yukawa_ref
        )
        .expect("write");
    }
    writeln!(txt).expect("write");

    writeln!(txt, "[quark_ratio_lanes]").expect("write");
    for row in &report.quark_ratios {
        writeln!(
            txt,
            "{}: pred={:.9}  ref={:.9}  rel_err={:.6}%",
            row.name,
            row.pred,
            row.ref_from_masses,
            row.rel_err_vs_ref * 100.0
        )
        .expect("write");
    }
    writeln!(txt).expect("write");

    writeln!(txt, "[quark_ratio_cycle_consistency]").expect("write");
    for c in &report.ratio_consistency_cycles {
        writeln!(
            txt,
            "{}: lhs={:.9} rhs={:.9} mismatch={:.6}%",
            c.name,
            c.lhs,
            c.rhs,
            c.rel_mismatch * 100.0
        )
        .expect("write");
    }
    writeln!(txt).expect("write");

    writeln!(txt, "[quark_masses_one_anchor_d_lstsq]").expect("write");
    for row in &report.quarks_one_anchor_d_lstsq {
        writeln!(
            txt,
            "{}: m_pred={:.6} MeV  m_ref={:.6} MeV  rel_err={:.6}%  y_pred={:.12e}  y_ref={:.12e}",
            row.name,
            row.mass_pred_mev,
            row.mass_ref_mev,
            row.rel_err * 100.0,
            row.yukawa_pred,
            row.yukawa_ref
        )
        .expect("write");
    }
    writeln!(txt).expect("write");

    writeln!(txt, "[neutrino_lane]").expect("write");
    writeln!(
        txt,
        "m1={:.12e} eV  m2={:.12e} eV  m3={:.12e} eV  sum={:.12e} eV",
        report.neutrino_lane.m1_ev,
        report.neutrino_lane.m2_ev,
        report.neutrino_lane.m3_ev,
        report.neutrino_lane.sum_ev
    )
    .expect("write");
    writeln!(
        txt,
        "y1={:.12e} y2={:.12e} y3={:.12e}  split_32_over_21={:.6}",
        report.neutrino_lane.y1,
        report.neutrino_lane.y2,
        report.neutrino_lane.y3,
        report.neutrino_lane.splitting_ratio_32_over_21
    )
    .expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[summary]").expect("write");
    writeln!(
        txt,
        "lepton_max_rel_err_pct = {:.6}",
        report.summary.lepton_max_rel_err_pct
    )
    .expect("write");
    writeln!(
        txt,
        "quark_ratio_max_rel_err_pct = {:.6}",
        report.summary.quark_ratio_max_rel_err_pct
    )
    .expect("write");
    writeln!(
        txt,
        "quark_cycle_max_mismatch_pct = {:.6}",
        report.summary.quark_cycle_max_mismatch_pct
    )
    .expect("write");
    writeln!(
        txt,
        "quark_lstsq_max_rel_err_pct = {:.6}",
        report.summary.quark_lstsq_max_rel_err_pct
    )
    .expect("write");
    writeln!(txt, "dof_note = {}", report.summary.dof_note).expect("write");

    let json = serde_json::to_string_pretty(&report).expect("serialize");
    fs::write(&json_path, json).expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "summary: lepton_max={:.3}% quark_ratio_max={:.3}% quark_cycle_max={:.3}% quark_lstsq_max={:.3}%",
        report.summary.lepton_max_rel_err_pct,
        report.summary.quark_ratio_max_rel_err_pct,
        report.summary.quark_cycle_max_mismatch_pct,
        report.summary.quark_lstsq_max_rel_err_pct,
    );
}
