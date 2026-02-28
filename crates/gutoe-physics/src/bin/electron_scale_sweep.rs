use gutoe_physics::constants::{
    ALPHA_LEADING_ORDER, CLIFFORD_STATE_COUNT_STRUCTURAL, DARK_TO_VISIBLE_COUNT_RATIO,
    DARK_TO_VISIBLE_GEOMETRIC_RATIO, LAMBDA_QG, PLANCK_MASS, VISIBLE_STATE_COUNT_STRUCTURAL,
};
use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

const ELECTRON_MASS_MEV_OBS: f64 = 0.510_998_950;
const PROTON_MASS_MEV_OBS: f64 = 938.272_088_16;
const KG_TO_MEV: f64 = 5.609_588_603e29;
const C_INF: f64 = 67.0 / 66.0;
const MP_ME_STRUCT: f64 = 1836.0;

#[derive(Clone)]
struct Candidate {
    expr: String,
    f: f64,
    op_count: u32,
    m_pred_mev: f64,
    rel_err: f64,
    proton_pred_mev: f64,
    proton_rel_err: f64,
}

fn main() {
    let out_dir = std::env::var("GUTOE_ELECTRON_SWEEP_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/electron_scale_sweep".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let planck_mev = PLANCK_MASS * KG_TO_MEV;
    let target_ratio = ELECTRON_MASS_MEV_OBS / planck_mev;

    let alpha = ALPHA_LEADING_ORDER;
    let lambda = LAMBDA_QG;
    let geom_ratio = DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let c_inf = C_INF;
    let dv = DARK_TO_VISIBLE_COUNT_RATIO;
    let grade1_over_clifford = 4.0 / CLIFFORD_STATE_COUNT_STRUCTURAL;
    let one_over_visible = 1.0 / VISIBLE_STATE_COUNT_STRUCTURAL;

    let mut cands: Vec<Candidate> = Vec::new();

    // Hand-written hypotheses first (easy to audit).
    let manual: [(&str, f64, u32); 8] = [
        ("alpha^10", alpha.powi(10), 1),
        ("alpha^10 * lambda", alpha.powi(10) * lambda, 3),
        ("alpha^10 * dark_to_visible(5/11)", alpha.powi(10) * dv, 3),
        ("alpha^10 * grade1/clifford(1/4)", alpha.powi(10) * grade1_over_clifford, 4),
        (
            "alpha^10 * dark_geom(60/11)^(-1)",
            alpha.powi(10) * geom_ratio.powi(-1),
            4,
        ),
        ("alpha^10 * c_inf", alpha.powi(10) * c_inf, 3),
        ("alpha^10 * c_inf * (1/visible)", alpha.powi(10) * c_inf * one_over_visible, 5),
        ("alpha^9 * lambda^2 * c_inf", alpha.powi(9) * lambda.powi(2) * c_inf, 5),
    ];
    for (expr, f, op_count) in manual {
        let m_pred = planck_mev * f;
        let rel_err = (m_pred - ELECTRON_MASS_MEV_OBS) / ELECTRON_MASS_MEV_OBS;
        let proton_pred = m_pred * MP_ME_STRUCT;
        let proton_rel_err = (proton_pred - PROTON_MASS_MEV_OBS) / PROTON_MASS_MEV_OBS;
        cands.push(Candidate {
            expr: expr.to_string(),
            f,
            op_count,
            m_pred_mev: m_pred,
            rel_err,
            proton_pred_mev: proton_pred,
            proton_rel_err,
        });
    }

    // Structured brute-force over shared primitives with small integer exponents.
    for p_alpha in 6..=14 {
        for p_lambda in 0..=8 {
            for p_geom in -4..=2 {
                for p_cinf in -2..=2 {
                    for p_dv in -3..=3 {
                        let f = alpha.powi(p_alpha)
                            * lambda.powi(p_lambda)
                            * geom_ratio.powi(p_geom)
                            * c_inf.powi(p_cinf)
                            * dv.powi(p_dv);
                        if !f.is_finite() || f <= 0.0 {
                            continue;
                        }
                        let expr = format!(
                            "alpha^{p_alpha} * lambda^{p_lambda} * (60/11)^{p_geom} * (67/66)^{p_cinf} * (5/11)^{p_dv}"
                        );
                        let nonzero = [p_alpha, p_lambda, p_geom, p_cinf, p_dv]
                            .into_iter()
                            .filter(|p| *p != 0)
                            .count() as u32;
                        let op_count = if nonzero == 0 { 0 } else { nonzero + (nonzero - 1) };
                        let m_pred = planck_mev * f;
                        let rel_err = (m_pred - ELECTRON_MASS_MEV_OBS) / ELECTRON_MASS_MEV_OBS;
                        let proton_pred = m_pred * MP_ME_STRUCT;
                        let proton_rel_err = (proton_pred - PROTON_MASS_MEV_OBS) / PROTON_MASS_MEV_OBS;
                        cands.push(Candidate {
                            expr,
                            f,
                            op_count,
                            m_pred_mev: m_pred,
                            rel_err,
                            proton_pred_mev: proton_pred,
                            proton_rel_err,
                        });
                    }
                }
            }
        }
    }

    // Kill gate: expression complexity.
    cands.retain(|c| c.op_count <= 8);

    cands.sort_by(|a, b| {
        a.rel_err
            .abs()
            .partial_cmp(&b.rel_err.abs())
            .unwrap_or(Ordering::Equal)
    });

    let top_n = 25usize.min(cands.len());
    let survivors: Vec<&Candidate> = cands
        .iter()
        .filter(|c| c.rel_err.abs() <= 0.05)
        .take(5)
        .collect();
    let top = &cands[..top_n];
    let best = &cands[0];

    let txt_path = out.join("electron_scale_sweep.txt");
    let json_path = out.join("electron_scale_sweep.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[inputs]").expect("write");
    writeln!(txt, "planck_mass_kg = {:.12e}", PLANCK_MASS).expect("write");
    writeln!(txt, "planck_mass_mev = {:.12e}", planck_mev).expect("write");
    writeln!(txt, "electron_mass_mev_obs = {:.12}", ELECTRON_MASS_MEV_OBS).expect("write");
    writeln!(txt, "target_ratio_me_over_mplanck = {:.12e}", target_ratio).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[best]").expect("write");
    writeln!(txt, "expr = {}", best.expr).expect("write");
    writeln!(txt, "op_count = {}", best.op_count).expect("write");
    writeln!(txt, "F = {:.12e}", best.f).expect("write");
    writeln!(txt, "predicted_electron_mev = {:.12}", best.m_pred_mev).expect("write");
    writeln!(txt, "relative_error = {:.12e}", best.rel_err).expect("write");
    writeln!(txt, "predicted_proton_mev_via_1836 = {:.12}", best.proton_pred_mev).expect("write");
    writeln!(txt, "proton_relative_error = {:.12e}", best.proton_rel_err).expect("write");
    writeln!(
        txt,
        "status = {}",
        if best.rel_err.abs() < 1.0e-3 {
            "promote_lean"
        } else if best.rel_err.abs() < 1.0e-2 {
            "promote_investigation"
        } else if best.rel_err.abs() > 5.0e-2 {
            "killed_accuracy"
        } else {
            "review"
        }
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[survivors_top5_under_5pct]").expect("write");
    writeln!(txt, "rank,expr,op_count,predicted_electron_mev,relative_error,predicted_proton_mev,proton_relative_error,status").expect("write");
    for (i, c) in survivors.iter().enumerate() {
        let status = if c.rel_err.abs() < 1.0e-3 {
            "promote_lean"
        } else if c.rel_err.abs() < 1.0e-2 {
            "promote_investigation"
        } else {
            "review"
        };
        writeln!(
            txt,
            "{},{},{},{:.12},{:.12e},{:.12},{:.12e},{}",
            i + 1,
            c.expr,
            c.op_count,
            c.m_pred_mev,
            c.rel_err,
            c.proton_pred_mev,
            c.proton_rel_err,
            status
        ).expect("write");
    }
    writeln!(txt).expect("write");
    writeln!(txt, "[top_{}]", top_n).expect("write");
    writeln!(txt, "rank,expr,op_count,F,predicted_electron_mev,relative_error").expect("write");
    for (i, c) in top.iter().enumerate() {
        writeln!(
            txt,
            "{},{},{},{:.12e},{:.12},{:.12e}",
            i + 1,
            c.expr,
            c.op_count,
            c.f,
            c.m_pred_mev,
            c.rel_err
        )
        .expect("write");
    }

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"inputs\": {{\"planck_mass_kg\": {:.12e}, \"planck_mass_mev\": {:.12e}, \"electron_mass_mev_obs\": {:.12}, \"target_ratio_me_over_mplanck\": {:.12e}}},\n  \"best\": {{\"expr\": \"{}\", \"op_count\": {}, \"F\": {:.12e}, \"predicted_electron_mev\": {:.12}, \"relative_error\": {:.12e}, \"predicted_proton_mev_via_1836\": {:.12}, \"proton_relative_error\": {:.12e}}},\n  \"survivors_top5_under_5pct\": [",
        PLANCK_MASS,
        planck_mev,
        ELECTRON_MASS_MEV_OBS,
        target_ratio,
        best.expr.replace('"', "'"),
        best.op_count,
        best.f,
        best.m_pred_mev,
        best.rel_err,
        best.proton_pred_mev,
        best.proton_rel_err
    )
    .expect("write");
    for (i, c) in survivors.iter().enumerate() {
        let comma = if i + 1 == survivors.len() { "" } else { "," };
        let status = if c.rel_err.abs() < 1.0e-3 {
            "promote_lean"
        } else if c.rel_err.abs() < 1.0e-2 {
            "promote_investigation"
        } else {
            "review"
        };
        writeln!(
            json,
            "    {{\"rank\": {}, \"expr\": \"{}\", \"op_count\": {}, \"F\": {:.12e}, \"predicted_electron_mev\": {:.12}, \"relative_error\": {:.12e}, \"predicted_proton_mev_via_1836\": {:.12}, \"proton_relative_error\": {:.12e}, \"status\": \"{}\"}}{}",
            i + 1,
            c.expr.replace('"', "'"),
            c.op_count,
            c.f,
            c.m_pred_mev,
            c.rel_err,
            c.proton_pred_mev,
            c.proton_rel_err,
            status,
            comma
        )
        .expect("write row");
    }
    writeln!(json, "  ],\n  \"top\": [").expect("write");
    for (i, c) in top.iter().enumerate() {
        let comma = if i + 1 == top.len() { "" } else { "," };
        writeln!(
            json,
            "    {{\"rank\": {}, \"expr\": \"{}\", \"op_count\": {}, \"F\": {:.12e}, \"predicted_electron_mev\": {:.12}, \"relative_error\": {:.12e}}}{}",
            i + 1,
            c.expr.replace('"', "'"),
            c.op_count,
            c.f,
            c.m_pred_mev,
            c.rel_err,
            comma
        )
        .expect("write row");
    }
    writeln!(json, "  ]\n}}").expect("write end");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "best: {} | m_e_pred={:.9} MeV | rel_err={:.3e}",
        best.expr, best.m_pred_mev, best.rel_err
    );
}
