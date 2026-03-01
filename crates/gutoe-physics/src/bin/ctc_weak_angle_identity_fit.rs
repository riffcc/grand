//! Weak-angle identity fitter against public multi-point anchors.
//!
//! Goal:
//! - Test candidate identities for sin^2(theta_W) against measured anchors.
//! - Quantify whether a correction/running term rescues the leading-order miss.
//!
//! Scope:
//! - Coarse weighted least-squares with publicly cited anchors.
//! - Not a full global electroweak fit.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

const MZ_GEV: f64 = 91.1876;
const A_BASE: f64 = 3.0 / 13.0;
const A_SHIFTED_RATIONAL: f64 = 3.0 / 13.0 + 1.0 / (13.0 * 13.0 * 13.0);
const ALPHA_INV_PDG_2025: f64 = 137.035_999_177;

#[derive(Clone, Copy)]
struct Point {
    name: &'static str,
    q_gev: f64,
    value: f64,
    sigma: f64,
    source: &'static str,
    note: &'static str,
}

#[derive(Clone, Copy)]
struct FitSummary {
    chi2: f64,
    dof: f64,
    red_chi2: f64,
    params: [f64; 2], // interpretation depends on model
    used_params: usize,
}

fn ln_ratio(q: f64) -> f64 {
    (MZ_GEV / q).ln()
}

fn eval_model(model: &str, p: &Point, params: [f64; 2]) -> f64 {
    let alpha = 1.0 / ALPHA_INV_PDG_2025;
    let delta1_fixed = alpha * 10.0_f64.ln() / (4.0 * std::f64::consts::PI);
    match model {
        // fixed candidates
        "base_fixed" => A_BASE,
        "rational_fixed" => A_SHIFTED_RATIONAL,
        // shifted candidates
        "base_plus_delta0" => A_BASE + params[0],
        "base_plus_delta0_delta1log" => A_BASE + params[0] + params[1] * ln_ratio(p.q_gev),
        "rational_plus_delta1log" => A_SHIFTED_RATIONAL + params[1] * ln_ratio(p.q_gev),
        "rational_plus_alpha_ln10_over_4pi_log_fixed" => {
            A_SHIFTED_RATIONAL + delta1_fixed * ln_ratio(p.q_gev)
        }
        _ => A_BASE,
    }
}

fn chi2_for(model: &str, params: [f64; 2], points: &[Point]) -> f64 {
    points
        .iter()
        .map(|pt| {
            let pred = eval_model(model, pt, params);
            let r = (pred - pt.value) / pt.sigma;
            r * r
        })
        .sum()
}

fn fit_base_plus_delta0(points: &[Point]) -> FitSummary {
    // Weighted mean of (y - A_BASE)
    let mut sw = 0.0;
    let mut swy = 0.0;
    for p in points {
        let w = 1.0 / (p.sigma * p.sigma);
        sw += w;
        swy += w * (p.value - A_BASE);
    }
    let delta0 = if sw > 0.0 { swy / sw } else { 0.0 };
    let params = [delta0, 0.0];
    let chi2 = chi2_for("base_plus_delta0", params, points);
    let dof = (points.len() as f64 - 1.0).max(1.0);
    FitSummary {
        chi2,
        dof,
        red_chi2: chi2 / dof,
        params,
        used_params: 1,
    }
}

fn fit_base_plus_delta0_delta1log(points: &[Point]) -> FitSummary {
    // Weighted linear fit of y-A_BASE = d0 + d1*x
    let mut s = 0.0;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;

    for p in points {
        let x = ln_ratio(p.q_gev);
        let y = p.value - A_BASE;
        let w = 1.0 / (p.sigma * p.sigma);
        s += w;
        sx += w * x;
        sy += w * y;
        sxx += w * x * x;
        sxy += w * x * y;
    }

    let det = s * sxx - sx * sx;
    let (d0, d1) = if det.abs() > 1e-30 {
        ((sy * sxx - sxy * sx) / det, (s * sxy - sx * sy) / det)
    } else {
        (0.0, 0.0)
    };
    let params = [d0, d1];
    let chi2 = chi2_for("base_plus_delta0_delta1log", params, points);
    let dof = (points.len() as f64 - 2.0).max(1.0);
    FitSummary {
        chi2,
        dof,
        red_chi2: chi2 / dof,
        params,
        used_params: 2,
    }
}

fn fit_rational_plus_delta1log(points: &[Point]) -> FitSummary {
    // Weighted linear fit y-A_SHIFTED = d1*x through known intercept.
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for p in points {
        let x = ln_ratio(p.q_gev);
        let y = p.value - A_SHIFTED_RATIONAL;
        let w = 1.0 / (p.sigma * p.sigma);
        sxx += w * x * x;
        sxy += w * x * y;
    }
    let d1 = if sxx > 0.0 { sxy / sxx } else { 0.0 };
    let params = [0.0, d1];
    let chi2 = chi2_for("rational_plus_delta1log", params, points);
    let dof = (points.len() as f64 - 1.0).max(1.0);
    FitSummary {
        chi2,
        dof,
        red_chi2: chi2 / dof,
        params,
        used_params: 1,
    }
}

fn fit_fixed(model: &str, points: &[Point]) -> FitSummary {
    let params = [0.0, 0.0];
    let chi2 = chi2_for(model, params, points);
    let dof = (points.len() as f64).max(1.0);
    FitSummary {
        chi2,
        dof,
        red_chi2: chi2 / dof,
        params,
        used_params: 0,
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_WEAK_ANGLE_FIT_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_weak_angle_identity_fit".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Public anchors used in this coarse multi-point fit:
    // 1) PDG MS-bar angle at MZ (table value used in repo context)
    // 2) LHCb 2024 effective leptonic weak mixing angle at Z-scale
    // 3) SLAC E158 low-Q measurement
    let points = vec![
        Point {
            name: "PDG_MSbar_at_MZ",
            q_gev: MZ_GEV,
            value: 0.23122,
            sigma: 0.00006,
            source: "https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf",
            note: "table value used in previous lanes",
        },
        Point {
            name: "LHCb_2024_sin2thetaeff_l",
            q_gev: MZ_GEV,
            value: 0.23147,
            sigma: (0.00044_f64.powi(2) + 0.00005_f64.powi(2) + 0.00023_f64.powi(2)).sqrt(),
            source: "https://doi.org/10.1007/JHEP12(2024)026",
            note: "combined stat+syst+model in quadrature",
        },
        Point {
            name: "CMS_2025_sin2thetaeff_l",
            q_gev: MZ_GEV,
            value: 0.23152,
            sigma: 0.00031,
            source: "https://doi.org/10.1016/j.physletb.2025.139526",
            note: "published hadron-collider measurement",
        },
        Point {
            name: "Tevatron_RunII_combination",
            q_gev: MZ_GEV,
            value: 0.23148,
            sigma: 0.00033,
            source: "https://doi.org/10.1103/PhysRevD.97.112007",
            note: "CDF+D0 combined result",
        },
        Point {
            name: "ATLAS_2015_sin2thetaeff_l",
            q_gev: MZ_GEV,
            value: 0.2308,
            sigma: (0.0005_f64.powi(2) + 0.0006_f64.powi(2) + 0.0009_f64.powi(2)).sqrt(),
            source: "https://doi.org/10.1007/JHEP09(2015)049",
            note: "combined stat+syst+PDF in quadrature",
        },
        Point {
            name: "LEP_SLD_combined_2005",
            q_gev: MZ_GEV,
            value: 0.23153,
            sigma: 0.00016,
            source: "https://doi.org/10.1016/j.physrep.2005.12.006",
            note: "final LEP/SLD combined effective leptonic weak mixing angle",
        },
        Point {
            name: "CDF_RunII_muon_plus_electron_combined",
            q_gev: MZ_GEV,
            value: 0.23221,
            sigma: 0.00046,
            source: "https://doi.org/10.1103/PhysRevD.93.112016",
            note: "CDF combination quoted in PRD 93, 112016",
        },
        Point {
            name: "D0_RunII_electron",
            q_gev: MZ_GEV,
            value: 0.23147,
            sigma: 0.00047,
            source: "https://doi.org/10.1103/PhysRevLett.115.041801",
            note: "D0 extraction from A_FB near Z pole",
        },
        Point {
            name: "SLAC_E158_lowQ",
            q_gev: 0.026_f64.sqrt(),
            value: 0.2397,
            sigma: (0.0010_f64.powi(2) + 0.0008_f64.powi(2)).sqrt(),
            source: "https://www.slac.stanford.edu/exp/e158/plots/results.html",
            note: "Run I-III final result",
        },
    ];

    let fit_base = fit_fixed("base_fixed", &points);
    let fit_rat = fit_fixed("rational_fixed", &points);
    let fit_shift = fit_base_plus_delta0(&points);
    let fit_shift_run = fit_base_plus_delta0_delta1log(&points);
    let fit_rat_run = fit_rational_plus_delta1log(&points);
    let fit_rat_alpha_fixed = fit_fixed("rational_plus_alpha_ln10_over_4pi_log_fixed", &points);

    let models = vec![
        ("base_fixed", fit_base),
        ("rational_fixed", fit_rat),
        ("base_plus_delta0", fit_shift),
        ("base_plus_delta0_delta1log", fit_shift_run),
        ("rational_plus_delta1log", fit_rat_run),
        ("rational_plus_alpha_ln10_over_4pi_log_fixed", fit_rat_alpha_fixed),
    ];

    let mut model_rows = Vec::new();
    for (name, f) in &models {
        let mut pulls = Vec::new();
        for p in &points {
            let pred = eval_model(name, p, f.params);
            let pull = (pred - p.value) / p.sigma;
            pulls.push(json!({
              "point": p.name,
              "pred": pred,
              "obs": p.value,
              "sigma": p.sigma,
              "pull_sigma": pull
            }));
        }
        model_rows.push(json!({
          "model": name,
          "chi2": f.chi2,
          "dof": f.dof,
          "reduced_chi2": f.red_chi2,
          "used_params": f.used_params,
          "params": {
            "delta0": f.params[0],
            "delta1_log": f.params[1]
          },
          "pulls": pulls
        }));
    }

    let best = models
        .iter()
        .min_by(|a, b| a.1.red_chi2.total_cmp(&b.1.red_chi2))
        .map(|(n, f)| json!({"model": n, "reduced_chi2": f.red_chi2}))
        .unwrap_or_else(|| json!({"model":"none","reduced_chi2": f64::NAN}));

    let payload = json!({
      "scope": "coarse multi-point weak-angle identity fit",
      "assumptions": {
        "base_identity": "3/13",
        "rational_candidate": "3/13 + 1/13^3",
        "running_basis": "ln(MZ/Q)",
        "alpha_inverse": ALPHA_INV_PDG_2025,
        "delta1_fixed_candidate": "alpha * ln(10) / (4*pi)"
      },
      "points": points.iter().map(|p| json!({
        "name": p.name,
        "Q_GeV": p.q_gev,
        "value": p.value,
        "sigma": p.sigma,
        "source": p.source,
        "note": p.note
      })).collect::<Vec<_>>(),
      "model_fits": model_rows,
      "best_by_reduced_chi2": best
    });

    let txt_path = out.join("ctc_weak_angle_identity_fit.txt");
    let json_path = out.join("ctc_weak_angle_identity_fit.json");

    let mut txt = String::new();
    txt.push_str("[ctc_weak_angle_identity_fit]\n");
    txt.push_str("coarse multi-point weak-angle identity fit\n\n");
    txt.push_str(&format!("A_base = {:.12}\n", A_BASE));
    txt.push_str(&format!("A_rational = {:.12}\n", A_SHIFTED_RATIONAL));
    txt.push_str(&format!("points = {}\n", points.len()));
    for (name, f) in &models {
        txt.push_str(&format!(
            "{}: chi2={:.6}, dof={:.1}, red_chi2={:.6}, delta0={:.9e}, delta1={:.9e}\n",
            name, f.chi2, f.dof, f.red_chi2, f.params[0], f.params[1]
        ));
    }

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
