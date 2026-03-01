//! Execute next-three falsification tests:
//! - T07: DESI void scalar construction (coarse, literature-null extraction)
//! - T11: weak-angle chronological holdout
//! - T13: weak-angle scheme-coherent fit
//!
//! Scope:
//! - Fast public-data falsification pass.
//! - Coarse assumptions are explicit in output.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

const MZ_GEV: f64 = 91.1876;
const ALPHA_INV: f64 = 137.035_999_177;
const A_RATIONAL: f64 = 3.0 / 13.0 + 1.0 / (13.0 * 13.0 * 13.0);

#[derive(Clone, Copy)]
enum Scheme {
    Msbar,
    Effective,
}

#[derive(Clone, Copy)]
struct EwPoint {
    name: &'static str,
    year: i32,
    q_gev: f64,
    value: f64,
    sigma: f64,
    scheme: Scheme,
}

fn pred_zero_free(q_gev: f64) -> f64 {
    let alpha = 1.0 / ALPHA_INV;
    let k = alpha * 10.0_f64.ln() / (4.0 * std::f64::consts::PI);
    A_RATIONAL + k * (MZ_GEV / q_gev).ln()
}

fn reduced_chi2(points: &[EwPoint], values: &[f64], sigmas: &[f64]) -> f64 {
    let mut chi2 = 0.0_f64;
    for i in 0..points.len() {
        let pred = pred_zero_free(points[i].q_gev);
        let z = (pred - values[i]) / sigmas[i];
        chi2 += z * z;
    }
    let dof = (points.len() as f64).max(1.0);
    chi2 / dof
}

fn main() {
    let out_dir = std::env::var("GUTOE_NEXT3_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_next3_t07_t11_t13".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Same anchors as the expanded weak fit lane.
    let points = vec![
        EwPoint {
            name: "PDG_MSbar_at_MZ",
            year: 2025,
            q_gev: MZ_GEV,
            value: 0.23122,
            sigma: 0.00006,
            scheme: Scheme::Msbar,
        },
        EwPoint {
            name: "LHCb_2024_sin2thetaeff_l",
            year: 2024,
            q_gev: MZ_GEV,
            value: 0.23147,
            sigma: (0.00044_f64.powi(2) + 0.00005_f64.powi(2) + 0.00023_f64.powi(2)).sqrt(),
            scheme: Scheme::Effective,
        },
        EwPoint {
            name: "CMS_2025_sin2thetaeff_l",
            year: 2025,
            q_gev: MZ_GEV,
            value: 0.23152,
            sigma: 0.00031,
            scheme: Scheme::Effective,
        },
        EwPoint {
            name: "Tevatron_RunII_combination",
            year: 2018,
            q_gev: MZ_GEV,
            value: 0.23148,
            sigma: 0.00033,
            scheme: Scheme::Effective,
        },
        EwPoint {
            name: "ATLAS_2015_sin2thetaeff_l",
            year: 2015,
            q_gev: MZ_GEV,
            value: 0.2308,
            sigma: (0.0005_f64.powi(2) + 0.0006_f64.powi(2) + 0.0009_f64.powi(2)).sqrt(),
            scheme: Scheme::Effective,
        },
        EwPoint {
            name: "LEP_SLD_combined_2005",
            year: 2005,
            q_gev: MZ_GEV,
            value: 0.23153,
            sigma: 0.00016,
            scheme: Scheme::Effective,
        },
        EwPoint {
            name: "CDF_RunII_muon_plus_electron_combined",
            year: 2016,
            q_gev: MZ_GEV,
            value: 0.23221,
            sigma: 0.00046,
            scheme: Scheme::Effective,
        },
        EwPoint {
            name: "D0_RunII_electron",
            year: 2015,
            q_gev: MZ_GEV,
            value: 0.23147,
            sigma: 0.00047,
            scheme: Scheme::Effective,
        },
        EwPoint {
            name: "SLAC_E158_lowQ",
            year: 2005,
            q_gev: 0.026_f64.sqrt(),
            value: 0.2397,
            sigma: (0.0010_f64.powi(2) + 0.0008_f64.powi(2)).sqrt(),
            scheme: Scheme::Msbar,
        },
    ];

    // T13: scheme-coherent conversion of effective-angle points to MSbar-like lane.
    // Conversion anchor from PDG Z-pole differences:
    //   s_eff - s_MS ~ 0.23154 - 0.23122 = 0.00032.
    // Coarse conversion uncertainty chosen as 0.00009.
    let delta_eff_to_ms = 0.00032_f64;
    let sigma_conv = 0.00009_f64;
    let mut values_sc = Vec::new();
    let mut sigmas_sc = Vec::new();
    for p in &points {
        match p.scheme {
            Scheme::Msbar => {
                values_sc.push(p.value);
                sigmas_sc.push(p.sigma);
            }
            Scheme::Effective => {
                values_sc.push(p.value - delta_eff_to_ms);
                sigmas_sc.push((p.sigma * p.sigma + sigma_conv * sigma_conv).sqrt());
            }
        }
    }
    let red_t13 = reduced_chi2(&points, &values_sc, &sigmas_sc);
    let t13_status = if red_t13 <= 1.5 {
        "PASS"
    } else if red_t13 > 3.0 {
        "FAIL"
    } else {
        "OPEN"
    };

    // T11: chronological holdout (post-2016 as holdout), no refit.
    let holdout: Vec<EwPoint> = points.iter().copied().filter(|p| p.year >= 2017).collect();
    let holdout_values: Vec<f64> = holdout.iter().map(|p| p.value).collect();
    let holdout_sigmas: Vec<f64> = holdout.iter().map(|p| p.sigma).collect();
    let red_t11 = reduced_chi2(&holdout, &holdout_values, &holdout_sigmas);
    let t11_status = if red_t11 <= 1.5 {
        "PASS"
    } else if red_t11 > 3.0 {
        "FAIL"
    } else {
        "OPEN"
    };

    // T07: DESI void scalar construction (coarse literature-null extraction).
    // Current published wording: "generally consistent void properties".
    // Convert this to scalar anomaly z=0 for kill/pass gate tracking.
    let z_void_scalar = 0.0_f64;
    let t07_status = if z_void_scalar >= 5.0 {
        "PASS"
    } else if z_void_scalar < 3.0 {
        "FAIL"
    } else {
        "OPEN"
    };

    let payload = json!({
      "scope": "next-three tests execution (T07, T11, T13)",
      "assumptions": {
        "t13_eff_to_ms_delta": delta_eff_to_ms,
        "t13_conversion_sigma": sigma_conv,
        "t11_holdout_year_threshold": 2017,
        "t07_void_scalar_method": "coarse literature-null extraction from DESIVAST consistency statement"
      },
      "tests": {
        "T13_scheme_coherent_fit": {
          "metric": "reduced_chi2_scheme_clean",
          "value": red_t13,
          "pass_if": "<=1.5",
          "kill_if": ">3.0",
          "status": t13_status
        },
        "T11_chronological_holdout": {
          "metric": "holdout_reduced_chi2",
          "value": red_t11,
          "pass_if": "<=1.5",
          "kill_if": ">3.0",
          "holdout_points": holdout.iter().map(|p| p.name).collect::<Vec<_>>(),
          "status": t11_status
        },
        "T07_void_scalar_lane": {
          "metric": "z_void_scalar",
          "value": z_void_scalar,
          "pass_if": ">=5.0",
          "kill_if": "<3.0",
          "status": t07_status,
          "source_note": "DESIVAST/adb559 abstract wording: generally consistent void properties"
        }
      }
    });

    let txt_path = out.join("ctc_next3_t07_t11_t13.txt");
    let json_path = out.join("ctc_next3_t07_t11_t13.json");

    let mut txt = String::new();
    txt.push_str("[ctc_next3_t07_t11_t13]\n");
    txt.push_str("next-three falsification tests\n\n");
    txt.push_str(&format!(
        "T13 reduced_chi2_scheme_clean = {:.6} status={}\n",
        red_t13, t13_status
    ));
    txt.push_str(&format!(
        "T11 holdout_reduced_chi2 = {:.6} status={}\n",
        red_t11, t11_status
    ));
    txt.push_str(&format!(
        "T07 z_void_scalar = {:.6} status={}\n",
        z_void_scalar, t07_status
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
