//! Dual-lane search requested by user:
//! Lane A: exact two-term coefficient extraction for alpha^-1 = 137 + 5a - b a^2
//! Lane B: shared higher-order structure across alpha^-1 and mp/me with grade-tied scaling.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

const ALPHA_INV_PHYS: f64 = 137.035_999_177;
const ALPHA: f64 = 1.0 / ALPHA_INV_PHYS;
const PI: f64 = std::f64::consts::PI;
const MP_ME_EXP: f64 = 1836.152_673_43;

const LO_ALPHA_INV: f64 = 137.0;
const A1: f64 = 5.0; // grade levels

fn ppm(pred: f64, truth: f64) -> f64 {
    (pred - truth).abs() / truth * 1.0e6
}

fn main() {
    let out_dir = std::env::var("GUTOE_ALPHA_DUAL_LANE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_alpha_dual_lane_search".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let lo_mp = 6.0 * PI.powi(5);

    // Residuals after linear 5a term.
    let d_alpha = ALPHA_INV_PHYS - (LO_ALPHA_INV + A1 * ALPHA);
    let d_mp = MP_ME_EXP - (lo_mp + A1 * ALPHA);

    // Exact two-term coefficients (cubic term absent).
    let b_alpha_exact = -d_alpha / (ALPHA * ALPHA);
    let b_mp_exact = -d_mp / (ALPHA * ALPHA);

    // Lane A: compact structural approximants near b_alpha_exact.
    let mut approximants = Vec::new();
    let denoms = [
        7_i32, 12, 13, 16, 17, 24, 27, 32, 36, 48, 64, 81, 96, 108, 128, 144,
    ];

    for &q in &denoms {
        for p in 1_i32..=32_i32 {
            let b = 9.0 + (p as f64) / (q as f64);
            approximants.push((format!("9 + {}/{}", p, q), b));
        }
    }

    approximants.sort_by(|a, b| {
        (a.1 - b_alpha_exact)
            .abs()
            .total_cmp(&(b.1 - b_alpha_exact).abs())
    });
    let top_a = approximants.into_iter().take(8).collect::<Vec<_>>();

    // Lane B: shared higher-order structure
    // alpha^-1 = 137 + 5a - b a^2 + c a^3
    // mp/me    = 6pi^5 + 5a - B a^2 + C a^3
    // with B=g b, C=g c (grade-tied scale factor g).

    // If both equations are exact under this tie, g is fixed by residual ratio.
    let g_exact = d_mp / d_alpha;

    // For any chosen b, c is fixed by alpha equation.
    let solve_c = |b: f64| -> f64 { (d_alpha + b * ALPHA * ALPHA) / (ALPHA * ALPHA * ALPHA) };

    let scenarios = vec![
        ("b=9, g=4", 9.0_f64, 4.0_f64),
        ("b=9+5/32, g=4", 9.0 + 5.0 / 32.0, 4.0_f64),
        ("b=b_alpha_exact, g=4", b_alpha_exact, 4.0_f64),
        ("b=9, g=g_exact", 9.0_f64, g_exact),
        ("b=9+5/32, g=g_exact", 9.0 + 5.0 / 32.0, g_exact),
    ];

    let mut scenario_rows = Vec::new();
    for (name, b, g) in scenarios {
        let c = solve_c(b);
        let b_mp = g * b;
        let c_mp = g * c;

        let alpha_pred = LO_ALPHA_INV + A1 * ALPHA - b * ALPHA * ALPHA + c * ALPHA * ALPHA * ALPHA;
        let mp_pred = lo_mp + A1 * ALPHA - b_mp * ALPHA * ALPHA + c_mp * ALPHA * ALPHA * ALPHA;

        scenario_rows.push(json!({
          "name": name,
          "b": b,
          "c": c,
          "g": g,
          "B": b_mp,
          "C": c_mp,
          "alpha_pred": alpha_pred,
          "alpha_ppm_error": ppm(alpha_pred, ALPHA_INV_PHYS),
          "mp_pred": mp_pred,
          "mp_ppm_error": ppm(mp_pred, MP_ME_EXP)
        }));
    }

    let payload = json!({
      "scope": "dual-lane coefficient search",
      "laneA_exact_two_term": {
        "equation": "alpha^-1 = 137 + 5a - b a^2",
        "b_alpha_exact": b_alpha_exact,
        "equation_mp": "mp/me = 6pi^5 + 5a - B a^2",
        "B_mp_exact": b_mp_exact,
        "ratio_B_over_b": b_mp_exact / b_alpha_exact,
        "top_structural_approximants_for_b": top_a.iter().map(|(name, b)| {
          json!({
            "expr": name,
            "value": b,
            "abs_error": (b - b_alpha_exact).abs()
          })
        }).collect::<Vec<_>>()
      },
      "laneB_shared_higher_order": {
        "equations": [
          "alpha^-1 = 137 + 5a - b a^2 + c a^3",
          "mp/me    = 6pi^5 + 5a - B a^2 + C a^3"
        ],
        "tie": "B=g*b, C=g*c",
        "g_exact_from_residual_ratio": g_exact,
        "scenarios": scenario_rows
      }
    });

    let txt_path = out.join("ctc_alpha_dual_lane_search.txt");
    let json_path = out.join("ctc_alpha_dual_lane_search.json");

    let mut txt = String::new();
    txt.push_str("[ctc_alpha_dual_lane_search]\n");
    txt.push_str("Lane A + Lane B requested dual-lane run\n\n");

    txt.push_str("[Lane A exact two-term]\n");
    txt.push_str(&format!("b_alpha_exact = {:.12}\n", b_alpha_exact));
    txt.push_str(&format!("B_mp_exact    = {:.12}\n", b_mp_exact));
    txt.push_str(&format!("B/b ratio      = {:.12}\n", b_mp_exact / b_alpha_exact));
    txt.push_str("top approximants near b_alpha:\n");
    for (name, b) in top_a {
        txt.push_str(&format!(
            "  {} = {:.12}  abs_err={:.6e}\n",
            name,
            b,
            (b - b_alpha_exact).abs()
        ));
    }

    txt.push_str("\n[Lane B shared higher-order]\n");
    txt.push_str(&format!("g_exact = {:.12}\n", g_exact));
    txt.push_str("scenarios (alpha ppm, mp ppm):\n");
    if let Some(arr) = payload
        .get("laneB_shared_higher_order")
        .and_then(|v| v.get("scenarios"))
        .and_then(|v| v.as_array())
    {
        for row in arr {
            let name = row.get("name").and_then(|v| v.as_str()).unwrap_or("?");
            let b = row.get("b").and_then(|v| v.as_f64()).unwrap_or(f64::NAN);
            let c = row.get("c").and_then(|v| v.as_f64()).unwrap_or(f64::NAN);
            let g = row.get("g").and_then(|v| v.as_f64()).unwrap_or(f64::NAN);
            let alpha_ppm = row
                .get("alpha_ppm_error")
                .and_then(|v| v.as_f64())
                .unwrap_or(f64::NAN);
            let mp_ppm = row
                .get("mp_ppm_error")
                .and_then(|v| v.as_f64())
                .unwrap_or(f64::NAN);
            txt.push_str(&format!(
                "  {}: b={:.9}, c={:.9}, g={:.9} -> alpha_ppm={:.6}, mp_ppm={:.6}\n",
                name, b, c, g, alpha_ppm, mp_ppm
            ));
        }
    }

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
