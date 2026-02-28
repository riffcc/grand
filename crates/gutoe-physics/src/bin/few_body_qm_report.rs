//! Few-body QM light-nucleus report (A<=16 focus).

use gutoe_physics::{scan_nuclear_chart, ScanConfig};
use serde_json::json;
use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

#[derive(Clone, Copy)]
struct Spot {
    z: u16,
    n: u16,
    a: u16,
    b_per_a_obs: f64,
}

// AME2020 spot checks spanning few-body and near-transition region.
const SPOTS: &[Spot] = &[
    Spot {
        z: 1,
        n: 1,
        a: 2,
        b_per_a_obs: 1.112,
    }, // H-2
    Spot {
        z: 1,
        n: 2,
        a: 3,
        b_per_a_obs: 2.827,
    }, // H-3
    Spot {
        z: 2,
        n: 1,
        a: 3,
        b_per_a_obs: 2.573,
    }, // He-3
    Spot {
        z: 2,
        n: 2,
        a: 4,
        b_per_a_obs: 7.074,
    }, // He-4
    Spot {
        z: 3,
        n: 4,
        a: 7,
        b_per_a_obs: 5.606,
    }, // Li-7
    Spot {
        z: 4,
        n: 5,
        a: 9,
        b_per_a_obs: 6.463,
    }, // Be-9
    Spot {
        z: 5,
        n: 6,
        a: 11,
        b_per_a_obs: 6.928,
    }, // B-11
    Spot {
        z: 6,
        n: 6,
        a: 12,
        b_per_a_obs: 7.680,
    }, // C-12
    Spot {
        z: 7,
        n: 7,
        a: 14,
        b_per_a_obs: 7.476,
    }, // N-14
    Spot {
        z: 8,
        n: 8,
        a: 16,
        b_per_a_obs: 7.976,
    }, // O-16
];

fn main() {
    let out_dir = std::env::var("GUTOE_FEW_BODY_OUT").unwrap_or_else(|_| "/tmp/nuclear_chart".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let cfg = ScanConfig {
        z_min: 1,
        z_max: 16,
        n_min: 1,
        n_max: 20,
        ..ScanConfig::default()
    };
    let records = scan_nuclear_chart(cfg);
    let mut map = BTreeMap::new();
    for r in &records {
        map.insert((r.z, r.n), r.binding_per_nucleon_mev);
    }

    let mut rows = Vec::new();
    let mut sae = 0.0;
    let mut sse = 0.0;
    let mut matched = 0usize;
    for s in SPOTS {
        if let Some(&pred) = map.get(&(s.z, s.n)) {
            let resid = pred - s.b_per_a_obs;
            sae += resid.abs();
            sse += resid * resid;
            matched += 1;
            rows.push(json!({
                "z": s.z,
                "n": s.n,
                "a": s.a,
                "b_per_a_obs_mev": s.b_per_a_obs,
                "b_per_a_pred_mev": pred,
                "residual_mev_per_a": resid
            }));
        }
    }
    let mae = if matched > 0 { sae / matched as f64 } else { f64::NAN };
    let rmse = if matched > 0 {
        (sse / matched as f64).sqrt()
    } else {
        f64::NAN
    };

    let report = json!({
        "meta": {
            "lane": "few_body_qm_light_nuclei",
            "scan_bounds": {"z_min": 1, "z_max": 16, "n_min": 1, "n_max": 20},
            "spot_rows": SPOTS.len()
        },
        "summary": {
            "matched_rows": matched,
            "mae_mev_per_a": mae,
            "rmse_mev_per_a": rmse
        },
        "spots": rows
    });

    let mut csv = String::from("Z,N,A,b_per_a_obs_mev,b_per_a_pred_mev,residual_mev_per_a\n");
    for row in &rows {
        csv.push_str(&format!(
            "{},{},{},{:.6},{:.6},{:.6}\n",
            row["z"].as_u64().unwrap_or(0),
            row["n"].as_u64().unwrap_or(0),
            row["a"].as_u64().unwrap_or(0),
            row["b_per_a_obs_mev"].as_f64().unwrap_or(0.0),
            row["b_per_a_pred_mev"].as_f64().unwrap_or(0.0),
            row["residual_mev_per_a"].as_f64().unwrap_or(0.0)
        ));
    }

    let txt = format!(
        "[few_body_qm]\nmatched_rows = {}\nmae_mev_per_a = {:.6}\nrmse_mev_per_a = {:.6}\n",
        matched, mae, rmse
    );

    let txt_path = out.join("few_body_qm_report.txt");
    let json_path = out.join("few_body_qm_report.json");
    let csv_path = out.join("few_body_qm_report.csv");
    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&report).expect("json"),
    )
    .expect("write json");
    fs::write(&csv_path, csv).expect("write csv");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("wrote {}", csv_path.display());
    println!(
        "few_body_qm: matched={}, mae={:.4} MeV/A, rmse={:.4} MeV/A",
        matched, mae, rmse
    );
}
