//! Topological defect bundle probe.
//!
//! Measures distance reduction when a compact defect creates a bridge edge in a
//! 1D projected manifold model.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;
use std::fs;
use std::path::PathBuf;

fn abs(x: f64) -> f64 {
    x.abs()
}

fn base_distance(a: f64, b: f64) -> f64 {
    abs(b - a)
}

fn defect_distance(a: f64, b: f64, l: f64, r: f64) -> f64 {
    let direct = base_distance(a, b);
    let via_lr = abs(a - l) + abs(b - r);
    let via_rl = abs(a - r) + abs(b - l);
    direct.min(via_lr.min(via_rl))
}

fn main() {
    let out_dir = std::env::var("GUTOE_DEFECT_BUNDLE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/topological_defect_bundle_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let samples: usize = std::env::var("GUTOE_DEFECT_BUNDLE_SAMPLES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(100_000);

    // Deterministic canonical witness.
    let witness_a = -1000.0;
    let witness_b = 1500.0;
    let witness_l = -10.0;
    let witness_r = 10.0;
    let witness_base = base_distance(witness_a, witness_b);
    let witness_defect = defect_distance(witness_a, witness_b, witness_l, witness_r);

    // Random scan over endpoints and local bridge inside compact support [-R, R].
    let r_support = 50.0;
    let mut rng = StdRng::seed_from_u64(0xD3F3_2026);

    let mut best_ratio = f64::INFINITY;
    let mut best = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut improved_count = 0usize;

    for _ in 0..samples {
        let a = rng.gen_range(-2000.0..2000.0);
        let b = rng.gen_range(-2000.0..2000.0);

        let l = rng.gen_range(-r_support..r_support);
        let r = rng.gen_range(-r_support..r_support);

        let d0 = base_distance(a, b);
        if d0 <= 1e-12 {
            continue;
        }
        let d1 = defect_distance(a, b, l, r);

        if d1 + 1e-12 < d0 {
            improved_count += 1;
        }

        let ratio = d1 / d0;
        if ratio < best_ratio {
            best_ratio = ratio;
            best = (a, b, l, r, d0, d1);
        }
    }

    let payload = json!({
        "model": {
            "base_distance": "|b-a|",
            "defect_distance": "min(direct, via bridge l<->r)",
            "compact_support": [-r_support, r_support]
        },
        "witness": {
            "a": witness_a,
            "b": witness_b,
            "l": witness_l,
            "r": witness_r,
            "base": witness_base,
            "defect": witness_defect,
            "ratio": witness_defect / witness_base
        },
        "random_scan": {
            "samples": samples,
            "improved_count": improved_count,
            "improved_fraction": improved_count as f64 / samples as f64,
            "best_ratio": best_ratio,
            "best_case": {
                "a": best.0,
                "b": best.1,
                "l": best.2,
                "r": best.3,
                "base": best.4,
                "defect": best.5
            }
        }
    });

    let txt_path = out.join("topological_defect_bundle_probe.txt");
    let json_path = out.join("topological_defect_bundle_probe.json");

    let mut txt = String::new();
    txt.push_str("[topological_defect_bundle_probe]\n");
    txt.push_str(&format!(
        "witness: base={:.12e}, defect={:.12e}, ratio={:.12e}\n",
        witness_base,
        witness_defect,
        witness_defect / witness_base
    ));
    txt.push_str(&format!(
        "random: samples={}, improved_fraction={:.12e}, best_ratio={:.12e}\n",
        samples,
        improved_count as f64 / samples as f64,
        best_ratio
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
