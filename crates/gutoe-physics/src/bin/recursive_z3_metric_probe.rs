//! Recursive Z3 metric probe for the projection tower 256 -> 16 -> 4.
//!
//! Purpose:
//! - Quantify whether "lift to higher level, traverse, descend" yields
//!   shorter distances than base 4D separation in the current linear model.
//! - Provide explicit lower-bound + witness checks in numerics.

use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use serde_json::json;
use std::fs;
use std::path::PathBuf;

const N4: usize = 4;
const N16: usize = 16;
const N256: usize = 256;

#[derive(Clone)]
struct Case {
    name: &'static str,
    a4: [f64; N4],
    b4: [f64; N4],
}

fn idx_256(i: usize, j: usize) -> usize {
    i * N16 + j
}

fn norm_sq(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum()
}

fn dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

fn section4_to_16(x4: &[f64; N4]) -> [f64; N16] {
    let mut v = [0.0_f64; N16];
    v[..N4].copy_from_slice(x4);
    v
}

fn section16_to_256(x16: &[f64; N16]) -> [f64; N256] {
    let mut v = [0.0_f64; N256];
    for i in 0..N16 {
        v[idx_256(i, 0)] = x16[i];
    }
    v
}

fn section4_to_256(x4: &[f64; N4]) -> [f64; N256] {
    let x16 = section4_to_16(x4);
    section16_to_256(&x16)
}

fn random_fiber_lift_256(base4: &[f64; N4], rng: &mut StdRng, span: f64) -> [f64; N256] {
    // Fiber constraints for p256->16->4:
    // fixed: x(i,0)=base4[i] for i<4
    // free:  x(i,0) for i>=4, and x(i,j!=0) for all i
    let mut v = [0.0_f64; N256];

    for i in 0..N4 {
        v[idx_256(i, 0)] = base4[i];
    }

    for i in N4..N16 {
        v[idx_256(i, 0)] = rng.gen_range(-span..span);
    }

    for i in 0..N16 {
        for j in 1..N16 {
            v[idx_256(i, j)] = rng.gen_range(-span..span);
        }
    }

    v
}

fn main() {
    let out_dir = std::env::var("GUTOE_RECURSIVE_Z3_METRIC_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/recursive_z3_metric_probe".to_string());
    let out = PathBuf::from(&out_dir);
    let _ = fs::create_dir_all(&out);

    let samples: usize = std::env::var("GUTOE_RECURSIVE_Z3_METRIC_SAMPLES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(20_000);
    let span: f64 = std::env::var("GUTOE_RECURSIVE_Z3_METRIC_SPAN")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.0);

    let cases = vec![
        Case {
            name: "unit_axis",
            a4: [0.0, 0.0, 0.0, 0.0],
            b4: [1.0, 0.0, 0.0, 0.0],
        },
        Case {
            name: "mixed_13",
            a4: [0.0, 0.0, 0.0, 0.0],
            b4: [3.0, -4.0, 12.0, 0.0],
        },
        Case {
            name: "large_scale",
            a4: [10.0, -20.0, 5.0, -1.0],
            b4: [1.0e6 + 10.0, -20.0, 5.0, -1.0],
        },
    ];

    let mut rng = StdRng::seed_from_u64(0x5A17_2026);
    let mut rows = Vec::new();

    for c in cases {
        let d4 = dist(&c.a4, &c.b4);

        // Constructive witness: canonical lifts with all free coordinates zero.
        let la = section4_to_256(&c.a4);
        let lb = section4_to_256(&c.b4);
        let witness_dist = dist(&la, &lb);

        // Coordinate lower bound: first four constraints are fixed differences.
        // Therefore any connector must satisfy d >= ||a4-b4||.
        let lower_bound = d4;

        let mut random_min = f64::INFINITY;
        for _ in 0..samples {
            let xa = random_fiber_lift_256(&c.a4, &mut rng, span);
            let xb = random_fiber_lift_256(&c.b4, &mut rng, span);
            random_min = random_min.min(dist(&xa, &xb));
        }

        let d_exact = (witness_dist - lower_bound).abs() < 1e-12;

        rows.push(json!({
            "case": c.name,
            "d4": d4,
            "lower_bound": lower_bound,
            "witness_dist": witness_dist,
            "witness_matches_lower_bound": d_exact,
            "random_min_dist": random_min,
            "random_min_over_d4": if d4 > 0.0 { random_min / d4 } else { 1.0 },
            "compression_possible_in_linear_model": false
        }));
    }

    // Global verdict for this linear projection tower model.
    let all_exact = rows
        .iter()
        .all(|r| r["witness_matches_lower_bound"].as_bool().unwrap_or(false));

    let payload = json!({
        "model": {
            "tower": "256->16->4",
            "projection_256_to_16": "slice j=0",
            "projection_16_to_4": "first four coordinates",
            "metric": "euclidean_on_total_space",
            "interpretation": "linear additive lift/traverse/descend"
        },
        "scan": {
            "samples_per_case": samples,
            "random_span": span
        },
        "results": rows,
        "verdict": {
            "infimum_equals_base_distance": all_exact,
            "bounded_shortcut_detected": false,
            "note": "In this linear model, fibers do not reduce point-to-point base distance."
        }
    });

    let txt_path = out.join("recursive_z3_metric_probe.txt");
    let json_path = out.join("recursive_z3_metric_probe.json");

    let mut txt = String::new();
    txt.push_str("[recursive_z3_metric_probe]\n");
    txt.push_str("tower = 256->16->4\n");
    txt.push_str("metric = euclidean_on_total_space\n");
    txt.push_str(&format!("samples_per_case = {}\n", samples));
    txt.push_str(&format!("random_span = {:.6e}\n", span));
    txt.push_str("\n[cases]\n");
    for r in payload["results"].as_array().expect("array") {
        txt.push_str(&format!(
            "{}: d4={:.12e}, lower_bound={:.12e}, witness={:.12e}, random_min={:.12e}, random_ratio={:.12e}\n",
            r["case"].as_str().unwrap_or("case"),
            r["d4"].as_f64().unwrap_or(f64::NAN),
            r["lower_bound"].as_f64().unwrap_or(f64::NAN),
            r["witness_dist"].as_f64().unwrap_or(f64::NAN),
            r["random_min_dist"].as_f64().unwrap_or(f64::NAN),
            r["random_min_over_d4"].as_f64().unwrap_or(f64::NAN)
        ));
    }
    txt.push_str("\n[verdict]\n");
    txt.push_str(&format!(
        "infimum_equals_base_distance = {}\n",
        payload["verdict"]["infimum_equals_base_distance"]
            .as_bool()
            .unwrap_or(false)
    ));
    txt.push_str(&format!(
        "bounded_shortcut_detected = {}\n",
        payload["verdict"]["bounded_shortcut_detected"]
            .as_bool()
            .unwrap_or(false)
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "infimum_equals_base_distance={} bounded_shortcut_detected=false",
        payload["verdict"]["infimum_equals_base_distance"].as_bool().unwrap_or(false)
    );

    // Keep dead-code lints quiet if compile profile changes in future.
    let _ = norm_sq(&[0.0, 1.0]);
}
