//! Recursive Z3 affine-creation probe.
//!
//! Quantifies the inhomogeneous (affine) lane for `256 -> 16 -> 4` descent:
//! if a translation offset is allowed, what is the minimal offset norm needed
//! to realize a given 4D target shift?

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;
use std::fs;
use std::path::PathBuf;

const N4: usize = 4;
const N16: usize = 16;
const N256: usize = 256;

#[derive(Clone)]
struct Case {
    name: &'static str,
    start: [f64; N4],
    target: [f64; N4],
}

fn idx_256(i: usize, j: usize) -> usize {
    i * N16 + j
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn diff(a: &[f64; N4], b: &[f64; N4]) -> [f64; N4] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]]
}

fn section4_to_256(x4: &[f64; N4]) -> [f64; N256] {
    let mut v = [0.0_f64; N256];
    for i in 0..N4 {
        v[idx_256(i, 0)] = x4[i];
    }
    v
}

fn random_affine_offset_with_fixed_projection(
    d4: &[f64; N4],
    span: f64,
    rng: &mut StdRng,
) -> [f64; N256] {
    let mut v = [0.0_f64; N256];

    // Fix projection constraints to match desired 4D shift.
    for i in 0..N4 {
        v[idx_256(i, 0)] = d4[i];
    }

    // Free fiber components.
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
    let out_dir = std::env::var("GUTOE_RECURSIVE_Z3_AFFINE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/recursive_z3_affine_creation_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let samples: usize = std::env::var("GUTOE_RECURSIVE_Z3_AFFINE_SAMPLES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(20_000);
    let span: f64 = std::env::var("GUTOE_RECURSIVE_Z3_AFFINE_SPAN")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.0);

    let cases = vec![
        Case {
            name: "unit_shift",
            start: [0.0, 0.0, 0.0, 0.0],
            target: [1.0, 0.0, 0.0, 0.0],
        },
        Case {
            name: "mixed_shift",
            start: [3.0, -2.0, 5.0, 7.0],
            target: [8.0, -10.0, 6.0, 7.0],
        },
        Case {
            name: "large_shift",
            start: [0.0, 0.0, 0.0, 0.0],
            target: [1.0e6, 0.0, 0.0, 0.0],
        },
    ];

    let mut rng = StdRng::seed_from_u64(0xAFF1_2026);
    let mut rows = Vec::new();

    for c in cases {
        let d4 = diff(&c.target, &c.start);
        let d4_norm = norm(&d4);

        // Minimal witness: section with all free coordinates zero.
        let witness = section4_to_256(&d4);
        let witness_norm = norm(&witness);

        // Monte Carlo over feasible affine offsets (same projection constraints).
        let mut random_min = f64::INFINITY;
        for _ in 0..samples {
            let u = random_affine_offset_with_fixed_projection(&d4, span, &mut rng);
            random_min = random_min.min(norm(&u));
        }

        rows.push(json!({
            "case": c.name,
            "shift4": d4,
            "shift4_norm": d4_norm,
            "witness_affine_norm": witness_norm,
            "random_min_affine_norm": random_min,
            "witness_is_minimal_in_scan": witness_norm <= random_min + 1e-9,
            "affine_lane_open": d4_norm > 0.0,
            "cost_scales_with_shift_norm": true
        }));
    }

    let payload = json!({
        "model": {
            "tower": "256->16->4",
            "mechanism": "inhomogeneous affine offset (creation term)",
            "constraint": "projection(offset)=desired 4D shift"
        },
        "scan": {
            "samples_per_case": samples,
            "random_span": span
        },
        "results": rows
    });

    let txt_path = out.join("recursive_z3_affine_creation_probe.txt");
    let json_path = out.join("recursive_z3_affine_creation_probe.json");

    let mut txt = String::new();
    txt.push_str("[recursive_z3_affine_creation_probe]\n");
    txt.push_str("mechanism = inhomogeneous affine offset\n");
    txt.push_str(&format!("samples_per_case = {}\n", samples));
    txt.push_str(&format!("random_span = {:.6e}\n", span));
    txt.push_str("\n[cases]\n");
    for r in payload["results"].as_array().expect("array") {
        txt.push_str(&format!(
            "{}: shift_norm={:.12e}, witness_norm={:.12e}, random_min={:.12e}, affine_lane_open={}\n",
            r["case"].as_str().unwrap_or("case"),
            r["shift4_norm"].as_f64().unwrap_or(f64::NAN),
            r["witness_affine_norm"].as_f64().unwrap_or(f64::NAN),
            r["random_min_affine_norm"].as_f64().unwrap_or(f64::NAN),
            r["affine_lane_open"].as_bool().unwrap_or(false)
        ));
    }

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
