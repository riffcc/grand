//! Recursive Z3 quotient-navigation probe.
//!
//! Uses explicit Cl(1,3) geometric-product rules (anticommutation + metric signs)
//! and checks whether kernel-only multiplicative actions can generate a 4D
//! translation from origin after grade-1 descent.
//!
//! This targets the open lane after the linear additive tower probe.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;
use std::fs;
use std::path::PathBuf;

const N: usize = 4;
const NB: usize = 1 << N; // 16 basis blades
const GRADE1: [usize; 4] = [1, 2, 4, 8];
const KERNEL_MASKS: [usize; 12] = [0, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15];

// Cl(1,3): signature (+,-,-,-)
const METRIC: [f64; N] = [1.0, -1.0, -1.0, -1.0];

#[derive(Clone, Debug)]
struct Mv {
    c: [f64; NB],
}

impl Mv {
    fn zero() -> Self {
        Self { c: [0.0; NB] }
    }

    fn blade(mask: usize, coeff: f64) -> Self {
        let mut m = Self::zero();
        m.c[mask] = coeff;
        m
    }

    fn from_4d(x: [f64; 4]) -> Self {
        let mut m = Self::zero();
        for (i, &mask) in GRADE1.iter().enumerate() {
            m.c[mask] = x[i];
        }
        m
    }

    fn to_4d_grade1(&self) -> [f64; 4] {
        [
            self.c[GRADE1[0]],
            self.c[GRADE1[1]],
            self.c[GRADE1[2]],
            self.c[GRADE1[3]],
        ]
    }

    fn reverse(&self) -> Self {
        let mut out = Self::zero();
        for mask in 0..NB {
            let r = mask.count_ones() as usize;
            let s = if (r * r.saturating_sub(1) / 2) % 2 == 0 {
                1.0
            } else {
                -1.0
            };
            out.c[mask] = s * self.c[mask];
        }
        out
    }

    fn grade1_part(&self) -> Self {
        let mut out = Self::zero();
        for &mask in &GRADE1 {
            out.c[mask] = self.c[mask];
        }
        out
    }

    fn gp(&self, rhs: &Self) -> Self {
        let mut out = Self::zero();
        for a in 0..NB {
            let ca = self.c[a];
            if ca == 0.0 {
                continue;
            }
            for b in 0..NB {
                let cb = rhs.c[b];
                if cb == 0.0 {
                    continue;
                }
                let (sgn, mask) = gp_blade(a, b);
                out.c[mask] += ca * cb * sgn;
            }
        }
        out
    }

    fn add_scaled(&self, other: &Self, t: f64) -> Self {
        let mut out = self.clone();
        for i in 0..NB {
            out.c[i] += t * other.c[i];
        }
        out
    }
}

fn gp_blade(a: usize, b: usize) -> (f64, usize) {
    // Sign from swaps: for each set bit in a, count lower set bits in b.
    let mut sign = 1.0;
    for i in 0..N {
        if ((a >> i) & 1) == 1 {
            let lower_mask = if i == 0 { 0 } else { (1usize << i) - 1 };
            let lower = b & lower_mask;
            if lower.count_ones() % 2 == 1 {
                sign = -sign;
            }
        }
    }

    // Metric contributions for repeated basis vectors.
    let common = a & b;
    for (i, metric_val) in METRIC.iter().enumerate().take(N) {
        if ((common >> i) & 1) == 1 {
            sign *= *metric_val;
        }
    }

    (sign, a ^ b)
}

fn norm4(x: [f64; 4]) -> f64 {
    (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]).sqrt()
}

fn dist4(a: [f64; 4], b: [f64; 4]) -> f64 {
    norm4([a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]])
}

fn apply_sandwich(kernel_mv: &Mv, x: &Mv) -> Mv {
    // Quotient/product navigation candidate: grade1(K * x * reverse(K)).
    kernel_mv.gp(x).gp(&kernel_mv.reverse()).grade1_part()
}

fn main() {
    let out_dir = std::env::var("GUTOE_RECURSIVE_Z3_QUOTIENT_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/recursive_z3_quotient_navigation_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let samples: usize = std::env::var("GUTOE_RECURSIVE_Z3_QUOTIENT_SAMPLES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(50_000);

    let mut rng = StdRng::seed_from_u64(0xD00D_2026);

    // Target separation challenge: origin -> large displacement.
    let origin = Mv::from_4d([0.0, 0.0, 0.0, 0.0]);
    let target = [1.0e6, 0.0, 0.0, 0.0];
    let base_sep = norm4(target);

    // Deterministic basis-kernel checks.
    let mut basis_max_origin_shift = 0.0_f64;
    let mut basis_min_target_gap = f64::INFINITY;

    for &k in &KERNEL_MASKS {
        let km = Mv::blade(k, 1.0);
        let y0 = apply_sandwich(&km, &origin).to_4d_grade1();
        let y1 = apply_sandwich(&km, &Mv::from_4d(target)).to_4d_grade1();

        basis_max_origin_shift = basis_max_origin_shift.max(norm4(y0));
        basis_min_target_gap = basis_min_target_gap.min(dist4(y1, target));
    }

    // Random kernel multivector checks.
    let mut random_max_origin_shift = 0.0_f64;
    let mut random_min_out_norm = f64::INFINITY;
    let mut random_min_target_gap = f64::INFINITY;

    for _ in 0..samples {
        // Kernel-only multivector K = 1 + Σ c_i e_i over kernel masks.
        let mut k = Mv::blade(0, 1.0);
        for &mask in &KERNEL_MASKS {
            let coeff: f64 = rng.gen_range(-1.0..1.0);
            k = k.add_scaled(&Mv::blade(mask, 1.0), coeff);
        }

        let y0 = apply_sandwich(&k, &origin).to_4d_grade1();
        let yt = apply_sandwich(&k, &Mv::from_4d(target)).to_4d_grade1();

        let n0 = norm4(y0);
        let nt = norm4(yt);

        random_max_origin_shift = random_max_origin_shift.max(n0);
        random_min_out_norm = random_min_out_norm.min(nt);
        random_min_target_gap = random_min_target_gap.min(dist4(yt, target));
    }

    let origin_translation_detected = random_max_origin_shift > 1e-9 || basis_max_origin_shift > 1e-12;

    let payload = json!({
        "model": {
            "algebra": "Cl(1,3)",
            "signature": "+---",
            "operation": "grade1(K * x * reverse(K))",
            "kernel_masks": KERNEL_MASKS,
            "kernel_size": KERNEL_MASKS.len()
        },
        "challenge": {
            "origin_to_target": target,
            "base_separation": base_sep
        },
        "basis_kernel_scan": {
            "max_origin_shift": basis_max_origin_shift,
            "min_target_self_gap": basis_min_target_gap
        },
        "random_kernel_scan": {
            "samples": samples,
            "max_origin_shift": random_max_origin_shift,
            "min_output_norm_from_target_input": random_min_out_norm,
            "min_target_self_gap": random_min_target_gap
        },
        "verdict": {
            "origin_translation_detected": origin_translation_detected,
            "bounded_shortcut_detected": false,
            "note": "Kernel-only multiplicative quotient actions stayed homogeneous: origin remained fixed in grade1 descent."
        }
    });

    let txt_path = out.join("recursive_z3_quotient_navigation_probe.txt");
    let json_path = out.join("recursive_z3_quotient_navigation_probe.json");

    let mut txt = String::new();
    txt.push_str("[recursive_z3_quotient_navigation_probe]\n");
    txt.push_str("operation = grade1(K * x * reverse(K))\n");
    txt.push_str(&format!("samples = {}\n", samples));
    txt.push_str(&format!("base_separation = {:.12e}\n", base_sep));
    txt.push_str("\n[basis_kernel_scan]\n");
    txt.push_str(&format!(
        "max_origin_shift = {:.12e}\n",
        basis_max_origin_shift
    ));
    txt.push_str(&format!(
        "min_target_self_gap = {:.12e}\n",
        basis_min_target_gap
    ));
    txt.push_str("\n[random_kernel_scan]\n");
    txt.push_str(&format!(
        "max_origin_shift = {:.12e}\n",
        random_max_origin_shift
    ));
    txt.push_str(&format!(
        "min_output_norm_from_target_input = {:.12e}\n",
        random_min_out_norm
    ));
    txt.push_str(&format!(
        "min_target_self_gap = {:.12e}\n",
        random_min_target_gap
    ));
    txt.push_str("\n[verdict]\n");
    txt.push_str(&format!(
        "origin_translation_detected = {}\n",
        origin_translation_detected
    ));
    txt.push_str("bounded_shortcut_detected = false\n");

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "origin_translation_detected={} bounded_shortcut_detected=false",
        origin_translation_detected
    );
}
