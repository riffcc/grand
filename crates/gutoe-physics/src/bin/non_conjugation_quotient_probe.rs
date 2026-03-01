//! Non-conjugation quotient-product probe.
//!
//! Scans maps of the form:
//!   F_{L,R}(x) = grade1( L * X(x) * R )
//! with `R` independent of `reverse(L)` (non-conjugation lane).
//!
//! Reports distance-compression behavior and rank (injectivity proxy) on grade-1.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;
use std::fs;
use std::path::PathBuf;

const N: usize = 4;
const NB: usize = 1 << N;
const GRADE1_MASKS: [usize; 4] = [1, 2, 4, 8];
const KERNEL_MASKS: [usize; 12] = [0, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15];
const METRIC: [f64; N] = [1.0, -1.0, -1.0, -1.0];

#[derive(Clone)]
struct Mv {
    c: [f64; NB],
}

impl Mv {
    fn zero() -> Self {
        Self { c: [0.0; NB] }
    }

    fn blade(mask: usize, coeff: f64) -> Self {
        let mut out = Self::zero();
        out.c[mask] = coeff;
        out
    }

    fn from_vec4(v: [f64; 4]) -> Self {
        let mut out = Self::zero();
        for (i, &m) in GRADE1_MASKS.iter().enumerate() {
            out.c[m] = v[i];
        }
        out
    }

    fn to_vec4_grade1(&self) -> [f64; 4] {
        [
            self.c[GRADE1_MASKS[0]],
            self.c[GRADE1_MASKS[1]],
            self.c[GRADE1_MASKS[2]],
            self.c[GRADE1_MASKS[3]],
        ]
    }

    fn norm_coeff(&self) -> f64 {
        self.c.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    fn normalized(mut self) -> Self {
        let n = self.norm_coeff();
        if n > 0.0 {
            for i in 0..NB {
                self.c[i] /= n;
            }
        }
        self
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
    let common = a & b;
    for (i, metric_val) in METRIC.iter().enumerate().take(N) {
        if ((common >> i) & 1) == 1 {
            sign *= *metric_val;
        }
    }
    (sign, a ^ b)
}

fn norm4(v: [f64; 4]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]).sqrt()
}

fn mat_rank(mut a: [[f64; 4]; 4], eps: f64) -> usize {
    let mut rank = 0usize;
    let mut row = 0usize;
    for col in 0..4 {
        let mut pivot = row;
        let mut best = a[pivot][col].abs();
        for (r, rr) in a.iter().enumerate().skip(row + 1) {
            let v = rr[col].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best <= eps {
            continue;
        }
        a.swap(row, pivot);
        let diag = a[row][col];
        for c in col..4 {
            a[row][c] /= diag;
        }
        for r in 0..4 {
            if r == row {
                continue;
            }
            let factor = a[r][col];
            if factor.abs() <= eps {
                continue;
            }
            for c in col..4 {
                a[r][c] -= factor * a[row][c];
            }
        }
        rank += 1;
        row += 1;
        if row == 4 {
            break;
        }
    }
    rank
}

fn map_vec4_nonconj(l: &Mv, r: &Mv, x: [f64; 4]) -> [f64; 4] {
    let x_mv = Mv::from_vec4(x);
    l.gp(&x_mv).gp(r).to_vec4_grade1()
}

fn map_matrix_nonconj(l: &Mv, r: &Mv) -> [[f64; 4]; 4] {
    let e = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    let mut m = [[0.0; 4]; 4];
    for (j, ej) in e.iter().enumerate() {
        let y = map_vec4_nonconj(l, r, *ej);
        for i in 0..4 {
            m[i][j] = y[i];
        }
    }
    m
}

fn main() {
    let out_dir = std::env::var("GUTOE_NONCONJ_QUOTIENT_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/non_conjugation_quotient_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let target = [1.0e6, 0.0, 0.0, 0.0];
    let base = norm4(target);

    let mut best_kernel_any = f64::INFINITY;
    let mut best_kernel_rank4 = f64::INFINITY;
    let mut best_kernel_any_pair = (0usize, 0usize, 0usize);
    let mut best_kernel_rank4_pair = (0usize, 0usize, 0usize);

    // Exhaustive basis-blade kernel-only scan.
    for &lm in &KERNEL_MASKS {
        for &rm in &KERNEL_MASKS {
            let l = Mv::blade(lm, 1.0);
            let r = Mv::blade(rm, 1.0);
            let y = map_vec4_nonconj(&l, &r, target);
            let ratio = norm4(y) / base;
            let rank = mat_rank(map_matrix_nonconj(&l, &r), 1e-10);

            if ratio < best_kernel_any {
                best_kernel_any = ratio;
                best_kernel_any_pair = (lm, rm, rank);
            }
            if rank == 4 && ratio < best_kernel_rank4 {
                best_kernel_rank4 = ratio;
                best_kernel_rank4_pair = (lm, rm, rank);
            }
        }
    }

    // Random dense non-conjugation multivectors, normalized (unit coefficient norm).
    let samples: usize = std::env::var("GUTOE_NONCONJ_RANDOM_SAMPLES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(50_000);
    let mut rng = StdRng::seed_from_u64(0xBEEF_2026);

    let mut random_best_any = f64::INFINITY;
    let mut random_best_rank4 = f64::INFINITY;

    for _ in 0..samples {
        let mut l = Mv::zero();
        let mut r = Mv::zero();

        // Build from kernel masks only to stay in requested lane.
        for &m in &KERNEL_MASKS {
            l = l.add_scaled(&Mv::blade(m, 1.0), rng.gen_range(-1.0..1.0));
            r = r.add_scaled(&Mv::blade(m, 1.0), rng.gen_range(-1.0..1.0));
        }
        l = l.normalized();
        r = r.normalized();

        let y = map_vec4_nonconj(&l, &r, target);
        let ratio = norm4(y) / base;
        let rank = mat_rank(map_matrix_nonconj(&l, &r), 1e-9);

        if ratio < random_best_any {
            random_best_any = ratio;
        }
        if rank == 4 && ratio < random_best_rank4 {
            random_best_rank4 = ratio;
        }
    }

    let payload = json!({
        "model": {
            "operation": "grade1(L * X * R)",
            "lane": "non-conjugation quotient products",
            "domain": "kernel-restricted factors"
        },
        "target_norm": base,
        "basis_kernel_scan": {
            "best_ratio_any_rank": best_kernel_any,
            "best_pair_any_rank": {
                "L_mask": best_kernel_any_pair.0,
                "R_mask": best_kernel_any_pair.1,
                "rank": best_kernel_any_pair.2
            },
            "best_ratio_rank4": best_kernel_rank4,
            "best_pair_rank4": {
                "L_mask": best_kernel_rank4_pair.0,
                "R_mask": best_kernel_rank4_pair.1,
                "rank": best_kernel_rank4_pair.2
            }
        },
        "random_dense_scan": {
            "samples": samples,
            "best_ratio_any_rank": random_best_any,
            "best_ratio_rank4": random_best_rank4
        }
    });

    let txt_path = out.join("non_conjugation_quotient_probe.txt");
    let json_path = out.join("non_conjugation_quotient_probe.json");

    let mut txt = String::new();
    txt.push_str("[non_conjugation_quotient_probe]\n");
    txt.push_str("operation = grade1(L*X*R)\n");
    txt.push_str(&format!("target_norm = {:.12e}\n", base));
    txt.push_str(&format!(
        "basis_best_any = {:.12e} (L={},R={},rank={})\n",
        best_kernel_any, best_kernel_any_pair.0, best_kernel_any_pair.1, best_kernel_any_pair.2
    ));
    txt.push_str(&format!(
        "basis_best_rank4 = {:.12e} (L={},R={},rank={})\n",
        best_kernel_rank4, best_kernel_rank4_pair.0, best_kernel_rank4_pair.1, best_kernel_rank4_pair.2
    ));
    txt.push_str(&format!(
        "random_best_any = {:.12e}, random_best_rank4 = {:.12e}, samples={}\n",
        random_best_any, random_best_rank4, samples
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
