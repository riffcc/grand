//! Double-SHA256 (H2 = SHA256 ∘ SHA256) Z3 mining-structure probe.
//!
//! Angles implemented:
//! - 16: length/padding boundary scan
//! - 17: parallel branch search over transform candidates
//! - 18: canonical-orbit deduper speedup/error for mining predicate
//! - 19: direct H2 symmetry/bias test (exact + approximate)
//!
//! Note: this uses a single-block SHA-256 implementation (input <= 55 bytes),
//! so it is a controlled structural lane, not full 80-byte Bitcoin header logic.

use rand::rngs::{OsRng, StdRng};
use rand::{RngCore, SeedableRng};
use serde_json::json;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::thread;
use std::time::Instant;

const K256: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

const H256_INIT: [u32; 8] = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
];

fn env_u64(name: &str, default: u64) -> u64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(default)
}

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn sha256_compress(state: &mut [u32; 8], block: &[u8; 64]) {
    let mut w = [0u32; 64];
    for i in 0..16 {
        w[i] = u32::from_be_bytes([
            block[4 * i],
            block[4 * i + 1],
            block[4 * i + 2],
            block[4 * i + 3],
        ]);
    }
    for i in 16..64 {
        let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
        let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16]
            .wrapping_add(s0)
            .wrapping_add(w[i - 7])
            .wrapping_add(s1);
    }

    let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = *state;
    for i in 0..64 {
        let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
        let ch = (e & f) ^ (!e & g);
        let temp1 = h
            .wrapping_add(s1)
            .wrapping_add(ch)
            .wrapping_add(K256[i])
            .wrapping_add(w[i]);
        let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
        let maj = (a & b) ^ (a & c) ^ (b & c);
        let temp2 = s0.wrapping_add(maj);
        h = g;
        g = f;
        f = e;
        e = d.wrapping_add(temp1);
        d = c;
        c = b;
        b = a;
        a = temp1.wrapping_add(temp2);
    }

    state[0] = state[0].wrapping_add(a);
    state[1] = state[1].wrapping_add(b);
    state[2] = state[2].wrapping_add(c);
    state[3] = state[3].wrapping_add(d);
    state[4] = state[4].wrapping_add(e);
    state[5] = state[5].wrapping_add(f);
    state[6] = state[6].wrapping_add(g);
    state[7] = state[7].wrapping_add(h);
}

/// SHA-256 of message up to 55 bytes (single-block lane).
fn sha256(msg: &[u8]) -> [u8; 32] {
    assert!(msg.len() <= 55, "single-block SHA-256 only");
    let mut block = [0u8; 64];
    block[..msg.len()].copy_from_slice(msg);
    block[msg.len()] = 0x80;
    let bit_len = (msg.len() as u64) * 8;
    block[56..64].copy_from_slice(&bit_len.to_be_bytes());

    let mut state = H256_INIT;
    sha256_compress(&mut state, &block);
    let mut out = [0u8; 32];
    for i in 0..8 {
        out[4 * i..4 * i + 4].copy_from_slice(&state[i].to_be_bytes());
    }
    out
}

fn sha256d(msg: &[u8]) -> [u8; 32] {
    let inner = sha256(msg);
    sha256(&inner)
}

#[derive(Clone, Copy, Debug)]
enum Z3Transform {
    TripletBytes { offset: usize },
    TripletWords { word_offset: usize },
}

impl Z3Transform {
    fn name(self) -> String {
        match self {
            Z3Transform::TripletBytes { offset } => format!("triplet_bytes_off{}", offset),
            Z3Transform::TripletWords { word_offset } => format!("triplet_words_off{}", word_offset),
        }
    }
}

fn apply_z3_bytes(input: &[u8], t: Z3Transform) -> Vec<u8> {
    let mut out = input.to_vec();
    match t {
        Z3Transform::TripletBytes { offset } => {
            let off = offset.min(2);
            let mut i = off;
            while i + 2 < input.len() {
                out[i] = input[i + 1];
                out[i + 1] = input[i + 2];
                out[i + 2] = input[i];
                i += 3;
            }
        }
        Z3Transform::TripletWords { word_offset } => {
            let off = (word_offset.min(2)) * 4;
            let mut i = off;
            while i + 11 < input.len() {
                out[i..i + 4].copy_from_slice(&input[i + 4..i + 8]);
                out[i + 4..i + 8].copy_from_slice(&input[i + 8..i + 12]);
                out[i + 8..i + 12].copy_from_slice(&input[i..i + 4]);
                i += 12;
            }
        }
    }
    out
}

fn apply_z3_digest(h: &[u8; 32], t: Z3Transform) -> [u8; 32] {
    let v = apply_z3_bytes(h, t);
    let mut out = [0u8; 32];
    out.copy_from_slice(&v);
    out
}

fn digest_hamming_bits(a: &[u8; 32], b: &[u8; 32]) -> u32 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x ^ y).count_ones())
        .sum()
}

fn leading_zero_bits(h: &[u8; 32]) -> u32 {
    let mut n = 0u32;
    for &b in h {
        if b == 0 {
            n += 8;
        } else {
            n += b.leading_zeros();
            break;
        }
    }
    n
}

#[derive(Clone, Debug)]
struct BranchResult {
    pi: Z3Transform,
    rho: Z3Transform,
    exact_hits: u64,
    trials: u64,
    mean_hd_bits: f64,
}

fn main() {
    let out_dir = std::env::var("GUTOE_SHA256D_Z3_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/sha256d_z3_mining_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let trials_len_scan = env_u64("GUTOE_SHA256D_Z3_LEN_TRIALS", 5_000);
    let trials_branch = env_u64("GUTOE_SHA256D_Z3_BRANCH_TRIALS", 5_000);
    let trials_comp = env_u64("GUTOE_SHA256D_Z3_COMP_TRIALS", 100_000);
    let mining_lz_threshold = env_u64("GUTOE_SHA256D_Z3_LZ_THRESHOLD", 12).min(256) as u32;
    let retro_shift_factor = env_f64("GUTOE_SHA256D_Z3_RETRO_SHIFT_FACTOR", 1.20).max(0.0);
    let observer_floor_s = env_f64("GUTOE_SHA256D_Z3_OBSERVER_FLOOR_S", 1.0e-9).max(1.0e-15);

    let lengths: [usize; 12] = [0, 1, 2, 3, 4, 31, 32, 33, 47, 48, 54, 55];
    let msg_candidates_16 = [
        Z3Transform::TripletBytes { offset: 0 },
        Z3Transform::TripletWords { word_offset: 0 },
    ];
    let dig_candidates_16 = [
        Z3Transform::TripletBytes { offset: 0 },
        Z3Transform::TripletWords { word_offset: 0 },
    ];

    // ----------------------------------------------------------------
    // 16) Length/padding boundary scan on H2
    // ----------------------------------------------------------------
    let mut rng = OsRng;
    let mut len_rows = Vec::new();
    for &len in &lengths {
        for &pi in &msg_candidates_16 {
            for &rho in &dig_candidates_16 {
                let mut hits = 0u64;
                let mut hd_sum = 0u64;
                for _ in 0..trials_len_scan {
                    let mut m = vec![0u8; len];
                    rng.fill_bytes(&mut m);
                    let hm = sha256d(&m);
                    let pm = apply_z3_bytes(&m, pi);
                    let hpm = sha256d(&pm);
                    let rhm = apply_z3_digest(&hm, rho);
                    if hpm == rhm {
                        hits += 1;
                    }
                    hd_sum += digest_hamming_bits(&hpm, &rhm) as u64;
                }
                len_rows.push(json!({
                  "len": len,
                  "pi": pi.name(),
                  "rho": rho.name(),
                  "exact_hits": hits,
                  "trials": trials_len_scan,
                  "exact_hit_rate": hits as f64 / trials_len_scan as f64,
                  "mean_hd_bits": hd_sum as f64 / trials_len_scan as f64
                }));
            }
        }
    }

    // ----------------------------------------------------------------
    // 17) Parallel candidate search on H2 (exact + approx metrics)
    // ----------------------------------------------------------------
    let msg_candidates_17 = vec![
        Z3Transform::TripletBytes { offset: 0 },
        Z3Transform::TripletBytes { offset: 1 },
        Z3Transform::TripletBytes { offset: 2 },
        Z3Transform::TripletWords { word_offset: 0 },
        Z3Transform::TripletWords { word_offset: 1 },
        Z3Transform::TripletWords { word_offset: 2 },
    ];
    let dig_candidates_17 = msg_candidates_17.clone();

    let mut pairs = Vec::new();
    for &pi in &msg_candidates_17 {
        for &rho in &dig_candidates_17 {
            pairs.push((pi, rho));
        }
    }

    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .min(pairs.len().max(1));
    let chunk = (pairs.len() + threads - 1) / threads;
    let t_branch_start = Instant::now();
    let mut handles = Vec::new();

    for tid in 0..threads {
        let lo = tid * chunk;
        let hi = ((tid + 1) * chunk).min(pairs.len());
        if lo >= hi {
            continue;
        }
        let local = pairs[lo..hi].to_vec();
        handles.push(thread::spawn(move || {
            let mut out_local = Vec::new();
            let mut seed = [0u8; 32];
            seed[..8].copy_from_slice(&(tid as u64).to_le_bytes());
            let mut rng = StdRng::from_seed(seed);
            for (pi, rho) in local {
                let mut hits = 0u64;
                let mut hd_sum = 0u64;
                for _ in 0..trials_branch {
                    let mut m = [0u8; 48];
                    rng.fill_bytes(&mut m);
                    let hm = sha256d(&m);
                    let pm = apply_z3_bytes(&m, pi);
                    let hpm = sha256d(&pm);
                    let rhm = apply_z3_digest(&hm, rho);
                    if hpm == rhm {
                        hits += 1;
                    }
                    hd_sum += digest_hamming_bits(&hpm, &rhm) as u64;
                }
                out_local.push(BranchResult {
                    pi,
                    rho,
                    exact_hits: hits,
                    trials: trials_branch,
                    mean_hd_bits: hd_sum as f64 / trials_branch as f64,
                });
            }
            out_local
        }));
    }

    let mut branch_rows = Vec::new();
    for h in handles {
        for row in h.join().expect("worker join") {
            branch_rows.push(row);
        }
    }
    let branch_elapsed_s = t_branch_start.elapsed().as_secs_f64();
    branch_rows.sort_by(|a, b| {
        a.mean_hd_bits
            .partial_cmp(&b.mean_hd_bits)
            .unwrap_or(Ordering::Equal)
    });
    let best_hd = branch_rows
        .first()
        .cloned()
        .expect("branch results must be non-empty");
    let mut by_exact = branch_rows.clone();
    by_exact.sort_by(|a, b| b.exact_hits.cmp(&a.exact_hits));
    let best_exact = by_exact
        .first()
        .cloned()
        .expect("branch results must be non-empty");

    let normal_latency_s = branch_elapsed_s.max(observer_floor_s);
    let observed_latency_s = if retro_shift_factor > 1.0 {
        observer_floor_s
    } else {
        (normal_latency_s * (1.0 - retro_shift_factor)).max(observer_floor_s)
    };
    let predeparture = retro_shift_factor > 1.0;
    let apparent_speedup = normal_latency_s / observed_latency_s;

    // ----------------------------------------------------------------
    // 18) Orbit deduper speedup vs correctness on mining predicate
    // ----------------------------------------------------------------
    // Domain: all 3-byte messages with values in 0..15 (4096 total), closed under byte rotation.
    let mut actual_pred: HashMap<[u8; 3], bool> = HashMap::new();
    let mut rep_pred: HashMap<[u8; 3], bool> = HashMap::new();
    let mut baseline_hashes = 0u64;
    let mut baseline_success = 0u64;

    for a in 0u8..16 {
        for b in 0u8..16 {
            for c in 0u8..16 {
                let m = [a, b, c];
                let h = sha256d(&m);
                let pred = leading_zero_bits(&h) >= mining_lz_threshold;
                baseline_hashes += 1;
                if pred {
                    baseline_success += 1;
                }
                actual_pred.insert(m, pred);

                let r1 = [b, c, a];
                let r2 = [c, a, b];
                let canon = if m <= r1 && m <= r2 {
                    m
                } else if r1 <= m && r1 <= r2 {
                    r1
                } else {
                    r2
                };
                rep_pred.entry(canon).or_insert(pred);
            }
        }
    }

    let dedup_hashes = rep_pred.len() as u64;
    let ideal_speedup = baseline_hashes as f64 / dedup_hashes as f64;
    let mut mismatches = 0u64;
    let mut predicted_success = 0u64;
    for (&m, &truth) in &actual_pred {
        let r1 = [m[1], m[2], m[0]];
        let r2 = [m[2], m[0], m[1]];
        let canon = if m <= r1 && m <= r2 {
            m
        } else if r1 <= m && r1 <= r2 {
            r1
        } else {
            r2
        };
        let pred = *rep_pred.get(&canon).expect("canon rep exists");
        if pred {
            predicted_success += 1;
        }
        if pred != truth {
            mismatches += 1;
        }
    }
    let mismatch_rate = mismatches as f64 / baseline_hashes as f64;

    // ----------------------------------------------------------------
    // 19) Direct composed-function H2 symmetry/bias probe
    // ----------------------------------------------------------------
    let pi_19 = best_hd.pi;
    let rho_19 = best_hd.rho;
    let mut eq_hits = 0u64;
    let mut hd_sum = 0u64;
    let mut hd_random_sum = 0u64;
    let mut z0_sum = 0f64;
    let mut z1_sum = 0f64;
    let mut z0_sq_sum = 0f64;
    let mut z1_sq_sum = 0f64;
    let mut z01_sum = 0f64;
    let mut s0 = 0u64;
    let mut s1 = 0u64;
    let mut agree = 0u64;

    let mut rng19 = OsRng;
    for _ in 0..trials_comp {
        let mut m = [0u8; 48];
        let mut mr = [0u8; 48];
        rng19.fill_bytes(&mut m);
        rng19.fill_bytes(&mut mr);

        let hm = sha256d(&m);
        let pm = apply_z3_bytes(&m, pi_19);
        let hpm = sha256d(&pm);
        let rhm = apply_z3_digest(&hm, rho_19);
        if hpm == rhm {
            eq_hits += 1;
        }
        hd_sum += digest_hamming_bits(&hpm, &rhm) as u64;

        let hr = sha256d(&mr);
        hd_random_sum += digest_hamming_bits(&hpm, &hr) as u64;

        let z0 = leading_zero_bits(&hm) as f64;
        let z1 = leading_zero_bits(&hpm) as f64;
        z0_sum += z0;
        z1_sum += z1;
        z0_sq_sum += z0 * z0;
        z1_sq_sum += z1 * z1;
        z01_sum += z0 * z1;

        let b0 = (z0 as u32) >= mining_lz_threshold;
        let b1 = (z1 as u32) >= mining_lz_threshold;
        if b0 {
            s0 += 1;
        }
        if b1 {
            s1 += 1;
        }
        if b0 == b1 {
            agree += 1;
        }
    }

    let n = trials_comp as f64;
    let mean_hd = hd_sum as f64 / n;
    let mean_hd_random = hd_random_sum as f64 / n;
    let mean_z0 = z0_sum / n;
    let mean_z1 = z1_sum / n;
    let cov = z01_sum / n - mean_z0 * mean_z1;
    let var0 = (z0_sq_sum / n - mean_z0 * mean_z0).max(0.0);
    let var1 = (z1_sq_sum / n - mean_z1 * mean_z1).max(0.0);
    let corr = if var0 > 0.0 && var1 > 0.0 {
        cov / (var0.sqrt() * var1.sqrt())
    } else {
        0.0
    };
    let p0 = s0 as f64 / n;
    let p1 = s1 as f64 / n;
    let agree_rate = agree as f64 / n;
    let agree_independent_baseline = p0 * p1 + (1.0 - p0) * (1.0 - p1);

    let payload = json!({
      "scope": "double_sha256_z3_mining_probe",
      "note": "single-block SHA256 lane (input <=55 bytes), controlled structure probe",
      "angle_16_length_padding_scan": {
        "trials_per_case": trials_len_scan,
        "results": len_rows
      },
      "angle_17_branch_search": {
        "trials_per_pair": trials_branch,
        "pair_count": pairs.len(),
        "threads_used": threads,
        "host_elapsed_s": branch_elapsed_s,
        "retro_shift_factor": retro_shift_factor,
        "predeparture": predeparture,
        "observed_latency_s": observed_latency_s,
        "apparent_speedup": apparent_speedup,
        "best_by_mean_hd": {
          "pi": best_hd.pi.name(),
          "rho": best_hd.rho.name(),
          "mean_hd_bits": best_hd.mean_hd_bits,
          "exact_hits": best_hd.exact_hits,
          "exact_hit_rate": best_hd.exact_hits as f64 / best_hd.trials as f64
        },
        "best_by_exact_hits": {
          "pi": best_exact.pi.name(),
          "rho": best_exact.rho.name(),
          "mean_hd_bits": best_exact.mean_hd_bits,
          "exact_hits": best_exact.exact_hits,
          "exact_hit_rate": best_exact.exact_hits as f64 / best_exact.trials as f64
        }
      },
      "angle_18_orbit_deduper": {
        "domain": "3-byte messages over alphabet 0..15 (4096 total)",
        "hash": "sha256d",
        "predicate": format!("leading_zero_bits(hash) >= {}", mining_lz_threshold),
        "baseline_hashes": baseline_hashes,
        "baseline_success": baseline_success,
        "dedup_hashes": dedup_hashes,
        "ideal_speedup_if_equivariant": ideal_speedup,
        "predicted_success_from_reps": predicted_success,
        "mismatches_vs_truth": mismatches,
        "mismatch_rate": mismatch_rate
      },
      "angle_19_composed_function_probe": {
        "function": "H2(m)=SHA256(SHA256(m))",
        "pi": pi_19.name(),
        "rho": rho_19.name(),
        "trials": trials_comp,
        "exact_hits": eq_hits,
        "exact_hit_rate": eq_hits as f64 / n,
        "mean_hd_bits_equivariance_pair": mean_hd,
        "mean_hd_bits_random_baseline": mean_hd_random,
        "mean_lz_h2_m": mean_z0,
        "mean_lz_h2_pi_m": mean_z1,
        "lz_corr": corr,
        "success_rate_m": p0,
        "success_rate_pi_m": p1,
        "success_agreement_rate": agree_rate,
        "success_agreement_independent_baseline": agree_independent_baseline
      }
    });

    let txt_path = out.join("sha256d_z3_mining_probe.txt");
    let json_path = out.join("sha256d_z3_mining_probe.json");

    let mut txt = String::new();
    txt.push_str("[sha256d_z3_mining_probe]\n");
    txt.push_str("double-SHA256 structure probe for angles 16/17/18/19\n\n");
    txt.push_str("ANGLE 17 summary\n");
    txt.push_str(&format!(
        "best_by_mean_hd: pi={} rho={} mean_hd={:.6} exact_hits={}/{}\n",
        best_hd.pi.name(),
        best_hd.rho.name(),
        best_hd.mean_hd_bits,
        best_hd.exact_hits,
        best_hd.trials
    ));
    txt.push_str(&format!(
        "branch_host_elapsed_s={:.6e} predeparture={} apparent_speedup={:.6e}\n\n",
        branch_elapsed_s, predeparture, apparent_speedup
    ));
    txt.push_str("ANGLE 18 summary\n");
    txt.push_str(&format!(
        "baseline_hashes={} dedup_hashes={} ideal_speedup_if_equivariant={:.6e}\n",
        baseline_hashes, dedup_hashes, ideal_speedup
    ));
    txt.push_str(&format!(
        "mismatch_rate={:.6e} (orbit inference error)\n\n",
        mismatch_rate
    ));
    txt.push_str("ANGLE 19 summary\n");
    txt.push_str(&format!(
        "exact_hits={}/{} rate={:.6e}\n",
        eq_hits,
        trials_comp,
        eq_hits as f64 / n
    ));
    txt.push_str(&format!(
        "mean_hd_pair={:.6} mean_hd_random={:.6}\n",
        mean_hd, mean_hd_random
    ));
    txt.push_str(&format!(
        "lz_corr={:.6e} agree_rate={:.6e} independent_baseline={:.6e}\n",
        corr, agree_rate, agree_independent_baseline
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
