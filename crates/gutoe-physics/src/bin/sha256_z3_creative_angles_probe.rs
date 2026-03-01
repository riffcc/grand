//! SHA-256 creative Z3 angles probe (16, 17, 18, 19).
//!
//! 16) Length/padding boundary scan.
//! 17) Retrocompute-like branch search over transform candidates.
//! 18) Canonical-orbit deduper speedup + error measurement.
//! 19) Wrapped construction H'(m)=decode(H(encode(m))) with induced Z3 behavior.

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

fn is_order3_on_sample(t: Z3Transform, sample: &[u8]) -> bool {
    let p1 = apply_z3_bytes(sample, t);
    let p2 = apply_z3_bytes(&p1, t);
    let p3 = apply_z3_bytes(&p2, t);
    p3 == sample
}

fn canonical_rotation3_block(a: u8, b: u8, c: u8) -> (u8, u8, u8) {
    let r0 = [a, b, c];
    let r1 = [b, c, a];
    let r2 = [c, a, b];
    match r0.cmp(&r1) {
        Ordering::Less | Ordering::Equal => {
            if r0 <= r2 {
                (r0[0], r0[1], r0[2])
            } else {
                (r2[0], r2[1], r2[2])
            }
        }
        Ordering::Greater => {
            if r1 <= r2 {
                (r1[0], r1[1], r1[2])
            } else {
                (r2[0], r2[1], r2[2])
            }
        }
    }
}

fn encode_z3_triplet_canonical(msg: &[u8]) -> Vec<u8> {
    let mut out = msg.to_vec();
    let mut i = 0usize;
    while i + 2 < msg.len() {
        let (x, y, z) = canonical_rotation3_block(msg[i], msg[i + 1], msg[i + 2]);
        out[i] = x;
        out[i + 1] = y;
        out[i + 2] = z;
        i += 3;
    }
    out
}

fn decode_z3_triplet_canonical_digest(h: &[u8; 32]) -> [u8; 32] {
    let mut out = *h;
    let mut i = 0usize;
    while i + 2 < 30 {
        let (x, y, z) = canonical_rotation3_block(h[i], h[i + 1], h[i + 2]);
        out[i] = x;
        out[i + 1] = y;
        out[i + 2] = z;
        i += 3;
    }
    out
}

fn h_prime(msg: &[u8]) -> [u8; 32] {
    let encoded = encode_z3_triplet_canonical(msg);
    let h = sha256(&encoded);
    decode_z3_triplet_canonical_digest(&h)
}

#[derive(Clone, Debug)]
struct BranchResult {
    pi: Z3Transform,
    rho: Z3Transform,
    hits: u64,
    trials: u64,
}

fn main() {
    let out_dir = std::env::var("GUTOE_SHA256_Z3_CREATIVE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/sha256_z3_creative_angles_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let trials_len_scan = env_u64("GUTOE_SHA256_Z3_LEN_TRIALS", 5_000);
    let trials_branch = env_u64("GUTOE_SHA256_Z3_BRANCH_TRIALS", 20_000);
    let trials_wrap = env_u64("GUTOE_SHA256_Z3_WRAP_TRIALS", 20_000);
    let retro_shift_factor = env_f64("GUTOE_SHA256_Z3_RETRO_SHIFT_FACTOR", 1.20).max(0.0);
    let observer_floor_s = env_f64("GUTOE_SHA256_Z3_OBSERVER_FLOOR_S", 1.0e-9).max(1.0e-15);

    let msg_candidates_16 = [
        Z3Transform::TripletBytes { offset: 0 },
        Z3Transform::TripletWords { word_offset: 0 },
    ];
    let dig_candidates_16 = [
        Z3Transform::TripletBytes { offset: 0 },
        Z3Transform::TripletWords { word_offset: 0 },
    ];
    let lengths: [usize; 12] = [0, 1, 2, 3, 4, 31, 32, 33, 47, 48, 54, 55];

    // ----------------------------------------------------------------
    // 16) Length/padding boundary scan
    // ----------------------------------------------------------------
    let mut rng = OsRng;
    let mut len_rows = Vec::new();
    for &len in &lengths {
        for &pi in &msg_candidates_16 {
            for &rho in &dig_candidates_16 {
                let mut hits = 0_u64;
                for _ in 0..trials_len_scan {
                    let mut m = vec![0u8; len];
                    rng.fill_bytes(&mut m);
                    let hm = sha256(&m);
                    let pm = apply_z3_bytes(&m, pi);
                    let hpm = sha256(&pm);
                    let rhm = apply_z3_digest(&hm, rho);
                    if hpm == rhm {
                        hits += 1;
                    }
                }
                len_rows.push(json!({
                  "len": len,
                  "pi": pi.name(),
                  "rho": rho.name(),
                  "hits": hits,
                  "trials": trials_len_scan,
                  "hit_rate": hits as f64 / trials_len_scan as f64
                }));
            }
        }
    }

    // ----------------------------------------------------------------
    // 17) Retrocompute-like branch search over transform candidates
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

    let t_branch_start = Instant::now();
    let chunk = (pairs.len() + threads - 1) / threads;
    let mut handles = Vec::new();

    for tid in 0..threads {
        let lo = tid * chunk;
        let hi = ((tid + 1) * chunk).min(pairs.len());
        if lo >= hi {
            continue;
        }
        let slice = pairs[lo..hi].to_vec();
        handles.push(thread::spawn(move || {
            let mut out_local = Vec::new();
            let mut seed = [0u8; 32];
            seed[..8].copy_from_slice(&(tid as u64).to_le_bytes());
            let mut rng = StdRng::from_seed(seed);
            for (pi, rho) in slice {
                let mut hits = 0_u64;
                for _ in 0..trials_branch {
                    let mut m = [0u8; 48];
                    rng.fill_bytes(&mut m);
                    let hm = sha256(&m);
                    let pm = apply_z3_bytes(&m, pi);
                    let hpm = sha256(&pm);
                    let rhm = apply_z3_digest(&hm, rho);
                    if hpm == rhm {
                        hits += 1;
                    }
                }
                out_local.push(BranchResult {
                    pi,
                    rho,
                    hits,
                    trials: trials_branch,
                });
            }
            out_local
        }));
    }

    let mut branch_rows = Vec::new();
    for h in handles {
        let rows = h.join().expect("branch worker should join");
        for r in rows {
            branch_rows.push(r);
        }
    }
    let branch_elapsed_s = t_branch_start.elapsed().as_secs_f64();
    branch_rows.sort_by(|a, b| b.hits.cmp(&a.hits));
    let best = branch_rows
        .first()
        .cloned()
        .expect("branch search must have at least one result");
    let best_hit_rate = best.hits as f64 / best.trials as f64;

    let normal_latency_s = branch_elapsed_s.max(observer_floor_s);
    let observed_latency_s = if retro_shift_factor > 1.0 {
        observer_floor_s
    } else {
        (normal_latency_s * (1.0 - retro_shift_factor)).max(observer_floor_s)
    };
    let predeparture = retro_shift_factor > 1.0;
    let apparent_speedup = normal_latency_s / observed_latency_s;

    // ----------------------------------------------------------------
    // 18) Canonical orbit deduper speedup + validity
    // ----------------------------------------------------------------
    // Domain: all 3-byte messages with bytes in 0..15 (16^3 = 4096).
    let mut baseline_true = 0_u64;
    let mut baseline_hashes = 0_u64;
    let mut orbit_rep_prop: HashMap<[u8; 3], bool> = HashMap::new();
    let mut orbit_size: HashMap<[u8; 3], u64> = HashMap::new();
    let mut actual_prop: HashMap<[u8; 3], bool> = HashMap::new();

    for a in 0u8..16 {
        for b in 0u8..16 {
            for c in 0u8..16 {
                let m = [a, b, c];
                let h = sha256(&m);
                let prop = h[0] == 0; // toy predicate
                baseline_hashes += 1;
                if prop {
                    baseline_true += 1;
                }
                actual_prop.insert(m, prop);

                let r1 = [b, c, a];
                let r2 = [c, a, b];
                let canon = if m <= r1 && m <= r2 {
                    m
                } else if r1 <= m && r1 <= r2 {
                    r1
                } else {
                    r2
                };

                orbit_size.entry(canon).and_modify(|v| *v += 1).or_insert(1);
                orbit_rep_prop.entry(canon).or_insert(prop);
            }
        }
    }

    let dedup_hashes = orbit_rep_prop.len() as u64;
    let ideal_speedup = baseline_hashes as f64 / dedup_hashes as f64;

    let mut predicted_true = 0_u64;
    let mut mismatches = 0_u64;
    for a in 0u8..16 {
        for b in 0u8..16 {
            for c in 0u8..16 {
                let m = [a, b, c];
                let r1 = [b, c, a];
                let r2 = [c, a, b];
                let canon = if m <= r1 && m <= r2 {
                    m
                } else if r1 <= m && r1 <= r2 {
                    r1
                } else {
                    r2
                };
                let pred = *orbit_rep_prop.get(&canon).expect("canon rep should exist");
                let truth = *actual_prop.get(&m).expect("truth should exist");
                if pred {
                    predicted_true += 1;
                }
                if pred != truth {
                    mismatches += 1;
                }
            }
        }
    }
    let mismatch_rate = mismatches as f64 / baseline_hashes as f64;

    // ----------------------------------------------------------------
    // 19) Wrapped construction H'(m)=decode(H(encode(m)))
    // ----------------------------------------------------------------
    let mut wrap_hits = 0_u64;
    let mut wrap_rng = OsRng;
    let pi_wrap = Z3Transform::TripletBytes { offset: 0 };
    for _ in 0..trials_wrap {
        let mut m = [0u8; 48];
        wrap_rng.fill_bytes(&mut m);
        let p = apply_z3_bytes(&m, pi_wrap);
        let h1 = h_prime(&m);
        let h2 = h_prime(&p);
        if h1 == h2 {
            wrap_hits += 1;
        }
    }
    let wrap_hit_rate = wrap_hits as f64 / trials_wrap as f64;

    // Order-3 checks on canonical samples for reporting sanity.
    let sample48: [u8; 48] = core::array::from_fn(|i| i as u8);
    let sample32: [u8; 32] = core::array::from_fn(|i| (255 - i as u8));
    let mut order3_checks = Vec::new();
    for &t in &msg_candidates_17 {
        order3_checks.push(json!({
          "domain": "message",
          "transform": t.name(),
          "order3_identity_on_sample": is_order3_on_sample(t, &sample48)
        }));
    }
    for &t in &dig_candidates_17 {
        let p1 = apply_z3_digest(&sample32, t);
        let p2 = apply_z3_digest(&p1, t);
        let p3 = apply_z3_digest(&p2, t);
        order3_checks.push(json!({
          "domain": "digest",
          "transform": t.name(),
          "order3_identity_on_sample": p3 == sample32
        }));
    }

    let payload = json!({
      "scope": "creative angles 16/17/18/19",
      "order3_checks": order3_checks,
      "angle_16_length_padding_scan": {
        "trials_per_case": trials_len_scan,
        "lengths": lengths,
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
        "best_pair": {
          "pi": best.pi.name(),
          "rho": best.rho.name(),
          "hits": best.hits,
          "trials": best.trials,
          "hit_rate": best_hit_rate
        },
        "all_pairs": branch_rows.iter().map(|r| json!({
          "pi": r.pi.name(),
          "rho": r.rho.name(),
          "hits": r.hits,
          "trials": r.trials,
          "hit_rate": r.hits as f64 / r.trials as f64
        })).collect::<Vec<_>>()
      },
      "angle_18_orbit_deduper": {
        "domain": "3-byte messages over alphabet 0..15 (4096 total)",
        "predicate": "sha256(m)[0] == 0",
        "baseline_hashes": baseline_hashes,
        "baseline_true": baseline_true,
        "dedup_hashes": dedup_hashes,
        "ideal_speedup_if_equivariant": ideal_speedup,
        "predicted_true_from_reps": predicted_true,
        "mismatches_vs_truth": mismatches,
        "mismatch_rate": mismatch_rate
      },
      "angle_19_wrapped_construction": {
        "definition": "H'(m)=decode(sha256(encode(m)))",
        "encode": "triplet-wise canonical rotation",
        "decode": "digest triplet-wise canonical rotation",
        "test_action": pi_wrap.name(),
        "trials": trials_wrap,
        "invariance_hits": wrap_hits,
        "invariance_hit_rate": wrap_hit_rate
      }
    });

    let txt_path = out.join("sha256_z3_creative_angles_probe.txt");
    let json_path = out.join("sha256_z3_creative_angles_probe.json");

    let mut txt = String::new();
    txt.push_str("[sha256_z3_creative_angles_probe]\n");
    txt.push_str("angles: 16(length scan), 17(branch search), 18(orbit deduper), 19(wrapped H')\n\n");
    txt.push_str("ANGLE 16\n");
    txt.push_str(&format!("trials_per_case = {}\n", trials_len_scan));
    for row in payload["angle_16_length_padding_scan"]["results"]
        .as_array()
        .expect("array")
    {
        txt.push_str(&format!(
            "len={} pi={} rho={} hits={}/{} rate={:.6e}\n",
            row["len"].as_u64().unwrap_or_default(),
            row["pi"].as_str().unwrap_or(""),
            row["rho"].as_str().unwrap_or(""),
            row["hits"].as_u64().unwrap_or_default(),
            row["trials"].as_u64().unwrap_or_default(),
            row["hit_rate"].as_f64().unwrap_or_default()
        ));
    }

    txt.push_str("\nANGLE 17\n");
    txt.push_str(&format!("pairs={} threads={}\n", pairs.len(), threads));
    txt.push_str(&format!("host_elapsed_s={:.6e}\n", branch_elapsed_s));
    txt.push_str(&format!("predeparture={}\n", predeparture));
    txt.push_str(&format!("apparent_speedup={:.6e}\n", apparent_speedup));
    txt.push_str(&format!(
        "best_pair pi={} rho={} hits={}/{} rate={:.6e}\n",
        best.pi.name(),
        best.rho.name(),
        best.hits,
        best.trials,
        best_hit_rate
    ));

    txt.push_str("\nANGLE 18\n");
    txt.push_str(&format!("baseline_hashes={}\n", baseline_hashes));
    txt.push_str(&format!("dedup_hashes={}\n", dedup_hashes));
    txt.push_str(&format!("ideal_speedup_if_equivariant={:.6e}\n", ideal_speedup));
    txt.push_str(&format!("mismatch_rate={:.6e}\n", mismatch_rate));

    txt.push_str("\nANGLE 19\n");
    txt.push_str(&format!(
        "wrapped_invariance hits={}/{} rate={:.6e}\n",
        wrap_hits, trials_wrap, wrap_hit_rate
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

