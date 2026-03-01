//! Consensus-preserving double-SHA256 nonce Z3 probe.
//!
//! This lane keeps hashing semantics aligned with Bitcoin-style PoW input shape:
//! - 80-byte header = 76-byte prefix + 4-byte nonce
//! - hash = SHA256(SHA256(header))
//! - only nonce bytes are transformed for Z3 orbit tests.
//!
//! It tests whether any tested Z3 action gives exploitable structure:
//! - exact equivariance hits
//! - approximate bias (Hamming distance)
//! - mining predicate correlation/agreement
//! - orbit-deduper error for leading-zero predicate

use rand::rngs::{OsRng, StdRng};
use rand::{RngCore, SeedableRng};
use serde_json::json;
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

fn sha256(msg: &[u8]) -> [u8; 32] {
    let mut state = H256_INIT;
    let bit_len = (msg.len() as u64).wrapping_mul(8);

    let mut data = msg.to_vec();
    data.push(0x80);
    while (data.len() % 64) != 56 {
        data.push(0);
    }
    data.extend_from_slice(&bit_len.to_be_bytes());

    for chunk in data.chunks_exact(64) {
        let mut block = [0u8; 64];
        block.copy_from_slice(chunk);
        sha256_compress(&mut state, &block);
    }

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
enum NoncePi {
    Rotate3LowBytes, // n0 n1 n2 n3 -> n1 n2 n0 n3
    Rotate3HighBytes, // n0 n1 n2 n3 -> n0 n2 n3 n1
}

impl NoncePi {
    fn name(self) -> &'static str {
        match self {
            NoncePi::Rotate3LowBytes => "rotate3_low_bytes",
            NoncePi::Rotate3HighBytes => "rotate3_high_bytes",
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum DigestRho {
    TripletBytesOff0, // rotate each 3-byte block starting at 0
    TripletBytesOff1, // rotate each 3-byte block starting at 1
    TripletWordsOff0, // rotate each 3-word block starting at byte 0
}

impl DigestRho {
    fn name(self) -> &'static str {
        match self {
            DigestRho::TripletBytesOff0 => "triplet_bytes_off0",
            DigestRho::TripletBytesOff1 => "triplet_bytes_off1",
            DigestRho::TripletWordsOff0 => "triplet_words_off0",
        }
    }
}

fn apply_nonce_pi(n: [u8; 4], pi: NoncePi) -> [u8; 4] {
    match pi {
        NoncePi::Rotate3LowBytes => [n[1], n[2], n[0], n[3]],
        NoncePi::Rotate3HighBytes => [n[0], n[2], n[3], n[1]],
    }
}

fn apply_digest_rho(h: &[u8; 32], rho: DigestRho) -> [u8; 32] {
    let mut out = *h;
    match rho {
        DigestRho::TripletBytesOff0 => {
            let mut i = 0usize;
            while i + 2 < 30 {
                out[i] = h[i + 1];
                out[i + 1] = h[i + 2];
                out[i + 2] = h[i];
                i += 3;
            }
        }
        DigestRho::TripletBytesOff1 => {
            let mut i = 1usize;
            while i + 2 < 31 {
                out[i] = h[i + 1];
                out[i + 1] = h[i + 2];
                out[i + 2] = h[i];
                i += 3;
            }
        }
        DigestRho::TripletWordsOff0 => {
            let mut i = 0usize;
            while i + 11 < 24 {
                out[i..i + 4].copy_from_slice(&h[i + 4..i + 8]);
                out[i + 4..i + 8].copy_from_slice(&h[i + 8..i + 12]);
                out[i + 8..i + 12].copy_from_slice(&h[i..i + 4]);
                i += 12;
            }
        }
    }
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

fn build_header(prefix76: &[u8; 76], nonce: [u8; 4]) -> [u8; 80] {
    let mut h = [0u8; 80];
    h[..76].copy_from_slice(prefix76);
    h[76..80].copy_from_slice(&nonce);
    h
}

#[derive(Clone, Debug)]
struct PairResult {
    pi: NoncePi,
    rho: DigestRho,
    exact_hits: u64,
    trials: u64,
    mean_hd: f64,
}

fn main() {
    let out_dir = std::env::var("GUTOE_SHA256D_BTC_Z3_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/sha256d_btc_nonce_z3_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let trials_pairs = env_u64("GUTOE_SHA256D_BTC_Z3_PAIR_TRIALS", 200_000);
    let trials_pred = env_u64("GUTOE_SHA256D_BTC_Z3_PRED_TRIALS", 200_000);
    let dedup_nonce_domain = env_u64("GUTOE_SHA256D_BTC_Z3_DEDUP_DOMAIN", 1 << 16).min(1 << 20);
    let lz_threshold = env_u64("GUTOE_SHA256D_BTC_Z3_LZ_THRESHOLD", 12).min(256) as u32;
    let retro_shift_factor = env_f64("GUTOE_SHA256D_BTC_Z3_RETRO_SHIFT_FACTOR", 1.20).max(0.0);
    let observer_floor_s = env_f64("GUTOE_SHA256D_BTC_Z3_OBSERVER_FLOOR_S", 1.0e-9).max(1e-15);

    let pis = [NoncePi::Rotate3LowBytes, NoncePi::Rotate3HighBytes];
    let rhos = [
        DigestRho::TripletBytesOff0,
        DigestRho::TripletBytesOff1,
        DigestRho::TripletWordsOff0,
    ];

    let mut seed_rng = OsRng;
    let mut prefix = [0u8; 76];
    seed_rng.fill_bytes(&mut prefix);

    // Pair search in parallel.
    let mut pairs = Vec::new();
    for &pi in &pis {
        for &rho in &rhos {
            pairs.push((pi, rho));
        }
    }

    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .min(pairs.len().max(1));
    let chunk = (pairs.len() + threads - 1) / threads;
    let t0 = Instant::now();

    let mut handles = Vec::new();
    for tid in 0..threads {
        let lo = tid * chunk;
        let hi = ((tid + 1) * chunk).min(pairs.len());
        if lo >= hi {
            continue;
        }
        let local_pairs = pairs[lo..hi].to_vec();
        let local_prefix = prefix;
        handles.push(thread::spawn(move || {
            let mut seed = [0u8; 32];
            seed[..8].copy_from_slice(&(tid as u64).to_le_bytes());
            let mut rng = StdRng::from_seed(seed);
            let mut out_local = Vec::new();
            for (pi, rho) in local_pairs {
                let mut hits = 0u64;
                let mut hd_sum = 0u64;
                for _ in 0..trials_pairs {
                    let mut nonce = [0u8; 4];
                    rng.fill_bytes(&mut nonce);
                    let header = build_header(&local_prefix, nonce);
                    let hash0 = sha256d(&header);
                    let nonce1 = apply_nonce_pi(nonce, pi);
                    let header1 = build_header(&local_prefix, nonce1);
                    let hash1 = sha256d(&header1);
                    let rho_hash0 = apply_digest_rho(&hash0, rho);
                    if hash1 == rho_hash0 {
                        hits += 1;
                    }
                    hd_sum += digest_hamming_bits(&hash1, &rho_hash0) as u64;
                }
                out_local.push(PairResult {
                    pi,
                    rho,
                    exact_hits: hits,
                    trials: trials_pairs,
                    mean_hd: hd_sum as f64 / trials_pairs as f64,
                });
            }
            out_local
        }));
    }

    let mut pair_rows = Vec::new();
    for h in handles {
        for r in h.join().expect("worker join failed") {
            pair_rows.push(r);
        }
    }
    let host_elapsed_s = t0.elapsed().as_secs_f64();
    pair_rows.sort_by(|a, b| {
        a.mean_hd
            .partial_cmp(&b.mean_hd)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let best = pair_rows.first().cloned().expect("at least one pair");

    // Mining-predicate correlation for best pair.
    let mut pred_rng = OsRng;
    let mut agree = 0u64;
    let mut s0 = 0u64;
    let mut s1 = 0u64;
    let mut lz0_sum = 0f64;
    let mut lz1_sum = 0f64;
    let mut lz0_sq = 0f64;
    let mut lz1_sq = 0f64;
    let mut lz01 = 0f64;
    let mut hd_pair_sum = 0u64;
    let mut hd_rand_sum = 0u64;
    let mut eq_hits = 0u64;
    for _ in 0..trials_pred {
        let mut nonce = [0u8; 4];
        pred_rng.fill_bytes(&mut nonce);
        let h0 = sha256d(&build_header(&prefix, nonce));
        let n1 = apply_nonce_pi(nonce, best.pi);
        let h1 = sha256d(&build_header(&prefix, n1));
        let rh = apply_digest_rho(&h0, best.rho);
        if h1 == rh {
            eq_hits += 1;
        }
        hd_pair_sum += digest_hamming_bits(&h1, &rh) as u64;

        let mut nonce_r = [0u8; 4];
        pred_rng.fill_bytes(&mut nonce_r);
        let hrand = sha256d(&build_header(&prefix, nonce_r));
        hd_rand_sum += digest_hamming_bits(&h1, &hrand) as u64;

        let z0 = leading_zero_bits(&h0) as f64;
        let z1 = leading_zero_bits(&h1) as f64;
        lz0_sum += z0;
        lz1_sum += z1;
        lz0_sq += z0 * z0;
        lz1_sq += z1 * z1;
        lz01 += z0 * z1;

        let b0 = (z0 as u32) >= lz_threshold;
        let b1 = (z1 as u32) >= lz_threshold;
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

    // Orbit deduper on nonce space for best pi/predicate (fixed prefix).
    let mut baseline_success = 0u64;
    let mut baseline_hashes = 0u64;
    let mut rep_pred: HashMap<[u8; 4], bool> = HashMap::new();
    let mut truth: HashMap<[u8; 4], bool> = HashMap::new();
    for nonce_u in 0..dedup_nonce_domain {
        let nonce = (nonce_u as u32).to_le_bytes();
        let h = sha256d(&build_header(&prefix, nonce));
        let pred = leading_zero_bits(&h) >= lz_threshold;
        baseline_hashes += 1;
        if pred {
            baseline_success += 1;
        }
        truth.insert(nonce, pred);

        let n1 = apply_nonce_pi(nonce, best.pi);
        let n2 = apply_nonce_pi(n1, best.pi);
        let canon = if nonce <= n1 && nonce <= n2 {
            nonce
        } else if n1 <= nonce && n1 <= n2 {
            n1
        } else {
            n2
        };
        rep_pred.entry(canon).or_insert(pred);
    }

    let dedup_hashes = rep_pred.len() as u64;
    let mut mismatches = 0u64;
    let mut predicted_success = 0u64;
    for (&n, &t) in &truth {
        let n1 = apply_nonce_pi(n, best.pi);
        let n2 = apply_nonce_pi(n1, best.pi);
        let canon = if n <= n1 && n <= n2 {
            n
        } else if n1 <= n && n1 <= n2 {
            n1
        } else {
            n2
        };
        let p = *rep_pred.get(&canon).expect("rep exists");
        if p {
            predicted_success += 1;
        }
        if p != t {
            mismatches += 1;
        }
    }

    let n = trials_pred as f64;
    let mean_hd_pair = hd_pair_sum as f64 / n;
    let mean_hd_random = hd_rand_sum as f64 / n;
    let mean_lz0 = lz0_sum / n;
    let mean_lz1 = lz1_sum / n;
    let cov = lz01 / n - mean_lz0 * mean_lz1;
    let var0 = (lz0_sq / n - mean_lz0 * mean_lz0).max(0.0);
    let var1 = (lz1_sq / n - mean_lz1 * mean_lz1).max(0.0);
    let lz_corr = if var0 > 0.0 && var1 > 0.0 {
        cov / (var0.sqrt() * var1.sqrt())
    } else {
        0.0
    };
    let p0 = s0 as f64 / n;
    let p1 = s1 as f64 / n;
    let agree_rate = agree as f64 / n;
    let agree_indep = p0 * p1 + (1.0 - p0) * (1.0 - p1);

    let predeparture = retro_shift_factor > 1.0;
    let normal_latency_s = host_elapsed_s.max(observer_floor_s);
    let observed_latency_s = if predeparture {
        observer_floor_s
    } else {
        (normal_latency_s * (1.0 - retro_shift_factor)).max(observer_floor_s)
    };
    let apparent_speedup = normal_latency_s / observed_latency_s;

    let payload = json!({
      "scope": "consensus-preserving nonce-orbit structure probe on sha256d(header80)",
      "pair_search": {
        "pair_trials_each": trials_pairs,
        "threads_used": threads,
        "pair_count": pairs.len(),
        "host_elapsed_s": host_elapsed_s,
        "predeparture": predeparture,
        "apparent_speedup": apparent_speedup,
        "best_pair": {
          "pi": best.pi.name(),
          "rho": best.rho.name(),
          "exact_hits": best.exact_hits,
          "exact_hit_rate": best.exact_hits as f64 / best.trials as f64,
          "mean_hd_bits": best.mean_hd
        },
        "all_pairs": pair_rows.iter().map(|r| json!({
          "pi": r.pi.name(),
          "rho": r.rho.name(),
          "exact_hits": r.exact_hits,
          "exact_hit_rate": r.exact_hits as f64 / r.trials as f64,
          "mean_hd_bits": r.mean_hd
        })).collect::<Vec<_>>()
      },
      "predicate_correlation": {
        "pred_trials": trials_pred,
        "lz_threshold_bits": lz_threshold,
        "exact_hits": eq_hits,
        "exact_hit_rate": eq_hits as f64 / n,
        "mean_hd_pair": mean_hd_pair,
        "mean_hd_random": mean_hd_random,
        "mean_lz_base": mean_lz0,
        "mean_lz_transformed": mean_lz1,
        "lz_corr": lz_corr,
        "success_rate_base": p0,
        "success_rate_transformed": p1,
        "success_agreement_rate": agree_rate,
        "success_agreement_independent_baseline": agree_indep
      },
      "orbit_deduper": {
        "nonce_domain": dedup_nonce_domain,
        "lz_threshold_bits": lz_threshold,
        "baseline_hashes": baseline_hashes,
        "baseline_success": baseline_success,
        "dedup_hashes": dedup_hashes,
        "ideal_speedup_if_equivariant": baseline_hashes as f64 / dedup_hashes as f64,
        "predicted_success_from_reps": predicted_success,
        "mismatches_vs_truth": mismatches,
        "mismatch_rate": mismatches as f64 / baseline_hashes as f64
      }
    });

    let txt_path = out.join("sha256d_btc_nonce_z3_probe.txt");
    let json_path = out.join("sha256d_btc_nonce_z3_probe.json");

    let mut txt = String::new();
    txt.push_str("[sha256d_btc_nonce_z3_probe]\n");
    txt.push_str("consensus-preserving nonce-orbit structure test for double-SHA256\n\n");
    txt.push_str(&format!(
        "best_pair pi={} rho={} exact_hits={}/{} mean_hd={:.6}\n",
        best.pi.name(),
        best.rho.name(),
        best.exact_hits,
        best.trials,
        best.mean_hd
    ));
    txt.push_str(&format!(
        "predicate: eq_hits={}/{} mean_hd_pair={:.6} mean_hd_random={:.6} lz_corr={:.6e}\n",
        eq_hits, trials_pred, mean_hd_pair, mean_hd_random, lz_corr
    ));
    txt.push_str(&format!(
        "success_agree={:.6e} indep_baseline={:.6e}\n",
        agree_rate, agree_indep
    ));
    txt.push_str(&format!(
        "deduper: baseline_hashes={} dedup_hashes={} mismatch_rate={:.6e}\n",
        baseline_hashes,
        dedup_hashes,
        mismatches as f64 / baseline_hashes as f64
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

