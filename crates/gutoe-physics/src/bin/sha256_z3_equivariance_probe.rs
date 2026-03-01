//! SHA-256 Z3-equivariance probe.
//!
//! Tests empirical candidates of the form:
//!   H(pi(m)) == rho(H(m))
//! where pi and rho are explicit order-3 transforms.
//!
//! This is a structural probe, not a proof of impossibility.

use rand::rngs::OsRng;
use rand::RngCore;
use serde_json::json;
use std::fs;
use std::path::PathBuf;

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

/// SHA-256 of a message up to 55 bytes (single block).
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

fn msg_pi_triplet_bytes(m: &[u8; 48]) -> [u8; 48] {
    let mut out = *m;
    for i in (0..48).step_by(3) {
        out[i] = m[i + 1];
        out[i + 1] = m[i + 2];
        out[i + 2] = m[i];
    }
    out
}

fn msg_pi_triplet_words(m: &[u8; 48]) -> [u8; 48] {
    // 12 words (4 bytes each) -> 4 groups of 3 words, each rotated.
    let mut out = [0u8; 48];
    for g in 0..4 {
        let base = g * 12;
        // word0 <- word1, word1 <- word2, word2 <- word0
        out[base..base + 4].copy_from_slice(&m[base + 4..base + 8]);
        out[base + 4..base + 8].copy_from_slice(&m[base + 8..base + 12]);
        out[base + 8..base + 12].copy_from_slice(&m[base..base + 4]);
    }
    out
}

fn rho_triplet_bytes_32(h: &[u8; 32]) -> [u8; 32] {
    // 10 triplets rotated (30 bytes), last 2 fixed.
    let mut out = *h;
    for i in (0..30).step_by(3) {
        out[i] = h[i + 1];
        out[i + 1] = h[i + 2];
        out[i + 2] = h[i];
    }
    out[30] = h[30];
    out[31] = h[31];
    out
}

fn rho_triplet_words_32(h: &[u8; 32]) -> [u8; 32] {
    // 8 words -> two 3-word rotations plus two fixed words.
    let mut out = [0u8; 32];
    // words 0,1,2 rotate
    out[0..4].copy_from_slice(&h[4..8]);
    out[4..8].copy_from_slice(&h[8..12]);
    out[8..12].copy_from_slice(&h[0..4]);
    // words 3,4,5 rotate
    out[12..16].copy_from_slice(&h[16..20]);
    out[16..20].copy_from_slice(&h[20..24]);
    out[20..24].copy_from_slice(&h[12..16]);
    // words 6,7 fixed
    out[24..28].copy_from_slice(&h[24..28]);
    out[28..32].copy_from_slice(&h[28..32]);
    out
}

fn apply_pi(name: &str, m: &[u8; 48]) -> [u8; 48] {
    match name {
        "triplet_bytes" => msg_pi_triplet_bytes(m),
        "triplet_words" => msg_pi_triplet_words(m),
        _ => unreachable!("unknown pi"),
    }
}

fn apply_rho(name: &str, h: &[u8; 32]) -> [u8; 32] {
    match name {
        "triplet_bytes" => rho_triplet_bytes_32(h),
        "triplet_words" => rho_triplet_words_32(h),
        _ => unreachable!("unknown rho"),
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_SHA256_Z3_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/sha256_z3_equivariance_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let trials = env_u64("GUTOE_SHA256_Z3_TRIALS", 20_000).max(1);

    let pi_candidates = ["triplet_bytes", "triplet_words"];
    let rho_candidates = ["triplet_bytes", "triplet_words"];

    // Verify each transform is order-3 on a sample.
    let sample: [u8; 48] = core::array::from_fn(|i| i as u8);
    let sample_h = sha256(&sample);
    let mut order3_checks = Vec::new();
    for &pi in &pi_candidates {
        let p1 = apply_pi(pi, &sample);
        let p2 = apply_pi(pi, &p1);
        let p3 = apply_pi(pi, &p2);
        order3_checks.push(json!({
          "domain": "message",
          "name": pi,
          "order3_identity_on_sample": p3 == sample
        }));
    }
    for &rho in &rho_candidates {
        let r1 = apply_rho(rho, &sample_h);
        let r2 = apply_rho(rho, &r1);
        let r3 = apply_rho(rho, &r2);
        order3_checks.push(json!({
          "domain": "digest",
          "name": rho,
          "order3_identity_on_sample": r3 == sample_h
        }));
    }

    let mut rng = OsRng;
    let mut rows = Vec::new();
    let mut report_txt = String::new();
    report_txt.push_str("[sha256_z3_equivariance_probe]\n");
    report_txt.push_str("tests H(pi(m)) == rho(H(m)) for explicit order-3 transforms\n\n");
    report_txt.push_str(&format!("trials = {}\n", trials));
    report_txt.push_str("message_len = 48 (single-block SHA256 lane)\n\n");

    for &pi in &pi_candidates {
        for &rho in &rho_candidates {
            let mut hits = 0_u64;
            for _ in 0..trials {
                let mut m = [0u8; 48];
                rng.fill_bytes(&mut m);
                let hm = sha256(&m);
                let pm = apply_pi(pi, &m);
                let hpm = sha256(&pm);
                let rhm = apply_rho(rho, &hm);
                if hpm == rhm {
                    hits += 1;
                }
            }
            let hit_rate = hits as f64 / trials as f64;
            rows.push(json!({
              "pi": pi,
              "rho": rho,
              "hits": hits,
              "trials": trials,
              "hit_rate": hit_rate
            }));
            report_txt.push_str(&format!(
                "pair pi={} rho={}: hits={} / {} (rate {:.6e})\n",
                pi, rho, hits, trials, hit_rate
            ));
        }
    }

    report_txt.push_str("\nrandom_baseline_for_exact_256bit_match ~= 2^-256 (~8.64e-78)\n");
    report_txt.push_str(
        "zero hits at practical trial counts indicates no detectable nontrivial Z3-equivariance for tested transforms.\n",
    );

    let payload = json!({
      "trials": trials,
      "message_len_bytes": 48,
      "order3_checks": order3_checks,
      "results": rows,
      "baseline_exact_match_probability": "2^-256",
      "scope": "empirical candidate test only; not a proof of impossibility"
    });

    let txt_path = out.join("sha256_z3_equivariance_probe.txt");
    let json_path = out.join("sha256_z3_equivariance_probe.json");

    fs::write(&txt_path, report_txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}

