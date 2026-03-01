//! Commuting Z3-pair probe on consensus-preserving Bitcoin-style nonce lane.
//!
//! Goal:
//! - enumerate all nontrivial order-3 permutations of 4 nonce bytes,
//! - find commuting pairs,
//! - test whether any commuting pair can support useful joint equivariance
//!   for H2 = SHA256(SHA256(header80)).
//!
//! This directly tests the "recursive 3^k elimination" hinge for k=2
//! in the strict nonce-only action family.

use rand::rngs::OsRng;
use rand::RngCore;
use serde_json::json;
use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Perm4 {
    idx: [usize; 4],
}

impl Perm4 {
    fn identity() -> Self {
        Self { idx: [0, 1, 2, 3] }
    }

    fn apply(self, n: [u8; 4]) -> [u8; 4] {
        [n[self.idx[0]], n[self.idx[1]], n[self.idx[2]], n[self.idx[3]]]
    }

    fn compose(self, other: Self) -> Self {
        // self ∘ other
        let o = other.idx;
        let s = self.idx;
        Self {
            idx: [o[s[0]], o[s[1]], o[s[2]], o[s[3]]],
        }
    }

    fn pow(self, k: u32) -> Self {
        let mut out = Self::identity();
        for _ in 0..k {
            out = out.compose(self);
        }
        out
    }

    fn is_order3_nontrivial(self) -> bool {
        self != Self::identity() && self.pow(3) == Self::identity()
    }

    fn commutes(self, other: Self) -> bool {
        self.compose(other) == other.compose(self)
    }

    fn cyclic_subgroup(self) -> BTreeSet<Self> {
        let mut s = BTreeSet::new();
        s.insert(Self::identity());
        s.insert(self);
        s.insert(self.pow(2));
        s
    }

    fn name(self) -> String {
        format!(
            "[{} {} {} {}]",
            self.idx[0], self.idx[1], self.idx[2], self.idx[3]
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum DigestRho {
    TripletBytesOff0,
    TripletBytesOff1,
    TripletWordsOff0,
}

impl DigestRho {
    fn all() -> [Self; 3] {
        [
            Self::TripletBytesOff0,
            Self::TripletBytesOff1,
            Self::TripletWordsOff0,
        ]
    }

    fn name(self) -> &'static str {
        match self {
            Self::TripletBytesOff0 => "triplet_bytes_off0",
            Self::TripletBytesOff1 => "triplet_bytes_off1",
            Self::TripletWordsOff0 => "triplet_words_off0",
        }
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

fn digest_hd(a: &[u8; 32], b: &[u8; 32]) -> u32 {
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

fn all_perms4() -> Vec<Perm4> {
    let mut out = Vec::new();
    for a in 0..4 {
        for b in 0..4 {
            if b == a {
                continue;
            }
            for c in 0..4 {
                if c == a || c == b {
                    continue;
                }
                for d in 0..4 {
                    if d == a || d == b || d == c {
                        continue;
                    }
                    out.push(Perm4 { idx: [a, b, c, d] });
                }
            }
        }
    }
    out
}

fn generated_group_size(gens: &[Perm4]) -> usize {
    let mut seen = HashSet::new();
    let mut q = VecDeque::new();
    let id = Perm4::identity();
    seen.insert(id);
    q.push_back(id);
    while let Some(cur) = q.pop_front() {
        for &g in gens {
            let nxt = cur.compose(g);
            if seen.insert(nxt) {
                q.push_back(nxt);
            }
        }
    }
    seen.len()
}

fn main() {
    let out_dir = std::env::var("GUTOE_SHA256D_BTC_COMMUTE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/sha256d_btc_commuting_z3_pair_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let trials_single = env_u64("GUTOE_SHA256D_BTC_COMMUTE_SINGLE_TRIALS", 80_000);
    let trials_joint = env_u64("GUTOE_SHA256D_BTC_COMMUTE_JOINT_TRIALS", 80_000);
    let lz_threshold = env_u64("GUTOE_SHA256D_BTC_COMMUTE_LZ_THRESHOLD", 12).min(256) as u32;

    let mut rng = OsRng;
    let mut prefix = [0u8; 76];
    rng.fill_bytes(&mut prefix);

    let order3 = all_perms4()
        .into_iter()
        .filter(|p| p.is_order3_nontrivial())
        .collect::<Vec<_>>();

    let mut commuting_pairs = Vec::new();
    let mut independent_pairs = Vec::new();
    for i in 0..order3.len() {
        for j in (i + 1)..order3.len() {
            let p = order3[i];
            let q = order3[j];
            if p.commutes(q) {
                let gp = p.cyclic_subgroup();
                let gq = q.cyclic_subgroup();
                let same_subgroup = gp == gq;
                let gsize = generated_group_size(&[p, q]);
                let independent = !same_subgroup && gsize >= 9;
                commuting_pairs.push((p, q, same_subgroup, gsize));
                if independent {
                    independent_pairs.push((p, q, gsize));
                }
            }
        }
    }

    // For each order-3 permutation, pick best single rho by mean HD.
    let mut best_rho_for_perm: HashMap<Perm4, (DigestRho, f64, u64)> = HashMap::new();
    for &p in &order3 {
        let mut best = (DigestRho::TripletBytesOff0, f64::INFINITY, 0u64);
        for rho in DigestRho::all() {
            let mut hd_sum = 0u64;
            let mut hits = 0u64;
            for _ in 0..trials_single {
                let mut nonce = [0u8; 4];
                rng.fill_bytes(&mut nonce);
                let h0 = sha256d(&build_header(&prefix, nonce));
                let h1 = sha256d(&build_header(&prefix, p.apply(nonce)));
                let rh = apply_digest_rho(&h0, rho);
                if h1 == rh {
                    hits += 1;
                }
                hd_sum += digest_hd(&h1, &rh) as u64;
            }
            let mhd = hd_sum as f64 / trials_single as f64;
            if mhd < best.1 {
                best = (rho, mhd, hits);
            }
        }
        best_rho_for_perm.insert(p, best);
    }

    // Joint test for commuting pairs using composed best rhos.
    let mut joint_rows = Vec::new();
    for &(p, q, same_subgroup, gsize) in &commuting_pairs {
        let rho_p = best_rho_for_perm[&p].0;
        let rho_q = best_rho_for_perm[&q].0;
        let mut eq_hits = 0u64;
        let mut hd_sum = 0u64;
        let mut agree = 0u64;
        let mut s0 = 0u64;
        let mut s1 = 0u64;
        for _ in 0..trials_joint {
            let mut nonce = [0u8; 4];
            rng.fill_bytes(&mut nonce);
            let h0 = sha256d(&build_header(&prefix, nonce));
            let n2 = p.apply(q.apply(nonce));
            let h2 = sha256d(&build_header(&prefix, n2));
            let rhs = apply_digest_rho(&apply_digest_rho(&h0, rho_q), rho_p);
            if h2 == rhs {
                eq_hits += 1;
            }
            hd_sum += digest_hd(&h2, &rhs) as u64;
            let b0 = leading_zero_bits(&h0) >= lz_threshold;
            let b1 = leading_zero_bits(&h2) >= lz_threshold;
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
        let n = trials_joint as f64;
        let p0 = s0 as f64 / n;
        let p1 = s1 as f64 / n;
        let agree_rate = agree as f64 / n;
        let indep = p0 * p1 + (1.0 - p0) * (1.0 - p1);
        joint_rows.push(json!({
          "p": p.name(),
          "q": q.name(),
          "same_subgroup": same_subgroup,
          "generated_group_size": gsize,
          "rho_p": rho_p.name(),
          "rho_q": rho_q.name(),
          "eq_hits": eq_hits,
          "eq_hit_rate": eq_hits as f64 / n,
          "mean_hd_bits": hd_sum as f64 / n,
          "success_agreement_rate": agree_rate,
          "success_independent_baseline": indep
        }));
    }

    let best_joint = joint_rows
        .iter()
        .min_by(|a, b| {
            a["mean_hd_bits"]
                .as_f64()
                .unwrap_or(f64::INFINITY)
                .partial_cmp(&b["mean_hd_bits"].as_f64().unwrap_or(f64::INFINITY))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .cloned();

    let payload = json!({
      "scope": "commuting-pair test for recursive 3^k hypothesis on consensus nonce lane",
      "order3_nonce_permutations_count": order3.len(),
      "order3_nonce_permutations": order3.iter().map(|p| p.name()).collect::<Vec<_>>(),
      "commuting_pairs_count": commuting_pairs.len(),
      "commuting_pairs_independent_count": independent_pairs.len(),
      "independent_pairs": independent_pairs.iter().map(|(p,q,g)| json!({
        "p": p.name(),
        "q": q.name(),
        "generated_group_size": g
      })).collect::<Vec<_>>(),
      "best_single_rho_per_perm": order3.iter().map(|p| {
        let (rho,mhd,hits) = best_rho_for_perm[p];
        json!({
          "perm": p.name(),
          "rho": rho.name(),
          "single_mean_hd_bits": mhd,
          "single_exact_hits": hits,
          "single_trials": trials_single
        })
      }).collect::<Vec<_>>(),
      "joint_rows": joint_rows,
      "best_joint": best_joint,
      "trials_single": trials_single,
      "trials_joint": trials_joint,
      "lz_threshold": lz_threshold
    });

    let txt_path = out.join("sha256d_btc_commuting_z3_pair_probe.txt");
    let json_path = out.join("sha256d_btc_commuting_z3_pair_probe.json");

    let mut txt = String::new();
    txt.push_str("[sha256d_btc_commuting_z3_pair_probe]\n");
    txt.push_str("recursive 3^k hinge test: two commuting nontrivial order-3 nonce actions\n\n");
    txt.push_str(&format!(
        "order3_nonce_permutations_count={}\n",
        order3.len()
    ));
    txt.push_str(&format!("commuting_pairs_count={}\n", commuting_pairs.len()));
    txt.push_str(&format!(
        "commuting_pairs_independent_count={}\n",
        independent_pairs.len()
    ));
    if let Some(best) = &payload["best_joint"].as_object() {
        txt.push_str("\nbest_joint_pair:\n");
        txt.push_str(&format!("p={}\n", best["p"].as_str().unwrap_or("")));
        txt.push_str(&format!("q={}\n", best["q"].as_str().unwrap_or("")));
        txt.push_str(&format!(
            "generated_group_size={}\n",
            best["generated_group_size"].as_u64().unwrap_or_default()
        ));
        txt.push_str(&format!(
            "eq_hits={}\n",
            best["eq_hits"].as_u64().unwrap_or_default()
        ));
        txt.push_str(&format!(
            "eq_hit_rate={:.6e}\n",
            best["eq_hit_rate"].as_f64().unwrap_or_default()
        ));
        txt.push_str(&format!(
            "mean_hd_bits={:.6}\n",
            best["mean_hd_bits"].as_f64().unwrap_or_default()
        ));
        txt.push_str(&format!(
            "success_agreement_rate={:.6e}\n",
            best["success_agreement_rate"].as_f64().unwrap_or_default()
        ));
        txt.push_str(&format!(
            "success_independent_baseline={:.6e}\n",
            best["success_independent_baseline"]
                .as_f64()
                .unwrap_or_default()
        ));
    }

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
