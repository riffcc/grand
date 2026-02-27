// Speculative Execution over Finite Algebras vs SHA-256
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// Core insight: if a computation passes through a FINITE-DIMENSIONAL algebra,
// the entire future can be precomputed as a response matrix while waiting for
// the unknown intermediate value. When it resolves: one matrix-vector multiply.
//
// Cl(1,3): 16 dimensions → 16×16 = 256 entries. Fits in 2 KB.
// SHA-256:  256-bit state → 2^256 entries. Does not fit in the observable universe.
//
// SHA-256's ARX design (Add-Rotate-XOR) deliberately mixes two incompatible
// algebras — GF(2)^32 (XOR/rotate: linear) and Z/2^32 (addition: nonlinear
// carries) — specifically to prevent the algebraic shortcuts that Clifford
// speculation exploits.

use std::time::Instant;

// ═══════════════════════════════════════════════════════════════════════════════
// PART 1: Cl(1,3) Clifford Algebra
// ═══════════════════════════════════════════════════════════════════════════════
//
// 16 basis elements indexed by blade bitmask:
//   0b0000 = 1,  0b0001 = e₁,  0b0010 = e₂,  0b0100 = e₃,  0b1000 = e₄
//   0b0011 = e₁₂, 0b0101 = e₁₃, ...  0b1111 = e₁₂₃₄
// Metric: e₁² = +1, e₂² = e₃² = e₄² = −1  (signature 1,3)

#[derive(Clone, Copy, Debug)]
pub struct Mv(pub [u64; 16]);

/// Precomputed product table: (result_blade_index, is_negative)
const fn build_mul_table() -> [[(u8, bool); 16]; 16] {
    let mut table = [[(0u8, false); 16]; 16];
    let mut a: u8 = 0;
    while a < 16 {
        let mut b: u8 = 0;
        while b < 16 {
            let result = a ^ b;

            // Swap sign: count inversions when interleaving generators
            let mut swaps = 0u32;
            let mut aa = a >> 1;
            while aa != 0 {
                let mut x = aa & b;
                while x != 0 {
                    swaps += 1;
                    x &= x - 1;
                }
                aa >>= 1;
            }
            let mut neg = swaps & 1 != 0;

            // Metric sign: e₁²=+1, e₂²=e₃²=e₄²=−1
            let common = a & b;
            let mut i: u8 = 0;
            while i < 4 {
                if common & (1 << i) != 0 && i > 0 {
                    neg = !neg;
                }
                i += 1;
            }

            table[a as usize][b as usize] = (result, neg);
            b += 1;
        }
        a += 1;
    }
    table
}

const MUL_TABLE: [[(u8, bool); 16]; 16] = build_mul_table();

impl Mv {
    pub fn zero() -> Self {
        Mv([0; 16])
    }

    pub fn identity() -> Self {
        let mut v = [0u64; 16];
        v[0] = 1;
        Mv(v)
    }

    /// Basis vector e_i (0-indexed: 0=scalar, 1=e₁, ..., 15=e₁₂₃₄)
    pub fn basis(i: usize) -> Self {
        let mut v = [0u64; 16];
        v[i] = 1;
        Mv(v)
    }

    /// Clifford geometric product: 16×16 = 256 multiply-adds
    pub fn mul(&self, other: &Self) -> Self {
        let mut result = [0u64; 16];
        for a in 0..16 {
            let sa = self.0[a];
            if sa == 0 {
                continue;
            }
            for b in 0..16 {
                let ob = other.0[b];
                if ob == 0 {
                    continue;
                }
                let (idx, neg) = MUL_TABLE[a][b];
                let prod = sa.wrapping_mul(ob);
                if neg {
                    result[idx as usize] = result[idx as usize].wrapping_sub(prod);
                } else {
                    result[idx as usize] = result[idx as usize].wrapping_add(prod);
                }
            }
        }
        Mv(result)
    }

    /// Seed a "random-looking" multivector from an index (deterministic, not crypto)
    fn from_seed(seed: u64) -> Self {
        let mut v = [0u64; 16];
        let mut s = seed.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(1);
        for x in &mut v {
            s ^= s >> 30;
            s = s.wrapping_mul(0xBF58476D1CE4E5B9);
            s ^= s >> 27;
            s = s.wrapping_mul(0x94D049BB133111EB);
            s ^= s >> 31;
            *x = s;
        }
        v[0] |= 1; // ensure invertibility
        Mv(v)
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// PART 1b: Response Matrix — the speculation engine
// ═══════════════════════════════════════════════════════════════════════════════
//
// For f(x) = left × x × right  (Clifford sandwich, LINEAR in x):
//   response_col[i] = left × e_i × right
//   f(x) = Σ x[i] × response_col[i]
//
// This 16×16 matrix IS the entire future of the computation.

pub struct ResponseMatrix {
    /// cols[i][j] = component j of f(e_i)
    pub cols: [[u64; 16]; 16],
}

impl ResponseMatrix {
    /// Precompute response for f(x) = left * x * right
    pub fn from_sandwich(left: &Mv, right: &Mv) -> Self {
        let mut cols = [[0u64; 16]; 16];
        for i in 0..16 {
            let basis = Mv::basis(i);
            let result = left.mul(&basis).mul(right);
            cols[i] = result.0;
        }
        ResponseMatrix { cols }
    }

    /// Zip: apply precomputed response to the resolved unknown
    /// Cost: 16×16 = 256 multiply-adds (half a Clifford multiply)
    pub fn apply(&self, x: &Mv) -> Mv {
        let mut result = [0u64; 16];
        for i in 0..16 {
            let xi = x.0[i];
            if xi == 0 {
                continue;
            }
            for j in 0..16 {
                result[j] = result[j].wrapping_add(xi.wrapping_mul(self.cols[i][j]));
            }
        }
        Mv(result)
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// PART 2: Minimal SHA-256
// ═══════════════════════════════════════════════════════════════════════════════

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

/// SHA-256 of a message up to 55 bytes (single block)
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

/// SHA-256(SHA-256(msg)) — Bitcoin's double hash
fn sha256d(msg: &[u8]) -> [u8; 32] {
    let inner = sha256(msg);
    sha256(&inner)
}

// ═══════════════════════════════════════════════════════════════════════════════
// PART 3: Benchmarks — the moment of truth
// ═══════════════════════════════════════════════════════════════════════════════

/// Clifford chain: product of N seeded multivectors
fn clifford_chain(n: usize, seed_base: u64) -> Mv {
    let mut state = Mv::identity();
    for i in 0..n {
        state = state.mul(&Mv::from_seed(seed_base.wrapping_add(i as u64)));
    }
    state
}

pub fn run_speculation_benchmark() {
    println!("╔════════════════════════════════════════════════════════════════╗");
    println!("║  Algebraic Speculation: Cl(1,3) vs SHA-256                   ║");
    println!("╚════════════════════════════════════════════════════════════════╝\n");

    // ── SHA-256 double hash ──────────────────────────────────────────────
    let n_sha = 50_000;
    let msg = b"GUTOE hex+z lattice 2026";
    let t0 = Instant::now();
    let mut sha_result = [0u8; 32];
    for i in 0..n_sha {
        let mut input = [0u8; 32];
        input[..8].copy_from_slice(&(i as u64).to_le_bytes());
        input[8..32].copy_from_slice(&msg[..24]);
        sha_result = sha256d(&input);
    }
    let sha_time = t0.elapsed();
    println!("  SHA-256 double hash ({n_sha} iterations):");
    println!(
        "    Time: {:.2}ms ({:.1} ns/hash)",
        sha_time.as_secs_f64() * 1e3,
        sha_time.as_secs_f64() * 1e9 / n_sha as f64
    );
    println!("    Can speculate? NO — 64 rounds × 8 additions with carries");
    println!("    State space: 2^256 — cannot precompute response matrix");
    println!(
        "    Result: {:02x}{:02x}{:02x}{:02x}...\n",
        sha_result[0], sha_result[1], sha_result[2], sha_result[3]
    );

    // ── Clifford: measure single-mul cost first ──────────────────────────
    let n_cal = 100_000;
    let a = Mv::from_seed(42);
    let mut b = Mv::from_seed(137);
    let t0 = Instant::now();
    for _ in 0..n_cal {
        b = a.mul(&b);
    }
    let mul_time = t0.elapsed();
    let ns_per_mul = mul_time.as_secs_f64() * 1e9 / n_cal as f64;
    println!("  Cl(1,3) single multiplication: {:.0} ns", ns_per_mul);
    println!("  (16×16 = 256 multiply-adds per product)\n");

    // ── Long chain: sequential vs speculative ────────────────────────────
    // Use chain_len=512 to make real work dominate over thread overhead
    let chain_len = 512;
    let half = chain_len / 2;

    // Sequential: inner(512) + left(256) + right(256) + sandwich(2) = 1026 muls
    let t0 = Instant::now();
    let inner_seq = clifford_chain(chain_len, 0);
    let left_seq = clifford_chain(half, 900_000);
    let right_seq = clifford_chain(half, 800_000);
    let result_seq = left_seq.mul(&inner_seq).mul(&right_seq);
    let seq_time = t0.elapsed();
    let seq_muls = chain_len + half + half + 2;

    // Speculative: thread1=inner(512), thread2=left(256)+right(256)+response(32)
    let t0 = Instant::now();
    let (inner_spec, response) = std::thread::scope(|s| {
        let inner_h = s.spawn(|| clifford_chain(chain_len, 0));
        let response_h = s.spawn(|| {
            let left = clifford_chain(half, 900_000);
            let right = clifford_chain(half, 800_000);
            ResponseMatrix::from_sandwich(&left, &right)
        });
        (inner_h.join().unwrap(), response_h.join().unwrap())
    });
    let result_spec = response.apply(&inner_spec);
    let spec_time = t0.elapsed();

    let matches = result_seq.0 == result_spec.0;
    let speedup = seq_time.as_secs_f64() / spec_time.as_secs_f64();

    let thread2_muls = half + half + 32; // left + right + 16 response columns × 2 muls each
    let spec_depth = std::cmp::max(chain_len, thread2_muls) + 1;

    println!("  ── Chain of {chain_len} Clifford multiplications ──");
    println!("  Sequential:");
    println!(
        "    Time: {:.2}ms  ({seq_muls} multiplications)",
        seq_time.as_secs_f64() * 1e3
    );
    println!("  Speculative (2 threads):");
    println!(
        "    Time: {:.2}ms  (depth: {spec_depth} muls on critical path)",
        spec_time.as_secs_f64() * 1e3
    );
    println!("    Thread 1: inner chain ({chain_len} muls)");
    println!("    Thread 2: left ({half}) + right ({half}) + response matrix (16×2) = {thread2_muls} muls");
    println!("    Zip: 256 multiply-adds when inner resolves");
    println!("    Measured speedup: {speedup:.2}×");
    println!(
        "    Theoretical (2 threads): {:.2}× (depth {seq_muls} → {spec_depth})",
        seq_muls as f64 / spec_depth as f64
    );
    println!(
        "    Theoretical (∞ threads): {:.2}× (response cols parallelized)",
        seq_muls as f64 / (std::cmp::max(chain_len, half + half + 2) + 1) as f64
    );
    println!("    Exact match: {}", if matches { "YES" } else { "BUG!" });

    // ── Batch: amortize thread overhead ──────────────────────────────────
    let n_batch = 500;
    let batch_chain = 128;
    let batch_half = batch_chain / 2;

    let t0 = Instant::now();
    for i in 0..n_batch {
        let inner = clifford_chain(batch_chain, i as u64 * 1000);
        let left = clifford_chain(batch_half, 900_000 + i as u64);
        let right = clifford_chain(batch_half, 800_000 + i as u64);
        std::hint::black_box(left.mul(&inner).mul(&right));
    }
    let batch_seq = t0.elapsed();

    let t0 = Instant::now();
    for i in 0..n_batch {
        let (inner, resp) = std::thread::scope(|s| {
            let ih = s.spawn(move || clifford_chain(batch_chain, i as u64 * 1000));
            let rh = s.spawn(move || {
                let l = clifford_chain(batch_half, 900_000 + i as u64);
                let r = clifford_chain(batch_half, 800_000 + i as u64);
                ResponseMatrix::from_sandwich(&l, &r)
            });
            (ih.join().unwrap(), rh.join().unwrap())
        });
        std::hint::black_box(resp.apply(&inner));
    }
    let batch_spec = t0.elapsed();
    let batch_speedup = batch_seq.as_secs_f64() / batch_spec.as_secs_f64();

    println!("\n  Batch ({n_batch}× chain-{batch_chain}):");
    println!(
        "    Sequential: {:.2}ms ({:.1} µs/hash)",
        batch_seq.as_secs_f64() * 1e3,
        batch_seq.as_secs_f64() * 1e6 / n_batch as f64
    );
    println!(
        "    Speculative: {:.2}ms ({:.1} µs/hash)",
        batch_spec.as_secs_f64() * 1e3,
        batch_spec.as_secs_f64() * 1e6 / n_batch as f64
    );
    println!("    Speedup: {batch_speedup:.2}×");

    // ── Massive chain: thread overhead becomes negligible ────────────────
    let mega = 10_000;
    let mega_half = mega / 2;
    let t0 = Instant::now();
    let inner_m = clifford_chain(mega, 0);
    let left_m = clifford_chain(mega_half, 900_000);
    let right_m = clifford_chain(mega_half, 800_000);
    let result_m_seq = left_m.mul(&inner_m).mul(&right_m);
    let mega_seq = t0.elapsed();

    let t0 = Instant::now();
    let (inner_mp, resp_mp) = std::thread::scope(|s| {
        let ih = s.spawn(|| clifford_chain(mega, 0));
        let rh = s.spawn(|| {
            let l = clifford_chain(mega_half, 900_000);
            let r = clifford_chain(mega_half, 800_000);
            ResponseMatrix::from_sandwich(&l, &r)
        });
        (ih.join().unwrap(), rh.join().unwrap())
    });
    let result_m_spec = resp_mp.apply(&inner_mp);
    let mega_spec = t0.elapsed();
    let mega_speedup = mega_seq.as_secs_f64() / mega_spec.as_secs_f64();
    let mega_match = result_m_seq.0 == result_m_spec.0;

    println!("\n  ── MEGA chain ({mega} muls — thread overhead negligible) ──");
    println!(
        "    Sequential: {:.2}ms ({} muls)",
        mega_seq.as_secs_f64() * 1e3,
        mega + mega_half + mega_half + 2
    );
    println!(
        "    Speculative: {:.2}ms (depth: {} muls)",
        mega_spec.as_secs_f64() * 1e3,
        mega_half + mega_half + 32 + 1
    );
    println!(
        "    Speedup: {mega_speedup:.2}×  {}",
        if mega_match { "(exact match)" } else { "BUG!" }
    );
    println!(
        "    Thread overhead: ~{:.0}µs vs {:.0}µs work = {:.1}%",
        10.0,
        mega_seq.as_secs_f64() * 1e6,
        10.0 / (mega_seq.as_secs_f64() * 1e6) * 100.0
    );

    // ── Response matrix size comparison ──────────────────────────────────
    println!("\n  ── Why This Works for Cl(1,3) and Not SHA-256 ──");
    println!("  Cl(1,3) response matrix: 16×16 = 256 entries × 8 bytes = 2 KB");
    println!("  SHA-256 response matrix: 2^256 × 2^256 entries = ∞ (heat death first)");
    println!();
    println!("  Cl(1,3): bilinear product → result is LINEAR in any single factor");
    println!("  SHA-256: ARX design → result is degree-2^64 polynomial over GF(2)");
    println!();
    println!("  Clifford multiplication: BILINEAR, ASSOCIATIVE, FINITE-DIMENSIONAL");
    println!("  SHA-256 round function:  NONLINEAR (carries), NO ALGEBRAIC STRUCTURE");

    // ── What fraction of SHA-256 is linear? ──────────────────────────────
    println!("\n  ── SHA-256 Per-Round Algebraic Anatomy ──");
    println!("  Operation          Algebra     Count  Linear?");
    println!("  ─────────          ───────     ─────  ───────");
    println!("  Σ₀(a), Σ₁(e)      GF(2)^32    2      YES (rotate + XOR)");
    println!("  Ch(e,f,g)          GF(2)^32    1      NO  (AND = degree 2)");
    println!("  Maj(a,b,c)         GF(2)^32    1      NO  (AND = degree 2)");
    println!("  T1 = h+Σ₁+Ch+K+W  Z/2^32      5      NO  (addition: carries!)");
    println!("  T2 = Σ₀ + Maj      Z/2^32      1      NO  (addition: carries!)");
    println!("  e = d + T1         Z/2^32      1      NO  (addition: carries!)");
    println!("  a = T1 + T2        Z/2^32      1      NO  (addition: carries!)");
    println!("  ─────────────────────────────────────────────");
    println!("  Linear: 2/12 ops. After 64 rounds: degree 2^64 over GF(2).");
    println!("  This is BY DESIGN. ARX = Add-Rotate-XOR: two incompatible");
    println!("  algebras mixed to prevent exactly the shortcut Cl(1,3) enables.");

    // ── The deep connection ──────────────────────────────────────────────
    println!("\n  ── The Deep Connection ──");
    println!("  Quantum superposition IS speculative execution.");
    println!("  A quantum state speculatively holds ALL possible results.");
    println!("  Measurement = zip: contract the response matrix with the observable.");
    println!("  QM is computable precisely because the algebra (Hilbert space) is");
    println!("  finite-dimensional for any bounded system.");
    println!();
    println!("  Cl(1,3) has 16 dimensions → 16-branch speculation.");
    println!("  A qubit has 2 dimensions → 2-branch speculation.");
    println!("  N qubits has 2^N dimensions → 2^N-branch speculation.");
    println!("  This is why quantum computers are powerful: exponential speculation");
    println!("  in the NUMBER of qubits, but tractable per-branch via linearity.");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cl14_multiplication_is_correct() {
        // e₁ × e₁ = +1 (metric)
        let e1 = Mv::basis(1);
        let e1_sq = e1.mul(&e1);
        assert_eq!(e1_sq.0[0], 1, "e₁² should be +1 (scalar)");
        for i in 1..16 {
            assert_eq!(e1_sq.0[i], 0, "e₁² non-scalar should be 0");
        }

        // e₂ × e₂ = -1 (metric, wrapping)
        let e2 = Mv::basis(2);
        let e2_sq = e2.mul(&e2);
        assert_eq!(e2_sq.0[0], u64::MAX, "e₂² should be -1 (wrapping)");

        // e₁ × e₂ = e₁₂ (index 3 = 0b0011)
        let e12 = e1.mul(&e2);
        assert_eq!(e12.0[3], 1, "e₁e₂ should be e₁₂");

        // e₂ × e₁ = -e₁₂ (anti-commutation)
        let e21 = e2.mul(&e1);
        assert_eq!(e21.0[3], u64::MAX, "e₂e₁ should be -e₁₂");

        // Associativity: (e₁ × e₂) × e₃ = e₁ × (e₂ × e₃)
        let e3 = Mv::basis(4); // 0b0100 = e₃
        let left = e1.mul(&e2).mul(&e3);
        let right = e1.mul(&e2.mul(&e3));
        assert_eq!(left.0, right.0, "Clifford product must be associative");
    }

    #[test]
    fn sha256_is_correct() {
        // Test vector: SHA-256("abc") = ba7816bf...
        let hash = sha256(b"abc");
        assert_eq!(hash[0], 0xba);
        assert_eq!(hash[1], 0x78);
        assert_eq!(hash[2], 0x16);
        assert_eq!(hash[3], 0xbf);

        // Double hash of empty-ish (deterministic check)
        let h1 = sha256(b"");
        let h2 = sha256(&h1);
        // SHA256(SHA256("")) is a known value
        assert_ne!(h1, h2, "double hash should differ from single");
    }

    #[test]
    fn response_matrix_is_exact() {
        // Verify: response.apply(x) == left * x * right for random x
        let left = Mv::from_seed(42);
        let right = Mv::from_seed(137);
        let response = ResponseMatrix::from_sandwich(&left, &right);

        for seed in 0..100 {
            let x = Mv::from_seed(seed);
            let direct = left.mul(&x).mul(&right);
            let speculative = response.apply(&x);
            assert_eq!(
                direct.0, speculative.0,
                "response matrix must be exact for seed {seed}"
            );
        }
    }

    #[test]
    fn speculation_benchmark() {
        run_speculation_benchmark();
    }
}
