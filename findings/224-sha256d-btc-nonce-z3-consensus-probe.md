# Finding 224 — SHA256d Bitcoin-Style Nonce Z3 Probe (Consensus-Preserving)

## Summary

Added and ran a strict lane that keeps PoW input shape intact:

- 80-byte header (`76-byte prefix + 4-byte nonce`)
- hash = `SHA256(SHA256(header))`
- only nonce bytes transformed by tested order-3 actions

Binary:

- `crates/gutoe-physics/src/bin/sha256d_btc_nonce_z3_probe.rs`

Outputs:

- `/tmp/bh_renders/sha256d_btc_nonce_z3_probe/sha256d_btc_nonce_z3_probe.txt`
- `/tmp/bh_renders/sha256d_btc_nonce_z3_probe/sha256d_btc_nonce_z3_probe.json`

## Results

Pair search (`200,000` trials per pair):

- Best pair: `pi=rotate3_high_bytes`, `rho=triplet_words_off0`
- `exact_hits = 0 / 200000`
- `mean_hd_bits = 127.99703` (near random 128)

Predicate correlation (`200,000` trials, threshold `leading_zero_bits >= 12`):

- `exact_hits = 0 / 200000`
- `mean_hd_pair = 127.95695`
- `mean_hd_random = 128.005825`
- `lz_corr = -5.1819579e-3`
- `success_agreement_rate = 0.999535`
- `success_agreement_independent_baseline = 0.9995351081`

Orbit deduper (nonce domain `65536`, same threshold):

- `baseline_hashes = 65536`
- `dedup_hashes = 65536`
- `ideal_speedup_if_equivariant = 1.0`
- `mismatch_rate = 0.0`

## Interpretation

No detectable exploitable Z3 structure was observed in this consensus-preserving
nonce lane at current sample sizes. Exact equivariance is absent; approximate
metrics remain at random-baseline levels.
