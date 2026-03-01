# Finding 223 — Double-SHA256 Z3 Mining Structure Probe (Angles 16/17/18/19)

## Summary

Implemented and ran a dedicated H2 probe where:

- `H2(m) = SHA256(SHA256(m))`

Binary:

- `crates/gutoe-physics/src/bin/sha256d_z3_mining_probe.rs`

Outputs:

- `/tmp/bh_renders/sha256d_z3_mining_probe/sha256d_z3_mining_probe.txt`
- `/tmp/bh_renders/sha256d_z3_mining_probe/sha256d_z3_mining_probe.json`

## Angle 16 — Length/Padding Boundary Scan

- Tested lengths: `0,1,2,3,4,31,32,33,47,48,54,55`
- Tested nontrivial order-3 transform pairs (byte/word cyclic variants).
- Exact equivariance hits: `0` across all tested pairs/cases.
- For operational lengths (`>=3`), mean Hamming distance stayed near random baseline (~128 bits).

## Angle 17 — Parallel Branch Search

- Candidate pairs: `36`
- Trials per pair: `5000`
- Threads used: `16`
- Best pair by mean-HD:
  - `pi=triplet_bytes_off1`
  - `rho=triplet_words_off0`
  - `mean_hd_bits=127.76`
  - `exact_hits=0/5000`
- Retro-observed wrapper (sim semantics):
  - `host_elapsed_s=3.44423667e-1`
  - `predeparture=true`
  - `apparent_speedup=3.44423667e8`

## Angle 18 — Orbit Deduper

Domain: 3-byte messages over alphabet `0..15` (`4096` points), mining predicate:

- `leading_zero_bits(sha256d(m)) >= 12`

Results:

- `baseline_hashes=4096`
- `dedup_hashes=1376`
- `ideal_speedup_if_equivariant=2.976744...`
- Actual inference mismatch vs truth: `2/4096`
- `mismatch_rate=4.8828125e-4`

Interpretation: dedup can reduce evaluations, but nonzero mismatch means orbit inference is not exact for H2 under tested transforms.

## Angle 19 — Direct H2 Bias Test

Using best pair from angle 17:

- Exact hits: `0/100000`
- `mean_hd_pair=127.98522`
- Random baseline: `mean_hd_random=127.96009`
- Leading-zero correlation:
  - `lz_corr=5.29539e-4`
  - success agreement `0.99949`
  - independent baseline `0.9994901276`

Interpretation: no detectable exploitable Z3 bias in double-SHA256 under tested actions at this sample scale.

## Scope note

This is a controlled structural probe (single-block SHA lane), not a full mining stack implementation.
