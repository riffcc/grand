# Finding 222 — SHA256 Z3 Equivariance Probe

## Summary

Added and ran a direct SHA256 equivariance test for candidate order-3 transforms:

- `crates/gutoe-physics/src/bin/sha256_z3_equivariance_probe.rs`

Tested empirical condition:

- `H(pi(m)) == rho(H(m))`

with nontrivial order-3 `pi` and `rho` maps.

## Tested transforms

Message transforms (`pi`, order 3):

- `triplet_bytes` (rotate each 3-byte block)
- `triplet_words` (rotate each 3-word block over 48-byte message)

Digest transforms (`rho`, order 3):

- `triplet_bytes` over first 30 bytes (last 2 fixed)
- `triplet_words` over first 6 words (last 2 fixed)

All order-3 checks passed (`T^3 = identity`) on sample vectors.

## Run

```bash
cargo run -p gutoe-physics --bin sha256_z3_equivariance_probe
```

Outputs:

- `/tmp/bh_renders/sha256_z3_equivariance_probe/sha256_z3_equivariance_probe.txt`
- `/tmp/bh_renders/sha256_z3_equivariance_probe/sha256_z3_equivariance_probe.json`

## Results

### Default run (20,000 trials each pair)

- `pi=triplet_bytes, rho=triplet_bytes`: `0 / 20000`
- `pi=triplet_bytes, rho=triplet_words`: `0 / 20000`
- `pi=triplet_words, rho=triplet_bytes`: `0 / 20000`
- `pi=triplet_words, rho=triplet_words`: `0 / 20000`

No detectable nontrivial Z3-equivariance signal for these candidates.

### High-stat rerun (200,000 trials each pair)

- `pi=triplet_bytes, rho=triplet_bytes`: `0 / 200000`
- `pi=triplet_bytes, rho=triplet_words`: `0 / 200000`
- `pi=triplet_words, rho=triplet_bytes`: `0 / 200000`
- `pi=triplet_words, rho=triplet_words`: `0 / 200000`

## Scope note

This is an empirical candidate test, not a formal impossibility proof for all
possible group actions.
