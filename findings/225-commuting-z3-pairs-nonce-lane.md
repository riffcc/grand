# Finding 225 — Commuting Z3 Pair Test on Bitcoin-Style Nonce Lane

## Summary

Implemented and ran a direct test of the recursive `3^k` hinge in the strict
consensus-preserving nonce-action family.

Binary:

- `crates/gutoe-physics/src/bin/sha256d_btc_commuting_z3_pair_probe.rs`

Outputs:

- `/tmp/bh_renders/sha256d_btc_commuting_z3_pair_probe/sha256d_btc_commuting_z3_pair_probe.txt`
- `/tmp/bh_renders/sha256d_btc_commuting_z3_pair_probe/sha256d_btc_commuting_z3_pair_probe.json`

## Core result

- Nontrivial order-3 nonce permutations found: `8`
- Commuting pairs found: `4`
- **Independent commuting pairs found: `0`**
- All commuting pairs are in the same cyclic subgroup (`generated_group_size = 3`)

So in this nonce-only permutation family there is no `k=2` (`9x`) route.

## Joint equivariance metrics (best row)

- `eq_hits = 0`
- `eq_hit_rate = 0.0`
- `generated_group_size = 3`
- `same_subgroup = true`

No exact joint equivariance signal appears in this lane.

## Interpretation

The recursive elimination step requires two independent commuting Z3 actions.
That condition fails in the tested consensus-preserving nonce permutation family.
