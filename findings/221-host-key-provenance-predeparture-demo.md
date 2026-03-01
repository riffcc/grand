# Finding 221 — Host Key Provenance + In-Sim Predeparture Crack (Faithful Bridge)

## Summary

Implemented and validated a faithful host-to-sim lane where a real host key is
sampled from OS entropy, injected into simulated reality as public challenge
material, cracked in forward compute, and observed predeparture in simulated
coordinate time.

Binary:

- `crates/gutoe-physics/src/bin/ctc_hostkey_sim_predeparture_demo.rs`

## Important fix

The earlier toy-ECC variant could recover an equivalent scalar in a small group,
which broke strict `key_match`/commit verification despite a valid DLP solve.

This lane now uses a unique hash-preimage challenge (`blake3_preimage`) so
recovery is exact.

## Run

```bash
cargo run -p gutoe-physics --bin ctc_hostkey_sim_predeparture_demo
```

Outputs:

- `/tmp/bh_renders/ctc_hostkey_sim_predeparture_demo/ctc_hostkey_sim_predeparture_demo.txt`
- `/tmp/bh_renders/ctc_hostkey_sim_predeparture_demo/ctc_hostkey_sim_predeparture_demo.json`

## Default-run highlights

- `challenge_type = blake3_preimage`
- `key_match = true`
- `challenge_verified = true`
- `commitment_verified = true`
- `guesses_used = 382620`
- `host_crack_time_s = 4.190687e-1`
- `predeparture = true`
- `predeparture_margin_s = 8.381373e-2`
- `apparent_speedup = 4.190687e8`

## Interpretation

- Host provenance is now strict and auditable (`key_match` + commitment true).
- Compute still occurs forward on host hardware.
- Simulated observer sees predeparture availability via the retro-shift model.
