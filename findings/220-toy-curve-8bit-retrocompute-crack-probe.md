# Finding 220 — Toy Curve 8-bit Retrocompute Crack Probe

## Summary

Implemented and executed a toy ECC discrete-log crack benchmark for an 8-bit
key space, with a retrocompute observability overlay.

Binary:

- `crates/gutoe-physics/src/bin/ctc_toy_curve_key_crack_probe.rs`

## Run

```bash
cargo run -p gutoe-physics --bin ctc_toy_curve_key_crack_probe
```

Outputs:

- `/tmp/bh_renders/ctc_toy_curve_key_crack_probe/ctc_toy_curve_key_crack_probe.txt`
- `/tmp/bh_renders/ctc_toy_curve_key_crack_probe/ctc_toy_curve_key_crack_probe.json`

## Default-run metrics

- `trials = 20000`
- `solved = 20000`
- `success_rate = 1.0`
- `avg_guesses = 128.4696` (expected near midpoint of 8-bit keyspace)
- `normal_avg_time_s = 6.080978e-5`
- `retro_observed_avg_latency_s = 1.0e-9` (instrument floor)
- `apparent_speedup = 6.080978e4`

## 32-bit extension run

Using `GUTOE_TOY_CURVE_KEY_BITS=32` (estimated-scaling mode):

- `mode = estimated_scaling`
- `keyspace_size = 4294967296`
- `avg_guesses = 2147483648`
- `normal_avg_time_s = 3.429730e3` (~57.2 minutes average)
- `retro_observed_avg_latency_s = 3.429730e-3`
- `apparent_speedup = 1.0e6`

## Notes

- This is strictly a toy 8-bit lane for sanity checks.
- Internal compute work is unchanged; speedup is in observer-facing latency in
  the retrocompute model.
