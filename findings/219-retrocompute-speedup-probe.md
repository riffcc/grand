# Finding 219 — Retrocompute Speedup Probe (Finite Hardware, Large Apparent Speed)

## Summary

Added a dedicated retrocompute simulation lane that quantifies observer-facing
latency collapse when results are returned predeparture.

Binary:

- `crates/gutoe-physics/src/bin/ctc_retrocompute_speedup_probe.rs`

## Model

- Internal compute time: `task_flops / hardware_flops_per_s`
- Predeparture shift: `r * compute_time` with `r in [0,1]`
- External latency: `max(compute_time * (1-r), instrumentation_floor)`
- Apparent speedup: `compute_time / external_latency`
- Retry closure:
  `p_eventual = 1 - (1 - p_single)^retry_depth`

## Run

```bash
cargo run -p gutoe-physics --bin ctc_retrocompute_speedup_probe
```

Outputs:

- `/tmp/bh_renders/ctc_retrocompute_speedup_probe/ctc_retrocompute_speedup_probe.txt`
- `/tmp/bh_renders/ctc_retrocompute_speedup_probe/ctc_retrocompute_speedup_probe.json`

## Default-run highlights

- `compute_time_s_internal = 1.000000e5`
- `observed_latency_s_external = 1.000000e-1`
- `apparent_speedup = 1.000000e6`
- `p_single_pass = 0.12`
- `retry_depth = 100`
- `eventual_success_prob = 0.999997192840`
- `effective_speedup_with_retries = 8.333333e6`

## Notes

- This lane does not claim a physical CTC engine; it quantifies the operational
  consequence if the predeparture channel exists.
- The included near-unity `r` sweep shows unbounded trend in apparent speedup as
  external latency approaches the instrumentation floor.

