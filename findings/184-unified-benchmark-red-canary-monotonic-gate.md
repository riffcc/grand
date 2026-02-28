# 184 — Unified Benchmark Red Canary (Monotonic Non-Increase)

## Goal
Prevent silent regressions in `elements_with_any_red` while allowing aggressive experimentation.

## What was added
File:
- `crates/gutoe-physics/src/bin/element_unified_external_benchmark.rs`

New behavior:
- Persistent best-state file for red-count baseline.
- Canary comparison on every run:
  - `current_red <= previous_best` => pass
  - `current_red > previous_best` => regression
- Strict mode default is **on**.
- Benchmark now emits canary details in text summary and dedicated canary JSON.

## New outputs
- `element_unified_external_benchmark.best_red_count`
- `element_unified_external_benchmark_canary.json`
- Extra lines in `element_unified_external_benchmark.txt`:
  - enabled/strict/state path
  - previous best/current/best after
  - improved/regressed/pass flags

## Environment controls
- `GUTOE_BENCH_RED_CANARY_ENABLED` (default `true`)
- `GUTOE_BENCH_RED_CANARY_STRICT` (default `true`)
- `GUTOE_BENCH_RED_CANARY_RESET` (default `false`)
- `GUTOE_BENCH_RED_CANARY_STATE` (optional custom baseline path)

## Validation
1. Baseline init pass:
- No prior state => baseline initialized to current red count.

2. Stable rerun pass:
- Prior best = 90, current = 90 => pass.

3. Forced regression fail:
- Injected baseline `80`, current `90` => process exits non-zero with:
  - `red canary regression: elements_with_any_red increased from 80 to 90`

## Current benchmark baseline (after rollback)
- `phase_accuracy = 1.000000`
- `ionization_mae_ev = 0.307058`
- `density_mae_g_cm3 = 2.829566`
- `melting_mae_k = 448.703068`
- `boiling_mae_k = 973.868170`
- `elements_with_any_red = 90`

This establishes a hard floor: future passes cannot land with `elements_with_any_red > 90` unless canary strict is intentionally disabled.
