# Finding 022: Parity Dashboard Column-Safe CSV Parsing

Date: 2026-02-24
Issue: GRAND-179

## Problem

`bh_parity_dashboard` used fixed column indices and could misinterpret CSV schema
changes. In particular, when no `gpu_ms` column existed, it could accidentally
parse another column as timing data.

## Change

Updated `crates/gutoe-gpu/src/bin/bh_parity_dashboard.rs` to parse by header name,
not hardcoded index:

- `mad` by header key
- `transfer_delta_parity_abs` by header key
- `gpu_ms` only when the header explicitly contains `gpu_ms`

## Validation

- `cargo run -p gutoe-gpu --bin bh_parity_dashboard`
  - compiles/runs cleanly
  - still exits `2` with current header-only inputs (expected)

## Result

Dashboard aggregation is now schema-safe and won’t silently contaminate perf
metrics when CSV layouts evolve.
