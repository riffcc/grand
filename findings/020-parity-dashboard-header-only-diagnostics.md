# Finding 020: Parity Dashboard Header-Only Diagnostics

Date: 2026-02-24
Issue: GRAND-179

## Problem

`transfer_parity_*.csv` could exist with only a header row (no data), which made
parity dashboards look "generated" while hiding that no valid parity samples were captured.

## Change

Updated `crates/gutoe-gpu/src/bin/bh_parity_dashboard.rs` to:

- track header-only parity CSV files,
- include a dedicated diagnostics section in `parity_dashboard.md`,
- exit with status `2` when no real parity rows are present.

This turns silent false-green behavior into explicit diagnostic failure.

## Validation

- `cargo run -p gutoe-gpu --bin bh_parity_dashboard`
  - produced `/tmp/bh_renders/parity_dashboard.md`
  - exited with code `2` when only header rows existed
  - listed `transfer_parity_m87star.csv` under diagnostics

## Result

Parity collection remains blocked on non-empty CUDA/ROCm captures, but the dashboard now reports that state honestly and machine-detectably.
