# Finding 021: Transfer Parity Backend-Tagged CSV Outputs

Date: 2026-02-24
Issue: GRAND-179

## Problem

Transfer parity reports were written to a single filename pattern:

- `transfer_parity_<view>.csv`

This can cause backend runs to overwrite each other and weakens multi-backend
(CUDA vs ROCm) comparison workflows.

## Change

Updated `run_transfer_parity_report` in `crates/gutoe-gpu/src/bin/bh_render.rs`
to include backend tags in output filenames.

New pattern:

- `transfer_parity_<view>_<backend>.csv`

Backend tag resolution:

- `BH_BACKEND_TAG` env var if provided
- else inferred from compile features (`cuda`, `rocm`, `multi`, fallback `gpu`)

## Validation

- `cargo check -p gutoe-gpu --bin bh_render` ✅

## Impact

- Prevents accidental overwrite across backend runs.
- Aligns with `bh_parity_dashboard` backend grouping semantics.
- Improves reproducibility for GRAND-179 parity/perf closure.
