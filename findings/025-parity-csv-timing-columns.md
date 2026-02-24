# Finding 025: Parity CSV Timing Columns + Dashboard Speedup

Date: 2026-02-24
Issue: GRAND-179

## Change

Extended parity-report data model for performance diagnostics:

- `bh_render transfer_parity` now writes extra columns:
  - `backend`
  - `gpu_ms`
  - `cpu_ms`
- `bh_parity_dashboard` now aggregates:
  - mean `gpu_ms`
  - mean `cpu_ms`
  - speedup `cpu/gpu`

## Validation

- `cargo check -p gutoe-gpu --bin bh_render --bin bh_parity_dashboard` ✅
- `cargo run -p gutoe-gpu --bin bh_parity_dashboard` ✅ (still exits `2` on header-only data, by design)

## Notes

CUDA/ROCm kernel-side compilation for the updated tracer interface requires
backend-enabled builds/runs on GPU hosts; this sandbox run validated host-side
Rust paths only.
