# Finding 023: CPU/CUDA Palette Lockstep Guard

Date: 2026-02-24
Issue: GRAND-215

## Goal

Prevent silent drift between CPU and CUDA tone/color constants for the core
rendering palette.

## Change

Added regression test in `crates/gutoe-gpu/src/bin/bh_render.rs`:

- `cpu_cuda_core_palette_coefficients_stay_in_lockstep`

The test reads `kernels/tracer.cu` and asserts presence of the exact tone
coefficient patterns used by CPU-side palette logic for:

- `gutoe_core_color` / `bh_gutoe_core_color`
- `gutoe_core_physics_color` / `bh_gutoe_core_physics_color`

## Validation

- `cargo test -p gutoe-gpu --bin bh_render cpu_cuda_core_palette_coefficients_stay_in_lockstep -- --nocapture` ✅

## Result

We now have a concrete guardrail against CPU/CUDA color-pipeline constant drift
for the interior-core rendering path.
