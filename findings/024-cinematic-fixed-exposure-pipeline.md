# Finding 024: Cinematic Path + Fixed Exposure Pipeline

Date: 2026-02-24
Issue: GRAND-175

## What shipped

Implemented a fixed-exposure path for cinematic rendering and wired it into
`interstellar_spin`.

### 1) Fixed exposure override support

- CPU path:
  - Added `band_weight_with_exposure(...)` in `crates/gutoe-gpu/src/synchrotron.rs`.
  - `bh_render` now uses `BH_FIXED_EXPOSURE` (if set) via `fixed_exposure_override()`.
- GPU path:
  - Extended CUDA/HIP tracer interfaces (`tracer.cu`) so disk color uses
    `fixed_exposure` when `>= 0`.
  - FFI now passes fixed exposure from `bh_render` into `gutoe_render_bh`.

### 2) Cinematic sequence defaults

`render_interstellar_spin(...)` now forces a stable cinematic pipeline for the
sequence:

- `BH_SPECTRUM=optical`
- `BH_FIXED_EXPOSURE=1.4`

and restores previous env values afterwards.

## Validation

- `cargo test -p gutoe-gpu fixed_exposure_override_changes_band_weight -- --nocapture` ✅
- `cargo check -p gutoe-gpu --bin bh_render` ✅
- `cargo run -p gutoe-gpu --bin bh_render -- interstellar_spin 1 320x180` ✅
  - produced `/tmp/bh_renders/sgr_astar__spin_0000.png`
  - run log confirms `spectrum=optical`

## Result

Cinematic renders now run with physically stable camera choreography and a fixed
exposure path, reducing frame-to-frame exposure drift and making demos more
comparable.
