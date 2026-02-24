# Finding 007: Wilson Loop Confinement Sweep

Date: 2026-02-24  
Related issue: GRAND-96 (Confinement / Wilson loops)

## What was run

- Binary: `cargo run --release -p gutoe-em --bin wilson_sweep`
- Config: `24x24x1`, `seeds=6`, `therm=300`, `meas=140`
- Output CSV: `/tmp/bh_renders/confinement_wilson_sweep.csv`

## Result summary

Across increasing `beta` (weakening confinement), we observe:

- Plaquette mean rises monotonically: `0.099 -> 0.762`
- Effective potential proxy `V3 = -ln <W_triangle>` falls monotonically: `2.391 -> 0.270`

Representative points:

- `beta=0.20`: plaquette `0.09896±0.00215`, `V3=2.39101±0.07226`
- `beta=1.00`: plaquette `0.43299±0.00219`, `V3=0.84605±0.05643`
- `beta=3.00`: plaquette `0.76208±0.00145`, `V3=0.27009±0.01777`

## Interpretation

- Low-`beta` (strong-coupling) regime exhibits larger `V3`, consistent with stronger confinement / area-law behavior.
- High-`beta` regime trends toward deconfinement signatures (larger plaquette, lower `V3`).
- This quantitatively extends the existing in-progress confinement work with a reproducible multi-seed sweep artifact.

