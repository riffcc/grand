# Finding 097 — GRAND-344 Differential Envelope Cross-Check

Date: 2026-02-27
Status: COMPLETE (negative result, retained)

## Summary
Implemented and tested a **differential** (not absolute) CMB envelope operator intended to reconcile microphysics-derived damping scale with CLASS baseline damping, then propagated it across TT/TE/EE with a single operator.

New files:
- `crates/gutoe-physics/src/cmb_differential.rs`
- `crates/gutoe-physics/src/bin/cmb_differential_crosscheck.rs`

Updated:
- `crates/gutoe-physics/src/lib.rs`

## Operator
Let
- `ell_diff_struct` from microphysics diffusion derivation,
- `ell_diff_class` from CLASS high-`ell` tail estimate,
- and gate at `ell_transition ~ ell_peak1`.

Apply:
- `F(ell) = exp( -ell^2 * (1/ell_struct^2 - 1/ell_class^2) * gate(ell) )`
- to all channels TT/TE/EE identically.

## Run result
From:
- `GUTOE_CLASS_BIN=/tmp/class_public/class cargo run -q -p gutoe-physics --bin cmb_differential_crosscheck`

Observed scales:
- `ell_diff_struct = 1598.6`
- `ell_diff_class = 1389.1`

Because `ell_diff_struct > ell_diff_class`, the differential factor amplifies high-`ell` instead of suppressing it.

Fit impact (worse in all channels):
- TT full red: `1.269 -> 15.720`
- TE full red: `1.139 -> 1.467`
- EE full red: `1.067 -> 1.206`
- EE band `1200..1600` mean `|pull|`: `0.721 -> 2.071`

## Interpretation
This differential pass is falsified in current form.

It still provides useful diagnostic information:
- Current microphysics diffusion lane predicts a weaker damping envelope than CLASS effective damping in this multipole band.
- Therefore this operator direction cannot be used as a refinement without revisiting either:
  - diffusion-scale derivation details, or
  - the mapping from microphysics scales into effective `ell`-space envelope.

Retained as a reproducible negative ablation, not promoted into the forward lane.
