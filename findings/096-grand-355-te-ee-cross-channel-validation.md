# Finding 096 — GRAND-355 TE/EE Cross-Channel Validation

Date: 2026-02-27
Status: COMPLETE

## Summary
Implemented cross-channel validation across TT/TE/EE using the same structurally derived `tau_reio` (no channel-specific fitting).

New file:
- `crates/gutoe-physics/src/bin/cmb_te_ee_crosscheck.rs`

Updated:
- `crates/gutoe-physics/src/cmb_class.rs` (generic CLASS column parser)
- `crates/gutoe-physics/data/README.md`

Data snapshots added:
- `crates/gutoe-physics/data/COM_PowerSpect_CMB-TE-full_R3.01.txt`
- `crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-full_R3.01.txt`
- `crates/gutoe-physics/data/COM_PowerSpect_CMB-TE-binned_R3.02.txt`
- `crates/gutoe-physics/data/COM_PowerSpect_CMB-EE-binned_R3.02.txt`

## Structural setup
- Assumed lane: `tau_reio = 0.054`
- Derived lane: `tau_reio = 0.067531` from structural reionization timing (`z_reion = 9.035`)

## Results (full-spectrum reduced chi2)
From:
- `GUTOE_CLASS_BIN=/tmp/class_public/class cargo run -q -p gutoe-physics --bin cmb_te_ee_crosscheck`

- TT full red: `1.607 -> 1.269`
- TE full red: `1.192 -> 1.139`
- EE full red: `1.095 -> 1.067`

All three channels improve under the same derived parameter change.

## Interpretation
This is the intended no-overfit behavior:
- one upstream structural refinement (`tau_reio`) propagates through TT/TE/EE,
- no per-channel retuning required,
- polarization channels improve in the same direction as TT.

This validates the “derive first, cross-check second” order for the CMB lane.
