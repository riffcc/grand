# 029 — Strong-CP Structural Gate (GRAND-125)

Status: theorem-level structural closure, with explicit scope boundary.

## What is now proven

- In `lean/Gutoe/StrongCP.lean`:
  - `cp_odd_sector_imbalance_zero`
  - `theta_qcd_structural_zero`
  - `neutron_edm_from_structural_theta_zero`
  - `neutron_edm_structural_within_catalog_bound`
  - `theta_zero_of_proportional_cp_odd_source`

The core structural chain is:

1. Cl(1,3) Lorentz bivector split is exactly `3 + 3` (rotations + boosts).
2. CP-odd structural source is modeled as that grade-2 imbalance.
3. Imbalance is zero, so structural `theta_QCD` is zero.
4. The bridge EDM estimate `|d_n| ≈ 2.4e-16 * |theta_QCD| e·cm` is therefore zero.

## Runtime parity status

- `crates/gutoe-physics/src/bin/theorem_parity.rs` now includes:
  - `theta_qcd_structural = 0`
  - `neutron_edm_from_theta_structural = 0`
- Current run: parity rows all within tolerance.

## Scope boundary (explicit)

This is a **structural theorem closure** for the current GUTOE formalization.
It is not yet a full nonperturbative QCD path-integral derivation of the physical
vacuum angle from first principles.

That deeper derivation remains a separate milestone.
