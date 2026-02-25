# Finding 043 — Lane-Emden n=0 Baseline (Lean/Rust parity)

## Scope
Start GRAND-281 by replacing pure proxy-only ignition structure with an exact Lane-Emden baseline that is analytically checkable.

## Lean
- Added exact n=0 profile and derivatives in `lean/Gutoe/StellarFusion.lean`:
  - `laneEmdenThetaN0`
  - `laneEmdenThetaN0Prime`
  - `laneEmdenThetaN0PrimePrime`
- Added multiplied-form residual:
  - `laneEmdenResidualN0 ξ = ξ² θ'' + 2ξ θ' + ξ²`
- Proved:
  - `lane_emden_residual_n0_zero` (exact algebraic zero for all `ξ`).

## Rust parity
- Added parity functions in `crates/gutoe-physics/src/equations.rs`:
  - `lane_emden_theta_n0`
  - `lane_emden_theta_n0_prime`
  - `lane_emden_theta_n0_prime_prime`
  - `lane_emden_residual_n0`
- Added tests:
  - `lane_emden_n0_profile_matches_closed_form_values`
  - `lane_emden_n0_residual_is_numerically_zero`

## Verification
- `lake build Gutoe.StellarFusion` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-physics lane_emden_n0 -- --nocapture` ✅

## Note
This is a baseline kill step: exact finite-checkable Lane-Emden structure is now in-tree and parity-verified. Remaining work is general-index ODE bridge and boundary-regular handling near the origin.
