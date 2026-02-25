# Finding 048 — Flux-Derivative Window Envelope for General Lane-Emden Branches

## Scope
Push GRAND-284 beyond pure derivative-sign assumptions by proving a finite-window envelope bound from the Lane-Emden divergence-form flux derivative.

## Lean (`lean/Gutoe/StellarFusion.lean`)
Added theorem:
- `lane_emden_average_theta_le_one_of_flux_deriv_negative_on_window`

This theorem proves sampled-average envelope (`avg theta <= 1`) on `[0,a]` from:
- positivity of `theta` on `(0,a)`
- continuity/differentiability assumptions
- divergence-form derivative identity for flux `flux(xi)=xi^2*theta'(xi)`
- regular center condition `theta'(0)=0`
- sample points constrained to `[0,a]`

Mechanism:
- derive `deriv flux < 0`
- obtain `StrictAntiOn flux`
- infer `theta' < 0` on `(0,a)`
- derive `AntitoneOn theta`
- conclude finite sampled average bound

## Rust parity (`crates/gutoe-physics/src/equations.rs`)
Added monotonic-envelope witness utilities over sampled profiles:
- `lane_emden_profile_all_nonnegative_xi`
- `lane_emden_profile_is_nonincreasing`
- `lane_emden_envelope_le_one_from_monotone_profile`

Added tests:
- `lane_emden_n1_profile_is_nonincreasing_on_nonnegative_window`
- `lane_emden_n3_profile_is_nonincreasing_on_nonnegative_window`

## Verification
- `lake build Gutoe.StellarFusion` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-physics lane_emden_n -- --nocapture` ✅

## Remaining gap
To fully close GRAND-284, we still need to derive the flux-derivative hypothesis itself directly from the full Lane-Emden ODE solution contract for a broad `n` class, not assume it as input.
