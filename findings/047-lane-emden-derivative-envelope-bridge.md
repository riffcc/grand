# Finding 047 — Derivative-Sign Envelope Bridge for Lane-Emden Profiles

## Scope
Advance GRAND-284 by replacing a raw `avg(theta) <= 1` assumption with a theorem-level envelope route derived from monotonicity assumptions on the ODE branch.

## Lean (`lean/Gutoe/StellarFusion.lean`)
Added:
- `lane_emden_average_theta_le_one_of_sample_bound`
- `lane_emden_theta_le_one_of_deriv_nonpos_on_nonneg`
- `lane_emden_average_theta_le_one_of_deriv_nonpos_on_nonneg`
- `polytropic_ignition_from_lane_emden_profile_deriv_nonpos`

What this buys us:
- If a Lane-Emden branch is continuous on `ξ>=0`, differentiable on `ξ>0`, has nonpositive derivative on `ξ>0`, and satisfies `θ(0)=1`, then sampled profile averages are proven `<=1`.
- That bound now feeds directly into the ignition bridge theorem, reducing ad-hoc envelope assumptions.

## Rust parity (`crates/gutoe-physics/src/equations.rs`)
Added:
- `lane_emden_profile_all_nonnegative_xi`
- `lane_emden_profile_is_nonincreasing`
- `lane_emden_envelope_le_one_from_monotone_profile`

Tests added:
- `lane_emden_n1_profile_is_nonincreasing_on_nonnegative_window`
- `lane_emden_n3_profile_is_nonincreasing_on_nonnegative_window`

## Verification
- `lake build Gutoe.StellarFusion` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-physics lane_emden_n -- --nocapture` ✅

## Remaining gap
Theorems now consume derivative-sign assumptions cleanly, but we still need to *derive* nonpositivity (`θ' <= 0`) from full Lane-Emden ODE dynamics for broad `n` classes (currently tracked by GRAND-284).
