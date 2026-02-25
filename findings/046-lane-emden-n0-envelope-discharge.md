# Finding 046 — n=0 Lane-Emden Envelope Bound Discharged

## Scope
Discharge the profile-envelope assumption (`avg(theta) <= 1`) for the exact n=0 Lane-Emden branch and wire it into the ignition bridge.

## Lean (`lean/Gutoe/StellarFusion.lean`)
Added:
- `lane_emden_theta_n0_le_one` (pointwise envelope)
- `lane_emden_average_theta_n0_le_one` (sampled average envelope)
- `polytropic_ignition_from_lane_emden_n0_profile`

Result:
- For n=0, profile-based ignition no longer requires an explicit external `avg(theta) <= 1` hypothesis.

## Rust parity (`crates/gutoe-physics/src/equations.rs`)
Added:
- `polytropic_ignition_condition_from_lane_emden_n0_profile`
- tests:
  - `lane_emden_n0_profile_average_theta_is_bounded_by_one`
  - `lane_emden_n0_profile_specialized_ignition_bridge_is_usable`

## Verification
- `lake build Gutoe.StellarFusion` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-physics lane_emden_n0_profile -- --nocapture` ✅

## Remaining
General-index envelope proof (beyond n=0) is still open and tracked in GRAND-284.
