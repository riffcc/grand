# Finding 045 — Lane-Emden Profile-to-Ignition Bridge

## Scope
Advance GRAND-281 from standalone ODE scaffolding to a bridge that uses sampled Lane-Emden profiles directly in ignition criteria.

## Lean (`lean/Gutoe/StellarFusion.lean`)
Added:
- `laneEmdenAverageTheta`
- `laneEmdenProfileCompression`
- `proxy_compression_from_profile_threshold`
- `polytropic_ignition_from_lane_emden_profile`

Key result:
- If profile-weighted compression clears threshold and sampled profile obeys envelope bound `avg(theta) <= 1`, then base proxy compression also clears threshold, yielding ignition via existing theorem.

## Rust (`crates/gutoe-physics/src/equations.rs`)
Added:
- `lane_emden_average_theta_from_profile`
- `lane_emden_profile_weighted_compression`
- `polytropic_ignition_condition_from_lane_emden_profile`

This now consumes RK4 Lane-Emden trajectories in ignition checks.

## Tests
Added and passing:
- `lane_emden_profile_weighted_ignition_condition_tracks_threshold`
- `lane_emden_profile_ignition_rejects_avg_theta_above_one`

## Verification
- `lake build Gutoe.StellarFusion` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-physics lane_emden_profile -- --nocapture` ✅

## Remaining Gap
The envelope assumption `avg(theta) <= 1` is explicit and honest. Removing it requires an ODE-level monotonic/envelope proof for regular Lane-Emden branches (tracked as GRAND-284).
