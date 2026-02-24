# 013 — Unitarity, Causality, and Hawking Sign Audit

Issues targeted:
- GRAND-116 (`Unitarity of time evolution`)
- GRAND-117 (`Causality / light cone structure from lattice`)
- GRAND-91 (`Hawking temperature correction`) adjudication + fix

## GRAND-116: Unitarity evidence

Code anchors:
- `crates/gutoe-qm/src/lib.rs`
- `crates/gutoe-qm/src/gates.rs`

Executed:
- `cargo test -p gutoe-qm --lib -- --nocapture`

Result:
- `13 passed; 0 failed`.
- Explicit norm-preservation tests for hop/phase/Z3 gates pass.
- Interference tests also preserve coherence/visibility in unitary evolution.

Conclusion:
- Runtime quantum evolution gates satisfy unitarity checks in current test suite.

## GRAND-117: Causality/light-cone evidence

Code anchors:
- `crates/gutoe-physics/src/equations.rs` (`WaveEquation::group_velocity`)
- `lean/Gutoe/LorentzInvariance.lean`

Executed:
- `cargo test -p gutoe-physics group_velocity_bounded_by_phase_velocity -- --nocapture`
- `cargo test -p gutoe-physics group_velocity_is_positive -- --nocapture`

Result:
- both targeted causality tests pass.
- test assertions explicitly enforce `v_g <= v` and `v_g >= 0`.

Conclusion:
- Current lattice-wave runtime checks satisfy causality constraints used by the model.

## GRAND-91: Hawking sign adjudication (resolved)

Conflict:
- `lean/Gutoe/GravityMetric.lean` proves `hawking_temp_gt_gr` (hotter).
- `lean/Gutoe/HawkingCorrection.lean` proves `gutoe_hawking_cooler` (cooler).

Decision:
- Canonical branch for GRAND-91 is the subluminal-dispersion sign (`cooler`).
- `GravityMetric` was updated to match `HawkingCorrection`:
  - `hawking_temp` uses `1 - lambda_qg * (...)^2`
  - theorem renamed to `hawking_temp_lt_gr`
  - master theorem updated accordingly.
- Rust parity updated in `crates/gutoe-gpu/src/metric.rs`:
  - Hawking correction sign set negative
  - test renamed to `hawking_temperature_is_below_gr`

Validation:
- `lake build Gutoe.GravityMetric` passed.
- `lake build Gutoe` passed.
- `cargo test -p gutoe-gpu hawking_temperature_is_below_gr -- --nocapture` passed.
