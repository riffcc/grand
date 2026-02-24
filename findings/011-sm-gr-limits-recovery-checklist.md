# 011 — Correct SM + GR Limits Recovered (Checklist)

Status: verification checklist and evidence index for GRAND-123.

## Standard Model side (structural limits)

1. Electroweak structure
- `sin^2(theta_W) = 3/13`: `lean/Gutoe/Z3Uniqueness.lean`, `lean/Gutoe/GaugeConstants.lean`.
- `mZ/mW = sqrt(13/10)`: `lean/Gutoe/GaugeConstants.lean`.

2. Gauge algebra size
- Total generators = 12: `lean/Gutoe/GaugeGroupSM.lean`, `lean/Gutoe/GaugeConstants.lean`.

3. Alpha consistency split
- Structural LO: `alpha_LO = 1/137`.
- Runtime measured: `alpha ~= 1/137.036`.
- Reference: `crates/gutoe-physics/src/constants.rs`, findings `002-alpha-consistency-rust-lean.md`.

## GR side (classical/Schwarzschild/Kerr limits)

1. Schwarzschild/classical recovery
- Classical limit statements in `lean/Gutoe/GravityMetric.lean`.
- Runtime metric checks in `crates/gutoe-gpu/src/metric.rs`.

2. Kerr consistency limits
- Horizon/static-limit behavior in `lean/Gutoe/KerrGeometry.lean`.
- Camera constants and projected limits in `lean/Gutoe/KerrCameraStability.lean`.
- Runtime reference: `crates/gutoe-gpu/src/kerr.rs` and tracer tests.

3. Transfer positivity/boundedness
- Proven constraints: `lean/Gutoe/SynchrotronTransfer.lean`.
- Runtime transfer tests: `crates/gutoe-gpu/src/transfer.rs`.

## Evidence quality gates

- Lean proofs compile (`lake build Gutoe`).
- Runtime parity tests pass for mapped coefficients and limits.
- Any mismatch gets filed as a parity issue before claiming limit recovery.

## Current gap summary

- Several high-level limits are structurally proven, but full observational closure still requires ongoing tickets (e.g., GRAND-217/218).
- This issue should be marked complete only as a checklist/evidence-index milestone, not as full phenomenology closure.
