# Finding 015: Rust/Lean Constant Parity Hardening

Date: 2026-02-24

## Summary

We tightened the runtime parity guardrails around the shared structural constants used by Lean proofs and Rust simulation code.

## Changes

- `crates/gutoe-core/src/states.rs`
  - Set `constants::LAMBDA_QG` to the structural value `1/12`.
  - Updated the associated unit test to enforce `|λ_QG - 1/12| < 1e-15`.
- `crates/gutoe-lattice/src/vector_rails.rs`
  - Replaced hardcoded test values (`0.084372`) with `gutoe_core::constants::LAMBDA_QG`.
- `crates/gutoe-physics/src/equations.rs`
  - Updated quadratic λ-mass-scaling test to use shared `LAMBDA_QG` as the Rust baseline.
- `crates/gutoe-physics/src/bin/theorem_parity.rs`
  - Added explicit `lambda_qg_core` parity row so the parity report now checks both:
    - `gutoe_physics::constants::LAMBDA_QG`
    - `gutoe_core::constants::LAMBDA_QG`

## Validation

- `lake build Gutoe` ✅
- `cargo test -p gutoe-core lambda_qg_matches_claimed_value -- --nocapture` ✅
- `cargo test -p gutoe-lattice test_dispersion_relation -- --nocapture` ✅
- `cargo test -p gutoe-lattice test_critical_wave_number -- --nocapture` ✅
- `cargo test -p gutoe-lattice test_stability_check -- --nocapture` ✅
- `cargo test -p gutoe-physics quark_mass_scales_as_lambda_qg_squared -- --nocapture` ✅
- `cargo run -p gutoe-physics --bin theorem_parity` ✅

Output artifact:
- `/tmp/bh_renders/theorem_runtime_parity.csv`

## Result

The parity report shows all structural constants in tolerance, with both core and physics stacks pinned to `λ_QG = 1/12`.
