# GRAND-363 — PoC Push: Alpha Correction Lane, 3-Mode G Bridge, PMNS θ23 Triage

## Scope

Single-session push on three fronts:

1. Canonicalize the `0.036` alpha-correction lane in reproducible CI artifacts.
2. Implement third-mode `G` bridge using muon/phase-derived electron input.
3. Triage PMNS `θ23` tension as perturbative-fixable vs hard mismatch.

## 1) Alpha correction lane (canonicalized)

Updated unified report:

- `crates/gutoe-physics/src/bin/alpha_web_ci_report.rs`

Now emits explicit block:

- `delta_target = α_inv_physical - 137`
- `delta_first_order_5alpha`
- `delta_second_order_5alpha_minus_9alpha2`
- first/second absolute errors
- `second_order_improves` boolean (gated)

Run:

- `cargo run -q -p gutoe-physics --bin alpha_web_ci_report`

Artifact:

- `/tmp/bh_renders/alpha_web_ci_report/alpha_web_ci_report.txt`
- `/tmp/bh_renders/alpha_web_ci_report/alpha_web_ci_report.json`

Observed:

- `delta_target = 3.5999084e-2`
- first-order `5α = 3.6486763e-2` (abs err `4.8768e-4`)
- second-order `5α - 9α² = 3.6007501e-2` (abs err `8.4167e-6`)
- improvement factor ~58x
- `ci_gate passes_all=true`

## 2) G bridge third mode (implemented)

Updated:

- `crates/gutoe-physics/src/bin/g_bridge_report.rs`

Added mode:

- `mode_muon_phase_alpha2_electron`
- phase model: `δ = 3π/4 - 5α*(13/12) - (15/16)α²`

Run:

- `cargo run -q -p gutoe-physics --bin g_bridge_report`

Artifact:

- `/tmp/bh_renders/g_bridge_report/g_bridge_report.txt`
- `/tmp/bh_renders/g_bridge_report/g_bridge_report.json`

Observed `G` relative errors:

- measured-electron mode: `-9.699e-4`
- proton-anchor structural mode: `-1.136e-3`
- muon-phase mode: `-1.720e-3`

This closes GRAND-345 implementation scope for a third independent mode.

## 3) PMNS θ23 triage (result: perturbative-fixable)

Current direct structural value from existing flavor lane:

- `θ23_direct = 49.106605°` (target `49.0°`, `+0.1066°`)

Model triage with minimal correction ansatz:

- `sin²(θ23) = 4/7 - c α²`

Solved target coefficient:

- `c_target ≈ 34.591`

Short structural-rational candidates from existing primitives already near target:

- `c = 137/4 = 34.25` → `θ23 = 49.001051°` (`+0.001051°`)
- `c = 67/2 = 33.5` → `θ23 = 49.003362°` (`+0.003362°`)
- `c = 33` → `θ23 = 49.004902°` (`+0.004902°`)

Interpretation:

- This does **not** look like a hard-theory conflict.
- It behaves like a small second-order correction lane (same pattern as other sectors).

## Additional hard gate added

In `crates/gutoe-em/src/alpha.rs`:

- test `structural_alpha_identity_and_lane_regression_gate`

Asserts:

- exact identity `triangular(2^4)+1=137`
- structural-alpha lepton lane remains in bounded error band.

Run:

- `cargo test -p gutoe-em structural_alpha_identity_and_lane_regression_gate -- --nocapture`

Result: pass.

## Status summary

Tonight’s push converted three previously “open” PoC items into executable artifacts/gates:

- alpha residual lane is now explicit and gated,
- G has three-mode bridge reporting,
- PMNS θ23 classified as likely perturbative refinement, not a structural break.
