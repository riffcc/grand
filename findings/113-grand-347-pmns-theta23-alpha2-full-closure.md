# GRAND-347 — PMNS `θ23` α² Lane Full Closure (Lean + Rust Parity)

## Reality

The optional PMNS correction lane is now structurally aligned with Lean and gated in runtime CI.

Baseline direct lane:

- `sin²(θ23) = 4/7`
- `θ23 = 49.106605°` (residual to 49.0° target: `0.106605°`)

Corrected lane (default structural coefficient):

- `sin²(θ23) = 4/7 - c α²`
- `c = 137/4` (structural default)
- `α = 1/137` (structural lane parity)
- `θ23 = 49.000995°` (residual: `0.000995°`)

Improvement:

- residual reduction factor: `~107x`

## Rust lane changes

Files:

- `crates/gutoe-em/src/flavor.rs`
- `crates/gutoe-em/src/lib.rs`
- `crates/gutoe-em/src/bin/flavor_mix_report.rs`
- `crates/gutoe-em/src/bin/flavor_ci_gate.rs`

What changed:

1. PMNS correction now uses structural alpha parity by default

- corrected lane computes `α` from `ALPHA_INVERSE_STRUCTURAL` (`1/137`)
- shared default coefficient exported as:
  - `PMNS_THETA23_ALPHA2_COEFF_STRUCTURAL = 137/4`

2. Explicit theta23 correction helpers added

- `pmns_theta23_sq_direct()`
- `pmns_theta23_sq_alpha2_corrected(c_alpha2)`

3. Runtime report now emits direct vs corrected improvement block

- `pmns.theta23_alpha2_improvement` in JSON with:
  - direct residual
  - corrected residual
  - improvement factor
  - `improves_10x` boolean

4. CI gate keeps correction optional but enforces hard 10x improvement when enabled

- `pmns_theta23_improvement.pass` must be true

## Lean lane changes

File:

- `lean/Gutoe/FlavorMixing.lean`

What changed:

- added shared coefficient definition from primitive:
  - `pmnsTheta23Alpha2CoeffQ := alphaInverse 4 / 4`
- theorem:
  - `pmns_theta23_alpha2_coeff_eq` (`= 137/4`)
- corrected lane now written via shared coefficient term
- added formal 10x improvement theorem:
  - `pmns_theta23_corrected_improves_10x`

Existing closure theorems retained:

- `pmns_theta23_void_term`
- `pmns_theta23_corrected_closed_form`
- `pmns_theta23_corrected_closer_than_direct`

## Verification

Commands run:

- `cargo run -q -p gutoe-em --bin flavor_mix_report`
- `cargo run -q -p gutoe-em --bin flavor_ci_gate`
- `cargo test -q -p gutoe-em pmns_theta23_alpha2 -- --nocapture`
- `cd lean && lake build Gutoe`

Key outputs:

- flavor report JSON:
  - `direct_abs_residual_deg = 0.106605350869`
  - `corrected_abs_residual_deg = 0.000995373260`
  - `improves_10x = true`
- flavor CI gate:
  - `overall_pass = true`
  - `pmns_theta23_improvement.pass = true`
- targeted Rust tests:
  - `2 passed, 0 failed`
- Lean build:
  - `Build completed successfully (8109 jobs)`

## Acceptance check

- one-command repro includes direct vs corrected PMNS in one artifact: **done**
- `θ23` residual improved by >=10x without breaking other observables: **done**
- Lean parity landed with explicit coefficient + improvement theorem: **done**
