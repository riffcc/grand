# 144 - Everyday Extremes + Lean Constraint Tightening

Date: 2026-02-28  
Scope: Implement four new everyday-physics lanes and convert their gate windows into Lean constraints.

## Rust Lane Additions

- Module:
  - `crates/gutoe-physics/src/everyday_extremes.rs`
- Report binary:
  - `crates/gutoe-physics/src/bin/everyday_extremes_report.rs`
- CI gate binary:
  - `crates/gutoe-physics/src/bin/everyday_extremes_ci_gate.rs`
- Library wiring:
  - `crates/gutoe-physics/src/lib.rs`

## Returned Physics Numbers

From `/tmp/bh_renders/everyday_extremes/everyday_extremes_report.txt`:

- Ice slipperiness (fixed pressure/speed sweep):
  - `μ(-20°C) = 0.2564`
  - `μ(-2°C) = 0.1151`
  - drop `Δμ = 0.1413`
- Popcorn:
  - pop temperature `155.0°C`
  - internal pressure `1.090705 MPa`
  - rupture threshold `1.066667 MPa`
  - expansion ratio `23.766`
- Raindrop shape:
  - optimal stable diameter `5.900 mm`
  - aspect ratio `0.6731`
  - Weber number `8.7778` (stable regime)
- Mpemba:
  - hot freeze time `199.047 min`
  - cold freeze time `274.865 min`
  - hot faster by `75.818 min`
  - sweep fraction hot-faster `16/16 = 1.0`

CI gate:
- `/tmp/bh_renders/everyday_extremes_ci_gate.json`
- `overall_pass = true`

## Lean Constraint Tightening

Added:
- `lean/Gutoe/EverydayExtremes.lean`
- `lean/lakefile.lean` root registration (`Gutoe.EverydayExtremes`)

Formalized gate constraints as theorem-level inequalities:

- `alpha_leading_order_q` (ties lane to `alphaInverse 4 = 137`)
- `ice_slipperiness_drop_positive` (`Δμ > 0.1`)
- `popcorn_temperature_gate` (`150 ≤ T_pop ≤ 220`)
- `popcorn_expansion_gate` (`expansion ≥ 10`)
- `raindrop_diameter_window` (`2.5 ≤ d_opt ≤ 6`)
- `raindrop_aspect_window` (`0.6 ≤ aspect ≤ 0.8`)
- `mpemba_default_ordering` (`t_hot < t_cold`)
- `mpemba_sweep_fraction_gate` (`fraction ≥ 0.2`)
- Aggregated as `everyday_extremes_constraint_spine`

This converts runtime behavior into explicit Lean gate invariants so regressions become theorem-visible, not just report-visible.

## Verification

- `cargo run -q -p gutoe-physics --bin everyday_extremes_report` passed.
- `cargo run -q -p gutoe-physics --bin everyday_extremes_ci_gate` passed (`overall_pass=true`).
- `cargo test -q -p gutoe-physics --lib everyday_extremes` passed (4/4).
- `cd lean && lake build Gutoe.EverydayExtremes` passed.
- `cd lean && lake build Gutoe` passed (warnings only, no errors).
