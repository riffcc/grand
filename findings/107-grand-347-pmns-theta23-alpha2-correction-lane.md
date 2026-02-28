# GRAND-347 — PMNS θ23 alpha² Correction Lane (Implemented)

## Scope

Implement optional second-order PMNS correction lane from structural ansatz:

- `sin²(θ23) = 4/7 - c α²`
- default coefficient `c = 137/4` (configurable)

Goal: reduce direct PMNS `θ23` residual by >=10x while preserving envelope and CP witness gates.

## Code changes

1. New corrected PMNS observable function

- `crates/gutoe-em/src/flavor.rs`
- Added:
  - `pmns_from_clifford_theta23_alpha2(c_alpha2: f64) -> MixingObservables`

2. Public export

- `crates/gutoe-em/src/lib.rs`
- Re-exported:
  - `pmns_from_clifford_theta23_alpha2`

3. Report lane wiring

- `crates/gutoe-em/src/bin/flavor_mix_report.rs`
- Added corrected PMNS block and JSON object:
  - `pmns.direct_theta23_alpha2_corrected`
- Config env:
  - `GUTOE_PMNS_TH23_ALPHA2_C` (default `137/4`)

4. CI gate wiring

- `crates/gutoe-em/src/bin/flavor_ci_gate.rs`
- Added corrected PMNS envelope+CP checks
- Added hard improvement gate:
  - `|Δθ23_corrected| <= |Δθ23_direct| / 10`
- JSON now includes:
  - `pmns_theta23_improvement` summary block

## Verification

Commands:

- `cargo run -q -p gutoe-em --bin flavor_mix_report`
- `cargo run -q -p gutoe-em --bin flavor_ci_gate`

Observed direct vs corrected:

- direct PMNS `θ23 = 49.106605°` (residual `+0.106605°`)
- corrected PMNS `θ23 = 49.001051°` with `c=137/4` (residual `+0.001051°`)

Improvement:

- factor ≈ `101.4x` (passes 10x requirement)

Gate status:

- CKM direct/texture: pass
- PMNS direct/texture: pass
- PMNS corrected: pass
- PMNS θ23 improvement gate: pass

Artifacts:

- `/tmp/bh_renders/flavor_mix_report.txt`
- `/tmp/bh_renders/flavor_mix_report.json`
- `/tmp/bh_renders/flavor_ci_gate.json`

## Notes

- This correction lane is additive/optional and does not replace legacy direct PMNS output.
- The coefficient remains configurable for future structural/Lean refinement.
