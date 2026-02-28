# 163 — GRAND-291: Texture Phase-Fidelity Closure (CKM/PMNS delta drift)

Date: 2026-02-28

## Scope
Reduced texture-branch CP-phase drift (relative to direct Cl(1,3) branch) without introducing continuous fit knobs.

Files changed:
- `crates/gutoe-em/src/flavor.rs`

## Method
Reworked texture phase placement as a finite structural convention search:
- CKM texture: finite sign lattice over 5 phase placements (`2^5` candidates).
- PMNS texture: finite sign lattice over 6 phase placements (`2^6` candidates).
- Candidate selection objective:
  - minimize phase drift to the direct branch,
  - tie-break with angle/Jarlskog alignment,
  - require PDG envelope pass.

No continuous parameter fitting, no empirical delta injection.

## Results
From `cargo run -q -p gutoe-em --bin flavor_mix_report`:

Direct branch:
- CKM delta: `68.130°`
- PMNS delta: `198.435°`

Texture branch (after closure):
- CKM delta: `67.544°`
- PMNS delta: `198.753°`

Texture drift vs direct:
- CKM: `-0.586°` (was `-3.291°`)
- PMNS: `+0.318°` (was `-5.727°`)

## Regression guards
Added test:
- `texture_phase_conventions_keep_delta_drift_subdegree`

Guard thresholds:
- `|Δδ_ckm(texture-direct)| <= 1.0°`
- `|Δδ_pmns(texture-direct)| <= 1.0°`

## Verification
Commands run:

```bash
cargo test -q -p gutoe-em flavor
cargo run -q -p gutoe-em --bin flavor_mix_report
cargo run -q -p gutoe-em --bin flavor_ci_gate
```

Observed:
- flavor tests: `23 passed, 0 failed`
- flavor CI gate: `overall_pass=true`
