# GRAND-290 + GRAND-291 (pass 1): Flavor trend persistence + neutrino splitting-ratio tightening

Date: 2026-02-28

## Scope landed

1. GRAND-290 trend persistence artifacts
- Added run-to-run trend files for flavor lanes:
  - `/tmp/bh_renders/flavor_mix_trend.csv`
  - `/tmp/bh_renders/flavor_ci_gate_trend.csv`
- Added versioned snapshots per run:
  - `/tmp/bh_renders/flavor_mix_report.<timestamp>.txt`
  - `/tmp/bh_renders/flavor_mix_report.<timestamp>.json`
  - `/tmp/bh_renders/flavor_ci_gate.<timestamp>.json`
- Header/schema guard added: legacy trend files are auto-archived when schema changes.

2. GRAND-291 root-lane tightening (neutrino split ratio)
- Added shared structural absolute-neutrino transduction in `gutoe-em::flavor`:
  - `neutrino_hierarchy_exponent_structural()`
  - `neutrino_absolute_masses_from_texture()`
- Structural exponent introduced from existing constants:
  - `p = α^{-1} / (|grade1| + |grade2|) = 137/10 = 13.7`
- Updated consumers to use shared lane:
  - `neutrino_tiny_mass_report`
  - `neutrino_ci_gate`
  - `neutrino_oscillation_ci_gate`
  - `ew_flavor_coupled_ci_gate`
  - `blind_prediction_register`
  - `blind_prediction_register_ci_gate`

## Verification run

Commands:

```bash
cargo test -q -p gutoe-em flavor
cargo run -q -p gutoe-em --bin flavor_mix_report
cargo run -q -p gutoe-em --bin flavor_ci_gate
cargo run -q -p gutoe-physics --bin ew_flavor_coupled_ci_gate
cargo run -q -p gutoe-em --bin neutrino_oscillation_ci_gate || true
cargo run -q -p gutoe-physics --bin blind_prediction_register
cargo run -q -p gutoe-physics --bin blind_prediction_register_ci_gate
```

Results:
- `gutoe-em flavor` tests: **19 passed, 0 failed**.
- Flavor trend files + snapshots emitted successfully.
- Coupled EW+flavor gate now passes:
  - `overall_pass = true`
  - `sin2_mz_bridge = 0.231195465518`
  - `splitting_ratio = 32.678997019738`
  - `splitting_ratio_rel_err = 3.150621926727e-3` (0.315%)
- Texture delta drift (GRAND-291 original symptom) remains:
  - CKM texture delta drift: `-3.291117673703 deg`
  - PMNS texture delta drift: `-5.726999545924 deg`

## Honest status

- GRAND-290 objective (artifact persistence + trend append) is implemented.
- GRAND-291 has **partial closure**:
  - The coupled EW+flavor splitting-ratio bottleneck is now tightened and gate-passing.
  - Texture-phase delta drift itself is still unresolved and needs another pass on phase placement/diagonalization conventions.
- Absolute oscillation-scale magnitudes remain open (`neutrino_oscillation_ci_gate` still fails on `dm21`/`dm32` amplitude, despite corrected ratio).

