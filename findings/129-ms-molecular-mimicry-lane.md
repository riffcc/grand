# MS Molecular Mimicry Lane — First Mechanistic Slice

Date: 2026-02-28  
Status: Implemented (report + CI gate)

## Goal
Build a reduced-order mechanistic lane for multiple sclerosis framed as a molecular
binding/misrecognition problem:
- TCR-like interface vs self epitope (MBP-like),
- TCR-like interface vs mimic epitope,
- tolerance-threshold excess and mimicry overlap,
- therapy-effect proxies and targeted-blocker feasibility.

## What shipped

- Module:
  - `crates/gutoe-physics/src/ms_autoimmunity.rs`
- Binaries:
  - `crates/gutoe-physics/src/bin/ms_autoimmunity_report.rs`
  - `crates/gutoe-physics/src/bin/ms_autoimmunity_ci_gate.rs`
- Wiring:
  - `crates/gutoe-physics/src/lib.rs` (module export)

## Executed output (default lane)

Command:
- `cargo run -q -p gutoe-physics --bin ms_autoimmunity_report`

Result:
- `self_binding_kj_mol = -29.909`
- `mimic_binding_kj_mol = -30.732`
- `mimicry_gap_kj_mol = 0.823`
- `tolerance_threshold_kj_mol = -29.000`
- `activation_excess_kj_mol = 0.909`
- `misrecognition_risk_index = 0.221`

Therapy proxy outputs:
- Combined drive reduction (ocrelizumab-like + natalizumab-like): `0.809`
- Residual drive index: `0.042`

Targeted blocker proxy:
- Required occupancy for tolerance restoration with safety buffer: `0.564`
- Feasible at default concentration: `true`

Artifacts:
- `/tmp/bh_renders/ms_autoimmunity/ms_autoimmunity_report.txt`
- `/tmp/bh_renders/ms_autoimmunity/ms_autoimmunity_report.json`

## CI gate

Command:
- `cargo run -q -p gutoe-physics --bin ms_autoimmunity_ci_gate`

Default windows:
- `mimicry_gap_kj_mol <= 2.0`
- `0.2 <= activation_excess_kj_mol <= 3.0`
- `therapy_relative_drive_reduction_fraction >= 0.5`
- `targeted blocker required occupancy <= 1.0` and feasible

Result:
- `overall_pass = true`
- Artifact: `/tmp/bh_renders/ms_autoimmunity_ci_gate.json`

## Honesty statement

This lane is a reduced mechanistic transduction, not full immunology/QM closure and
not a clinical dosing or treatment recommendation tool. It is intended to quantify
interface-margin logic and intervention leverage in a transparent, auditable form.
