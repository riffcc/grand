# GRAND (Cardiovascular lane) — Atorvastatin/HMGCR Binding Transduction

Date: 2026-02-28  
Status: Implemented (report + CI gate + residual decomposition)

## Goal
Open a cardiovascular-first molecular lane by anchoring one concrete, high-value benchmark:

- Drug: atorvastatin
- Target: HMG-CoA reductase (HMGCR)
- Observable: binding free energy from `Ki`

The lane computes:
1. exact thermodynamic benchmark `Ki -> ΔG`,
2. a first-principles QED electrostatic floor from `α ħ c / (ε r)`,
3. explicit residual stabilization still required beyond the floor,
4. a term-by-term non-electrostatic residual breakdown.

## What shipped

- New module:
  - `crates/gutoe-physics/src/cardiovascular_binding.rs`
- New binaries:
  - `crates/gutoe-physics/src/bin/cardiovascular_binding_report.rs`
  - `crates/gutoe-physics/src/bin/cardiovascular_binding_ci_gate.rs`
- Wiring:
  - `crates/gutoe-physics/src/lib.rs` (module export)
  - `crates/gutoe-physics/src/bin/global_gate_report.rs` (gate integration)

## Core formulas

- Standard-state binding free energy:
  - `ΔG = R T ln(Ki / 1 M)`
- QED contact energy per pair:
  - `E_pair = -(q1 q2) α ħ c / (ε r)`
  - Converted to kJ/mol via Avogadro scaling.

## Executed outputs (default lane)

Command:
- `cargo run -q -p gutoe-physics --bin cardiovascular_binding_report`

Observed:
- `Ki = 8.0 nM`
- `T = 298.15 K`
- `ΔG_exp = -46.217 kJ/mol`
- `QED ionic floor = -16.544 kJ/mol`
- `QED H-bond floor = -19.967 kJ/mol`
- `QED floor total = -36.512 kJ/mol`
- `Residual required = -9.706 kJ/mol`
- `Explained fraction = 0.790`

Residual decomposition (kJ/mol):
- Hydrophobic stabilization: `-10.350`
- Aromatic packing stabilization: `-2.750`
- Water-release stabilization: `-3.150`
- Conformational entropy penalty: `+4.680`
- Polar desolvation penalty: `+1.350`
- Ligand strain penalty: `+0.600`
- Modeled residual total: `-9.620`
- Residual closure error: `+0.086`

Artifacts:
- `/tmp/bh_renders/cardiovascular_binding/cardiovascular_binding_report.txt`
- `/tmp/bh_renders/cardiovascular_binding/cardiovascular_binding_report.json`

## CI gate

Command:
- `cargo run -q -p gutoe-physics --bin cardiovascular_binding_ci_gate`

Default windows:
- `40 <= |ΔG_exp| <= 55 kJ/mol`
- `|QED_floor| >= 20 kJ/mol`
- `explained_fraction >= 0.50`
- `|residual| <= 35 kJ/mol` and residual stabilizing (negative)
- `|residual_closure_error| <= 3 kJ/mol`

Result:
- `overall_pass = true`
- Artifact: `/tmp/bh_renders/cardiovascular_binding_ci_gate.json`

## Verification

- Unit tests:
  - `cargo test -q -p gutoe-physics --lib cardiovascular_binding`
  - Result: `3 passed, 0 failed`
- Global gate binary compiles with new lane:
  - `cargo check -q -p gutoe-physics --bin global_gate_report`

## Honesty statement

This is not full molecular-QM closure yet.  
It is a disciplined bridge:
- exact experimental thermodynamic anchor,
- first-principles electrostatic floor from the QED chain,
- explicit residual term that quantifies what remains to be derived,
- explicit non-electrostatic decomposition with an auditable closure error
  against the residual target.
