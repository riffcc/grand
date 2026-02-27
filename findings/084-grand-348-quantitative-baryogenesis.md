# Finding 084: GRAND-348 Quantitative Baryogenesis Gate

Date: 2026-02-27  
Status: GRAND-348 initial quantitative lane landed

## Goal

Turn the baryogenesis lane into a falsifiable quantitative check against the
baryon-to-photon target `η_B ≈ 6.12e-10`, using only already-derived primitives.

## What landed

### Rust lane

New module:

- `crates/gutoe-physics/src/baryogenesis.rs`

New binaries:

- `crates/gutoe-physics/src/bin/baryogenesis_report.rs`
- `crates/gutoe-physics/src/bin/baryogenesis_ci_gate.rs`

Core quantitative expression:

- `η_B = J * α_LO² * (5/11) * (1 - λ_QG) * (486/485)`

Where all factors are already shared:

- `J` from CKM texture diagonalization (`ckm_from_textures`)
- `α_LO = 1/137`
- `5/11` from dark/visible Clifford split
- `λ_QG = 1/12`
- `486/485` finite-mode rescale

Sakharov check wiring:

- CP violation: from `cp_violation_witness(...)`
- B-violation channel witness: weak non-Abelian lane present (`sin²θ_W` physical)
- non-equilibrium witness: structural survival factor in `(0, 1)`

Hard gate:

- `cargo run -q -p gutoe-physics --bin baryogenesis_ci_gate`
- emits `/tmp/bh_renders/baryogenesis_ci_gate.json`
- exits nonzero on failure

### Lean parity lane

New module:

- `lean/Gutoe/Baryogenesis.lean`

Wired into build:

- `lean/lakefile.lean` now includes `Gutoe.Baryogenesis`

Theorems (no `sorry`):

- `baryo_micro_mode_count_eq_486`
- `baryo_micro_finite_rescale_eq`
- `baryo_nonequilibrium_survival_eq`
- `baryo_nonequilibrium_survival_bounds`
- `baryogenesis_prefactor_eq`
- `baryogenesis_prefactor_pos`
- `eta_baryon_structural_pos`
- `eta_baryon_structural_from_shared_primitives`

## Quantitative result

From `/tmp/bh_renders/baryogenesis_report.txt`:

- `η_B(pred) = 6.301488e-10`
- `η_B(obs)  = 6.120000e-10`
- relative error `= 0.0297` (2.97%)
- default gate window `η_rel_error_max = 0.15` → pass
- Sakharov checks in this lane: all pass

## Honest boundary

This is a **first quantitative closure lane**, not full cosmological
baryogenesis dynamics yet.

Still open for GRAND-348 hard closure:

1. Replace structural B-violation witness with explicit sphaleron-rate lane
   coupled to temperature history.
2. Replace static non-equilibrium witness with Boltzmann/Friedmann freeze-out
   evolution.
3. Tighten the current 15% gate window once dynamical washout is in place.

## Build sanity

- `cargo check -p gutoe-physics --bin baryogenesis_report --bin baryogenesis_ci_gate` ✅
- `cargo test -p gutoe-physics baryogenesis -- --nocapture` ✅
- `cargo run -q -p gutoe-physics --bin baryogenesis_report` ✅
- `cargo run -q -p gutoe-physics --bin baryogenesis_ci_gate` ✅
- `cd lean && lake build Gutoe.Baryogenesis` ✅
- `cd lean && lake build Gutoe` ✅

No new `sorry`.
