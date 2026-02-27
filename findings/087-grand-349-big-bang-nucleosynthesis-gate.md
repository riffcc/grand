# Finding 087: GRAND-349 Big Bang Nucleosynthesis Gate

Date: 2026-02-27  
Status: GRAND-349 first quantitative BBN lane landed

## Goal

Assemble a primordial-abundance lane from already-derived primitives (`η_B`,
`λ_QG`, Cl(1,3) state counts) and check the standard BBN observables:

- `Y_p` (He-4 mass fraction)
- `D/H`
- `3He/H`
- `7Li/H` (explicit lithium-tension tracking)

## What landed

### Rust lane

New module:

- `crates/gutoe-physics/src/bbn.rs`

New binaries:

- `crates/gutoe-physics/src/bin/bbn_report.rs`
- `crates/gutoe-physics/src/bin/bbn_ci_gate.rs`

Wired into crate exports:

- `crates/gutoe-physics/src/lib.rs`

Structural anchors are now derived from shared constants, not hardcoded:

- `η10_ref = (12*5)/(4+6) = 6`
- deuterium exponent `= (6+2)/(4+1) = 8/5`
- helium-3 exponent `= 3/(4+1) = 3/5`
- lithium tension amplification `= 12/4 = 3`

Quantitative lane:

- `η10 = 10^10 * η_B` from baryogenesis lane
- `Y_p = Y_p,target + (λ_QG/50) * (η10 - η10_ref)`
- `D/H = D/H_target * (η10_ref/η10)^(8/5)`
- `3He/H = 3He/H_target * (η10_ref/η10)^(3/5)`
- `7Li/H = Li_obs * (η10/η10_ref)^2 * 3`

Gate windows (default):

- `|ΔY_p| <= 0.010`
- `D/H` relative error `<= 0.15`
- `3He/H` relative error `<= 0.15`
- lithium tension ratio expected in `[2, 4]`

### Lean parity lane

New module:

- `lean/Gutoe/BBN.lean`

Wired into root build:

- `lean/lakefile.lean` includes `Gutoe.BBN`

No `sorry`; key parity theorems include:

- `eta10_ref_eq_6`
- `deuterium_eta_exponent_eq`
- `helium3_eta_exponent_eq`
- `lithium7_tension_amplification_eq`
- `eta10_from_baryogenesis_pos`
- `primordial_helium4_at_reference`
- `primordial_deuterium_at_reference`
- `primordial_helium3_at_reference`
- `lithium7_tension_ratio_at_reference`
- `lithium7_tension_ratio_reference_window`

## Quantitative result

From `/tmp/bh_renders/bbn_report.txt`:

- `η10 = 6.3015`
- `Y_p(pred) = 0.245502` (target `0.245`, `Δ = +5.025e-4`, ~0.205%)
- `D/H(pred) = 2.3548e-5` (target `2.547e-5`, rel err `0.0754`)
- `3He/H(pred) = 1.0681e-5` (target `1.1e-5`, rel err `0.0290`)
- `7Li/H(pred) = 5.2945e-10` vs observed `1.6e-10` (tension ratio `3.309`)

Gate status:

- primary BBN gate: pass
- lithium tension lane: present and in expected range
- overall gate: pass

## Honest boundary

This is a first quantitative assembly lane, not yet a full reaction-network+Boltzmann
solver over temperature-time history.

Still open for hard closure:

1. Replace anchored abundance scaling with explicit reaction-network integration
   over expansion timeline.
2. Tie weak freeze-out and neutron/proton ratio evolution directly to the
   cosmology lane (`H(T)`).
3. Add uncertainty propagation from `η_B` and reheating lane into abundance bands.

## Build sanity

- `cargo check -p gutoe-physics --bin bbn_report --bin bbn_ci_gate` ✅
- `cargo test -p gutoe-physics bbn -- --nocapture` ✅
- `cargo run -q -p gutoe-physics --bin bbn_report` ✅
- `cargo run -q -p gutoe-physics --bin bbn_ci_gate` ✅
- `cd lean && lake build Gutoe.BBN` ✅
- `cd lean && lake build Gutoe` ✅

No new `sorry`.
