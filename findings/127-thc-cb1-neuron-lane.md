# THC-CB1 Neuron Lane — Focused Simulation Slice

Date: 2026-02-28  
Status: Implemented (report + CI gate)

## Goal
Run a tightly scoped lane for:
1. THC binding thermodynamics at CB1,
2. decomposition into QED floor vs non-electrostatic residual,
3. direct neuronal response curve (occupancy -> release/firing suppression).

## What shipped

- Module:
  - `crates/gutoe-physics/src/neuro_thc.rs`
- Binaries:
  - `crates/gutoe-physics/src/bin/thc_cb1_neuron_report.rs`
  - `crates/gutoe-physics/src/bin/thc_cb1_neuron_ci_gate.rs`
- Wiring:
  - `crates/gutoe-physics/src/lib.rs` (module export)

## Default run outputs

Command:
- `cargo run -q -p gutoe-physics --bin thc_cb1_neuron_report`

Binding decomposition (kJ/mol):
- `ΔG_exp = -43.927`
- `QED floor total = -5.239`
  - H-bond floor: `-3.706`
  - Polar floor: `-1.533`
- `Residual required = -38.688`
- `Explained fraction by QED floor = 0.119`

Residual breakdown (kJ/mol):
- Hydrophobic stabilization: `-36.400`
- Aromatic packing stabilization: `-4.650`
- Water-release stabilization: `-3.360`
- Conformational entropy penalty: `+3.300`
- Polar desolvation penalty: `+0.770`
- Ligand strain penalty: `+0.750`
- Modeled residual total: `-39.590`
- Residual closure error: `-0.902`

Neuronal response (default coupling) at selected concentrations:
- 10 nM:
  - occupancy `0.200`
  - release probability `0.321`
  - firing rate `7.604 Hz`
- 30 nM:
  - occupancy `0.429`
  - release probability `0.288`
  - firing rate `7.151 Hz`
- 100 nM:
  - occupancy `0.714`
  - release probability `0.247`
  - firing rate `6.586 Hz`

Artifacts:
- `/tmp/bh_renders/thc_cb1_neuron/thc_cb1_neuron_report.txt`
- `/tmp/bh_renders/thc_cb1_neuron/thc_cb1_neuron_report.json`
- `/tmp/bh_renders/thc_cb1_neuron/thc_cb1_neuron_sweep.csv`

## CI gate

Command:
- `cargo run -q -p gutoe-physics --bin thc_cb1_neuron_ci_gate`

Default windows:
- `38 <= |ΔG_exp| <= 50` kJ/mol
- `3 <= |QED_floor| <= 15` kJ/mol
- `|residual_closure_error| <= 5` kJ/mol
- monotone occupancy and monotone suppression curves

Result:
- `overall_pass = true`
- Artifact: `/tmp/bh_renders/thc_cb1_neuron_ci_gate.json`

## Honesty statement

This lane is a reduced biophysical transduction, not full receptor-level MD/QM closure.
It is explicitly decomposed and auditable: thermodynamic target, QED floor, residual
terms, and neuron-level response curve.
