# GRAND Rung: Unified Element Algebra Table

Date: 2026-02-28

## What Landed

A single element-level artifact was added to `mass_periodic_report` that fuses:

- nuclear stability classification and shell/beta-derived statistics,
- atomic SCF observables,
- thermodynamic properties from the coupled chemistry lane,
- coupling diagnostics and deltas vs proxy chemistry.

New outputs:

- `/tmp/nuclear_chart/element_unified_algebra_table.csv`
- `/tmp/nuclear_chart/element_unified_algebra_summary.json`

## Implementation

File updated:

- `crates/gutoe-physics/src/bin/mass_periodic_report.rs`

Key additions:

- Imports for atomic SCF + coupled thermo (`predict_atomic_scf`, `predict_element_thermo`, `predict_element_thermo_coupled_with_diagnostics`), and phase/state/constant helpers.
- `representative_mass_number_for_z(...)` helper:
  - prefers observed-stable isotopes (`Z <= 94`) nearest to the structural mass target,
  - falls back to best predicted stable-like isotope,
  - then falls back to structural target if needed.
- Unified table writer with one row per element (`Z_min..Z_max`) including:
  - stable-like presence/stats,
  - SCF energies and descriptors,
  - coupled thermo state/mass/volume/transition properties,
  - coupling deltas and diagnostics.
- Unified summary JSON with saturation/state/delta metrics.
- CLI output now explicitly prints unified artifact paths.

## Verification

Commands run:

- `cargo check -q -p gutoe-physics --bin mass_periodic_report` (pass)
- `cargo run -q -p gutoe-physics --bin mass_periodic_report` (pass)

Representative output metrics (`Z<=140`, coupled mode on):

- rows: `140`
- density clamp saturation count: `39`
- state counts @ 298 K, 1 bar: `solid=119, liquid=5, gas=16`
- mean abs delta vs proxy:
  - density: `5.445865417 g/cm^3`
  - melting: `153.384001838 K`
  - boiling: `405.155853969 K`

## Honest Status

This is a cross-domain fused artifact from one pipeline, not final high-accuracy materials closure. The immediate next tightening step is calibration of the remaining high-density tail and transition-temperature spread while preserving the coupled-mode gains and ambient phase realism.
