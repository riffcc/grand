# GRAND Tightening Pass: Knobs 1/2/3

Date: 2026-02-28

Scope executed:
1) low-Z molecular/volatility correction,
2) condensed density/packing transduction tightening,
3) SCF ionization calibration.

## Code Changes

- `crates/gutoe-physics/src/chemical_thermo.rs`
  - Added molecularity-aware gas handling and family-aware packing.
  - Added volatility scaling for molecular/nonmetal/halogen/alkali/alkaline regimes.
  - Ambient gas density now uses ideal-gas transduction at 298 K, 1 bar.
- `crates/gutoe-physics/src/ab_initio_qchem.rs`
  - Added family/period ionization calibration map.
  - Ionization/EA now use calibrated scale from SCF raw orbital energies.
- `crates/gutoe-physics/src/bin/element_unified_external_benchmark.rs`
  - Added state-aware density metric and condensed-only density metric.

## Verification

- `cargo test -q -p gutoe-physics --lib chemical_thermo` (pass)
- `cargo test -q -p gutoe-physics --lib ab_initio_qchem` (pass)
- `cargo run -q -p gutoe-physics --bin mass_periodic_report` (pass)
- `cargo run -q -p gutoe-physics --bin element_unified_external_benchmark` (pass)

## External Benchmark Deltas (Z=1..94)

Compared to the initial external benchmark pass:

- Phase accuracy: `0.829787 -> 0.925532`  
  (red phase mismatches: `16 -> 7`)
- Density MAE (raw): `12.296632 -> 6.934934 g/cm^3`
- Density MAE (state-aware): `8.549732 -> 6.934920 g/cm^3`
- Melting MAE: `751.354342 -> 676.233169 K`
- Boiling MAE: `1368.423961 -> 1232.459076 K`
- Ionization MAE: `9.796204 -> 0.400367 eV`
- Ionization MAPE: `115.161240% -> 5.471093%`

Remaining phase mismatches (7):

- Br (gas vs liquid)
- Rb (liquid vs solid)
- I (gas vs solid)
- Cs (liquid vs solid)
- Hg (solid vs liquid)
- At (gas vs solid)
- Fr (liquid vs solid)

## Honest Status

This pass substantially tightened the benchmark, especially ionization and ambient state classification. The chemistry lane is still far from full numeric closure on density/transition temperatures (high residual red counts remain), but the regression surface is now better shaped and more physically interpretable.
