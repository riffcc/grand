# GRAND Rung: Chemical Thermodynamics Coupled to Atomic SCF

Date: 2026-02-28

## What Landed

Implemented a native coupling from the atomic SCF lane into the periodic thermodynamics lane:

- Added `predict_element_thermo_coupled_with_diagnostics(z, a)` and `predict_element_thermo_coupled(z, a)` in `crates/gutoe-physics/src/chemical_thermo.rs`.
- Added `CoupledThermoDiagnostics` to expose SCF-side constraints per element:
  - SCF iterations/residual
  - SCF valence electron count
  - SCF atomic radius
  - SCF IE / EA / Mulliken electronegativity / hardness
  - coupled radius and frontier cohesive proxy
- Refactored shared thermodynamic assembly into `assemble_element_thermo(...)`.
- Added/kept tests for both proxy and coupled lanes (`cargo test --lib chemical_thermo`: 8 passed).

Report wiring (`crates/gutoe-physics/src/bin/chemical_thermo_report.rs`):

- New env toggle: `GUTOE_CHEM_THERMO_COUPLED` (default `true`).
- Emits extended CSV/JSON diagnostics and deltas vs proxy.
- Adds STP columns and summary delta metrics.
- Supports dynamic `z_max` (default 140).

## Tuned Coupling Behavior

Initial coupling overcorrected volatility (ambient gases collapsed to 0). Tuned with:

- family-specific frontier blending weights,
- bounded cohesive shift (`0.70x .. 1.35x` baseline),
- bounded radius shift (`0.70x .. 1.35x` baseline),
- softer hardness/valence gates.

## Verification

- `cargo check -q -p gutoe-physics --bin chemical_thermo_report`: pass
- `cargo test -q -p gutoe-physics --lib chemical_thermo`: pass (8/8)
- `GUTOE_CHEM_THERMO_COUPLED=1 cargo run -q -p gutoe-physics --bin chemical_thermo_report`: pass

Output snapshots:

- Proxy baseline: `/tmp/bh_renders/chemical_thermodynamics/chemical_thermo_report_proxy.{txt,json,csv}`
- Tuned coupled: `/tmp/bh_renders/chemical_thermodynamics/chemical_thermo_report_coupled_tuned.{txt,json,csv}`

## A/B (Proxy vs Tuned Coupled, Z <= 140)

- Density clamp saturation (`density >= 39.9 g/cm^3`):
  - proxy: 65
  - tuned coupled: 39
- Mean density (g/cm^3):
  - proxy: 23.2887834
  - tuned coupled: 17.8660862143
- Ambient state counts at 298 K:
  - proxy: solid=118, liquid=6, gas=16
  - tuned coupled: solid=119, liquid=5, gas=16
- Mean absolute shift vs proxy:
  - density: 5.4226971857 g/cm^3
  - melting: 153.3840017571 K
  - boiling: 405.1558541643 K

## Honest Status

This is a tighter, cross-lane constrained thermodynamics model (SCF + bulk transduction), not full many-body condensed-matter closure. It improves clamp behavior while preserving ambient gas count. Remaining gap: reduce residual density saturation further and calibrate transition-temperature spread without introducing fitted element tables.
