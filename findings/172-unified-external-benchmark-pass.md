# GRAND Rung: Unified External Benchmark Pass

Date: 2026-02-28

## What Landed

Added an external benchmark scorer for the unified element artifact.

New files:

- `crates/gutoe-physics/src/bin/element_unified_external_benchmark.rs`
- `crates/gutoe-physics/data/periodic_pubchem_reference.csv`

Outputs:

- `/tmp/nuclear_chart/element_unified_external_benchmark.csv`
- `/tmp/nuclear_chart/element_unified_external_benchmark.json`
- `/tmp/nuclear_chart/element_unified_external_benchmark.txt`

## Scope

Benchmark compares unified table predictions to external periodic references over configurable `Z` range (default `1..94`):

- phase at 298 K
- density (g/cm^3)
- melting temperature (K)
- boiling temperature (K)
- first ionization energy (eV)

It emits per-element error rows and red-light flags with aggregate MAE/MAPE and counts.

## Latest Metrics (Z=1..94)

From `/tmp/nuclear_chart/element_unified_external_benchmark.txt`:

- phase accuracy: `0.829787` (n=94, red=16)
- density MAE: `12.296632 g/cm^3` (MAPE `300750.520964%`, n=93, red=83)
- density MAE condensed-only (solid+liquid refs): `9.696637 g/cm^3` (MAPE `162.248486%`, n=82)
- melting MAE: `751.354342 K` (MAPE `238.594202%`, n=94, red=82)
- boiling MAE: `1368.423961 K` (MAPE `356.830202%`, n=91, red=72)
- ionization MAE: `9.796204 eV` (MAPE `115.161240%`, n=94, red=85)
- elements with any red: `94`

## Interpretation

The benchmark pass is now fully lit and explicit. Current chemistry/scf surfaces are still far from external numeric closure under this strict reference comparison. This gives a concrete red-light roadmap for tightening:

1. low-Z phase and volatility (H/alkali/noble-gas regime)
2. condensed-phase density scale calibration
3. transition temperatures (melt/boil transduction)
4. SCF ionization-energy scaling

## Notes

Reference source in this pass is PubChem periodic table JSON export (stored locally as CSV). If needed, this benchmark lane can be pointed at CRC/NIST-curated reference tables via `GUTOE_REFERENCE_TABLE` without changing scorer logic.
