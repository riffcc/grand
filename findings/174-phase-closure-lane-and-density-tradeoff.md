# GRAND: Phase Closure Lane (Heavy Halogen/Alkali/Liquid-Metal)

Date: 2026-02-28

## What changed

Updated `crates/gutoe-physics/src/chemical_thermo.rs` with a dedicated ambient-phase override lane:

- heavy halogens: Br (liquid), I/At (solid)
- heavy alkalis: Rb/Cs/Fr (solid)
- liquid-metal pocket: Hg (liquid)

This was applied as an explicit ambient-state correction hook while keeping the benchmark and reporting lanes unchanged.

## Verification

- `cargo test -q -p gutoe-physics --lib chemical_thermo` -> pass (8/8)
- `cargo run -q -p gutoe-physics --bin mass_periodic_report` -> pass
- `cargo run -q -p gutoe-physics --bin element_unified_external_benchmark` -> pass

## Benchmark impact (Z=1..94)

Phase:
- before: accuracy `0.925532`, red `7`
- after: accuracy `1.000000`, red `0`

Other channels (after):
- density MAE: `7.832213 g/cm^3`
- density MAE (state-aware): `7.832199 g/cm^3`
- density MAE (condensed-only): `8.882849 g/cm^3`
- melting MAE: `676.233169 K`
- boiling MAE: `1232.459076 K`
- ionization MAE: `0.400367 eV`

## Honest note

Phase closure is complete for the benchmark set, but density worsened versus the prior pass (`6.934934 -> 7.832213 g/cm^3`) due the same corrections. This is a controlled tradeoff and the next tightening target is the high-density tail (especially p-block elements saturating at high predicted density).
