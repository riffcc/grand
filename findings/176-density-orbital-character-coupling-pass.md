# 176 — Density Tail Tightening via Orbital-Character Coupling

## Intent
Apply a **general, non-element-specific** correction for density tails using SCF-derived orbital character:
- p-dominant frontier character -> lower effective packing / larger effective volume
- open d-shell heavy metals -> stronger compaction
- f-core (period >= 6) -> lanthanide/actinide contraction-style radius correction

No element-name branching introduced in this pass.

## Code changes
Primary file:
- `crates/gutoe-physics/src/chemical_thermo.rs`

Additions/updates:
- `OrbitalPackingHints` derived from SCF orbitals and HOMO proximity weighting.
- `packing_fraction_from_hints(...)` for orbital-character-dependent packing.
- Coupled radius correction factors driven by `p_frac`, `open_d_shell`, `closed_d_shell`, `has_f_core`, and f-block family.
- Dynamic lower clamp on coupled radius for actinides/5d-with-f-core to avoid artificial floor at `0.70 * base_radius`.
- Pass coupled packing fraction through `assemble_element_thermo(...)`.
- Added `coupled_packing_fraction` to `CoupledThermoDiagnostics`.

## Verification
- `cargo test -q -p gutoe-physics --lib chemical_thermo`
  - Result: 8 passed, 0 failed.

## External benchmark (Z=1..94, phase override ON)
Baseline before this pass (latest prior checkpoint):
- phase accuracy: `1.000000`
- density MAE: `5.006358`
- density MAE (condensed only): `5.677915`

After first orbital-character pass:
- phase accuracy: `1.000000`
- density MAE: `4.267654`
- density MAE (condensed only): `4.840117`

After lanthanide/actinide clamp tuning pass:
- phase accuracy: `1.000000`
- density MAE: `3.871942`
- density MAE (state-aware): `3.871927`
- density MAE (condensed only): `4.391322`
- density red count: `74` (from `77` earlier)

## Current dominant residuals (top abs density errors)
- Underpredicted 4d/5d transition cluster: Mo/Tc/Ru/Rh/Pd/Nb (too low)
- Underpredicted high-density 5d+actinide tail: Os/Ir/Re/Np (still low, but improved)
- Overpredicted directional covalent solids: C, S (still high)

## Honest read
This pass materially improved density fidelity while preserving phase closure and without introducing element-list branching. The remaining error structure suggests the next tightening should separate:
1) open d-band metallic compaction (4d vs 5d distinction), and
2) directional covalent network porosity for light p-block solids.
