# 167 — Stable-Presence + Chemical Extraction (Z<=140)

Date: 2026-02-28

## Scope
Produce a unified extraction table with all currently derivable per-element chemical details (mass/volume/density/phase/thermodynamics), joined to stable-presence classifier outputs, and extend projection beyond Z=94.

## Code updates
- `crates/gutoe-physics/src/bin/mass_periodic_report.rs`
  - Added `stable_like_isotope_aspects.csv` export (full isotope-level gate + weak-channel aspects for every predicted stable-like nuclide).
  - Added `element_stable_presence_aspects.csv` export (per-element stable-presence summary including best isotope and weak-Q margin stats).
- `crates/gutoe-physics/src/bin/chemical_thermo_report.rs`
  - Added configurable `GUTOE_CHEM_Z_MAX` (default 140).
  - Added explicit STP phase field (`state_273k_1bar_stp`).
  - Added mass/volume detail fields:
    - `atomic_mass_u`
    - `mass_per_atom_kg`
    - `volume_per_atom_ang3`
    - `molar_volume_stp_l_mol_if_gas`

## Commands run
```bash
cargo run -q -p gutoe-physics --bin mass_periodic_report
GUTOE_CHEM_Z_MAX=140 cargo run -q -p gutoe-physics --bin chemical_thermo_report
```

## Output artifacts
- `/tmp/nuclear_chart/stable_like_isotope_aspects.csv`
- `/tmp/nuclear_chart/element_stable_presence_aspects.csv`
- `/tmp/bh_renders/chemical_thermodynamics/chemical_thermo_report.csv`
- `/tmp/bh_renders/chemical_thermodynamics/chemical_thermo_report.json`
- `/tmp/nuclear_chart/element_chem_nuclear_full.csv` (merged)
- `/tmp/nuclear_chart/stable_presence_elements_chem_detail.csv` (merged, predicted stable-presence only)
- `/tmp/nuclear_chart/post94_chem_projection.csv` (merged, Z>94)

## Key counts
- Elements modeled (chemical lane): 140
- Predicted stable-presence elements: 80
- Predicted stable-presence elements for Z>94: 0
- Rows in merged full table: 140 (+ header)
- Rows in stable-presence-only table: 80 (+ header)
- Rows in post-94 table: 46 (+ header)

## Notes
- Chemical lane is proxy thermodynamics; values are explicit transduction outputs, not ab-initio quantum chemistry.
- Nuclear stable-presence beyond Z=94 remains false across Z=95..140 under current classifier gates.
