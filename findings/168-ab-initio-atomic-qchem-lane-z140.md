# 168 — Ab-initio Atomic Q-Chem Lane (Z<=140)

Date: 2026-02-28

## What landed
- New module: `crates/gutoe-physics/src/ab_initio_qchem.rs`
  - Self-consistent atomic SCF solver (spherical effective potential + screening).
  - Madelung orbital filling up to heavy-Z range.
  - Frontier descriptors derived from solved orbitals:
    - HOMO/LUMO energies
    - ionization-energy proxy
    - electron-affinity proxy
    - Mulliken electronegativity
    - hardness/softness
    - atomic/covalent radius proxies
    - polarizability proxy
- New report binary: `crates/gutoe-physics/src/bin/ab_initio_qchem_report.rs`
  - Emits element and orbital reports for configurable `Z` max (`GUTOE_ABINITIO_Z_MAX`, default 140).
- Library registration:
  - `crates/gutoe-physics/src/lib.rs` exports `ab_initio_qchem`.

## Verification
- Build check:
  - `cargo check -q -p gutoe-physics --bin ab_initio_qchem_report`
- Unit tests (library-only lane tests):
  - `cargo test -q -p gutoe-physics --lib ab_initio_qchem`
  - Result: 3 passed, 0 failed.

## Run
```bash
GUTOE_ABINITIO_Z_MAX=140 cargo run -q -p gutoe-physics --bin ab_initio_qchem_report
```

## Outputs
- `/tmp/bh_renders/ab_initio_qchem/ab_initio_qchem_report.csv`
- `/tmp/bh_renders/ab_initio_qchem/ab_initio_qchem_orbitals.csv`
- `/tmp/bh_renders/ab_initio_qchem/ab_initio_qchem_report.json`
- `/tmp/bh_renders/ab_initio_qchem/ab_initio_qchem_report.txt`

## Summary metrics (report JSON)
- `elements_modeled = 140`
- `mean_ionization_energy_ev = 22.6520374828`
- `mean_electronegativity_ev = 13.1200504491`
- `mean_atomic_radius_pm = 434.858087255`

## Integrated master table
- Three-way merged table (ab-initio + chemical thermo + stable-presence):
  - `/tmp/nuclear_chart/ab_initio_chemistry_master_z140.csv`
  - Rows: 141 (header + 140 elements).

## Notes
- This lane is atomic SCF first-principles style; it is not yet full molecular ab-initio chemistry (no multi-center basis, no post-HF correlation, no explicit molecular geometry optimization).
- It establishes the electron-structure foundation needed to replace remaining reduced-order chemistry transductions lane by lane.
