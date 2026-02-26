# gutoe-physics data snapshots

## `sparc_massmodels_2016c_baryon.csv`

Source:
- SPARC database (Lelli, McGaugh, Schombert 2016c)
- URL: `https://astroweb.cwru.edu/SPARC/MassModels_Lelli2016c.mrt`

Snapshot notes:
- Extracted on `2026-02-26`
- Rows retained: entries with `Vobs > 0` and positive baryonic speed
- Derived column:
  - `v_baryon_kms = sqrt(v_gas_kms^2 + v_disk_kms^2 + v_bulge_kms^2)`

This snapshot is used by `dark_matter_falsification.rs` for GRAND-346
dataset-backed rotation/lensing-proxy scoring and falsification gates.
