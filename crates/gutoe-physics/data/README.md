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

## `COM_PowerSpect_CMB-TT-binned_R3.01.txt` (GRAND-355)

Source target:
- Planck 2018 TT binned spectrum
- URL:
  - `https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/COM_PowerSpect_CMB-TT-binned_R3.01.txt`

Schema (whitespace):
- `l Dl -dDl +dDl BestFit`

Usage:
- `cmb_class_report` reads this path from `GUTOE_PLANCK_TT`.
- Default path is:
  - `crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt`

Notes:
- This lane compares full-shape TT residuals over `ell=2..2500` against
  CLASS output (CAMB-format `D_ell` in micro-Kelvin squared).
- The parser also accepts CSV files with:
  - `ell,d_ell_tt_uk2,sigma_uk2`

## `COM_PowerSpect_CMB-TT-full_R3.01.txt` (GRAND-355)

Source target:
- Planck 2018 TT full (unbinned) spectrum
- URL:
  - `https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/COM_PowerSpect_CMB-TT-full_R3.01.txt`

Schema (whitespace):
- `l Dl -dDl +dDl`

Usage:
- Used by:
  - `cmb_three_way_compare`
  - `cmb_likelihood_scan` (when `GUTOE_PLANCK_TT` points to the full file)
- Range used by the current CLASS harness:
  - `ell = 2..2500`

Notes:
- Three-way consistency check compares:
  - prediction vs binned
  - prediction vs full
  - binned vs rebinned-from-full
