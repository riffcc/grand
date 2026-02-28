# GRAND-363 (Gaia mode) — Gaia DR3 Life-Map Ingestion

Date: 2026-02-28
Status: Implemented (streaming ingest + map + target list)

## Goal
Replace synthetic-star-only mapping with a real Gaia DR3 ingestion lane that evaluates each star through:
1. derived habitability filter,
2. derived entropy progression stage timeline,
3. local Kauffman closure check scaled by measured metallicity.

## What shipped

### New binary
- `crates/gutoe-physics/src/bin/gaia_dr3_life_map.rs`

### Shared physics wiring
- `crates/gutoe-physics/src/galactic_life_map.rs`
  - Exported shared functions used by Gaia lane:
    - `derive_thresholds_and_multipliers`
    - `classify_stage`
    - `stage_entropy_multiplier`
    - `main_sequence_lifetime_gyr`
    - `habitability_score`
    - `is_habitable`
    - `infer_component_from_position`
- `crates/gutoe-physics/src/lib.rs`
  - Added module/export for `galactic_life_map`.

## Gaia ingest behavior
- Streaming CSV ingest (constant-memory, line-by-line processing).
- Accepts either:
  - direct Cartesian positions (`x_ly,y_ly,z_ly` or `x_pc,y_pc,z_pc`), or
  - sky coordinates (`ra,dec,parallax` or `distance_pc`) with ICRS->Galactic transform.
- Uses measured/derived stellar fields:
  - age (`age_gyr` / `age_years` / `log_age_years`)
  - metallicity (`metallicity` / `mh_gspphot` / `feh`)
  - mass (`mass_solar`) or Teff proxy (`teff`) if mass missing.
- Applies local Kauffman gate per star:
  - `N*p_local = N*p_baseline * sqrt(Z/Z_solar)`
  - with `Z_solar = 0.0142` and threshold `N*p >= 1`.
- Emits nearest signal targets as concrete star IDs + coordinates.

## Artifacts emitted
Default output dir: `/tmp/bh_renders/gaia_life_map`

- `gaia_life_map.png`
- `gaia_life_report.txt`
- `gaia_life_report.json`
- `gaia_signal_targets.csv`
- optional (if `GUTOE_GAIA_WRITE_SIGNAL_CATALOG=1`): `gaia_signal_catalog.csv`

## Runtime controls
- `GUTOE_GAIA_DR3_CSV` (required): path to Gaia CSV.
- `GUTOE_GAIA_MAX_ROWS` (optional): row cap for smoke runs.
- `GUTOE_GAIA_RENDER_SAMPLE` (optional): reservoir size for map rendering.
- `GUTOE_GAIA_NEAREST_K` (optional): number of nearest signal targets to keep.
- `GUTOE_GAIA_WRITE_SIGNAL_CATALOG` (optional): `1/true` to write all signal rows.

## Smoke verification
Run on a 10-row mock Gaia CSV:

- Command:
  - `GUTOE_GAIA_DR3_CSV=/tmp/gaia_mock.csv GUTOE_GAIA_MAX_ROWS=1000 GUTOE_GAIA_RENDER_SAMPLE=10000 cargo run -q -p gutoe-physics --bin gaia_dr3_life_map`
- Result:
  - `seen=10`, `used=10`, `habitable_now=7`, `signal_now=0`
  - all expected artifacts were written.

## Real Gaia sample render (executed)

Fetched real Gaia DR3 sample slices from ESA TAP and rendered immediately:

- Source sample: `/tmp/gaia_dr3_sample_ra_bins.csv` (18,000 real rows)
- Render command:
  - `GUTOE_GAIA_DR3_CSV=/tmp/gaia_dr3_sample_ra_bins.csv GUTOE_GAIA_RENDER_SAMPLE=18000 GUTOE_GAIA_NEAREST_K=128 cargo run -q -p gutoe-physics --bin gaia_dr3_life_map`
- Output:
  - `/tmp/bh_renders/gaia_life_map/gaia_life_map.png`
  - `/tmp/bh_renders/gaia_life_map/gaia_signal_targets.csv`
  - `/tmp/bh_renders/gaia_life_map/gaia_life_report.json`
- Run summary:
  - `rows_seen=18000`
  - `rows_used=18000`
  - `habitable_count_present=8247`
  - `signal_count_present=156`
  - `signal_fraction_present=0.008666666666666666`

## Readiness statement
The ingestion lane is now implemented and can process Gaia DR3 exports immediately.
Remaining operational step is data logistics: providing the full Gaia DR3 CSV export (or chunked exports) with the required columns.
