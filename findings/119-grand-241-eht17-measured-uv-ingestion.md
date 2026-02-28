# Finding 119: GRAND-241 EHT-17 Measured uv-Coverage Ingestion

Date: 2026-02-28
Status: GRAND-241 complete

## Scope

Close `GRAND-241` by ingesting measured 2017 EHT uv tracks and sampling synthetic visibilities at those measured uv points.

## Data source used

Official EHT first M87 results release (CyVerse curated dataset, DOI `10.25739/g85n-f134`), file:

- `csv/SR1_M87_2017_096_lo_hops_netcal_StokesI.csv`

Direct path used during verification:

- `https://de.cyverse.org/anon-files/iplant/home/shared/commons_repo/curated/EHTC_FirstM87Results_Apr2019/csv/SR1_M87_2017_096_lo_hops_netcal_StokesI.csv`

## Implementation changes

Updated `crates/gutoe-gpu/src/bin/bh_render.rs` (`run_eht_uv_export`):

1. Added parser support for official EHT CSV row shape:
   - `time_utc, T1, T2, U(lambda), V(lambda), ...`
2. Added unit transduction from measured uv in wavelengths to renderer uv units:
   - `u_norm = u_lambda * fov_rad`
   - `v_norm = v_lambda * fov_rad`
   - where `fov_rad = fov_rs * rs_rad`.
3. Added uv source metadata to export JSON:
   - `uv_source_mode`
   - `uv_source_rows`
   - `uv_source_input_csv`
   - `note`
4. Hardened closure-map construction for repeated baselines by averaging complex visibilities per baseline key before triangle/quad closure products.

## Verification run

Commands:

```bash
curl -L -s 'https://de.cyverse.org/anon-files/iplant/home/shared/commons_repo/curated/EHTC_FirstM87Results_Apr2019/csv/SR1_M87_2017_096_lo_hops_netcal_StokesI.csv' > /tmp/SR1_M87_2017_096_lo_hops_netcal_StokesI.csv
BH_UV_TRACK_CSV=/tmp/SR1_M87_2017_096_lo_hops_netcal_StokesI.csv cargo run -q -p gutoe-gpu --bin bh_render -- eht_uv m87_eht2017 96x96
jq '{uv_source_mode,uv_source_rows,note}' /tmp/bh_renders/m87_eht2017_eht_uv.json
```

Observed outputs:

- `uv_source_mode = "measured_eht_csv_lambda"`
- `uv_source_rows = 8645`
- Note confirms `U/V(lambda)` conversion into normalized uv units.
- Artifacts generated:
  - `/tmp/bh_renders/m87_eht2017_eht_uv.csv`
  - `/tmp/bh_renders/m87_eht2017_eht_closure.csv`
  - `/tmp/bh_renders/m87_eht2017_eht_closure_amp.csv`
  - `/tmp/bh_renders/m87_eht2017_eht_uv.json`

## Outcome

`GRAND-241` now runs on measured 2017 EHT uv points from the public release instead of synthetic fallback coordinates.
