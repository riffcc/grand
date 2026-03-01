# Finding 232 — Executed Next3 (T07, T11, T13)

## Summary

Executed the three requested tests:

- `T07` DESI void scalar anomaly lane
- `T11` chronological weak-angle holdout
- `T13` scheme-coherent weak-angle fit

Binary:

- `crates/gutoe-physics/src/bin/ctc_next3_t07_t11_t13.rs`

Outputs:

- `/tmp/bh_renders/ctc_next3_t07_t11_t13/ctc_next3_t07_t11_t13.txt`
- `/tmp/bh_renders/ctc_next3_t07_t11_t13/ctc_next3_t07_t11_t13.json`

## Results

1. `T13` **PASS**
   - metric: `reduced_chi2_scheme_clean`
   - value: `0.278100`
   - threshold: pass if `<= 1.5`

2. `T11` **PASS**
   - metric: `holdout_reduced_chi2`
   - value: `0.439208`
   - threshold: pass if `<= 1.5`

3. `T07` **FAIL**
   - metric: `z_void_scalar`
   - value: `0.000000`
   - threshold: fail if `< 3.0`

## Harness Integration

`ctc_falsification_20_harness` now ingests the Next3 artifact and updates test
statuses automatically.

Updated global scoreboard:

- `PASS = 5`
- `FAIL = 2`
- `OPEN = 13`
- `TOTAL = 20`

Output:

- `/tmp/bh_renders/ctc_falsification_20/ctc_falsification_20.json`
