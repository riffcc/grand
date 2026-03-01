# Finding 246 — Fixed `s²` Structural-Fraction Test (`13/12`, `3/8`)

## Hypothesis
Use structural fixed amplitudes:
- `s²_up = 2 + 13/12 = 3.083333...`
- `s²_down = 2 + 3/8 = 2.375`

while keeping fitted sector phases `δ_u`, `δ_d` (and `M_u`, `M_d`) from the existing Z3 extraction.

Question: does this collapse the down-sector residual (`mb/ms` lane)?

## Baseline (fitted `s`)
- `s²_up = 3.094038784751`
- `s²_down = 2.390533909547`
- `ms/md = 19.914346895075` (err vs 19: `4.812%`)
- `mb/ms = 44.946236559140` (err vs 51.434343...: `12.614%`)
- Down-only RMS-log closure: `0.100972056315`

## Fixed-`s²` test result
- `ms/md = 18.784999593029` (err vs 19: `1.132%`)  ← improves
- `mb/ms = 44.294318499813` (err vs 51.434343...: `13.882%`)  ← worsens
- Down-only RMS-log closure: `0.105982833397`  ← worse overall

## Structural-ratio impact (full 7-ratio check)
Max relative error against structural ratio targets:
- Baseline: `14.370%`
- Fixed-`s²`: `26.875%`

Largest degradations are cross/up-linked ratios (`m_u/m_d`, `m_c/m_u`).

## Verdict
This hypothesis is **not** the fix. It does not collapse the down-sector residual and worsens overall structural-ratio consistency.

## Notes
- Cycle equalities computed from any single internally generated mass set are tautological; the discriminating metric is mismatch vs structural targets.
