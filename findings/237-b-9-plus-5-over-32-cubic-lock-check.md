# 237 — Locked Candidate Check: `b = 9 + 5/32` and cubic closure

## Scope
Lock the two-loop coefficient candidate
- `b = 9 + 5/32 = 9.15625`

Then:
1. Solve exact cubic coefficient `c` from CODATA closure of
   `alpha^-1 = 137 + 5alpha - b alpha^2 + c alpha^3`.
2. Test whether linked mass-ratio series with shared scaling
   `B = 4b`, `C = 4c`
   closes `mp/me` simultaneously.

## Exact values (non-rounded)
Using:
- `alpha^-1_phys = 137.035999177`
- `alpha = 1/alpha^-1_phys`
- `mp/me_phys = 1836.15267343`

Solved cubic coefficient:
- `c_exact = -0.007996757405957876626533392673425759...`

### Alpha closure
By construction with this `c`:
- `alpha^-1_pred = alpha^-1_phys` (exact at arithmetic precision used).

### mp/me simultaneous check with `B=4b`, `C=4c`
Predicted:
- `mp/me_pred = 1836.152645131223748200994173331643...`

Residual:
- `Delta(mp/me) = -2.8298776251799005826668356736367e-5`
- `|Delta| = 0.0154119952339997208 ppm`
- `|Delta| = 15.4119952339997208 ppb`

## Joint-closure scale factor implied by data
If same `b,c` are retained but scale factor is free in
`B = g b`, `C = g c`, exact joint closure requires:
- `g_exact = 3.941961445565373788901249789775448...`
- `g_exact - 4 = -0.0580385544346262111...`

## Conclusion
- Locked candidate `b = 9 + 5/32` is extremely strong for the alpha lane.
- Applying strict grade-1 multiplier `4` at cubic order does **not** exactly close `mp/me`;
  it leaves a small but nonzero residual (`15.412 ppb`).
- The data-preferred shared scaling is `g ≈ 3.9419614456`, very near 4 but not equal.

This isolates the remaining problem cleanly: explain the `g` offset structurally, or add next-order terms with grade-constrained coefficients.
