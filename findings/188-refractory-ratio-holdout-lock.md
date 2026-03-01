# 188 — Refractory Ratio Holdout Lock (4d↔5d)

## Scope
- Run the actual bidirectional transition-row holdout under ratio lock:
  - 4d-fit -> 5d-validate (period-5 transition train, period-6 transition holdout)
  - 5d-fit -> 4d-validate (period-6 transition train, period-5 transition holdout)
- Candidate space: rational pairs with denominator <= 20 and exact lock `g_f / g_v = 12/7`.
- Verify exact structural ratio arithmetic for current structural pair:
  - `(3/5) / (7/20) = 12/7`.

## Runner
- Binary: `refractory_ratio_holdout`
- Output: `/tmp/nuclear_chart/refractory_ratio_holdout_latest/refractory_ratio_holdout.txt`
- Candidate count under lock: 15.

## Bidirectional Holdout Results
Best (pure thermal MAE-K objective under lock):
- `g_f = 12/13`, `g_v = 7/13`
- 4d-fit -> 5d-validate:
  - train p5 thermal MAE = `582.527047400 K`
  - validate p6 thermal MAE = `884.229500611 K`
- 5d-fit -> 4d-validate:
  - train p6 thermal MAE = `884.229500611 K`
  - validate p5 thermal MAE = `582.527047400 K`
- Consensus winner across both directions and minimax: `true`.

Structural pair check:
- `g_f = 3/5`, `g_v = 7/20`
- exact ratio lock: `(3/5)/(7/20) = 12/7` (true)
- holdout thermal MAE:
  - p5 = `585.776618350 K`
  - p6 = `884.230203944 K`
- delta vs pure-MAE winner:
  - p5: `+3.249570950 K` (`+0.557840%`)
  - p6: `+0.000703333 K` (`+0.0000795%`)

Interpretation:
- Under this objective, there is a near-flat basin along the ratio-locked family.
- `3/5, 7/20` remains exact-ratio-consistent and very close to the best MAE point, especially on p6.

## Lean Lock
Added ratio theorem (ratio-level claim, not coefficient-level lock):
- File: `lean/Gutoe/ThermalEntropyClosure.lean`
- New definitions/theorems:
  - `gaugeGeneratorCount`
  - `spatialVectorCount`
  - `spatialBivectorCount`
  - `oddParityBasisCount`
  - `refractorySuppressionRatioQ`
  - `refractory_suppression_ratio_eq_12_over_7`

Key closure:
- `refractorySuppressionRatioQ = 12/7`
- with `12` from gauge generator count and `7` from odd-parity basis proxy (`3 + 3 + 1`).

## Verification
- `lake build Gutoe.ThermalEntropyClosure` passes.
- `lake build Gutoe` passes.
- Holdout report generated at the output path above.
