# 207 — Even-subalgebra suppression (`1/2`) and critical offset

## New theorem module
- `lean/Gutoe/EvenSubalgebraSuppression.lean`

Key closures:
- even-grade count = `8`
- odd-grade count = `8`
- suppression = `8/16 = 1/2`
- canonical gain lane with `(branching, merge) = (2, 1)` gives:
  - `G_eff = 2 * 1 * (1/2) = 1` (exact)
- for uncapped measured gain `1.9992 = 2499/1250`:
  - required unit-gain suppression is exactly `1250/2499`
  - offset from half is exactly `1/4998`

Build status:
- `lake build Gutoe.EvenSubalgebraSuppression` passed
- `lake build Gutoe` passed (warnings only)

## New probe
- `ctc_even_subalgebra_suppression_probe`

Outputs:
- `/tmp/bh_renders/ctc_even_subalgebra_suppression_probe/ctc_even_subalgebra_suppression_probe.txt`
- `/tmp/bh_renders/ctc_even_subalgebra_suppression_probe/ctc_even_subalgebra_suppression_probe.json`

Measured outputs (default):
- `suppression_even = 0.5`
- canonical (`branching=2`, `merge=1`): `G_eff = 1.0`
- with `G_uncapped = 1.9992`:
  - `suppression_for_unit_gain = 0.500200080032...`
  - offset above half: `+2.000800320128e-4`

Cross-check (Z3+void lane in same probe):
- set `branching=3`, `merge=3/16`:
  - `G_eff = 3 * (3/16) * (1/2) = 9/32 = 0.28125`
  - so even-filter `1/2` alone does **not** close that lane.

## Interpretation
- The `1/2` suppression claim is structurally exact from Cl(1,3) grading.
- It closes exactly if the active gain lane is `2 × 1 × suppression`.
- The observed `0.500200...` is consistent with a tiny correction above `1/2`
  when matching against the measured `1.9992` uncapped gain.
