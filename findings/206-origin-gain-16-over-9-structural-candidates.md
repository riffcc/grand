# 206 — Origin gain `16/9` has exact structural decompositions

## New probe
- `ctc_gain_ratio_structural_probe`

Outputs:
- `/tmp/bh_renders/ctc_gain_ratio_structural_probe/ctc_gain_ratio_structural_probe.txt`
- `/tmp/bh_renders/ctc_gain_ratio_structural_probe/ctc_gain_ratio_structural_probe.json`

## Constraint tested
For the closure equation with structural branching and void merge:
- `1 = branching * merge_fraction * (eta * infra_gain)`
- `branching = 3` (Z3 order)
- `merge_fraction = 3/16` (void lane)

Required product:
- `eta * infra_gain = 1 / ((3)*(3/16)) = 16/9`

## Exact structural matches for `16/9`
The probe found exact (delta `= 0`) expressions from existing rails:

1. `16/9 = basis_16 / z3^2`
   - interpretation: full Clifford state count over generation-orbit square.
2. `16/9 = grade1^2 / z3^2`
   - interpretation: spacetime-vector square over generation-orbit square.
3. `16/9 = basis_16 / (ew_sum - z3_fixed_grade1)`
   - since `ew_sum = grade1 + grade2 = 10` and `z3_fixed_grade1 = 1`, denominator is `9`.

## Near misses
- `sin2w / lambda_h = 300/169 = 1.7751479` (very close; `|delta| ~ 2.63e-3`)
- `1 + cos2w = 23/13 = 1.7692308`
- `24/13 = 1.8461538`

## Interpretation
- `16/9` is not an arbitrary tuned value in this lane; it has multiple exact
  decompositions using established Cl(1,3)/Z3 counts.
- This does **not** yet uniquely determine a split between `eta` and `infra_gain`;
  it only fixes their product.

## One exact split candidate (`eta <= 1`)
A low-complexity structural sweep found:

- `eta = 2/3` (example count form: `grade1/grade2 = 4/6`)
- `infra_gain = 8/3` (example count form: `basis_16/grade2 = 16/6`)
- Product: `(2/3) * (8/3) = 16/9` (exact)

Then with `branching=3`, `merge=3/16`:
- `G_eff = 3 * (3/16) * (2/3) * (8/3) = 1` (exact knife-edge)
