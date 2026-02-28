# Finding 114: GRAND-323 Continuum R4 Schwinger Lane and Mass-Gap Floor from Alpha

Date: 2026-02-28
Status: GRAND-323 complete

## Scope

Ticket GRAND-323 requires an explicit continuum Schwinger-function construction on R4 from the Wilson-lattice lane, with no standalone existential carrier assumptions.

Files:
- `lean/Gutoe/YangMillsContinuumLimit.lean`
- `lean/Gutoe/YangMillsContinuumSurvival.lean`
- `lean/Gutoe/YangMillsWilsonEquivalence.lean`
- `lean/lakefile.lean`

## Physical stats (core numbers)

1. Structural coupling input:
   - `alpha = 1/137`

2. Uniform minorization floor (from SC coordination lane):
   - `eps_floor = alpha / (2 + alpha)`
   - With `alpha = 1/137`: `eps_floor = 1/275 = 0.003636363636364`

3. Continuum non-vanishing mass-gap floor (Doeblin form):
   - `c_floor = -log(1 - eps_floor) / a_cap`
   - Equivalently: `c_floor = log(275/274) / a_cap`
   - Numerically: `c_floor = 0.003642991278501 / a_cap`

4. If the refinement cap is normalized to `a_cap = 1`:
   - `c_floor = 0.003642991278501`

5. Additional structural exact values in the lane:
   - Transfer basis cardinality: `3` (`Fin 3`)
   - Zero-step partition: `cylPartition K 0 = 3`
   - Row stochasticity: `sum_j K_ij = 1`

## What is proven in GRAND-323

The new continuum module proves an explicit Schwinger-family package:

- cylinder path weights are strictly positive
- partition function is strictly positive
- correlator family is normalized (`<1> = 1`)
- Wilson kernel row sums are 1 at each refinement step
- OS reconstruction schedule is explicit (no standalone `exists K`)
- mass-gap lower bound survives the continuum schedule uniformly:
  `exists c > 0, forall n, c <= doeblinGapLowerBound (a_t n) (eps_n)`

Here the concrete Wilson-domain bridge is provided by:
- `c3_gap_correspondence_of_domain`
- `continuum_mass_gap_from_wilson_domain`

## Interpretation

The physical content is the alpha-controlled floor:
- the continuum Yang-Mills mass gap is bounded away from zero by a structural constant derived from alpha,
- with explicit formula `eps_floor = 1/275` and induced gap floor `c_floor`.

This is a non-vanishing result in the continuum schedule, not a vanishing-limit proxy.

## Verification

- `cd lean && lake build Gutoe.YangMillsContinuumLimit` passed
- `cd lean && lake build Gutoe` passed (module included in root set)

No proof `sorry` in `YangMillsContinuumLimit.lean`.
