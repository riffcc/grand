# Finding 072: GRAND-301 Full State-Space Lift Progress

Date: 2026-02-26
Status: GRAND-301 in progress (hardening pass complete)

## New module

- `lean/Gutoe/YangMillsFullStateLift.lean`

## What landed

### 1) Formal reduced->full mapping layer

- `liftObservable`
- `liftDiagonalOperator`
- `liftObservable_fiber_constant`

This gives an explicit Lean map from reduced observables to full-state operators.

### 2) Gap-preservation proof obligations

- `FullStateLiftGapObligations`
- `full_gap_positive_of_lift_obligations`

This isolates exactly what must be true to transfer strict gap positivity from
reduced-lane control into full-state-lane control.

### 3) Counterexample checks (failure modes)

- `NoExtraUnitMode`
- `identityKernel4_not_noExtraUnitMode`
- `mass_gap_zero_at_unit_ratio`
- `reduced_positive_full_zero_counterexample`

This documents a concrete failure pattern: reduced gap can be positive while full
gap is zero if extra unit modes are not excluded.

### 4) Obligation necessity is now explicit (not implicit)

- `FullStateLiftGapObligationsNoRatioDom`
- `FullStateLiftGapObligationsNoReducedBound`
- `no_ratio_dom_rule_false`
- `no_reduced_ratio_bound_rule_false`
- `lift_obligation_families_are_independently_necessary`

These prove each reduced↔full coupling obligation is load-bearing. Dropping
either one admits a concrete zero-gap counterexample.

### 5) Wilson/center construction now feeds the obligation lane directly

- `WilsonCenterLiftSpectralHypotheses`
- `reduced_ratio_bound_of_wilson_closed_form`
- `lift_schedule_obligations_of_wilson_center_schedule`
- `full_gap_positive_all_steps_of_wilson_center_schedule`

This closes the "plumbing" gap: the epsilon schedule is no longer hand-wired.
It is derived from the actual Wilson/center row-total construction via
`minorizationEps (wilsonRowTotalsSchedule ...)`.

### 6) Seam attack: reduced ratio bound now derivable from mode data

- `wilsonReducedKernel`
- `WilsonReducedModeHypotheses`
- `reduced_ratio_bound_of_wilson_reduced_modes`
- `WilsonCenterFullGapModeHypotheses`
- `full_gap_positive_all_steps_of_wilson_center_modes`

This is the high-impact reduction:
we can now obtain the reduced-lane bound from concrete zero-sum eigenmodes +
Doeblin decomposition on the Wilson-induced kernel, rather than assuming the
reduced ratio inequality directly.

### 7) Seam attack follow-through: ratio-dominance now derivable from mode-dominance

- `WilsonModeDominanceHypotheses`
- `ratio_dom_of_wilson_mode_dominance`
- `full_gap_positive_all_steps_of_wilson_center_modes_dominance`

This removes direct `ratio_dom` from the end-to-end theorem path. The remaining
load-bearing assumption is now explicit and narrower: reduced-mode dominance
over zero-sum modes on the Wilson-induced reduced kernel.

### 8) Max seam closure: full-lane cap now derivable from full-mode data only

- `full_gap_positive_all_steps_of_wilson_center_from_full_modes_only`

This theorem removes both:
- any separate reduced-lane spectral package, and
- any explicit reduced→full ratio-dominance hypothesis.

For each refinement step it derives the full-lane ratio cap directly from:
- a concrete eigenmode on `wilsonReducedKernel`,
- structural Doeblin control from Wilson row totals,
- principal normalization `lambda0F = 1`.

The remaining non-structural seam is now isolated to one hypothesis family:
realization of the targeted full-mode on the Wilson-induced reduced kernel.

### 9) Identified-mode specialization landed

- `full_gap_positive_all_steps_of_wilson_center_identified_mode`

If we identify the full-lane subdominant mode directly with the concrete
reduced-kernel mode (`lambda0 = 1`), strict positivity now follows from one
mode package (`WilsonReducedModeHypotheses`) plus structural Doeblin control.
This removes an extra bookkeeping layer in the closure lane.

### 10) Structural no-extra-unit-mode theorem (production kernel family)

- `wilson_no_extra_unit_mode_zero_sum`

This discharges the old "extra unit mode" concern on the Wilson-induced reduced
kernel itself: for `alpha>0`, no nontrivial zero-sum eigenmode can have
eigenvalue `1`. The result is purely structural (Doeblin + SC row totals), not
a numerical certificate.

### 11) Unconditional full-state Doeblin gap certificate

- `full_doeblin_gap_positive_all_steps_of_wilson_center`

This gives a fully structural, zero-assumption positivity certificate per step
for the Wilson/center lane:
`doeblinGapLowerBound(a_t n, eps_n) > 0` from only `a_t>0` and `alpha>0`.

Interpretation: even when an explicit real full-mode witness is not yet
constructed, strict positive full-state gap lower bounds are already proven at
every refinement step.

## Cl(1,3) anchor

- `reduced_transfer_basis_dim_eq_three` reuses transfer-basis dimension from
  the Z3/Cl(1,3) lane.

## Build sanity

- `lake build Gutoe.YangMillsFullStateLift` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.

## Remaining for GRAND-301 closure

- Discharge full-mode realization constructively on the production kernel family
  (replace `hFullMode` with a structural spectral realization theorem).
- Lift the reduced-kernel no-extra-unit result into the final full-state
  realization statement without reintroducing spectral assumptions.
