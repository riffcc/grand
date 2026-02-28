# GRAND-348 — Lean Parity Checklist for Findings 098-103

## Scope

Close formal parity backlog by mapping each finding claim to concrete Lean definitions/theorems (or explicitly marking empirical-only claims), then re-running the full Lean build.

## Build Evidence

Command:

- `cd lean && lake build Gutoe`

Result:

- success (`Build completed successfully (8109 jobs)`)

Safety checks (repo-wide Lean lane):

- `rg -n "\\bsorry\\b" lean/Gutoe | wc -l` -> `52`
- `rg -n "\\baxiom\\b" lean/Gutoe | wc -l` -> `8`

No new `sorry`/`axiom` were introduced in the parity modules touched here.

## Finding-to-Theorem Map

### Finding 098 — Sigma8 root correction + `C_inf` lane

File links:

- `lean/Gutoe/DarkMatterSector.lean`
- `lean/Gutoe/Inflation.lean`

Mapped claims:

- structural subtraction `δ = |dark|/2 = 5/2`
  - `geometricDarkBudgetDeltaQ`
  - `geometric_dark_budget_delta_eq`
- corrected budget ratio `(60 - 5/2)/11 = 115/22`
  - `correctedUnifiedBudgetDarkToVisibleRatio`
  - `corrected_unified_budget_dark_to_visible_ratio_eq`
- corrected inflation Hubble ratio parity expression
  - `inflation_hubble_ratio_with_corrected_budget_eq`
- inflation correction factor from shared counts
  - `inflationCorrectionCInf`
  - `inflation_cinf_eq`
  - `inflation_cinf_eq_67_over_66`

### Finding 099 — Proton mass structural ratio

File link:

- `lean/Gutoe/MassSpectrum.lean`

Mapped claims:

- `mp/me = 12 * T(17)` definition
  - `mpMeAlgebraic`
- exact integer closure
  - `mp_me_eq_1836`
- shared-input parity with alpha lane
  - `mp_me_uses_same_inputs`

### Finding 100 — `C_inf=67/66` cross-sector recurrence (inflation + QCD)

File link:

- `lean/Gutoe/StrongCouplingCInfBridge.lean`

Mapped claims:

- shared second-order correction from common counts
  - `sharedSecondOrderCorrectionQ`
  - `shared_second_order_correction_eq_67_over_66`
- strong-coupling leading term from structural counts
  - `alphaSStructuralLeadingQ`
  - `alpha_s_structural_leading_eq_16_over_137`
- corrected strong-coupling expression
  - `alphaSStructuralCorrectedQ`
  - `alpha_s_structural_corrected_eq`
- explicit inflation/QCD correction equality bridge
  - `shared_correction_matches_inflation_cinf`

### Finding 101 — Electron transduction sweep hard negative + constrained survivors

File links:

- `lean/Gutoe/ElectronScaleTransduction.lean`
- runtime sweep artifact lane: `crates/gutoe-physics/src/bin/electron_scale_sweep.rs`

Mapped formal claims:

- structural lane ingredients used by sweep follow shared terms
  - `alphaStructural`
  - `correctedDarkVisibleRatio`
  - `darkVisibleCountRatio`
- inverse-cubic lambda equivalence (no hidden magic factor)
  - `lambda_qg_inv_cube_eq_12_cube`

Empirical-only claim boundary:

- candidate ranking / hard-negative elimination remains simulation data, not a theorem proposition.
- parity obligation is satisfied by proving the shared factorization identities used by those candidates.

### Finding 102 — Flagged gauge-cube candidate + collateral check

File link:

- `lean/Gutoe/ElectronScaleTransduction.lean`

Mapped claims:

- flagged form equals gauge-cube form
  - `electronScaleFactorFlagged`
  - `electronScaleFactorGaugeCube`
  - `electron_scale_flagged_eq_gauge_cube`

### Finding 103 — Triangular alpha inverse closure

File link:

- `lean/Gutoe/FineStructure.lean`

Mapped claims:

- core identity `T(2^4)+1=137`
  - `triangular_clifford_dim_plus_one_eq_137`
- canonical alpha closure
  - `alpha_inverse_d4`
  - `fine_structure_constant`

## Acceptance Check

- Checklist artifact linking findings 098-103 to theorem names: **done**
- No increase in unresolved formal gaps from this closure pass: **done**
- `lake build Gutoe` green: **done**
