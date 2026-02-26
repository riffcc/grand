# Finding 055 — Signature `√2` + Micro-Mode Closure Candidate for Λ (Conjectural, GRAND-293/295)

## Status
Conjectural bridge (not final proof).

`GRAND-92` established the structural baseline:
- `Λ_struct = (13/100)^137 / l_P^2`
- `Λ_struct / Λ_obs = 1.411306927...`

This residual is numerically close to `√2 = 1.414213562...`.

## Numerical check
From `/tmp/bh_renders/lambda_cosmological_report.json`:
- `ratio_struct_over_obs = 1.411306927110`
- `sqrt2 = 1.414213562373`
- `k_signature = sqrt(6/3) = 1.414213562373`
- `k_required = Λ_struct / Λ_obs = 1.411306927110`
- `k_required / k_signature = 0.997944698495`
- relative gap vs `√2`: `0.20553%`
- `residual_over_sqrt2 = 0.997944698495`

Applying a `1/√2` correction gives:
- `Λ_sig = Λ_struct / √2 = 1.103327932868e-52`
- `Λ_sig / Λ_obs = 0.997944698495`
- relative error vs observation: `0.20553%`

Applying the GRAND-295 micro finite-mode factor
`k_micro = N_micro/(N_micro-1) = 486/485` with
`N_micro = ewsbScaleFactor + |grade-2| = 480 + 6 = 486` gives:
- `Λ_full = Λ_sig * (486/485) = 1.105602561045e-52`
- `Λ_full / Λ_obs = 1.000002316430`
- relative error vs observation: `2.316e-6` (~0.0002316%)

## Candidate chain
Conjectural corrected formula:

`Λ = (13/100)^137 / (√2 · l_P^2) · (486/485)`

where:
- `13/100` from `Gutoe.EWSBHiggs.higgs_quartic_eq_13_100`
- `137` from `Gutoe.FineStructure.alpha_inverse_d4`
- `√2` from the exact structural split factor
  `sqrt(|grade2| / |timelike-spacelike bivectors|) = sqrt(6/3) = sqrt(2)` (proved in Lean)
- `486/485` from finite-mode correction
  `N_micro/(N_micro - N_fixed)` with
  `N_micro = 480 + 6` and `N_fixed = 1` (unique Z3-fixed grade-1 generator)
- `l_P` shared Planck scale anchor

## Lean/Rust support added
- Lean:
  - `lean/Gutoe/CosmologicalConstant.lean`
  - `lorentzSignatureNormalization`
  - `lorentz_signature_normalization_eq_sqrt2`
  - `bivector_signature_split_3_3`
  - `lorentzSignatureNormalizationFromParity`
  - `lorentz_signature_normalization_from_parity_eq_sqrt2`
  - `lorentz_signature_normalization_eq_from_parity`
  - `lambdaCosmologicalSignatureCandidate`
  - `lambdaCosmologicalSignatureFromSplit`
  - `z3FixedGrade1Count`
  - `microModeCount`
  - `microFiniteModeRescale`
  - `lambdaCosmologicalFullCandidate`
  - `lambda_signature_from_split_eq_candidate`
  - `lambda_cosmological_signature_candidate_eq`
  - `lambda_cosmological_signature_candidate_pos`
  - `lambda_cosmological_full_candidate_eq`
  - `lambda_cosmological_full_candidate_pos`
- Rust:
  - `crates/gutoe-physics/src/constants.rs`
  - `lambda_cosmological_suppression()`
  - `lambda_cosmological_structural()`
  - `lambda_cosmological_signature_candidate()`
  - `lambda_micro_mode_count()`
  - `lambda_micro_finite_mode_rescale()`
  - `lambda_cosmological_full_candidate()`
  - test `test_lambda_cosmological_signature_candidate_is_close_to_observed`
  - `crates/gutoe-physics/src/bin/lambda_cosmological_report.rs` now emits both structural and signature-candidate branches

## What remains for closure
- `GRAND-293` (done): `√2` factor derived from Cl(1,3) signature split.
- `GRAND-295`: Lean/Rust parity now includes the finite-mode rescale path (`486/485`).
  Remaining closure criterion is interpretive, not numerical:
  derive/justify this finite-mode normalization as the unique physically admissible
  continuum-compatible correction (no hidden free parameter).
