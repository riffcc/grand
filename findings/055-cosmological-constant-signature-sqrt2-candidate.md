# Finding 055 — Signature `√2` Candidate for Λ Residual (Conjectural, GRAND-293)

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

## Candidate chain
Conjectural corrected formula:

`Λ = (13/100)^137 / (√2 · l_P^2)`

where:
- `13/100` from `Gutoe.EWSBHiggs.higgs_quartic_eq_13_100`
- `137` from `Gutoe.FineStructure.alpha_inverse_d4`
- `√2` from the exact structural split factor
  `sqrt(|grade2| / |timelike-spacelike bivectors|) = sqrt(6/3) = sqrt(2)` (proved in Lean)
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
  - `lambda_signature_from_split_eq_candidate`
  - `lambda_cosmological_signature_candidate_eq`
  - `lambda_cosmological_signature_candidate_pos`
- Rust:
  - `crates/gutoe-physics/src/constants.rs`
  - `lambda_cosmological_suppression()`
  - `lambda_cosmological_structural()`
  - `lambda_cosmological_signature_candidate()`
  - test `test_lambda_cosmological_signature_candidate_is_close_to_observed`
  - `crates/gutoe-physics/src/bin/lambda_cosmological_report.rs` now emits both structural and signature-candidate branches

## What remains for closure
- `GRAND-293` (now closed): `√2` factor derived from Cl(1,3) signature split.
- Remaining micro-residual tracked in `GRAND-295`:
  - `k_micro = k_required / k_signature = 0.997944698495` (~0.2055%).
  - target is a no-knob structural correction for this final factor.
