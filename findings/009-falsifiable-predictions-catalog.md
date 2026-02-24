# 009 — Falsifiable Predictions Catalog (GUTOE)

Status: draft-complete catalog for issue GRAND-124.

## Scope

This is the explicit kill-list for GUTOE: observations that would falsify core claims.

## Predictions and Falsifiers

1. Electroweak angle from Z3 orbit structure
- Claim: `sin^2(theta_W) = 3/13` (tree-level structural value).
- Basis: `lean/Gutoe/Z3Uniqueness.lean`, `lean/Gutoe/GaugeConstants.lean`, `crates/gutoe-em/src/weak.rs`.
- Falsifier: if precision-corrected structure cannot be embedded as perturbative corrections around the fixed structural value, model fails.

2. Gauge generator count from grade decomposition
- Claim: `dim(U(1))+dim(SU(2))+dim(SU(3)) = 12`.
- Basis: `lean/Gutoe/GaugeGroupSM.lean`, `lean/Gutoe/GaugeConstants.lean`.
- Falsifier: if Clifford-grade map cannot produce observed gauge algebra without ad-hoc additions.

3. Fine-structure leading order from triangular structure
- Claim: `alpha_LO^-1 = 137` with runtime offset interpreted as higher-order/QED.
- Basis: `lean/Gutoe/FineStructure.lean`, `crates/gutoe-physics/src/constants.rs`.
- Falsifier: if correction hierarchy cannot reproduce measured `alpha^-1 ~= 137.036` without free fit knobs.

4. Lattice dispersion coefficient
- Claim: `lambda_qg = 1/12` at leading order.
- Basis: `lean/Gutoe/LambdaQG.lean`, `lean/Gutoe/DispersionRelation.lean`, `crates/gutoe-physics/src/constants.rs`.
- Falsifier: if required sign/magnitude for stable causal dispersion is inconsistent with proofs + runtime parity.

5. Black-hole observables (GUTOE vs GR)
- Claim: measurable deviations in high-resolution horizon observables (ring substructure, transfer behavior) must be internally consistent and bounded.
- Basis: `lean/Gutoe/GravityMetric.lean`, `lean/Gutoe/SynchrotronTransfer.lean`, renderer/runtime parity checks.
- Falsifier: if parity harnesses fail or inferred observables demand contradictions with proved constraints.

6. Upcoming explicit high-risk tests
- Photon-ring coherent interference and lattice diffraction signatures.
- Basis: GRAND-217, GRAND-218.
- Falsifier: null/incompatible signal where model predicts resolvable structure at stated scales.

## Acceptance Artifacts

- This file establishes explicit falsification gates.
- Follow-on work should bind each row to a reproducible dataset/test command and confidence envelope.
