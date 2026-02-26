# 009 — Falsifiable Predictions Catalog (GUTOE)

Status: quantitative kill-criteria catalog for issue GRAND-124.

## Scope

This is the explicit kill-list for GUTOE: observations that would falsify core claims.

## Predictions and Numeric Kill Criteria

All criteria below are intentionally binary: pass window or dead.

1. `sin^2(theta_W)` structural anchor and EW corrected value
- Structural claim: `sin^2(theta_W)_LO = 3/13 = 0.2307692308` (exact).
- Measured target (EW scale): `sin^2(theta_W)_exp ~= 0.23122`.
- Required correction magnitude at EW scale: `delta_EW ~= 4.51e-4`.
- Kill criterion: if the no-free-parameter correction pipeline predicts
  `sin^2(theta_W)` outside `[0.23100, 0.23140]`, GUTOE fails this gate.

2. `m_Z/m_W` ratio from Weinberg structure
- Structural claim: `(m_Z/m_W)^2 = 13/10`, so `m_Z/m_W_LO = sqrt(13/10)`.
- Measured target: `m_Z/m_W_exp ~= 1.1345`.
- Kill criterion: corrected prediction must lie in `[1.1335, 1.1355]`.
  Outside this band, fail.

3. Fine-structure constant bridge (Lean LO -> runtime measured)
- Structural claim: `alpha_LO^-1 = 137` (exact Lean theorem).
- Measured target: `alpha_exp^-1 ~= 137.036`.
- Kill criterion: correction chain must produce `alpha^-1` in
  `[137.034, 137.038]` with no fitted free parameter.
  Any value outside this band fails.

4. Lattice dispersion coefficient
- Claim: `lambda_qg = 1/12` (exact leading-order structural coefficient).
- Kill criterion: any sign flip (`lambda_qg <= 0`) or coefficient outside
  `1/12 ± 1e-12` in theorem/runtime parity checks fails.

5. Gauge-sector counting
- Claim: total gauge generators = `12`.
- Kill criterion: if low-energy phenomenology requires extra unscreened gauge
  generators not derivable from the Cl(1,3) decomposition, fail.

6. Strong-CP / neutron-EDM gate
- Runtime bridge claim: `theta_qcd` is structurally pinned to zero in the
  current Cl(1,3) dynamics map (`theta_qcd = 0`).
- Chiral estimate bridge: `|d_n| ~= 2.4e-16 * |theta_qcd| e*cm`.
- Kill criteria (numeric):
  - `|theta_qcd| <= 4.2e-11`
  - `|d_n| <= 1.0e-26 e*cm`
  Any violation fails this gate.

7. Null-result consistency gates
- Proton decay: predicted lifetime must satisfy `tau_p >= 1e34 years`.
- Fifth force/LV: predicted signal strengths must stay below configured bounds.
- Kill criterion: any bound violation fails.

8. Yang-Mills mass-gap structural gate (SU(3) sector)
- Structural theorem path: Doeblin minorization + stochastic spectral contraction gives
  `m_gap >= -log(1-eps)/a_t > 0` in the transfer-basis construction.
- Kill criterion (structural): if Lean cannot discharge positivity of the gap from
  row-stochastic/Doeblin hypotheses (without an external `|mu| ≤ 1` axiom), this gate fails.
- Status:
  - Theorem A (structural gap existence): closed.
  - Theorem B (continuum-survival bridge): open.
  - Theorem C (Wilson-action equivalence bridge): open.

## Basis Theorems / Code Anchors

- `lean/Gutoe/Z3Uniqueness.lean`
- `lean/Gutoe/GaugeConstants.lean`
- `lean/Gutoe/FineStructure.lean`
- `lean/Gutoe/LambdaQG.lean`
- `lean/Gutoe/DispersionRelation.lean`
- `lean/Gutoe/FalsifiabilityCatalog.lean`
- `lean/Gutoe/YangMillsStructuralGap.lean`
- `lean/Gutoe/YangMillsMassGap.lean`
- `crates/gutoe-physics/src/constants.rs`

## Acceptance Artifacts

- This file now defines numeric fail windows.
- Runtime gates are now wired in Rust:
  - `crates/gutoe-physics/src/falsifiability.rs`
  - `evaluate_structural()` enforces structural/gauge/`lambda_qg` checks.
  - `evaluate_with_corrected(...)` enforces numeric EW-scale windows.
- Current status is explicit and honest:
  - Structural gates pass in tests.
  - Corrected EW gates are still marked as a bridge gap until the full no-free-parameter correction chain is fully wired.
