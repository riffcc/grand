# 030 — Strong CP Path-Integral Guardrails (GRAND-267)

Status: in progress; structural skeleton formalized, full nonperturbative closure still open.

## What was added now

- New Lean module: `lean/Gutoe/StrongCPPathIntegral.lean`
- Added finite-sector path-integral skeleton tied to Cl(1,3) sector counts:
  - `N₊ = |magneticTriplet| = 3`
  - `N₋ = |emTriplet| = 3`
  - `Z_im(θ) = (N₊ - N₋) sin θ`
  - therefore `Z_im(θ)=0` for all `θ` in this skeleton.
- Strengthened route-2 reflection channel:
  - `zRe_even_theta_reflection`: real channel is even under `θ ↦ -θ`
  - `zIm_odd_theta_reflection`: imaginary channel is odd under `θ ↦ -θ`
  - `zComplex_even_theta_reflection`: complex partition coefficient is reflection-even in the balanced Cl(1,3) skeleton
- Added route-2 exclusion lemma in the principal branch:
  - if candidates are constrained to `{0, π}` and vacuum weight is positive,
    then `θ = 0` (since `Z_re(π) = -6`).

## Why this matters

This goes beyond the previous classical proxy-only statement by encoding
vacuum-sector phase-channel behavior in Lean and making the `π` branch
exclusion explicit under a positivity condition.

## What is still missing (the real GRAND-267 core)

To claim a full nonperturbative Strong-CP dissolution, one of these must be
proved in the Cl(1,3) framework (not just assumed):

1. No nontrivial topological sectors are physically accessible (`Q ≠ 0` absent).
2. Or `Z(θ)=Z(-θ)` plus a model-internal argument that excludes `π` without
   importing ad hoc assumptions.
3. Or nonzero-`Q` sectors are suppressed in a continuum-stable way.

The finite-sector skeleton is a bridge artifact, not the final theorem.
