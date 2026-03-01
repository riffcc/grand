# Finding 229 — Weak-Angle Identity Fit (Multi-Point)

## Summary

Ran a coarse multi-point fit to test whether the leading-order identity
`sin²(theta_W) = 3/13` survives public anchors, and whether a corrected form is
recoverable.

Binary:

- `crates/gutoe-physics/src/bin/ctc_weak_angle_identity_fit.rs`

Outputs:

- `/tmp/bh_renders/ctc_weak_angle_identity_fit/ctc_weak_angle_identity_fit.txt`
- `/tmp/bh_renders/ctc_weak_angle_identity_fit/ctc_weak_angle_identity_fit.json`

## Data Points Used

1. `PDG_MSbar_at_MZ`: `0.23122 ± 0.00006`
2. `LHCb_2024_sin2thetaeff_l`: `0.23147 ± 0.000499` (combined)
3. `SLAC_E158_lowQ`: `0.2397 ± 0.001281`

## Models Tested

1. `base_fixed`: `3/13`
2. `rational_fixed`: `3/13 + 1/13^3`
3. `base_plus_delta0`: `3/13 + δ0`
4. `base_plus_delta0_delta1log`: `3/13 + δ0 + δ1 ln(M_Z/Q)`
5. `rational_plus_delta1log`: `(3/13 + 1/13^3) + δ1 ln(M_Z/Q)`

## Fit Results

- `base_fixed`: `red_chi2 = 35.68` (hard fail)
- `rational_fixed`: `red_chi2 = 14.68` (still hard fail)
- `base_plus_delta0`: `red_chi2 = 21.98` (M_Z corrected, low-Q still fails)
- `base_plus_delta0_delta1log`: `red_chi2 = 0.247`
- `rational_plus_delta1log`: `red_chi2 = 0.124` (best in this set)

Best-fit parameters:

- `δ0 ≈ 4.54e-4`
- `δ1 ≈ 1.34e-3` (log-running coefficient)

## Key Pulls

For `base_fixed = 3/13`:

- PDG `M_Z`: `-7.51σ`
- LHCb `M_Z`: `-1.40σ`
- E158 low-Q: `-6.97σ`

For `rational_plus_delta1log` (best coarse fit):

- PDG `M_Z`: `+0.073σ`
- LHCb `M_Z`: `-0.492σ`
- E158 low-Q: `~0σ`

## Interpretation

- Leading-order mapping `3/13 -> sin²(theta_W)` is falsified by these anchors.
- A corrected identity with a small additive shift and log-running term
  reconciles this coarse set.
- This does **not** close a full electroweak global fit; it is a focused
  three-anchor sanity fit.
