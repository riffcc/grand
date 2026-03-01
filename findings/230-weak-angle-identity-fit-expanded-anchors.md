# Finding 230 — Weak-Angle Identity Fit (Expanded Anchors)

## Summary

Expanded the weak-angle identity test from 3 to 6 anchors and re-ran the same
model family, including the zero-free-parameter candidate:

`sin²(theta_W)(Q) = 508/2197 + [alpha*ln(10)/(4*pi)]*ln(M_Z/Q)`

Binary:

- `crates/gutoe-physics/src/bin/ctc_weak_angle_identity_fit.rs`

Outputs:

- `/tmp/bh_renders/ctc_weak_angle_identity_fit/ctc_weak_angle_identity_fit.txt`
- `/tmp/bh_renders/ctc_weak_angle_identity_fit/ctc_weak_angle_identity_fit.json`

## Anchor Set (6 points)

At `Q = M_Z`:

1. PDG MS-bar value
2. LHCb 2024 measurement
3. CMS 2025 measurement
4. Tevatron Run-II combination
5. ATLAS 2015 measurement

Low-Q:

6. SLAC E158

## Results

### Hard fail for naive identity

- `3/13`: `reduced chi2 = 19.59`

### Hard fail for static rational-only correction

- `3/13 + 1/13^3`: `reduced chi2 = 7.61`

### Surviving corrected-running forms

- `3/13 + delta0 + delta1*ln(M_Z/Q)`: `reduced chi2 = 0.451`
- `(3/13 + 1/13^3) + delta1*ln(M_Z/Q)`: `reduced chi2 = 0.377`
- `(3/13 + 1/13^3) + [alpha ln(10)/(4pi)] ln(M_Z/Q)`:
  - `chi2 = 1.88367`
  - `dof = 6`
  - `reduced chi2 = 0.31395`

## Pulls for Zero-Free-Parameter Candidate

1. PDG `M_Z`: `+0.073σ`
2. LHCb 2024: `-0.492σ`
3. CMS 2025: `-0.954σ`
4. Tevatron combo: `-0.775σ`
5. ATLAS 2015: `+0.356σ`
6. E158 low-Q: `-0.001σ`

All six anchors remain within `|pull| < 1σ` in this coarse fit.

## Interpretation

- The leading-order `3/13` mapping is robustly falsified.
- The corrected identity with fixed `alpha ln(10)/(4pi)` running remains
  consistent across an expanded multi-experiment anchor set.
