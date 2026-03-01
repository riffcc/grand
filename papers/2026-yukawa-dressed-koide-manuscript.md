# Dressed Koide Structure Across SM Sectors from Full-Dynamics GUTOE Lanes
## Draft Manuscript (Working)
### 2026-03-01

## Abstract
We report a constrained decomposition of Yukawa-sector structure in a Cl(1,3)-based GUTOE runtime lane with coupled full dynamics (2-loop QCD running, running QED, flavor-sensitive mass flow, threshold activation, and Yukawa self terms). Three empirical relations survive hardening: (i) a down-sector plateau shift relative to the lepton baseline consistent with a Casimir ratio (`s²_down - 2 ≈ 4/9` at sub-percent level), (ii) near-perfect derivative-level correlation between up/down split growth and top-Yukawa strength (`dΔ_ud/dlnμ` vs `y_t²`, `R²≈0.9998`), and (iii) monotonic UV approach of geometric ladder Koide parameter toward triplet value (`s²(L_g) → 3` from below). A quadratic gap model in `α_s` has coefficients that numerically lock to rational group-theory candidates (`23/12`, `-98/9`) within O(10^-4) relative error on this lane. We present these as falsifiable structural candidates, not final proofs.

## 1. Setup
### 1.1 Runtime lane
Data comes from:
- `crates/gutoe-em/src/bin/yukawa_full_dynamics_scan.rs`
- scan anchors: `μ = {m_t, 1e4, 1e8, 1e12, 1e16, 1e18, 1e19} GeV`
- outputs under `/tmp/bh_renders/`.

The lane integrates:
- coupled `α_s(μ), α_em(μ)` RK4 transport,
- 2-loop QCD beta structure (active-flavor thresholds),
- flavor-sensitive mass flow (`QCD + QED + Yukawa self`),
- extraction of ladder/isospin modes:
  - `L_g = [sqrt(m_u m_d), sqrt(m_c m_s), sqrt(m_t m_b)]`
  - `S_g = [0.5 ln(m_u/m_d), 0.5 ln(m_c/m_s), 0.5 ln(m_t/m_b)]`.

### 1.2 Koide parameterization
For a mass triplet, define:
- `K = Σm_i / (Σ√m_i)^2`
- `s²_eff = 6K - 2`.

This maps `K` to the Z3-harmonic amplitude parameter used throughout the lane.

## 2. Main Results
### 2.1 Down-sector plateau (Casimir candidate)
From hardened fit:
- `mean(s²_down) = 2.440686475`
- `Δ_QCD := s²_down - 2 = 0.440686475`.

Compared to:
- `4/9 = C_F/N_c = (4/3)/3 = 0.444444444`
- mismatch: `-0.003757970` (`-0.846%`).

This is a strong structural candidate that the down-sector dressing is Casimir-structured.

### 2.2 Up/down split growth tracks top-Yukawa channel
Define:
- `Δ_ud = s²_up - s²_down`.

Segment derivative regression:
- `d(Δ_ud)/dlnμ` vs midpoint `y_t²`
- fit `R² = 0.999786`, `RMSE ~ 6.6e-6`.

Interpretation: split growth is dominantly an anomalous-dimension-like Yukawa channel, not a static offset.

### 2.3 Geometric ladder approaches triplet UV value
Across the scan:
- `s²(L_g): 2.9175 -> 2.9684`
- `3 - s²(L_g)` decreases monotonically,
- fixed-`s²=3` closure RMS improves monotonically.

This supports a UV triplet emergence in the geometric mean channel.

## 3. Gap Law in Strong Coupling
Define:
- `gap_lg = 3 - s²(L_g)`.

### 3.1 Baseline forms
- origin linear: `gap = c1 α_s`, `c1=0.9976`, but broad under hardening.
- affine linear: nonzero intercept present.
- quadratic no-intercept: substantial residual reduction.

### 3.2 Rational-lock candidate
Free quadratic fit:
- `gap = c1 α_s + c2 α_s²`
- `c1=1.916340718`, `c2=-10.888164740`.

Structural candidate:
- `c1=23/12`, `c2=-98/9`.

Forced-vs-free fit quality:
- `RMSE_free = 0.003042030`
- `RMSE_forced = 0.003042064`
- `ΔRMSE ≈ 3.35e-8` (negligible at current depth).

This is consistent with a beta-structured coefficient interpretation, but requires densified scans for lock confirmation.

## 4. Hardening Summary
Hardening artifacts:
- `/tmp/bh_renders/yukawa_mode_decomp_hardening.txt`
- `/tmp/bh_renders/yukawa_mode_decomp_hardening.json`
- findings 252, 253, 254.

Stable lanes:
- `Δ_QCD` plateau value is tightly stable under LOO and bootstrap.
- derivative Yukawa lane remains robust under coupling-scale perturbations.

Less stable lane:
- single-parameter origin fit coefficient (`c1`) remains subset/scheme sensitive due to residual intercept and short anchor set.

## 5. Neutrino Lane Status
Current in-engine neutrino tiny-mass lane:
- predicts tiny normal hierarchy but fails oscillation-splitting gate magnitudes in current implementation.

Independent constraint check:
- imposing `K_ν = 1/2` with standard oscillation splittings admits physically valid NO and IO mass solutions.

So the neutrino lane is open and testable, not closed.

## 6. Interpretation Boundary
What is supported now:
- Yukawa structure can be decomposed into UV-representation anchor + vacuum/Casimir dressing + Yukawa-driven split evolution.

What is not yet proven:
- exact coefficient-lock of all gap terms as universal group invariants at finite scan depth,
- full derivation of neutrino sector and absolute mass anchor inside one closed structural pipeline.

## 7. Falsification Path
Next high-value tests:
1. Scan densification across thresholds (especially around flavor transitions).
2. Scheme conversion envelope with explicit uncertainty propagation.
3. Joint fit including third-order `α_s` term with theory priors on coefficient families.
4. Neutrino lane closure using `K_ν=1/2` prior + oscillation constraints + cosmology bounds.

## 8. Conclusion
This run does not finish the Yukawa problem, but it materially compresses it: the dominant structure is now expressed in three sharpened equations with robust empirical support and explicit failure envelopes. The remaining work is no longer unconstrained fitting; it is targeted closure on a narrow residual manifold.
