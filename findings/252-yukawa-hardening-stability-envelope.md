# Finding 252 — Yukawa Hardening: Stability Envelope

## Scope
Performed a hardening pass on the decomposition outputs using:

1. Leave-one-scale-out (LOO) stability.
2. Bootstrap uncertainty bands.
3. Perturbation sensitivity checks (subset, coupling scheme, and down-shape perturbations).

Inputs:
- `/tmp/bh_renders/yukawa_full_dynamics_scan.csv`

Outputs:
- `/tmp/bh_renders/yukawa_mode_decomp_hardening.txt`
- `/tmp/bh_renders/yukawa_mode_decomp_hardening.json`

## Baseline metrics

- `Δ_QCD = mean(s²_down) - 2 = 0.440686475`
- `c1 (origin fit, gap_lg = c1 * alpha_s) = 0.997641667`
- derivative lane:
  - `d(Δ_ud)/dlnμ` vs `y_t²` slope = `0.003362222`
  - `R² = 0.999785557`

## 1) LOO stability

Across 7 folds (dropping each scale once):

- `Δ_QCD` remained in:
  - `[0.440683679, 0.440690587]`
- derivative `R²` remained high:
  - `[0.996142, 0.999793]`
- `c1` ranged:
  - `[0.952488, 1.298987]`
  - with only the lowest-scale holdout producing the largest upward excursion.

Interpretation:
- Casimir-like plateau and derivative correlation are robust.
- `c1` is sensitive to low-scale anchor inclusion (expected with short 7-point span and residual intercept).

## 2) Bootstrap bands (20,000 samples)

95% intervals:

- `Δ_QCD`: `[0.440674971, 0.440696896]`
- `c1 (origin fit)`: `[0.833385515, 1.665761266]`
- derivative slope: `[0.003292160, 0.003517813]`
- derivative `R²`: `[0.999076548, 0.999993742]`

Interpretation:
- `Δ_QCD` and derivative lane are tightly constrained.
- `c1` remains broad under nonparametric resampling due to short dataset and intercept leakage.

## 3) Perturbation sensitivity

### 3a. High-μ subsets

- Full:
  - `Δ_QCD = 0.440686475`, `c1 = 0.997641667`
- `μ >= 1e8`:
  - `Δ_QCD = 0.440694899`, `c1 = 1.633646149`
- `μ >= 1e12`:
  - `Δ_QCD = 0.440698371`, `c1 = 1.735435638`
- `μ >= 1e16`:
  - `Δ_QCD = 0.440700965`, `c1 = 1.733325616`

Interpretation:
- `Δ_QCD` is stable.
- `c1` inflates on truncated high-μ subsets, consistent with fitting through-origin while a residual intercept is present.

### 3b. Alpha_s scheme-style transforms

Scenarios (`alpha_s -> a*alpha_s + b`) preserved `R²=0.879701` in free affine fit.

Origin-fit `c1`:
- `a=0.95`: `1.050149123`
- `a=1.05`: `0.950134921`
- `b=+0.002`: `0.973873642`
- `b=-0.002`: `1.021548680`

Interpretation:
- `c1 ~ 1` remains in moderate scheme perturbations.

### 3c. Top-Yukawa scheme scaling

Under `y_t -> k y_t`:
- derivative `R²` stayed fixed at `0.999786`.
- slope rescaled as expected (`~1/k²`), confirming structural correlation and normalization sensitivity separation.

### 3d. Down-shape perturbation

Applied shape tilt to down masses:
- `m_d -> m_d(1+eps)`, `m_s` fixed, `m_b -> m_b(1-eps)`, `eps ∈ {±0.01, ±0.02}`.

Results:
- `Δ_QCD` shifted between `0.426866881` and `0.454232271`.
- derivative lane unchanged (`R²=0.999786`).

Interpretation:
- Casimir-proximity claim is sensitive to down-sector mass-shape systematics (as expected).
- derivative Yukawa lane remains orthogonal and robust.

## Overall verdict

Hardening confirms:

- **Strong**: derivative Yukawa-integral structure (`R² ~ 0.9998`) is robust to all tested perturbations.
- **Strong**: down-sector plateau exists and is numerically stable (`Δ_QCD ~ 0.44069`).
- **Qualified**: `c1 ~ 1` for `gap_lg ~ alpha_s` is plausible at full range but not yet a tightly constrained universal constant under truncation/bootstrap.

This places the decomposition on a stable empirical footing while identifying exactly where additional data depth (more scan points / refined residual modeling) is needed.
