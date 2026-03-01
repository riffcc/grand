# Finding 267 — RH Scaling/Plateau Diagnostic (Multi-Resolution)

Date: 2026-03-01  
Scope: Decide whether current RH structural lane is mainly under-resolved or map-mismatched by measuring resolution law, floor candidate, and residual curvature.

## Added

New runtime binary:

- `crates/gutoe-physics/src/bin/riemann_nail_scaling_plateau_report.rs`

Artifacts written by this lane:

- `/tmp/bh_renders/riemann_nail_scaling_plateau_report.json`
- `/tmp/bh_renders/riemann_nail_scaling_plateau_report.txt`

## Run

```bash
GUTOE_RIEMANN_REF_PATH=/tmp/bh_renders/zeta_zeros_first_1000_odlyzko.txt \
GUTOE_RIEMANN_NS=512,768,1024,1536,2048 \
cargo run -q -p gutoe-physics --bin riemann_nail_scaling_plateau_report
```

## Protocol

- Train/Hold/Freeze remains identical to hardened lane:
  - train: first 40 zeros
  - hold: next 40
  - freeze: next 40
- Same branch-locked quadratic map capacity:
  - objective on canonical core: `hold + freeze + 0.05*train + 0.001*|c|`
- This report adds:
  - scaling fit `err(n) = A * n^{-p}`
  - floor fit `err(n) = A * n^{-p} + C`
  - residual-curvature probe on first 500 points

## Results

Core hold+freeze MAPE by resolution:

- `n=512`:  `3.6221e-2`
- `n=768`:  `2.1739e-2`
- `n=1024`: `1.5540e-2`
- `n=1536`: `9.9077e-3`
- `n=2048`: `7.2690e-3`

Long unseen holdout:

- `121..500` MAPE: `0.1307 -> 0.0200` (monotone improvement)
- `121..1000` MAPE (where available): `0.1243 -> 0.0608`

Branch identity:

- Best `start` stayed at `0` for every `n` in this scan.
- Interpretation: improvement is not coming from branch drift.

Scaling law:

- No-floor fit `err = A n^{-p}`:
  - `p = 1.15459`
  - `SSE = 1.7546e-3`
- Floor fit `err = A n^{-p} + C`:
  - `C = 1.7160e-3`
  - `p = 1.31207`
  - `SSE = 1.5557e-4`
  - `SSE gain vs no-floor = 91.13%`

Residual-curvature trend (first 500 points, quadratic coeff in log-target space):

- `n=768`:  `-5.204e-2`
- `n=1024`: `-3.428e-2`
- `n=1536`: `-1.507e-2`
- `n=2048`: `-6.972e-3`

Magnitude decreases with resolution, consistent with under-resolution dominating current mismatch.

## Honest Interpretation

- RH is still **not proven**.
- This diagnostic strengthens the empirical lane:
  - strong monotone error decay with `n`,
  - no branch-hopping artifact,
  - residual curvature shrinking as resolution increases.
- A nonzero floor candidate appears in the current finite scan, but with only 5 points it is not yet decisive; higher `n` is required to classify true asymptote vs finite-range artifact.

## Next concrete step

Run same diagnostic with one or two larger matrices (e.g. `n=3072`, `n=4096` if feasible) to test whether fitted `C` stabilizes or collapses.

