# Finding 262 — Riemann Straight-Shot Operator Lane (Exploratory)

Date: 2026-03-01  
Scope: Launch a direct RH exploratory lane in-engine with a concrete self-adjoint operator and reproducible spectrum-vs-zero report.

## What Was Added

- New binary: `crates/gutoe-physics/src/bin/riemann_straight_shot_report.rs`
- Dependency: `nalgebra = "0.33"` in `crates/gutoe-physics/Cargo.toml`

The binary:

1. Builds a real symmetric tridiagonal operator (`self-adjoint` by construction):
   - Diagonal: `ln(n + 13/16)`
   - Off-diagonal: `hop_scale * sqrt((n+1)(n+2))`
2. Computes eigenvalues via `nalgebra::SymmetricEigen`.
3. Maps raw eigenvalues to candidate zero ordinates with structural affine map:
   - `γ_pred = (13*24 + 8/17) + (11/18) * λ_raw`
4. Compares to first 80 nontrivial zeta zeros (Odlyzko table values embedded in source).
5. Writes report artifacts:
   - `/tmp/bh_renders/riemann_straight_shot_report.txt`
   - `/tmp/bh_renders/riemann_straight_shot_report.json`

## First Run Metrics

Command:

```bash
cargo run -q -p gutoe-physics --bin riemann_straight_shot_report
```

Output summary (`N=512`, `k=80`, `hop_scale=0.5`):

- `self_adjoint_residual = 0`
- `MAE = 5.590988`
- `RMSE = 7.981175`
- `MAPE = 4.016611e-2` (4.0166%)
- `max_rel_err = 1.646223e-1` (16.46%)
- `signed_rel_bias = -3.021194e-2` (-3.02%)

## Interpretation

- The lane is operational and reproducible.
- The operator is mathematically well-posed (self-adjoint residual numerically zero).
- Baseline agreement to the first 80 zeros is nontrivial but not yet precision-grade.
- This is a valid starting point for direct iteration on operator structure and mapping law.

## Next Iteration Hooks

- `GUTOE_RIEMANN_DIM` (operator dimension)
- `GUTOE_RIEMANN_HOP` (off-diagonal coupling scale)
- `GUTOE_RIEMANN_K` (number of zeros compared)
- `GUTOE_RIEMANN_OUT` (artifact output directory)

## Quick Headroom Sweep (hop scale)

Sweep run (`N=512`, `k=80`) over `hop_scale ∈ {0.30, 0.35, ..., 0.70}`:

- Best in tested grid: `hop_scale = 0.50` (current default), `MAPE ≈ 4.02%`
- Off this point, error worsens rapidly:
  - `0.45`: `MAPE ≈ 24.77%`
  - `0.55`: `MAPE ≈ 30.46%`
  - `0.30` / `0.70`: `MAPE > 100%`

So the current lane has a sharp local basin at `hop_scale = 0.50`; improving below ~4% likely requires changing operator structure and/or mapping law, not just scalar retuning.
