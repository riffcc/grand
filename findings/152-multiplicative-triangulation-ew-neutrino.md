# GRAND — Multiplicative Triangulation (EW + Neutrino)

Date: 2026-02-28

## Scope

Implemented a concrete triangulation lane that forces latent multiplicative factors from independent anchors:

- ratio anchor: `Δm²32 / Δm²21`
- absolute anchors: `Δm²21`, `Δm²32`
- cross-lane anchor: `sin²θ_W(M_Z)`

Added Lean formal spine for product+ratio identifiability and multiplicative EW factorization.

## Files

- `crates/gutoe-physics/src/bin/triangulate_params.rs`
- `lean/Gutoe/MultiplicativeTriangulation.lean`
- `lean/lakefile.lean` (added new Lean root)

## Runtime Outputs

From `/tmp/bh_renders/triangulate_params_report.txt`:

- `p_ratio = 13.688110433760`
- `kappa_geo = 34.660950672521`
- `ew_coeff_required = 8.460487692308`
- `ratio_fit_rel_err = 1.768048732288e-11`
- `kappa_vs_structural_rel = 5.354507623295e0`
- `ew_coeff_delta_rel = 5.756096153845e-2`

Structural references:

- `p_structural = 13.7`
- `kappa_structural = 60/11 = 5.454545...`
- `ew_coeff_structural = d/2 = 8`

Interpretation:

- Pattern constraint (`Δm²` ratio) is satisfiable near the existing exponent.
- Absolute neutrino scale requires a much larger multiplicative factor than current structural `60/11`.
- EW shift needs a modest uplift over `d/2`.

## Lean Validation

`lake build Gutoe` passes (build green after integrating `Gutoe.MultiplicativeTriangulation`).

Key theorem content:

- recovery of positive latent factors from product+ratio anchors
- identifiability of `(θ₁, θ₂)` under equal anchors
- EW shift multiplicative form and closed numeric form

## Honesty

This lane exposes parameter tension; it does not claim closure.

- It is a forcing diagnostic: where one shared multiplicative chain fits, and where scale drift remains.
- No per-lane fallback/retune paths were introduced.
