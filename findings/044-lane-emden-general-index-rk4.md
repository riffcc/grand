# Finding 044 — General-Index Lane-Emden ODE Scaffold (Lean + Rust Parity)

## Scope
Advance GRAND-281 from a fixed n=0 identity to a reusable integer-index Lane-Emden solution schema with executable ODE support.

## Lean additions (`lean/Gutoe/StellarFusion.lean`)
- Added generalized residual:
  - `laneEmdenResidualNat (n : ℕ) (ξ θ θ' θ'' : ℝ)`
- Added regularity and solution predicates:
  - `laneEmdenRegularOrigin`
  - `LaneEmdenSolutionNat`
- Added theorems:
  - `lane_emden_n0_regular_origin`
  - `lane_emden_residual_nat_n0_zero`
  - `lane_emden_n0_solution`

This keeps the n=0 exact witness and lifts it into a general-index formulation.

## Rust additions (`crates/gutoe-physics/src/equations.rs`)
- Added parity primitives:
  - `lane_emden_residual_nat`
  - `lane_emden_regular_origin`
- Added RK4 integrator for integer index:
  - `lane_emden_integrate_rk4_nat`
  - Uses regular-center expansion seed (`theta≈1-xi^2/6`, `theta'≈-xi/3`) to avoid ξ=0 singularity.

## Tests added
- `lane_emden_residual_nat_reduces_to_n0_residual`
- `lane_emden_regular_origin_accepts_exact_center`
- `lane_emden_rk4_n0_matches_closed_form_near_xi1`
- `lane_emden_rk4_n1_matches_sinxi_over_xi_near_xi1`

## Verification
- `lake build Gutoe.StellarFusion` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-physics lane_emden_ -- --nocapture` ✅

## Impact
This is the bridge from symbolic baseline to executable Lane-Emden dynamics:
- Lean now has a non-vacuous general-index schema.
- Rust now has a deterministic ODE path that can be linked to ignition constraints in follow-up work.
