# Finding 016: Geodesic3D Lean/Runtime Parity Closure

Date: 2026-02-24
Issue: GRAND-224

## What was verified

Lean theorem scaffold in `lean/Gutoe/Geodesic3DProjection.lean` is green and non-vacuous for the 3D pinhole reduction layer:

- `rayNormSq_eval`
- `impactRadius_even_beta`
- `rayDir_unit_normSq`
- `rayDir_z_positive`
- `kerrXi_beta_invariant`
- `kerrEta_equatorial_from_ray`

## Runtime parity wiring

Added explicit Rust tests in `crates/gutoe-gpu/src/kerr.rs` mirroring Lean Kerr image-constant statements:

- `image_plane_xi_is_beta_invariant`
- `image_plane_eta_equatorial_is_beta_squared`

Existing parity tests in `crates/gutoe-gpu/src/geodesic3d.rs` remained green:

- `lean_ray_norm_sq_eval_parity`
- `lean_impact_radius_even_beta_parity`
- `lean_ray_dir_unit_norm_parity`

## Validation commands

- `lake build Gutoe/Geodesic3DProjection.lean` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-gpu image_plane_xi_is_beta_invariant -- --nocapture` ✅
- `cargo test -p gutoe-gpu image_plane_eta_equatorial_is_beta_squared -- --nocapture` ✅
- `cargo test -p gutoe-gpu geodesic3d::tests::lean_ -- --nocapture` ✅

## Result

The 3D projection invariants and Kerr image-constant reductions are now formally represented in Lean and exercised by concrete runtime tests, preserving Lean↔Rust parity at the camera/ray reduction boundary.
