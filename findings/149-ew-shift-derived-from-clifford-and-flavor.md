# 149 - EW M_Z Shift Derived from Clifford + Flavor Structure

Superseded by finding 150 for the frozen publishable state.

## Change
Replaced fixed EW bridge constant (`+4.51e-4`) with a derived structural
expression coupled to neutrino texture anisotropy.

File:
- `crates/gutoe-physics/src/dynamics_map.rs`

## New expression
`sin²(theta_W)(M_Z) = sin²(theta_W)_structural + Δsin²`

with
`Δsin² = α² * (d/2 + ρ_tex/2)`

where:
- `α = 1/137` (structural leading value)
- `d = 16` (Cl(1,3) dimension)
- `ρ_tex = |t3 - t2| / |t2 - t1|` from neutrino texture eigenvalue spacings
  (this is a texture anisotropy ratio, not the oscillation `Δm²` ratio)

No observed `sin²(theta_W)(M_Z)` value is used as an input to this shift.

## Current numbers (from coupled gate artifact)
From `/tmp/bh_renders/ew_flavor_coupled_ci_gate.json`:

- `sin2_structural = 0.230769230769`
- `alpha = 0.007299270073`
- `clifford_half_dim = 8.0`
- `flavor_anisotropy = 0.904222214321`
- `shift_coeff = 8.452111107161`
- `delta_sin2 = 4.503229318110e-4`
- `sin2_mz_bridge = 0.231219553701`
- `mz_abs_err = 4.462989582343e-7` vs target `0.23122`

## Status
- EW side remains in-window (`ew_ok=true`).
- Flavor splitting ratio remains far from target (`ratio_ok=false`).

This confirms the closure bottleneck is in flavor mass splitting structure,
while EW bridge is now produced from internal structural terms rather than a
fixed additive constant.
