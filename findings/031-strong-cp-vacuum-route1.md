# 031 — Strong CP Route-1 Vacuum Sector Gate (GRAND-267)

Status: theorem-level route-1 scaffold implemented; final nonperturbative closure still open.

## New Lean artifacts

- `lean/Gutoe/StrongCPVacuum.lean`
  - `continuous_to_fundamental_group_constant`
  - `based_continuous_to_fundamental_group_zero`
  - `based_fundamental_sector_unique`
  - `theta_unphysical_of_onlyZeroSupport`
  - `onlyZeroSupport_of_inherits_zero`
  - `cl13_theta_unphysical`
  - `cl13_route1_scale_stable`

## What this establishes

1. Using a discrete fundamental gauge carrier (`Fin 3` / Z₃ skeleton), continuous
   maps from preconnected domains are constant.
2. With fixed basepoint value, based maps collapse to the trivial map.
3. Therefore, in the route-1 support model (`Q=0` only), the CP-odd partition
   channel is identically zero and θ becomes unphysical in the finite-sector sum.
4. A separate inheritance theorem is explicit: if an effective support set is
   inherited from fundamental `{0}`, nonzero sectors cannot be populated.

## Scope boundary

This is not yet a full proof that emergent SU(3) in the full nonperturbative
theory cannot repopulate nontrivial winding sectors. That remaining step is the
core open part of GRAND-267.
