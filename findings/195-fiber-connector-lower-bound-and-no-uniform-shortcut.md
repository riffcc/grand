# 195 — Fiber Connector Geometry: Lower Bound and No Uniform Shortcut

## Prompt target
Given fibers over two spacetime points under the 16->4 projection, determine whether
connector lengths in total space can be bounded independently of base separation.

## Lean file
- `lean/Gutoe/ProjectionFibers.lean`

## New theorems

- `fiberAt`
- `abs_coord_le_norm`
- `fiber_coord0`
- `fiber_connector_lower_bound_axis0`
- `no_uniform_fiber_connector_bound`

## Formal result (current natural model)
For any base points `x,y : Fin 4 -> ℝ` and any states `v ∈ fiberAt x`, `w ∈ fiberAt y`:

`|y 0 - x 0| <= ‖w - v‖`.

Consequently, for every bound `B`, there exist fibers such that **all** connectors between them
have norm `> B`.

So in this lane there is no distance-independent global shortcut bound.

## Interpretation vs requested trichotomy
This rules out the "bounded independent of D" outcome in the present formalization.
It supports the rigid side (at least axis-wise linear lower bound with base separation).

## Verification
- `lake build Gutoe.ProjectionFibers` ✅
- `lake build Gutoe` ✅ (`8140 jobs`, warning-only)

## Honest boundary
This theorem gives a rigorous lower bound and a no-uniform-bound result.
It does not yet prove an exact infimum formula over all fiber pairs under a full
Minkowski+fiber connection model.
