# 209 — Containment scope: law-level yes, state-level no (current projection lane)

## New theorem module
- `lean/Gutoe/ContainmentScope.lean`

## Proven
1. Every fiber has a canonical basepoint (`grade1Section x`) that projects to `x`.
2. Fiber membership is exactly:
   - `v ∈ fiberAt x` iff `v - fiberBase x ∈ ker(grade1Projection)`.
3. Every fiber is an affine translate of the same shared kernel.
4. Every fiber is nonempty.
5. Fibers over different base points are related by basepoint-delta translation.
6. Shared kernel rank remains `12` globally.

Interpretation:
- The algebraic law structure is uniform everywhere in this lane:
  each local fiber is "the same shape" (same kernel model) shifted in basepoint.

## Also proven (hard boundary)
7. No full state reconstructor from `grade1Projection` can exist:
   - there is no `rec` with `rec (grade1Projection v) = v` for all `v`.
   - reason: `grade1Projection` is non-injective.

Interpretation:
- Strong "state-level everything-from-everywhere" is not available in the
  current linear projection lane.

## Build status
- `lake build Gutoe.ContainmentScope` passed.
- `lake build Gutoe` passed (warnings only).
