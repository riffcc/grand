# 193 — CTC Classification: Global Boundary Condition vs Local Tool

## Why
After finding 192 (finite-local, unbounded-global Escher staircase), the unresolved question was:
- Is the staircase only a landscape/topology boundary condition (Path A)?
- Or does the current formal model already support compact local creation (Path B)?

## Lean additions
File: `lean/Gutoe/CTCLegality.lean`

Added:
- `nontrivialTimeShiftAtX`
- `nontrivial_shift_at_every_x`
- `sameOnTimeCylinder_not_compactly_supported`
- `current_model_global_boundary_condition`

## Result
Within the current `sameOnTimeCylinder` model:
- Nontrivial identified time-shifts exist at **every** spatial `x`.
- There is **no finite radius** outside which such nontrivial identifications vanish.

So the formal status is:
- Current model = **global boundary condition** (Path A shape).
- Current model does **not** provide a compact-support local creation mechanism (Path B not yet modeled).

## Verification
- `lake build Gutoe.CTCLegality` ✅
- `lake build Gutoe` ✅ (`8139 jobs`, zero errors)

## Honest boundary
This is a classification of the current formalization, not a no-go theorem against all possible dynamic-topology physics.
