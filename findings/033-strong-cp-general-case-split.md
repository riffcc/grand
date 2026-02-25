# 033 — Strong CP General-Case Split (GRAND-267 Follow-Through)

Status: formal split theorem added and Lean-verified.

## New Lean module

- `lean/Gutoe/StrongCPGeneralCases.lean`

Key theorems:

- `theta_phase_at_zero`
- `theta_phase_nontrivial_of_nonzero_charge`
- `theta_physical_of_nonzero_sector`
- `gutoe_theta_unphysical_on_emergent_image`
- `strong_cp_general_case_split`

## What this locks down

This module formalizes both global regimes side-by-side:

1. **GUTOE route-1 emergent-image regime**  
   On the concrete `Z₃ → SU(3)` emergent image with anchored topological charge,
   `θ` is unphysical (`θ`-phase equality for all `θ1, θ2`).

2. **Standard nonzero-sector regime**  
   If any nonzero integer topological sector is physically accessible, then
   `θ` is physical: there exist `θ1, θ2` with distinct phase factors.

The split is explicit in `strong_cp_general_case_split`.

## Interpretation

This does not claim "all QCD has θ = 0".  
It proves a clean formal dichotomy:

- route-1 GUTOE physical-image sector ⇒ `θ` unphysical,
- nonzero topological-sector accessibility ⇒ `θ` physically observable.

That is the exact boundary needed to keep the claim honest and falsifiable.
