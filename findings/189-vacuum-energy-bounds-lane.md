# 189 — Vacuum Energy Bounds Lane (Ford-Roman + Casimir + Cl(1,3) EW Proxy)

Date: 2026-02-28

## Scope
Added a new Lean research lane:
- `lean/Gutoe/VacuumEnergyBounds.lean`

Goal: formalize the constraint skeleton around local negative-energy proposals, without introducing vacuous closure claims.

## What Was Added

### 1) Structural constants carried into the lane
From existing Cl(1,3) modules:
- `voidFractionQ = 3/16`
- `higgsQuartic = 13/100`
- `vevOverProton = 40/153`

New theoremized proxy:
- `ewsbBarrierProxyQ = λ_H * (v/mp)^4 / 4`
- closed form: `(13/100) * ((40/153)^4) / 4`
- positivity theorem: `0 < ewsbBarrierProxyQ`

### 2) Ford-Roman-style inequality formalization
Defined:
- `FordRomanBound (rhoNeg tau qeiK) : Prop := |rhoNeg| * tau^4 ≤ qeiK`

Proved:
- `ford_roman_no_durable_window`:
  if target magnitude and minimum dwell time imply budget exceedance (`qeiK < rhoTarget * tauMin^4`), a durable window is impossible.

Interpretation: this captures the standard “borrow more, hold shorter” constraint logic as a formal no-go gate.

### 3) Casimir geometric suppression formalization
Defined:
- `casimirMagnitude (kappa a) = kappa / a^4`

Proved:
- `casimir_magnitude_antitone_in_gap`:
  for `kappa ≥ 0`, increasing gap `a` decreases attainable magnitude.
- `casimir_no_go_from_min_gap`:
  if the minimum achievable gap already undershoots a target, all larger gaps undershoot too.

Interpretation: this is the rigorous finite-gap `a^-4` bottleneck statement.

## Build Verification
- `lake build Gutoe.VacuumEnergyBounds` ✅
- `lake build Gutoe` ✅ (full library build)

No `sorry` added.

## Honesty Boundary (Explicit)
This lane proves bound structure, not a propulsion mechanism.

Still open and intentionally not claimed:
- whether a *physical* local void-orientation field exists beyond basis relabeling,
- whether any such field can evade the Ford-Roman + finite-gap suppression gates,
- any macroscopic FTL-enabling conclusion.

## Why This Matters
This turns an otherwise open-ended speculative discussion into machine-checkable constraints:
- if a proposed channel is only a basis relabeling, it has no physical effect,
- if it is a physical local field, it must satisfy hard duration/magnitude and geometric suppression bounds.

That gives a clean decision surface for future proposals and prevents silent drift into unconstrained claims.
