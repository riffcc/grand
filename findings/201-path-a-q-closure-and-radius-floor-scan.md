# 201 — Path-A `q` closure from budget + radius-floor scan

## What changed
- Added Lean closure module: `Gutoe.CTCPathAQClosure`
  - `loopGainFromBudget = budget / structuralCreationThreshold`
  - `effectiveArrivalFromBudget`
  - proofs that:
    - threshold budget gives `q_eff = 1`,
    - over-threshold budget gives `q_eff > 1`,
    - dynamic gate with positive radius gives `q_eff ≥ 1`.
- Added to `lean/lakefile.lean` roots:
  - `Gutoe.CTCPathAQClosure`
- Updated Rust bin: `ctc_path_a_time_travel_sim`
  - derives `q_eff` from SI budget/threshold:
    - `q_eff = budget_j / (kappa*(3/16)*|R|*|T|)`
  - keeps integer `q` sweep as diagnostic comparison lane.

## Build status
- `lake build Gutoe.CTCPathAQClosure` passes.
- `lake build Gutoe` passes.
- `cargo run -q -p gutoe-physics --bin ctc_path_a_time_travel_sim` passes.

## Operational check (derived mode)
- At `budget = threshold`: `q_eff = 1.0`, no pre-departure.
- At `budget = 2*threshold`: `q_eff = 2.0`, pre-departure appears (`first_pre_n = 1`).

## Radius-floor scan (current SI calibration)
Using:
- `kappa = 5.645474097135e37 J/(m*s)`
- `threshold_j = kappa*(3/16)*R*T`
- with a timelike minimum-loop model `T_min ≈ 2R/c`

Then:
- `threshold_j_min(R) = (2*kappa*(3/16)/c) * R^2`
- coefficient:
  - `2*kappa*(3/16)/c = 7.061727971907903e28 J/m^2`

Selected points:
- `R=1e-19 m -> E_min ≈ 7.06e-10 J`
- `R=8e-19 m -> E_min ≈ 4.52e-8 J`
- `R=1e-18 m -> E_min ≈ 7.06e-8 J`
- `R=1e-16 m -> E_min ≈ 7.06e-4 J`
- `R=1e-15 m -> E_min ≈ 7.06e-2 J`

Derived radii for target energy:
- `E=0.1 J -> R ≈ 1.19e-15 m`
- `E=1.5e-10 J (proton rest-energy) -> R ≈ 4.61e-20 m`

## Honest boundary
This confirms `q` is no longer free in this lane (`q_eff` is budget-closed).

The radius scan result depends on the chosen loop-time floor model (`T_min ≈ 2R/c`).
If a stronger physical lower bound on `T` is imposed, the floor shifts accordingly.
