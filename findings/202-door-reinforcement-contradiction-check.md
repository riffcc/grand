# 202 — Door reinforcement contradiction check (closed vs open cycle)

## Request
"Try it and prove it anyways."

This pass adds a dynamic door-state bookkeeping lane and tests the exact
self-funding contradiction condition.

## New Lean module
- `Gutoe.CTCDoorEnergyBookkeeping`

Definitions:
- `LoopConservation`:
  - `Ein + Eprev = Eout + Enext + Export + Loss`
- `ClosedPacketCycle`:
  - `Ein = Eout`
- `NoDoorDrawdown`:
  - `Eprev ≤ Enext`

Theorems:
- `closed_cycle_no_positive_export`
  - under closed packet cycle + no drawdown + nonnegative loss, `Export ≤ 0`
- `export_equals_inflow_plus_drawdown_minus_loss`
  - exact rearrangement:
  - `Export = (Ein - Eout) + (Eprev - Enext) - Loss`
- `positive_export_requires_inflow_or_drawdown`
  - if `Export > 0` and `Loss ≥ 0`, then either net inflow (`Ein > Eout`) or drawdown (`Enext < Eprev`)

Status:
- `lake build Gutoe.CTCDoorEnergyBookkeeping` passes
- full `lake build Gutoe` passes

## New Rust bin
- `ctc_door_reinforcement_probe`

Outputs:
- `/tmp/bh_renders/ctc_door_reinforcement_probe/ctc_door_reinforcement_probe.txt`
- `/tmp/bh_renders/ctc_door_reinforcement_probe/ctc_door_reinforcement_probe.json`

### Result A — contradiction guard (closed-cycle sweep)
- Guard tested:
  - `Ein=Eout`, `Enext>=Eprev`, `Loss>=0`
- Sweep result:
  - `positive_export_count = 0`
  - `max_export_j = -0.0`

So no positive export survives the closed/no-drawdown/nonnegative-loss regime.

### Result B — open-flux run
- Positive export appears when packet-energy inflow is present.
- Conservation residual is numerically zero (machine precision).
- `theorem_guard_violations = 0`

Interpretation:
- The model allows export with throughput, but not "free" export in a closed cycle.
- Any positive export traces to net inflow and/or door drawdown, exactly as the theorem states.

## Additional fix
- Removed artificial period clamps that were hiding quark-scale tests:
  - `ctc_path_a_time_travel_sim`: `T >= 1e-30` (zero-guard only)
  - `ctc_door_reinforcement_probe`: `T >= 1e-30` (zero-guard only)
- Added `n_required_predeparture` reporting in `ctc_path_a_time_travel_sim`.

Quark-scale check (`R≈8.04e-19 m`, `T≈5.37e-27 s`):
- threshold remains tiny (`~4.57e-8 J`)
- but required loops for predeparture are enormous (`~1.55e21` at `q=2`)

So tiny patch-energy alone is not sufficient; the loop-count/time-shift scaling is the choke point.
