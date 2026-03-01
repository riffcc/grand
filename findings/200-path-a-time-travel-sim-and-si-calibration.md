# 200 — Path-A simulation + SI calibration for creation threshold

## Delivered
You asked for:
1. SI calibration of the structural threshold.
2. A non-toy time-travel simulation for Path A.

Both are now wired and executed.

## New Lean parity module
- `Gutoe.CTCPathAEffectiveArrival`
  - `effectiveArrival`, `coverArrival`
  - `effective_equals_cover_minus_shift`
  - `q_one_no_coordinate_gain`
  - `q_gt_one_predeparture_possible`

This keeps the Rust formula lane tied to theorem statements.

## SI calibration (new Rust bin)
- Bin: `dynamic_topology_creation_si_calibration`
- Output:
  - `/tmp/bh_renders/dynamic_topology_creation_si_calibration/dynamic_topology_creation_si_calibration.txt`
  - `/tmp/bh_renders/dynamic_topology_creation_si_calibration/dynamic_topology_creation_si_calibration.json`

Calibration model:
- Gate units: `(3/16)|R||T|`
- SI mapping via derived worldsheet energy:
  - `E_sheet = sigma_rear * (2πR) * (cT)`
  - `sigma_rear` from structural EW gradient tension (`fc=3/16`, `v=245.3 GeV`, `Δθ=π`, `l=ħc/v`)

Computed constants:
- `sigma_front = 5.619544567498e28 J/m²`
- `sigma_rear  = 5.619544567498e27 J/m²`
- `kappa = 5.645474097135e37 J/(m·s)`

Human-scale examples (this calibration):
- `R=1 m, T=1 s -> E ≈ 1.0585e37 J`
- `R=10 m, T=60 s -> E ≈ 6.3512e39 J`
- `R=20 m, T=50 s -> E ≈ 1.0585e40 J`

## Path-A time-travel simulation (new Rust bin)
- Bin: `ctc_path_a_time_travel_sim`
- Output:
  - `/tmp/bh_renders/ctc_path_a_time_travel_sim/ctc_path_a_time_travel_sim.txt`
  - `/tmp/bh_renders/ctc_path_a_time_travel_sim/ctc_path_a_time_travel_sim.json`

Model:
- Timelike start->node and node->goal segments at `βc`
- `n` local loops of period `T`
- Branch-shift quanta per loop `q`
- Effective arrival formula:
  - `t_eff = dt_in + dt_out + n(1-q)T`

Default-run result:
- `q=1`: no coordinate gain, no pre-departure
- `q>=2`: pre-departure appears (`first_pre_n=1`), while local segments remain timelike

Interpretation:
- In this explicit Path-A branch-shift model, pre-departure coordinate arrival is operationally reachable for `q>1`.
- This is a model result. Physical mechanism for `q>1` branch-shift realization remains open.

## Build state
- New modules compiled.
- Full `lake build Gutoe` passes.
