# 190 — FTL Exhaustive Measure Lane + Lean Wall-Speed Closure

Date: 2026-02-28

## What was done

### Lean (proof-first)
1. Extended `VacuumEnergyBounds.lean` with explicit Higgs-orientation wall energetics:
- `higgsOrientationGradientEnergy = (fc * v^2 * A * Δθ^2)/(2*L)`
- inverse-thickness theorem (`E ∝ 1/L`), linear-in-area theorem, antitone-in-thickness theorem.

2. Added rear-face 1/10 factor kinematics (from `VoidRearFace`):
- `rearFaceSuppressionR = rearCostFactor = 1/10`
- rear wall tension theorem: `σ_rear = σ_front / 10`
- Lorentz-factor theorem: `γ_rear = 1 + 10*(E/A)/σ_front`
- monotonic theorem: for nonnegative drive, `γ_rear >= γ_front`.

3. Added finite-universe mechanism exhaustion module:
- `FTLExhaustion.lean`
- explicit mechanism enum + blocker classes + theorem
  `exhaustive_no_go_declared_universe : ∀ m, ¬ feasibleUnderCurrentLane m`
  over the declared mechanism universe.

Build status:
- `lake build Gutoe.VacuumEnergyBounds` ✅
- `lake build Gutoe.FTLExhaustion` ✅
- `lake build Gutoe` ✅

### Rust (measure-first)
Added:
- `crates/gutoe-physics/src/bin/ftl_exhaustive_measure.rs`

This bin computes from equations (no categorical hardcoded no-go):
1. Casimir sweep vs Higgs-restoration density.
2. Higgs wall-surf tension for front/rear faces across thickness and Δθ sweeps.
3. Required areal drive for target β values.
4. Finite-energy β sweep with analytic causal margins (`1 - β²` from gamma), to avoid float-rounding ambiguity.

Output files:
- `/tmp/bh_renders/ftl_exhaustive_measure/ftl_exhaustive_measure.txt`
- `/tmp/bh_renders/ftl_exhaustive_measure/ftl_exhaustive_measure.csv`
- `/tmp/bh_renders/ftl_exhaustive_measure/ftl_exhaustive_measure.json`

## Key measured results

From `ftl_exhaustive_measure.txt`:
- `higgs_restoration_density_j_m3 = 2.453717368943e45`
- `casimir_max_density_j_m3 = 4.333752574826e32` (at `gap=1e-15 m` in this sweep)
- `casimir_to_higgs_ratio_max = 1.766198760166e-13`
- `casimir_deficit_orders_of_magnitude = 12.752960`

Rear-face wall channel:
- `rear_face_factor = 0.1`
- finite-sweep causal margins:
  - `front_min_1_minus_beta_sq_from_gamma = 9.999999980000e-19`
  - `rear_min_1_minus_beta_sq_from_gamma = 9.999999998000e-21`
- `front_ftl_detected = false`
- `rear_ftl_detected = false`

Interpretation:
- Rear-face suppression measurably improves required areal drive (10x on tension and thus on required `E/A` for fixed target β).
- In this relativistic wall-surf model, finite-energy sweeps approach `β -> 1` but do not cross `β > 1`.

## Honesty boundary
- This lane measures the continuous wall-energy model and proves finite-universe closure over declared mechanism classes.
- It does **not** yet formalize a local discrete automorphism actuator that bypasses field-gradient dynamics.
- If such an actuator is proposed, it needs a new explicit state variable + locality law + energy functional before it can be measured/proved.
