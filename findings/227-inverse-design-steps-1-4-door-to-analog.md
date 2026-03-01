# Finding 227 — Inverse-Design Steps 1-4 (Door -> Analog)

## Summary

Implemented a dedicated inverse-design lane that executes the requested
pipeline directly from the current campaign output:

1. Minimal door invariants.
2. Sensitivity and ablation.
3. Dimensionless target windows.
4. Candidate physical analog mappings.

Binary:

- `crates/gutoe-physics/src/bin/ctc_inverse_design_steps_1_4.rs`

## Run

```bash
cargo run -p gutoe-physics --bin ctc_inverse_design_steps_1_4
```

Default input:

- `/tmp/bh_renders/ctc_50y_campaign_fast/ctc_50y_campaign_fast.json`

Outputs:

- `/tmp/bh_renders/ctc_inverse_design_steps_1_4/ctc_inverse_design_steps_1_4.txt`
- `/tmp/bh_renders/ctc_inverse_design_steps_1_4/ctc_inverse_design_steps_1_4.json`

## Step 1 — Minimal Invariants (Extracted)

From current baseline row:

- `beta = 0.8` (timelike local lane holds)
- `s = 0.1` -> `chi = u/c = 10`
- `q_eff = 3.63e19`
- `mu = n*(q_eff-1) = 4.35e21` (for `n=120`)
- `g = branching * f_void * eta * infra = 1`
- `|g-1| = 0`

Gate checks:

- `local_timelike = true`
- `coordinate_superluminal = true`
- `budget_gate_open = true`
- `predeparture = true`

## Step 2 — Sensitivity + Ablation

Local sensitivity (`±5%`) confirms structure:

- `s` only changes `chi = 1/s` (coordinate channel strength).
- `{kappa, f_void, radius, period}` all enter via threshold and therefore scale
  `q_eff` and predeparture margin inversely.
- `budget_per_door` and `n_loops` scale predeparture margin directly.
- `beta` is primarily a timelike/safety margin knob.

Key ablation results:

- `s -> 1` kills coordinate superluminal (`viable = false`).
- `beta -> 1` kills local timelike condition (`viable = false`).
- `budget -> 0.9 * threshold` closes gate (`viable = false`).
- `q_eff = 1` is break-even and not predeparture (`viable = false`).
- `n = 0` removes loop accumulation (`viable = false`).

This cleanly separates required conditions from optional ones.

## Step 3 — Dimensionless Target Windows

Grid sweep (`beta x s x q x n = 600` points):

- passing: `400/600` (`0.667`)
- robust subset (extra margin cut): `400/600` (`0.667`)

Robust windows:

- `beta_local in [0.6, 0.95]`
- `s in [0.05, 0.4]`
- `q_eff in [1.1, 2.0]`
- `n in [20, 200]`

Recommended Pi target point (best objective in robust set):

- `Pi1 = q_eff = 2.0`
- `Pi2 = chi = u/c = 20.0`
- `Pi3 = beta = 0.6`
- `Pi4 = mu = n(q-1) = 200`
- `Pi5 = g = 1.0`
- `Pi6 = |g-1| = 0.0`
- `Pi7 = stability_margin = 0.4`

PDG anchor lock (2025 constants embedded in probe):

- `alpha^-1`: model `137` vs PDG `137.035999177` (`-0.02627%`)
- `sin^2(theta_W)(M_Z)_MS`: model `3/13 = 0.230769...` vs PDG `0.23122`
  (`-0.19495%`)
- `N_nu`: model `3` vs PDG `2.9963` (absolute `+0.0037`)

These anchors are exported under:

- `step_3_dimensionless_targets.pdg_ratio_anchors_2025`

## Step 4 — Physical Analog Classes (Mapped)

Generated analog mappings parameterized by the same Pi groups:

1. Phase-locked photonic delay mesh.
2. Superconducting resonator loop (flux/phase domain).
3. RF ring with digital-twin feedback controller.
4. Analog-gravity metamaterial loop.

Each class includes:

- Pi-to-control equation mapping,
- control knobs,
- lab observables,
- target Pi vector.

## Notes

- This lane is inverse-design from simulation constraints to measurable analog
  control objectives.
- It does **not** claim a physical engine exists.
- It does provide a falsifiable, parameterized bridge from in-engine theorem
  structure to real-lab analog experiments.
