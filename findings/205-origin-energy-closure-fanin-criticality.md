# 205 — Origin energy closure requires near-critical fan-in gain

## New probe
- `ctc_origin_energy_closure_probe`

Outputs:
- `/tmp/bh_renders/ctc_origin_energy_closure_probe/ctc_origin_energy_closure_probe.txt`
- `/tmp/bh_renders/ctc_origin_energy_closure_probe/ctc_origin_energy_closure_probe.json`

## Question
Given seed contribution `E_seed` and target origin energy `E_target`, what
effective fan-in gain is required under finite causal depth?

Definitions:
- `E_seed` (default quark-scale threshold): `4.569939e-8 J`
- `E_target` (default observable-universe order): `1e69 J`
- ratio: `E_target / E_seed ≈ 2.188e76`
- causal depth estimate:
  - `K_max = floor(horizon / period)`
  - default `period = 5.366854e-27 s`, `horizon = 4.354e17 s`
  - `K_max ≈ 8.113e43`

Effective gain:
- `b_eff = branching * merge_fraction * eta * infra_gain`

## Main result
With structural merge default `merge_fraction = 3/16` and `branching = 2`:
- `b_eff ≈ 0.375`
- `finite_horizon_reaches_target = false`

Required minimum for finite-horizon closure:
- `b_min ≈ 1.0000000000000102`

So the origin-energy closure lane is **critical**:
- below 1: cannot reach target at defaults,
- just above 1: target becomes reachable within finite depth,
- much above 1: runaway fan-in.

## Cross-checks
If merge is unconstrained (`merge_fraction = 1`, same other defaults):
- `b_eff ≈ 2.0`
- finite-horizon closure is reachable (but in runaway-prone regime).

User suppression candidate (`G_uncapped * (13/16) * (1/10)`):
- `G_uncapped = 1.9992`
- `b_eff = 1.9992 * 0.8125 * 0.1 = 0.162435`
- `finite_horizon_reaches_target = false`

If interpreted as infinite geometric closure with `B < 1`:
- `E_origin = E_seed / (1 - b_eff)`
- with `E_seed = 4.569939e-8 J`, `E_origin ≈ 5.45622e-8 J` (about `1.19x` seed)
- i.e. strongly convergent, no large amplification.

## Infinite closed geometric note
For an infinite closed geometric sum with `B < 1`, required epsilon is:
- `epsilon = E_seed/E_target ≈ 4.57e-77`
- i.e., `B = 1 - epsilon` (numerically rounds to `1.0` in f64).

Interpretation:
- If the origin is a fan-in closure, the algebra needs a mechanism that pins
  effective gain near criticality (`b_eff ~ 1`) while preventing runaway.
