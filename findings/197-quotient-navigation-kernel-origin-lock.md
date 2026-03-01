# 197 — Quotient navigation probe: kernel multiplicative actions keep origin fixed

## Goal
Try the next open lane: full quotient Cl(1,3) multiplication (anticommutation + metric)
for recursive lift/traverse/descend, and test whether kernel-only multiplicative routing
can generate a 4D point translation from origin.

## Lean lock (new)
Module: `Gutoe.RecursiveNavigationNoTranslation`

Added structural theorems for the current `256 -> 16 -> 4` descent map:
- `descended_linear_step_preserves_origin`
- `no_linear_origin_to_nonzero_target`
- `linear_origin_reachability_iff_zero`

These close the linear/multiplicative homogeneous lane in proof form.

## Rust measurement lane (new)
Bin: `recursive_z3_quotient_navigation_probe`

Mechanism tested:
- Explicit Cl(1,3) geometric product on basis blades (signature `+---`)
- Action: `grade1(K * x * reverse(K))`
- `K` restricted to kernel masks `{0,3,5,6,7,9,10,11,12,13,14,15}`

Outputs:
- `/tmp/bh_renders/recursive_z3_quotient_navigation_probe/recursive_z3_quotient_navigation_probe.txt`
- `/tmp/bh_renders/recursive_z3_quotient_navigation_probe/recursive_z3_quotient_navigation_probe.json`

## Measured results
- Challenge separation: `|target - origin| = 1e6`
- Basis-kernel scan:
  - `max_origin_shift = 0.0`
- Random kernel multivector scan (`50,000` samples):
  - `max_origin_shift = 0.0`
  - `origin_translation_detected = false`

## Conclusion
In this quotient-product kernel-only multiplicative lane,
origin stays fixed after grade-1 descent. No 4D point translation was detected.

This does **not** close nonlinear/topological creation lanes; it closes the tested
homogeneous multiplicative mechanism.
