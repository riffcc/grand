# 196 — Recursive Z3 metric probe (linear tower): no bounded point shortcut

## Scope
Test the hypothesis: "lift to higher level (256D), traverse, descend" may reduce effective 4D point-to-point distance.

We measured this in the current **linear projection tower model**:
- `256 -> 16`: slice `j=0`
- `16 -> 4`: first four coordinates
- total-space metric: Euclidean

## Lean status
`Gutoe.RecursiveZ3Tower` now compiles cleanly (no `sorry`) and locks:
- `z3_16_order3`
- `z3_256_index_order3`
- `proj256to16_surjective`
- `proj256to16_kernel_finrank = 240`
- tower profile with `ker(16->4)=12`, total hidden dims `252`

## Rust measurement lane
New bin: `recursive_z3_metric_probe`.
Outputs:
- `/tmp/bh_renders/recursive_z3_metric_probe/recursive_z3_metric_probe.txt`
- `/tmp/bh_renders/recursive_z3_metric_probe/recursive_z3_metric_probe.json`

### Results (20k random samples/case)
- `unit_axis`: `d4 = 1.0`, witness `= 1.0`, random min `= 11.208669...`
- `mixed_13`: `d4 = 13.0`, witness `= 13.0`, random min `= 17.121798...`
- `large_scale`: `d4 = 1e6`, witness `= 1e6`, random min `= 1_000_000.000061...`

Global verdict:
- `infimum_equals_base_distance = true`
- `bounded_shortcut_detected = false`

## Interpretation
In this linear additive tower model, lifting to 256D and traversing fibers does **not** compress base point distance: infimum connector length equals base 4D separation.

This does **not** close the full recursive hypothesis; it closes only the linear model. The remaining open lane is quotient/anticommutation-induced navigation (multiplicative, not additive).
