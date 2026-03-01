# 210 — Topological access operator: active lane opens, linear lane remains closed

## New theorem module
- `lean/Gutoe/TopologicalAccessOperator.lean`

## Wiring change
- Added `Gutoe.TopologicalAccessOperator` to `lean/lakefile.lean` roots.

## Proven in Lean
1. **Controllable entry+exit witness (conditional on dynamic gate)**:
   - If `dynamicCreationGate budget radius period` holds, there exist events
     with local nontrivial identification shift (`b.t ≠ a.t`) and an affine
     witness that hits any chosen descended 4D target.
2. **Linear no-go still holds**:
   - For nonzero descended targets, purely linear origin mapping remains blocked.
3. **Bypass is explicit and scoped**:
   - For nonzero targets, the linear lane fails while the topological/affine
     lane has a witness when the dynamic gate is open.
4. **Closed-cycle bookkeeping guard preserved**:
   - Under `Ein = Eout`, `Enext ≥ Eprev`, and `Loss ≥ 0`, positive net export is
     forbidden (`Export ≤ 0`).
5. **Quotient bridge kinematic witness**:
   - Endpoint identification gives strict shortcut in the defect-distance model.

Interpretation:
- Theorems separate **passive reconstruction no-go** from **active operator
  existence**.
- Access can be modeled through quotient/topological operations without
  invalidating the linear no-go or conservation guard.

## Runtime probes (Rust)
- `dynamic_topology_creation_probe`:
  - under-budget: `gate=false`, operational flags suppressed.
  - at-budget/over-budget: `gate=true` with toy effective-superluminal and
    pre-departure coordinate arrival flags.
- `non_conjugation_quotient_probe`:
  - strong compression exists in non-conjugation maps, but rank-preserving basis
    optimum is identity-like in the scanned basis lane; dense random scan found
    low-ratio maps.
- `topological_defect_bundle_probe`:
  - strict shortcut witness and broad random improvement fraction (~0.50064).
- `ctc_door_reinforcement_probe` + `ctc_bootstrap_fixedpoint_probe`:
  - closed-cycle conservation guard remained intact; no positive-export
    violations under theorem guard assumptions.

## Build status
- `lake build Gutoe.TopologicalAccessOperator` passed.
- `lake build Gutoe` passed (warnings only).
