# 203 — Bootstrap fixed-point + temporal-loan stability criterion

## What was tested
User hypothesis:
- "Send energy back to create the same door" (self-consistent bootstrap loop),
- then compound capability by investing early and producing more future energy.

## Formalized (Lean)
New module:
- `Gutoe.CTCBootstrapFixedPoint`

Core map:
- `E_past = η * E_future - loss`
- fixed-point condition `E* = ηE* - loss`

Key theorem:
- under `η ≤ 1`, `loss ≥ 0`, and `E* > 0`,
  positive fixed point implies exactly:
  - `η = 1` and `loss = 0`.

Interpretation:
- Closed-cycle positive bootstrap is only possible in the ideal lossless limit.
- Any nonzero loss with non-amplifying return kills positive closed fixed points.

Build:
- `lake build Gutoe.CTCBootstrapFixedPoint` passes.
- full `lake build Gutoe` passes.

## Numerical probe (Rust)
New bin:
- `ctc_bootstrap_fixedpoint_probe`

Outputs:
- `/tmp/bh_renders/ctc_bootstrap_fixedpoint_probe/ctc_bootstrap_fixedpoint_probe.txt`
- `/tmp/bh_renders/ctc_bootstrap_fixedpoint_probe/ctc_bootstrap_fixedpoint_probe.json`

Default (quark-scale threshold, near-unity eta, nonzero loss, zero inflow):
- closed fixed point: negative (not feasible)
- threshold bootstrap: not feasible

Ideal closed case (`η=1`, `loss=0`, `inflow=0`):
- continuum of fixed points (self-consistent circulation class)
- no net export implied (matches finding 202 constraints)

## Temporal-loan (compounding) criterion
With capability feedback reduced to an affine recurrence:
- `L_{k+1} = A + B * L_k`
- `B` is effective compounding gain (`B ~ η * sendback_fraction * capability_return`).

Regimes:
- `B < 1`: convergent finite temporal loan.
- `B = 1`: marginal; requires ideal lossless closure to avoid collapse.
- `B > 1`: runaway growth unless bounded by hard caps (resource/flux/topology limits).

So yes, "loaning from the future" can compound, but only within this stability structure.
