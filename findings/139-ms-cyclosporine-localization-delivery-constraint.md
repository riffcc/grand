# 139 — MS Cyclosporine: Delivery-Constraint Framing (Localization Follow-up)

## Why this note exists
The historical signal and the current simulation now align on one point:
- binding/immune-mechanism signal exists,
- but broad systemic exposure is the limiting factor.

This lane reframes the question from "does the molecule bind?" to
"can effect-site efficacy be sustained while systemic toxicity risk stays bounded?"

## Re-run verification
Executed:
- `cargo run -q -p gutoe-physics --bin ms_localized_dual_compartment`
- `cargo run -q -p gutoe-physics --bin ms_tolerance_induction_dynamics`
- `cargo run -q -p gutoe-physics --bin ms_boundary_shift_combined_sweep`
- `cargo test -q -p gutoe-physics --lib ms_boundary_shift`

Artifacts:
- `/tmp/bh_renders/ms_localized_dual_compartment/ms_localized_dual_compartment.json`
- `/tmp/bh_renders/ms_tolerance_induction_dynamics/ms_tolerance_induction_dynamics.json`
- `/tmp/bh_renders/ms_boundary_shift_combined_sweep/ms_boundary_shift_combined_sweep.json`

## Core readout

### 1) Split-gate localization model passes
Default conceptual setting (`localization_factor=0.60`, `transduction_efficiency=0.30`):
- `overall_pass = true`
- efficacy gate:
  - ARR reduction vs standard (2y): `0.1702`
  - lesion reduction vs standard (10y): `0.4562`
- systemic gate:
  - `p50=115.94 ng/mL`, `p95=250.35 ng/mL`
  - `P(>renal caution)=0.00930`
  - `P(>renal high)=0.00084`
  - `P(>neuro caution)=0.00178`

Interpretation: once efficacy-site and systemic gates are treated explicitly, a plausible
in-model region exists where both pass.

### 2) Dynamic tolerance shift improves over no-tolerance control
Under the same candidate assumptions:
- ARR reduction vs no-tolerance: `0.0477`
- lesion reduction vs no-tolerance: `0.1926`
- disability index: `0.0214 -> 0.0173`

Interpretation: boundary-shift/tolerance terms reduce required suppression pressure and
improve outcomes in the reduced-order dynamics.

### 3) Combined sweep finds stable pass regions
Top-ranked pass row from combined sweep:
- `transduction_efficiency=0.15`
- `localization_factor=0.50`
- `tolerance_shift=0.60 kJ/mol`
- ARR reduction (2y): `0.2088`
- lesion reduction (10y): `0.5576`
- `P(>renal high)=0.00020`
- `overall_pass=true`

Interpretation: in this scaffold, localization plus boundary-shift terms can recover both
safety and efficacy without requiring high systemic exposure.

## What this does and does not claim
- Claims supported: model-level evidence that delivery/localization constraints are the
  dominant limiter after target-mechanism matching.
- Not claimed: clinical dosing guidance or patient-level efficacy/safety predictions.

This remains a translational simulation lane and should be treated as hypothesis generation
for subsequent PK/PD and clinical-design work.
