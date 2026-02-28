# 135 — MS Cyclosporine Safety Gate

## Scope
Evaluate PK-bridge exposure distribution against configurable simulation safety windows.

This is a model gate for translational triage, not clinical decision support.

## Gate configuration (ng/mL)
- Target zone: `[80, 300]`
- Renal caution threshold: `350`
- Renal high threshold: `500`
- Neuro caution threshold: `450`

Probability thresholds:
- `P(> renal caution) <= 0.15`
- `P(> renal high) <= 0.05`
- `P(> neuro caution) <= 0.08`
- `P(in target zone) >= 0.50`

## Results
- `overall_pass = true`
- `P(in target zone) = 0.79608`
- `P(> renal caution) = 0.10206`
- `P(> renal high) = 0.02136`
- `P(> neuro caution) = 0.03576`
- PK center: `p50 = 193.23 ng/mL`, `p95 = 417.26 ng/mL`

## Artifact paths
- `/tmp/bh_renders/ms_cyclosporine_safety_gate/ms_cyclosporine_safety_gate.txt`
- `/tmp/bh_renders/ms_cyclosporine_safety_gate/ms_cyclosporine_safety_gate.json`
