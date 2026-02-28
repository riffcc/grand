# 131 — MS Macrocycle Application Simulation

## Scope
Direct reduced-order application simulation for candidate:
`macrocycle_A__c20nM__buf3`.

This lane propagates molecular mimicry drive and blocker energetics into a 10-year expected-value disease trajectory (baseline vs therapy scenarios).

## Candidate and blocker transduction
- Candidate: `macrocycle_A__c20nM__buf3`
- Concentration: `20 nM`
- Target Ki: `3 nM`
- Off-target Ki: `120 nM`
- Required shift: `1.2092 kJ/mol`
- Achieved shift: `2.6087 kJ/mol`
- Effective shift (transduction efficiency 0.30): `0.7826 kJ/mol`
- Efficacy margin: `+1.3995 kJ/mol`
- Target occupancy: `0.8696`
- Off-target occupancy: `0.1429`
- Feasible: `true`

## 10-year simulation outputs
- Baseline drive index: `0.22153`
- Standard therapy residual drive: `0.04238`
- Macrocycle-only drive: `0.04885`
- Macrocycle + standard drive: `0.00935`

Annualized relapse rate (per year):
- Baseline: `0.49612`
- Standard: `0.26521`
- Macrocycle-only: `0.27363`
- Macrocycle + standard: `0.22213`

Final lesion index (10y):
- Baseline: `1.34495`
- Standard: `0.37530`
- Macrocycle-only: `0.40791`
- Macrocycle + standard: `0.21156`

Final disability index (10y):
- Baseline: `0.13201`
- Standard: `0.03874`
- Macrocycle-only: `0.04203`
- Macrocycle + standard: `0.02202`

Derived reductions vs baseline:
- Lesion reduction (macrocycle-only): `69.67%`
- Lesion reduction (macrocycle + standard): `84.27%`
- Relapse reduction (macrocycle-only): `44.84%`
- Relapse reduction (macrocycle + standard): `55.23%`

## Artifacts
- `/tmp/bh_renders/ms_macrocycle_application/ms_macrocycle_application_report.txt`
- `/tmp/bh_renders/ms_macrocycle_application/ms_macrocycle_application_report.json`
- `/tmp/bh_renders/ms_macrocycle_application/ms_macrocycle_application_summary.csv`
- `/tmp/bh_renders/ms_macrocycle_application/ms_macrocycle_application_trajectory.csv`

## Notes
- This is a reduced-order simulation lane for comparative dynamics, not clinical guidance.
- Calibration update: lesion repair now scales with current lesion burden, preventing unrealistic full-clearance saturation under low-drive scenarios.
