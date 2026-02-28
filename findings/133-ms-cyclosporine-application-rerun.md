# 133 — MS Application Rerun with Cyclosporine Parameters

## Scope
Rerun the 10-year MS application simulation using real-world cyclosporine Ki profile
against the same lane assumptions used for `macrocycle_A__c20nM__buf3`.

## Inputs
Cyclosporine run:
- `GUTOE_MS_CANDIDATE_LABEL=cyclosporine__c20nM__buf3`
- `GUTOE_MS_CANDIDATE_CONC_NM=20`
- `GUTOE_MS_CANDIDATE_TARGET_KI_NM=2.64`
- `GUTOE_MS_CANDIDATE_OFFTARGET_KI_NM=200`
- `GUTOE_MS_CANDIDATE_MAX_SHIFT_KJ_MOL=3.0`
- `GUTOE_MS_CANDIDATE_SAFETY_BUFFER_KJ_MOL=0.3`

Baseline comparison run:
- `macrocycle_A__c20nM__buf3` (default candidate settings)

## Cyclosporine lane output
Blocker-level:
- target occupancy: `0.883392`
- off-target occupancy: `0.090909`
- achieved shift: `2.650177 kJ/mol`
- effective shift (efficiency 0.30): `0.795053 kJ/mol`
- efficacy margin: `+1.440945 kJ/mol`
- feasible: `true`

10-year outcomes:
- lesion reduction (candidate-only): `72.02%`
- lesion reduction (candidate + standard): `84.71%`
- relapse reduction (candidate-only): `46.49%`
- relapse reduction (candidate + standard): `55.54%`
- final disability index (candidate + standard): `0.021419`

## Delta vs macrocycle_A__c20nM__buf3
- target occupancy: `+0.013827`
- off-target occupancy: `-0.051948`
- efficacy margin: `+0.041481 kJ/mol`
- lesion reduction (candidate-only): `+2.35 pp`
- lesion reduction (candidate + standard): `+0.44 pp`
- relapse reduction (candidate-only): `+1.65 pp`
- relapse reduction (candidate + standard): `+0.32 pp`
- final disability index (candidate + standard): `-0.000604`

## Artifacts
Baseline:
- `/tmp/bh_renders/ms_macrocycle_application_baseline/ms_macrocycle_application_report.json`

Cyclosporine rerun:
- `/tmp/bh_renders/ms_macrocycle_application_cyclosporine/ms_macrocycle_application_report.txt`
- `/tmp/bh_renders/ms_macrocycle_application_cyclosporine/ms_macrocycle_application_report.json`
- `/tmp/bh_renders/ms_macrocycle_application_cyclosporine/ms_macrocycle_application_summary.csv`
- `/tmp/bh_renders/ms_macrocycle_application_cyclosporine/ms_macrocycle_application_trajectory.csv`

## Notes
- This is a reduced-order mechanistic simulation comparison, not clinical guidance.
- Ki mapping to this lane remains an in-silico transduction step, not a direct trial-outcome predictor.
