# 136 — MS Cyclosporine Trial Hypothesis Pack

## Scope
Create a falsifiable trial-hypothesis artifact by combining:
1. Mechanistic efficacy lane (`cyclosporine__c20nM__buf3`)
2. PK bridge exposure window
3. Safety-gate probabilities

This is a simulation hypothesis pack, not treatment advice.

## Candidate and mechanism
- Candidate: `cyclosporine__c20nM__buf3`
- Target/off-target Ki: `2.64 / 200 nM`
- Target/off-target occupancy @20 nM: `0.8834 / 0.0909`
- Efficacy margin: `+1.4409 kJ/mol`

## Efficacy projection
2-year projection:
- ARR standard proxy: `0.2741`
- ARR combo (standard + cyclosporine): `0.2279`
- Relative ARR reduction (combo vs standard): `16.86%`

10-year projection:
- Lesion index: `0.3753 -> 0.2057` (standard -> combo)
- Disability index: `0.0387 -> 0.0214` (standard -> combo)
- Lesion reduction vs standard at 10y: `45.19%`

## PK/safety linkage
- Recommended exposure window (p25..p75): `[141.0, 264.8] ng/mL`
- Safety gate: `pass=true`
- `P(> renal caution)=0.1021`, `P(> renal high)=0.0214`, `P(> neuro caution)=0.0358`

## Power sizing (Poisson ARR approximation)
- alpha (two-sided): `0.05`
- power: `0.80`
- follow-up: `2 years`
- estimated `n per arm ≈ 921`
- estimated `n total ≈ 1842`

## Artifact paths
- `/tmp/bh_renders/ms_cyclosporine_trial_hypothesis/ms_cyclosporine_trial_hypothesis.txt`
- `/tmp/bh_renders/ms_cyclosporine_trial_hypothesis/ms_cyclosporine_trial_hypothesis.json`
- `/tmp/bh_renders/ms_cyclosporine_trial_hypothesis/ms_cyclosporine_trial_key_metrics.csv`
