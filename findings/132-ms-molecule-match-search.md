# 132 — MS Molecule Match Search

## Scope
Find the closest molecule proxy match to the validated MS application profile from:
`macrocycle_A__c20nM__buf3`.

The match objective uses the same reduced-order lane and minimizes profile distance over:
- target occupancy
- off-target occupancy
- efficacy margin
- combo lesion reduction
- combo relapse reduction
- combo disability endpoint

## Reference profile
Derived from `macrocycle_A__c20nM__buf3` in current calibration:
- target occupancy: `0.869565217`
- off-target occupancy: `0.142857143`
- efficacy margin: `1.399463709 kJ/mol`
- combo lesion reduction: `0.842696763`
- combo relapse reduction: `0.552270676`
- combo disability index: `0.022023777`

## Search result
- Feasible screened candidates: `171`
- Best match: `macrocycle_A__c20nM__buf3`
  - profile distance: `0.000000000`
  - profile match score: `1.000000000`
- Best alternative (non-identical): `macrocycle_A__c20nM__buf4`
  - profile distance: `0.166666667`
  - profile match score: `0.857142857`

## Top non-identical scaffold alternatives
1. `macrocycle_C__c30nM__buf5`
   - profile distance: `0.251521727`
   - target/off-target occupancy: `0.882352941 / 0.142857143`
   - efficacy margin: `1.414297469 kJ/mol`
   - combo lesion/relapse reduction: `0.856877012 / 0.562557877`
2. `macrocycle_D__c20nM__buf3`
   - profile distance: `0.277481377`
   - target/off-target occupancy: `0.884955752 / 0.160000000`
   - efficacy margin: `1.357139738 kJ/mol`
   - combo lesion/relapse reduction: `0.839454948 / 0.549924009`
3. `spiro_B__c35nM__buf5`
   - profile distance: `0.352735383`
   - target/off-target occupancy: `0.892857143 / 0.127272727`
   - efficacy margin: `1.358625200 kJ/mol`
   - combo lesion/relapse reduction: `0.853645112 / 0.560210038`

## Artifacts
- `/tmp/bh_renders/ms_molecule_match_search/ms_molecule_match_search.txt`
- `/tmp/bh_renders/ms_molecule_match_search/ms_molecule_match_search.csv`
- `/tmp/bh_renders/ms_molecule_match_search/ms_molecule_match_search.json`

## Notes
- This lane is in-silico ranking only, not clinical guidance.
- Validation command used: `cargo run -q -p gutoe-physics --bin ms_molecule_match_search`.
