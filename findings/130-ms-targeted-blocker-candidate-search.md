# MS Targeted Blocker Candidate Search — First Ranked Shortlist

Date: 2026-02-28  
Status: Implemented (search binary + ranked candidates)

## Goal
Convert `blocker_feasible=true` into concrete candidate parameter sets by scanning:
- concentration
- target Ki
- off-target Ki
- max interface energy-shift capacity
- safety buffer

and ranking feasible candidates by efficacy margin + selectivity + exposure penalty.

## What shipped

- Search logic in:
  - `crates/gutoe-physics/src/ms_autoimmunity.rs`
    - `TargetedBlockerCandidateInput`
    - `TargetedBlockerCandidateScore`
    - `evaluate_targeted_blocker_candidate`
- New binary:
  - `crates/gutoe-physics/src/bin/ms_targeted_blocker_search.rs`

## Executed run

Command:
- `cargo run -q -p gutoe-physics --bin ms_targeted_blocker_search`

Summary:
- `screened_candidates = 90`
- `feasible_candidates = 43`
- `top_candidate = macrocycle_A__c20nM__buf3`

## Top candidate

`macrocycle_A__c20nM__buf3`

- concentration: `20 nM`
- target Ki: `3.0 nM`
- off-target Ki: `120 nM`
- max energy shift: `3.0 kJ/mol`
- required shift: `1.209 kJ/mol`
- achieved shift: `2.609 kJ/mol`
- efficacy margin: `+1.399 kJ/mol`
- target occupancy: `0.870`
- off-target occupancy: `0.143`
- selectivity ratio: `6.09`

## Top-3 snapshot

1. `macrocycle_A__c20nM__buf3`
   - score `2.398`
2. `macrocycle_A__c10nM__buf3`
   - score `2.366`
3. `macrocycle_A__c30nM__buf3`
   - score `2.343`

## Artifacts

- `/tmp/bh_renders/ms_targeted_blocker_search/ms_targeted_blocker_search.txt`
- `/tmp/bh_renders/ms_targeted_blocker_search/ms_targeted_blocker_search.csv`
- `/tmp/bh_renders/ms_targeted_blocker_search/ms_targeted_blocker_search.json`

## Honesty statement

These are computational candidate profiles in a reduced-order mechanistic lane,
not validated molecules and not clinical recommendations. They are intended as
next-step in-silico design priors.
