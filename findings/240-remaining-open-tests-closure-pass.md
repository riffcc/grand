# 240 — Remaining-Open Tests Closure Pass (T05/T08/T09/T10/T12/T14/T15/T16/T17/T18/T19/T20)

## Scope
Executed a dedicated closure pass for the previously open set:
- `T05, T08, T09, T10, T12, T14, T15, T16, T17, T18, T19, T20`

Runner:
- `crates/gutoe-physics/src/bin/ctc_remaining12_runner.rs`

Independent replication step for `T16`:
- reran `ctc_public_door_joint_fit` with
  `GUTOE_DOOR_JOINT_OUT=/tmp/bh_renders/ctc_public_door_joint_fit_rerun`

Artifacts:
- `/tmp/bh_renders/ctc_remaining12_runner/ctc_remaining12_runner.txt`
- `/tmp/bh_renders/ctc_remaining12_runner/ctc_remaining12_runner.json`
- `/tmp/bh_renders/ctc_public_door_joint_fit_rerun/ctc_public_door_joint_fit.json`

## Results
- `T05` OPEN (blocked): isotropic CB lane not quantified in current joint artifact.
- `T08` OPEN (blocked): no quantified void+lensing scalar pair.
- `T09` FAIL: EHT coarse residual lane gives `z=0` (<3 sigma threshold).
- `T10` OPEN (blocked): no separate M87*/SgrA* parameter rows in coarse artifact.
- `T12` PASS: `max_delta_pred_sigma_LOO=0.000000` (usable folds: 8).
- `T14` PASS: slope lock `z=0.000982`.
- `T15` PASS: topology gain lock `|G-1|=0` from topology overdetermination artifact.
- `T16` PASS: independent rerun reproduces null (`door_detected=false` in both runs).
- `T17` OPEN (blocked): no photonic Pi-target lab artifact in current runset.
- `T18` OPEN (blocked): no superconducting lock-threshold artifact in current runset.
- `T19` OPEN (blocked): no RF-PLL causal-order artifact in current runset.
- `T20` OPEN (blocked): no quantified cross-lane posterior overlap artifact.

## Battery impact
After importing this runner into the main harness:
- `PASS=11`, `FAIL=4`, `OPEN=7`, `TOTAL=22`

Main artifact:
- `/tmp/bh_renders/ctc_falsification_20/ctc_falsification_20.txt`

## Clarification
`T06` is the anisotropic cosmic birefringence amplitude lane (`A_CB/sigma`), not the weak-angle slope-lock lane.
Weak-angle slope-lock is `T14`, and it now passes with high margin.
