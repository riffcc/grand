# 235 — Alpha Correction Priority: Structural Coefficient Probe

## Scope
Priority lane requested: structural justification for
- `alpha^-1 = 137 + 5alpha - 9alpha^2`

and linked mass-ratio correction lane
- `mp/me = 6pi^5 + 5alpha - 36alpha^2`

Runner:
- `crates/gutoe-physics/src/bin/ctc_alpha_correction_priority_probe.rs`

Artifacts:
- `/tmp/bh_renders/ctc_alpha_correction_priority_probe/ctc_alpha_correction_priority_probe.txt`
- `/tmp/bh_renders/ctc_alpha_correction_priority_probe/ctc_alpha_correction_priority_probe.json`

## Results
With fixed linear coefficient `a1=5` (grade-level count):
- implied exact quadratic coefficient from CODATA:
  - `b_alpha_exact = 9.156308355`
- implied exact quadratic coefficient from `mp/me`:
  - `b_mp_exact = 36.093814512`
- ratio:
  - `b_mp_exact / b_alpha_exact = 3.941913136` (near 4)

Structural candidates tested:
- `9 = z3^2 = 3^2`
- `9 = grade2 + z3 = 6 + 3`
- `36 = grade1 * z3^2 = 4*9`
- `36 = grade2^2 = 6^2`

Numerical lane checks:
- `137 + 5alpha - 9alpha^2 = 137.036007500632` (0.061 ppm)
- `6pi^5 + 5alpha - 36alpha^2 = 1836.152678425750` (0.003 ppm)

## Interpretation (strict)
- This probe gives a concrete structural candidate family for coefficients `5, 9, 36`.
- It does **not** yet prove uniqueness of the `9` decomposition.
- It materially narrows the correction hunt from unconstrained fitting to integer-structured candidates tied to existing Cl(1,3)/Z3 counts.
