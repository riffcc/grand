# GRAND-359 — Electron Absolute-Scale Sweep (Hard Negative + Constrained Survivors)

## Target

Derive absolute electron scale from Planck scale via short structural maps:

- `m_e = M_Pl * F(shared dimensionless primitives)`
- target ratio: `m_e / M_Pl = 4.185462873156e-23`

Kill criteria used:

- operation budget: `<= 8`
- `|rel_err| > 5%` => kill
- `|rel_err| < 1%` => promote investigation
- `|rel_err| < 0.1%` => promote Lean

## Hard negative (critical)

Previous best candidate (0.254% error):

- `F = α^11 * (60/11)^2 * (5/11) * (66/67)`

is structurally invalid because it relies on superseded ratio `60/11`.

One-point substitution with corrected Lean-verified ratio:

- replace `(60/11)` with `(115/22)`
- prediction shifts from `0.509699 MeV` to `0.468109 MeV`
- error degrades from `0.254%` to `8.393%`

Conclusion: this candidate is killed as a false positive.

## Sweep implementation

Binary:

- `crates/gutoe-physics/src/bin/electron_scale_sweep.rs`

Artifacts:

- `/tmp/bh_renders/electron_scale_sweep/electron_scale_sweep.txt`
- `/tmp/bh_renders/electron_scale_sweep/electron_scale_sweep.json`

## Constrained results

### A) Short-form subset (`<= 6` ops)

Best survivor:

- `F = α^10 * (67/66)^2 * (5/11)^3`
- `m_e,pred = 0.507300 MeV`
- error `-0.724%` (investigation-grade)

No `<0.1%` candidate in this subset.

### B) Corrected-ratio product family (forcing `R = 115/22` participation)

Allowing signed exponents produced:

- `F = α^13 * R^3 * (67/66)^1 * λ_QG^-3`
- `m_e,pred = 0.510751 MeV`
- error `-0.0485%` (<0.1%)

But this uses strong inverse-`λ_QG` amplification; status is **structurally
suspicious** until justified by an independent theorem.

### C) Corrected-ratio products with nonnegative exponents only

Best survivor:

- `F = α^11 * R^3 * (5/11)^3`
- `m_e,pred = 0.513225 MeV`
- error `+0.4356%` (investigation-grade)

No `<0.1%` candidate under nonnegative-exponent discipline.

## Honest interpretation

- Confirmed hard negative: naive short multiplicative map using superseded `60/11`
  does not survive consistency.
- Viable investigation-grade survivors exist (~0.4–0.7%).
- A numerically excellent candidate exists (<0.1%) but currently depends on
  inverse `λ_QG` powers and is not yet trusted as structural physics.

## Status

- Finding type: **hard negative + constrained follow-up**
- Mountain remains open.
- Next closure step: derive or reject inverse-`λ_QG` participation by theorem
  before promoting any `<0.1%` candidate.
