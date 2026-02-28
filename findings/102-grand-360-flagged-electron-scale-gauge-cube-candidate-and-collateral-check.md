# GRAND-360 — Flagged Electron Scale Candidate (`λ_QG^-3 = 12^3`) + Collateral Check

## Candidate

Flagged transduction candidate from the constrained corrected-ratio sweep:

- `m_e = M_Pl * F_flagged`
- `F_flagged = α^13 * (115/22)^3 * (67/66) * λ_QG^-3`

Using `λ_QG = 1/12`, this is algebraically equivalent to:

- `F_flagged = α^13 * (115/22)^3 * (67/66) * 12^3`

Numerical fit quality (from sweep):

- `m_e,pred = 0.510751 MeV`
- relative error `-4.852e-4` (`0.0485%`)

This is below the `<0.1%` numerical threshold, but remains **flagged** pending
mechanistic justification for inverse-cubic `λ_QG` dependence.

## Lean status

Added theorem module:

- `lean/Gutoe/ElectronScaleTransduction.lean`

Key proved statements:

- `alphaStructural = 1/137`
- corrected ratio lane term is `115/22`
- `lambda_qg^(-3) = (12:ℝ)^3`
- flagged and gauge-cube forms are equal:
  `electronScaleFactorFlagged = electronScaleFactorGaugeCube`

Build:

- `lake build Gutoe` passes.

## Collateral checks (no regressions observed)

Checked existing validated lanes after adding the candidate theorem/scaffold.

1. Proton lane

- `cargo run --release -p gutoe-physics --bin proton_mass_report`
- unchanged:
  - `mp/me = 1836`
  - `m_p,pred = 938.194072 MeV`
  - relative error `-8.315e-5`

2. CMB + sigma8 lane

- `GUTOE_CLASS_BIN=/tmp/class_public/class cargo run --release -p gutoe-physics --bin cmb_full_derived_report`
- unchanged:
  - `TT = 1.263`
  - `TE = 1.109`
  - `EE = 1.064`
  - `sigma8 = 0.811221`

## Honest boundary

- The candidate is **numerically strong** and now **formalized as a theorem lane
  scaffold**.
- It is not yet accepted as final derivation because inverse `λ_QG` powers may
  be amplifying a small UV parameter without a completed physical bridge.
- Next acceptance gate: derive (or falsify) why electron scale transduction
  should include exactly cubic inverse `λ_QG` dependence.
