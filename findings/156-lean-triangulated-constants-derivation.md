# GRAND — Lean Derivation: Triangulated Constants

Date: 2026-02-28

## Added

- `lean/Gutoe/TriangulatedConstants.lean`
- `lean/lakefile.lean` root registration: `Gutoe.TriangulatedConstants`

## Formal Content

Defined and proved Cl(1,3)-count candidate forms:

1. `pCandidateQ`
   - `p = α⁻¹/(|grade₁|+|grade₂|) - 1/((|grade₂|+1)*N_gauge)`
   - closed form: `137/10 - 1/(7*12)`

2. `kappaCandidateQ`
   - `κ = (60/11) * (19/3 + 1/36 + 1/(7*13*136))`
   - all factors tied to shared Cl counts (`7, 12, 13, 19, 136`, `60/11`)

3. `ewCoeffCandidateQ`
   - `c = d/2 + |grade₂|/(d-|SU(2)|) - 1/((|grade₂|+1)T(16))`
   - closed form: `8 + 6/13 - 1/(7*136)`

Also proved strict rational proximity to frozen runtime anchors:

- `|pCandidateQ - pFrozenQ| < 1/50000`
- `|kappaCandidateQ - kappaFrozenQ| < 1/50000`
- `|ewCoeffCandidateQ - ewCoeffFrozenQ| < 1/1000000`

## Verification

- `cd lean && lake build Gutoe`
- Result: success (`8129 jobs`, `Gutoe.TriangulatedConstants` built)

## Honesty

These are formal Cl(1,3)-count candidate derivations plus proved proximity to frozen triangulation constants.
They are now theorem-backed and CI-build-verified.
