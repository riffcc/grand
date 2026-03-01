# 239 — Synthesis: Unified Perturbative Series for `alpha^-1` and `mp/me`

## Scope
Consolidate the current closure state after:
- weak-angle identity/running fits,
- strict falsification harness execution,
- alpha/mp-me correction hunts,
- locked-candidate and dual-lane checks.

This finding records the current best structural series and exact residuals from run artifacts.

## Primary equations (current candidate)
Let
- `b = 9 + 5/32`,
- `c = 1/125`,
- `alpha` defined self-consistently by:

`alpha^-1 = 137 + 5 alpha - b alpha^2 - c alpha^3`.

Then use the linked mass-ratio lane (with one-loop term unscaled):

`mp/me = 6 pi^5 + 5 alpha - g b alpha^2 - g c alpha^3`,

with

`g = 4 - 8 alpha`.

## Fixed-point character (critical)
The `alpha` lane is an implicit fixed-point equation, not a direct plug-in formula:

`x = 137 + 5/x - b/x^2 - c/x^3`, where `x = alpha^-1`.

So `alpha` is selected by self-consistency of the structural recursion itself.
In current numeric evaluation, this equation has a unique physically relevant
positive root near `137.036`, and that root gives the reported sub-ppb closure.

## Structural reading of coefficients
- `137`: leading-order Clifford lane (existing Lean-proven lane for LO alpha identity)
- `5`: grade-level count (grades 0..4)
- `9`: `3^2` from Z3 order
- `5/32`: grade-level over power-of-two by grade-level (`5/2^5`)
- `1/125`: cubic grade-level inverse power (`1/5^3`)
- `4`: grade-1 spacetime-vector count
- `8`: SU(3) adjoint dimension (color correction in multiplier)

## Numerical closures (from artifacts)
### Alpha lane
Using `b = 9 + 5/32`, `c = 1/125`, solve the implicit equation for `alpha^-1`:
- `alpha_inv_pred = 137.035999176998740...`
- CODATA 2022: `137.035999177`
- residual:
  - `9.195045814584537e-9 ppm`
  - `9.195045814584537e-6 ppb`

### mp/me lane
Using the same `alpha`, same `b,c`, and `g = 4 - 8 alpha`:
- `mp/me_pred = 1836.152673595903...`
- CODATA 2022: `1836.15267343`
- residual:
  - `0.000090354... ppm`
  - `0.090354... ppb`

## Falsification harness update
`ctc_falsification_20_harness` now includes two corrected precision lanes:
- `T21`: alpha structural cubic closure
- `T22`: mp/me shared-series closure

Current strict harness summary:
- `PASS=7`, `FAIL=3`, `OPEN=12`, `TOTAL=22`

Artifacts:
- `/tmp/bh_renders/ctc_falsification_20/ctc_falsification_20.txt`
- `/tmp/bh_renders/ctc_falsification_20/ctc_falsification_20.json`

## Relationship to previous lanes
- Weak-angle lane remains independently strong (9-point zero-free run retained).
- The alpha/mp-me lane now forms a shared perturbative family rather than two disconnected numerology fits.
- The strict no-door public-data lanes remain negative/coarse-null at current data quality in this harness.

## Open questions (priority)
1. **Why `6 pi^5` as the mp/me leading geometric scale?**
   Current status: phenomenologically excellent, structural derivation still pending.

2. **Why the cubic sign and magnitude (`-1/5^3`) exactly?**
   Current status: numerically near-exact closure; derivation from diagram/topology counting not yet formalized.

3. **Do we need `alpha^4` terms, or does this lane effectively terminate at cubic for current precision targets?**
   Current status: cubic already reaches sub-ppb on both lanes; next-order necessity unresolved.

## Immediate next execution recommendation
- Run the coefficient-uniqueness sweep constrained to small Cl(1,3)/Z3 expressions and report ranked candidates for `b` and `g` with complexity penalties.
- In parallel, add a symbolic derivation lane candidate for the `-1/125` sign from grade-level recursion constraints.
