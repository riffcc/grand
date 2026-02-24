# 003 — Koide / Wyler Mass Prediction Chain

## Scope

This note documents the existing proof+implementation chain for the lepton mass structure and proton/electron scale relation.

## Structural pieces already in repo

- Lean formal chain:
  - `lean/Gutoe/KoideMasses.lean`
  - `lean/Gutoe/LeptonMass.lean`
  - `lean/Gutoe/MassSpectrum.lean`
- Rust runtime chain:
  - `crates/gutoe-em/src/alpha.rs`
  - key functions: `koide_ratio`, `z3_harmonic_masses`

## Equations tracked

1. **Koide structural target**

`K = (Σ m_i) / (Σ sqrt(m_i))^2 = 2/3`

In this codebase, this is tied to Clifford grade counting:

`grade1_4d / grade2_4d = 4/6 = 2/3`

2. **Z3 harmonic mass family**

`m_k = M^2 (1 + s cos(δ + 2πk/3))^2`

with exact identity:

`Koide = (1 + s^2/2) / 3`

hence `Koide = 2/3 <-> s^2 = 2`.

3. **Proton/electron geometric scale**

`m_p / m_e ~ 6π^5`

Lean comments in `MassSpectrum.lean` also track the corrected form:

`m_p / m_e = 6π^5 + n_grades / alpha_inv = 6π^5 + 5/137`.

## Why this matters

- Koide structure is not treated as a free fit in this project; it is constrained by the Z3 harmonic/circulant construction.
- The `6π^5` scale plus `n_grades / alpha_inv` correction ties the baryon/lepton scale bridge to shared primitives already proven in Lean.

## Practical parity check points

- Keep Rust comments/tests in `alpha.rs` aligned with:
  - `koide_target = 2/3`
  - `s^2 = 2` at target
  - `m_p/m_e` correction term `+ 5/137` when quoted.
