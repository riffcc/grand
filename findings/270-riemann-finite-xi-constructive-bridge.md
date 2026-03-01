# Finding 270 — RH Finite-Xi Constructive Bridge in Lean

Date: 2026-03-01  
Scope: Add a fully constructive Lean case where `Xi ↔ Spec` is proved exactly (not assumed), then derive RH for that model.

## Added

New module:

- `lean/Gutoe/RiemannFiniteXiModel.lean`

Wired into roots:

- `lean/lakefile.lean` includes `Gutoe.RiemannFiniteXiModel`.

## Constructive objects

- `finiteSpecSet : Finset ℝ → Set ℝ`
- `XiFinite : Finset ℝ → (ℂ → ℂ)` with
  - `XiFinite spec s = ∏ t∈spec (s - (1/2 + i t))`

## Proven theorems

- `XiFinite_zero_of_mem`
  - every listed spectral ordinate gives a zero.
- `XiFinite_zero_iff_exists`
  - exact finite zero characterization.
- `finiteXi_spectralBridge`
  - exact `SpectralBridge` for `XiFinite`.
- `rh_XiFinite`
  - RH-for-`XiFinite` follows by exact bridge.

## Why this matters

This is the first RH-lane module where bridge is established constructively from an explicit `Xi` definition, not just carried as a hypothesis.

It does **not** prove RH for the real completed zeta function, but it converts part of the bridge program into executable theorem content:

- “if the function is built from spectral factors, bridge and RH follow” is now formal and compiled.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe
```

Result: **passes** (`8167` jobs, warnings only).

