# Finding 256 — Neutrino Endgame Lane Split (No-Fit vs Triangulated)

Date: 2026-03-01  
Runner: `cargo run -q -p gutoe-em --bin yukawa_neutrino_endgame_report`

Artifacts:
- `/tmp/bh_renders/yukawa_neutrino_endgame_report.txt`
- `/tmp/bh_renders/yukawa_neutrino_endgame_report.json`

## Summary

We executed the remaining neutrino lane as two explicit tracks:

1. **Structural no-fit lane** (`neutrino_absolute_masses_from_texture`)
2. **Oscillation-triangulated lane** (`triangulate_neutrino_from_splittings`)

Both tracks keep the same texture ordering and hierarchy prediction; the split isolates where closure is structural vs forced.

## Core Results

### Structural no-fit lane
- Hierarchy: `normal` (pass)
- Tiny-mass bounds: `m3=7.9045e-3 eV`, `sum=9.2666e-3 eV` (pass)
- Splitting ratio: `Δm32/Δm21 = 32.6790` vs target `32.5764` (relative `+0.315%`, pass)
- Absolute splittings:
  - `Δm21 = 1.8552e-6` vs `7.53e-5` (relative `-97.536%`, fail)
  - `Δm32 = 6.0626e-5` vs `2.453e-3` (relative `-97.528%`, fail)
- Koide diagnostic:
  - `Kν = 0.5854167427`
  - `s²ν = 1.5125004564`
  - vs `Kν=1/2`: relative `+17.083%`

### Oscillation-triangulated lane
- `p_triangulated = 13.6881104338`
- `kappa_geo = 34.6973960555`
- Ratio closure: relative error `1.77e-11` (pass)
- Absolute splitting closure:
  - `Δm21` relative error `-8.84e-12`
  - `Δm32` relative error `+8.84e-12`
  - (both pass by construction)
- Koide diagnostic remains near structural lane:
  - `Kν = 0.58528601995`
  - `s²ν = 1.5117161197`
  - vs `Kν=1/2`: relative `+17.057%`

## Interpretation

- The **structural no-fit lane** correctly captures hierarchy, tinyness, and splitting **ratio**, but not absolute splitting scale.
- The **triangulated lane** closes oscillation targets numerically and exposes the needed normalization scale (`kappa_geo`), but this is a forced diagnostic lane, not a zero-parameter proof.
- Current in-engine data does **not** support `Kν = 1/2` in either lane.

