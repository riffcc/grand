# 015 — M87* Campaign Parameter Lock (EHT 2017 ↔ GUTOE Mapping)

Status: locked runtime defaults for synthetic-observation workflows.

## Locked Observer/Source Parameters

- Source: M87*
- Mass: `M = 6.5e9 M_sun`
- Distance: `D = 16.8 Mpc`
- Inclination: `theta_obs = 17 deg`
- Position angle: `PA = 288 deg` (stored as camera azimuth override in renderer workflows)
- Frequency anchor: `nu_obs = 230 GHz`

## GUTOE Runtime Mapping

- Schwarzschild radius:
  - `r_s = 2GM/c^2`
  - for locked `(M, D)`, angular scale is `r_s ~= 7.638217 uas`
- Camera FOV used by the EHT-aligned preset:
  - `fov_rs = 8.5`
  - `fov_uas = fov_rs * r_s_uas ~= 64.924849 uas`
- Pixel angular scale:
  - `pixel_uas = fov_uas / width`
  - at `512x512`: `pixel_uas ~= 0.126806 uas`
- GR shadow diameter reference:
  - `theta_shadow,GR = 2*sqrt(27)*GM/(c^2 D) ~= 39.689342 uas`

These values are emitted in `eht_uv` summary JSON for auditability.

## Code Anchors

- EHT-aligned preset view:
  - `crates/gutoe-gpu/src/bin/bh_render.rs` (`m87_eht2017`)
- Angular audit and parameter export:
  - `crates/gutoe-gpu/src/bin/bh_render.rs` (`run_eht_uv_export`)

## Primary Reference Set

- EHT Collaboration (2019), M87* first image paper set (Papers I–VI), source geometry and campaign conventions.
- EHT Collaboration parameter conventions used in M87* imaging/model libraries:
  mass, distance, inclination, position angle, and 230 GHz observing band.
