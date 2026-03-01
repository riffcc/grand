# Finding 257 — Absolute-Scale Endgame Closure (Lattice Branch)

Date: 2026-03-01  
Runner: `cargo run -q -p gutoe-em --bin yukawa_absolute_scale_endgame_report`

Artifacts:
- `/tmp/bh_renders/yukawa_absolute_scale_endgame_report.txt`
- `/tmp/bh_renders/yukawa_absolute_scale_endgame_report.json`

## Summary

Executed the second remaining lane: absolute scale closure from the proton/lattice branch against Fermi/PDG anchors.

All report gates were set at `1%` relative tolerance.

## Results

- Electron anchor (`mp / 1836`):
  - predicted `0.5110414194 MeV`
  - reference `0.5109989500 MeV`
  - relative error `+8.311e-5` (`+0.00831%`) — pass

- VEV closure (lattice vs Fermi):
  - `v_lattice = 245.2998813 GeV`
  - `v_fermi = 246.2196508 GeV`
  - relative error `-3.736e-3` (`-0.3736%`) — pass

- Lattice-branch masses:
  - `mW = 80.01348 GeV` vs `80.377` (rel `-0.452%`) — pass
  - `mZ = 91.22941 GeV` vs `91.1876` (rel `+0.0458%`) — pass
  - `mH = 125.07889 GeV` vs `125.25` (rel `-0.1366%`) — pass

- Fermi-branch masses:
  - `mW` rel `-0.0790%` — pass
  - `mZ` rel `+0.4210%` — pass
  - `mH` rel `+0.2378%` — pass

## Gate Verdict

`overall_pass = true`

This closes the absolute-scale endgame lane at the configured 1% envelope for:
- electron anchor,
- VEV mapping,
- and W/Z/H mass-sector outputs in both lattice and Fermi branches.

