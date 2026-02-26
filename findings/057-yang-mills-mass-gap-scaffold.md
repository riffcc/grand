# Finding 057 — Yang-Mills Mass Gap Scaffold (GRAND-297/298)

## Status
Scaffold implemented; finite-volume spectral-gap lane is operational.

## What was built
- New module: `crates/gutoe-physics/src/mass_gap.rs`
  - dense symmetric transfer-matrix representation
  - symmetry / nonnegativity checks
  - conservative Gershgorin lower bound on minimum eigenvalue
  - largest/second eigenvalue extraction (power + deflated power)
  - finite-volume gap observable:
    - `m_gap = -(1/a_t) ln(λ1/λ0)`
  - conservative lower-bound extraction from residual intervals:
    - with `λ0 ∈ [λ0-r0, λ0+r0]`, `λ1 ∈ [λ1-r1, λ1+r1]`
    - `m_gap ≥ -(1/a_t) ln((λ1+r1)/(λ0-r0))` when admissible
  - volume-trend helpers:
    - monotone non-increasing check with volume
    - continuum-stability interval intersection band
- New report binary: `crates/gutoe-physics/src/bin/ym_mass_gap_report.rs`
  - writes:
    - `/tmp/bh_renders/ym_mass_gap_report.txt`
    - `/tmp/bh_renders/ym_mass_gap_report.json`

## Current diagnostic output
From the synthetic finite-volume sequence:
- transfer matrices are symmetric and entrywise nonnegative
- extracted gap trend is monotone non-increasing with volume
- explicit finite-volume gap estimates are produced for each `L`
- conservative lower-bound extraction path is exercised end-to-end

## Why this matters
This closes the mechanical core for the mass-gap lane:
1. define transfer-matrix spectral observable
2. compute finite-volume gap estimates
3. attach explicit lower-bound/error machinery
4. check finite-volume trend behavior prior to continuum claims

## Remaining closure steps
- replace synthetic toy matrices with SU(3)/Z₃-orbit-derived correlator or transfer data
- pin physically justified uncertainty model (stat + systematic)
- produce continuum-stability result on real generated ensembles
