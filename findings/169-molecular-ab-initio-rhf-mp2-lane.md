# 169 — Molecular Ab-initio Lane (RHF + MP2)

Date: 2026-02-28

## Scope
Add a molecular quantum-chemistry lane beyond atomic SCF:
- explicit Gaussian AO integrals,
- closed-shell RHF self-consistent field,
- post-HF MP2 correlation correction,
- benchmark report over a small molecule panel.

## Code landed
- New module:
  - `crates/gutoe-physics/src/molecular_ab_initio.rs`
- New report binary:
  - `crates/gutoe-physics/src/bin/molecular_ab_initio_report.rs`
- Library exports updated:
  - `crates/gutoe-physics/src/lib.rs`

## Model details
- AO basis: compact s-type Gaussian primitives (multiple radial functions per atom).
- Integrals: analytic one-electron (`S`, `T`, `V`) and two-electron (`(ij|kl)`) integrals.
- SCF: generalized eigen solve (`F C = S C ε`) with symmetric orthogonalization.
- Correlation: closed-shell MP2 energy correction from MO-transformed two-electron terms.
- Additional observables:
  - HOMO/LUMO energies and gap,
  - Mulliken charges,
  - dipole magnitude (Debye),
  - orbital energy spectrum.

## Commands
```bash
cargo check -q -p gutoe-physics --bin molecular_ab_initio_report
cargo run -q -p gutoe-physics --bin molecular_ab_initio_report
cargo test -q -p gutoe-physics --lib molecular_ab_initio
```

## Outputs
- `/tmp/bh_renders/molecular_ab_initio/molecular_ab_initio_report.csv`
- `/tmp/bh_renders/molecular_ab_initio/molecular_ab_initio_report.json`
- `/tmp/bh_renders/molecular_ab_initio/molecular_ab_initio_report.txt`

## Benchmarks
Solved 8/8:
- H2, LiH, HF, H2O, NH3, CH4, N2, CO2

## Honest limits
- Current basis is compact radial s-primitives only (no explicit p/d angular basis), so this is a forward molecular ab-initio scaffold, not production-grade quantum chemistry accuracy yet.
- RHF lane is closed-shell only (odd-electron open-shell species excluded in this pass).
