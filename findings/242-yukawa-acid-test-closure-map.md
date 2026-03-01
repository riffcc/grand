# 242 — Yukawa Acid Test Closure Map

## Scope
Run a single integrated "Yukawa acid test" over current GUTOE lanes to separate:
- what is already constrained by Cl(1,3) structure,
- what is one-anchor reconstructable,
- what is still missing physics (non-integrable ratio cycles).

Binary:
- `crates/gutoe-em/src/bin/yukawa_acid_test_report.rs`

Artifacts:
- `/tmp/bh_renders/yukawa_acid_test_report.txt`
- `/tmp/bh_renders/yukawa_acid_test_report.json`

## Setup
- `v` from Fermi relation: `v = 246.219650794 GeV`.
- Yukawas computed via `y_f = sqrt(2) m_f / v`.
- Structural lanes used directly from existing modules:
  - charged-lepton lane from `mp/me` anchor + Koide phase lane,
  - quark-ratio lane from shared Clifford/Z3 ratio formulas,
  - neutrino absolute scale lane from texture + structural suppression.

## Results

### 1) Charged leptons (structural lane) — strong
- `e`: `0.0083%` relative error
- `mu`: `0.2312%` relative error
- `tau`: `0.2138%` relative error

Interpretation: charged-lepton Yukawa hierarchy is already tightly constrained by the current structural lane.

### 2) Quark ratios vs coarse PDG-like references — mixed
- best lanes near exact:
  - `m_t/m_c`: `0.023%`
  - `m_c/m_u`: `0.274%`
  - `m_t/m_b`: `1.11%`
  - `m_u/m_d`: `1.74%`
- weaker lanes:
  - `m_s/m_d`: `4.59%`
  - `m_c/m_s`: `12.56%`
  - `m_b/m_s`: `14.44%`

Interpretation: quark hierarchy shape is directionally correct, but not uniformly closed.

### 3) Quark cycle consistency (acid test) — fails integrability
Cycle mismatches:
- `(mc/mu)*(mu/md)` vs `(mc/ms)*(ms/md)`: `21.63%`
- `(mt/mc)*(mc/ms)` vs `(mt/mb)*(mb/ms)`: `24.45%`
- full chained cycle: `8.11%`

Interpretation: current ratio set is overdetermined but not exactly integrable into a single consistent absolute mass map. This is the current Yukawa closure gap.

### 4) One-anchor quark reconstruction (least squares, anchored at `m_d`)
Max absolute-mass relative error: `6.145%`.

Representative outputs:
- `u`: `1.61%`
- `s`: `1.34%`
- `c`: `5.12%`
- `b`: `6.14%`
- `t`: `0.90%`

Interpretation: with one absolute anchor and current ratio equations, quark Yukawas are reconstructable to few-percent level, but not yet exact/zero-parameter closure.

### 5) Neutrino Yukawa scale lane
- `m3 ~ 7.90e-3 eV` gives `y_nu3 ~ 4.54e-14`.
- `sum(m_nu) ~ 9.27e-3 eV`.

Interpretation: neutrino lane stays tiny and hierarchical as expected from current structural suppression path.

## Headline
The acid test is now quantified:
- **Lepton Yukawas: already tight** (`<0.25%` max error in this lane).
- **Quark Yukawas: close but not yet integrable** (cycle inconsistency up to `~24%`).
- **One-anchor quark closure exists** at few-percent, signaling a missing correction term rather than total model failure.

## Next step (physics work, not bookkeeping)
Add the missing correction to quark-ratio cycles (likely shared higher-order term) and rerun this exact report until cycle mismatch is driven toward zero.
