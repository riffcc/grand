# Finding 266 — RH Lean Proof Scaffold + Ablation Hardening

Date: 2026-03-01  
Scope: Execute RH proof-program kickoff end-to-end: formal Lean scaffold + fairness ablations.

## Added (Lean)

New Lean modules:

- `lean/Gutoe/RiemannCore.lean`
- `lean/Gutoe/RiemannSelfAdjoint.lean`
- `lean/Gutoe/RiemannBridge.lean`
- `lean/Gutoe/RiemannCounting.lean`
- `lean/Gutoe/RiemannRHClosure.lean`

Wired into `lean/lakefile.lean` roots:

- `Gutoe.RiemannCore`
- `Gutoe.RiemannSelfAdjoint`
- `Gutoe.RiemannBridge`
- `Gutoe.RiemannCounting`
- `Gutoe.RiemannRHClosure`

Build check:

```bash
cd lean
lake build Gutoe
```

Status: **passes** (with pre-existing repo linter warnings unrelated to this lane).

## Formal Content Delivered

This is a real reduction scaffold (no `sorry`, no fake closure):

1. **Core RH predicate (abstract `Xi`)**
   - `RiemannHypothesisXi : (ℂ → ℂ) → Prop`
   - critical-line embedding `t ↦ 1/2 + i t`
2. **Reduction theorem**
   - `rh_of_zero_parameterization`
   - If every zero is parameterized as `1/2 + i t`, RH-for-`Xi` follows.
3. **Finite structural operator lane**
   - structural tridiagonal matrix defined from shared constants
   - symmetry/self-adjoint proxy theorem:
     `structuralRiemannMatrix_finiteSelfAdjoint`
4. **Exact bridge theorem**
   - `bridge_implies_rh`
   - exact spectral bridge implies RH-for-`Xi`
5. **Counting algebra**
   - finite `countUpTo`
   - monotonicity + match transport theorems
6. **Program pack**
   - `RHProgramAssumptions`
   - closure theorem `rh_of_program_assumptions`

Interpretation: this closes the *proof plumbing* layer, not RH itself.

## Added (Runtime Ablation)

New binary:

- `crates/gutoe-physics/src/bin/riemann_nail_ablation_report.rs`

Run:

```bash
GUTOE_RIEMANN_REF_PATH=/tmp/bh_renders/zeta_zeros_first_1000_odlyzko.txt \
GUTOE_RIEMANN_NS=512,1024,2048 \
cargo run -q -p gutoe-physics --bin riemann_nail_ablation_report
```

Artifacts:

- `/tmp/bh_renders/riemann_nail_ablation_report.txt`
- `/tmp/bh_renders/riemann_nail_ablation_report.json`

## Ablation Results (equal model complexity)

Protocol: same quadratic-map capacity, same objective, same branch scan.
Compare:

- truth spectrum (structural operator eigenvalues)
- scrambled spectrum control
- linear surrogate control

Core hold+freeze error (`41..120`) and relative gain:

- `n=512`
  - truth: `3.622e-02`
  - scrambled: `1.196e+00`
  - linear: `5.218e-01`
  - truth gain vs scrambled: `96.97%`
  - truth gain vs linear: `93.06%`

- `n=1024`
  - truth: `1.554e-02`
  - scrambled: `1.196e+00`
  - linear: `5.218e-01`
  - truth gain vs scrambled: `98.70%`
  - truth gain vs linear: `97.02%`

- `n=2048`
  - truth: `7.269e-03`
  - scrambled: `1.191e+00`
  - linear: `4.682e-02`
  - truth gain vs scrambled: `99.39%`
  - truth gain vs linear: `84.47%`

Long holdout (`121..500`) at `n=2048`:

- truth: `2.003e-02`
- scrambled: `8.472e-01`
- linear: `6.784e-01`

This is a hard fairness win for the structural spectrum lane.

## Honest Status

- RH is **not yet proven**.
- What is now true:
  - theorem-stack plumbing exists in Lean and compiles in the full library;
  - runtime lane survives stronger equal-capacity ablations.
- Remaining decisive work is still the analytic identity bridge (`Xi ↔ spectrum(H)`) at theorem level.

