# CLAUDE.md — Project Context

## Project Constitution

**Read [VICE_RULES.md](./VICE_RULES.md) before writing any Lean code.**

## Lean Workflow

**Use `/lean-prover` when working on Lean proofs.** This skill provides Lean-specific guidance and will help fix compilation errors.

To activate: type `/lean-prover` in the conversation.

---

## VICE Rules Summary

VICE = "Verification-Integrated Constraint Enforcement". Every theorem must:
1. State a non-trivial proposition (no `True`, no tautologies)
2. Use shared Lean definitions (no hardcoded literals when named terms exist)
3. Trace physical claims to Cl(1,3) structure
4. Be honest about gaps (`sorry` with a Plane issue > fake proof of `True`)

Hollow theorems that compile but prove nothing are worse than `sorry` — they create false confidence.

---

## What is GUTOE?

GUTOE (Grand Unifying Theory of Everything) is a zero-free-parameter physics theory derived entirely from the Clifford algebra Cl(1,3).

### Core Claims
- α⁻¹ = 137 from T(16) + 1 (triangular number)
- sin²θ_W = 3/13 from Z₃ orbit structure
- n_gen = 3 from |Z₃|
- SU(3)×SU(2)×U(1) from grade decomposition
- Simple cubic lattice from spatial bivector counting

### Architecture
- **Lean proofs** (`lean/Gutoe/`): Formal proofs of combinatorial claims
- **Rust simulations** (`crates/gutoe-*/`): Monte Carlo, GPU solvers
- **Shared primitives**: `grade1_4d`, `grade2_4d`, `magneticTriplet`, `z3_4d`

### Key Files
| File | Content |
|------|---------|
| `Z3Uniqueness.lean` | Z₃ forces 1 lepton + 3 quarks |
| `FineStructure.lean` | α⁻¹ = 137 |
| `LatticeGeometry.lean` | SC lattice from bivectors |
| `GaugeConstants.lean` | β₀, Weinberg angle |
| `GaugeGroupSM.lean` | SU(3)×SU(2)×U(1) |

---

## TODO

- **Load `/mnt/riffcastle/gaia_dr3.mvl` and calculate the entire galaxy** 🌌
  Gaia DR3 full source catalog. Use `gaia_dr3_life_map.rs` as the starting point.
  Goal: run GUTOE Great Filter + life probability across every resolved star.

---

## Communication Style

- Be direct and blunt about quality
- Flag violations of VICE rules immediately
- No `sorry` without a corresponding Plane issue
- Prioritize constraint propagation over feature addition
