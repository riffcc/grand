# VICE_RULES.md — Project Constitution

**VICE = Verification-Integrated Constraint Enforcement**

Every Lean theorem must tighten constraints, not add freedom. This is the project constitution.

---

## 1. Zero Parameter Discipline

- No new free parameters unless forced and explicitly named
- Every physical quantity must trace to Cl(1,3) combinatorics via Lean definitions
- The single input is `m_p = 938.272 MeV`. Everything else is derived
- If you add a parameter, justify why Cl(1,3) doesn't determine it

---

## 2. No Vacuous Theorems

A theorem must have a **non-trivial proposition**. `True` is never acceptable.

**The GaugeConstants.lean violation:**
```lean
-- BAD: proves nothing
theorem wilson_loop_area_law : True := by trivial

-- BAD: tautology
theorem u1_hypercharge : 1 = 1 := by rfl

-- BAD: wrong statement (simplified away the physics)
theorem charge_sum_per_generation : (1 : ℚ) + (-1) = 0 := by norm_num
```

**What to do instead:**
```lean
-- GOOD: honest gap with clear statement
theorem confinement_from_area_law
    (σ : ℝ) (hσ : σ > 0)
    (area_law : ∀ C, wilson_expectation C = Real.exp (-σ * area C)) :
    Potential V(r) ∼ σ * r := sorry

-- GOOD: actual charge sum from Z₃ structure
theorem charge_sum_per_generation :
    3 * (2/3 : ℚ) + 3 * (-1/3 : ℚ) + (-1 : ℚ) + (0 : ℚ) = 0 := by norm_num
```

If you can't prove it yet, use `sorry` — but state the actual claim. A fake proof of `True` is worse than an honest `sorry`.

---

## 3. Bridges = Theorems, Not Comments

Every cross-module connection must be a Lean theorem or explicitly tagged as an open gap.

**Bad:**
```lean
-- This connects to InstantonMass...
-- (no actual theorem)
```

**Good:**
```lean
-- Actual bridge theorem
theorem clifford_forces_instanton_scale
    (h : Z3SymmetryBroken) :
    InstantonScale = magneticTriplet.card * LatticeScale := by sorry
```

If a connection is only informal, create a `sorry`'d theorem and file a Plane issue.

---

## 4. Shared Primitives Must Be Shared

Key objects must flow through ALL modules as the same Lean term:

| Name | Definition | Location |
|------|------------|----------|
| `cliffordDim` | `2^4 = 16` | DimensionalStructure |
| `magneticTriplet` | `{7, 11, 13}` | Z3Uniqueness |
| `grade1_4d` | `{2, 3, 5, 9}` | DimensionalStructure |
| `grade2_4d` | `{4, 6, 7, 10, 11, 13}` | Z3Uniqueness |
| `triangularNumber` | `n(n+1)/2` | FineStructure |
| `z3_4d` | ℕ → ℕ permutation | DimensionalStructure |

**Rule:** Import, don't redefine. No hardcoded `16` when `cliffordDim` exists.

---

## 5. Axiom Hygiene

Every `axiom` or `postulate` must be tagged with:
- **Justification**: Why can't this be proved?
- **Degree-of-freedom cost**: What does assuming this buy us?

`sorry` is a temporary axiom. It must have a corresponding Plane issue.

---

## 6. The Minimax Rule

Weaker models (tagged 🤖 minimax-safe) may only work on items that are:
- Pure algebra / finite combinatorics
- Provable with `norm_num`, `decide`, `native_decide`, `ring`, `omega`
- NOT involving real analysis, limits, or continuous functions

**Minimax output MUST be reviewed** against these VICE rules before items are closed.

A theorem that compiles but proves something trivially different from what was asked is a **VICE violation**.

---

## 7. CI Gate

- `lake build` must pass (zero errors)
- `grep sorry` count must not increase without a corresponding Plane issue
- `grep axiom` count must not increase without justification in PR
- All new theorems must import from existing shared primitives where applicable

---

## Exhibit A: What We Fixed

GaugeConstants.lean had 3 vacuous theorems:

1. `wilson_loop_area_law : True` — removed
2. `u1_hypercharge : 1 = 1` — removed
3. `charge_sum_per_generation : (1:ℚ) + (-1) = 0` — replaced with actual physics

The solid theorems (`beta_zero`, `cos_sq_theta_w`, `mZ_over_mW_sq`, `total_gauge_bosons`) connect to real infrastructure.

---

**Violation of these rules is worse than `sorry` — it creates false confidence.**
