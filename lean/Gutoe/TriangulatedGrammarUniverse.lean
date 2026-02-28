import Mathlib
import Gutoe.TriangulatedConstants
import Gutoe.TriangulatedClosureUniqueness

namespace Gutoe.TriangulatedGrammarUniverse

open Gutoe.TriangulatedConstants
open Gutoe.TriangulatedClosureUniqueness

/-!
Universe-level closure:

1. Define widened finite supergrammars for `p`, `kappa`, and `ew` candidate families.
2. Prove each near-solution set is a singleton in that supergrammar.
3. Lift uniqueness to *every* subgrammar via subset inheritance.
-/

def pTolQ : ℚ := 1 / 50000
def kappaTolQ : ℚ := 1 / 50000
def ewTolQ : ℚ := 1 / 1000000

def zBand : Finset ℤ := Finset.Icc (-8) 8
def zBandBase : Finset ℤ := Finset.Icc 0 2

/-- Supergrammar for `p` expressions with two correction channels. -/
def pUniverse : Finset (ℤ × ℤ) := zBand.product zBand

def pExprQ (t : ℤ × ℤ) : ℚ :=
  let a := t.1
  let b := t.2
  (137 / 10 : ℚ) + (a : ℚ) * (1 / (7 * 12 : ℚ)) + (b : ℚ) * (1 / (7 * 13 * 136 : ℚ))

def pCandidateTuple : ℤ × ℤ := (-1, 0)

theorem p_expr_candidate_eq :
    pExprQ pCandidateTuple = pCandidateQ := by
  unfold pExprQ pCandidateTuple
  rw [p_candidate_closed_form]
  ring

def pNearSet : Finset (ℤ × ℤ) :=
  pUniverse.filter (fun t => |pExprQ t - pFrozenQ| < pTolQ)

theorem p_near_set_singleton :
    pNearSet = ({pCandidateTuple} : Finset (ℤ × ℤ)) := by
  native_decide

theorem p_candidate_near_frozen :
    |pExprQ pCandidateTuple - pFrozenQ| < pTolQ := by
  have hmem : pCandidateTuple ∈ pNearSet := by
    rw [p_near_set_singleton]
    simp
  exact (Finset.mem_filter.mp hmem).2

theorem p_universe_unique_near
    {t : ℤ × ℤ}
    (ht : t ∈ pUniverse)
    (hnear : |pExprQ t - pFrozenQ| < pTolQ) :
    t = pCandidateTuple := by
  have hmem : t ∈ pNearSet := Finset.mem_filter.mpr ⟨ht, hnear⟩
  rw [p_near_set_singleton] at hmem
  simpa using hmem

theorem p_any_subgrammar_unique
    (G : Finset (ℤ × ℤ))
    (hsub : G ⊆ pUniverse)
    {t : ℤ × ℤ}
    (ht : t ∈ G)
    (hnear : |pExprQ t - pFrozenQ| < pTolQ) :
    t = pCandidateTuple :=
  p_universe_unique_near (hsub ht) hnear

theorem p_any_subgrammar_has_solution
    (G : Finset (ℤ × ℤ))
    (_hsub : G ⊆ pUniverse)
    (hCand : pCandidateTuple ∈ G) :
    ∃ t ∈ G, |pExprQ t - pFrozenQ| < pTolQ := by
  refine ⟨pCandidateTuple, hCand, ?_⟩
  exact p_candidate_near_frozen

/-- Supergrammar for `kappa` expressions with signed coefficients. -/
def kappaUniverse : Finset (ℤ × (ℤ × ℤ)) :=
  zBand.product (zBand.product zBand)

def kappaExprQ (t : ℤ × (ℤ × ℤ)) : ℚ :=
  let a := t.1
  let b := t.2.1
  let c := t.2.2
  (60 / 11 : ℚ)
    * ((a : ℚ) * (19 / 3 : ℚ)
      + (b : ℚ) * (1 / 36 : ℚ)
      + (c : ℚ) * (1 / (7 * 13 * 136 : ℚ)))

def kappaCandidateTuple : ℤ × (ℤ × ℤ) := (1, (1, 1))

theorem kappa_expr_candidate_eq :
    kappaExprQ kappaCandidateTuple = kappaCandidateQ := by
  unfold kappaExprQ kappaCandidateTuple
  rw [kappa_candidate_closed_form]
  ring

def kappaNearSet : Finset (ℤ × (ℤ × ℤ)) :=
  kappaUniverse.filter (fun t => |kappaExprQ t - kappaFrozenQ| < kappaTolQ)

theorem kappa_near_set_singleton :
    kappaNearSet = ({kappaCandidateTuple} : Finset (ℤ × (ℤ × ℤ))) := by
  native_decide

theorem kappa_candidate_near_frozen :
    |kappaExprQ kappaCandidateTuple - kappaFrozenQ| < kappaTolQ := by
  have hmem : kappaCandidateTuple ∈ kappaNearSet := by
    rw [kappa_near_set_singleton]
    simp
  exact (Finset.mem_filter.mp hmem).2

theorem kappa_universe_unique_near
    {t : ℤ × (ℤ × ℤ)}
    (ht : t ∈ kappaUniverse)
    (hnear : |kappaExprQ t - kappaFrozenQ| < kappaTolQ) :
    t = kappaCandidateTuple := by
  have hmem : t ∈ kappaNearSet := Finset.mem_filter.mpr ⟨ht, hnear⟩
  rw [kappa_near_set_singleton] at hmem
  simpa using hmem

theorem kappa_any_subgrammar_unique
    (G : Finset (ℤ × (ℤ × ℤ)))
    (hsub : G ⊆ kappaUniverse)
    {t : ℤ × (ℤ × ℤ)}
    (ht : t ∈ G)
    (hnear : |kappaExprQ t - kappaFrozenQ| < kappaTolQ) :
    t = kappaCandidateTuple :=
  kappa_universe_unique_near (hsub ht) hnear

theorem kappa_any_subgrammar_has_solution
    (G : Finset (ℤ × (ℤ × ℤ)))
    (_hsub : G ⊆ kappaUniverse)
    (hCand : kappaCandidateTuple ∈ G) :
    ∃ t ∈ G, |kappaExprQ t - kappaFrozenQ| < kappaTolQ := by
  refine ⟨kappaCandidateTuple, hCand, ?_⟩
  exact kappa_candidate_near_frozen

/-- Supergrammar for EW coefficient expressions. -/
def ewUniverse : Finset (ℤ × (ℤ × ℤ)) :=
  zBandBase.product (zBand.product zBand)

def ewExprQ (t : ℤ × (ℤ × ℤ)) : ℚ :=
  let m := t.1
  let b := t.2.1
  let c := t.2.2
  (m : ℚ) * (8 : ℚ) + (b : ℚ) * (6 / 13 : ℚ) + (c : ℚ) * (1 / (7 * 136 : ℚ))

def ewCandidateTuple : ℤ × (ℤ × ℤ) := (1, (1, -1))

theorem ew_expr_candidate_eq :
    ewExprQ ewCandidateTuple = ewCoeffCandidateQ := by
  unfold ewExprQ ewCandidateTuple
  rw [ew_coeff_candidate_closed_form]
  ring

def ewNearSet : Finset (ℤ × (ℤ × ℤ)) :=
  ewUniverse.filter (fun t => |ewExprQ t - ewCoeffFrozenQ| < ewTolQ)

theorem ew_near_set_singleton :
    ewNearSet = ({ewCandidateTuple} : Finset (ℤ × (ℤ × ℤ))) := by
  native_decide

theorem ew_candidate_near_frozen :
    |ewExprQ ewCandidateTuple - ewCoeffFrozenQ| < ewTolQ := by
  have hmem : ewCandidateTuple ∈ ewNearSet := by
    rw [ew_near_set_singleton]
    simp
  exact (Finset.mem_filter.mp hmem).2

theorem ew_universe_unique_near
    {t : ℤ × (ℤ × ℤ)}
    (ht : t ∈ ewUniverse)
    (hnear : |ewExprQ t - ewCoeffFrozenQ| < ewTolQ) :
    t = ewCandidateTuple := by
  have hmem : t ∈ ewNearSet := Finset.mem_filter.mpr ⟨ht, hnear⟩
  rw [ew_near_set_singleton] at hmem
  simpa using hmem

theorem ew_any_subgrammar_unique
    (G : Finset (ℤ × (ℤ × ℤ)))
    (hsub : G ⊆ ewUniverse)
    {t : ℤ × (ℤ × ℤ)}
    (ht : t ∈ G)
    (hnear : |ewExprQ t - ewCoeffFrozenQ| < ewTolQ) :
    t = ewCandidateTuple :=
  ew_universe_unique_near (hsub ht) hnear

theorem ew_any_subgrammar_has_solution
    (G : Finset (ℤ × (ℤ × ℤ)))
    (_hsub : G ⊆ ewUniverse)
    (hCand : ewCandidateTuple ∈ G) :
    ∃ t ∈ G, |ewExprQ t - ewCoeffFrozenQ| < ewTolQ := by
  refine ⟨ewCandidateTuple, hCand, ?_⟩
  exact ew_candidate_near_frozen

/-- Global closure theorem:
    any subgrammar of the supergrammar families has a unique near solution
    (when it contains the candidate tuple). -/
theorem all_subgrammars_unique_if_candidate_included :
    (∀ (Gp : Finset (ℤ × ℤ)),
        Gp ⊆ pUniverse → pCandidateTuple ∈ Gp →
        (∃ t ∈ Gp, |pExprQ t - pFrozenQ| < pTolQ) ∧
        (∀ t ∈ Gp, |pExprQ t - pFrozenQ| < pTolQ → t = pCandidateTuple)) ∧
    (∀ (Gk : Finset (ℤ × (ℤ × ℤ))),
        Gk ⊆ kappaUniverse → kappaCandidateTuple ∈ Gk →
        (∃ t ∈ Gk, |kappaExprQ t - kappaFrozenQ| < kappaTolQ) ∧
        (∀ t ∈ Gk, |kappaExprQ t - kappaFrozenQ| < kappaTolQ → t = kappaCandidateTuple)) ∧
    (∀ (Ge : Finset (ℤ × (ℤ × ℤ))),
        Ge ⊆ ewUniverse → ewCandidateTuple ∈ Ge →
        (∃ t ∈ Ge, |ewExprQ t - ewCoeffFrozenQ| < ewTolQ) ∧
        (∀ t ∈ Ge, |ewExprQ t - ewCoeffFrozenQ| < ewTolQ → t = ewCandidateTuple)) := by
  constructor
  · intro Gp hsub hCand
    constructor
    · exact p_any_subgrammar_has_solution Gp hsub hCand
    · intro t ht hnear
      exact p_any_subgrammar_unique Gp hsub ht hnear
  · constructor
    · intro Gk hsub hCand
      constructor
      · exact kappa_any_subgrammar_has_solution Gk hsub hCand
      · intro t ht hnear
        exact kappa_any_subgrammar_unique Gk hsub ht hnear
    · intro Ge hsub hCand
      constructor
      · exact ew_any_subgrammar_has_solution Ge hsub hCand
      · intro t ht hnear
        exact ew_any_subgrammar_unique Ge hsub ht hnear

/-!
Expanded ("for the road") supergrammar closure:

Same structural operators, strictly wider coefficient domains.
This strengthens uniqueness claims by proving the singleton near-solution
survives a much larger search region.
-/

set_option maxRecDepth 200000

def zBandXL : Finset ℤ := Finset.Icc (-32) 32
def zBandBaseXL : Finset ℤ := Finset.Icc 0 5
def zBandEwXL : Finset ℤ := Finset.Icc (-48) 48

def pUniverseXL : Finset (ℤ × ℤ) := zBandXL.product zBandXL

def pNearSetXL : Finset (ℤ × ℤ) :=
  pUniverseXL.filter (fun t => |pExprQ t - pFrozenQ| < pTolQ)

theorem p_near_set_xl_singleton :
    pNearSetXL = ({pCandidateTuple} : Finset (ℤ × ℤ)) := by
  native_decide

theorem p_xl_universe_unique_near
    {t : ℤ × ℤ}
    (ht : t ∈ pUniverseXL)
    (hnear : |pExprQ t - pFrozenQ| < pTolQ) :
    t = pCandidateTuple := by
  have hmem : t ∈ pNearSetXL := Finset.mem_filter.mpr ⟨ht, hnear⟩
  rw [p_near_set_xl_singleton] at hmem
  simpa using hmem

theorem p_xl_any_subgrammar_unique
    (G : Finset (ℤ × ℤ))
    (hsub : G ⊆ pUniverseXL)
    {t : ℤ × ℤ}
    (ht : t ∈ G)
    (hnear : |pExprQ t - pFrozenQ| < pTolQ) :
    t = pCandidateTuple :=
  p_xl_universe_unique_near (hsub ht) hnear

theorem p_xl_any_subgrammar_has_solution
    (G : Finset (ℤ × ℤ))
    (_hsub : G ⊆ pUniverseXL)
    (hCand : pCandidateTuple ∈ G) :
    ∃ t ∈ G, |pExprQ t - pFrozenQ| < pTolQ := by
  refine ⟨pCandidateTuple, hCand, ?_⟩
  exact p_candidate_near_frozen

def kappaUniverseXL : Finset (ℤ × (ℤ × ℤ)) :=
  zBandXL.product (zBandXL.product zBandXL)

def kappaNearSetXL : Finset (ℤ × (ℤ × ℤ)) :=
  kappaUniverseXL.filter (fun t => |kappaExprQ t - kappaFrozenQ| < kappaTolQ)

theorem kappa_near_set_xl_singleton :
    kappaNearSetXL = ({kappaCandidateTuple} : Finset (ℤ × (ℤ × ℤ))) := by
  native_decide

theorem kappa_xl_universe_unique_near
    {t : ℤ × (ℤ × ℤ)}
    (ht : t ∈ kappaUniverseXL)
    (hnear : |kappaExprQ t - kappaFrozenQ| < kappaTolQ) :
    t = kappaCandidateTuple := by
  have hmem : t ∈ kappaNearSetXL := Finset.mem_filter.mpr ⟨ht, hnear⟩
  rw [kappa_near_set_xl_singleton] at hmem
  simpa using hmem

theorem kappa_xl_any_subgrammar_unique
    (G : Finset (ℤ × (ℤ × ℤ)))
    (hsub : G ⊆ kappaUniverseXL)
    {t : ℤ × (ℤ × ℤ)}
    (ht : t ∈ G)
    (hnear : |kappaExprQ t - kappaFrozenQ| < kappaTolQ) :
    t = kappaCandidateTuple :=
  kappa_xl_universe_unique_near (hsub ht) hnear

theorem kappa_xl_any_subgrammar_has_solution
    (G : Finset (ℤ × (ℤ × ℤ)))
    (_hsub : G ⊆ kappaUniverseXL)
    (hCand : kappaCandidateTuple ∈ G) :
    ∃ t ∈ G, |kappaExprQ t - kappaFrozenQ| < kappaTolQ := by
  refine ⟨kappaCandidateTuple, hCand, ?_⟩
  exact kappa_candidate_near_frozen

def ewUniverseXL : Finset (ℤ × (ℤ × ℤ)) :=
  zBandBaseXL.product (zBandEwXL.product zBandEwXL)

def ewNearSetXL : Finset (ℤ × (ℤ × ℤ)) :=
  ewUniverseXL.filter (fun t => |ewExprQ t - ewCoeffFrozenQ| < ewTolQ)

theorem ew_near_set_xl_singleton :
    ewNearSetXL = ({ewCandidateTuple} : Finset (ℤ × (ℤ × ℤ))) := by
  native_decide

theorem ew_xl_universe_unique_near
    {t : ℤ × (ℤ × ℤ)}
    (ht : t ∈ ewUniverseXL)
    (hnear : |ewExprQ t - ewCoeffFrozenQ| < ewTolQ) :
    t = ewCandidateTuple := by
  have hmem : t ∈ ewNearSetXL := Finset.mem_filter.mpr ⟨ht, hnear⟩
  rw [ew_near_set_xl_singleton] at hmem
  simpa using hmem

theorem ew_xl_any_subgrammar_unique
    (G : Finset (ℤ × (ℤ × ℤ)))
    (hsub : G ⊆ ewUniverseXL)
    {t : ℤ × (ℤ × ℤ)}
    (ht : t ∈ G)
    (hnear : |ewExprQ t - ewCoeffFrozenQ| < ewTolQ) :
    t = ewCandidateTuple :=
  ew_xl_universe_unique_near (hsub ht) hnear

theorem ew_xl_any_subgrammar_has_solution
    (G : Finset (ℤ × (ℤ × ℤ)))
    (_hsub : G ⊆ ewUniverseXL)
    (hCand : ewCandidateTuple ∈ G) :
    ∃ t ∈ G, |ewExprQ t - ewCoeffFrozenQ| < ewTolQ := by
  refine ⟨ewCandidateTuple, hCand, ?_⟩
  exact ew_candidate_near_frozen

theorem legacy_universes_embed_in_xl :
    pUniverse ⊆ pUniverseXL ∧
    kappaUniverse ⊆ kappaUniverseXL ∧
    ewUniverse ⊆ ewUniverseXL := by
  native_decide

theorem all_subgrammars_unique_if_candidate_included_xl :
    (∀ (Gp : Finset (ℤ × ℤ)),
        Gp ⊆ pUniverseXL → pCandidateTuple ∈ Gp →
        (∃ t ∈ Gp, |pExprQ t - pFrozenQ| < pTolQ) ∧
        (∀ t ∈ Gp, |pExprQ t - pFrozenQ| < pTolQ → t = pCandidateTuple)) ∧
    (∀ (Gk : Finset (ℤ × (ℤ × ℤ))),
        Gk ⊆ kappaUniverseXL → kappaCandidateTuple ∈ Gk →
        (∃ t ∈ Gk, |kappaExprQ t - kappaFrozenQ| < kappaTolQ) ∧
        (∀ t ∈ Gk, |kappaExprQ t - kappaFrozenQ| < kappaTolQ → t = kappaCandidateTuple)) ∧
    (∀ (Ge : Finset (ℤ × (ℤ × ℤ))),
        Ge ⊆ ewUniverseXL → ewCandidateTuple ∈ Ge →
        (∃ t ∈ Ge, |ewExprQ t - ewCoeffFrozenQ| < ewTolQ) ∧
        (∀ t ∈ Ge, |ewExprQ t - ewCoeffFrozenQ| < ewTolQ → t = ewCandidateTuple)) := by
  constructor
  · intro Gp hsub hCand
    constructor
    · exact p_xl_any_subgrammar_has_solution Gp hsub hCand
    · intro t ht hnear
      exact p_xl_any_subgrammar_unique Gp hsub ht hnear
  · constructor
    · intro Gk hsub hCand
      constructor
      · exact kappa_xl_any_subgrammar_has_solution Gk hsub hCand
      · intro t ht hnear
        exact kappa_xl_any_subgrammar_unique Gk hsub ht hnear
    · intro Ge hsub hCand
      constructor
      · exact ew_xl_any_subgrammar_has_solution Ge hsub hCand
      · intro t ht hnear
        exact ew_xl_any_subgrammar_unique Ge hsub ht hnear

/-!
Infinite-space elimination by index ranges:

Instead of searching tuples directly, collapse each affine grammar to a single
integer lattice index. Near-target constraints then force that index into a
tiny finite interval.
-/

def pIndex (a b : ℤ) : ℤ := 442 * a + 3 * b

theorem pExprQ_eq_index (a b : ℤ) :
    pExprQ (a, b) = (137 / 10 : ℚ) + (pIndex a b : ℚ) / 37128 := by
  change (137 / 10 : ℚ) + (a : ℚ) * (1 / (7 * 12 : ℚ)) + (b : ℚ) * (1 / (7 * 13 * 136 : ℚ)) =
    (137 / 10 : ℚ) + ((442 * a + 3 * b : ℤ) : ℚ) / 37128
  have hcast : (((442 * a + 3 * b : ℤ) : ℚ)) = (a : ℚ) * 442 + (b : ℚ) * 3 := by
    norm_num [Int.cast_add, Int.cast_mul, add_comm, add_left_comm, add_assoc, mul_comm, mul_left_comm, mul_assoc]
  rw [hcast]
  ring

theorem p_near_forces_index_window
    {a b : ℤ}
    (hnear : |pExprQ (a, b) - pFrozenQ| < pTolQ) :
    pIndex a b = -442 ∨ pIndex a b = -441 := by
  let n : ℤ := pIndex a b
  have hnear' : |(n : ℚ) / 37128 - (pFrozenQ - (137 / 10 : ℚ))| < pTolQ := by
    simpa [n, pExprQ_eq_index, sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hnear
  have hlt := abs_lt.mp hnear'
  have hlow : (pFrozenQ - (137 / 10 : ℚ)) - pTolQ < (n : ℚ) / 37128 := by
    linarith [hlt.1]
  have hhigh : (n : ℚ) / 37128 < (pFrozenQ - (137 / 10 : ℚ)) + pTolQ := by
    linarith [hlt.2]
  have hLBconst : (-443 : ℚ) / 37128 < (pFrozenQ - (137 / 10 : ℚ)) - pTolQ := by
    native_decide
  have hUBconst : (pFrozenQ - (137 / 10 : ℚ)) + pTolQ < (-440 : ℚ) / 37128 := by
    native_decide
  have hge : -442 ≤ n := by
    by_contra hneg
    have hle : n ≤ -443 := by omega
    have hqle : (n : ℚ) ≤ (-443 : ℚ) := by exact_mod_cast hle
    have hqdiv : (n : ℚ) / 37128 ≤ (-443 : ℚ) / 37128 := by
      nlinarith
    linarith
  have hle : n ≤ -441 := by
    by_contra hpos
    have hge' : -440 ≤ n := by omega
    have hqge : (-440 : ℚ) ≤ (n : ℚ) := by exact_mod_cast hge'
    have hqdiv : (-440 : ℚ) / 37128 ≤ (n : ℚ) / 37128 := by
      nlinarith
    linarith
  have : n = -442 ∨ n = -441 := by omega
  simpa [n] using this

def kappaIndex (a b c : ℤ) : ℤ := 3527160 * a + 15470 * b + 45 * c

theorem kappaExprQ_eq_index (a b c : ℤ) :
    kappaExprQ (a, (b, c)) = (kappaIndex a b c : ℚ) / 102102 := by
  change (60 / 11 : ℚ) * ((a : ℚ) * (19 / 3 : ℚ) + (b : ℚ) * (1 / 36 : ℚ) + (c : ℚ) * (1 / (7 * 13 * 136 : ℚ))) =
    (((3527160 * a + 15470 * b + 45 * c : ℤ) : ℚ) / 102102)
  have hcast : (((3527160 * a + 15470 * b + 45 * c : ℤ) : ℚ)) =
      (a : ℚ) * 3527160 + (b : ℚ) * 15470 + (c : ℚ) * 45 := by
    norm_num [Int.cast_add, Int.cast_mul, add_comm, add_left_comm, add_assoc, mul_comm, mul_left_comm, mul_assoc]
  rw [hcast]
  ring

theorem kappa_near_forces_index_window
    {a b c : ℤ}
    (hnear : |kappaExprQ (a, (b, c)) - kappaFrozenQ| < kappaTolQ) :
    3542672 ≤ kappaIndex a b c ∧ kappaIndex a b c ≤ 3542675 := by
  let n : ℤ := kappaIndex a b c
  have hnear' : |(n : ℚ) / 102102 - kappaFrozenQ| < kappaTolQ := by
    simpa [n, kappaExprQ_eq_index] using hnear
  have hlt := abs_lt.mp hnear'
  have hlow : kappaFrozenQ - kappaTolQ < (n : ℚ) / 102102 := by
    linarith [hlt.1]
  have hhigh : (n : ℚ) / 102102 < kappaFrozenQ + kappaTolQ := by
    linarith [hlt.2]
  have hLBconst : (3542671 : ℚ) / 102102 < kappaFrozenQ - kappaTolQ := by
    native_decide
  have hUBconst : kappaFrozenQ + kappaTolQ < (3542676 : ℚ) / 102102 := by
    native_decide
  have hge : 3542672 ≤ n := by
    by_contra hneg
    have hle : n ≤ 3542671 := by omega
    have hqle : (n : ℚ) ≤ (3542671 : ℚ) := by exact_mod_cast hle
    have hqdiv : (n : ℚ) / 102102 ≤ (3542671 : ℚ) / 102102 := by
      nlinarith
    linarith
  have hle : n ≤ 3542675 := by
    by_contra hpos
    have hge' : 3542676 ≤ n := by omega
    have hqge : (3542676 : ℚ) ≤ (n : ℚ) := by exact_mod_cast hge'
    have hqdiv : (3542676 : ℚ) / 102102 ≤ (n : ℚ) / 102102 := by
      nlinarith
    linarith
  simpa [n] using And.intro hge hle

def ewIndex (m b c : ℤ) : ℤ := 99008 * m + 5712 * b + 13 * c

theorem ewExprQ_eq_index (m b c : ℤ) :
    ewExprQ (m, (b, c)) = (ewIndex m b c : ℚ) / 12376 := by
  change (m : ℚ) * (8 : ℚ) + (b : ℚ) * (6 / 13 : ℚ) + (c : ℚ) * (1 / (7 * 136 : ℚ)) =
    (((99008 * m + 5712 * b + 13 * c : ℤ) : ℚ) / 12376)
  have hcast : (((99008 * m + 5712 * b + 13 * c : ℤ) : ℚ)) =
      (m : ℚ) * 99008 + (b : ℚ) * 5712 + (c : ℚ) * 13 := by
    norm_num [Int.cast_add, Int.cast_mul, add_comm, add_left_comm, add_assoc, mul_comm, mul_left_comm, mul_assoc]
  rw [hcast]
  ring

theorem ew_near_forces_index_exact
    {m b c : ℤ}
    (hnear : |ewExprQ (m, (b, c)) - ewCoeffFrozenQ| < ewTolQ) :
    ewIndex m b c = 104707 := by
  let n : ℤ := ewIndex m b c
  have hnear' : |(n : ℚ) / 12376 - ewCoeffFrozenQ| < ewTolQ := by
    simpa [n, ewExprQ_eq_index] using hnear
  have hlt := abs_lt.mp hnear'
  have hlow : ewCoeffFrozenQ - ewTolQ < (n : ℚ) / 12376 := by
    linarith [hlt.1]
  have hhigh : (n : ℚ) / 12376 < ewCoeffFrozenQ + ewTolQ := by
    linarith [hlt.2]
  have hLBconst : (104706 : ℚ) / 12376 < ewCoeffFrozenQ - ewTolQ := by
    native_decide
  have hUBconst : ewCoeffFrozenQ + ewTolQ < (104708 : ℚ) / 12376 := by
    native_decide
  have hge : 104707 ≤ n := by
    by_contra hneg
    have hle : n ≤ 104706 := by omega
    have hqle : (n : ℚ) ≤ (104706 : ℚ) := by exact_mod_cast hle
    have hqdiv : (n : ℚ) / 12376 ≤ (104706 : ℚ) / 12376 := by
      nlinarith
    linarith
  have hle : n ≤ 104707 := by
    by_contra hpos
    have hge' : 104708 ≤ n := by omega
    have hqge : (104708 : ℚ) ≤ (n : ℚ) := by exact_mod_cast hge'
    have hqdiv : (104708 : ℚ) / 12376 ≤ (n : ℚ) / 12376 := by
      nlinarith
    linarith
  omega

/-!
Canonical selection from elimination windows:

Use the index-window certificates plus XL structural bounds to recover unique
candidate tuples.
-/

def pIndexWindowSetXL : Finset (ℤ × ℤ) :=
  pUniverseXL.filter (fun t => pIndex t.1 t.2 = -442 ∨ pIndex t.1 t.2 = -441)

theorem p_index_window_xl_singleton :
    pIndexWindowSetXL = ({pCandidateTuple} : Finset (ℤ × ℤ)) := by
  native_decide

theorem p_xl_bounds_and_index_select_candidate
    {a b : ℤ}
    (ha : a ∈ zBandXL)
    (hb : b ∈ zBandXL)
    (hidx : pIndex a b = -442 ∨ pIndex a b = -441) :
    (a, b) = pCandidateTuple := by
  have hmem : (a, b) ∈ pIndexWindowSetXL := by
    exact Finset.mem_filter.mpr ⟨Finset.mem_product.mpr ⟨ha, hb⟩, hidx⟩
  rw [p_index_window_xl_singleton] at hmem
  simpa using hmem

theorem p_infinite_elimination_xl_complete
    {a b : ℤ}
    (ha : a ∈ zBandXL)
    (hb : b ∈ zBandXL)
    (hnear : |pExprQ (a, b) - pFrozenQ| < pTolQ) :
    (a, b) = pCandidateTuple := by
  exact p_xl_bounds_and_index_select_candidate ha hb (p_near_forces_index_window hnear)

def kappaIndexWindowSetXL : Finset (ℤ × (ℤ × ℤ)) :=
  kappaUniverseXL.filter (fun t =>
    3542672 ≤ kappaIndex t.1 t.2.1 t.2.2 ∧ kappaIndex t.1 t.2.1 t.2.2 ≤ 3542675)

theorem kappa_index_window_xl_singleton :
    kappaIndexWindowSetXL = ({kappaCandidateTuple} : Finset (ℤ × (ℤ × ℤ))) := by
  native_decide

theorem kappa_xl_bounds_and_index_select_candidate
    {a b c : ℤ}
    (ha : a ∈ zBandXL)
    (hb : b ∈ zBandXL)
    (hc : c ∈ zBandXL)
    (hidx : 3542672 ≤ kappaIndex a b c ∧ kappaIndex a b c ≤ 3542675) :
    (a, (b, c)) = kappaCandidateTuple := by
  have hmem : (a, (b, c)) ∈ kappaIndexWindowSetXL := by
    refine Finset.mem_filter.mpr ?_
    exact ⟨Finset.mem_product.mpr ⟨ha, Finset.mem_product.mpr ⟨hb, hc⟩⟩, hidx⟩
  rw [kappa_index_window_xl_singleton] at hmem
  simpa using hmem

theorem kappa_infinite_elimination_xl_complete
    {a b c : ℤ}
    (ha : a ∈ zBandXL)
    (hb : b ∈ zBandXL)
    (hc : c ∈ zBandXL)
    (hnear : |kappaExprQ (a, (b, c)) - kappaFrozenQ| < kappaTolQ) :
    (a, (b, c)) = kappaCandidateTuple := by
  exact kappa_xl_bounds_and_index_select_candidate ha hb hc
    (kappa_near_forces_index_window hnear)

def ewIndexExactSetXL : Finset (ℤ × (ℤ × ℤ)) :=
  ewUniverseXL.filter (fun t => ewIndex t.1 t.2.1 t.2.2 = 104707)

theorem ew_index_exact_xl_singleton :
    ewIndexExactSetXL = ({ewCandidateTuple} : Finset (ℤ × (ℤ × ℤ))) := by
  native_decide

theorem ew_xl_bounds_and_index_select_candidate
    {m b c : ℤ}
    (hm : m ∈ zBandBaseXL)
    (hb : b ∈ zBandEwXL)
    (hc : c ∈ zBandEwXL)
    (hidx : ewIndex m b c = 104707) :
    (m, (b, c)) = ewCandidateTuple := by
  have hmem : (m, (b, c)) ∈ ewIndexExactSetXL := by
    refine Finset.mem_filter.mpr ?_
    exact ⟨Finset.mem_product.mpr ⟨hm, Finset.mem_product.mpr ⟨hb, hc⟩⟩, hidx⟩
  rw [ew_index_exact_xl_singleton] at hmem
  simpa using hmem

theorem ew_infinite_elimination_xl_complete
    {m b c : ℤ}
    (hm : m ∈ zBandBaseXL)
    (hb : b ∈ zBandEwXL)
    (hc : c ∈ zBandEwXL)
    (hnear : |ewExprQ (m, (b, c)) - ewCoeffFrozenQ| < ewTolQ) :
    (m, (b, c)) = ewCandidateTuple := by
  exact ew_xl_bounds_and_index_select_candidate hm hb hc
    (ew_near_forces_index_exact hnear)

end Gutoe.TriangulatedGrammarUniverse
