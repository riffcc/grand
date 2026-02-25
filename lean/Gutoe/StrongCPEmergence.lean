/- 
 * GUTOE — Strong CP Emergence Bridge (GRAND-267)
 *
 * This module formalizes the key non-repopulation bridge:
 *
 * If effective-sector topological charge is a pullback of a fundamental
 * Z3-carrier charge, and fundamental continuous maps are constant (from the
 * route-1 vacuum theorem), then emergent sectors on that image cannot carry
 * nonzero topological charge.
 *
 * In that case, theta-dependent phases are identically unity on the emergent
 * image, so theta is unphysical there.
 -/

import Mathlib
import Gutoe.StrongCPVacuum

namespace Gutoe.StrongCPEmergence

open Gutoe.StrongCPVacuum

/-- A charge functional on the emergent image cannot repopulate nonzero sectors
    when its source space is subsingleton and a single base state is normalized
    to zero charge. This captures the "no topology created from nothing" step
    without requiring an explicit pointwise coarse-graining formula. -/
theorem no_repopulation_from_subsingleton_source
    {S X Eff : Type}
    [Subsingleton S]
    [TopologicalSpace X]
    [TopologicalSpace Eff]
    (toEff : S → C(X, Eff))
    (qEff : C(X, Eff) → ℤ)
    (s0 : S)
    (hbase : qEff (toEff s0) = 0) :
    ∀ s : S, qEff (toEff s) = 0 := by
  intro s
  have hs : s = s0 := Subsingleton.elim s s0
  rw [hs]
  exact hbase

/-- Canonical based fundamental map: constant map to `0 : Fin 3`. -/
def basedFundamentalZero {X : Type} [TopologicalSpace X] (x0 : X) :
    {f : C(X, FundamentalGaugeGroup) // f x0 = 0} :=
  ⟨ContinuousMap.const X 0, by simp⟩

/-- Route-1 no-repopulation for based sectors:
    since based fundamental maps into the discrete Z3 carrier are subsingleton,
    any emergent charge functional is pinned by its value on one base state. -/
theorem no_repopulation_from_based_route1
    {X Eff : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    [TopologicalSpace Eff]
    (x0 : X)
    (toEff : {f : C(X, FundamentalGaugeGroup) // f x0 = 0} → C(X, Eff))
    (qEff : C(X, Eff) → ℤ)
    (hbase : qEff (toEff (basedFundamentalZero x0)) = 0) :
    ∀ fb : {f : C(X, FundamentalGaugeGroup) // f x0 = 0}, qEff (toEff fb) = 0 := by
  have hs : Subsingleton {f : C(X, FundamentalGaugeGroup) // f x0 = 0} :=
    based_fundamental_sector_unique x0
  letI : Subsingleton {f : C(X, FundamentalGaugeGroup) // f x0 = 0} := hs
  exact no_repopulation_from_subsingleton_source toEff qEff (basedFundamentalZero x0) hbase

/-- Homotopy-invariance route:
    if effective charge is homotopy-invariant on the emergent image and the
    emergence map preserves homotopy classes, then route-1 constantness of
    fundamental maps forces zero emergent charge. -/
theorem no_repopulation_of_homotopy_invariance
    {X Eff : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    [TopologicalSpace Eff]
    (toEff : C(X, FundamentalGaugeGroup) → C(X, Eff))
    (qEff : C(X, Eff) → ℤ)
    (hPreserve :
      ∀ f g : C(X, FundamentalGaugeGroup),
        f.Homotopic g → (toEff f).Homotopic (toEff g))
    (hInv :
      ∀ F G : C(X, Eff), F.Homotopic G → qEff F = qEff G)
    (hConstZero :
      ∀ g0 : FundamentalGaugeGroup, qEff (toEff (ContinuousMap.const X g0)) = 0) :
    ∀ f : C(X, FundamentalGaugeGroup), qEff (toEff f) = 0 := by
  intro f
  rcases continuous_to_fundamental_group_constant f with ⟨g0, hg0⟩
  have heq : f = ContinuousMap.const X g0 := by
    ext x
    simpa using congrArg (fun y : FundamentalGaugeGroup => (y : ℕ)) (hg0 x)
  have hhomFund : f.Homotopic (ContinuousMap.const X g0) := by
    simpa [heq] using (ContinuousMap.Homotopic.refl (ContinuousMap.const X g0))
  have hhomEff : (toEff f).Homotopic (toEff (ContinuousMap.const X g0)) :=
    hPreserve f (ContinuousMap.const X g0) hhomFund
  calc
    qEff (toEff f) = qEff (toEff (ContinuousMap.const X g0)) := hInv _ _ hhomEff
    _ = 0 := hConstZero g0

/-- Route-1 no-repopulation theorem:
    if emergent charge is pulled back from fundamental Z3-carrier maps and
    constant fundamental maps have zero charge, then emergent-image charge is zero. -/
theorem no_repopulation_on_emergent_image
    {X Eff : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    [TopologicalSpace Eff]
    (toEff : C(X, FundamentalGaugeGroup) → C(X, Eff))
    (qEff : C(X, Eff) → ℤ)
    (qFund : C(X, FundamentalGaugeGroup) → ℤ)
    (hqPullback : ∀ f, qEff (toEff f) = qFund f)
    (hqFundConstZero :
      ∀ f, (∃ g0 : FundamentalGaugeGroup, ∀ x : X, f x = g0) → qFund f = 0) :
    ∀ f : C(X, FundamentalGaugeGroup), qEff (toEff f) = 0 := by
  intro f
  rcases continuous_to_fundamental_group_constant f with ⟨g0, hg0⟩
  calc
    qEff (toEff f) = qFund f := hqPullback f
    _ = 0 := hqFundConstZero f ⟨g0, hg0⟩

/-- Theta-phase factor for integer topological charge. -/
noncomputable def thetaPhase (theta : ℝ) (q : ℤ) : ℂ :=
  Complex.exp (((theta * (q : ℝ)) : ℂ) * Complex.I)

/-- Zero charge gives a trivial theta phase. -/
theorem theta_phase_unity_of_zero_charge (theta : ℝ) :
    thetaPhase theta 0 = 1 := by
  unfold thetaPhase
  norm_num

/-- Emergent-image theta phases are unity under the no-repopulation hypotheses. -/
theorem theta_phase_unity_on_emergent_image
    {X Eff : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    [TopologicalSpace Eff]
    (toEff : C(X, FundamentalGaugeGroup) → C(X, Eff))
    (qEff : C(X, Eff) → ℤ)
    (qFund : C(X, FundamentalGaugeGroup) → ℤ)
    (hqPullback : ∀ f, qEff (toEff f) = qFund f)
    (hqFundConstZero :
      ∀ f, (∃ g0 : FundamentalGaugeGroup, ∀ x : X, f x = g0) → qFund f = 0) :
    ∀ (f : C(X, FundamentalGaugeGroup)) (theta : ℝ),
      thetaPhase theta (qEff (toEff f)) = 1 := by
  intro f theta
  have hq0 : qEff (toEff f) = 0 :=
    no_repopulation_on_emergent_image toEff qEff qFund hqPullback hqFundConstZero f
  rw [hq0]
  exact theta_phase_unity_of_zero_charge theta

/-- Coarse-graining no-creation theorem:
    if coarse-graining is pointwise (`CG f = φ ∘ f`) and effective charge
    vanishes on constant fields, then no nonzero effective sector can be
    created from fundamental Z3-carrier fields on a preconnected domain. -/
theorem coarse_grain_cannot_create_nontrivial_sector
    {X Eff : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    [TopologicalSpace Eff]
    (CG : C(X, FundamentalGaugeGroup) → C(X, Eff))
    (phi : C(FundamentalGaugeGroup, Eff))
    (hCG : ∀ f, CG f = phi.comp f)
    (qEff : C(X, Eff) → ℤ)
    (hqConstZero : ∀ e0 : Eff, qEff (ContinuousMap.const X e0) = 0) :
    ∀ f : C(X, FundamentalGaugeGroup), qEff (CG f) = 0 := by
  intro f
  rcases continuous_to_fundamental_group_constant f with ⟨g0, hg0⟩
  have hconst : CG f = ContinuousMap.const X (phi g0) := by
    rw [hCG f]
    ext x
    simpa [hg0 x]
  rw [hconst]
  exact hqConstZero (phi g0)

/-- Theta-phase corollary of coarse-grain no-creation. -/
theorem theta_phase_unity_of_coarse_grain_no_creation
    {X Eff : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    [TopologicalSpace Eff]
    (CG : C(X, FundamentalGaugeGroup) → C(X, Eff))
    (phi : C(FundamentalGaugeGroup, Eff))
    (hCG : ∀ f, CG f = phi.comp f)
    (qEff : C(X, Eff) → ℤ)
    (hqConstZero : ∀ e0 : Eff, qEff (ContinuousMap.const X e0) = 0) :
    ∀ (f : C(X, FundamentalGaugeGroup)) (theta : ℝ),
      thetaPhase theta (qEff (CG f)) = 1 := by
  intro f theta
  have hq0 : qEff (CG f) = 0 :=
    coarse_grain_cannot_create_nontrivial_sector CG phi hCG qEff hqConstZero f
  rw [hq0]
  exact theta_phase_unity_of_zero_charge theta

/-- SU(3)-matrix specialization of the coarse-grain no-repopulation theorem. -/
theorem su3_effective_no_repopulation
    {X : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    (CG : C(X, FundamentalGaugeGroup) → C(X, Matrix (Fin 3) (Fin 3) ℂ))
    (phi : C(FundamentalGaugeGroup, Matrix (Fin 3) (Fin 3) ℂ))
    (hCG : ∀ f, CG f = phi.comp f)
    (qEff : C(X, Matrix (Fin 3) (Fin 3) ℂ) → ℤ)
    (hqConstZero : ∀ e0 : Matrix (Fin 3) (Fin 3) ℂ, qEff (ContinuousMap.const X e0) = 0) :
    ∀ f : C(X, FundamentalGaugeGroup), qEff (CG f) = 0 := by
  intro f
  exact coarse_grain_cannot_create_nontrivial_sector CG phi hCG qEff hqConstZero f

/-- SU(3)-matrix specialization of theta-phase collapse on emergent image. -/
theorem su3_effective_theta_phase_unity
    {X : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    (CG : C(X, FundamentalGaugeGroup) → C(X, Matrix (Fin 3) (Fin 3) ℂ))
    (phi : C(FundamentalGaugeGroup, Matrix (Fin 3) (Fin 3) ℂ))
    (hCG : ∀ f, CG f = phi.comp f)
    (qEff : C(X, Matrix (Fin 3) (Fin 3) ℂ) → ℤ)
    (hqConstZero : ∀ e0 : Matrix (Fin 3) (Fin 3) ℂ, qEff (ContinuousMap.const X e0) = 0) :
    ∀ (f : C(X, FundamentalGaugeGroup)) (theta : ℝ),
      thetaPhase theta (qEff (CG f)) = 1 := by
  intro f theta
  exact theta_phase_unity_of_coarse_grain_no_creation CG phi hCG qEff hqConstZero f theta

end Gutoe.StrongCPEmergence
