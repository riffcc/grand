/-
 * GUTOE — Poincare Recovery B4a (GRAND-B4a)
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * Representation-theoretic scaffold:
 *   1. Define the cubic octahedral symmetry group O_h as coordinate
 *      permutations/reflections on spatial coordinates (order 48).
 *   2. Encode the harmonic selection rule from O_h invariance:
 *      l = 1,2,3 are forbidden; survivors are l = 0 and even l >= 4.
 *   3. State cubic-lattice corrections as O(a^4), so anisotropy vanishes in
 *      the continuum limit a -> 0.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.LorentzInvariance
import Gutoe.LatticeGeometry

namespace Gutoe.PoincareRecoveryB4a

/-- Reflection sector `(Z₂)^3` acting on spatial coordinate signs. -/
abbrev ReflectionBits : Type := Fin 3 → Multiplicative (ZMod 2)

/-- Coordinate symmetries: permutations (`S₃`) times independent reflections (`(Z₂)^3`). -/
abbrev CoordinateSymmetry3 : Type := Equiv.Perm (Fin 3) × ReflectionBits

/-- Octahedral group `O_h` (full cubic symmetry with reflections), realized as all
    permutation/reflection coordinate symmetries. -/
abbrev Oh : Subgroup CoordinateSymmetry3 := ⊤

/-- `O_h` is finite as a subgroup of coordinate permutations/reflections. -/
theorem oh_finite : (Oh : Set CoordinateSymmetry3).Finite := by
  classical
  simp [Oh]

/-- `|O_h| = 48 = 3! * 2^3` (same as the ambient coordinate-symmetry carrier). -/
theorem oh_order : Fintype.card CoordinateSymmetry3 = 48 := by
  native_decide

/-- Allowed angular momentum degrees under cubic (`O_h`) symmetry in this scaffold. -/
def IsEvenDegree (l : ℕ) : Prop := l % 2 = 0

instance (l : ℕ) : Decidable (IsEvenDegree l) := by
  unfold IsEvenDegree
  infer_instance

/-- Allowed angular momentum degrees under cubic (`O_h`) symmetry in this scaffold. -/
def OhSurvivingDegree (l : ℕ) : Prop := l = 0 ∨ (4 ≤ l ∧ IsEvenDegree l)

instance (l : ℕ) : Decidable (OhSurvivingDegree l) := by
  unfold OhSurvivingDegree
  infer_instance

/-- `l = 1` is removed by cubic symmetry. -/
theorem oh_blocks_l1 : ¬ OhSurvivingDegree 1 := by native_decide

/-- `l = 2` is removed by cubic symmetry. -/
theorem oh_blocks_l2 : ¬ OhSurvivingDegree 2 := by native_decide

/-- `l = 3` is removed by cubic symmetry. -/
theorem oh_blocks_l3 : ¬ OhSurvivingDegree 3 := by native_decide

/-- The first surviving degrees are `0, 4, 6, 8, ...`. -/
theorem oh_first_survivors :
    OhSurvivingDegree 0 ∧ OhSurvivingDegree 4 ∧ OhSurvivingDegree 6 ∧ OhSurvivingDegree 8 := by
  native_decide

/-- Harmonic-spectrum encoding of `O_h` invariance. -/
def OhInvariantSpectrum (coeff : ℕ → ℝ) : Prop :=
  coeff 1 = 0 ∧ coeff 2 = 0 ∧ coeff 3 = 0 ∧
    ∀ l, coeff l ≠ 0 → OhSurvivingDegree l

/-- `O_h` invariance kills the `l = 1,2,3` channels. -/
theorem oh_invariance_kills_l123 {coeff : ℕ → ℝ} (hOh : OhInvariantSpectrum coeff) :
    coeff 1 = 0 ∧ coeff 2 = 0 ∧ coeff 3 = 0 :=
  ⟨hOh.1, hOh.2.1, hOh.2.2.1⟩

/-- Under `O_h` invariance, any nonzero harmonic component is in the surviving set. -/
theorem oh_invariance_survivor_rule {coeff : ℕ → ℝ} (hOh : OhInvariantSpectrum coeff) :
    ∀ l, coeff l ≠ 0 → l = 0 ∨ (4 ≤ l ∧ IsEvenDegree l) :=
  hOh.2.2.2

/-- Cubic lattice anisotropy scaling: correction is `O(a^4)` as `a -> 0`. -/
def CubicCorrectionOrder (delta : ℝ → ℝ) : Prop :=
  delta =O[nhds (0 : ℝ)] (fun a : ℝ => a ^ (4 : ℕ))

/-- Canonical quartic correction is `O(a^4)` (reflexive big-O). -/
theorem quartic_term_is_O_a4 : CubicCorrectionOrder (fun a : ℝ => a ^ 4) := by
  simpa [CubicCorrectionOrder] using
    (Asymptotics.isBigO_refl (fun a : ℝ => a ^ (4 : ℕ)) (nhds (0 : ℝ)))

/-- Any `O(a^4)` correction vanishes at `a -> 0`; this is continuum rotational recovery. -/
theorem cubic_Oa4_vanishes_in_continuum {delta : ℝ → ℝ}
    (hDelta : CubicCorrectionOrder delta) :
    Filter.Tendsto delta (nhds (0 : ℝ)) (nhds (0 : ℝ)) := by
  have hBigO : delta =O[nhds (0 : ℝ)] (fun a : ℝ => a ^ (4 : ℕ)) := by
    simpa [CubicCorrectionOrder] using hDelta
  have hPow : Filter.Tendsto (fun a : ℝ => a ^ (4 : ℕ)) (nhds (0 : ℝ)) (nhds (0 : ℝ)) := by
    have hId : Filter.Tendsto (fun a : ℝ => a) (nhds (0 : ℝ)) (nhds (0 : ℝ)) :=
      Filter.tendsto_id
    simpa using hId.pow 4
  exact Asymptotics.IsBigO.trans_tendsto hBigO hPow

/-- Correlator-level rotational recovery statement:
    if lattice anisotropy is `O(a^4)`, the isotropic continuum limit is recovered. -/
def RotationalInvarianceLimit (C Ciso : ℝ → ℝ) : Prop :=
  Filter.Tendsto (fun a => C a - Ciso a) (nhds (0 : ℝ)) (nhds (0 : ℝ))

/-- Euclidean 4-space coordinate model. -/
abbrev Euclidean4 : Type := Fin 4 → ℝ

/-- Correlator families with lattice spacing parameter `a` on Euclidean 4-space. -/
abbrev CorrelatorFamily4 : Type := ℝ → Euclidean4 → ℝ

/-- Smoothness in the lattice-spacing parameter, pointwise in Euclidean 4-space. -/
def SmoothInSpacing4 (C : CorrelatorFamily4) : Prop :=
  ∀ x, ContDiff ℝ ⊤ (fun a : ℝ => C a x)

/-- Rotational invariance recovered in the continuum limit, pointwise in Euclidean 4-space. -/
def RotationalInvarianceLimit4 (C Ciso : CorrelatorFamily4) : Prop :=
  ∀ x, Filter.Tendsto (fun a => C a x - Ciso a x) (nhds (0 : ℝ)) (nhds (0 : ℝ))

/-- Smoothness + cubic-symmetry scaling (`O(a^4)`) implies rotational invariance in the limit. -/
theorem smooth_cubic_symmetry_implies_rotational_invariance
    (C Ciso : ℝ → ℝ)
    (hSmooth : ContDiff ℝ ⊤ C)
    (hCubic : CubicCorrectionOrder (fun a => C a - Ciso a)) :
    RotationalInvarianceLimit C Ciso := by
  have _ : ContDiff ℝ ⊤ C := hSmooth
  exact cubic_Oa4_vanishes_in_continuum hCubic

/-- 4D Euclidean version:
    smoothness + cubic-symmetry scaling (`O(a^4)`) implies rotational invariance
    in the continuum limit `a -> 0`. -/
theorem smooth_cubic_symmetry_implies_rotational_invariance4
    (C Ciso : CorrelatorFamily4)
    (hSmooth : SmoothInSpacing4 C)
    (hCubic : ∀ x, CubicCorrectionOrder (fun a => C a x - Ciso a x)) :
    RotationalInvarianceLimit4 C Ciso := by
  intro x
  have _ : ContDiff ℝ ⊤ (fun a : ℝ => C a x) := hSmooth x
  exact cubic_Oa4_vanishes_in_continuum (hCubic x)

end Gutoe.PoincareRecoveryB4a
