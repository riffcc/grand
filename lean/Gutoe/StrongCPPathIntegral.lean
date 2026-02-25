/-
 * GUTOE — Strong CP Path-Integral Skeleton (GRAND-267)
 *
 * This module adds a finite-sector path-integral skeleton tied directly to
 * the Cl(1,3) bivector split used in StrongCP:
 *
 * - + sector multiplicity := |magneticTriplet| = 3
 * - - sector multiplicity := |emTriplet| = 3
 *
 * Using the standard ±θ phase decomposition, the CP-odd channel is the
 * imaginary-part coefficient:
 *
 *   Z_im(θ) = (N₊ - N₋) * sin θ
 *
 * Since N₊ = N₋ in Cl(1,3), Z_im(θ)=0 for all θ in this skeleton.
 *
 * Scope note:
 * This is a structural finite-sector model, not yet the full QCD path integral.
 -/

import Mathlib
import Gutoe.LorentzInvariance
import Gutoe.StrongCP

namespace Gutoe.StrongCPPathIntegral

open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness
open Gutoe.LorentzInvariance
open Gutoe.StrongCP

/-- Positive CP-sector multiplicity in the skeleton model. -/
def plusSectorMultiplicity : ℤ := (magneticTriplet.card : ℤ)

/-- Negative CP-sector multiplicity in the skeleton model. -/
def minusSectorMultiplicity : ℤ := (emTriplet.card : ℤ)

/-- Cl(1,3) gives equal ± sector multiplicities (3 and 3). -/
theorem sector_multiplicities_equal :
    plusSectorMultiplicity = minusSectorMultiplicity := by
  unfold plusSectorMultiplicity minusSectorMultiplicity
  rcases lorentz_algebra_decomposition with ⟨_, _, hmag, hem⟩
  rw [hmag, hem]

/-- Explicit multiplicity values from the same decomposition. -/
theorem sector_multiplicities_are_three :
    plusSectorMultiplicity = 3 ∧ minusSectorMultiplicity = 3 := by
  unfold plusSectorMultiplicity minusSectorMultiplicity
  rcases lorentz_algebra_decomposition with ⟨_, _, hmag, hem⟩
  exact ⟨by simpa [hmag], by simpa [hem]⟩

/-- Real-channel coefficient in the finite-sector skeleton. -/
noncomputable def zRe (theta : ℝ) : ℝ :=
  (plusSectorMultiplicity : ℝ) * Real.cos theta +
  (minusSectorMultiplicity : ℝ) * Real.cos theta

/-- Imaginary-channel coefficient in the finite-sector skeleton. -/
noncomputable def zIm (theta : ℝ) : ℝ :=
  ((plusSectorMultiplicity - minusSectorMultiplicity : ℤ) : ℝ) * Real.sin theta

/-- CP-odd channel factorizes as multiplicity imbalance times `sin θ`. -/
theorem zIm_factorized (theta : ℝ) :
    zIm theta = ((plusSectorMultiplicity - minusSectorMultiplicity : ℤ) : ℝ) * Real.sin theta := by
  rfl

/-- The imbalance term is exactly the StrongCP structural source. -/
theorem cp_odd_source_matches_strongcp :
    (plusSectorMultiplicity - minusSectorMultiplicity : ℤ) = cpOddSectorImbalance := by
  unfold plusSectorMultiplicity minusSectorMultiplicity cpOddSectorImbalance
  ring

/-- Cl(1,3) sector balance kills the CP-odd channel for every θ. -/
theorem zIm_zero_all_theta (theta : ℝ) : zIm theta = 0 := by
  rw [zIm_factorized]
  rw [cp_odd_source_matches_strongcp]
  rw [cp_odd_sector_imbalance_zero]
  norm_num

/-- Real channel is even under `θ ↦ -θ`. -/
theorem zRe_even_theta_reflection (theta : ℝ) : zRe (-theta) = zRe theta := by
  unfold zRe
  simp [Real.cos_neg]

/-- Imaginary channel is odd under `θ ↦ -θ`. -/
theorem zIm_odd_theta_reflection (theta : ℝ) : zIm (-theta) = -zIm theta := by
  rw [zIm_factorized, zIm_factorized, Real.sin_neg]
  ring

/-- Complex skeleton partition coefficient. -/
noncomputable def zComplex (theta : ℝ) : ℂ :=
  ((zRe theta : ℝ) : ℂ) + (((zIm theta : ℝ) : ℂ) * Complex.I)

/-- In the Cl(1,3) balanced skeleton, the complex coefficient is `θ`-reflection even. -/
theorem zComplex_even_theta_reflection (theta : ℝ) :
    zComplex (-theta) = zComplex theta := by
  unfold zComplex
  rw [zRe_even_theta_reflection theta]
  rw [zIm_zero_all_theta (-theta), zIm_zero_all_theta theta]

/-- Equivalent expression using the structural θ proxy from StrongCP. -/
theorem zIm_as_structural_theta (theta : ℝ) :
    zIm theta = thetaQcdStructural * Real.sin theta := by
  rw [zIm_factorized]
  rw [cp_odd_source_matches_strongcp]
  unfold thetaQcdStructural
  ring

/-- Real channel collapses to `6*cos θ` from the 3+3 split. -/
theorem zRe_eq_six_cos (theta : ℝ) : zRe theta = 6 * Real.cos theta := by
  unfold zRe
  rcases sector_multiplicities_are_three with ⟨hplus, hminus⟩
  rw [hplus, hminus]
  ring

/-- At θ=0 the skeleton partition coefficient is strictly positive. -/
theorem zRe_theta_zero_positive : zRe 0 = 6 := by
  rw [zRe_eq_six_cos]
  norm_num

/-- At `θ = π`, the real channel is negative in this skeleton. -/
theorem zRe_theta_pi_negative : zRe Real.pi = -6 := by
  rw [zRe_eq_six_cos]
  norm_num

/-- Route-2 candidate set used in the principal-branch exclusion step. -/
def route2Candidate (theta : ℝ) : Prop := theta = 0 ∨ theta = Real.pi

/-- Nonnegative-weight variant of route-2 exclusion:
    if candidates are `{0, π}` and `zRe(θ) ≥ 0`, then `θ = 0`. -/
theorem theta_zero_of_discrete_candidates_and_nonnegative_weight
    (theta : ℝ)
    (hdisc : route2Candidate theta)
    (hnonneg : 0 ≤ zRe theta) :
    theta = 0 := by
  rcases hdisc with h0 | hpi
  · exact h0
  · exfalso
    have hneg : zRe theta = -6 := by simpa [hpi] using zRe_theta_pi_negative
    rw [hneg] at hnonneg
    linarith

/-- Route-2 exclusion step in the principal branch:
    if a candidate is constrained to `{0, π}` and the vacuum weight is positive,
    then `θ = 0` (the `π` branch is excluded). -/
theorem theta_zero_of_discrete_candidates_and_positive_weight
    (theta : ℝ)
    (hdisc : theta = 0 ∨ theta = Real.pi)
    (hpos : 0 < zRe theta) :
    theta = 0 := by
  rcases hdisc with h0 | hpi
  · exact h0
  · exfalso
    have hneg : zRe theta = -6 := by simpa [hpi] using zRe_theta_pi_negative
    rw [hneg] at hpos
    linarith

/-- Structural `θ_QCD` is a route-2 candidate (in fact the `0` branch). -/
theorem theta_qcd_structural_route2_candidate :
    route2Candidate thetaQcdStructural := by
  left
  exact theta_qcd_structural_zero

/-- Structural route-2 closure, phrased through the nonnegative-weight gate. -/
theorem theta_qcd_structural_route2_closed :
    thetaQcdStructural = 0 := by
  exact
    theta_zero_of_discrete_candidates_and_nonnegative_weight
      thetaQcdStructural
      theta_qcd_structural_route2_candidate
      (by
        rw [theta_qcd_structural_zero, zRe_theta_zero_positive]
        norm_num)

end Gutoe.StrongCPPathIntegral
