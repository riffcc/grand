import Mathlib
import Gutoe.FineStructure
import Gutoe.GaugeConstants

namespace Gutoe.MultiplicativeTriangulation

open Gutoe.FineStructure
open Gutoe.GaugeConstants

/-!
Multiplicative triangulation core:

If positive latent factors `(θ₁, θ₂)` are observed through both a product anchor
and a ratio anchor, then the factors are uniquely recoverable.

This is the formal spine behind log-space triangulation used in the runtime
flavor/EW lane.
-/

/-- Product/ratio map for positive latent factors. -/
def productAnchor (θ₁ θ₂ : ℝ) : ℝ := θ₁ * θ₂

/-- Ratio anchor paired with `productAnchor`. -/
noncomputable def ratioAnchor (θ₁ θ₂ : ℝ) : ℝ := θ₁ / θ₂

/-- Recover `θ₁` from product+ratio anchors (positive branch). -/
theorem theta1_from_product_ratio
    (θ₁ θ₂ : ℝ) (h1 : 0 < θ₁) (h2 : 0 < θ₂) :
    θ₁ = Real.sqrt ((productAnchor θ₁ θ₂) * (ratioAnchor θ₁ θ₂)) := by
  have hsq : (productAnchor θ₁ θ₂) * (ratioAnchor θ₁ θ₂) = θ₁ ^ 2 := by
    unfold productAnchor ratioAnchor
    field_simp [h2.ne']
  calc
    θ₁ = Real.sqrt (θ₁ ^ 2) := by
      rw [Real.sqrt_sq_eq_abs, abs_of_pos h1]
    _ = Real.sqrt ((productAnchor θ₁ θ₂) * (ratioAnchor θ₁ θ₂)) := by rw [hsq]

/-- Recover `θ₂` from product+ratio anchors (positive branch). -/
theorem theta2_from_product_ratio
    (θ₁ θ₂ : ℝ) (h1 : 0 < θ₁) (h2 : 0 < θ₂) :
    θ₂ = Real.sqrt ((productAnchor θ₁ θ₂) / (ratioAnchor θ₁ θ₂)) := by
  have hsq : (productAnchor θ₁ θ₂) / (ratioAnchor θ₁ θ₂) = θ₂ ^ 2 := by
    unfold productAnchor ratioAnchor
    field_simp [h1.ne', h2.ne']
  calc
    θ₂ = Real.sqrt (θ₂ ^ 2) := by
      rw [Real.sqrt_sq_eq_abs, abs_of_pos h2]
    _ = Real.sqrt ((productAnchor θ₁ θ₂) / (ratioAnchor θ₁ θ₂)) := by rw [hsq]

/-- Identifiability: same product and ratio anchors force the same positive factors. -/
theorem product_ratio_identifiable
    (θ₁ θ₂ θ₁' θ₂' : ℝ)
    (h1 : 0 < θ₁) (h2 : 0 < θ₂)
    (h1' : 0 < θ₁') (h2' : 0 < θ₂')
    (hprod : productAnchor θ₁ θ₂ = productAnchor θ₁' θ₂')
    (hratio : ratioAnchor θ₁ θ₂ = ratioAnchor θ₁' θ₂') :
    θ₁ = θ₁' ∧ θ₂ = θ₂' := by
  constructor
  · calc
      θ₁ = Real.sqrt ((productAnchor θ₁ θ₂) * (ratioAnchor θ₁ θ₂)) :=
        theta1_from_product_ratio θ₁ θ₂ h1 h2
      _ = Real.sqrt ((productAnchor θ₁' θ₂') * (ratioAnchor θ₁' θ₂')) := by
        rw [hprod, hratio]
      _ = θ₁' := by
        simpa using (theta1_from_product_ratio θ₁' θ₂' h1' h2').symm
  · calc
      θ₂ = Real.sqrt ((productAnchor θ₁ θ₂) / (ratioAnchor θ₁ θ₂)) :=
        theta2_from_product_ratio θ₁ θ₂ h1 h2
      _ = Real.sqrt ((productAnchor θ₁' θ₂') / (ratioAnchor θ₁' θ₂')) := by
        rw [hprod, hratio]
      _ = θ₂' := by
        simpa using (theta2_from_product_ratio θ₁' θ₂' h1' h2').symm

/-- Structural EW shift is multiplicative in `α²` and `(d/2)`. -/
theorem weinberg_shift_multiplicative_form :
    weinbergMZStructuralShiftQ =
      (1 / (alphaInverse 4 : ℚ)) ^ 2 * ((2 ^ 4 : ℚ) / 2) := by
  rfl

/-- Numeric closure for the structural EW shift factorization. -/
theorem weinberg_shift_multiplicative_closed :
    (1 / (alphaInverse 4 : ℚ)) ^ 2 * ((2 ^ 4 : ℚ) / 2) = 8 / (137 ^ 2 : ℚ) := by
  rw [alpha_inverse_d4]
  norm_num

end Gutoe.MultiplicativeTriangulation
