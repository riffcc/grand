/- 
 * GUTOE — Strong CP General-Case Split
 *
 * This module formalizes the two global regimes discussed in GRAND-267:
 *
 * (A) GUTOE route-1 emergent-image regime:
 *     θ is unphysical (all θ phases collapse on the physical image).
 *
 * (B) Standard-QCD nonzero-topological-sector regime:
 *     θ is physical (there exist θ choices with distinct phase factors).
 -/

import Mathlib
import Gutoe.StrongCPEmergence
import Gutoe.StrongCPPathIntegral

namespace Gutoe.StrongCPGeneralCases

open Gutoe.StrongCPEmergence
open Gutoe.StrongCPVacuum

/-- For any integer charge sector, `θ = 0` gives a trivial phase factor. -/
theorem theta_phase_at_zero (q : ℤ) : thetaPhase 0 q = 1 := by
  unfold thetaPhase
  simp

/-- Any nonzero integer topological charge yields a nontrivial phase for
    `θ = π / q` (the phase is `-1`, hence not `1`). -/
theorem theta_phase_nontrivial_of_nonzero_charge
    (q : ℤ) (hq : q ≠ 0) :
    ∃ theta : ℝ, thetaPhase theta q ≠ 1 := by
  let theta : ℝ := Real.pi / (q : ℝ)
  have hqR : (q : ℝ) ≠ 0 := by exact_mod_cast hq
  refine ⟨theta, ?_⟩
  have hmul : theta * (q : ℝ) = Real.pi := by
    calc
      theta * (q : ℝ) = (Real.pi / (q : ℝ)) * (q : ℝ) := by rfl
      _ = Real.pi := by field_simp [hqR]
  have hphase : thetaPhase theta q = (-1 : ℂ) := by
    unfold thetaPhase
    have harg : (theta : ℂ) * (((q : ℝ) : ℂ)) = (Real.pi : ℂ) := by
      exact_mod_cast hmul
    rw [Complex.exp_mul_I]
    change
      Complex.cos ((theta : ℂ) * (((q : ℝ) : ℂ))) +
      Complex.sin ((theta : ℂ) * (((q : ℝ) : ℂ))) * Complex.I = -1
    rw [harg]
    simp
  intro hEq
  have hneg : (-1 : ℂ) ≠ 1 := by norm_num
  exact hneg (hphase.symm.trans hEq)

/-- In any regime where at least one nonzero topological sector is physically
    accessible, `θ` is a physical parameter (phase factors can differ). -/
theorem theta_physical_of_nonzero_sector
    (q : ℤ) (hq : q ≠ 0) :
    ∃ theta1 theta2 : ℝ, thetaPhase theta1 q ≠ thetaPhase theta2 q := by
  rcases theta_phase_nontrivial_of_nonzero_charge q hq with ⟨theta, htheta⟩
  refine ⟨theta, 0, ?_⟩
  intro hEq
  apply htheta
  calc
    thetaPhase theta q = thetaPhase 0 q := hEq
    _ = 1 := theta_phase_at_zero q

/-- GUTOE route-1 general-case statement:
    on the concrete emergent image, all `θ` choices are physically equivalent. -/
theorem gutoe_theta_unphysical_on_emergent_image
    {X : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    (x0 : X)
    (qClass : HomotopyClass X Su3Matrix → ℤ) :
    ∀ (f : C(X, FundamentalGaugeGroup)) (theta1 theta2 : ℝ),
      thetaPhase theta1 (qEffFromClassAnchored x0 qClass (z3ToSu3.comp f)) =
      thetaPhase theta2 (qEffFromClassAnchored x0 qClass (z3ToSu3.comp f)) :=
by
  intro f theta1 theta2
  exact theta_unphysical_concrete x0 qClass f theta1 theta2

/-- Global two-case split (formal):
    (A) emergent-image GUTOE case: `θ` unphysical;
    (B) nonzero-sector case: `θ` physical. -/
theorem strong_cp_general_case_split
    {X : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    (x0 : X)
    (qClass : HomotopyClass X Su3Matrix → ℤ) :
    (∀ (f : C(X, FundamentalGaugeGroup)) (theta1 theta2 : ℝ),
      thetaPhase theta1 (qEffFromClassAnchored x0 qClass (z3ToSu3.comp f)) =
      thetaPhase theta2 (qEffFromClassAnchored x0 qClass (z3ToSu3.comp f))) ∧
    (∀ q : ℤ, q ≠ 0 → ∃ theta1 theta2 : ℝ,
      thetaPhase theta1 q ≠ thetaPhase theta2 q) := by
  constructor
  · exact gutoe_theta_unphysical_on_emergent_image x0 qClass
  · intro q hq
    exact theta_physical_of_nonzero_sector q hq

end Gutoe.StrongCPGeneralCases
