/- 
 * GUTOE — One-Generation Standard Model Anomaly Gates
 *
 * Exact rational anomaly sums over the canonical Weyl registry.
 * These are non-perturbative consistency gates: they either cancel or fail.
-/

import Mathlib
import Gutoe.SM.Rep

namespace Gutoe.SM.Anomalies

open Gutoe.SM.Rep
open scoped BigOperators

/-- SU(3)^2 · U(1)_Y contribution from one Weyl species. -/
def su3su3u1Term (f : WeylSpecies) : ℚ :=
  hypercharge f * dynkinSU3 (colorRep f) * (weakMultiplicitySpecies f : ℚ)

/-- SU(2)^2 · U(1)_Y contribution from one Weyl species. -/
def su2su2u1Term (f : WeylSpecies) : ℚ :=
  hypercharge f * dynkinSU2 (weakRep f) * (colorMultiplicitySpecies f : ℚ)

/-- U(1)_Y^3 contribution from one Weyl species. -/
def u1cubedTerm (f : WeylSpecies) : ℚ :=
  (hypercharge f) ^ 3 * (colorMultiplicitySpecies f : ℚ) * (weakMultiplicitySpecies f : ℚ)

/-- grav^2 · U(1)_Y contribution from one Weyl species. -/
def gravgravu1Term (f : WeylSpecies) : ℚ :=
  hypercharge f * (colorMultiplicitySpecies f : ℚ) * (weakMultiplicitySpecies f : ℚ)

/-- Total SU(3)^2 · U(1)_Y anomaly for one generation. -/
def anomalySU3SU3U1 : ℚ := Finset.sum oneGeneration (fun f => su3su3u1Term f)

/-- Total SU(2)^2 · U(1)_Y anomaly for one generation. -/
def anomalySU2SU2U1 : ℚ := Finset.sum oneGeneration (fun f => su2su2u1Term f)

/-- Total U(1)_Y^3 anomaly for one generation. -/
def anomalyU1Cubed : ℚ := Finset.sum oneGeneration (fun f => u1cubedTerm f)

/-- Total grav^2 · U(1)_Y anomaly for one generation. -/
def anomalyGravGravU1 : ℚ := Finset.sum oneGeneration (fun f => gravgravu1Term f)

/-- One-generation SU(3)^2 · U(1)_Y anomaly cancellation. -/
theorem anomaly_su3su3u1_cancels : anomalySU3SU3U1 = 0 := by
  native_decide

/-- One-generation SU(2)^2 · U(1)_Y anomaly cancellation. -/
theorem anomaly_su2su2u1_cancels : anomalySU2SU2U1 = 0 := by
  native_decide

/-- One-generation U(1)_Y^3 anomaly cancellation. -/
theorem anomaly_u1cubed_cancels : anomalyU1Cubed = 0 := by
  native_decide

/-- One-generation grav^2 · U(1)_Y anomaly cancellation. -/
theorem anomaly_gravgravu1_cancels : anomalyGravGravU1 = 0 := by
  native_decide

/-- Number of SU(2) doublet copies (including color multiplicity). -/
def su2DoubletCopies : ℕ :=
  Finset.sum oneGeneration (fun f => if weakRep f = .doublet then colorMultiplicitySpecies f else 0)

/-- Explicit value of SU(2) doublet copies for one generation. -/
theorem su2_doublet_copies_eq_four : su2DoubletCopies = 4 := by
  native_decide

/-- Witten global SU(2) anomaly cancellation condition: even number of doublets. -/
theorem witten_global_su2_anomaly_cancelled : Even su2DoubletCopies := by
  rw [su2_doublet_copies_eq_four]
  decide

end Gutoe.SM.Anomalies
