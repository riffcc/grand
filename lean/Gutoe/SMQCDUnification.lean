/- 
 * GUTOE — SM × QCD Unification Gate
 *
 * Formal bridge that combines:
 *   - one-generation Standard Model closure,
 *   - QCD structural sector gates (color, beta sign, strong-CP structural pin),
 *   - and the two-case Strong-CP general split.
 -/

import Mathlib
import Gutoe.SM.Closure
import Gutoe.GaugeConstants
import Gutoe.GaugeGroupSU3
import Gutoe.AsymptoticFreedomEntropy
import Gutoe.StrongCP
import Gutoe.StrongCPGeneralCases
import Gutoe.ChiralSymmetryBreaking

namespace Gutoe.SMQCDUnification

open Gutoe.SM.Closure
open Gutoe.SM.Anomalies
open Gutoe.GaugeConstants
open Gutoe.GaugeGroupSU3
open Gutoe.AsymptoticFreedomEntropy
open Gutoe.StrongCP
open Gutoe.StrongCPGeneralCases
open Gutoe.StrongCPEmergence
open Gutoe.StrongCPVacuum
open Gutoe.ChiralSymmetryBreaking

/-- QCD structural core gate derived from the Cl(1,3) chain. -/
def qcdCoreGate : Prop :=
  quarkOrbit.card = 3 ∧
  quarkOrbit.card ^ 2 - 1 = 8 ∧
  0 < beta0Clifford ∧
  thetaQcdStructural = 0

/-- QCD structural core gate is satisfied by existing theorem chain. -/
theorem qcd_core_gate_holds : qcdCoreGate := by
  refine ⟨quarkOrbit_card, quarks_predict_gluon_count, beta0_clifford_pos, theta_qcd_structural_zero⟩

/-- Unified SM × QCD structural bundle (algebraic closure + QCD gates). -/
def smQcdUnifiedStructural : Prop :=
  CanonicalRepConstraints ∧
  anomalySU3SU3U1 = 0 ∧
  anomalySU2SU2U1 = 0 ∧
  anomalyU1Cubed = 0 ∧
  anomalyGravGravU1 = 0 ∧
  Even su2DoubletCopies ∧
  qcdCoreGate ∧
  ((3 ^ 2 - 1) + (2 ^ 2 - 1) + 1 = 12)

/-- End-to-end structural unification theorem for SM × QCD. -/
theorem sm_qcd_unified_structural_holds : smQcdUnifiedStructural := by
  rcases sm_consistency_complete with
    ⟨hrep, hsu3, hsu2, hu1, hgrav, hwitten, _hC, _hQup, _hQdown, _hQnu, _hQe, _huniq⟩
  exact ⟨hrep, hsu3, hsu2, hu1, hgrav, hwitten, qcd_core_gate_holds, total_gauge_bosons⟩

/-- Full unification bundle including the Strong-CP two-case split:
    (A) physical GUTOE emergent image => `θ` unphysical,
    (B) any nonzero topological sector => `θ` physical. -/
def smQcdGeneralCaseBundle
    {X : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    (x0 : X)
    (qClass : HomotopyClass X Su3Matrix → ℤ) : Prop :=
  smQcdUnifiedStructural ∧
  (∀ (f : C(X, FundamentalGaugeGroup)) (theta1 theta2 : ℝ),
    thetaPhase theta1 (qEffFromClassAnchored x0 qClass (z3ToSu3.comp f)) =
    thetaPhase theta2 (qEffFromClassAnchored x0 qClass (z3ToSu3.comp f))) ∧
  (∀ q : ℤ, q ≠ 0 → ∃ theta1 theta2 : ℝ,
    thetaPhase theta1 q ≠ thetaPhase theta2 q)

/-- Master theorem: SM closure, QCD structural gate, and Strong-CP general split
    coexist in one formal package. -/
theorem sm_qcd_general_case_bundle_holds
    {X : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    (x0 : X)
    (qClass : HomotopyClass X Su3Matrix → ℤ) :
    smQcdGeneralCaseBundle x0 qClass := by
  have hsplit := strong_cp_general_case_split x0 qClass
  exact ⟨sm_qcd_unified_structural_holds, hsplit.1, hsplit.2⟩

/-- Unified QCD structural closure including the GRAND-126 chiral gate. -/
def qcdCoreWithChiralGate : Prop :=
  qcdCoreGate ∧
  quarkCondensateProxy < 0 ∧
  0 < pionMassSqProxy ∧
  pseudoGoldstoneRatio = (1 : ℝ) / 137 ∧
  pionMassSqFromExplicitBreaking 0 = 0 ∧
  0 < confinementChiralLinkStrength

/-- Existing QCD core gate and chiral-symmetry gate hold simultaneously. -/
theorem qcd_core_with_chiral_gate_holds : qcdCoreWithChiralGate := by
  exact ⟨qcd_core_gate_holds, chiral_symmetry_breaking_gate⟩

end Gutoe.SMQCDUnification
