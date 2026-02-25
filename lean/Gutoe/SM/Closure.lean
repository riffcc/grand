/- 
 * GUTOE — Standard Model Closure Gate
 *
 * End-to-end consistency theorem for one generation:
 * reps + anomalies + Yukawas + normalization force the canonical registry.
-/

import Mathlib
import Gutoe.SM.Rep
import Gutoe.SM.Anomalies
import Gutoe.SM.HyperchargeBridge

namespace Gutoe.SM.Closure

open Gutoe.SM.Rep
open Gutoe.SM.Anomalies
open Gutoe.SM.HyperchargeBridge

/-- Canonical one-generation SM representation assignments. -/
def CanonicalRepConstraints : Prop :=
  colorRep .qL = .triplet ∧ weakRep .qL = .doublet ∧
  colorRep .uRc = .antiTriplet ∧ weakRep .uRc = .singlet ∧
  colorRep .dRc = .antiTriplet ∧ weakRep .dRc = .singlet ∧
  colorRep .lL = .singlet ∧ weakRep .lL = .doublet ∧
  colorRep .eRc = .singlet ∧ weakRep .eRc = .singlet ∧
  colorRep .nuRc = .singlet ∧ weakRep .nuRc = .singlet

/--
Full hypercharge constraints package `C` used for the one-generation closure gate.
Reps are tracked separately by `CanonicalRepConstraints`.
-/
def ConstraintsC (q u d ℓ e h : ℚ) : Prop :=
  u + d = -2 * q ∧
  ℓ = -3 * q ∧
  6 * q + 3 * u + 3 * d + 2 * ℓ + e = 0 ∧
  6 * q ^ 3 + 3 * u ^ 3 + 3 * d ^ 3 + 2 * ℓ ^ 3 + e ^ 3 = 0 ∧
  q + h + u = 0 ∧
  q - h + d = 0 ∧
  ℓ - h + e = 0 ∧
  q = hypercharge .qL

/-- Canonical hypercharges satisfy Yukawa gauge-invariance constraints. -/
theorem canonical_yukawa_gauge_invariance :
    hypercharge .qL + YH + hypercharge .uRc = 0 ∧
    hypercharge .qL - YH + hypercharge .dRc = 0 ∧
    hypercharge .lL - YH + hypercharge .eRc = 0 := by
  norm_num [hypercharge, YqL, YuRc, YdRc, YlL, YeRc, YH]

/-- Canonical registry satisfies representation constraints. -/
theorem canonical_rep_constraints_hold : CanonicalRepConstraints := by
  simp [CanonicalRepConstraints, colorRep, weakRep]

/-- Canonical registry satisfies full constraints package `C`. -/
theorem canonical_constraintsC_hold :
    ConstraintsC (hypercharge .qL) (hypercharge .uRc) (hypercharge .dRc)
      (hypercharge .lL) (hypercharge .eRc) YH := by
  rcases canonical_assignment_satisfies_constraints with
    ⟨hsu3, hsu2, hgrav, hu1cube, _hq_nonzero, _hqnorm, _hud⟩
  rcases canonical_yukawa_gauge_invariance with ⟨hyukU, hyukD, hyukE⟩
  refine ⟨hsu3, hsu2, hgrav, hu1cube, hyukU, hyukD, hyukE, rfl⟩

/--
Uniqueness: any one-generation hypercharges satisfying full package `C`
must be the canonical registry assignment.
-/
theorem constraintsC_force_canonical_registry
    {q u d ℓ e h : ℚ}
    (hC : ConstraintsC q u d ℓ e h) :
    u = hypercharge .uRc ∧ d = hypercharge .dRc ∧
      ℓ = hypercharge .lL ∧ e = hypercharge .eRc ∧ h = YH := by
  rcases hC with ⟨hsu3, hsu2, hgrav, hu1cube, hyukU, hyukD, hyukE, hqnorm⟩
  exact canonical_hypercharge_unique_under_constraints_C_full
    hsu3 hsu2 hgrav hu1cube hyukU hyukD hyukE hqnorm

/-- Electric charges `Q = T3 + Y` for left-chiral doublet components. -/
def QupL : ℚ := (1 / 2 : ℚ) + hypercharge .qL
def QdownL : ℚ := (-1 / 2 : ℚ) + hypercharge .qL
def QnuL : ℚ := (1 / 2 : ℚ) + hypercharge .lL
def QeL : ℚ := (-1 / 2 : ℚ) + hypercharge .lL

/-- Canonical electric charges of left-chiral doublet components. -/
theorem canonical_electric_charges :
    QupL = 2 / 3 ∧ QdownL = -1 / 3 ∧ QnuL = 0 ∧ QeL = -1 := by
  norm_num [QupL, QdownL, QnuL, QeL, hypercharge, YqL, YlL]

/--
End-to-end one-generation SM closure gate:
reps are canonical, all local anomalies cancel, Witten global anomaly cancels,
canonical assignment satisfies full `C`, electric charges are canonical, and
`C` uniquely forces the canonical charged-sector hypercharges + Higgs charge.
-/
theorem sm_consistency_complete :
    CanonicalRepConstraints ∧
    anomalySU3SU3U1 = 0 ∧
    anomalySU2SU2U1 = 0 ∧
    anomalyU1Cubed = 0 ∧
    anomalyGravGravU1 = 0 ∧
    Even su2DoubletCopies ∧
    ConstraintsC (hypercharge .qL) (hypercharge .uRc) (hypercharge .dRc)
      (hypercharge .lL) (hypercharge .eRc) YH ∧
    QupL = 2 / 3 ∧ QdownL = -1 / 3 ∧ QnuL = 0 ∧ QeL = -1 ∧
    (∀ q u d ℓ e h, ConstraintsC q u d ℓ e h →
      u = hypercharge .uRc ∧ d = hypercharge .dRc ∧
        ℓ = hypercharge .lL ∧ e = hypercharge .eRc ∧ h = YH) := by
  refine ⟨canonical_rep_constraints_hold,
    anomaly_su3su3u1_cancels,
    anomaly_su2su2u1_cancels,
    anomaly_u1cubed_cancels,
    anomaly_gravgravu1_cancels,
    witten_global_su2_anomaly_cancelled,
    canonical_constraintsC_hold,
    ?_⟩
  rcases canonical_electric_charges with ⟨hup, hdown, hnu, he⟩
  refine ⟨hup, hdown, hnu, he, ?_⟩
  intro q u d ℓ e h hC
  exact constraintsC_force_canonical_registry hC

end Gutoe.SM.Closure
