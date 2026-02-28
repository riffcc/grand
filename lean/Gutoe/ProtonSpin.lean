import Mathlib
import Gutoe.DarkMatterSector
import Gutoe.GaugeConstants
import Gutoe.DimensionalStructure

namespace Gutoe.ProtonSpin

open Gutoe.DarkMatterSector
open Gutoe.GaugeConstants
open Gutoe.DimensionalStructure

/-- Quark-spin channel count from the finite dark-sector split:
    `|dark candidates| = 5`. -/
def quarkSpinChannelsQ : ℚ :=
  (darkSectorCandidates.card : ℚ)

/-- Gluon-spin channel count from SU(3):
    `3^2 - 1 = 8`. -/
def gluonSpinChannelsQ : ℚ :=
  ((3 ^ 2 - 1 : ℕ) : ℚ)

/-- Orbital-spin channel count from grade-1 directions:
    `|grade1_4d| = 4`. -/
def orbitalSpinChannelsQ : ℚ :=
  (grade1_4d.card : ℚ)

/-- Total proton-spin channel count in this structural lane:
    `5 + 8 + 4 = 17 = 16 + 1`. -/
def protonSpinChannelTotalQ : ℚ :=
  quarkSpinChannelsQ + gluonSpinChannelsQ + orbitalSpinChannelsQ

theorem quark_spin_channels_eq_5 :
    quarkSpinChannelsQ = 5 := by
  unfold quarkSpinChannelsQ
  rcases visible_dark_state_count_split with ⟨_, hDark, _, _⟩
  rw [hDark]
  norm_num

theorem gluon_spin_channels_eq_8 :
    gluonSpinChannelsQ = 8 := by
  unfold gluonSpinChannelsQ
  have h : (3 ^ 2 - 1 : ℕ) = 8 := su3_gluons
  rw [h]
  norm_num

theorem orbital_spin_channels_eq_4 :
    orbitalSpinChannelsQ = 4 := by
  unfold orbitalSpinChannelsQ
  rw [grade1_state_count_eq]
  norm_num

theorem proton_spin_channel_total_eq_17 :
    protonSpinChannelTotalQ = 17 := by
  unfold protonSpinChannelTotalQ
  rw [quark_spin_channels_eq_5, gluon_spin_channels_eq_8, orbital_spin_channels_eq_4]
  norm_num

theorem proton_spin_channel_total_eq_clifford_plus_identity :
    protonSpinChannelTotalQ = ((2 ^ 4 + 1 : ℕ) : ℚ) := by
  rw [proton_spin_channel_total_eq_17]
  norm_num

/-- Fraction of proton spin carried by quark helicity in this lane. -/
def quarkSpinFractionQ : ℚ :=
  quarkSpinChannelsQ / protonSpinChannelTotalQ

/-- Fraction of proton spin carried by gluon helicity in this lane. -/
def gluonSpinFractionQ : ℚ :=
  gluonSpinChannelsQ / protonSpinChannelTotalQ

/-- Fraction of proton spin carried by orbital angular momentum in this lane. -/
def orbitalSpinFractionQ : ℚ :=
  orbitalSpinChannelsQ / protonSpinChannelTotalQ

theorem quark_spin_fraction_eq_5_over_17 :
    quarkSpinFractionQ = 5 / 17 := by
  unfold quarkSpinFractionQ
  rw [quark_spin_channels_eq_5, proton_spin_channel_total_eq_17]

theorem gluon_spin_fraction_eq_8_over_17 :
    gluonSpinFractionQ = 8 / 17 := by
  unfold gluonSpinFractionQ
  rw [gluon_spin_channels_eq_8, proton_spin_channel_total_eq_17]

theorem orbital_spin_fraction_eq_4_over_17 :
    orbitalSpinFractionQ = 4 / 17 := by
  unfold orbitalSpinFractionQ
  rw [orbital_spin_channels_eq_4, proton_spin_channel_total_eq_17]

theorem proton_spin_fraction_partition :
    quarkSpinFractionQ + gluonSpinFractionQ + orbitalSpinFractionQ = 1 := by
  rw [quark_spin_fraction_eq_5_over_17, gluon_spin_fraction_eq_8_over_17,
    orbital_spin_fraction_eq_4_over_17]
  norm_num

/-- Proton total spin normalization (in units of `ℏ`). -/
def protonTotalSpinQ : ℚ := 1 / 2

/-- Quark angular-momentum share (`ℏ` units). -/
def quarkAngularMomentumQ : ℚ := protonTotalSpinQ * quarkSpinFractionQ

/-- Gluon angular-momentum share (`ℏ` units). -/
def gluonAngularMomentumQ : ℚ := protonTotalSpinQ * gluonSpinFractionQ

/-- Orbital angular-momentum share (`ℏ` units). -/
def orbitalAngularMomentumQ : ℚ := protonTotalSpinQ * orbitalSpinFractionQ

theorem angular_momentum_partition :
    quarkAngularMomentumQ + gluonAngularMomentumQ + orbitalAngularMomentumQ = protonTotalSpinQ := by
  unfold quarkAngularMomentumQ gluonAngularMomentumQ orbitalAngularMomentumQ protonTotalSpinQ
  rw [quark_spin_fraction_eq_5_over_17, gluon_spin_fraction_eq_8_over_17,
    orbital_spin_fraction_eq_4_over_17]
  norm_num

/-- First-pass phenomenology windows (broad ranges, no fit):
    quark 25–35%, gluon 35–55%, orbital 15–35%. -/
theorem proton_spin_fractions_within_broad_windows :
    (1 / 4 : ℚ) ≤ quarkSpinFractionQ ∧ quarkSpinFractionQ ≤ (7 / 20 : ℚ) ∧
    (7 / 20 : ℚ) ≤ gluonSpinFractionQ ∧ gluonSpinFractionQ ≤ (11 / 20 : ℚ) ∧
    (3 / 20 : ℚ) ≤ orbitalSpinFractionQ ∧ orbitalSpinFractionQ ≤ (7 / 20 : ℚ) := by
  rw [quark_spin_fraction_eq_5_over_17, gluon_spin_fraction_eq_8_over_17,
    orbital_spin_fraction_eq_4_over_17]
  native_decide

end Gutoe.ProtonSpin
