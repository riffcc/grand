import Mathlib
import Gutoe.DimensionalStructure
import Gutoe.DarkMatterSector

namespace Gutoe.AntimatterProduction

open Gutoe.DimensionalStructure
open Gutoe.DarkMatterSector

/-- Baseline chain efficiency used in the antimatter production lane:
    `10^-4 * 10^-3 * 10^-2 * 10^-1 * 10^0 = 10^-10`. -/
def antimatterChainBaselineQ : ℚ :=
  (1 / (10 ^ 4 : ℚ))
    * (1 / (10 ^ 3 : ℚ))
    * (1 / (10 ^ 2 : ℚ))
    * (1 / (10 ^ 1 : ℚ))
    * 1

theorem antimatter_chain_baseline_eq_1e10 :
    antimatterChainBaselineQ = 1 / (10 ^ 10 : ℚ) := by
  unfold antimatterChainBaselineQ
  norm_num

/-- Structural 2-OM upgrade multiplier from shared Cl(1,3) counts:
    `(dark states)^2 * (grade-1 states) = 5^2 * 4 = 100`. -/
def twoOmUpgradeMultiplierQ : ℚ :=
  (darkSectorCandidates.card : ℚ)
    * (darkSectorCandidates.card : ℚ)
    * (grade1_4d.card : ℚ)

theorem two_om_upgrade_multiplier_eq_100 :
    twoOmUpgradeMultiplierQ = 100 := by
  unfold twoOmUpgradeMultiplierQ
  rcases visible_dark_state_count_split with ⟨_, hDark, _, _⟩
  rw [hDark, grade1_state_count_eq]
  norm_num

/-- Upgraded chain efficiency after the structural 2-OM multiplier. -/
def antimatterChainUpgradedQ : ℚ :=
  antimatterChainBaselineQ * twoOmUpgradeMultiplierQ

theorem antimatter_chain_upgraded_eq_1e8 :
    antimatterChainUpgradedQ = 1 / (10 ^ 8 : ℚ) := by
  unfold antimatterChainUpgradedQ
  rw [antimatter_chain_baseline_eq_1e10, two_om_upgrade_multiplier_eq_100]
  norm_num

/-- Exact gain from the baseline to the upgraded chain is two orders of magnitude. -/
theorem antimatter_two_om_gain_exact :
    antimatterChainUpgradedQ / antimatterChainBaselineQ = 100 := by
  unfold antimatterChainUpgradedQ
  have hbase : antimatterChainBaselineQ ≠ 0 := by
    rw [antimatter_chain_baseline_eq_1e10]
    norm_num
  field_simp [hbase]
  exact two_om_upgrade_multiplier_eq_100

/-- Rest-output power proxy used by the antimatter comparison lane:
    `P_rest = η_net * P_beam`. -/
def antimatterRestOutputPowerQ (etaNetQ beamPowerQ : ℚ) : ℚ :=
  etaNetQ * beamPowerQ

/-- Annihilation-equivalent power proxy used for density comparison:
    `P_ann,eq = 2 * P_rest`. -/
def antimatterAnnihilationEquivalentPowerQ (etaNetQ beamPowerQ : ℚ) : ℚ :=
  2 * antimatterRestOutputPowerQ etaNetQ beamPowerQ

/-- By construction of the comparison lane, annihilation-equivalent power is
exactly two times the rest-output antimatter stream. -/
theorem antimatter_annihilation_equivalent_is_two_x
    (etaNetQ beamPowerQ : ℚ) :
    antimatterAnnihilationEquivalentPowerQ etaNetQ beamPowerQ =
      2 * antimatterRestOutputPowerQ etaNetQ beamPowerQ := by
  rfl

/-- For nonzero efficiency and beam power, the annihilation-equivalent to rest-output
power ratio is exactly 2. -/
theorem antimatter_annihilation_to_rest_ratio
    (etaNetQ beamPowerQ : ℚ)
    (heta : etaNetQ ≠ 0)
    (hbeam : beamPowerQ ≠ 0) :
    antimatterAnnihilationEquivalentPowerQ etaNetQ beamPowerQ
      / antimatterRestOutputPowerQ etaNetQ beamPowerQ = 2 := by
  unfold antimatterAnnihilationEquivalentPowerQ antimatterRestOutputPowerQ
  have hmul : etaNetQ * beamPowerQ ≠ 0 := mul_ne_zero heta hbeam
  field_simp [hmul]

/-- Specialization of the `2x` annihilation map to the structurally upgraded chain. -/
theorem upgraded_chain_annihilation_to_rest_ratio
    (beamPowerQ : ℚ) (hbeam : beamPowerQ ≠ 0) :
    antimatterAnnihilationEquivalentPowerQ antimatterChainUpgradedQ beamPowerQ
      / antimatterRestOutputPowerQ antimatterChainUpgradedQ beamPowerQ = 2 := by
  have hη : antimatterChainUpgradedQ ≠ 0 := by
    rw [antimatter_chain_upgraded_eq_1e8]
    norm_num
  exact antimatter_annihilation_to_rest_ratio antimatterChainUpgradedQ beamPowerQ hη hbeam

end Gutoe.AntimatterProduction
