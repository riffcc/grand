import Mathlib
import Gutoe.FineStructure
import Gutoe.DarkMatterSector
import Gutoe.AbiogenesisThreshold

namespace Gutoe.LifeProgressionEntropy

open Gutoe.FineStructure
open Gutoe.DarkMatterSector
open Gutoe.AbiogenesisThreshold

/-- Bare-rock dissipative baseline in normalized units. -/
def stage0BareRockQ : ℚ := 1

/-- Prebiotic incremental gain from the α lane: `11/137`. -/
def prebioticGainQ : ℚ := 11 / (alphaInverse 4 : ℚ)

theorem prebiotic_gain_eq :
    prebioticGainQ = 11 / 137 := by
  unfold prebioticGainQ
  rw [alpha_inverse_d4]
  norm_num

/-- Autocatalytic incremental gain from the abiogenesis closure margin. -/
def autocatalyticGainQ : ℚ := abiogenesisClosureMarginQ

theorem autocatalytic_gain_eq :
    autocatalyticGainQ = 3473 / 9727 := by
  unfold autocatalyticGainQ
  exact abiogenesis_closure_margin_eq

/-- Photosynthetic incremental gain from geometric dark fraction: `60/71`. -/
def photosyntheticGainQ : ℚ := geometricDarkFractionOfMatter

theorem photosynthetic_gain_eq :
    photosyntheticGainQ = 60 / 71 := by
  unfold photosyntheticGainQ
  exact geometric_dark_fraction_of_matter_eq

/-- Multicellular incremental gain from finite dark/visible split: `5/11`. -/
def multicellularGainQ : ℚ := darkToVisibleCountRatio

theorem multicellular_gain_eq :
    multicellularGainQ = 5 / 11 := by
  unfold multicellularGainQ
  exact dark_to_visible_count_ratio_eq

/-- Intelligence incremental gain from geometric dark/visible ratio: `60/11`. -/
def intelligenceGainQ : ℚ := geometricDarkToVisibleRatio

theorem intelligence_gain_eq :
    intelligenceGainQ = 60 / 11 := by
  unfold intelligenceGainQ
  exact geometric_dark_to_visible_ratio_eq

/-- Stage totals in normalized entropy-production units. -/
def stage1PrebioticQ : ℚ := stage0BareRockQ + prebioticGainQ
def stage2AutocatalyticQ : ℚ := stage1PrebioticQ + autocatalyticGainQ
def stage3PhotosyntheticQ : ℚ := stage2AutocatalyticQ + photosyntheticGainQ
def stage4MulticellularQ : ℚ := stage3PhotosyntheticQ + multicellularGainQ
def stage5IntelligenceQ : ℚ := stage4MulticellularQ + intelligenceGainQ

theorem stage_totals_closed_form :
    stage0BareRockQ = 1 ∧
    stage1PrebioticQ = 148 / 137 ∧
    stage2AutocatalyticQ = 13981 / 9727 ∧
    stage3PhotosyntheticQ = 22201 / 9727 ∧
    stage4MulticellularQ = 292846 / 106997 ∧
    stage5IntelligenceQ = 876466 / 106997 := by
  have h0 : stage0BareRockQ = 1 := by
    norm_num [stage0BareRockQ]
  have h1 : stage1PrebioticQ = 148 / 137 := by
    unfold stage1PrebioticQ stage0BareRockQ
    rw [prebiotic_gain_eq]
    norm_num
  have h2 : stage2AutocatalyticQ = 13981 / 9727 := by
    unfold stage2AutocatalyticQ
    rw [h1, autocatalytic_gain_eq]
    norm_num
  have h3 : stage3PhotosyntheticQ = 22201 / 9727 := by
    unfold stage3PhotosyntheticQ
    rw [h2, photosynthetic_gain_eq]
    norm_num
  have h4 : stage4MulticellularQ = 292846 / 106997 := by
    unfold stage4MulticellularQ
    rw [h3, multicellular_gain_eq]
    norm_num
  have h5 : stage5IntelligenceQ = 876466 / 106997 := by
    unfold stage5IntelligenceQ
    rw [h4, intelligence_gain_eq]
    norm_num
  constructor
  · exact h0
  constructor
  · exact h1
  constructor
  · exact h2
  constructor
  · exact h3
  constructor
  · exact h4
  · exact h5

/-- Strict progression: each stage raises normalized entropy-production capacity. -/
theorem strict_entropy_progression :
    stage0BareRockQ < stage1PrebioticQ ∧
    stage1PrebioticQ < stage2AutocatalyticQ ∧
    stage2AutocatalyticQ < stage3PhotosyntheticQ ∧
    stage3PhotosyntheticQ < stage4MulticellularQ ∧
    stage4MulticellularQ < stage5IntelligenceQ := by
  constructor
  · unfold stage1PrebioticQ stage0BareRockQ
    rw [prebiotic_gain_eq]
    norm_num
  constructor
  · unfold stage2AutocatalyticQ stage1PrebioticQ
    rw [autocatalytic_gain_eq]
    norm_num [prebiotic_gain_eq]
  constructor
  · unfold stage3PhotosyntheticQ stage2AutocatalyticQ
    rw [photosynthetic_gain_eq]
    norm_num [autocatalytic_gain_eq, prebiotic_gain_eq]
  constructor
  · unfold stage4MulticellularQ stage3PhotosyntheticQ
    rw [multicellular_gain_eq]
    norm_num [photosynthetic_gain_eq, autocatalytic_gain_eq, prebiotic_gain_eq]
  · unfold stage5IntelligenceQ stage4MulticellularQ
    rw [intelligence_gain_eq]
    norm_num [multicellular_gain_eq, photosynthetic_gain_eq, autocatalytic_gain_eq, prebiotic_gain_eq]

/-- Intelligence step dominates all prior incremental gains. -/
theorem intelligence_step_dominates :
    intelligenceGainQ > photosyntheticGainQ ∧
    intelligenceGainQ > multicellularGainQ ∧
    intelligenceGainQ > autocatalyticGainQ ∧
    intelligenceGainQ > prebioticGainQ := by
  constructor
  · rw [intelligence_gain_eq, photosynthetic_gain_eq]
    norm_num
  constructor
  · rw [intelligence_gain_eq, multicellular_gain_eq]
    norm_num
  constructor
  · rw [intelligence_gain_eq, autocatalytic_gain_eq]
    norm_num
  · rw [intelligence_gain_eq, prebiotic_gain_eq]
    norm_num

/-- Composite formal gate for the progression lane. -/
theorem life_progression_entropy_gate :
    stage0BareRockQ < stage1PrebioticQ ∧
    stage1PrebioticQ < stage2AutocatalyticQ ∧
    stage2AutocatalyticQ < stage3PhotosyntheticQ ∧
    stage3PhotosyntheticQ < stage4MulticellularQ ∧
    stage4MulticellularQ < stage5IntelligenceQ ∧
    intelligenceGainQ > photosyntheticGainQ ∧
    intelligenceGainQ > multicellularGainQ ∧
    intelligenceGainQ > autocatalyticGainQ ∧
    intelligenceGainQ > prebioticGainQ := by
  rcases strict_entropy_progression with ⟨h01, h12, h23, h34, h45⟩
  rcases intelligence_step_dominates with ⟨hi3, hi4, hi2, hi1⟩
  exact ⟨h01, h12, h23, h34, h45, hi3, hi4, hi2, hi1⟩

end Gutoe.LifeProgressionEntropy
