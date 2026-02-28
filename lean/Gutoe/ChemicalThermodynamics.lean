import Mathlib

namespace Gutoe.ChemicalThermodynamics

/-- Fusion entropy scale (Richards-like rule proxy) in J/(mol*K) units. -/
def fusionEntropyQ : ℚ := 10

/-- Vaporization entropy scale (Trouton-like rule proxy) in J/(mol*K) units. -/
def vaporEntropyQ : ℚ := 85

/-- Reference pressure anchor used by the runtime lane (1 atm, Pa units). -/
def pressureRefQ : ℚ := 101325

/-- Transition-family latent fusion fraction of cohesive energy in the runtime lane. -/
def transitionFusionFractionQ : ℚ := 7 / 200

/-- Transition-family latent vaporization fraction of cohesive energy in the runtime lane. -/
def transitionVaporFractionQ : ℚ := 41 / 50

/-- Reduced melting scale factor (cohesive multiplier divided by entropy scale). -/
def reducedMeltingScaleQ : ℚ := transitionFusionFractionQ / fusionEntropyQ

/-- Reduced boiling scale factor (cohesive multiplier divided by entropy scale). -/
def reducedBoilingScaleQ : ℚ := transitionVaporFractionQ / vaporEntropyQ

/-- Reduced liquid offset for phase selection at temperature ratio `T/Tm`. -/
def reducedLiquidOffsetQ (t tm : ℚ) : ℚ := transitionFusionFractionQ * (1 - t / tm)

/-- Reduced vapor offset (before pressure term) at temperature ratio `T/Tb`. -/
def reducedVaporOffsetQ (t tb : ℚ) : ℚ := transitionVaporFractionQ * (1 - t / tb)

theorem transition_fusion_fraction_eq :
    transitionFusionFractionQ = 7 / 200 := by
  rfl

theorem transition_vapor_fraction_eq :
    transitionVaporFractionQ = 41 / 50 := by
  rfl

theorem fusion_entropy_eq :
    fusionEntropyQ = 10 := by
  rfl

theorem vapor_entropy_eq :
    vaporEntropyQ = 85 := by
  rfl

theorem pressure_ref_eq :
    pressureRefQ = 101325 := by
  rfl

theorem reduced_melting_scale_eq :
    reducedMeltingScaleQ = 7 / 2000 := by
  unfold reducedMeltingScaleQ transitionFusionFractionQ fusionEntropyQ
  norm_num

theorem reduced_boiling_scale_eq :
    reducedBoilingScaleQ = 41 / 4250 := by
  unfold reducedBoilingScaleQ transitionVaporFractionQ vaporEntropyQ
  norm_num

theorem reduced_boiling_scale_gt_melting_scale :
    reducedBoilingScaleQ > reducedMeltingScaleQ := by
  rw [reduced_boiling_scale_eq, reduced_melting_scale_eq]
  norm_num

/-- The reduced vapor offset is zero at the boiling-point anchor. -/
theorem reduced_vapor_offset_at_boiling (tb : ℚ) (htb : tb ≠ 0) :
    reducedVaporOffsetQ tb tb = 0 := by
  unfold reducedVaporOffsetQ transitionVaporFractionQ
  have hdiv : tb / tb = 1 := by exact div_self htb
  rw [hdiv]
  ring

/-- The reduced liquid offset is zero at the melting-point anchor. -/
theorem reduced_liquid_offset_at_melting (tm : ℚ) (htm : tm ≠ 0) :
    reducedLiquidOffsetQ tm tm = 0 := by
  unfold reducedLiquidOffsetQ transitionFusionFractionQ
  have hdiv : tm / tm = 1 := by exact div_self htm
  rw [hdiv]
  ring

/-- For positive cohesive energy, transition-family boiling scale exceeds melting scale
    under the lane's fixed latent-fraction and entropy constants. -/
theorem transition_boiling_scale_above_melting
    {cohesive : ℚ} (hcoh : 0 < cohesive) :
    cohesive * reducedBoilingScaleQ > cohesive * reducedMeltingScaleQ := by
  have hscale : reducedBoilingScaleQ > reducedMeltingScaleQ :=
    reduced_boiling_scale_gt_melting_scale
  nlinarith [hcoh, hscale]

end Gutoe.ChemicalThermodynamics
