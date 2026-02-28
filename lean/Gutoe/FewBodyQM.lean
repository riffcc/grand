import Mathlib
import Gutoe.NuclearFirstPrinciples

namespace Gutoe.FewBodyQM

open Gutoe.NuclearFirstPrinciples

/-- Light-nucleus few-body lane ceiling: A <= 16 from Cl(1,3) dimension. -/
def fewBodyMaxAQ : ℚ := cliffordDimQ

/-- Two-body attraction scale used by the few-body variational lane. -/
def fewBodyAttractiveDepthQ : ℚ := (3 * cliffordDimQ + gaugeTotalQ / 2) / (su2GeneratorsQ + u1GeneratorsQ)

/-- Short-range core repulsion scale for the few-body lane. -/
def fewBodyRepulsiveCoreQ : ℚ := semfAAQ

/-- Pair saturation scale from gauge-triplet ratio. -/
def fewBodyPairSaturationScaleQ : ℚ := gaugeTotalQ / (su2GeneratorsQ + u1GeneratorsQ)

/-- Tensor support scale from color generators. -/
def fewBodyTensorScaleQ : ℚ := su3GeneratorsQ / su2GeneratorsQ

/-- Compact three-body scale in the few-body lane. -/
def fewBodyThreeBodyScaleQ : ℚ := gaugeTotalQ

/-- Net depth margin (attraction minus repulsive core) for light-nucleus binding. -/
def fewBodyDepthMarginQ : ℚ := fewBodyAttractiveDepthQ - fewBodyRepulsiveCoreQ

theorem few_body_max_a_eq_16 :
    fewBodyMaxAQ = 16 := by
  unfold fewBodyMaxAQ
  exact clifford_dim_q_eq_16

theorem few_body_attractive_depth_eq_27_over_2 :
    fewBodyAttractiveDepthQ = 27 / 2 := by
  unfold fewBodyAttractiveDepthQ
  rw [clifford_dim_q_eq_16, gauge_total_q_eq_12, su2_generators_q_eq_3, u1GeneratorsQ]
  norm_num

theorem few_body_repulsive_core_eq_23 :
    fewBodyRepulsiveCoreQ = 23 := by
  unfold fewBodyRepulsiveCoreQ
  exact semf_a_a_eq_23

theorem few_body_pair_saturation_scale_eq_3 :
    fewBodyPairSaturationScaleQ = 3 := by
  unfold fewBodyPairSaturationScaleQ
  rw [gauge_total_q_eq_12, su2_generators_q_eq_3, u1GeneratorsQ]
  norm_num

theorem few_body_tensor_scale_eq_8_over_3 :
    fewBodyTensorScaleQ = 8 / 3 := by
  unfold fewBodyTensorScaleQ
  rw [su3_generators_q_eq_8, su2_generators_q_eq_3]

theorem few_body_three_body_scale_eq_12 :
    fewBodyThreeBodyScaleQ = 12 := by
  unfold fewBodyThreeBodyScaleQ
  exact gauge_total_q_eq_12

theorem few_body_depth_margin_eq_neg_19_over_2 :
    fewBodyDepthMarginQ = -19 / 2 := by
  unfold fewBodyDepthMarginQ
  rw [few_body_attractive_depth_eq_27_over_2, few_body_repulsive_core_eq_23]
  norm_num

theorem few_body_core_exceeds_pair_depth :
    fewBodyDepthMarginQ < 0 := by
  rw [few_body_depth_margin_eq_neg_19_over_2]
  norm_num

end Gutoe.FewBodyQM
