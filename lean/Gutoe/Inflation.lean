import Mathlib
import Gutoe.DarkMatterSector
import Gutoe.GravityMetric

namespace Gutoe.Inflation

open Gutoe.DarkMatterSector
open Gutoe.GravityMetric

/-- Structural inflation e-fold count from shared Clifford dark-sector counts:
    `N = 12 * 5 = 60`. -/
noncomputable def inflationEfoldCount : ℝ :=
  (geometricDarkAmplificationQ : ℝ) * (darkSectorCandidates.card : ℝ)

/-- Slow-roll epsilon in the structural plateau lane. -/
noncomputable def slowRollEpsilon (N : ℝ) : ℝ :=
  3 / (4 * N ^ 2)

/-- Slow-roll eta in the structural plateau lane. -/
noncomputable def slowRollEta (N : ℝ) : ℝ :=
  -1 / N

/-- Scalar spectral index from first-order slow-roll observables. -/
noncomputable def scalarSpectralIndex (N : ℝ) : ℝ :=
  1 - 6 * slowRollEpsilon N + 2 * slowRollEta N

/-- Tensor-to-scalar ratio from epsilon. -/
noncomputable def tensorToScalarRatio (N : ℝ) : ℝ :=
  16 * slowRollEpsilon N

/-- End-of-inflation proxy where `ε = 1` in this lane. -/
noncomputable def inflationEndN : ℝ :=
  Real.sqrt 3 / 2

/-- Total expansion factor `exp(N)`. -/
noncomputable def inflationExpansionFactor : ℝ :=
  Real.exp inflationEfoldCount

theorem dark_sector_card_eq_five :
    darkSectorCandidates.card = 5 := by
  rcases visible_dark_state_count_split with ⟨_, hDark, _, _⟩
  exact hDark

theorem inflation_efolds_eq_60 :
    inflationEfoldCount = 60 := by
  unfold inflationEfoldCount
  rw [geometric_dark_amplification_eq, dark_sector_card_eq_five]
  norm_num

theorem epsilon_structural_eq :
    slowRollEpsilon inflationEfoldCount = (1 : ℝ) / 4800 := by
  rw [inflation_efolds_eq_60]
  norm_num [slowRollEpsilon]

theorem eta_structural_eq :
    slowRollEta inflationEfoldCount = -(1 : ℝ) / 60 := by
  rw [inflation_efolds_eq_60]
  norm_num [slowRollEta]

theorem ns_structural_eq :
    scalarSpectralIndex inflationEfoldCount = (2317 : ℝ) / 2400 := by
  rw [inflation_efolds_eq_60]
  norm_num [scalarSpectralIndex, slowRollEpsilon, slowRollEta]

theorem r_structural_eq :
    tensorToScalarRatio inflationEfoldCount = (1 : ℝ) / 300 := by
  rw [inflation_efolds_eq_60]
  norm_num [tensorToScalarRatio, slowRollEpsilon]

theorem ns_in_observational_window :
    (0.955 : ℝ) ≤ scalarSpectralIndex inflationEfoldCount ∧
    scalarSpectralIndex inflationEfoldCount ≤ 0.975 := by
  rw [ns_structural_eq]
  constructor <;> norm_num

theorem r_below_current_upper_bound :
    tensorToScalarRatio inflationEfoldCount < 0.06 := by
  rw [r_structural_eq]
  norm_num

theorem graceful_exit_condition :
    slowRollEpsilon inflationEndN = 1 := by
  unfold slowRollEpsilon inflationEndN
  have hsq : (Real.sqrt 3) ^ 2 = 3 := by
    have h3 : (0 : ℝ) ≤ 3 := by positivity
    simpa using Real.sq_sqrt h3
  calc
    3 / (4 * (Real.sqrt 3 / 2) ^ 2)
        = 3 / (4 * ((Real.sqrt 3) ^ 2 / 4)) := by ring
    _ = 3 / ((Real.sqrt 3) ^ 2) := by ring
    _ = 3 / 3 := by rw [hsq]
    _ = 1 := by norm_num

theorem expansion_factor_gt_zero :
    0 < inflationExpansionFactor := by
  unfold inflationExpansionFactor
  positivity

end Gutoe.Inflation
