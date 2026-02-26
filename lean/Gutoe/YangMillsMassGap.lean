/-
 * GUTOE — Yang-Mills Mass Gap Validation Gates
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * Formal validation slice for GRAND-296/297/298:
 *   - Cl(1,3)→Z3 basis cardinality gate (3-state transfer basis)
 *   - positivity gate: λ₀ > λ₁ > 0 ⇒ m_gap > 0
 *   - concrete GEVP datapoint checks from `ym_mass_gap_report`
 *   - finite-volume monotone trend gate on reported gaps
 *
 * No `sorry`.
-/

import Mathlib
import Gutoe.Z3Uniqueness
import Gutoe.YangMillsStructuralGap

namespace Gutoe.YangMillsMassGap

open Real
open Gutoe.Z3Uniqueness
open Gutoe.YangMillsStructuralGap

/-- Transfer basis dimension is fixed by the Cl(1,3) Z₃ magnetic orbit count. -/
def transferBasisDim : ℕ := magneticTriplet.card

/-- Cl(1,3) forces a 3-state transfer basis via `magneticTriplet`. -/
theorem transfer_basis_dim_eq_three : transferBasisDim = 3 := by
  unfold transferBasisDim
  exact su2_dim

/-- Spectral-gap estimator from transfer-matrix eigenvalue ratio. -/
noncomputable def massGapFromEigenRatio (a_t lambda0 lambda1 : ℝ) : ℝ :=
  -(Real.log (lambda1 / lambda0)) / a_t

/-- Explicit Doeblin-induced lower bound on the mass gap. -/
noncomputable def doeblinGapLowerBound (a_t eps : ℝ) : ℝ :=
  -(Real.log (1 - eps)) / a_t

/-- Positivity gate: `0 < λ₁ < λ₀` and `a_t>0` implies positive mass gap. -/
theorem mass_gap_positive_of_eigen_ratio
    {a_t lambda0 lambda1 : ℝ}
    (ha : 0 < a_t)
    (h0 : 0 < lambda0)
    (h1 : 0 < lambda1)
    (h10 : lambda1 < lambda0) :
    0 < massGapFromEigenRatio a_t lambda0 lambda1 := by
  unfold massGapFromEigenRatio
  have hratio_pos : 0 < lambda1 / lambda0 := by
    exact div_pos h1 h0
  have hratio_lt_one : lambda1 / lambda0 < 1 := by
    have hne0 : lambda0 ≠ 0 := ne_of_gt h0
    field_simp [hne0]
    exact h10
  have hlog_lt : Real.log (lambda1 / lambda0) < Real.log 1 := by
    exact Real.log_lt_log hratio_pos hratio_lt_one
  have hlog_neg : Real.log (lambda1 / lambda0) < 0 := by
    simpa using hlog_lt
  have hnum : 0 < -Real.log (lambda1 / lambda0) := by
    linarith
  exact div_pos hnum ha

/-- If `λ₁/λ₀ ≤ 1-ε`, the mass gap is bounded below by the Doeblin expression. -/
theorem mass_gap_ge_doeblin_bound
    {a_t lambda0 lambda1 eps : ℝ}
    (ha : 0 < a_t)
    (h0 : 0 < lambda0)
    (h1 : 0 < lambda1)
    (heps : eps < 1)
    (hratio : lambda1 / lambda0 ≤ 1 - eps) :
    doeblinGapLowerBound a_t eps ≤ massGapFromEigenRatio a_t lambda0 lambda1 := by
  unfold doeblinGapLowerBound massGapFromEigenRatio
  have honeMinusPos : 0 < 1 - eps := by linarith
  have hratioPos : 0 < lambda1 / lambda0 := div_pos h1 h0
  have hlog_le : Real.log (lambda1 / lambda0) ≤ Real.log (1 - eps) :=
    Real.log_le_log hratioPos hratio
  have hneg : -Real.log (1 - eps) ≤ -Real.log (lambda1 / lambda0) := by linarith
  exact div_le_div_of_nonneg_right hneg (le_of_lt ha)

/-- The Doeblin lower bound is strictly positive when `0<ε<1`. -/
theorem doeblin_bound_positive
    {a_t eps : ℝ}
    (ha : 0 < a_t)
    (heps0 : 0 < eps)
    (heps1 : eps < 1) :
    0 < doeblinGapLowerBound a_t eps := by
  unfold doeblinGapLowerBound
  have honeMinusPos : 0 < 1 - eps := by linarith
  have honeMinusLtOne : 1 - eps < 1 := by linarith
  have hlogNeg : Real.log (1 - eps) < 0 := by
    have : Real.log (1 - eps) < Real.log 1 := Real.log_lt_log honeMinusPos honeMinusLtOne
    simpa using this
  have hnum : 0 < -Real.log (1 - eps) := by linarith
  exact div_pos hnum ha

/-- Explicit positive-gap corollary from a Doeblin-style subdominant ratio bound. -/
theorem mass_gap_positive_of_doeblin_ratio
    {a_t lambda0 lambda1 eps : ℝ}
    (ha : 0 < a_t)
    (h0 : 0 < lambda0)
    (h1 : 0 < lambda1)
    (heps0 : 0 < eps)
    (heps1 : eps < 1)
    (hratio : lambda1 / lambda0 ≤ 1 - eps) :
    0 < massGapFromEigenRatio a_t lambda0 lambda1 := by
  have hbound :
      doeblinGapLowerBound a_t eps ≤ massGapFromEigenRatio a_t lambda0 lambda1 :=
    mass_gap_ge_doeblin_bound ha h0 h1 heps1 hratio
  have hpos : 0 < doeblinGapLowerBound a_t eps := doeblin_bound_positive ha heps0 heps1
  exact lt_of_lt_of_le hpos hbound

/-- Structural bridge: Doeblin decomposition on a zero-sum eigenmode implies positive mass gap. -/
theorem mass_gap_positive_of_doeblin_mode
    {a_t eps lam mu : ℝ}
    (P R : Matrix (Fin 3) (Fin 3) ℝ)
    (v : Fin 3 → ℝ)
    (ha : 0 < a_t)
    (heps0 : 0 < eps)
    (heps1 : eps < 1)
    (hlam0 : 0 < lam)
    (hdecomp : ∀ i j, P i j = eps * uniformKernel i j + (1 - eps) * R i j)
    (hEigP : P.mulVec v = lam • v)
    (hEigR : R.mulVec v = mu • v)
    (hRstoch : R ∈ Matrix.rowStochastic ℝ (Fin 3))
    (hsum0 : v 0 + v 1 + v 2 = 0)
    (hvne : v ≠ 0) :
    0 < massGapFromEigenRatio a_t 1 lam := by
  have hlam_le : lam ≤ 1 - eps :=
    eigenvalue_le_one_sub_eps_of_decomposition_stochastic
      P R eps lam mu v heps1 hlam0.le hdecomp hEigP hEigR hRstoch hsum0 hvne
  have hratio : lam / (1 : ℝ) ≤ 1 - eps := by simpa using hlam_le
  exact mass_gap_positive_of_doeblin_ratio
    (a_t := a_t) (lambda0 := 1) (lambda1 := lam) (eps := eps)
    ha (by norm_num) hlam0 heps0 heps1 hratio

/-- GEVP eigenvalue estimates from `/tmp/bh_renders/ym_mass_gap_report`. -/
noncomputable def lambda0L6  : ℝ := (4942970855113019 : ℚ) / 5000000000000000
noncomputable def lambda1L6  : ℝ := (9712073993584441 : ℚ) / 10000000000000000
noncomputable def lambda0L8  : ℝ := (9910176117914437 : ℚ) / 10000000000000000
noncomputable def lambda1L8  : ℝ := (4868132483369487 : ℚ) / 5000000000000000
noncomputable def lambda0L10 : ℝ := (4971113258305147 : ℚ) / 5000000000000000
noncomputable def lambda1L10 : ℝ := (1955130885848443 : ℚ) / 2000000000000000
noncomputable def lambda0L12 : ℝ := (621596346623220 : ℚ) / 625000000000000
noncomputable def lambda1L12 : ℝ := (1958670042118565 : ℚ) / 2000000000000000

noncomputable def gapL6  : ℝ := massGapFromEigenRatio 1 lambda0L6  lambda1L6
noncomputable def gapL8  : ℝ := massGapFromEigenRatio 1 lambda0L8  lambda1L8
noncomputable def gapL10 : ℝ := massGapFromEigenRatio 1 lambda0L10 lambda1L10
noncomputable def gapL12 : ℝ := massGapFromEigenRatio 1 lambda0L12 lambda1L12

/-- Reported GEVP gap estimates from `ym_mass_gap_report` (Rust artifact values). -/
noncomputable def gapReportedL6  : ℝ := (887193280298489 : ℚ) / 50000000000000000
noncomputable def gapReportedL8  : ℝ := (885227472609799 : ℚ) / 50000000000000000
noncomputable def gapReportedL10 : ℝ := (16895938306527147 : ℚ) / 1000000000000000000
noncomputable def gapReportedL12 : ℝ := (30841519000448253 : ℚ) / 2000000000000000000

/-- Numerical ordering gate for reported GEVP eigenvalues. -/
theorem gevp_eigenvalue_ordering :
    0 < lambda1L6 ∧ lambda1L6 < lambda0L6 ∧
    0 < lambda1L8 ∧ lambda1L8 < lambda0L8 ∧
    0 < lambda1L10 ∧ lambda1L10 < lambda0L10 ∧
    0 < lambda1L12 ∧ lambda1L12 < lambda0L12 := by
  norm_num [lambda0L6, lambda1L6, lambda0L8, lambda1L8,
    lambda0L10, lambda1L10, lambda0L12, lambda1L12]

/-- Every reported finite-volume GEVP point has strictly positive gap. -/
theorem gevp_gap_positive_all_volumes :
    0 < gapL6 ∧ 0 < gapL8 ∧ 0 < gapL10 ∧ 0 < gapL12 := by
  rcases gevp_eigenvalue_ordering with
    ⟨h1p6, hlt6, h1p8, hlt8, h1p10, hlt10, h1p12, hlt12⟩
  constructor
  · exact mass_gap_positive_of_eigen_ratio (by norm_num) (by norm_num [lambda0L6]) h1p6 hlt6
  constructor
  · exact mass_gap_positive_of_eigen_ratio (by norm_num) (by norm_num [lambda0L8]) h1p8 hlt8
  constructor
  · exact mass_gap_positive_of_eigen_ratio (by norm_num) (by norm_num [lambda0L10]) h1p10 hlt10
  · exact mass_gap_positive_of_eigen_ratio (by norm_num) (by norm_num [lambda0L12]) h1p12 hlt12

/-- Reported finite-volume trend gate (GEVP lane): non-increasing with volume. -/
theorem gevp_gap_monotone_nonincreasing :
    gapReportedL8 ≤ gapReportedL6 ∧
    gapReportedL10 ≤ gapReportedL8 ∧
    gapReportedL12 ≤ gapReportedL10 := by
  norm_num [gapReportedL6, gapReportedL8, gapReportedL10, gapReportedL12]

end Gutoe.YangMillsMassGap
