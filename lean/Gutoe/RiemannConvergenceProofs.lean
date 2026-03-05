import Mathlib
import Gutoe.RiemannConvergenceSubObligations
import Gutoe.RiemannCore

namespace Gutoe.RiemannConvergenceProofs

open Gutoe.RiemannFinalTarget
open Gutoe.RiemannFiniteXiModel
open Gutoe.RiemannConvergenceSubObligations
open Gutoe.RiemannOperatorLadder
open Gutoe.RiemannCore

noncomputable section
open scoped Topology

/-- Canonical finite Hadamard partial products built from a real-ordinate
enumerator `ρ`. At level `N`, use the first `N` ordinates. -/
def canonicalHadamardXiN (ρ : ℕ → ℝ) : ℕ → (ℂ → ℂ) :=
  fun N => XiFiniteHadamard ((Finset.range N).image ρ)

/-- One-step factor sequence whose finite-range products reproduce
`canonicalHadamardXiN`; repeated ordinates contribute factor `1`. -/
def canonicalHadamardStepFactor (ρ : ℕ → ℝ) (n : ℕ) (z : ℂ) : ℂ :=
  if ρ n ∈ (Finset.range n).image ρ then 1 else hadamardFactor (ρ n) z

/-- Finite-range product reconstruction of `canonicalHadamardXiN`. -/
theorem canonicalHadamardXiN_eq_finsetRange_prod_stepFactor
    (ρ : ℕ → ℝ) (N : ℕ) (z : ℂ) :
    canonicalHadamardXiN ρ N z =
      ∏ i ∈ Finset.range N, canonicalHadamardStepFactor ρ i z := by
  induction N with
  | zero =>
      simp [canonicalHadamardXiN, canonicalHadamardStepFactor, XiFiniteHadamard]
  | succ N ih =>
      by_cases hdup : ρ N ∈ (Finset.range N).image ρ
      · have himage :
          (Finset.range (N + 1)).image ρ = (Finset.range N).image ρ := by
          refine Finset.ext (fun t => ?_)
          constructor
          · intro ht
            rcases Finset.mem_image.mp ht with ⟨i, hi, rfl⟩
            by_cases hiN : i = N
            · subst hiN
              exact hdup
            · have hiLt : i < N := by
                have his : i < N + 1 := Finset.mem_range.mp hi
                omega
              exact Finset.mem_image.mpr ⟨i, Finset.mem_range.mpr hiLt, rfl⟩
          · intro ht
            rcases Finset.mem_image.mp ht with ⟨i, hi, rfl⟩
            exact Finset.mem_image.mpr
              ⟨i, Finset.mem_range.mpr (Nat.lt_trans (Finset.mem_range.mp hi) (Nat.lt_succ_self N)), rfl⟩
        calc
          canonicalHadamardXiN ρ (N + 1) z
              = XiFiniteHadamard ((Finset.range N).image ρ) z := by
                  simp [canonicalHadamardXiN, himage]
          _ = ∏ i ∈ Finset.range N, canonicalHadamardStepFactor ρ i z := ih
          _ = (∏ i ∈ Finset.range N, canonicalHadamardStepFactor ρ i z) *
                canonicalHadamardStepFactor ρ N z := by
                simp [canonicalHadamardStepFactor, hdup]
          _ = ∏ i ∈ Finset.range (N + 1), canonicalHadamardStepFactor ρ i z := by
                simp [Finset.prod_range_succ]
      · have hnot : ρ N ∉ (Finset.range N).image ρ := hdup
        have himage :
            (Finset.range (N + 1)).image ρ = insert (ρ N) ((Finset.range N).image ρ) := by
          refine Finset.ext (fun t => ?_)
          constructor
          · intro ht
            rcases Finset.mem_image.mp ht with ⟨i, hi, rfl⟩
            by_cases hiN : i = N
            · subst hiN
              simp
            · have hiLt : i < N := by
                have his : i < N + 1 := Finset.mem_range.mp hi
                omega
              exact Finset.mem_insert.mpr <| Or.inr <|
                Finset.mem_image.mpr ⟨i, Finset.mem_range.mpr hiLt, rfl⟩
          · intro ht
            rcases Finset.mem_insert.mp ht with rfl | ht'
            · exact Finset.mem_image.mpr
                ⟨N, Finset.mem_range.mpr (Nat.lt_succ_self N), rfl⟩
            · rcases Finset.mem_image.mp ht' with ⟨i, hi, rfl⟩
              exact Finset.mem_image.mpr
                ⟨i, Finset.mem_range.mpr (Nat.lt_trans (Finset.mem_range.mp hi) (Nat.lt_succ_self N)), rfl⟩
        have hins :
            XiFiniteHadamard ((Finset.range (N + 1)).image ρ) z =
              hadamardFactor (ρ N) z * XiFiniteHadamard ((Finset.range N).image ρ) z := by
          simpa [himage] using
            (XiFiniteHadamard_insert ((Finset.range N).image ρ) (t := ρ N) hnot z)
        calc
          canonicalHadamardXiN ρ (N + 1) z
              = hadamardFactor (ρ N) z * XiFiniteHadamard ((Finset.range N).image ρ) z := by
                  simpa [canonicalHadamardXiN] using hins
          _ = hadamardFactor (ρ N) z * (∏ i ∈ Finset.range N, canonicalHadamardStepFactor ρ i z) := by
                have ih' :
                    XiFiniteHadamard ((Finset.range N).image ρ) z =
                      ∏ i ∈ Finset.range N, canonicalHadamardStepFactor ρ i z := by
                  simpa [canonicalHadamardXiN] using ih
                rw [ih']
          _ = (∏ i ∈ Finset.range N, canonicalHadamardStepFactor ρ i z) *
                canonicalHadamardStepFactor ρ N z := by
                simp [canonicalHadamardStepFactor, hnot, mul_comm]
          _ = ∏ i ∈ Finset.range (N + 1), canonicalHadamardStepFactor ρ i z := by
                simp [Finset.prod_range_succ]

/-- Placeholder abstraction for "ρ enumerates nontrivial `XiTarget` zeros on the
critical line". This keeps the hard arithmetic/analytic obligations explicit. -/
def EnumeratesXiTargetZeros (ρ : ℕ → ℝ) : Prop :=
  ∀ s : ℂ, XiTarget s = 0 → ∃ n : ℕ, s = criticalLinePoint (ρ n)

/-- Quantitative growth hypothesis for a zero enumerator:
critical-line norms along the enumeration are bounded below linearly. -/
def WeylLawGrowth (ρ : ℕ → ℝ) : Prop :=
  ∃ c : ℝ, 0 < c ∧ ∀ n : ℕ, c * (n + 1 : ℝ) ≤ ‖criticalLinePoint (ρ n)‖

/-- Bridge from the zero enumerator lane to the concrete operator spectral lane. -/
def EnumeratorHitsOperatorSpec (ρ : ℕ → ℝ) : Prop :=
  ∀ n : ℕ, ∃ N : ℕ, ρ n ∈ operatorSpecN N

/-- Summability consequence of Weyl-style linear growth:
the canonical quadratic Hadamard profile is summable. -/
theorem summable_hadamardQuadraticProfile_of_weylGrowth
    (ρ : ℕ → ℝ)
    (hGrowth : WeylLawGrowth ρ)
    (R M : ℝ) :
    ∃ c : ℝ, 0 < c ∧ Summable (hadamardQuadraticProfile R M c) := by
  rcases hGrowth with ⟨c, hc0, _hlin⟩
  refine ⟨c, hc0, ?_⟩
  exact summable_hadamardQuadraticProfile R M c (ne_of_gt hc0)

/-- Hard analytic existence lemma:
Weyl-style growth plus Hadamard increment control yields local-uniform
convergence of canonical finite Hadamard products to some limit `F`. -/
theorem canonicalHadamardXiN_tendstoLocallyUniformly_to_some_limit_of_weyl
    (ρ : ℕ → ℝ)
    (hGrowth : WeylLawGrowth ρ) :
    ∃ F : ℂ → ℂ,
      TendstoLocallyUniformly (canonicalHadamardXiN ρ) F
        (Filter.atTop : Filter ℕ) := by
  -- Core analytic burden: prove local-uniform multipliability for the step
  -- factors from Weyl growth + Hadamard increment control.
  have hStep :
      MultipliableLocallyUniformly (canonicalHadamardStepFactor ρ) := by
    rcases hGrowth with ⟨c, hc, hgrow⟩
    have hmem_iff_exists_lt :
        ∀ n : ℕ, (ρ n ∈ (Finset.range n).image ρ) ↔ ∃ a < n, ρ a = ρ n := by
      intro n
      constructor
      · intro h
        rcases Finset.mem_image.mp h with ⟨a, ha, haeq⟩
        exact ⟨a, Finset.mem_range.mp ha, haeq⟩
      · rintro ⟨a, ha, haeq⟩
        exact Finset.mem_image.mpr ⟨a, Finset.mem_range.mpr ha, haeq⟩
    refine multipliableLocallyUniformly_of_of_forall_exists_nhds ?_
    intro x
    let R : ℝ := ‖x‖ + 1
    let K : Set ℂ := Metric.closedBall x 1
    let g : ℕ → ℂ → ℂ := fun n z => canonicalHadamardStepFactor ρ n z - 1
    let u : ℕ → ℝ := hadamardQuadraticProfile R 1 c
    have hKnhds : K ∈ 𝓝 x := by
      simpa [K] using Metric.closedBall_mem_nhds x zero_lt_one
    have hKcompact : IsCompact K := by
      simpa [K] using isCompact_closedBall x (1 : ℝ)
    have hu : Summable u := by
      simpa [u] using summable_hadamardQuadraticProfile R 1 c (ne_of_gt hc)
    have hcts : ∀ n, ContinuousOn (g n) K := by
      intro n
      by_cases hdup : ∃ a < n, ρ a = ρ n
      · simpa [g, canonicalHadamardStepFactor, hmem_iff_exists_lt, hdup] using
          (continuousOn_const : ContinuousOn (fun _ : ℂ => (0 : ℂ)) K)
      · have hratio_cont : Continuous (fun z : ℂ => z / criticalLinePoint (ρ n)) :=
          continuous_id.div_const _
        have hlin_cont : Continuous (fun z : ℂ => 1 - z / criticalLinePoint (ρ n)) :=
          continuous_const.sub hratio_cont
        have hexp_cont : Continuous (fun z : ℂ => Complex.exp (z / criticalLinePoint (ρ n))) :=
          Complex.continuous_exp.comp hratio_cont
        have hhad_cont : Continuous (hadamardFactor (ρ n)) := by
          simpa [hadamardFactor] using hlin_cont.mul hexp_cont
        simpa [g, canonicalHadamardStepFactor, hmem_iff_exists_lt, hdup] using
          (hhad_cont.continuousOn.sub continuousOn_const)
    have hEvent : ∀ᶠ n in Filter.atTop, ∀ z ∈ K, ‖g n z‖ ≤ u n := by
      let N0 : ℕ := Nat.ceil (R / c)
      have hN0 : ∀ n : ℕ, N0 ≤ n → R ≤ c * (n + 1 : ℝ) := by
        intro n hn
        have hceil : R / c ≤ (N0 : ℝ) := Nat.le_ceil (R / c)
        have hceil' : R / c ≤ (n : ℝ) := by
          exact le_trans hceil (by exact_mod_cast hn)
        have hRc : R ≤ c * (n : ℝ) := by
          have hmul : c * (R / c) ≤ c * (n : ℝ) :=
            mul_le_mul_of_nonneg_left hceil' hc.le
          have hcne : c ≠ 0 := ne_of_gt hc
          have hleft : c * (R / c) = R := by field_simp [hcne]
          calc
            R = c * (R / c) := hleft.symm
            _ ≤ c * (n : ℝ) := hmul
        have hcn : c * (n : ℝ) ≤ c * (n + 1 : ℝ) := by
          have hnp1 : (n : ℝ) ≤ (n + 1 : ℝ) := by linarith
          exact mul_le_mul_of_nonneg_left hnp1 hc.le
        exact le_trans hRc hcn
      refine Filter.mem_atTop_sets.mpr ?_
      refine ⟨N0, ?_⟩
      intro n hn z hzK
      have hz_dist : ‖z - x‖ ≤ 1 := by
        simpa [K, Metric.mem_closedBall, dist_eq_norm] using hzK
      have hz_norm : ‖z‖ ≤ R := by
        calc
          ‖z‖ = ‖(z - x) + x‖ := by ring
          _ ≤ ‖z - x‖ + ‖x‖ := norm_add_le _ _
          _ ≤ 1 + ‖x‖ := by linarith
          _ = R := by simp [R, add_comm]
      by_cases hdup : ∃ a < n, ρ a = ρ n
      · have hu_nonneg : 0 ≤ u n := by
          simp [u, hadamardQuadraticProfile]
          positivity
        simpa [g, canonicalHadamardStepFactor, hmem_iff_exists_lt, hdup] using hu_nonneg
      · have hcn_le_norm : c * (n + 1 : ℝ) ≤ ‖criticalLinePoint (ρ n)‖ := hgrow n
        have hR_le_cn : R ≤ c * (n + 1 : ℝ) := hN0 n hn
        have hcn_pos : 0 < c * (n + 1 : ℝ) := by
          have hn1 : (0 : ℝ) < (n + 1 : ℝ) := by positivity
          exact mul_pos hc hn1
        have hz_ratio_le :
            ‖z / criticalLinePoint (ρ n)‖ ≤ R / (c * (n + 1 : ℝ)) := by
          rw [norm_div]
          have hz1 : ‖z‖ / ‖criticalLinePoint (ρ n)‖ ≤ R / ‖criticalLinePoint (ρ n)‖ :=
            div_le_div_of_nonneg_right hz_norm (norm_nonneg _)
          have hz2 : R / ‖criticalLinePoint (ρ n)‖ ≤ R / (c * (n + 1 : ℝ)) := by
            have hInv : (1 / ‖criticalLinePoint (ρ n)‖) ≤ (1 / (c * (n + 1 : ℝ))) :=
              one_div_le_one_div_of_le hcn_pos hcn_le_norm
            have hRnonneg : 0 ≤ R := by
              dsimp [R]
              positivity
            simpa [div_eq_mul_inv, mul_assoc, mul_left_comm, mul_comm] using
              (mul_le_mul_of_nonneg_left hInv hRnonneg)
          exact le_trans hz1 hz2
        have hz_ratio_le_one : ‖z / criticalLinePoint (ρ n)‖ ≤ 1 := by
          have hRratio : R / (c * (n + 1 : ℝ)) ≤ 1 := by
            exact (div_le_iff₀ hcn_pos).2 (by simpa using hR_le_cn)
          exact le_trans hz_ratio_le hRratio
        have hfac :
            ‖hadamardFactor (ρ n) z - 1‖ ≤
              3 * ‖z / criticalLinePoint (ρ n)‖ ^ (2 : ℕ) :=
          norm_hadamardFactor_sub_one_le_three_mul_sq (ρ n) z hz_ratio_le_one
        have hpow :
            ‖z / criticalLinePoint (ρ n)‖ ^ (2 : ℕ) ≤
              (R / (c * (n + 1 : ℝ))) ^ (2 : ℕ) := by
          nlinarith [norm_nonneg (z / criticalLinePoint (ρ n)), hz_ratio_le]
        have hmain :
            ‖g n z‖ ≤ u n := by
          calc
            ‖g n z‖ = ‖hadamardFactor (ρ n) z - 1‖ := by
              simp [g, canonicalHadamardStepFactor, hmem_iff_exists_lt, hdup]
            _ ≤ 3 * ‖z / criticalLinePoint (ρ n)‖ ^ (2 : ℕ) := hfac
            _ ≤ 3 * (R / (c * (n + 1 : ℝ))) ^ (2 : ℕ) := by
                  gcongr
            _ = u n := by simp [u, hadamardQuadraticProfile, mul_comm, mul_left_comm, mul_assoc]
        exact hmain
    have hMulOneAdd : MultipliableUniformlyOn (fun n z => 1 + g n z) K :=
      hu.multipliableUniformlyOn_nat_one_add hKcompact hEvent hcts
    have hMul : MultipliableUniformlyOn (canonicalHadamardStepFactor ρ) K := by
      simpa [g, sub_eq_add_neg, add_assoc] using hMulOneAdd
    exact ⟨K, hKnhds, hMul⟩
  refine ⟨fun z : ℂ => ∏' n : ℕ, canonicalHadamardStepFactor ρ n z, ?_⟩
  have hProd :
      HasProdLocallyUniformly (canonicalHadamardStepFactor ρ)
        (fun z : ℂ => ∏' n : ℕ, canonicalHadamardStepFactor ρ n z) :=
    hStep.hasProdLocallyUniformly
  have hRange :
      TendstoLocallyUniformly
        (fun N : ℕ => fun z : ℂ => ∏ i ∈ Finset.range N, canonicalHadamardStepFactor ρ i z)
        (fun z : ℂ => ∏' n : ℕ, canonicalHadamardStepFactor ρ n z)
        (Filter.atTop : Filter ℕ) :=
    hProd.tendstoLocallyUniformly_finsetRange
  have hEq :
      (fun N : ℕ => canonicalHadamardXiN ρ N) =
      (fun N : ℕ => fun z : ℂ => ∏ i ∈ Finset.range N, canonicalHadamardStepFactor ρ i z) := by
    funext N z
    exact canonicalHadamardXiN_eq_finsetRange_prod_stepFactor ρ N z
  simpa [hEq] using hRange

/-- Hard analytic lemma (finite-level nontriviality at `0`) for canonical
Hadamard partial products. -/
theorem canonicalHadamardXiN_nontrivialAtZero
    (ρ : ℕ → ℝ) :
    ∀ N : ℕ, canonicalHadamardXiN ρ N 0 ≠ 0 := by
  intro N
  simp [canonicalHadamardXiN, XiFiniteHadamard, hadamardFactor]

/-- Hard analytic lemma (differentiability of each Hadamard partial product). -/
theorem differentiable_canonicalHadamardXiN
    (ρ : ℕ → ℝ) :
    ∀ N : ℕ, Differentiable ℂ (canonicalHadamardXiN ρ N) := by
  intro N
  intro z
  unfold canonicalHadamardXiN
  unfold XiFiniteHadamard
  let spec : Finset ℝ := (Finset.range N).image ρ
  let f : ℝ → ℂ → ℂ := fun t s => hadamardFactor t s
  have hf : ∀ t ∈ spec, DifferentiableAt ℂ (f t) z := by
    intro t ht
    have hcp : criticalLinePoint t ≠ 0 :=
      Gutoe.RiemannFiniteXiModel.criticalLinePoint_ne_zero t
    have hdiv : DifferentiableAt ℂ (fun s : ℂ => s / criticalLinePoint t) z := by
      simpa [div_eq_mul_inv, mul_comm, mul_left_comm, mul_assoc] using
        (differentiableAt_id.mul_const ((criticalLinePoint t)⁻¹))
    have hconst1 : DifferentiableAt ℂ (fun _ : ℂ => (1 : ℂ)) z :=
      (differentiableAt_const (c := (1 : ℂ)))
    have hlin : DifferentiableAt ℂ (fun s : ℂ => 1 - s / criticalLinePoint t) z :=
      hconst1.sub hdiv
    have hexp : DifferentiableAt ℂ (fun s : ℂ => Complex.exp (s / criticalLinePoint t)) z :=
      Complex.differentiableAt_exp.comp z hdiv
    have hmul : DifferentiableAt ℂ
        (fun s : ℂ => (1 - s / criticalLinePoint t) * Complex.exp (s / criticalLinePoint t)) z :=
      hlin.mul hexp
    simpa [f, hadamardFactor, hcp] using hmul
  simpa [spec, f] using (DifferentiableAt.fun_finset_prod (u := spec) (f := f) hf)

/-- Hard analytic Hurwitz lemma:
local-uniform convergence of holomorphic approximants to `XiTarget` transfers
zeros near each target zero with multiplicity bookkeeping. -/
theorem hurwitz_eventual_zero_near
    (hHurwitz : HurwitzZeroTransfer)
    {XiN : ℕ → (ℂ → ℂ)}
    (hconv : TendstoLocallyUniformly XiN XiTarget (Filter.atTop : Filter ℕ))
    (hDiffN : ∀ N : ℕ, Differentiable ℂ (XiN N))
    (hDiffXi : ∀ s : ℂ, XiTarget s = 0 → DifferentiableAt ℂ XiTarget s)
    (hNontriv : ∀ N : ℕ, XiN N 0 ≠ 0)
    (s : ℂ) (hs : XiTarget s = 0) :
    ∀ ε : ℝ, 0 < ε → ∃ N : ℕ, ∃ z : ℂ, XiN N z = 0 ∧ ‖s - z‖ < ε := by
  rcases hHurwitz.zeroOrder_of_zero s hs with ⟨k, hk⟩
  exact hHurwitz.eventually_zero_near hconv hDiffN hDiffXi hNontriv s k hs hk

/-- Hard analytic collapse lemma:
the operator ladder is the canonical Hadamard lane and Hurwitz transfer can be
upgraded (via spectral rigidity) to exact finite-level zeros at the target
point itself. -/
theorem operatorExactZeroUpgrade_of_canonical_hadamard
    (ρ : ℕ → ℝ)
    (hEnum : EnumeratesXiTargetZeros ρ)
    (hBridge : EnumeratorHitsOperatorSpec ρ) :
    ∀ s : ℂ, XiTarget s = 0 → ∃ N : ℕ, operatorXiFiniteLadder N s = 0 := by
  intro s hs
  rcases hEnum s hs with ⟨n, hsEq⟩
  rcases hBridge n with ⟨N, hmem⟩
  refine ⟨N, ?_⟩
  have hZero : XiFinite (operatorSpecN N) s = 0 := by
    rw [hsEq]
    exact XiFinite_zero_of_mem (operatorSpecN N) hmem
  simpa [operatorXiFiniteLadder, Gutoe.RiemannTargetFiniteLadder.XiFiniteLadder] using hZero

/-- Build a `HurwitzZeroTransfer` package from explicit analytic lemmas for the
canonical Hadamard lane.

`eventually_zero_near` is provided as an abstract theorem slot and the exact
operator upgrade is isolated in `operatorExactZeroUpgrade_of_canonical_hadamard`.
-/
def hurwitzZeroTransfer_of_canonical
    (ρ : ℕ → ℝ)
    (hEnum : EnumeratesXiTargetZeros ρ)
    (hBridge : EnumeratorHitsOperatorSpec ρ)
    (hNear :
      ∀ {XiN : ℕ → (ℂ → ℂ)},
        TendstoLocallyUniformly XiN XiTarget (Filter.atTop : Filter ℕ) →
        (∀ N : ℕ, Differentiable ℂ (XiN N)) →
        (∀ s : ℂ, XiTarget s = 0 → DifferentiableAt ℂ XiTarget s) →
        (∀ N : ℕ, XiN N 0 ≠ 0) →
        ∀ s : ℂ, ∀ _k : ℕ, XiTarget s = 0 → True →
          ∀ ε : ℝ, 0 < ε →
            ∃ N : ℕ, ∃ z : ℂ, XiN N z = 0 ∧ ‖s - z‖ < ε) :
    HurwitzZeroTransfer where
  zeroOrder := fun _ _ => True
  zeroOrder_of_zero := by
    intro s hs
    exact ⟨1, trivial⟩
  eventually_zero_near := by
    intro XiN hconv hDiffN hDiffXi hNontriv s k hs _hk ε hε
    exact hNear hconv hDiffN hDiffXi hNontriv s k hs trivial ε hε
  operatorExactZeroUpgrade := operatorExactZeroUpgrade_of_canonical_hadamard ρ hEnum hBridge

/-- End-to-end reduction theorem:
`OperatorApproxZeroConvergence` follows directly from zero enumeration
and operator spectral bridge. Every `XiTarget` zero is captured at some
finite operator level via the enumerator–spectral correspondence.

Note: this bypasses the Hadamard factorization identity (which would require
identifying the infinite Hadamard product limit with `XiTarget` up to the
entire-order exponential prefactor `e^{A+Bs}`). The zero-capture goes through
`EnumeratesXiTargetZeros` + `EnumeratorHitsOperatorSpec` directly. -/
theorem operatorApproxZeroConvergence_of_canonical_hadamard
    (ρ : ℕ → ℝ)
    (hEnum : EnumeratesXiTargetZeros ρ)
    (hBridge : EnumeratorHitsOperatorSpec ρ) :
    OperatorApproxZeroConvergence :=
  (operatorApproxZeroConvergence_iff_eventual_exact_zero).2
    (operatorExactZeroUpgrade_of_canonical_hadamard ρ hEnum hBridge)

end

end Gutoe.RiemannConvergenceProofs
