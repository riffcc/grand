import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannSelfAdjoint
import Gutoe.RiemannTargetFiniteLadder
import Gutoe.RiemannLimitBridge
import Gutoe.RiemannConvergenceTransfer
import Gutoe.RiemannFinalTarget
import Gutoe.RiemannFiniteXiModel

namespace Gutoe.RiemannOperatorLadder

open Gutoe.RiemannCore
open Gutoe.RiemannSelfAdjoint
open Gutoe.RiemannTargetFiniteLadder
open Gutoe.RiemannLimitBridge
open Gutoe.RiemannConvergenceTransfer
open Gutoe.RiemannFinalTarget
open Gutoe.RiemannFiniteXiModel

noncomputable section
open scoped Topology

/-- Complex-lifted structural matrix lane for spectral statements. -/
def structuralRiemannMatrixC (n : ℕ) : Matrix (Fin n) (Fin n) ℂ :=
  fun i j => (structuralRiemannMatrix n i j : ℂ)

/-- Endomorphism induced by the structural matrix on coordinate vectors. -/
def structuralRiemannEnd (n : ℕ) : Module.End ℂ (Fin n → ℂ) :=
  Matrix.mulVecLin (structuralRiemannMatrixC n)

/-- A real ordinate is spectral for level `n` if it is an eigenvalue of the
complex structural endomorphism at that level. -/
def ordinateIsEigenvalue (n : ℕ) (t : ℝ) : Prop :=
  (t : ℂ) ∈ spectrum ℂ (structuralRiemannMatrixC n)

/-- Direct operator-defined spectrum lane at level `N`:
real ordinates whose complex lift lies in the matrix spectrum. -/
def operatorSpecSet (N : ℕ) : Set ℝ :=
  fun t => (t : ℂ) ∈ spectrum ℂ (structuralRiemannMatrixC (N + 1))

theorem ordinateIsEigenvalue_iff_mem_operatorSpecSet
    (N : ℕ) (t : ℝ) :
    ordinateIsEigenvalue (N + 1) t ↔ t ∈ operatorSpecSet N := by
  rfl

/-- The structural complex matrix is Hermitian (real symmetric data lifted to `ℂ`). -/
theorem structuralRiemannMatrixC_isHermitian (n : ℕ) :
    (structuralRiemannMatrixC n).IsHermitian := by
  ext i j
  simp [Matrix.conjTranspose, structuralRiemannMatrixC]
  exact (structuralRiemannMatrix_isSymm n).apply i j

/-- Canonical real eigenvalue enumerator for the structural matrix at level `N+1`. -/
noncomputable def operatorEigenvalues (N : ℕ) : Fin (N + 1) → ℝ :=
  (structuralRiemannMatrixC_isHermitian (N + 1)).eigenvalues

/-- Concrete operator-derived finite ladder:
real eigenvalues of the structural Hermitian matrix at each level `N+1`. -/
noncomputable def operatorSpecN (N : ℕ) : Finset ℝ :=
  Finset.univ.image (operatorEigenvalues N)

/-- Every ordinate listed in `operatorSpecN N` lies in the real spectrum of
the structural matrix at level `N+1`. -/
theorem mem_operatorSpecN_implies_mem_real_spectrum
    {N : ℕ} {t : ℝ} (ht : t ∈ operatorSpecN N) :
    t ∈ spectrum ℝ (structuralRiemannMatrixC (N + 1)) := by
  rcases Finset.mem_image.mp ht with ⟨i, _hi, rfl⟩
  exact Matrix.IsHermitian.eigenvalues_mem_spectrum_real
    (hA := structuralRiemannMatrixC_isHermitian (N + 1)) i

/-- Exact membership characterization for the concrete operator ladder:
listed ordinates are exactly spectral ordinates of the structural matrix. -/
theorem mem_operatorSpecN_iff_ordinateIsEigenvalue
    (N : ℕ) (t : ℝ) :
    t ∈ operatorSpecN N ↔ ordinateIsEigenvalue (N + 1) t := by
  constructor
  · intro ht
    have hreal : t ∈ spectrum ℝ (structuralRiemannMatrixC (N + 1)) :=
      mem_operatorSpecN_implies_mem_real_spectrum ht
    exact (spectrum.algebraMap_mem_iff ℂ).2 hreal
  · intro ht
    have hreal : t ∈ spectrum ℝ (structuralRiemannMatrixC (N + 1)) :=
      (spectrum.algebraMap_mem_iff ℂ).1 ht
    have hrange :
        t ∈ Set.range ((structuralRiemannMatrixC_isHermitian (N + 1)).eigenvalues) := by
      simpa [Matrix.IsHermitian.spectrum_real_eq_range_eigenvalues
        (hA := structuralRiemannMatrixC_isHermitian (N + 1))] using hreal
    rcases hrange with ⟨i, rfl⟩
    exact Finset.mem_image.mpr ⟨i, Finset.mem_univ i, rfl⟩

/-- Concrete operator ladder membership is exactly operator-spectrum set membership. -/
theorem mem_operatorSpecN_iff_mem_operatorSpecSet
    (N : ℕ) (t : ℝ) :
    t ∈ operatorSpecN N ↔ t ∈ operatorSpecSet N := by
  exact (mem_operatorSpecN_iff_ordinateIsEigenvalue N t).trans
    (ordinateIsEigenvalue_iff_mem_operatorSpecSet N t)

/-- Operator-native finite ladder witness:
each finite level list is exactly the real ordinates detected as operator
eigenvalues for that level. -/
structure OperatorSpecLadder where
  specN : ℕ → Finset ℝ
  eigenvalue_exact :
    ∀ N : ℕ, ∀ t : ℝ, t ∈ specN N ↔ ordinateIsEigenvalue (N + 1) t

/-- Canonical operator ladder witness coming directly from `operatorSpecN`. -/
def canonicalOperatorSpecLadder : OperatorSpecLadder where
  specN := operatorSpecN
  eigenvalue_exact := mem_operatorSpecN_iff_ordinateIsEigenvalue

/-- Obligation-1 in operator-native form. -/
def OperatorNontrivialCapture (hOp : OperatorSpecLadder) : Prop :=
  RiemannNontrivialLadderZeroCapture hOp.specN

/-- Direct operator-set capture surface (no external reference data). -/
def OperatorSetNontrivialCapture : Prop :=
  ∀ s : ℂ, riemannZeta s = 0 →
    (¬ ∃ n : ℕ, s = -2 * (n + 1)) →
    s ≠ 1 →
    ∃ N : ℕ, ∃ t : ℝ, t ∈ operatorSpecSet N ∧ s = criticalLinePoint t

/-- Concrete finite-ladder capture obligation specialized to the
operator-constructed `operatorSpecN`. -/
def OperatorEnumeratedNontrivialCapture : Prop :=
  RiemannNontrivialLadderZeroCapture operatorSpecN

/-- Operator-native equivalence of endgame capture surfaces. -/
theorem operator_nontrivial_capture_iff_xiTarget_capture
    (hOp : OperatorSpecLadder) :
    OperatorNontrivialCapture hOp ↔ XiTargetLadderZeroCapture hOp.specN := by
  exact nontrivial_capture_iff_xiTarget_capture hOp.specN

/-- Operator-native RH closure:
once nontrivial-`ζ` capture is established for the operator ladder, RH follows. -/
theorem mathlibRH_of_operator_nontrivial_capture
    (hOp : OperatorSpecLadder)
    (hCap : OperatorNontrivialCapture hOp) :
    RiemannHypothesis := by
  exact mathlibRH_of_nontrivial_capture hOp.specN hCap

/-- RH closure directly from nontrivial capture over the concrete operator
eigenvalue ladder `operatorSpecN`. -/
theorem mathlibRH_of_operator_enumerated_nontrivial_capture
    (hCap : OperatorEnumeratedNontrivialCapture) :
    RiemannHypothesis := by
  exact mathlibRH_of_nontrivial_capture operatorSpecN hCap

/-- Concrete operator finite-product ladder used in the convergence lane. -/
def operatorXiFiniteLadder : ℕ → (ℂ → ℂ) :=
  XiFiniteLadder operatorSpecN

/-- Critical-line points are never zero. -/
theorem criticalLinePoint_ne_zero (t : ℝ) :
    criticalLinePoint t ≠ 0 := by
  intro h
  have hre : (criticalLinePoint t).re = 0 := by
    simpa [h] using congrArg Complex.re h
  have hhalf : (1 / 2 : ℝ) = 0 := by
    simpa [criticalLinePoint_re] using hre
  norm_num at hhalf

/-- The finite model `XiFinite spec` is complex-differentiable everywhere. -/
theorem differentiable_XiFinite
    (spec : Finset ℝ) :
    Differentiable ℂ (XiFinite spec) := by
  intro z
  unfold XiFinite
  let f : ℝ → ℂ → ℂ := fun t s => s - criticalLinePoint t
  have hf : ∀ t ∈ spec, DifferentiableAt ℂ (f t) z := by
    intro t ht
    simpa [f] using (differentiableAt_id.sub_const (criticalLinePoint t))
  simpa [f] using (DifferentiableAt.fun_finset_prod (u := spec) (f := f) hf)

/-- Each finite level of the concrete operator Xi ladder is differentiable. -/
theorem differentiable_operatorXiFiniteLadder (N : ℕ) :
    Differentiable ℂ (operatorXiFiniteLadder N) := by
  simpa [operatorXiFiniteLadder, XiFiniteLadder] using
    (differentiable_XiFinite (operatorSpecN N))

/-- The finite operator Xi products are nontrivial: they never vanish at `0`. -/
theorem xiFinite_operatorSpecN_zero_ne (N : ℕ) :
    XiFinite (operatorSpecN N) 0 ≠ 0 := by
  unfold XiFinite
  refine Finset.prod_ne_zero_iff.mpr ?_
  intro t ht
  simp [criticalLinePoint_ne_zero t]

/-- Constant tolerance profile for the concrete operator ladder. -/
def operatorTolConst (τ : ℝ) : ℕ → ℝ := fun _ => τ

/-- Local-uniform convergence implies nonzero-tolerance approximate zero capture
at every target zero. -/
theorem operatorApprox_of_locallyUniform_constTol
    (hconv : TendstoLocallyUniformly operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ))
    {τ : ℝ} (hτ : 0 < τ) :
    ApproxZeroConvergence XiTarget operatorXiFiniteLadder (operatorTolConst τ) := by
  intro s hsXi
  have hconvOn : TendstoLocallyUniformlyOn operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ) (Set.univ) := by
    simpa [tendstoLocallyUniformlyOn_univ] using hconv
  have hpt : Filter.Tendsto (fun N : ℕ => operatorXiFiniteLadder N s)
      (Filter.atTop : Filter ℕ) (𝓝 (XiTarget s)) :=
    hconvOn.tendsto_at (by simp)
  have hball : ∀ᶠ N : ℕ in (Filter.atTop : Filter ℕ),
      ‖operatorXiFiniteLadder N s - XiTarget s‖ < τ := by
    simpa [Metric.ball, dist_eq_norm] using
      (hpt (Metric.ball_mem_nhds (XiTarget s) hτ))
  rcases hball.exists_forall_of_atTop with ⟨N0, hN0⟩
  refine ⟨N0, ?_⟩
  have hlt : ‖operatorXiFiniteLadder N0 s - XiTarget s‖ < τ := hN0 N0 (le_rfl)
  have hlt0 : ‖operatorXiFiniteLadder N0 s‖ < τ := by
    simpa [hsXi] using hlt
  exact le_of_lt hlt0

/-- If local-uniform convergence is available and the constant tolerance level
is rigid, one gets zero-forward transfer for the operator Xi ladder. -/
theorem operatorZeroForward_of_locallyUniform_and_constRigidity
    (hconv : TendstoLocallyUniformly operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ))
    {τ : ℝ} (hτ : 0 < τ)
    (hRigid : SpectralRigidity operatorXiFiniteLadder (operatorTolConst τ)) :
    ZeroForwardTransfer XiTarget operatorXiFiniteLadder := by
  exact zeroForward_of_convergence_and_rigidity
    XiTarget
    operatorXiFiniteLadder
    (operatorTolConst τ)
    (operatorApprox_of_locallyUniform_constTol hconv hτ)
    hRigid

/-- Operator-ladder convergence obligation at zero tolerance. -/
def OperatorApproxZeroConvergence : Prop :=
  ApproxZeroConvergence
    XiTarget
    operatorXiFiniteLadder
    Gutoe.RiemannTargetFiniteLadder.tolZero

/-- For the concrete operator ladder, zero-tolerance approximation is exactly
finite-ladder `XiTarget` zero-capture. -/
theorem operatorApproxZero_iff_xiTarget_capture :
    OperatorApproxZeroConvergence ↔ XiTargetLadderZeroCapture operatorSpecN := by
  simpa [OperatorApproxZeroConvergence, operatorXiFiniteLadder] using
    (approxZero_tolZero_iff_zeroCapture operatorSpecN)

/-- Convergence closure in the operator lane:
zero-tolerance approximation implies nontrivial-`ζ` capture for `operatorSpecN`. -/
theorem operator_enumerated_capture_of_approxZero
    (hApprox : OperatorApproxZeroConvergence) :
    OperatorEnumeratedNontrivialCapture := by
  exact (nontrivial_capture_iff_xiTarget_capture operatorSpecN).2
    ((operatorApproxZero_iff_xiTarget_capture).1 hApprox)

/-- RH closure from operator-ladder convergence obligation. -/
theorem mathlibRH_of_operator_approxZero
    (hApprox : OperatorApproxZeroConvergence) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_enumerated_nontrivial_capture
    (operator_enumerated_capture_of_approxZero hApprox)

/-- Positive-threshold rigidity surface: there exists at least one strictly
positive tolerance level at which small values force exact finite-level zeros. -/
def OperatorPositiveRigidity : Prop :=
  ∃ τ : ℝ, 0 < τ ∧ SpectralRigidity operatorXiFiniteLadder (operatorTolConst τ)

/-- Local-uniform convergence plus any positive-threshold rigidity already
forces the zero-tolerance operator approximation obligation. -/
theorem operatorApproxZero_of_locallyUniform_and_positiveRigidity
    (hconv : TendstoLocallyUniformly operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ))
    (hPosRigid : OperatorPositiveRigidity) :
    OperatorApproxZeroConvergence := by
  rcases hPosRigid with ⟨τ, hτ, hRigid⟩
  intro s hsXi
  have hApproxτ :
      ApproxZeroConvergence XiTarget operatorXiFiniteLadder (operatorTolConst τ) :=
    operatorApprox_of_locallyUniform_constTol hconv hτ
  rcases hApproxτ s hsXi with ⟨N, hNτ⟩
  have hZero : operatorXiFiniteLadder N s = 0 := hRigid N s hNτ
  refine ⟨N, ?_⟩
  simpa [Gutoe.RiemannTargetFiniteLadder.tolZero, hZero]

/-- RH closure from the analytic surface
`local-uniform convergence + positive-threshold rigidity`. -/
theorem mathlibRH_of_locallyUniform_and_positiveRigidity
    (hconv : TendstoLocallyUniformly operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ))
    (hPosRigid : OperatorPositiveRigidity) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_approxZero
    (operatorApproxZero_of_locallyUniform_and_positiveRigidity hconv hPosRigid)

/-- The two operator capture surfaces are equivalent:
finite ladder capture (`operatorSpecN`) and direct spectrum-set capture (`operatorSpecSet`). -/
theorem operator_enumerated_capture_iff_operator_set_capture :
    OperatorEnumeratedNontrivialCapture ↔ OperatorSetNontrivialCapture := by
  constructor
  · intro hCap
    intro s hs htriv h1
    rcases hCap s hs htriv h1 with ⟨N, t, ht, hsEq⟩
    exact ⟨N, t, (mem_operatorSpecN_iff_mem_operatorSpecSet N t).1 ht, hsEq⟩
  · intro hCap
    intro s hs htriv h1
    rcases hCap s hs htriv h1 with ⟨N, t, ht, hsEq⟩
    exact ⟨N, t, (mem_operatorSpecN_iff_mem_operatorSpecSet N t).2 ht, hsEq⟩

/-- Endgame equivalence in operator-native form:
the convergence obligation is equivalent to direct operator-spectrum capture. -/
theorem operatorApproxZero_iff_operator_set_capture :
    OperatorApproxZeroConvergence ↔ OperatorSetNontrivialCapture := by
  exact (operatorApproxZero_iff_xiTarget_capture.trans
    (nontrivial_capture_iff_xiTarget_capture operatorSpecN).symm).trans
      operator_enumerated_capture_iff_operator_set_capture

/-- A finite operator-ladder witness immediately yields set-capture on the
explicit operator spectrum lane. -/
theorem operator_set_capture_of_operator_ladder
    (hOp : OperatorSpecLadder)
    (hCap : OperatorNontrivialCapture hOp) :
    OperatorSetNontrivialCapture := by
  intro s hs htriv h1
  rcases hCap s hs htriv h1 with ⟨N, t, ht, hsEq⟩
  refine ⟨N, t, ?_, hsEq⟩
  exact (ordinateIsEigenvalue_iff_mem_operatorSpecSet N t).1
    ((hOp.eigenvalue_exact N t).1 ht)

/-- RH closure from the direct operator spectrum set-capture surface. -/
theorem mathlibRH_of_operator_set_nontrivial_capture
    (hCap : OperatorSetNontrivialCapture) :
    RiemannHypothesis := by
  intro s hs htriv h1
  rcases hCap s hs htriv h1 with ⟨_N, t, _ht, hsEq⟩
  simpa [hsEq, onCriticalLine, criticalLinePoint_re] using (criticalLinePoint_re t)

/-- Approximate operator-spectrum capture:
every `XiTarget` zero can be approximated arbitrarily well by critical-line
points coming from finite operator spectra. -/
def OperatorApproximateCapture : Prop :=
  ∀ s : ℂ, XiTarget s = 0 → ∀ ε : ℝ, 0 < ε →
    ∃ N : ℕ, ∃ t : ℝ, t ∈ operatorSpecSet N ∧ ‖s - criticalLinePoint t‖ < ε

/-- Hurwitz-output surface specialized to the concrete operator lane:
every `XiTarget` zero can be approximated by an actual finite-level zero of the
operator finite product. -/
def OperatorHurwitzZeroApproxTransfer : Prop :=
  ∀ s : ℂ, XiTarget s = 0 → ∀ ε : ℝ, 0 < ε →
    ∃ N : ℕ, ∃ z : ℂ,
      XiFinite (operatorSpecN N) z = 0 ∧ ‖s - z‖ < ε

/-- The operator Hurwitz-output surface implies the internal approximate-capture
surface used by the RH closure theorem. -/
theorem operatorApproximateCapture_of_hurwitzTransfer
    (hHurwitz : OperatorHurwitzZeroApproxTransfer) :
    OperatorApproximateCapture := by
  intro s hsXi ε hε
  rcases hHurwitz s hsXi ε hε with ⟨N, z, hz0, hdist⟩
  rcases (XiFinite_zero_iff_exists (operatorSpecN N) z).1 hz0 with ⟨t, htN, hzEq⟩
  refine ⟨N, t, (mem_operatorSpecN_iff_mem_operatorSpecSet N t).1 htN, ?_⟩
  simpa [hzEq] using hdist

/-- The approximate-capture surface also yields the Hurwitz-output surface:
finite-level critical-line witnesses are concrete finite-product zeros. -/
theorem operatorHurwitzTransfer_of_operatorApproximateCapture
    (hApprox : OperatorApproximateCapture) :
    OperatorHurwitzZeroApproxTransfer := by
  intro s hsXi ε hε
  rcases hApprox s hsXi ε hε with ⟨N, t, htSet, hdist⟩
  have htN : t ∈ operatorSpecN N :=
    (mem_operatorSpecN_iff_mem_operatorSpecSet N t).2 htSet
  refine ⟨N, criticalLinePoint t, XiFinite_zero_of_mem (operatorSpecN N) htN, ?_⟩
  simpa using hdist

/-- Exact equivalence between the operator Hurwitz-output and approximate-capture
surfaces. This isolates the remaining analytic gap to either form. -/
theorem operatorHurwitz_iff_operatorApproximateCapture :
    OperatorHurwitzZeroApproxTransfer ↔ OperatorApproximateCapture := by
  constructor
  · exact operatorApproximateCapture_of_hurwitzTransfer
  · exact operatorHurwitzTransfer_of_operatorApproximateCapture

/-- RH closure from approximate operator-spectrum capture.
This is the exact interface needed for a future Hurwitz-style transfer lemma. -/
theorem mathlibRH_of_operator_approximate_capture
    (hApprox : OperatorApproximateCapture) :
    RiemannHypothesis := by
  intro s hs htriv h1
  have hsXi : XiTarget s = 0 := nontrivialZeroTransferToXiTarget s hs htriv h1
  have hsmall : ∀ ε : ℝ, 0 < ε → |s.re - (1 / 2 : ℝ)| < ε := by
    intro ε hε
    rcases hApprox s hsXi ε hε with ⟨N, t, _ht, hdist⟩
    have hreLe : |(s - criticalLinePoint t).re| ≤ ‖s - criticalLinePoint t‖ := by
      simpa using (Complex.abs_re_le_norm (s - criticalLinePoint t))
    have hreLt : |(s - criticalLinePoint t).re| < ε := lt_of_le_of_lt hreLe hdist
    simpa [criticalLinePoint_re, Complex.sub_re] using hreLt
  have hEq : s.re = (1 / 2 : ℝ) := by
    by_contra hne
    let δ : ℝ := |s.re - (1 / 2 : ℝ)| / 2
    have hδpos : 0 < δ := by
      unfold δ
      exact half_pos (abs_pos.mpr (sub_ne_zero.mpr hne))
    have hlt : |s.re - (1 / 2 : ℝ)| < δ := hsmall δ hδpos
    have hge : δ ≤ |s.re - (1 / 2 : ℝ)| := by
      unfold δ
      nlinarith [abs_nonneg (s.re - (1 / 2 : ℝ))]
    exact (not_lt_of_ge hge) hlt
  exact hEq

/-- RH closure from the operator Hurwitz-output obligation. -/
theorem mathlibRH_of_operator_hurwitz_zero_approx
    (hHurwitz : OperatorHurwitzZeroApproxTransfer) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_approximate_capture
    (operatorApproximateCapture_of_hurwitzTransfer hHurwitz)

/-- Zero-tolerance operator approximation already implies the Hurwitz-output
surface (with exact finite-level zeros at the target point itself). -/
theorem operatorHurwitz_of_operatorApproxZero
    (hApprox0 : OperatorApproxZeroConvergence) :
    OperatorHurwitzZeroApproxTransfer := by
  have hCapXi : XiTargetLadderZeroCapture operatorSpecN :=
    (operatorApproxZero_iff_xiTarget_capture).1 hApprox0
  intro s hsXi ε hε
  rcases hCapXi s hsXi with ⟨N, t, htN, hsEq⟩
  refine ⟨N, criticalLinePoint t, XiFinite_zero_of_mem (operatorSpecN N) htN, ?_⟩
  subst hsEq
  simpa using hε

/-- RH closure from the zero-tolerance operator approximation routed through
the Hurwitz-output surface. -/
theorem mathlibRH_of_operator_approxZero_via_hurwitz
    (hApprox0 : OperatorApproxZeroConvergence) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_hurwitz_zero_approx
    (operatorHurwitz_of_operatorApproxZero hApprox0)

end

end Gutoe.RiemannOperatorLadder
