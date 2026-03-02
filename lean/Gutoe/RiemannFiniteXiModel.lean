import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannBridge
import Gutoe.RiemannLimitBridge
import Gutoe.RiemannConvergenceTransfer

namespace Gutoe.RiemannFiniteXiModel

open Gutoe.RiemannCore
open Gutoe.RiemannBridge
open Gutoe.RiemannLimitBridge
open Gutoe.RiemannConvergenceTransfer

noncomputable section
open scoped BigOperators

/-- Finite spectral set as a `Set`. -/
def finiteSpecSet (spec : Finset ℝ) : Set ℝ := fun t => t ∈ spec

/-- Finite spectral Xi model:
    one zero factor per spectral ordinate. -/
def XiFinite (spec : Finset ℝ) : ℂ → ℂ :=
  fun s => Finset.prod spec (fun t => (s - criticalLinePoint t : ℂ))

/-- Critical-line points are never zero, so the Hadamard-normalized factors are
well-defined for every real ordinate. -/
theorem criticalLinePoint_ne_zero (t : ℝ) :
    criticalLinePoint t ≠ 0 := by
  intro h
  have hre : (criticalLinePoint t).re = 0 := by
    simpa [h] using congrArg Complex.re h
  have hhalf : (1 / 2 : ℝ) = 0 := by
    simpa [criticalLinePoint_re] using hre
  norm_num at hhalf

/-- Canonical genus-1 Hadamard factor attached to ordinate `t`. -/
def hadamardFactor (t : ℝ) (s : ℂ) : ℂ :=
  (1 - s / criticalLinePoint t) * Complex.exp (s / criticalLinePoint t)

/-- Finite Hadamard-normalized Xi model. This has the same zero set as `XiFinite`
but improved asymptotic behavior for infinite-level limits. -/
def XiFiniteHadamard (spec : Finset ℝ) : ℂ → ℂ :=
  fun s => Finset.prod spec (fun t => hadamardFactor t s)

/-- A canonical Hadamard factor vanishes exactly at its designated critical-line
point. -/
theorem hadamardFactor_eq_zero_iff
    (t : ℝ) (s : ℂ) :
    hadamardFactor t s = 0 ↔ s = criticalLinePoint t := by
  constructor
  · intro h
    have hmul := mul_eq_zero.mp h
    rcases hmul with hlin | hexp
    · have hneq : criticalLinePoint t ≠ 0 := criticalLinePoint_ne_zero t
      have hratio : (1 : ℂ) = s / criticalLinePoint t := by
        exact sub_eq_zero.mp hlin
      have hratio' : (criticalLinePoint t : ℂ) = s := by
        field_simp [hneq] at hratio
        simpa [mul_comm] using hratio
      simpa using hratio'.symm
    · exact (Complex.exp_ne_zero _ hexp).elim
  · intro hs
    subst hs
    unfold hadamardFactor
    have hneq : criticalLinePoint t ≠ 0 := criticalLinePoint_ne_zero t
    simp [hneq]

/-- One-step factorization in the finite Hadamard model when inserting a new
ordinate into the finite set. -/
theorem XiFiniteHadamard_insert
    (spec : Finset ℝ) {t : ℝ} (ht : t ∉ spec) (s : ℂ) :
    XiFiniteHadamard (insert t spec) s =
      hadamardFactor t s * XiFiniteHadamard spec s := by
  simp [XiFiniteHadamard, Finset.prod_insert, ht, hadamardFactor]

/-- Near-zero second-order control for a canonical Hadamard factor:
for `z = s / criticalLinePoint t` with `‖z‖ ≤ 1`,
`(1 - z)exp(z) - 1` is `O(‖z‖²)` with an explicit constant. -/
theorem norm_hadamardFactor_sub_one_le_three_mul_sq
    (t : ℝ) (s : ℂ)
    (hz : ‖s / criticalLinePoint t‖ ≤ 1) :
    ‖hadamardFactor t s - 1‖ ≤
      3 * ‖s / criticalLinePoint t‖ ^ 2 := by
  let z : ℂ := s / criticalLinePoint t
  have hz_nonneg : 0 ≤ ‖z‖ := norm_nonneg _
  have hz1 : ‖Complex.exp z - 1 - z‖ ≤ ‖z‖ ^ 2 :=
    Complex.norm_exp_sub_one_sub_id_le (by simpa [z] using hz)
  have hz2 : ‖Complex.exp z - 1‖ ≤ ‖z‖ + ‖z‖ ^ 2 := by
    calc
      ‖Complex.exp z - 1‖ = ‖(Complex.exp z - 1 - z) + z‖ := by ring
      _ ≤ ‖Complex.exp z - 1 - z‖ + ‖z‖ := norm_add_le _ _
      _ ≤ ‖z‖ ^ 2 + ‖z‖ := by gcongr
      _ = ‖z‖ + ‖z‖ ^ 2 := by ring
  have hcube_le_sq : ‖z‖ ^ 3 ≤ ‖z‖ ^ 2 := by
    have hz_le : ‖z‖ ≤ 1 := by simpa [z] using hz
    nlinarith [hz_nonneg, hz_le]
  calc
    ‖hadamardFactor t s - 1‖
        = ‖(Complex.exp z - 1 - z) - z * (Complex.exp z - 1)‖ := by
            simp [hadamardFactor, z]
            ring_nf
    _ ≤ ‖Complex.exp z - 1 - z‖ + ‖z * (Complex.exp z - 1)‖ := by
          simpa [sub_eq_add_neg] using norm_add_le (Complex.exp z - 1 - z) (-(z * (Complex.exp z - 1)))
    _ = ‖Complex.exp z - 1 - z‖ + ‖z‖ * ‖Complex.exp z - 1‖ := by simp [norm_mul]
    _ ≤ ‖z‖ ^ 2 + ‖z‖ * (‖z‖ + ‖z‖ ^ 2) := by gcongr
    _ = 2 * ‖z‖ ^ 2 + ‖z‖ ^ 3 := by ring
    _ ≤ 2 * ‖z‖ ^ 2 + ‖z‖ ^ 2 := by gcongr
    _ = 3 * ‖z‖ ^ 2 := by ring
    _ = 3 * ‖s / criticalLinePoint t‖ ^ 2 := by simp [z]

/-- One-step increment bound for the finite Hadamard products. -/
theorem norm_XiFiniteHadamard_insert_sub_le
    (spec : Finset ℝ) {t : ℝ} (ht : t ∉ spec) (s : ℂ) :
    ‖XiFiniteHadamard (insert t spec) s - XiFiniteHadamard spec s‖ ≤
      ‖hadamardFactor t s - 1‖ * ‖XiFiniteHadamard spec s‖ := by
  have hfactor := XiFiniteHadamard_insert spec ht s
  calc
    ‖XiFiniteHadamard (insert t spec) s - XiFiniteHadamard spec s‖
        = ‖(hadamardFactor t s - 1) * XiFiniteHadamard spec s‖ := by
            rw [hfactor]
            ring_nf
    _ = ‖hadamardFactor t s - 1‖ * ‖XiFiniteHadamard spec s‖ := by simp [norm_mul]
    _ ≤ ‖hadamardFactor t s - 1‖ * ‖XiFiniteHadamard spec s‖ := le_rfl

/-- One-step increment bound with explicit second-order Hadamard control
whenever `‖s / criticalLinePoint t‖ ≤ 1`. -/
theorem norm_XiFiniteHadamard_insert_sub_le_three_mul_sq
    (spec : Finset ℝ) {t : ℝ} (ht : t ∉ spec) (s : ℂ)
    (hz : ‖s / criticalLinePoint t‖ ≤ 1) :
    ‖XiFiniteHadamard (insert t spec) s - XiFiniteHadamard spec s‖ ≤
      (3 * ‖s / criticalLinePoint t‖ ^ 2) * ‖XiFiniteHadamard spec s‖ := by
  have hstep := norm_XiFiniteHadamard_insert_sub_le spec ht s
  have hfac : ‖hadamardFactor t s - 1‖ ≤ 3 * ‖s / criticalLinePoint t‖ ^ 2 :=
    norm_hadamardFactor_sub_one_le_three_mul_sq t s hz
  exact le_trans hstep (mul_le_mul_of_nonneg_right hfac (norm_nonneg _))

/-- Finite Hadamard products vanish exactly at listed critical-line points. -/
theorem XiFiniteHadamard_zero_iff_exists
    (spec : Finset ℝ) (s : ℂ) :
    XiFiniteHadamard spec s = 0 ↔ ∃ t : ℝ, t ∈ spec ∧ s = criticalLinePoint t := by
  classical
  constructor
  · intro hs
    unfold XiFiniteHadamard at hs
    rcases (Finset.prod_eq_zero_iff.mp hs) with ⟨t, ht, hfac⟩
    exact ⟨t, ht, (hadamardFactor_eq_zero_iff t s).1 hfac⟩
  · rintro ⟨t, ht, hsEq⟩
    subst hsEq
    unfold XiFiniteHadamard
    refine Finset.prod_eq_zero_iff.mpr ?_
    refine ⟨t, ht, ?_⟩
    exact (hadamardFactor_eq_zero_iff t (criticalLinePoint t)).2 rfl

/-- Basic growth envelope for finite Xi products:
the norm is bounded by the product of factor envelopes. -/
theorem norm_XiFinite_le_factorized_envelope
    (spec : Finset ℝ) (s : ℂ) :
    ‖XiFinite spec s‖ ≤
      Finset.prod spec (fun t => (‖s‖ + ‖criticalLinePoint t‖)) := by
  classical
  unfold XiFinite
  refine Finset.induction_on spec ?h0 ?hstep
  · simp
  · intro a spec ha ih
    have henvelope_nonneg : 0 ≤ ‖s‖ + ‖criticalLinePoint a‖ := by
      exact add_nonneg (norm_nonneg _) (norm_nonneg _)
    calc
      ‖Finset.prod (insert a spec) (fun t => (s - criticalLinePoint t : ℂ))‖
          = ‖s - criticalLinePoint a‖ *
              ‖Finset.prod spec (fun t => (s - criticalLinePoint t : ℂ))‖ := by
            simp [Finset.prod_insert, ha]
      _ ≤ (‖s‖ + ‖criticalLinePoint a‖) *
            ‖Finset.prod spec (fun t => (s - criticalLinePoint t : ℂ))‖ := by
            exact mul_le_mul_of_nonneg_right (norm_sub_le s (criticalLinePoint a))
              (norm_nonneg _)
      _ ≤ (‖s‖ + ‖criticalLinePoint a‖) *
            Finset.prod spec (fun t => (‖s‖ + ‖criticalLinePoint t‖)) := by
            exact mul_le_mul_of_nonneg_left ih henvelope_nonneg
      _ = Finset.prod (insert a spec) (fun t => (‖s‖ + ‖criticalLinePoint t‖)) := by
            simp [Finset.prod_insert, ha]

/-- Cardinality growth corollary for finite Xi products under a uniform ordinate bound. -/
theorem norm_XiFinite_le_pow_card_of_ordinate_bound
    (spec : Finset ℝ) (s : ℂ) (B : ℝ)
    (hB : ∀ t ∈ spec, ‖criticalLinePoint t‖ ≤ B) :
    ‖XiFinite spec s‖ ≤ (‖s‖ + B) ^ spec.card := by
  have hEnvelope :
      ‖XiFinite spec s‖ ≤
        Finset.prod spec (fun t => (‖s‖ + ‖criticalLinePoint t‖)) :=
    norm_XiFinite_le_factorized_envelope spec s
  have hProdLeAux :
      ∀ s' : Finset ℝ,
        (∀ t ∈ s', ‖criticalLinePoint t‖ ≤ B) →
        Finset.prod s' (fun t => (‖s‖ + ‖criticalLinePoint t‖))
          ≤ Finset.prod s' (fun _t => (‖s‖ + B)) := by
    intro s'
    refine Finset.induction_on s' ?h0 ?hstep
    · intro _hBs
      simp
    · intro a s' ha ih hBs
      have hBa : ‖criticalLinePoint a‖ ≤ B := hBs a (by simp [ha])
      have hBs' : ∀ t ∈ s', ‖criticalLinePoint t‖ ≤ B := by
        intro t ht
        exact hBs t (by simp [ht])
      have hfactor : ‖s‖ + ‖criticalLinePoint a‖ ≤ ‖s‖ + B := by
        linarith [hBa]
      have hconst_nonneg : 0 ≤ ‖s‖ + B := by
        have hleft_nonneg : 0 ≤ ‖s‖ + ‖criticalLinePoint a‖ := by
          exact add_nonneg (norm_nonneg _) (norm_nonneg _)
        exact le_trans hleft_nonneg hfactor
      have hprod_nonneg :
          0 ≤ Finset.prod s' (fun t => (‖s‖ + ‖criticalLinePoint t‖)) := by
        refine Finset.prod_nonneg ?_
        intro t ht
        exact add_nonneg (norm_nonneg _) (norm_nonneg _)
      calc
        Finset.prod (insert a s') (fun t => (‖s‖ + ‖criticalLinePoint t‖))
            = (‖s‖ + ‖criticalLinePoint a‖) *
                Finset.prod s' (fun t => (‖s‖ + ‖criticalLinePoint t‖)) := by
                  simp [Finset.prod_insert, ha]
        _ ≤ (‖s‖ + B) * Finset.prod s' (fun t => (‖s‖ + ‖criticalLinePoint t‖)) := by
              exact mul_le_mul_of_nonneg_right hfactor hprod_nonneg
        _ ≤ (‖s‖ + B) * Finset.prod s' (fun _t => (‖s‖ + B)) := by
              exact mul_le_mul_of_nonneg_left (ih hBs') hconst_nonneg
        _ = (‖s‖ + B) ^ (s'.card + 1) := by
              simpa [pow_succ, mul_comm, mul_left_comm, mul_assoc]
        _ = Finset.prod (insert a s') (fun _t => (‖s‖ + B)) := by
              simpa [ha]
  have hProdLe :
      Finset.prod spec (fun t => (‖s‖ + ‖criticalLinePoint t‖))
        ≤ Finset.prod spec (fun _t => (‖s‖ + B)) :=
    hProdLeAux spec hB
  calc
    ‖XiFinite spec s‖
        ≤ Finset.prod spec (fun t => (‖s‖ + ‖criticalLinePoint t‖)) := hEnvelope
    _ ≤ Finset.prod spec (fun _t => (‖s‖ + B)) := hProdLe
    _ = (‖s‖ + B) ^ spec.card := by simp

/-- One-step factorization when inserting a new ordinate into the finite product. -/
theorem XiFinite_insert
    (spec : Finset ℝ) {t : ℝ} (ht : t ∉ spec) (s : ℂ) :
    XiFinite (insert t spec) s = (s - criticalLinePoint t) * XiFinite spec s := by
  simp [XiFinite, Finset.prod_insert, ht]

/-- One-step increment bound for the finite Xi products. -/
theorem norm_XiFinite_insert_sub_le
    (spec : Finset ℝ) {t : ℝ} (ht : t ∉ spec) (s : ℂ) :
    ‖XiFinite (insert t spec) s - XiFinite spec s‖ ≤
      ‖s - criticalLinePoint t - 1‖ *
        Finset.prod spec (fun u => (‖s‖ + ‖criticalLinePoint u‖)) := by
  have hfactor := XiFinite_insert spec ht s
  calc
    ‖XiFinite (insert t spec) s - XiFinite spec s‖
        = ‖((s - criticalLinePoint t) - 1) * XiFinite spec s‖ := by
            rw [hfactor]
            ring_nf
    _ = ‖s - criticalLinePoint t - 1‖ * ‖XiFinite spec s‖ := by
          simp [sub_eq_add_neg, add_assoc]
    _ ≤ ‖s - criticalLinePoint t - 1‖ *
          Finset.prod spec (fun u => (‖s‖ + ‖criticalLinePoint u‖)) := by
          exact mul_le_mul_of_nonneg_left
            (norm_XiFinite_le_factorized_envelope spec s)
            (norm_nonneg _)

theorem XiFinite_zero_of_mem
    (spec : Finset ℝ) {t : ℝ}
    (ht : t ∈ spec) :
    XiFinite spec (criticalLinePoint t) = 0 := by
  classical
  unfold XiFinite
  refine Finset.prod_eq_zero_iff.mpr ?_
  exact ⟨t, ht, by simp [criticalLinePoint]⟩

theorem XiFinite_zero_iff_exists
    (spec : Finset ℝ) (s : ℂ) :
    XiFinite spec s = 0 ↔ ∃ t : ℝ, t ∈ spec ∧ s = criticalLinePoint t := by
  classical
  constructor
  · intro hs
    unfold XiFinite at hs
    rcases (Finset.prod_eq_zero_iff.mp hs) with ⟨t, ht, hfac⟩
    refine ⟨t, ht, ?_⟩
    exact sub_eq_zero.mp hfac
  · rintro ⟨t, ht, rfl⟩
    exact XiFinite_zero_of_mem spec ht

/-- The finite Hadamard model has the same zero set as `XiFinite`. -/
theorem XiFiniteHadamard_zero_iff_XiFinite_zero
    (spec : Finset ℝ) (s : ℂ) :
    XiFiniteHadamard spec s = 0 ↔ XiFinite spec s = 0 := by
  constructor
  · intro hs
    rcases (XiFiniteHadamard_zero_iff_exists spec s).1 hs with ⟨t, ht, hsEq⟩
    exact (XiFinite_zero_iff_exists spec s).2 ⟨t, ht, hsEq⟩
  · intro hs
    rcases (XiFinite_zero_iff_exists spec s).1 hs with ⟨t, ht, hsEq⟩
    exact (XiFiniteHadamard_zero_iff_exists spec s).2 ⟨t, ht, hsEq⟩

/-- Exact finite bridge theorem for the explicit finite Xi model. -/
theorem finiteXi_spectralBridge
    (spec : Finset ℝ) :
    SpectralBridge (XiFinite spec) (finiteSpecSet spec) := by
  intro s
  constructor
  · intro hs
    rcases (XiFinite_zero_iff_exists spec s).1 hs with ⟨t, ht, hsEq⟩
    exact ⟨t, ht, hsEq⟩
  · intro hs
    rcases hs with ⟨t, ht, hsEq⟩
    rw [hsEq]
    exact XiFinite_zero_of_mem spec ht

/-- RH holds for the explicit finite Xi model by exact bridge. -/
theorem rh_XiFinite
    (spec : Finset ℝ) :
    RiemannHypothesisXi (XiFinite spec) := by
  exact bridge_implies_rh (XiFinite spec) (finiteSpecSet spec) (finiteXi_spectralBridge spec)

/-- Constant finite-level family used to instantiate transfer contracts for `XiFinite`. -/
def XiFiniteConst (spec : Finset ℝ) : ℕ → (ℂ → ℂ) := fun _ => XiFinite spec

/-- Constant spectral ladder used with `XiFiniteConst`. -/
def specConst (spec : Finset ℝ) : ℕ → Finset ℝ := fun _ => spec

/-- Zero tolerance profile. -/
def tolZero : ℕ → ℝ := fun _ => 0

theorem finiteBridgeFamily_XiFiniteConst
    (spec : Finset ℝ) :
    FiniteBridgeFamily (XiFiniteConst spec) (specConst spec) := by
  intro N
  simpa [XiFiniteConst, specConst, levelSpecSet, finiteSpecSet]
    using finiteXi_spectralBridge spec

theorem zeroTol_tolZero : zeroTol tolZero := by
  intro N
  simp [tolZero]

theorem approxZero_XiFiniteConst
    (spec : Finset ℝ) :
    ApproxZeroConvergence (XiFinite spec) (XiFiniteConst spec) tolZero := by
  intro s hs
  refine ⟨0, ?_⟩
  simp [XiFiniteConst, tolZero, hs]

theorem rigidity_XiFiniteConst
    (spec : Finset ℝ) :
    SpectralRigidity (XiFiniteConst spec) tolZero := by
  intro N s hs
  have hnorm0 : ‖XiFiniteConst spec N s‖ = 0 := by
    exact le_antisymm hs (norm_nonneg _)
  exact norm_eq_zero.mp hnorm0

/-- Explicit convergence-transfer contract for `XiFinite` (all obligations discharged). -/
def XiFiniteConvergenceContract
    (spec : Finset ℝ) :
    RHConvergenceTransferContract (XiFinite spec) where
  XiN := XiFiniteConst spec
  specN := specConst spec
  finiteBridge := finiteBridgeFamily_XiFiniteConst spec
  tol := tolZero
  tolNonneg := zeroTol_tolZero
  approxZero := approxZero_XiFiniteConst spec
  rigidity := rigidity_XiFiniteConst spec

/-- RH for `XiFinite` via convergence-transfer contract path. -/
theorem rh_XiFinite_via_convergence_contract
    (spec : Finset ℝ) :
    RiemannHypothesisXi (XiFinite spec) := by
  exact rh_of_convergence_transfer_contract (XiFinite spec) (XiFiniteConvergenceContract spec)

end

end Gutoe.RiemannFiniteXiModel
