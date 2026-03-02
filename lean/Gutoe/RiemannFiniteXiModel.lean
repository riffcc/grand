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
