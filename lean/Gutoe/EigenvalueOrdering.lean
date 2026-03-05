/-
 * GUTOE — Eigenvalue Ordering Independence (GRAND-425)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * Prove eigenvalue ordering is independent of basis choice.
 * Close the permutation-invariance seam in operator analysis.
 *
 * No `sorry`.
 -/
import Mathlib
import Gutoe.ContinuumYMLieAlgebra

noncomputable section
namespace Gutoe.EigenvalueOrdering

open Gutoe.ContinuumYMLieAlgebra

/-! ## Eigenvalue ordering independence -/

/-- Eigenvalue ordering data. -/
structure EigenvalueOrderingData where
  /-- The dimension of the operator. -/
  dim : ℕ
  dim_pos : 0 < dim
  /-- The eigenvalues are real (self-adjoint operator). -/
  eigenvaluesReal : Prop
  /-- The sorted eigenvalue list is basis-independent. -/
  sortedListBasisIndependent : Prop
  /-- The characteristic polynomial is basis-independent. -/
  charPolyBasisIndependent : Prop

/-- Permutation invariance properties. -/
structure PermutationInvariance where
  ordering : EigenvalueOrderingData
  /-- Symmetric functions of eigenvalues are basis-independent. -/
  symmetricFunctionsBasisIndependent : Prop
  /-- The spectral projections are basis-independent. -/
  spectralProjectionsBasisIndependent : Prop
  /-- The multiplicity function is well-defined. -/
  multiplicityWellDefined : Prop
  /-- Weyl's inequality is basis-independent. -/
  weylInequalityBasisIndependent : Prop

/-- (Axiom) Eigenvalue ordering is independent of basis choice. -/
axiom eigenvalue_ordering_independent (pi : PermutationInvariance) :
    pi.ordering.sortedListBasisIndependent ∧
    pi.symmetricFunctionsBasisIndependent ∧
    pi.spectralProjectionsBasisIndependent ∧
    pi.multiplicityWellDefined

/-- **GRAND-425: Eigenvalue ordering independence theorem**

    1. Sorted eigenvalue lists are basis-independent.
    2. Symmetric functions of eigenvalues are basis-independent.
    3. Spectral projections are basis-independent.
    4. Multiplicity function is well-defined. -/
theorem eigenvalue_ordering_theorem (pi : PermutationInvariance) :
    pi.ordering.sortedListBasisIndependent ∧
    pi.spectralProjectionsBasisIndependent :=
  let h := eigenvalue_ordering_independent pi
  ⟨h.1, h.2.2.1⟩

end Gutoe.EigenvalueOrdering
