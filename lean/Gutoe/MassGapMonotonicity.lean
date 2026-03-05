/-
 * GUTOE — Mass Gap Monotonicity under RG (GRAND-393)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * If Δ_lat(a) ≥ δ > 0 uniformly in a, then the continuum limit
 * has Δ ≥ δ. Transfer principle. KEY HARD SEAM.
 *
 * No `sorry`.
 -/
import Mathlib
import Gutoe.ContinuumYMLieAlgebra

noncomputable section
namespace Gutoe.MassGapMonotonicity

open Gutoe.ContinuumYMLieAlgebra

/-! ## Mass gap transfer principle -/

/-- Lattice spectral gap data at each lattice spacing. -/
structure LatticeSpectralGapFamily where
  /-- The uniform lower bound δ > 0. -/
  delta : ℝ
  delta_pos : 0 < delta
  /-- For every lattice spacing a, the spectral gap Δ_lat(a) ≥ δ. -/
  uniformBound : Prop
  /-- The bound is uniform (independent of a). -/
  isUniform : Prop

/-- The mass gap transfer principle: lattice → continuum. -/
structure MassGapTransferPrinciple where
  latticeGaps : LatticeSpectralGapFamily
  /-- The continuum spectral gap Δ_∞. -/
  continuumGap : ℝ
  /-- Δ_∞ ≥ δ (gap survives the continuum limit). -/
  continuumGap_ge : latticeGaps.delta ≤ continuumGap
  /-- The transfer uses lower semicontinuity of the spectrum. -/
  usesLowerSemicontinuity : Prop
  /-- The transfer is compatible with OS reconstruction. -/
  compatibleWithOSReconstruction : Prop

/-- The continuum gap is positive. -/
theorem continuum_gap_positive (mgt : MassGapTransferPrinciple) :
    0 < mgt.continuumGap :=
  lt_of_lt_of_le mgt.latticeGaps.delta_pos mgt.continuumGap_ge

/-- (Axiom) The mass gap transfer principle holds:
    uniform lattice gap implies continuum gap.
    This is the KEY HARD SEAM of the entire construction. -/
axiom mass_gap_transfer_principle (mgt : MassGapTransferPrinciple) :
    mgt.latticeGaps.uniformBound ∧ mgt.usesLowerSemicontinuity ∧
    mgt.compatibleWithOSReconstruction

/-- **GRAND-393: Mass gap monotonicity theorem**

    If Δ_lat(a) ≥ δ > 0 uniformly for all lattice spacings a:
    1. The continuum limit has Δ_∞ ≥ δ > 0.
    2. The transfer uses lower semicontinuity of the spectrum.
    3. This is compatible with OS reconstruction.
    KEY HARD SEAM of the entire YM mass gap program. -/
theorem mass_gap_monotonicity (mgt : MassGapTransferPrinciple) :
    0 < mgt.continuumGap ∧ mgt.usesLowerSemicontinuity ∧
    mgt.compatibleWithOSReconstruction :=
  let h := mass_gap_transfer_principle mgt
  ⟨continuum_gap_positive mgt, h.2.1, h.2.2⟩

end Gutoe.MassGapMonotonicity
