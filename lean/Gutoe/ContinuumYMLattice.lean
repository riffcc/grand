/-
 * GUTOE — Continuum YM Hypercubic Lattice (GRAND-416)
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * Top-down formalization of the standard Wilson hypercubic lattice:
 *   Λ_a = aℤ⁴.
 -/

import Mathlib

namespace Gutoe.ContinuumYMLattice

open scoped BigOperators

/-- Integer lattice coordinates for sites (`ℤ⁴`). -/
abbrev LatticeSite := Fin 4 → ℤ

/-- Real embedding of a site into physical coordinates with spacing `a`.
    This realizes `x = a n` with `n ∈ ℤ⁴`. -/
def siteEmbedding (a : ℝ) (x : LatticeSite) : Fin 4 → ℝ :=
  fun μ => a * (x μ : ℝ)

/-- The geometric subset `aℤ⁴ ⊂ ℝ⁴`. -/
def scaledIntegerLattice (a : ℝ) : Set (Fin 4 → ℝ) :=
  {y | ∃ x : LatticeSite, y = siteEmbedding a x}

/-- Hypercubic lattice sites at spacing `a`.
    We keep the canonical coordinate model `ℤ⁴`; `siteEmbedding` maps it into `aℤ⁴`. -/
abbrev HypercubicLattice (a : ℝ) (ha : 0 < a) : Type := LatticeSite

/-- Positive spacing implies non-degenerate embedding scale. -/
lemma spacing_ne_zero {a : ℝ} (ha : 0 < a) : a ≠ 0 := by linarith

/-- The canonical embedding of lattice sites into `aℤ⁴` is injective. -/
theorem siteEmbedding_injective (a : ℝ) (ha : 0 < a) :
    Function.Injective (siteEmbedding a) := by
  intro x y hxy
  funext μ
  have hμ : a * (x μ : ℝ) = a * (y μ : ℝ) := by
    exact congrArg (fun f => f μ) hxy
  have hcast : (x μ : ℝ) = (y μ : ℝ) := by
    exact mul_left_cancel₀ (spacing_ne_zero ha) hμ
  exact Int.cast_inj.mp hcast

/-- The image of `siteEmbedding` is exactly `aℤ⁴`. -/
theorem range_siteEmbedding_eq_scaledIntegerLattice (a : ℝ) :
    Set.range (siteEmbedding a) = scaledIntegerLattice a := by
  ext y
  constructor
  · intro hy
    rcases hy with ⟨x, rfl⟩
    exact ⟨x, rfl⟩
  · intro hy
    rcases hy with ⟨x, rfl⟩
    exact ⟨x, rfl⟩

/-- Orientation of a nearest-neighbor lattice link. -/
inductive LinkOrientation
  | forward
  | backward
  deriving DecidableEq, Repr

/-- Unit basis vector `e_μ` in `ℤ⁴`. -/
def unitVec (μ : Fin 4) : LatticeSite :=
  fun ν => if ν = μ then 1 else 0

@[simp] theorem unitVec_same (μ : Fin 4) : unitVec μ μ = 1 := by
  simp [unitVec]

@[simp] theorem unitVec_ne (μ ν : Fin 4) (h : ν ≠ μ) : unitVec μ ν = 0 := by
  simp [unitVec, h]

/-- Signed nearest-neighbor step in direction `μ`. -/
def orientedStep (μ : Fin 4) : LinkOrientation → LatticeSite
  | .forward => unitVec μ
  | .backward => -unitVec μ

/-- Nearest-neighbor relation on sites in `ℤ⁴`. -/
def IsNearestNeighbor (x y : LatticeSite) : Prop :=
  ∃ μ : Fin 4, y = x + unitVec μ ∨ y = x - unitVec μ

/-- Oriented nearest-neighbor edge. -/
structure LatticeLink where
  src : LatticeSite
  dir : Fin 4
  orient : LinkOrientation
  deriving DecidableEq, Repr

/-- Target endpoint of an oriented link. -/
def LatticeLink.dst (ℓ : LatticeLink) : LatticeSite :=
  ℓ.src + orientedStep ℓ.dir ℓ.orient

/-- Every `LatticeLink` connects nearest neighbors. -/
theorem LatticeLink.isNearestNeighbor (ℓ : LatticeLink) :
    IsNearestNeighbor ℓ.src ℓ.dst := by
  rcases ℓ with ⟨src, dir, orient⟩
  cases orient with
  | forward =>
      refine ⟨dir, ?_⟩
      left
      rfl
  | backward =>
      refine ⟨dir, ?_⟩
      right
      simp [LatticeLink.dst, orientedStep, sub_eq_add_neg]

/-- Elementary oriented plaquette in the `(μ, ν)` plane with `μ ≠ ν`. -/
structure LatticePlaquette where
  base : LatticeSite
  μ : Fin 4
  ν : Fin 4
  hμν : μ ≠ ν
  deriving DecidableEq, Repr

/-- First edge of the plaquette boundary. -/
def LatticePlaquette.link1 (p : LatticePlaquette) : LatticeLink :=
  { src := p.base, dir := p.μ, orient := .forward }

/-- Second edge of the plaquette boundary. -/
def LatticePlaquette.link2 (p : LatticePlaquette) : LatticeLink :=
  { src := p.base + unitVec p.μ, dir := p.ν, orient := .forward }

/-- Third edge of the plaquette boundary. -/
def LatticePlaquette.link3 (p : LatticePlaquette) : LatticeLink :=
  { src := p.base + unitVec p.μ + unitVec p.ν, dir := p.μ, orient := .backward }

/-- Fourth edge of the plaquette boundary. -/
def LatticePlaquette.link4 (p : LatticePlaquette) : LatticeLink :=
  { src := p.base + unitVec p.ν, dir := p.ν, orient := .backward }

/-- The boundary links of an elementary plaquette as an ordered 4-tuple. -/
def LatticePlaquette.boundaryLinks (p : LatticePlaquette) : Fin 4 → LatticeLink
  | ⟨0, _⟩ => p.link1
  | ⟨1, _⟩ => p.link2
  | ⟨2, _⟩ => p.link3
  | ⟨3, _⟩ => p.link4

/-- A plaquette is a closed oriented unit square. -/
theorem LatticePlaquette.boundary_closes (p : LatticePlaquette) :
    p.link4.dst = p.base := by
  ext κ
  simp [LatticePlaquette.link4, LatticeLink.dst, orientedStep, unitVec, sub_eq_add_neg]

/-- Number of distinct plaquette planes through each site in 4D. -/
def plaquetteCountPerSite : ℕ := Nat.choose 4 2

/-- In 4D, each site belongs to exactly six elementary plaquette planes. -/
theorem plaquette_count_per_site : plaquetteCountPerSite = 6 := by
  decide

/-- Integer-valued `k`-chains with finite support. -/
abbrev SiteChain := LatticeSite →₀ ℤ
abbrev LinkChain := LatticeLink →₀ ℤ
abbrev PlaquetteChain := LatticePlaquette →₀ ℤ

/-- Boundary of an oriented 1-cell (link): `∂[x→y] = y - x`. -/
def boundaryLink (ℓ : LatticeLink) : SiteChain :=
  Finsupp.single ℓ.dst 1 - Finsupp.single ℓ.src 1

/-- Boundary of an elementary 2-cell (plaquette) as oriented sum of its four links. -/
def boundaryPlaquette (p : LatticePlaquette) : LinkChain :=
  Finsupp.single p.link1 1 +
  Finsupp.single p.link2 1 +
  Finsupp.single p.link3 1 +
  Finsupp.single p.link4 1

/-- Boundary operator on 1-chains: `∂₁ : C₁ → C₀`. -/
def boundary₁ (c : LinkChain) : SiteChain :=
  c.sum (fun ℓ n => n • boundaryLink ℓ)

/-- Boundary operator on 2-chains: `∂₂ : C₂ → C₁`. -/
def boundary₂ (c : PlaquetteChain) : LinkChain :=
  c.sum (fun p n => n • boundaryPlaquette p)

notation "∂₁" => boundary₁
notation "∂₂" => boundary₂

/-- `Λ_a` has the same underlying site topology as `ℤ⁴`. -/
def hypercubic_homeomorph_Z4 (a : ℝ) (ha : 0 < a) :
    HypercubicLattice a ha ≃ₜ (Fin 4 → ℤ) :=
  Homeomorph.refl _

/-- The hypercubic lattice is discrete (as for `ℤ⁴`). -/
theorem hypercubic_discrete (a : ℝ) (ha : 0 < a) :
    DiscreteTopology (HypercubicLattice a ha) := by
  infer_instance

/-- The full lattice has no topological boundary. -/
theorem hypercubic_no_boundary (a : ℝ) (ha : 0 < a) :
    frontier (Set.univ : Set (HypercubicLattice a ha)) = ∅ := by
  simp

end Gutoe.ContinuumYMLattice
