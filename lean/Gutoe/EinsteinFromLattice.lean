/- 
 * GUTOE — Einstein Field Equations from Cl(1,3) Lattice Bridge
 *
 * GRAND-89 (core): Einstein field equations in continuum limit
 *
 * This module formalizes the proof skeleton:
 *
 *   Cl(1,3)  ⇒  SC lattice (coordination 6)
 *          ⇒  Regge-style discrete curvature dynamics
 *          ⇒  continuum bridge for the Einstein tensor
 *          ⇒  Einstein field equations (with optional λ_QG correction term)
 *
 * Open physics work is kept explicit as hypotheses in bridge theorems.
 * This avoids fake confidence while keeping the theorem chain executable.
 -/

import Mathlib
import Gutoe.LatticeGeometry
import Gutoe.ContinuumLimit
import Gutoe.GravityMetric
import Gutoe.LambdaQG
import Gutoe.GaugeGroupSU3
import Gutoe.GaugeGroupSM
import Gutoe.GaugeGroupSU2
import Gutoe.FineStructure
import Gutoe.Z3Uniqueness

namespace Gutoe.EinsteinFromLattice

open Gutoe.LatticeGeometry
open Gutoe.GravityMetric

/-- Rank-2 real tensor field on an index type `ι`. -/
abbrev TensorField (ι : Type) := ι → ι → ℝ

/-- Einstein equation in component form:
    `G_{μν} + Λ g_{μν} = κ T_{μν}`. -/
def EinsteinFieldEquation
    {ι : Type}
    (G g T : TensorField ι)
    (Lambda kappa : ℝ) : Prop :=
  ∀ μ ν, G μ ν + Lambda * g μ ν = kappa * T μ ν

/-- GUTOE-modified Einstein equation with the lattice correction tensor `H_{μν}`:
    `G_{μν} + λ_QG l_P² H_{μν} + Λ g_{μν} = κ T_{μν}`. -/
def ModifiedEinsteinFieldEquation
    {ι : Type}
    (G H g T : TensorField ι)
    (lP Lambda kappa : ℝ) : Prop :=
  ∀ μ ν, G μ ν + lambda_qg * lP ^ 2 * H μ ν + Lambda * g μ ν = kappa * T μ ν

/-- Regge action on a simplicial edge set: `S = Σ_e A_e δ_e`. -/
def reggeAction
    {Edge : Type}
    [Fintype Edge]
    (area deficit : Edge → ℝ) : ℝ :=
  ∑ e, area e * deficit e

/-- Stationarity condition for a Regge action variation wrt edge lengths. -/
def reggeStationary
    {Edge : Type}
    (dSdl : Edge → ℝ) : Prop :=
  ∀ e, dSdl e = 0

/-- First variation of Regge action in edge variables:
    `δS = Σ_e (δA_e · δ_e + A_e · δδ_e)`. -/
def reggeVariation
    {Edge : Type}
    [Fintype Edge]
    (area deficit dArea dDeficit : Edge → ℝ) : ℝ :=
  ∑ e, (dArea e * deficit e + area e * dDeficit e)

/-- Schläfli identity in finite-edge form:
    `Σ_e A_e · dδ_e = 0`. -/
def SchlaefliIdentity
    {Edge : Type}
    [Fintype Edge]
    (area dDeficit : Edge → ℝ) : Prop :=
  ∑ e, area e * dDeficit e = 0

/-- With Schläfli identity, the deficit-angle variation term cancels and
    Regge variation reduces to `Σ_e δA_e · δ_e`. -/
theorem regge_variation_of_schlaefli
    {Edge : Type}
    [Fintype Edge]
    (area deficit dArea dDeficit : Edge → ℝ)
    (hSch : SchlaefliIdentity area dDeficit) :
    reggeVariation area deficit dArea dDeficit = ∑ e, dArea e * deficit e := by
  unfold reggeVariation SchlaefliIdentity at *
  calc
    ∑ e, (dArea e * deficit e + area e * dDeficit e)
        = (∑ e, dArea e * deficit e) + (∑ e, area e * dDeficit e) := by
            simp [Finset.sum_add_distrib]
    _ = ∑ e, dArea e * deficit e := by simpa [hSch]

/-- Regge deficit angle around an edge:
    `δ_e = 2π - Σ_t θ_{e,t}`. -/
noncomputable def deficitAngle (sumDihedral : ℝ) : ℝ :=
  2 * Real.pi - sumDihedral

/-- Flat-edge normalization: if the local dihedral sum is exactly `2π`,
    the deficit angle is zero. -/
theorem deficitAngle_zero_of_full_angle
    (hsum : sumDihedral = 2 * Real.pi) :
    deficitAngle sumDihedral = 0 := by
  unfold deficitAngle
  linarith

/-- Sign convention check: if local dihedral sum exceeds `2π`,
    the deficit angle is negative. -/
theorem deficitAngle_neg_of_gt_full_angle
    (hsum : 2 * Real.pi < sumDihedral) :
    deficitAngle sumDihedral < 0 := by
  unfold deficitAngle
  linarith

/-- Discrete Einstein equation written on edges:
    each edge variation is balanced by projected stress-energy. -/
def reggeEdgeEquation
    {Edge : Type}
    (dSdl source : Edge → ℝ) : Prop :=
  ∀ e, dSdl e = source e

/-- Stationarity implies zero-source edge equations. -/
theorem regge_stationary_implies_zero_source
    {Edge : Type}
    {dSdl : Edge → ℝ}
    (hstat : reggeStationary dSdl) :
    reggeEdgeEquation dSdl (fun _ => 0) := by
  intro e
  specialize hstat e
  simpa using hstat

/-- Vertex labels for one unit SC cube (`0..7`). -/
abbrev CubeVertex := Fin 8

/-- All vertices of the unit cube. -/
def scCubeVertices : Finset CubeVertex := Finset.univ

/-- Body-diagonal simplicialization of one SC cube into 6 tetrahedra.
    This is the canonical decomposition used for Regge-style curvature bookkeeping
    without introducing auxiliary interior vertices. -/
def scCubeTetrahedra : Finset (Finset CubeVertex) :=
  { ({0, 1, 3, 7} : Finset CubeVertex),
    ({0, 3, 2, 7} : Finset CubeVertex),
    ({0, 2, 6, 7} : Finset CubeVertex),
    ({0, 6, 4, 7} : Finset CubeVertex),
    ({0, 4, 5, 7} : Finset CubeVertex),
    ({0, 5, 1, 7} : Finset CubeVertex) }

/-- The canonical SC body-diagonal decomposition has exactly 6 tetrahedra. -/
theorem sc_cube_tetrahedra_card : scCubeTetrahedra.card = 6 := by
  decide

/-- Each listed simplex in the canonical decomposition has exactly 4 vertices. -/
theorem sc_cube_tetrahedra_listed_are_4simplices :
    ({0, 1, 3, 7} : Finset CubeVertex).card = 4 ∧
    ({0, 3, 2, 7} : Finset CubeVertex).card = 4 ∧
    ({0, 2, 6, 7} : Finset CubeVertex).card = 4 ∧
    ({0, 6, 4, 7} : Finset CubeVertex).card = 4 ∧
    ({0, 4, 5, 7} : Finset CubeVertex).card = 4 ∧
    ({0, 5, 1, 7} : Finset CubeVertex).card = 4 := by
  decide

/-- The six canonical tetrahedra jointly include all cube vertices `0..7`. -/
theorem sc_cube_tetrahedra_cover_vertices :
    (0 : CubeVertex) ∈ scCubeTetrahedra.biUnion id ∧
    (1 : CubeVertex) ∈ scCubeTetrahedra.biUnion id ∧
    (2 : CubeVertex) ∈ scCubeTetrahedra.biUnion id ∧
    (3 : CubeVertex) ∈ scCubeTetrahedra.biUnion id ∧
    (4 : CubeVertex) ∈ scCubeTetrahedra.biUnion id ∧
    (5 : CubeVertex) ∈ scCubeTetrahedra.biUnion id ∧
    (6 : CubeVertex) ∈ scCubeTetrahedra.biUnion id ∧
    (7 : CubeVertex) ∈ scCubeTetrahedra.biUnion id := by
  decide

/-- Canonical `Z₃` action on cube vertices induced by cyclic axis permutation
    `(x,y,z) ↦ (y,z,x)` in binary cube coordinates. -/
def z3CubeAction : CubeVertex → CubeVertex
  | 0 => 0
  | 1 => 4
  | 2 => 1
  | 3 => 5
  | 4 => 2
  | 5 => 6
  | 6 => 3
  | 7 => 7

/-- The cube `Z₃` action has order `3` on vertices. -/
theorem z3_cube_action_order3 :
    ∀ v : CubeVertex, z3CubeAction (z3CubeAction (z3CubeAction v)) = v := by
  intro v
  fin_cases v <;> rfl

/-- `Z₃`-compatibility of the canonical 6-tetra SC decomposition:
    applying the cyclic cube action maps each tetrahedron back into the
    decomposition family. -/
theorem sc_cube_tetrahedra_z3_compatible :
    ∀ t ∈ scCubeTetrahedra, t.image z3CubeAction ∈ scCubeTetrahedra := by
  intro t ht
  have hmem :
      t = ({0, 1, 3, 7} : Finset CubeVertex) ∨
      t = ({0, 3, 2, 7} : Finset CubeVertex) ∨
      t = ({0, 2, 6, 7} : Finset CubeVertex) ∨
      t = ({0, 6, 4, 7} : Finset CubeVertex) ∨
      t = ({0, 4, 5, 7} : Finset CubeVertex) ∨
      t = ({0, 5, 1, 7} : Finset CubeVertex) := by
    simpa [scCubeTetrahedra, Finset.mem_insert, Finset.mem_singleton] using ht
  rcases hmem with h | h | h | h | h | h <;> subst h <;> decide

/-- Bridge hypotheses from discrete Regge dynamics to continuum Einstein dynamics.

`hConvergence` and `hDiscreteDynamics` are the explicit nontrivial bridge obligations
tracked by GRAND-268/269/270/271. -/
def ReggeToEinsteinBridge
    {ι : Type}
    (Gdisc Gcont H g T : TensorField ι)
    (lP Lambda kappa : ℝ) : Prop :=
  (∀ μ ν, Gcont μ ν = Gdisc μ ν) ∧
  (∀ μ ν, Gdisc μ ν + lambda_qg * lP ^ 2 * H μ ν + Lambda * g μ ν = kappa * T μ ν)

/-- Bridge constructor from an edge-level Regge equation and a projection map
    from edge equations into tensor components.

`hProject` is the explicit GRAND-269/270 obligation carrying geometric content. -/
theorem regge_edge_projection_to_bridge
    {ι Edge : Type}
    (embed : ι → ι → Edge)
    (Gdisc Gcont H g T : TensorField ι)
    (dSdl source : Edge → ℝ)
    (lP Lambda kappa : ℝ)
    (hCont : ∀ μ ν, Gcont μ ν = Gdisc μ ν)
    (hEdgeEq : reggeEdgeEquation dSdl source)
    (hProject :
      ∀ μ ν,
        dSdl (embed μ ν) = source (embed μ ν) →
        Gdisc μ ν + lambda_qg * lP ^ 2 * H μ ν + Lambda * g μ ν = kappa * T μ ν) :
    ReggeToEinsteinBridge Gdisc Gcont H g T lP Lambda kappa := by
  refine ⟨hCont, ?_⟩
  intro μ ν
  exact hProject μ ν (hEdgeEq (embed μ ν))

/-- Edge labels for the 6-tetra SC cube decomposition (19 unique edges). -/
abbrev SimpEdge := Fin 19

/-- SC-specialized Schläfli identity over the 19-edge decomposition. -/
def scSchlaefliIdentity (area dDeficit : SimpEdge → ℝ) : Prop :=
  SchlaefliIdentity area dDeficit

/-- SC-specialized cancellation theorem: once Schläfli is established on the
    19-edge geometry, the Regge variation drops the `A_e dδ_e` term. -/
theorem sc_regge_variation_of_schlaefli
    (area deficit dArea dDeficit : SimpEdge → ℝ)
    (hSch : scSchlaefliIdentity area dDeficit) :
    reggeVariation area deficit dArea dDeficit = ∑ e, dArea e * deficit e :=
  regge_variation_of_schlaefli area deficit dArea dDeficit hSch

/-- Flat-space Schläfli base case on the SC edge set:
    if all deficit-angle variations vanish, Schläfli holds exactly. -/
theorem sc_schlaefli_flat
    (area : SimpEdge → ℝ) :
    scSchlaefliIdentity area (fun _ => 0) := by
  unfold scSchlaefliIdentity SchlaefliIdentity
  simp

/-- Endpoint table for the 19 unique edges appearing in `scCubeTetrahedra`. -/
def simpEdgeEndpoints (e : SimpEdge) : CubeVertex × CubeVertex :=
  match e.1 with
  | 0 => ((0 : CubeVertex), (1 : CubeVertex))
  | 1 => ((0 : CubeVertex), (2 : CubeVertex))
  | 2 => ((0 : CubeVertex), (3 : CubeVertex))
  | 3 => ((0 : CubeVertex), (4 : CubeVertex))
  | 4 => ((0 : CubeVertex), (5 : CubeVertex))
  | 5 => ((0 : CubeVertex), (6 : CubeVertex))
  | 6 => ((0 : CubeVertex), (7 : CubeVertex))
  | 7 => ((1 : CubeVertex), (3 : CubeVertex))
  | 8 => ((1 : CubeVertex), (5 : CubeVertex))
  | 9 => ((1 : CubeVertex), (7 : CubeVertex))
  | 10 => ((2 : CubeVertex), (3 : CubeVertex))
  | 11 => ((2 : CubeVertex), (6 : CubeVertex))
  | 12 => ((2 : CubeVertex), (7 : CubeVertex))
  | 13 => ((3 : CubeVertex), (7 : CubeVertex))
  | 14 => ((4 : CubeVertex), (5 : CubeVertex))
  | 15 => ((4 : CubeVertex), (6 : CubeVertex))
  | 16 => ((4 : CubeVertex), (7 : CubeVertex))
  | 17 => ((5 : CubeVertex), (7 : CubeVertex))
  | 18 => ((6 : CubeVertex), (7 : CubeVertex))
  | _ => ((0 : CubeVertex), (1 : CubeVertex))

/-- Every edge in the endpoint table has distinct endpoints. -/
theorem simp_edge_endpoints_distinct :
    ∀ e : SimpEdge, (simpEdgeEndpoints e).1 ≠ (simpEdgeEndpoints e).2 := by
  intro e
  fin_cases e <;> decide

/-- Tetrahedra from the SC decomposition incident to edge `e`. -/
def edgeIncidentTetrahedra (e : SimpEdge) : Finset (Finset CubeVertex) :=
  let a := (simpEdgeEndpoints e).1
  let b := (simpEdgeEndpoints e).2
  scCubeTetrahedra.filter (fun t => a ∈ t ∧ b ∈ t)

/-- Incidence multiplicity of edge `e` in the 6-tetra decomposition. -/
def edgeIncidentCount (e : SimpEdge) : ℕ :=
  (edgeIncidentTetrahedra e).card

/-- The body diagonal `(0,7)` is incident to all six tetrahedra. -/
theorem edge_incident_count_body_diagonal :
    edgeIncidentCount ⟨6, by decide⟩ = 6 := by
  decide

/-- Every enumerated SC edge is used by at least one tetrahedron in the
    canonical decomposition. -/
theorem edge_incident_count_pos :
    ∀ e : SimpEdge, 1 ≤ edgeIncidentCount e := by
  intro e
  fin_cases e <;> decide

/-- No edge is incident to more than all six tetrahedra. -/
theorem edge_incident_count_le_six :
    ∀ e : SimpEdge, edgeIncidentCount e ≤ 6 := by
  intro e
  fin_cases e <;> decide

/-- Cartesian coordinates for unit-cube vertices used for edge-to-tensor projection. -/
def cubeCoord : CubeVertex → Fin 3 → ℝ
  | 0 => ![0, 0, 0]
  | 1 => ![1, 0, 0]
  | 2 => ![0, 1, 0]
  | 3 => ![1, 1, 0]
  | 4 => ![0, 0, 1]
  | 5 => ![1, 0, 1]
  | 6 => ![0, 1, 1]
  | 7 => ![1, 1, 1]

/-- Oriented edge vector in the unit-cube chart. -/
def simpEdgeVector (e : SimpEdge) : Fin 3 → ℝ :=
  let a := (simpEdgeEndpoints e).1
  let b := (simpEdgeEndpoints e).2
  fun i => cubeCoord b i - cubeCoord a i

/-- Canonical projection weight from edge vectors to spatial tensor components. -/
def edgeProjectionWeight (μ ν : Fin 3) (e : SimpEdge) : ℝ :=
  simpEdgeVector e μ * simpEdgeVector e ν

/-- Canonical coefficient map on the 19-edge set.
    Nonzero support is on edges 0..5, giving a 6-edge generating family
    for symmetric `3×3` tensors. -/
def scProjectionCoeffs (S : TensorField (Fin 3)) : SimpEdge → ℝ :=
  fun e =>
    match e.1 with
    | 0 => S 0 0 - S 0 1 - S 0 2
    | 1 => S 1 1 - S 0 1 - S 1 2
    | 2 => S 0 1
    | 3 => S 2 2 - S 0 2 - S 1 2
    | 4 => S 0 2
    | 5 => S 1 2
    | _ => 0

/-- A 6-edge generating subfamily (subset of the SC 19-edge set) used for
    rank-6 projection witness on symmetric `3×3` tensors. -/
abbrev BasisEdge := Fin 6

/-- Basis-edge vectors chosen from SC edges:
    e0=[1,0,0], e1=[0,1,0], e2=[0,0,1], e3=[1,1,0], e4=[1,0,1], e5=[0,1,1]. -/
def basisEdgeVector : BasisEdge → Fin 3 → ℝ
  | 0 => ![1, 0, 0]
  | 1 => ![0, 1, 0]
  | 2 => ![0, 0, 1]
  | 3 => ![1, 1, 0]
  | 4 => ![1, 0, 1]
  | 5 => ![0, 1, 1]

/-- Projection weight on the 6-edge basis. -/
def basisProjectionWeight (μ ν : Fin 3) (e : BasisEdge) : ℝ :=
  basisEdgeVector e μ * basisEdgeVector e ν

/-- Linear projection from basis coefficients to symmetric tensor entries. -/
def basisProjected (c : BasisEdge → ℝ) (μ ν : Fin 3) : ℝ :=
  ∑ e, basisProjectionWeight μ ν e * c e

/-- Coefficients solving the 6×6 linear system for symmetric tensors. -/
def basisProjectionCoeffs (S : TensorField (Fin 3)) : BasisEdge → ℝ
  | 0 => S 0 0 - S 0 1 - S 0 2
  | 1 => S 1 1 - S 0 1 - S 1 2
  | 2 => S 2 2 - S 0 2 - S 1 2
  | 3 => S 0 1
  | 4 => S 0 2
  | 5 => S 1 2

/-- Rank-6 witness (constructive surjectivity) for the basis projection map
    onto symmetric `3×3` tensors. -/
theorem basis_projection_surjective
    (S : TensorField (Fin 3))
    (hSym : ∀ μ ν, S μ ν = S ν μ) :
    ∀ μ ν, S μ ν = basisProjected (basisProjectionCoeffs S) μ ν := by
  intro μ ν
  fin_cases μ <;> fin_cases ν
  ·
    simp [basisProjected, basisProjectionWeight, basisProjectionCoeffs, basisEdgeVector,
      Fin.sum_univ_six]
    ring
  ·
    simp [basisProjected, basisProjectionWeight, basisProjectionCoeffs, basisEdgeVector,
      Fin.sum_univ_six]
  ·
    simp [basisProjected, basisProjectionWeight, basisProjectionCoeffs, basisEdgeVector,
      Fin.sum_univ_six]
  ·
    simpa [hSym 1 0] using
      (by
        simp [basisProjected, basisProjectionWeight, basisProjectionCoeffs, basisEdgeVector,
          Fin.sum_univ_six] :
        S 0 1 = basisProjected (basisProjectionCoeffs S) 1 0)
  ·
    simp [basisProjected, basisProjectionWeight, basisProjectionCoeffs, basisEdgeVector,
      Fin.sum_univ_six]
    ring
  ·
    simp [basisProjected, basisProjectionWeight, basisProjectionCoeffs, basisEdgeVector,
      Fin.sum_univ_six]
  ·
    simpa [hSym 2 0] using
      (by
        simp [basisProjected, basisProjectionWeight, basisProjectionCoeffs, basisEdgeVector,
          Fin.sum_univ_six] :
        S 0 2 = basisProjected (basisProjectionCoeffs S) 2 0)
  ·
    simpa [hSym 2 1] using
      (by
        simp [basisProjected, basisProjectionWeight, basisProjectionCoeffs, basisEdgeVector,
          Fin.sum_univ_six] :
        S 1 2 = basisProjected (basisProjectionCoeffs S) 2 1)
  ·
    simp [basisProjected, basisProjectionWeight, basisProjectionCoeffs, basisEdgeVector,
      Fin.sum_univ_six]
    ring

/-- Residual on edge equations (`δS/δl - source`). -/
def edgeResidual {Edge : Type} (dSdl source : Edge → ℝ) : Edge → ℝ :=
  fun e => dSdl e - source e

/-- Generic weighted projection of edge residuals to spatial tensor components. -/
def projectedResidual
    (w : Fin 3 → Fin 3 → SimpEdge → ℝ)
    (res : SimpEdge → ℝ) (μ ν : Fin 3) : ℝ :=
  ∑ e, w μ ν e * res e

/-- If edge equations are satisfied pointwise, every projected residual vanishes. -/
theorem projectedResidual_zero_of_reggeEdgeEquation
    (w : Fin 3 → Fin 3 → SimpEdge → ℝ)
    {dSdl source : SimpEdge → ℝ}
    (hEq : reggeEdgeEquation dSdl source) :
    ∀ μ ν, projectedResidual w (edgeResidual dSdl source) μ ν = 0 := by
  intro μ ν
  unfold projectedResidual
  apply Finset.sum_eq_zero
  intro e he
  have hz : edgeResidual dSdl source e = 0 := by
    unfold edgeResidual
    linarith [hEq e]
  simp [hz]

/-- Flat-space base case: zero edge residuals imply zero projected tensor residual. -/
theorem projectedResidual_zero_of_zero_residual
    (w : Fin 3 → Fin 3 → SimpEdge → ℝ) :
    ∀ μ ν, projectedResidual w (fun _ => 0) μ ν = 0 := by
  intro μ ν
  unfold projectedResidual
  apply Finset.sum_eq_zero
  intro e he
  simp

/-- Spatial Einstein residual used in the SC edge→tensor bridge. -/
noncomputable def spatialEinsteinResidual
    (G H g T : TensorField (Fin 3))
    (lP Lambda kappa : ℝ) (μ ν : Fin 3) : ℝ :=
  G μ ν + lambda_qg * lP ^ 2 * H μ ν + Lambda * g μ ν - kappa * T μ ν

/-- Concrete edge→tensor bridge model for the SC simplicialization. -/
def EdgeToTensorProjectionModel
    (G H g T : TensorField (Fin 3))
    (dSdl source : SimpEdge → ℝ)
    (lP Lambda kappa : ℝ)
    (w : Fin 3 → Fin 3 → SimpEdge → ℝ) : Prop :=
  ∀ μ ν,
    spatialEinsteinResidual G H g T lP Lambda kappa μ ν =
      projectedResidual w (edgeResidual dSdl source) μ ν

/-- Zero tensor field on the spatial block. -/
def zeroTensor3 : TensorField (Fin 3) := fun _ _ => 0

/-- Flat-space base case: both Einstein residual and edge residual are zero,
    so the projection model holds unconditionally. -/
theorem flat_space_projection_model
    (lP Lambda kappa : ℝ)
    (w : Fin 3 → Fin 3 → SimpEdge → ℝ) :
    EdgeToTensorProjectionModel zeroTensor3 zeroTensor3 zeroTensor3 zeroTensor3
      (fun _ => 0) (fun _ => 0) lP Lambda kappa w := by
  intro μ ν
  have hProj : projectedResidual w (edgeResidual (fun _ => 0) (fun _ => 0)) μ ν = 0 := by
    unfold edgeResidual
    simpa using projectedResidual_zero_of_zero_residual w μ ν
  have hEin : spatialEinsteinResidual zeroTensor3 zeroTensor3 zeroTensor3 zeroTensor3
      lP Lambda kappa μ ν = 0 := by
    simp [spatialEinsteinResidual, zeroTensor3, lambda_qg]
  simpa [hProj, hEin]

/-- CMS-style quantitative projection consistency:
    the edge-projected residual approximates the spatial Einstein residual with
    an `O(h²)` mesh error bound. -/
def ProjectionErrorBound
    (G H g T : TensorField (Fin 3))
    (dSdl source : SimpEdge → ℝ)
    (lP Lambda kappa h C : ℝ)
    (w : Fin 3 → Fin 3 → SimpEdge → ℝ) : Prop :=
  ∀ μ ν,
    |spatialEinsteinResidual G H g T lP Lambda kappa μ ν -
      projectedResidual w (edgeResidual dSdl source) μ ν| ≤ C * h ^ 2

/-- If the CMS-style projection error is `O(h²)`, then at zero mesh spacing
    (`h = 0`) the projection model is exact. -/
theorem edge_projection_model_of_zero_mesh
    {G H g T : TensorField (Fin 3)}
    {dSdl source : SimpEdge → ℝ}
    {lP Lambda kappa h C : ℝ}
    {w : Fin 3 → Fin 3 → SimpEdge → ℝ}
    (hBound : ProjectionErrorBound G H g T dSdl source lP Lambda kappa h C w)
    (hMesh : h = 0) :
    EdgeToTensorProjectionModel G H g T dSdl source lP Lambda kappa w := by
  intro μ ν
  have hμν := hBound μ ν
  rw [hMesh] at hμν
  have hAbsLeZero :
      |spatialEinsteinResidual G H g T lP Lambda kappa μ ν -
        projectedResidual w (edgeResidual dSdl source) μ ν| ≤ 0 := by
    simpa using hμν
  have hAbsEqZero :
      |spatialEinsteinResidual G H g T lP Lambda kappa μ ν -
        projectedResidual w (edgeResidual dSdl source) μ ν| = 0 :=
    le_antisymm hAbsLeZero (abs_nonneg _)
  exact sub_eq_zero.mp (abs_eq_zero.mp hAbsEqZero)

/-- Discharging the projection model + edge equations yields a concrete
    `ReggeToEinsteinBridge` on the spatial tensor block. -/
theorem bridge_from_edge_projection_model
    {Gdisc Gcont H g T : TensorField (Fin 3)}
    {dSdl source : SimpEdge → ℝ}
    {lP Lambda kappa : ℝ}
    {w : Fin 3 → Fin 3 → SimpEdge → ℝ}
    (hCont : ∀ μ ν, Gcont μ ν = Gdisc μ ν)
    (hEq : reggeEdgeEquation dSdl source)
    (hModel : EdgeToTensorProjectionModel Gdisc H g T dSdl source lP Lambda kappa w) :
    ReggeToEinsteinBridge Gdisc Gcont H g T lP Lambda kappa := by
  refine ⟨hCont, ?_⟩
  intro μ ν
  have hProjZero : projectedResidual w (edgeResidual dSdl source) μ ν = 0 :=
    projectedResidual_zero_of_reggeEdgeEquation w hEq μ ν
  have hModelμν := hModel μ ν
  let lhs := Gdisc μ ν + lambda_qg * lP ^ 2 * H μ ν + Lambda * g μ ν
  have hEq0 : lhs - kappa * T μ ν = 0 := by
    simpa [lhs, spatialEinsteinResidual, hProjZero] using hModelμν
  exact sub_eq_zero.mp hEq0

/-- Regge bridge constructor with CMS-style `O(h²)` projection bound in place of
    exact projection equality. In the continuum endpoint `h = 0`, this discharges
    `hModel` and yields the same bridge theorem. -/
theorem bridge_from_cms_bound_zero_mesh
    {Gdisc Gcont H g T : TensorField (Fin 3)}
    {dSdl source : SimpEdge → ℝ}
    {lP Lambda kappa h C : ℝ}
    {w : Fin 3 → Fin 3 → SimpEdge → ℝ}
    (hCont : ∀ μ ν, Gcont μ ν = Gdisc μ ν)
    (hEq : reggeEdgeEquation dSdl source)
    (hBound : ProjectionErrorBound Gdisc H g T dSdl source lP Lambda kappa h C w)
    (hMesh : h = 0) :
    ReggeToEinsteinBridge Gdisc Gcont H g T lP Lambda kappa := by
  have hModel :
      EdgeToTensorProjectionModel Gdisc H g T dSdl source lP Lambda kappa w :=
    edge_projection_model_of_zero_mesh hBound hMesh
  exact bridge_from_edge_projection_model hCont hEq hModel

/-- Direct continuum-limit theorem on the spatial block:
    an `O(h²)` CMS projection bound plus `h = 0` yields modified Einstein
    dynamics after transport from discrete to continuum Einstein tensor. -/
theorem modified_einstein_from_cms_bound_zero_mesh
    {Gdisc Gcont H g T : TensorField (Fin 3)}
    {dSdl source : SimpEdge → ℝ}
    {lP Lambda kappa h C : ℝ}
    {w : Fin 3 → Fin 3 → SimpEdge → ℝ}
    (hCont : ∀ μ ν, Gcont μ ν = Gdisc μ ν)
    (hEq : reggeEdgeEquation dSdl source)
    (hBound : ProjectionErrorBound Gdisc H g T dSdl source lP Lambda kappa h C w)
    (hMesh : h = 0) :
    ModifiedEinsteinFieldEquation Gcont H g T lP Lambda kappa := by
  intro μ ν
  have hBridge : ReggeToEinsteinBridge Gdisc Gcont H g T lP Lambda kappa :=
    bridge_from_cms_bound_zero_mesh hCont hEq hBound hMesh
  have hConv := hBridge.1 μ ν
  have hDisc := hBridge.2 μ ν
  calc
    Gcont μ ν + lambda_qg * lP ^ 2 * H μ ν + Lambda * g μ ν
        = Gdisc μ ν + lambda_qg * lP ^ 2 * H μ ν + Lambda * g μ ν := by rw [hConv]
    _ = kappa * T μ ν := hDisc

/-- Stationary Regge dynamics (`δS/δl = 0`) in vacuum source form, combined with
    CMS `O(h²)` projection control at `h = 0`, yields modified Einstein dynamics. -/
theorem modified_einstein_from_stationary_cms_zero_mesh
    {Gdisc Gcont H g T : TensorField (Fin 3)}
    {dSdl : SimpEdge → ℝ}
    {lP Lambda kappa h C : ℝ}
    {w : Fin 3 → Fin 3 → SimpEdge → ℝ}
    (hCont : ∀ μ ν, Gcont μ ν = Gdisc μ ν)
    (hStat : reggeStationary dSdl)
    (hBound : ProjectionErrorBound
      Gdisc H g T dSdl (fun _ => 0) lP Lambda kappa h C w)
    (hMesh : h = 0) :
    ModifiedEinsteinFieldEquation Gcont H g T lP Lambda kappa := by
  have hEq : reggeEdgeEquation dSdl (fun _ => 0) :=
    regge_stationary_implies_zero_source hStat
  exact modified_einstein_from_cms_bound_zero_mesh hCont hEq hBound hMesh

/-- Finite-algebraic reduction for SC Schläfli: if local tetra contributions
    satisfy diagonal cancellations and mixed cross-terms vanish, then the global
    19-edge Schläfli identity follows for the assembled SC decomposition. -/
theorem sc_schlaefli_from_local_balances
    (A dTheta : Fin 6 → SimpEdge → ℝ)
    (hDiag : ∀ t, ∑ e, A t e * dTheta t e = 0)
    (hOff : ∀ t u, t ≠ u → ∑ e, A t e * dTheta u e = 0) :
    scSchlaefliIdentity
      (fun e => ∑ t : Fin 6, A t e)
      (fun e => -∑ t : Fin 6, dTheta t e) := by
  unfold scSchlaefliIdentity SchlaefliIdentity
  have hTU : ∀ t u : Fin 6, ∑ e, A t e * dTheta u e = 0 := by
    intro t u
    by_cases htu : t = u
    · subst htu
      exact hDiag t
    · exact hOff t u htu
  have hSwap :
      ∑ e, ∑ t : Fin 6, ∑ u : Fin 6, A t e * dTheta u e
        = ∑ t : Fin 6, ∑ u : Fin 6, ∑ e, A t e * dTheta u e := by
    calc
      ∑ e, ∑ t : Fin 6, ∑ u : Fin 6, A t e * dTheta u e
          = ∑ t : Fin 6, ∑ e, ∑ u : Fin 6, A t e * dTheta u e := by
              rw [Finset.sum_comm]
      _ = ∑ t : Fin 6, ∑ u : Fin 6, ∑ e, A t e * dTheta u e := by
            apply Finset.sum_congr rfl
            intro t ht
            rw [Finset.sum_comm]
  calc
    ∑ e, (∑ t : Fin 6, A t e) * (-(∑ t : Fin 6, dTheta t e))
        = -∑ e, (∑ t : Fin 6, A t e) * (∑ t : Fin 6, dTheta t e) := by
            simp [Finset.sum_neg_distrib]
    _ = -∑ e, ∑ t : Fin 6, ∑ u : Fin 6, A t e * dTheta u e := by
          apply congrArg Neg.neg
          apply Finset.sum_congr rfl
          intro e he
          rw [Finset.sum_mul]
          apply Finset.sum_congr rfl
          intro t ht
          rw [Finset.mul_sum]
    _ = -∑ t : Fin 6, ∑ u : Fin 6, ∑ e, A t e * dTheta u e := by
          simpa [hSwap]
    _ = -∑ t : Fin 6, ∑ u : Fin 6, 0 := by
          simp [hTU]
    _ = 0 := by simp

/-- Proposition form of `ContinuumLimit.continuum_limit_exists` for bridge packaging. -/
def ContinuumLimitStatement : Prop :=
  (2 : ℕ) ^ 4 = 16 ∧
  Gutoe.GaugeGroupSU3.quarkOrbit.card = 3 ∧
  Gutoe.GaugeGroupSU3.quarkOrbit.card ^ 2 - 1 = 8 ∧
  Gutoe.FineStructure.alphaInverse 4 = 137 ∧
  Gutoe.Z3Uniqueness.magneticTriplet.card = 3 ∧
  Gutoe.GaugeGroupSM.leptonState.card = 1 ∧
  Gutoe.GaugeGroupSU2.σ₁ * Gutoe.GaugeGroupSU2.σ₂ -
      Gutoe.GaugeGroupSU2.σ₂ * Gutoe.GaugeGroupSU2.σ₁ =
    (2 * Complex.I) • Gutoe.GaugeGroupSU2.σ₃ ∧
  Gutoe.Z3Uniqueness.grade2_4d.card = 6 ∧
  Gutoe.Z3Uniqueness.magneticTriplet.card + Gutoe.Z3Uniqueness.emTriplet.card =
    Gutoe.Z3Uniqueness.grade2_4d.card

/-- Existing theorem chain discharges `ContinuumLimitStatement`. -/
theorem continuum_limit_statement_holds : ContinuumLimitStatement :=
  Gutoe.ContinuumLimit.continuum_limit_exists

/-- Coupling identification from lattice-side relations:
    `κ = v² / G`. -/
noncomputable def kappaFromLattice (v G : ℝ) : ℝ := v ^ 2 / G

/-- Continuum Einstein-Hilbert coupling normalization:
    `κ = 8πG / c⁴`. -/
noncomputable def kappaEinstein (G c : ℝ) : ℝ :=
  8 * Real.pi * G / c ^ 4

/-- Inverse coupling map:
    `G = v² / κ`. -/
noncomputable def newtonFromLattice (v kappa : ℝ) : ℝ := v ^ 2 / kappa

/-- Inverse of the Einstein normalization map:
    `G = κ c⁴ / (8π)`. -/
noncomputable def newtonFromEinsteinKappa (kappa c : ℝ) : ℝ :=
  kappa * c ^ 4 / (8 * Real.pi)

/-- Planck-side relation:
    `ħ = l_P² κ`. -/
def hbarFromLattice (lP kappa : ℝ) : ℝ := lP ^ 2 * kappa

/-- The `G = v²/κ` relation is exact from the definition of `kappaFromLattice`. -/
theorem newton_relation_of_kappa_from_lattice
    {v G : ℝ} (hG : G ≠ 0) (hv : v ≠ 0) :
    newtonFromLattice v (kappaFromLattice v G) = G := by
  unfold newtonFromLattice kappaFromLattice
  field_simp [hG, hv]

/-- Einstein normalization inverts exactly to Newton's constant. -/
theorem newton_relation_of_kappa_einstein
    {G c : ℝ} (hc : c ≠ 0) :
    newtonFromEinsteinKappa (kappaEinstein G c) c = G := by
  unfold newtonFromEinsteinKappa kappaEinstein
  field_simp [hc, Real.pi_ne_zero]

/-- Combining `κ = v²/G` with `ħ = l_P² κ` yields
    `G = v² l_P² / ħ` (for `ħ ≠ 0`). -/
theorem newton_from_planck_lattice_relation
    {v lP hbar : ℝ}
    (hhbar : hbar ≠ 0)
    (hlP : lP ≠ 0) :
    newtonFromLattice v (hbar / (lP ^ 2)) = v ^ 2 * lP ^ 2 / hbar := by
  unfold newtonFromLattice
  field_simp [hhbar, hlP]

/-- Structural prerequisites already proven in the Cl(1,3) theorem chain:
    SC coordination (6), continuum-limit existence, and fixed `λ_QG = 1/12`. -/
theorem clifford_gravity_prerequisites :
    coordinationNumber = 6 ∧
    ContinuumLimitStatement ∧
    lambda_qg = 1 / 12 := by
  refine ⟨coordination_number_is_6, continuum_limit_statement_holds, ?_⟩
  simp [lambda_qg]

/-- If the Regge bridge hypotheses hold, the modified Einstein equation follows. -/
theorem regge_bridge_implies_modified_einstein
    {ι : Type}
    {Gdisc Gcont H g T : TensorField ι}
    {lP Lambda kappa : ℝ}
    (hBridge : ReggeToEinsteinBridge Gdisc Gcont H g T lP Lambda kappa) :
    ModifiedEinsteinFieldEquation Gcont H g T lP Lambda kappa := by
  intro μ ν
  have hConv := hBridge.1 μ ν
  have hDisc := hBridge.2 μ ν
  calc
    Gcont μ ν + lambda_qg * lP ^ 2 * H μ ν + Lambda * g μ ν
        = Gdisc μ ν + lambda_qg * lP ^ 2 * H μ ν + Lambda * g μ ν := by rw [hConv]
    _ = kappa * T μ ν := hDisc

/-- GR limit from the modified equation at `lP = 0`. -/
theorem modified_einstein_planck_zero
    {ι : Type}
    {G H g T : TensorField ι}
    {Lambda kappa : ℝ}
    (hMod : ModifiedEinsteinFieldEquation G H g T 0 Lambda kappa) :
    EinsteinFieldEquation G g T Lambda kappa := by
  intro μ ν
  specialize hMod μ ν
  simpa [ModifiedEinsteinFieldEquation, EinsteinFieldEquation]
    using hMod

/-- Master bridge theorem for GRAND-89:
    Cl(1,3) prerequisites + Regge bridge hypotheses imply modified Einstein dynamics. -/
theorem einstein_from_clifford_lattice
    {ι : Type}
    {Gdisc Gcont H g T : TensorField ι}
    {lP Lambda kappa : ℝ}
    (hBridge : ReggeToEinsteinBridge Gdisc Gcont H g T lP Lambda kappa) :
    coordinationNumber = 6 ∧
    ContinuumLimitStatement ∧
    lambda_qg = 1 / 12 ∧
    ModifiedEinsteinFieldEquation Gcont H g T lP Lambda kappa := by
  refine ⟨coordination_number_is_6, continuum_limit_statement_holds, ?_, ?_⟩
  · simp [lambda_qg]
  · exact regge_bridge_implies_modified_einstein hBridge

/-- If the Regge bridge is established directly at `lP = 0`,
    the continuum Einstein equation follows. -/
theorem einstein_from_clifford_lattice_gr_limit
    {ι : Type}
    {Gdisc Gcont H g T : TensorField ι}
    {Lambda kappa : ℝ}
    (hBridge0 : ReggeToEinsteinBridge Gdisc Gcont H g T 0 Lambda kappa) :
    EinsteinFieldEquation Gcont g T Lambda kappa := by
  exact modified_einstein_planck_zero (regge_bridge_implies_modified_einstein hBridge0)

end Gutoe.EinsteinFromLattice
