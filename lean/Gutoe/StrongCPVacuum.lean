/- 
 * GUTOE — Strong CP Vacuum-Sector Formalization (GRAND-267)
 *
 * Route-1 skeleton formalization:
 *   If only topological sector Q=0 is accessible, then θ is not a physical
 *   parameter of the partition function (θ-independent observables in this model).
 *
 * This module encodes that implication abstractly, then instantiates it with
 * the current Cl(1,3) structural support model from StrongCP.
 -/

import Mathlib
import Gutoe.StrongCP

namespace Gutoe.StrongCPVacuum

open Gutoe.StrongCP
open scoped BigOperators

/-- Sector support set over integer topological charges. -/
abbrev SectorSupport := Finset ℤ

/-- Fundamental gauge carrier in this path: the discrete Z₃ color orbit. -/
abbrev FundamentalGaugeGroup : Type := Fin 3

/-- Fundamental gauge carrier has exactly three states. -/
theorem fundamental_gauge_group_card : Fintype.card FundamentalGaugeGroup = 3 := by
  decide

/-- The fundamental gauge carrier is discrete. -/
theorem fundamental_gauge_group_discrete : DiscreteTopology FundamentalGaugeGroup := by
  infer_instance

/-- Any continuous map from a preconnected domain into the fundamental discrete
    gauge carrier is constant. This is the route-1 topology gate. -/
theorem continuous_to_fundamental_group_constant
    {X : Type} [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    (f : C(X, FundamentalGaugeGroup)) :
    ∃ g0 : FundamentalGaugeGroup, ∀ x : X, f x = g0 := by
  classical
  have himage : IsPreconnected (Set.range f) := by
    simpa [Set.range_comp] using (isPreconnected_univ.image f f.continuous.continuousOn)
  have hs : (Set.range f).Subsingleton := himage.subsingleton
  let x0 : X := Classical.choice (inferInstance : Nonempty X)
  refine ⟨f x0, ?_⟩
  intro x
  exact hs ⟨x, rfl⟩ ⟨x0, rfl⟩

/-- With a fixed basepoint value `0`, the continuous map is forced to the
    trivial constant map. -/
theorem based_continuous_to_fundamental_group_zero
    {X : Type} [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    (x0 : X) (f : C(X, FundamentalGaugeGroup)) (hbase : f x0 = 0) :
    ∀ x : X, f x = 0 := by
  rcases continuous_to_fundamental_group_constant f with ⟨g0, hg0⟩
  have hg0_zero : g0 = 0 := by
    calc
      g0 = f x0 := (hg0 x0).symm
      _ = 0 := hbase
  intro x
  calc
    f x = g0 := hg0 x
    _ = 0 := hg0_zero

/-- Based continuous maps into the discrete fundamental gauge carrier are
    unique (no nontrivial based topological sectors). -/
theorem based_fundamental_sector_unique
    {X : Type} [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    (x0 : X) :
    Subsingleton {f : C(X, FundamentalGaugeGroup) // f x0 = 0} := by
  classical
  refine ⟨?_⟩
  intro a b
  cases' a with fa hfa
  cases' b with fb hfb
  have hfa0 : ∀ x : X, fa x = 0 := based_continuous_to_fundamental_group_zero x0 fa hfa
  have hfb0 : ∀ x : X, fb x = 0 := based_continuous_to_fundamental_group_zero x0 fb hfb
  have hfun : fa = fb := by
    ext x
    rw [hfa0 x, hfb0 x]
  simp [hfun]

/-- Route-1 hypothesis: only `Q = 0` is accessible. -/
def onlyZeroSupport (S : SectorSupport) : Prop :=
  ∀ n : ℤ, n ∈ S → n = 0

/-- Real channel of a finite-sector partition sum. -/
noncomputable def zRe (S : SectorSupport) (w : ℤ → ℝ) (theta : ℝ) : ℝ :=
  Finset.sum S (fun n => w n * Real.cos (theta * (n : ℝ)))

/-- Imaginary channel of a finite-sector partition sum. -/
noncomputable def zIm (S : SectorSupport) (w : ℤ → ℝ) (theta : ℝ) : ℝ :=
  Finset.sum S (fun n => w n * Real.sin (theta * (n : ℝ)))

/-- If support is only `Q=0`, the CP-odd channel vanishes for all `θ`. -/
theorem zIm_zero_of_onlyZeroSupport
    (S : SectorSupport) (w : ℤ → ℝ) (hS : onlyZeroSupport S) :
    ∀ theta : ℝ, zIm S w theta = 0 := by
  intro theta
  unfold zIm
  refine Finset.sum_eq_zero ?_
  intro n hn
  have hn0 : n = 0 := hS n hn
  rw [hn0]
  norm_num

/-- If support is only `Q=0`, the real channel is θ-independent. -/
theorem zRe_theta_invariant_of_onlyZeroSupport
    (S : SectorSupport) (w : ℤ → ℝ) (hS : onlyZeroSupport S) :
    ∀ theta1 theta2 : ℝ, zRe S w theta1 = zRe S w theta2 := by
  intro theta1 theta2
  unfold zRe
  refine Finset.sum_congr rfl ?_
  intro n hn
  have hn0 : n = 0 := hS n hn
  rw [hn0]
  norm_num

/-- Route-1 consequence:
    with only `Q=0` support, both channels are θ-invariant (θ unphysical). -/
theorem theta_unphysical_of_onlyZeroSupport
    (S : SectorSupport) (w : ℤ → ℝ) (hS : onlyZeroSupport S) :
    ∀ theta1 theta2 : ℝ,
      zRe S w theta1 = zRe S w theta2 ∧ zIm S w theta1 = zIm S w theta2 := by
  intro theta1 theta2
  constructor
  · exact zRe_theta_invariant_of_onlyZeroSupport S w hS theta1 theta2
  · rw [zIm_zero_of_onlyZeroSupport S w hS theta1, zIm_zero_of_onlyZeroSupport S w hS theta2]

/-- If an effective support set inherits from the fundamental `{0}` sector
    support, then it cannot populate nonzero sectors. -/
theorem onlyZeroSupport_of_inherits_zero
    (S : SectorSupport)
    (hInherit : S ⊆ ({0} : SectorSupport)) :
    onlyZeroSupport S := by
  intro n hn
  have h0 : n ∈ ({0} : SectorSupport) := hInherit hn
  simpa using Finset.mem_singleton.mp h0

/-- Current Cl(1,3) structural sector support model: singleton at structural source. -/
def cl13Support : SectorSupport := {cpOddSectorImbalance}

/-- Cl(1,3) support model satisfies Route-1 hypothesis (`Q=0` only). -/
theorem cl13_support_only_zero : onlyZeroSupport cl13Support := by
  intro n hn
  have hnEq : n = cpOddSectorImbalance := by
    simpa [cl13Support] using Finset.mem_singleton.mp hn
  rw [hnEq, cp_odd_sector_imbalance_zero]

/-- In the current Cl(1,3) support model, `zIm` vanishes for all θ. -/
theorem cl13_zIm_zero (w : ℤ → ℝ) : ∀ theta : ℝ, zIm cl13Support w theta = 0 :=
  zIm_zero_of_onlyZeroSupport cl13Support w cl13_support_only_zero

/-- In the current Cl(1,3) support model, θ is unphysical in the finite-sector sum. -/
theorem cl13_theta_unphysical (w : ℤ → ℝ) :
    ∀ theta1 theta2 : ℝ,
      zRe cl13Support w theta1 = zRe cl13Support w theta2 ∧
      zIm cl13Support w theta1 = zIm cl13Support w theta2 :=
  theta_unphysical_of_onlyZeroSupport cl13Support w cl13_support_only_zero

/-- Scale-lift of the same support model (explicitly no scale dependence). -/
def cl13SupportAtScale (_a : ℝ) : SectorSupport := cl13Support

/-- The Route-1 support claim is scale-stable in this model. -/
theorem cl13_route1_scale_stable :
    ∀ a : ℝ, onlyZeroSupport (cl13SupportAtScale a) := by
  intro a
  simpa [cl13SupportAtScale] using cl13_support_only_zero

end Gutoe.StrongCPVacuum
