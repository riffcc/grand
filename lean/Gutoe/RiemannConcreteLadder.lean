import Mathlib
import Gutoe.RiemannTargetFiniteLadder
import Gutoe.RiemannFiniteXiModel
import Gutoe.RiemannCore

namespace Gutoe.RiemannConcreteLadder

open Gutoe.RiemannTargetFiniteLadder
open Gutoe.RiemannFiniteXiModel
open Gutoe.RiemannCore

noncomputable section

/-- Public 120-ordinate reference lane (from the shared physics dataset). -/
def referenceOrdinates : List ℝ :=
[
  14.134725142,
  21.022039639,
  25.010857580,
  30.424876126,
  32.935061588,
  37.586178159,
  40.918719012,
  43.327073281,
  48.005150881,
  49.773832478,
  52.970321478,
  56.446247697,
  59.347044003,
  60.831778525,
  65.112544048,
  67.079810529,
  69.546401711,
  72.067157674,
  75.704690699,
  77.144840069,
  79.337375020,
  82.910380854,
  84.735492981,
  87.425274613,
  88.809111208,
  92.491899271,
  94.651344041,
  95.870634228,
  98.831194218,
  101.317851006,
  103.725538040,
  105.446623052,
  107.168611184,
  111.029535543,
  111.874659177,
  114.320220915,
  116.226680321,
  118.790782866,
  121.370125002,
  122.946829294,
  124.256818554,
  127.516683880,
  129.578704200,
  131.087688531,
  133.497737203,
  134.756509753,
  138.116042055,
  139.736208952,
  141.123707404,
  143.111845808,
  146.000982487,
  147.422765343,
  150.053520421,
  150.925257612,
  153.024693811,
  156.112909294,
  157.597591818,
  158.849988171,
  161.188964138,
  163.030709687,
  165.537069188,
  167.184439978,
  169.094515416,
  169.911976479,
  173.411536520,
  174.754191523,
  176.441434298,
  178.377407776,
  179.916484020,
  182.207078484,
  184.874467848,
  185.598783678,
  187.228922584,
  189.416158656,
  192.026656361,
  193.079726604,
  195.265396680,
  196.876481841,
  198.015309676,
  201.264751944,
  202.493594514,
  204.189671803,
  205.394697202,
  207.906258888,
  209.576509717,
  211.690862595,
  213.347919360,
  214.547044783,
  216.169538508,
  219.067596349,
  220.714918839,
  221.430705555,
  224.007000255,
  224.983324670,
  227.421444280,
  229.337413306,
  231.250188700,
  231.987235253,
  233.693404179,
  236.524229666,
  237.769820481,
  239.555477573,
  241.049157796,
  242.823271934,
  244.070898497,
  247.136990075,
  248.101990060,
  249.573689645,
  251.014947795,
  253.069986748,
  255.306256455,
  256.380713694,
  258.610439492,
  259.874406990,
  260.805084505,
  263.573893905,
  265.557851839,
  266.614973782,
  267.921915083,
  269.970449024
]

/-- Total ordinate function used by the concrete ladder.
For indices beyond the table length, returns `0`. -/
def referenceOrdinate (n : ℕ) : ℝ :=
  referenceOrdinates.getD n 0

/-- Concrete finite-prefix ladder built from the reference ordinate table. -/
def referenceSpecN : ℕ → Finset ℝ :=
  prefixSpec referenceOrdinate

theorem referenceSpecN_nested (N : ℕ) :
    referenceSpecN N ⊆ referenceSpecN (N + 1) := by
  intro t ht
  unfold referenceSpecN prefixSpec at *
  rcases Finset.mem_image.mp ht with ⟨k, hk, rfl⟩
  refine Finset.mem_image.mpr ?_
  refine ⟨k, ?_, rfl⟩
  exact Finset.mem_range.mpr (Nat.lt_trans (Finset.mem_range.mp hk) (Nat.lt_succ_self (N + 1)))

/-- Every listed ordinate index up to `N` is present in the `N`th concrete prefix. -/
theorem mem_referenceSpecN_of_le {n N : ℕ} (hn : n ≤ N) :
    referenceOrdinate n ∈ referenceSpecN N := by
  unfold referenceSpecN prefixSpec
  refine Finset.mem_image.mpr ?_
  exact ⟨n, Finset.mem_range.mpr (Nat.lt_succ_of_le hn), rfl⟩

/-- Finite-level exactness on listed ordinates:
`XiFinite (referenceSpecN N)` vanishes at each indexed reference ordinate `n ≤ N`. -/
theorem xiFinite_referenceSpecN_zero_of_le {n N : ℕ} (hn : n ≤ N) :
    XiFinite (referenceSpecN N) (criticalLinePoint (referenceOrdinate n)) = 0 := by
  exact XiFinite_zero_of_mem (referenceSpecN N) (mem_referenceSpecN_of_le hn)

/-- Pointwise finite-ladder capture for each listed reference ordinate. -/
theorem reference_capture_of_index (n : ℕ) :
    ∃ N : ℕ, ∃ t : ℝ, t ∈ referenceSpecN N ∧
      criticalLinePoint (referenceOrdinate n) = criticalLinePoint t := by
  refine ⟨n, referenceOrdinate n, mem_referenceSpecN_of_le (Nat.le_refl n), rfl⟩

/-- In the concrete ladder, zero-tolerance convergence is equivalent to exact
target zero-capture (canonical reduction theorem specialized). -/
theorem approxZero_tolZero_iff_reference_capture :
    Gutoe.RiemannConvergenceTransfer.ApproxZeroConvergence
      Gutoe.RiemannFinalTarget.XiTarget (XiFiniteLadder referenceSpecN)
      Gutoe.RiemannTargetFiniteLadder.tolZero
      ↔ XiTargetLadderZeroCapture referenceSpecN :=
  approxZero_tolZero_iff_zeroCapture referenceSpecN

/-- Concrete-ladder nontrivial capture from a direct ordinate-enumeration hypothesis
on the reference ordinate function. -/
theorem reference_nontrivial_capture_of_ordinate_enumeration
    (hEnum : RiemannNontrivialZeroOrdinateEnumeration referenceOrdinate) :
    RiemannNontrivialLadderZeroCapture referenceSpecN := by
  simpa [referenceSpecN] using
    (ladder_capture_of_ordinate_enumeration referenceOrdinate hEnum)

/-- Concrete-ladder RH closure surface:
if the reference ladder captures all nontrivial `ζ` zeros, RH follows. -/
theorem mathlibRH_of_reference_nontrivial_capture
    (hCap : RiemannNontrivialLadderZeroCapture referenceSpecN) :
    RiemannHypothesis :=
  mathlibRH_of_nontrivial_capture referenceSpecN hCap

/-- Concrete-ladder RH closure from a direct ordinate-enumeration hypothesis
on the reference ordinate function. -/
theorem mathlibRH_of_reference_ordinate_enumeration
    (hEnum : RiemannNontrivialZeroOrdinateEnumeration referenceOrdinate) :
    RiemannHypothesis :=
  mathlibRH_of_reference_nontrivial_capture
    (reference_nontrivial_capture_of_ordinate_enumeration hEnum)

end

end Gutoe.RiemannConcreteLadder
