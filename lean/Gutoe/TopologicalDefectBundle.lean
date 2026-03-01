import Mathlib

/-!
GUTOE — Topological defect bundle lane scaffold

A minimal geometric bridge model:
- Base distance on a 1D projection is Euclidean (`|b-a|`).
- A defect can introduce an identification bridge between two points `l` and `r`.
- Path distance becomes the minimum of direct and bridge-assisted routes.

This captures "bundle structure changed by defect" at the kinematic level.
-/

namespace Gutoe.TopologicalDefectBundle

/-- Base 1D projected distance. -/
def baseDistance (a b : ℝ) : ℝ := |b - a|

/-- Defect-assisted distance with a bridge between `l` and `r`. -/
def defectDistance (a b l r : ℝ) : ℝ :=
  min (baseDistance a b)
    (min (|a - l| + |b - r|) (|a - r| + |b - l|))

/-- Defect distance is never worse than base distance. -/
theorem defect_distance_le_base (a b l r : ℝ) :
    defectDistance a b l r ≤ baseDistance a b := by
  unfold defectDistance
  exact min_le_left _ _

/-- If one bridge-assisted route is strictly shorter than base, defect distance
is strictly shorter than base. -/
theorem defect_strict_improvement_if_bridge_shorter
    (a b l r : ℝ)
    (hshort : |a - l| + |b - r| < baseDistance a b) :
    defectDistance a b l r < baseDistance a b := by
  unfold defectDistance
  have hmid : min (|a - l| + |b - r|) (|a - r| + |b - l|) < baseDistance a b := by
    exact lt_of_le_of_lt (min_le_left _ _) hshort
  exact min_lt_iff.2 (Or.inr hmid)

/-- Canonical witness: endpoints exactly on the defect bridge collapse to zero
defect distance while base distance stays positive when `l ≠ r`. -/
theorem bridge_endpoints_strict_shortcut (l r : ℝ) (hneq : l ≠ r) :
    defectDistance l r l r < baseDistance l r := by
  have hbase : 0 < baseDistance l r := by
    unfold baseDistance
    have h : r - l ≠ 0 := sub_ne_zero.mpr (Ne.symm hneq)
    exact abs_pos.mpr h
  have hshort : |l - l| + |r - r| < baseDistance l r := by
    simp [hbase]
  exact defect_strict_improvement_if_bridge_shorter l r l r hshort

end Gutoe.TopologicalDefectBundle
