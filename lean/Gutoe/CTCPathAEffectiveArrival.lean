import Mathlib

/-!
GUTOE — Path-A effective arrival algebra

For a path that:
- takes `dtIn` to access the identification channel,
- spends `n` loops of period `T`,
- takes `dtOut` to reach destination,
- applies `q` coordinate-time shift quanta per loop,

the effective arrival coordinate is
`dtIn + dtOut + n*(1-q)*T`.
-/

namespace Gutoe.CTCPathAEffectiveArrival

/-- Effective coordinate-time arrival for Path-A traversal model. -/
def effectiveArrival (dtIn dtOut T : ℝ) (n q : ℕ) : ℝ :=
  dtIn + dtOut + (n : ℝ) * (1 - (q : ℝ)) * T

/-- Cover-time (unwrapped) arrival in the same model. -/
def coverArrival (dtIn dtOut T : ℝ) (n : ℕ) : ℝ :=
  dtIn + dtOut + (n : ℝ) * T

/-- Effective/cover relation via loop shift quanta. -/
theorem effective_equals_cover_minus_shift
    (dtIn dtOut T : ℝ) (n q : ℕ) :
    effectiveArrival dtIn dtOut T n q =
      coverArrival dtIn dtOut T n - (n : ℝ) * (q : ℝ) * T := by
  unfold effectiveArrival coverArrival
  ring

/-- If `q = 1`, loops do not change effective coordinate arrival. -/
theorem q_one_no_coordinate_gain
    (dtIn dtOut T : ℝ) (n : ℕ) :
    effectiveArrival dtIn dtOut T n 1 = dtIn + dtOut := by
  unfold effectiveArrival
  ring

/-- For positive `T` and `q > 1`, sufficiently many loops force pre-departure
(`effectiveArrival < 0`) provided baseline access+egress is finite. -/
theorem q_gt_one_predeparture_possible
    (dtIn dtOut T : ℝ) (q : ℕ)
    (hT : 0 < T) (hq : 1 < q) :
    ∃ n : ℕ, effectiveArrival dtIn dtOut T n q < 0 := by
  let k : ℝ := ((q : ℝ) - 1) * T
  have hk : 0 < k := by
    have hqR : (1 : ℝ) < (q : ℝ) := by exact_mod_cast hq
    have hqPos : 0 < (q : ℝ) - 1 := by linarith
    exact mul_pos hqPos hT
  rcases exists_nat_gt ((dtIn + dtOut) / k) with ⟨n, hn⟩
  refine ⟨n, ?_⟩
  have hnk : dtIn + dtOut < (n : ℝ) * k := by
    exact (div_lt_iff₀ hk).1 hn
  have hqk : (1 - (q : ℝ)) * T = -k := by
    dsimp [k]
    ring
  calc
    effectiveArrival dtIn dtOut T n q
        = dtIn + dtOut + (n : ℝ) * ((1 - (q : ℝ)) * T) := by
            unfold effectiveArrival
            ring
    _ = dtIn + dtOut - (n : ℝ) * k := by rw [hqk]; ring
    _ < 0 := by linarith [hnk]

end Gutoe.CTCPathAEffectiveArrival
