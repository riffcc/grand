import Mathlib
import Gutoe.OverdeterminedTopologyRatios

/-!
GUTOE — CTC Retro-Retry Closure

This lane formalizes the retry closure claim:

- single-attempt success probability `p`,
- independent retry depth `n`,
- eventual success after `n` attempts:
  `S_n = 1 - (1 - p)^n`.

Key result:
if `0 < p < 1`, `S_n` is strictly increasing and tends to `1` as retries grow.
This captures "nonzero success + retry channel => eventual closure."
-/

namespace Gutoe.CTCRetroRetryClosure

open Gutoe.OverdeterminedTopologyRatios

/-- Eventual success probability after `n` independent retry attempts. -/
def retrySuccessProb (p : ℝ) (n : ℕ) : ℝ := 1 - (1 - p) ^ n

/-- Failure probability after `n` independent retry attempts. -/
def retryFailProb (p : ℝ) (n : ℕ) : ℝ := (1 - p) ^ n

theorem retry_partition (p : ℝ) (n : ℕ) :
    retrySuccessProb p n + retryFailProb p n = 1 := by
  unfold retrySuccessProb retryFailProb
  ring

theorem retry_success_zero (p : ℝ) :
    retrySuccessProb p 0 = 0 := by
  unfold retrySuccessProb
  norm_num

theorem retry_success_one (p : ℝ) :
    retrySuccessProb p 1 = p := by
  unfold retrySuccessProb
  ring

/-- Monotone retry closure for `0 ≤ p ≤ 1`. -/
theorem retry_success_monotone
    (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    Monotone (retrySuccessProb p) := by
  intro m n hmn
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hmn
  have hbase_nonneg : 0 ≤ 1 - p := by linarith
  have hbase_le_one : 1 - p ≤ 1 := by linarith
  have hpow_m_nonneg : 0 ≤ (1 - p) ^ m := by
    exact pow_nonneg hbase_nonneg m
  have hpow_k_le_one : (1 - p) ^ k ≤ 1 := by
    exact pow_le_one₀ hbase_nonneg hbase_le_one
  have hpow : (1 - p) ^ (m + k) ≤ (1 - p) ^ m := by
    calc
      (1 - p) ^ (m + k) = (1 - p) ^ m * (1 - p) ^ k := by rw [pow_add]
      _ ≤ (1 - p) ^ m * 1 := by
            exact mul_le_mul_of_nonneg_left hpow_k_le_one hpow_m_nonneg
      _ = (1 - p) ^ m := by ring
  unfold retrySuccessProb
  linarith

/-- Strict improvement per retry if `0 < p < 1`. -/
theorem retry_success_strict_step
    (p : ℝ) (n : ℕ)
    (hp0 : 0 < p) (hp1 : p < 1) :
    retrySuccessProb p n < retrySuccessProb p (n + 1) := by
  have hbase_pos : 0 < 1 - p := by linarith
  unfold retrySuccessProb
  have hpow_pos : 0 < (1 - p) ^ n := by exact pow_pos hbase_pos n
  have hmul_pos : 0 < p * (1 - p) ^ n := mul_pos hp0 hpow_pos
  have hstep :
      1 - (1 - p) ^ (n + 1) = (1 - (1 - p) ^ n) + p * (1 - p) ^ n := by
    rw [pow_succ]
    ring
  rw [hstep]
  linarith

/-- Retry closure tends to certainty for any non-degenerate success rate. -/
theorem retry_success_tendsto_one
    (p : ℝ) (hp0 : 0 < p) (hp1 : p < 1) :
    Filter.Tendsto (fun n : ℕ => retrySuccessProb p n) Filter.atTop (nhds 1) := by
  have hbase_nonneg : 0 ≤ 1 - p := by linarith
  have habs : |1 - p| < 1 := by
    rw [abs_of_nonneg hbase_nonneg]
    linarith
  have hpow0 : Filter.Tendsto (fun n : ℕ => (1 - p) ^ n) Filter.atTop (nhds 0) := by
    exact tendsto_pow_atTop_nhds_zero_of_abs_lt_one habs
  have hsub : Filter.Tendsto (fun n : ℕ => 1 - (1 - p) ^ n) Filter.atTop (nhds (1 - 0)) := by
    exact tendsto_const_nhds.sub hpow0
  simpa [retrySuccessProb] using hsub

/-- For any target below `1`, enough retries exceed it. -/
  theorem exists_retry_count_for_target
    (p target : ℝ)
    (hp0 : 0 < p) (hp1 : p < 1)
    (ht : target < 1) :
    ∃ N : ℕ, target < retrySuccessProb p N := by
  have hlim : Filter.Tendsto (fun n : ℕ => retrySuccessProb p n) Filter.atTop (nhds 1) :=
    retry_success_tendsto_one p hp0 hp1
  have hset : Set.Ioi target ∈ nhds (1 : ℝ) := Ioi_mem_nhds ht
  have hev : ∀ᶠ n in Filter.atTop, retrySuccessProb p n ∈ Set.Ioi target := hlim hset
  rcases (Filter.eventually_atTop.1 hev) with ⟨N, hN⟩
  refine ⟨N, ?_⟩
  exact hN N (le_rfl)

/-- Structural seed success from Cl(1,3) void split (`3/16` as a canonical
nonzero retry seed). -/
def structuralRetrySeed : ℝ := (voidStructuralQ : ℝ)

theorem structural_retry_seed_bounds :
    0 < structuralRetrySeed ∧ structuralRetrySeed < 1 := by
  unfold structuralRetrySeed
  rw [void_structural_eq_3_16]
  norm_num

/-- Structural closure corollary: for any target `< 1`, the canonical structural
seed eventually exceeds it under retries. -/
theorem structural_seed_eventual_closure
    (target : ℝ) (ht : target < 1) :
    ∃ N : ℕ, target < retrySuccessProb structuralRetrySeed N := by
  rcases structural_retry_seed_bounds with ⟨hp0, hp1⟩
  exact exists_retry_count_for_target structuralRetrySeed target hp0 hp1 ht

end Gutoe.CTCRetroRetryClosure
