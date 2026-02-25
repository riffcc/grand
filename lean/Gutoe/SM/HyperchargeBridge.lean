/- 
 * GUTOE — Hypercharge Bridge from Anomaly Constraints
 *
 * Algebraic bridge layer: anomaly equations constrain one-generation
 * hypercharges to the Standard-Model pattern up to normalization/sign choice.
-/

import Mathlib
import Gutoe.SM.Rep

namespace Gutoe.SM.HyperchargeBridge

open Gutoe.SM.Rep

/--
Anomaly constraints force quark-singlet hypercharges to two algebraic branches.
All fields are left-chiral Weyl variables:
- `q`  : quark doublet hypercharge
- `u`  : charge-conjugated up singlet hypercharge
- `d`  : charge-conjugated down singlet hypercharge
- `ℓ`  : lepton doublet hypercharge
- `e`  : charge-conjugated charged-lepton singlet hypercharge

Constraints:
- SU(3)^2·U(1): `u + d = -2 q`
- SU(2)^2·U(1): `ℓ = -3 q`
- grav^2·U(1): `6q + 3u + 3d + 2ℓ + e = 0`
- U(1)^3: `6q^3 + 3u^3 + 3d^3 + 2ℓ^3 + e^3 = 0`
-/
theorem hypercharge_branches_from_anomaly_constraints
    {q u d ℓ e : ℚ}
    (hsu3 : u + d = -2 * q)
    (hsu2 : ℓ = -3 * q)
    (hgrav : 6 * q + 3 * u + 3 * d + 2 * ℓ + e = 0)
    (hu1cube : 6 * q ^ 3 + 3 * u ^ 3 + 3 * d ^ 3 + 2 * ℓ ^ 3 + e ^ 3 = 0)
    (hq : q ≠ 0) :
    (u = 2 * q ∧ d = -4 * q ∧ e = 6 * q ∧ ℓ = -3 * q) ∨
    (u = -4 * q ∧ d = 2 * q ∧ e = 6 * q ∧ ℓ = -3 * q) := by
  have hℓ : ℓ = -3 * q := hsu2
  have hd : d = -2 * q - u := by linarith [hsu3]
  have he : e = 6 * q := by linarith [hgrav, hsu3, hsu2]

  have hpoly : (-18 : ℚ) * q * (u ^ 2 + 2 * q * u - 8 * q ^ 2) = 0 := by
    rw [hd, hℓ, he] at hu1cube
    ring_nf at hu1cube ⊢
    exact hu1cube

  have hq6 : (-18 : ℚ) * q ≠ 0 := by
    exact mul_ne_zero (by norm_num) hq

  have hquad : u ^ 2 + 2 * q * u - 8 * q ^ 2 = 0 := by
    exact (mul_eq_zero.mp hpoly).resolve_left hq6

  have hfact : (u - 2 * q) * (u + 4 * q) = 0 := by
    have hmul : (u - 2 * q) * (u + 4 * q) = u ^ 2 + 2 * q * u - 8 * q ^ 2 := by ring
    rw [hmul, hquad]

  rcases mul_eq_zero.mp hfact with hu2 | hu4
  · left
    have hu2' : u = 2 * q := by linarith
    refine ⟨hu2', ?_, he, hℓ⟩
    · linarith [hd, hu2']
  · right
    have hu4' : u = -4 * q := by linarith
    refine ⟨hu4', ?_, he, hℓ⟩
    linarith [hd, hu4']

/--
Choose the physical branch (`u < d`) and positive normalization (`q > 0`):
this fixes the Standard-Model sign branch.
-/
theorem physical_branch_selects_sm_sign
    {q u d ℓ e : ℚ}
    (hsu3 : u + d = -2 * q)
    (hsu2 : ℓ = -3 * q)
    (hgrav : 6 * q + 3 * u + 3 * d + 2 * ℓ + e = 0)
    (hu1cube : 6 * q ^ 3 + 3 * u ^ 3 + 3 * d ^ 3 + 2 * ℓ ^ 3 + e ^ 3 = 0)
    (hq : q ≠ 0)
    (hqpos : 0 < q)
    (hud : u < d) :
    u = -4 * q ∧ d = 2 * q ∧ e = 6 * q ∧ ℓ = -3 * q := by
  rcases hypercharge_branches_from_anomaly_constraints hsu3 hsu2 hgrav hu1cube hq with hA | hB
  · rcases hA with ⟨hu2, hd4, he6, hℓ⟩
    have : ¬(u < d) := by
      rw [hu2, hd4]
      nlinarith [hqpos]
    exact False.elim (this hud)
  · exact hB

/-- With canonical normalization `q = 1/6`, the physical branch is exactly SM. -/
theorem normalized_physical_branch_is_standard_model
    {q u d ℓ e : ℚ}
    (hsu3 : u + d = -2 * q)
    (hsu2 : ℓ = -3 * q)
    (hgrav : 6 * q + 3 * u + 3 * d + 2 * ℓ + e = 0)
    (hu1cube : 6 * q ^ 3 + 3 * u ^ 3 + 3 * d ^ 3 + 2 * ℓ ^ 3 + e ^ 3 = 0)
    (hq : q ≠ 0)
    (hqnorm : q = 1 / 6)
    (hud : u < d) :
    u = -2 / 3 ∧ d = 1 / 3 ∧ ℓ = -1 / 2 ∧ e = 1 := by
  have hqpos : 0 < q := by
    rw [hqnorm]
    norm_num
  rcases physical_branch_selects_sm_sign hsu3 hsu2 hgrav hu1cube hq hqpos hud with ⟨hu, hd, he, hℓ⟩
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [hu, hqnorm]
    ring
  · rw [hd, hqnorm]
    ring
  · rw [hℓ, hqnorm]
    ring
  · rw [he, hqnorm]
    ring

/-- Canonical registry hypercharges satisfy the charged-sector anomaly constraints. -/
theorem canonical_assignment_satisfies_constraints :
    hypercharge .uRc + hypercharge .dRc = -2 * hypercharge .qL ∧
    hypercharge .lL = -3 * hypercharge .qL ∧
    6 * hypercharge .qL + 3 * hypercharge .uRc + 3 * hypercharge .dRc +
      2 * hypercharge .lL + hypercharge .eRc = 0 ∧
    6 * (hypercharge .qL) ^ 3 + 3 * (hypercharge .uRc) ^ 3 + 3 * (hypercharge .dRc) ^ 3 +
      2 * (hypercharge .lL) ^ 3 + (hypercharge .eRc) ^ 3 = 0 ∧
    hypercharge .qL ≠ 0 ∧
    hypercharge .qL = 1 / 6 ∧
    hypercharge .uRc < hypercharge .dRc := by
  constructor
  · norm_num [hypercharge, YqL, YuRc, YdRc]
  constructor
  · norm_num [hypercharge, YqL, YlL]
  constructor
  · norm_num [hypercharge, YqL, YuRc, YdRc, YlL, YeRc]
  constructor
  · norm_num [hypercharge, YqL, YuRc, YdRc, YlL, YeRc]
  constructor
  · norm_num [hypercharge, YqL]
  constructor
  · rfl
  · norm_num [hypercharge, YuRc, YdRc]

/--
Uniqueness bridge: under anomaly constraints + physical branch (`u < d`) +
canonical normalization (`q = Y_qL`), the charged-sector hypercharges are
exactly the canonical one-generation registry values.
-/
theorem canonical_hypercharge_unique_under_anomalies
    {q u d ℓ e : ℚ}
    (hsu3 : u + d = -2 * q)
    (hsu2 : ℓ = -3 * q)
    (hgrav : 6 * q + 3 * u + 3 * d + 2 * ℓ + e = 0)
    (hu1cube : 6 * q ^ 3 + 3 * u ^ 3 + 3 * d ^ 3 + 2 * ℓ ^ 3 + e ^ 3 = 0)
    (hqnorm : q = hypercharge .qL)
    (hud : u < d) :
    u = hypercharge .uRc ∧ d = hypercharge .dRc ∧
      ℓ = hypercharge .lL ∧ e = hypercharge .eRc := by
  have hqnorm' : q = 1 / 6 := by simpa [hypercharge, YqL] using hqnorm
  have hqpos : 0 < q := by
    rw [hqnorm']
    norm_num
  have hq : q ≠ 0 := ne_of_gt hqpos
  rcases normalized_physical_branch_is_standard_model
      hsu3 hsu2 hgrav hu1cube hq hqnorm' hud with ⟨hu, hd, hℓ, he⟩
  refine ⟨?_, ?_, ?_, ?_⟩
  · simpa [hypercharge, YuRc] using hu
  · simpa [hypercharge, YdRc] using hd
  · simpa [hypercharge, YlL] using hℓ
  · simpa [hypercharge, YeRc] using he

end Gutoe.SM.HyperchargeBridge
