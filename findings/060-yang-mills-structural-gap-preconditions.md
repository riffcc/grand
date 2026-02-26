# Finding 060 — Yang-Mills Structural Gap Preconditions (Lean)

Date: 2026-02-26

Scope: GRAND-296/297/298 Theorem A track

## Added module

- `lean/Gutoe/YangMillsStructuralGap.lean`
- Added to roots in `lean/lakefile.lean`

## What is now proven (no `sorry`)

1. **Cl(1,3) basis lock**
   - `transfer_basis_dim_eq_three`
   - Transfer basis dimension is fixed at 3 via `magneticTriplet.card = 3`.

2. **Laplace smoothing positivity (structural, all volumes)**
   - `smooth_entry_pos`
   - `smoothed_transition_entry_pos`
   - For any `alpha > 0`, every smoothed transition entry is strictly positive.

3. **Primitive + irreducible from entrywise positivity**
   - `isPrimitive_of_entrywise_pos`
   - `isIrreducible_of_entrywise_pos`
   - For `3×3` real matrices, entrywise positivity is enough to prove primitivity (choose `k=1`) and irreducibility.

4. **Smoothed transition kernel is primitive/irreducible**
   - `smoothed_transition_isPrimitive`
   - `smoothed_transition_isIrreducible`
   - This discharges a key structural precondition for Perron-Frobenius style gap arguments.

5. **Gram/symmetric transfer proxy inherits positivity**
   - `gram_entry_pos_of_entrywise_pos`
   - `gram_isPrimitive_of_entrywise_pos`
   - `gram_isIrreducible_of_entrywise_pos`
   - If `S` is entrywise positive, then `S * Sᵀ` is entrywise positive, primitive, and irreducible.

## Why this matters

This upgrades Theorem A from purely numerical ordering checks to structural matrix-theory guarantees tied to construction:

- positivity from smoothing (`alpha>0`)
- primitivity/irreducibility of transfer objects from that positivity

What remains for full Theorem A closure is the explicit **spectral dominance** step (`λ₀ > |λ₁|`) from primitivity/irreducibility (Perron-Frobenius spectral theorem bridge).

## Build verification

- `cd lean && lake build Gutoe.YangMillsStructuralGap` ✅
- `cd lean && lake build Gutoe` ✅
