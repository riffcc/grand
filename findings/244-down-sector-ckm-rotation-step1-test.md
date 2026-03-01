# 244 — Down-Sector CKM Rotation Step-1 Test

## Scope
Run the exact requested Step-1 test:

`d_mass = V_CKM · d_weak  =>  d_weak = V_CKM† · d_mass`

using structural CKM from `ckm_from_clifford`, then re-evaluate down-sector closure.

Implementation:
- `crates/gutoe-em/src/bin/yukawa_down_ckm_rotate_test.rs`

Artifacts:
- `/tmp/bh_renders/yukawa_down_ckm_rotate_test.txt`
- `/tmp/bh_renders/yukawa_down_ckm_rotate_test.json`

## Method (fixed, no optimizer)
1. Use structural CKM angles/phases:
   - `theta12 = 13.262676°`
   - `theta23 = 2.388015°`
   - `theta13 = 0.210647°`
   - `delta = 68.130102°`
2. Build PDG-form CKM matrix `V`.
3. Treat down masses in amplitude space:
   - `a_mass = [sqrt(m_d), sqrt(m_s), sqrt(m_b)]`.
4. Rotate:
   - `a_weak = V† a_mass`.
5. Map back to masses:
   - `m_weak_i = |a_weak_i|^2`.
6. Compare pre/post against structural targets (`m_s/m_d`, `m_b/m_s`) and cross checks (`m_c/m_s`, `m_t/m_b`).

## Results

Mass basis (input):
- `[4.67, 93.0, 4180.0] MeV`
- `m_s/m_d = 19.9143`
- `m_b/m_s = 44.9462`
- down-closure RMS-log vs structural targets: `0.100972`

Rotated weak basis (`V†`):
- `[0.226588, 52.298824, 4225.144588] MeV`
- `m_s/m_d = 230.8104`
- `m_b/m_s = 80.7885`
- down-closure RMS-log vs structural targets: `1.794391`

Cross ratios:
- `m_c/m_s`: `13.6559 -> 24.2835` (worse)
- `m_t/m_b`: `41.3301 -> 40.8885` (slight drift, not improvement driver)

Boolean verdict:
- `closure_improved = false`

## Conclusion
In this direct amplitude-space realization of the hypothesis, Step-1 **fails**:
`V_CKM†` rotation does not collapse down-sector closure; it strongly degrades it.

This rules out the naive form:
- “apply CKM dagger directly to down mass amplitudes and read off corrected Z3 spectrum.”

## Practical implication
The next correction lane should not be another scalar tweak. It needs a different map between:
- Z3 harmonic parameters in weak space,
- CKM mixing structure,
- observed mass eigenvalues.

The failure here is useful: it removes one clean-but-wrong bridge and narrows the search.
