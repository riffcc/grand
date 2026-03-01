# 198 — Creation lanes tested: affine offset + compact local topology patch

## Prompted objective
Test the remaining creation lanes after linear and kernel-homogeneous no-shortcut closures.

## Lean additions
- `Gutoe.CreationLanes` (new)
  - `affine_origin_descends_as_offset`
  - `affine_origin_reaches_target_if_offset_projects`
  - `sameOnLocalPatch`
  - `local_patch_nontrivial_shift_exists`
  - `local_patch_no_nontrivial_shift_outside`

Interpretation:
- Affine/inhomegeneous lane is formally distinct from homogeneous linear lane.
- Compact-support local identification can be represented and constrained in-theory.

## Rust probes

### 1) Affine creation probe
Bin: `recursive_z3_affine_creation_probe`

Artifact:
- `/tmp/bh_renders/recursive_z3_affine_creation_probe/recursive_z3_affine_creation_probe.{txt,json}`

Result highlights:
- `unit_shift`: shift norm `1.0`, witness affine norm `1.0`
- `mixed_shift`: shift norm `9.486832980505`, witness affine norm `9.486832980505`
- `large_shift`: shift norm `1e6`, witness affine norm `1e6`

Across Monte Carlo feasible offsets, witness (section, free coords zero) remained minimal.

Conclusion:
- Affine lane is open by construction; cost floor scales with required 4D shift norm.

### 2) Compact local topology patch probe
Bin: `ctc_local_patch_creation_probe`

Artifact:
- `/tmp/bh_renders/ctc_local_patch_creation_probe/ctc_local_patch_creation_probe.{txt,json}`

Toy-model assumptions:
- Timelike local travel with `c=1`
- Identification active only for `|x| <= R`
- Each loop contributes effective coordinate shift `-T`

Result highlights:
- All tested cases found effective-superluminal and pre-departure coordinate arrivals
  when loop count `n` is sufficiently large.

Conclusion:
- In this compact-support identification toy model, creation lane is operationally open.
- This is a model result, not a physical realizability proof.

## Overall
Creation lanes now split cleanly:
- **Affine inhomogeneous lane:** mathematically open with explicit cost floor.
- **Local topological patch lane (toy quotient model):** operationally open in simulation.
- Physical instantiation mechanism remains open and is the next falsification target.
