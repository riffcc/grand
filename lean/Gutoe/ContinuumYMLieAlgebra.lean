/-
 * GUTOE — Continuum YM Lie-algebra scaffold
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND lane support:
 * Minimal compact-simple Lie group metadata used by continuum YM bundles.
 -/

import Mathlib

set_option autoImplicit false

namespace Gutoe

universe uG uLie

/-- Minimal metadata for a compact simple Lie group and its Lie algebra.

TODO (GRAND lane):
- refine compactness/simplicity as formal predicates,
- provide concrete instances (e.g. `SU n`, `SO n`),
- connect adjoint/coadjoint constructions used by YM bundles.
-/
structure CompactSimpleLieGroupData where
  G : Type uG
  𝔤 : Type uLie
  instGroup : Group G

attribute [instance] CompactSimpleLieGroupData.instGroup

end Gutoe
