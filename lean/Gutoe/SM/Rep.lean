/- 
 * GUTOE — Standard Model Chiral Representation Registry
 *
 * Canonical one-generation Weyl content used by anomaly and flavor proofs.
 * This is the single source of truth for representation bookkeeping.
-/

import Mathlib
import Gutoe.Z3Uniqueness

namespace Gutoe.SM.Rep

open Gutoe.Z3Uniqueness

/-- Color representation class for one Weyl species. -/
inductive ColorRep
  | singlet
  | triplet
  | antiTriplet
  deriving DecidableEq, Repr

/-- Weak representation class for one Weyl species. -/
inductive WeakRep
  | singlet
  | doublet
  deriving DecidableEq, Repr

/--
Canonical left-chiral one-generation Weyl species:
- `qL`   : quark doublet
- `uRc`  : charge-conjugated up singlet
- `dRc`  : charge-conjugated down singlet
- `lL`   : lepton doublet
- `eRc`  : charge-conjugated electron singlet
- `nuRc` : optional right-handed neutrino (hypercharge 0, anomaly-neutral)
-/
inductive WeylSpecies
  | qL
  | uRc
  | dRc
  | lL
  | eRc
  | nuRc
  deriving DecidableEq, Fintype, Repr

/-- Canonical species set for one generation. -/
def oneGeneration : Finset WeylSpecies := Finset.univ

/-- Color representation assignment. -/
def colorRep : WeylSpecies → ColorRep
  | .qL => .triplet
  | .uRc => .antiTriplet
  | .dRc => .antiTriplet
  | .lL => .singlet
  | .eRc => .singlet
  | .nuRc => .singlet

/-- Weak representation assignment. -/
def weakRep : WeylSpecies → WeakRep
  | .qL => .doublet
  | .uRc => .singlet
  | .dRc => .singlet
  | .lL => .doublet
  | .eRc => .singlet
  | .nuRc => .singlet

/-- Hypercharge assignment for left-chiral Weyl fields. -/
def YqL : ℚ := 1 / 6
def YuRc : ℚ := -2 / 3
def YdRc : ℚ := 1 / 3
def YlL : ℚ := -1 / 2
def YeRc : ℚ := 1
def YnuRc : ℚ := 0
def YH : ℚ := 1 / 2

/-- Hypercharge assignment for left-chiral Weyl fields. -/
def hypercharge : WeylSpecies → ℚ
  | .qL => YqL
  | .uRc => YuRc
  | .dRc => YdRc
  | .lL => YlL
  | .eRc => YeRc
  | .nuRc => YnuRc

/-- Color multiplicity from Cl(1,3) Z₃ orbit structure. -/
def colorMultiplicity : ColorRep → ℕ
  | .singlet => 1
  | .triplet => magneticTriplet.card
  | .antiTriplet => magneticTriplet.card

/-- Weak multiplicity for species counting. -/
def weakMultiplicity : WeakRep → ℕ
  | .singlet => 1
  | .doublet => 2

/-- Fundamental Dynkin index normalization for SU(3). -/
def dynkinSU3 : ColorRep → ℚ
  | .singlet => 0
  | .triplet => 1 / 2
  | .antiTriplet => 1 / 2

/-- Fundamental Dynkin index normalization for SU(2). -/
def dynkinSU2 : WeakRep → ℚ
  | .singlet => 0
  | .doublet => 1 / 2

/-- Color-copy multiplicity for a Weyl species. -/
def colorMultiplicitySpecies (f : WeylSpecies) : ℕ := colorMultiplicity (colorRep f)

/-- Weak-component multiplicity for a Weyl species. -/
def weakMultiplicitySpecies (f : WeylSpecies) : ℕ := weakMultiplicity (weakRep f)

/-- Shared Cl(1,3) fact: magnetic triplet has cardinality 3. -/
theorem magnetic_triplet_card_eq_three : magneticTriplet.card = 3 := by
  decide

/-- One-generation Weyl registry has six species. -/
theorem one_generation_card : oneGeneration.card = 6 := by
  native_decide

end Gutoe.SM.Rep
