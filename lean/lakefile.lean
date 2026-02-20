import Lake
open Lake DSL

package gutoe

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git"

lean_lib Gutoe where
  roots := #[`Gutoe.Basic, `Gutoe.HexStates, `Gutoe.GateProperties, `Gutoe.RealGates,
             `Gutoe.ParticleFormation, `Gutoe.DispersionRelation, `Gutoe.RailSpace,
             `Gutoe.Spacetime, `Gutoe.HawkingCorrection]
