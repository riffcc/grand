import Lake
open Lake DSL

package gutoe

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git"

lean_lib Gutoe where
  roots := #[`Gutoe.Basic, `Gutoe.HexStates, `Gutoe.GateProperties, `Gutoe.RealGates,
             `Gutoe.ParticleFormation, `Gutoe.DispersionRelation, `Gutoe.RailSpace,
             `Gutoe.Spacetime, `Gutoe.HawkingCorrection, `Gutoe.CliffordStructure,
             `Gutoe.BaryonPhysics, `Gutoe.HexFermions, `Gutoe.DynamicsProperties,
             `Gutoe.TopologyForces, `Gutoe.Conjectures, `Gutoe.HydrogenFormation,
             `Gutoe.FineStructure, `Gutoe.MassSpectrum]
