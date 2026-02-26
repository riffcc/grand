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
             `Gutoe.FineStructure, `Gutoe.MassSpectrum, `Gutoe.DimensionalStructure,
             `Gutoe.Z3Uniqueness, `Gutoe.BellInequality, `Gutoe.ThreeGenerations,
             `Gutoe.KoideMasses, `Gutoe.LeptonMass, `Gutoe.PerturbativeSymmetry,
             `Gutoe.FlavorMixing,
             `Gutoe.SignatureUniqueness, `Gutoe.LatticeGeometry, `Gutoe.InstantonMass,
             `Gutoe.Z3ForcedStructure,
             `Gutoe.GaugeGroupSU2, `Gutoe.GaugeGroupSU3, `Gutoe.GaugeGroupSM,
             `Gutoe.GaugeConstants, `Gutoe.Chirality, `Gutoe.ContinuumLimit,
             `Gutoe.LorentzInvariance, `Gutoe.GrandMasterTheorem,
             `Gutoe.GravityMetric,
             `Gutoe.EinsteinFromLattice,
             `Gutoe.KerrGeometry,
             `Gutoe.KerrCameraStability,
             `Gutoe.KerrTracerEquations,
             `Gutoe.KerrReferenceFrame,
             `Gutoe.Geodesic3DProjection,
             `Gutoe.SynchrotronGRMHD,
             `Gutoe.SynchrotronTransfer,
             `Gutoe.StrongCP,
             `Gutoe.StrongCPPathIntegral,
             `Gutoe.StrongCPVacuum,
             `Gutoe.StrongCPEmergence,
             `Gutoe.StrongCPGeneralCases,
             `Gutoe.SMQCDUnification,
             `Gutoe.EWSBHiggs,
             `Gutoe.CosmologicalConstant,
             `Gutoe.StellarFusion,
             `Gutoe.FalsifiabilityCatalog,
             `Gutoe.AsymptoticFreedomEntropy,
             `Gutoe.YangMillsMassGap,
             `Gutoe.YangMillsContinuumSurvival,
             `Gutoe.YangMillsStructuralGap,
             `Gutoe.YangMillsWilsonBridge,
             `Gutoe.SM.Rep,
             `Gutoe.SM.Anomalies,
             `Gutoe.SM.HyperchargeBridge,
             `Gutoe.SM.Closure,
             `Gutoe.LambdaQG]

lean_exe kerr_ref_frame where
  root := `Gutoe.KerrReferenceFrameGen
