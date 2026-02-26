// Debemus Server - Physics API Layer
// Copyright (C) 2026  Riff Labs
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

use axum::{routing::get, Json, Router};
use serde::Serialize;
use tower_http::cors::CorsLayer;
use tower_http::services::ServeDir;

// ── Shared JSON response envelope ────────────────────────────────────────────

#[derive(Serialize)]
struct Constant {
    name: &'static str,
    value: f64,
    derivation_method: &'static str,
    experimental_value: f64,
    error_percent: f64,
    lean_theorem: &'static str,
    proof_status: &'static str,
}

#[derive(Serialize)]
struct FlavorResponse {
    sector: &'static str,
    direct: gutoe_em::MixingObservables,
    texture: gutoe_em::MixingObservables,
    target: gutoe_em::MixingTargets,
    residuals_direct: gutoe_em::MixingResiduals,
    residuals_texture: gutoe_em::MixingResiduals,
}

#[derive(Serialize)]
struct LambdaResponse {
    lambda_base: LambdaStep,
    lorentz_factor: LambdaStep,
    lambda_signature: LambdaStep,
    micro_mode_count: LambdaStep,
    finite_mode_rescale: LambdaStep,
    lambda_full: LambdaStep,
    observed: LambdaStep,
}

#[derive(Serialize)]
struct LambdaStep {
    name: &'static str,
    value: f64,
    derivation: &'static str,
}

#[derive(Serialize)]
struct HubbleResponse {
    h0_km_s_mpc: f64,
    lambda_full: f64,
    omega_lambda: f64,
    omega_m0: f64,
    derivation_method: &'static str,
    experimental_value: f64,
    error_percent: f64,
    lean_theorem: &'static str,
    proof_status: &'static str,
}

// ── Z₃ orbit / Clifford algebra response types ──────────────────────────────

#[derive(Serialize)]
struct Z3Orbit {
    name: &'static str,
    orbit_type: &'static str,
    elements: Vec<&'static str>,
    physical_role: &'static str,
    grade: u32,
}

#[derive(Serialize)]
struct Z3OrbitsResponse {
    total_dimension: u32,
    singlet_count: u32,
    triplet_count: u32,
    orbits: Vec<Z3Orbit>,
    proof_status: &'static str,
    derivation: &'static str,
}

// ── Gauge generator response types ──────────────────────────────────────────

#[derive(Serialize)]
struct GaugeGroup {
    name: &'static str,
    dimension: u32,
    generators: Vec<&'static str>,
    source: &'static str,
    physical_role: &'static str,
}

#[derive(Serialize)]
struct GaugeGeneratorsResponse {
    groups: Vec<GaugeGroup>,
    total_generators: u32,
    charge_quantization: &'static str,
    anomaly_cancellation: &'static str,
    proof_status: &'static str,
    derivation: &'static str,
}

// ── Fermion content response types ──────────────────────────────────────────

#[derive(Serialize)]
struct SpinorDecomposition {
    name: &'static str,
    representation: &'static str,
    dimension: u32,
}

#[derive(Serialize)]
struct FermionContentResponse {
    generations: u32,
    generation_origin: &'static str,
    lepton_quark_distinction: &'static str,
    chirality_operator: &'static str,
    chirality_formula: &'static str,
    parity_violation: &'static str,
    spinor_decomposition: Vec<SpinorDecomposition>,
    total_spinor_dimension: u32,
    proof_status: &'static str,
    derivation: &'static str,
}

// ── Regge calculus / gravity response types ──────────────────────────────────

#[derive(Serialize)]
struct ReggeDecomposition {
    source: &'static str,
    tetrahedra: u32,
    edges: u32,
    coordination_number: u32,
}

#[derive(Serialize)]
struct ReggeGravityResponse {
    simplicial: ReggeDecomposition,
    deficit_angle_formula: &'static str,
    regge_action: &'static str,
    einstein_hilbert_convergence: &'static str,
    lambda_qg: f64,
    efe: &'static str,
    proof_status: &'static str,
    derivation: &'static str,
}

// ── Stellar fusion / PP-chain response types ────────────────────────────────

#[derive(Serialize)]
struct PPStep {
    step: u32,
    reaction: &'static str,
    description: &'static str,
}

#[derive(Serialize)]
struct SolarModelOutputs {
    core_temperature_k: f64,
    core_density_g_cm3: f64,
    specific_power_erg_g_s: f64,
    depletion_time_yr: f64,
}

#[derive(Serialize)]
struct PPChainResponse {
    steps: Vec<PPStep>,
    net_energy_mev: f64,
    net_energy_exact: &'static str,
    gamow_factor: &'static str,
    solar_model: SolarModelOutputs,
    alpha_connection: &'static str,
    proof_status: &'static str,
    derivation: &'static str,
}

// ── Yang-Mills mass gap / Doeblin response types ────────────────────────────

#[derive(Serialize)]
struct BridgeTheorem {
    index: u32,
    name: &'static str,
    statement: &'static str,
}

#[derive(Serialize)]
struct DoeblinMassGapResponse {
    alpha_inv: f64,
    doeblin_epsilon: f64,
    doeblin_formula: &'static str,
    bridge_theorems: Vec<BridgeTheorem>,
    mass_gap_lower_bound: &'static str,
    spectral_gap: &'static str,
    spectral_gap_numeric: f64,
    status: &'static str,
    plane_references: Vec<&'static str>,
    proof_status: &'static str,
    derivation: &'static str,
}

// ── Cosmological helpers (from frw_hz_report.rs) ─────────────────────────────

const METER_PER_MPC: f64 = 3.085_677_581_491_367e22;
const DEFAULT_OMEGA_M0: f64 = 0.315;
const DEFAULT_OMEGA_R0: f64 = 9.0e-5;
const DEFAULT_OMEGA_K0: f64 = 0.0;

fn s_inv_to_km_s_mpc(h0_s_inv: f64) -> f64 {
    h0_s_inv * METER_PER_MPC / 1_000.0
}

fn h0_from_lambda_and_omega_lambda(lambda: f64, omega_lambda: f64) -> f64 {
    let h0_s_inv =
        gutoe_physics::constants::C * (lambda / (3.0 * omega_lambda)).sqrt();
    s_inv_to_km_s_mpc(h0_s_inv)
}

// ── Endpoint handlers ────────────────────────────────────────────────────────

async fn health() -> &'static str {
    r#"{"status": "void", "coherence": 1.0, "veracity": 1.0}"#
}

async fn all_constants() -> Json<Vec<Constant>> {
    Json(vec![
        alpha_constant(),
        weinberg_constant(),
        mass_ratio_constant(),
        w_z_ratio_constant(),
        koide_constant(),
        higgs_quartic_constant(),
    ])
}

async fn alpha_endpoint() -> Json<Constant> {
    Json(alpha_constant())
}

async fn weinberg_endpoint() -> Json<Constant> {
    Json(weinberg_constant())
}

async fn mass_ratio_endpoint() -> Json<Constant> {
    Json(mass_ratio_constant())
}

async fn ckm_endpoint() -> Json<FlavorResponse> {
    let direct = gutoe_em::ckm_from_clifford();
    let texture = gutoe_em::ckm_from_textures();
    let target = gutoe_em::CKM_TARGET;
    Json(FlavorResponse {
        sector: "CKM",
        direct,
        texture,
        target,
        residuals_direct: gutoe_em::residuals(direct, target),
        residuals_texture: gutoe_em::residuals(texture, target),
    })
}

async fn pmns_endpoint() -> Json<FlavorResponse> {
    let direct = gutoe_em::pmns_from_clifford();
    let texture = gutoe_em::pmns_from_textures();
    let target = gutoe_em::PMNS_TARGET;
    Json(FlavorResponse {
        sector: "PMNS",
        direct,
        texture,
        target,
        residuals_direct: gutoe_em::residuals(direct, target),
        residuals_texture: gutoe_em::residuals(texture, target),
    })
}

async fn lambda_endpoint() -> Json<LambdaResponse> {
    use gutoe_physics::constants::*;

    let lambda_struct = lambda_cosmological_structural();
    let sig_factor = lorentz_signature_factor_from_bivector_split();
    let lambda_sig = lambda_cosmological_signature_candidate();
    let n_micro = lambda_micro_mode_count();
    let k_micro = lambda_micro_finite_mode_rescale();
    let lambda_full = lambda_cosmological_full_candidate();

    Json(LambdaResponse {
        lambda_base: LambdaStep {
            name: "Lambda_structural",
            value: lambda_struct,
            derivation: "lambda_H^(alpha_inv_LO) / l_P^2 = (13/100)^137 / l_P^2",
        },
        lorentz_factor: LambdaStep {
            name: "lorentz_signature_factor",
            value: sig_factor,
            derivation: "sqrt(bivector_total / bivector_timelike_spacelike) = sqrt(6/3) = sqrt(2)",
        },
        lambda_signature: LambdaStep {
            name: "Lambda_signature",
            value: lambda_sig,
            derivation: "Lambda_structural / sqrt(2)",
        },
        micro_mode_count: LambdaStep {
            name: "N_micro",
            value: n_micro,
            derivation: "ewsbScaleFactor + |grade-2| = 480 + 6 = 486 = 2 * 3^5",
        },
        finite_mode_rescale: LambdaStep {
            name: "k_micro",
            value: k_micro,
            derivation: "N_micro / (N_micro - 1) = 486/485",
        },
        lambda_full: LambdaStep {
            name: "Lambda_full",
            value: lambda_full,
            derivation: "Lambda_structural / sqrt(2) * 486/485",
        },
        observed: LambdaStep {
            name: "Lambda_observed",
            value: LAMBDA_COSMOLOGICAL_OBSERVED,
            derivation: "Planck 2018 + BAO (1/m^2)",
        },
    })
}

async fn hubble_endpoint() -> Json<HubbleResponse> {
    let lambda_full = gutoe_physics::constants::lambda_cosmological_full_candidate();
    let omega_lambda =
        1.0 - DEFAULT_OMEGA_M0 - DEFAULT_OMEGA_R0 - DEFAULT_OMEGA_K0;
    let h0 = h0_from_lambda_and_omega_lambda(lambda_full, omega_lambda);

    let experimental = 67.4_f64; // Planck 2018
    let error = ((h0 - experimental) / experimental).abs() * 100.0;

    Json(HubbleResponse {
        h0_km_s_mpc: h0,
        lambda_full,
        omega_lambda,
        omega_m0: DEFAULT_OMEGA_M0,
        derivation_method: "H0 = c * sqrt(Lambda_full / (3 * Omega_Lambda)), flat FRW with Omega_m = 0.315",
        experimental_value: experimental,
        error_percent: error,
        lean_theorem: "pending (FRW chain not yet in Lean)",
        proof_status: "computed",
    })
}

// ── Clifford Z₃ orbit structure ──────────────────────────────────────────────

async fn z3_orbits_endpoint() -> Json<Z3OrbitsResponse> {
    Json(Z3OrbitsResponse {
        total_dimension: 16,
        singlet_count: 4,
        triplet_count: 4,
        orbits: vec![
            Z3Orbit {
                name: "scalar singlet",
                orbit_type: "singlet",
                elements: vec!["1"],
                physical_role: "vacuum / identity",
                grade: 0,
            },
            Z3Orbit {
                name: "timelike singlet",
                orbit_type: "singlet",
                elements: vec!["\u{03b3}\u{2070}"],
                physical_role: "timelike vector / lepton seed",
                grade: 1,
            },
            Z3Orbit {
                name: "spatial volume singlet",
                orbit_type: "singlet",
                elements: vec!["\u{03b3}\u{00b9}\u{00b2}\u{00b3}"],
                physical_role: "spatial pseudoscalar",
                grade: 3,
            },
            Z3Orbit {
                name: "pseudoscalar singlet",
                orbit_type: "singlet",
                elements: vec!["\u{03b3}\u{2070}\u{00b9}\u{00b2}\u{00b3}"],
                physical_role: "chirality volume form",
                grade: 4,
            },
            Z3Orbit {
                name: "quark triplet",
                orbit_type: "triplet",
                elements: vec!["\u{03b3}\u{00b9}", "\u{03b3}\u{00b2}", "\u{03b3}\u{00b3}"],
                physical_role: "spatial vectors / color carriers",
                grade: 1,
            },
            Z3Orbit {
                name: "EM triplet",
                orbit_type: "triplet",
                elements: vec!["\u{03b3}\u{2070}\u{00b9}", "\u{03b3}\u{2070}\u{00b2}", "\u{03b3}\u{2070}\u{00b3}"],
                physical_role: "electric field / boosts",
                grade: 2,
            },
            Z3Orbit {
                name: "magnetic triplet",
                orbit_type: "triplet",
                elements: vec!["\u{03b3}\u{00b9}\u{00b2}", "\u{03b3}\u{00b2}\u{00b3}", "\u{03b3}\u{00b3}\u{00b9}"],
                physical_role: "magnetic field / SU(2) weak isospin generators",
                grade: 2,
            },
            Z3Orbit {
                name: "dual EM triplet",
                orbit_type: "triplet",
                elements: vec!["\u{03b3}\u{2070}\u{00b9}\u{00b2}", "\u{03b3}\u{2070}\u{00b2}\u{00b3}", "\u{03b3}\u{2070}\u{00b3}\u{00b9}"],
                physical_role: "Hodge duals of electric field",
                grade: 3,
            },
        ],
        proof_status: "proven",
        derivation: "Z_3 automorphism of Cl(1,3) partitions 16 basis elements into 4 singlets + 4 triplets; sin^2(theta_W) = 3/13 follows from |magnetic triplet|/(dim - |magnetic triplet|)",
    })
}

// ── Gauge generators from Clifford algebra ──────────────────────────────────

async fn gauge_generators_endpoint() -> Json<GaugeGeneratorsResponse> {
    Json(GaugeGeneratorsResponse {
        groups: vec![
            GaugeGroup {
                name: "SU(3)",
                dimension: 8,
                generators: vec![
                    "\u{03b3}\u{2070}\u{00b9}", "\u{03b3}\u{2070}\u{00b2}", "\u{03b3}\u{2070}\u{00b3}",
                    "\u{03b3}\u{00b9}\u{00b2}", "\u{03b3}\u{00b2}\u{00b3}", "\u{03b3}\u{00b3}\u{00b9}",
                    "even-subalgebra composite 1",
                    "even-subalgebra composite 2",
                ],
                source: "6 grade-2 bivectors + 2 composite generators from Cl+(1,3) even subalgebra",
                physical_role: "strong force / color confinement",
            },
            GaugeGroup {
                name: "SU(2)",
                dimension: 3,
                generators: vec![
                    "\u{03b3}\u{00b9}\u{00b2}", "\u{03b3}\u{00b2}\u{00b3}", "\u{03b3}\u{00b3}\u{00b9}",
                ],
                source: "magnetic triplet from Z_3 orbit decomposition",
                physical_role: "weak isospin",
            },
            GaugeGroup {
                name: "U(1)",
                dimension: 1,
                generators: vec![
                    "Y = linear combination of EM triplet",
                ],
                source: "linear combination of {gamma^{01}, gamma^{02}, gamma^{03}} = hypercharge direction",
                physical_role: "hypercharge",
            },
        ],
        total_generators: 12,
        charge_quantization: "Q = T_3 + Y/2",
        anomaly_cancellation: "automatic from algebraic structure",
        proof_status: "proven",
        derivation: "Gauge group 8+3+1=12 generators emerge from Cl(1,3) grade-2 bivector structure and Z_3 orbit decomposition; dim(SU(3))+dim(SU(2))+dim(U(1)) = 8+3+1 = 12",
    })
}

// ── Fermion content ─────────────────────────────────────────────────────────

async fn fermion_content_endpoint() -> Json<FermionContentResponse> {
    Json(FermionContentResponse {
        generations: 3,
        generation_origin: "Z_3 automorphism of Cl(1,3): three-fold cyclic symmetry generates exactly 3 fermion generations",
        lepton_quark_distinction: "singlet orbits -> leptons, triplet orbits -> quarks (color triplets under SU(3))",
        chirality_operator: "\u{03b3}\u{2075}",
        chirality_formula: "\u{03b3}\u{2075} = \u{03b3}\u{2070}\u{03b3}\u{00b9}\u{03b3}\u{00b2}\u{03b3}\u{00b3}",
        parity_violation: "built into algebra: left-handed and right-handed spinors live in different irreducible representations of Cl(1,3)",
        spinor_decomposition: vec![
            SpinorDecomposition {
                name: "Delta_plus",
                representation: "(4,1)",
                dimension: 4,
            },
            SpinorDecomposition {
                name: "Delta_minus",
                representation: "(1,4)",
                dimension: 4,
            },
        ],
        total_spinor_dimension: 8,
        proof_status: "structural",
        derivation: "Spinor space of Cl(1,3) decomposes as Delta+ direct-sum Delta- = (4,1) direct-sum (1,4) under chirality; Z_3 generates 3 copies (generations)",
    })
}

// ── Regge calculus / discrete gravity ────────────────────────────────────────

async fn regge_gravity_endpoint() -> Json<ReggeGravityResponse> {
    Json(ReggeGravityResponse {
        simplicial: ReggeDecomposition {
            source: "1 cube decomposed into simplices",
            tetrahedra: 6,
            edges: 19,
            coordination_number: 6,
        },
        deficit_angle_formula: "\u{03b4}_e = 2\u{03c0} - \u{03a3} \u{03b8}_t (sum of dihedral angles around edge e)",
        regge_action: "S_Regge = (1/8\u{03c0}G) \u{03a3}_e A_e \u{03b4}_e",
        einstein_hilbert_convergence: "Regge action converges to Einstein-Hilbert action S_EH = (1/16piG) integral(R sqrt(g) d^4x) in the continuum limit as lattice spacing -> 0",
        lambda_qg: gutoe_physics::constants::LAMBDA_QG,
        efe: "R_\u{03bc}\u{03bd} - \u{00bd}g_\u{03bc}\u{03bd}R + \u{039b}g_\u{03bc}\u{03bd} = 8\u{03c0}G T_\u{03bc}\u{03bd}",
        proof_status: "structural",
        derivation: "Regge calculus discretizes GR on simplicial complexes; lambda_QG = 1/12 from Planck-lattice dispersion k^4 coefficient; 1 cube -> 6 tetrahedra with 19 edges and coordination number 6",
    })
}

// ── Stellar fusion PP-chain ─────────────────────────────────────────────────

async fn pp_chain_endpoint() -> Json<PPChainResponse> {
    Json(PPChainResponse {
        steps: vec![
            PPStep {
                step: 1,
                reaction: "p + p \u{2192} d + e\u{207a} + \u{03bd}_e",
                description: "proton-proton fusion: weak interaction converts proton to neutron, producing deuterium, positron, and electron neutrino",
            },
            PPStep {
                step: 2,
                reaction: "p + d \u{2192} \u{00b3}He + \u{03b3}",
                description: "proton-deuterium fusion: strong interaction produces helium-3 and gamma ray",
            },
            PPStep {
                step: 3,
                reaction: "\u{00b3}He + \u{00b3}He \u{2192} \u{2074}He + 2p",
                description: "helium-3 fusion: two helium-3 nuclei fuse to produce helium-4 and two protons (PP-I termination)",
            },
        ],
        net_energy_mev: 26.732,
        net_energy_exact: "6683/250 MeV",
        gamow_factor: "P ~ exp(-b/sqrt(E)), where b = pi * alpha * sqrt(2 * m_reduced * c^2)",
        solar_model: SolarModelOutputs {
            core_temperature_k: 1.57e7,
            core_density_g_cm3: 150.0,
            specific_power_erg_g_s: 2.7e-4,
            depletion_time_yr: 5.0e10,
        },
        alpha_connection: "Gamow penetration factor depends on alpha = 1/137; tunneling probability P ~ exp(-2*pi*alpha*sqrt(m_reduced*c^2/(2*E))), making fusion rate exquisitely sensitive to the fine-structure constant",
        proof_status: "computed",
        derivation: "PP-I chain: 4p -> He-4 + 2e+ + 2nu_e + 2gamma; Q = 26.732 MeV from mass deficit; single-zone solar model with Lane-Emden polytrope",
    })
}

// ── Yang-Mills mass gap via Doeblin coupling ────────────────────────────────

async fn doeblin_massgap_endpoint() -> Json<DoeblinMassGapResponse> {
    let alpha_inv = gutoe_physics::constants::ALPHA_INV_LEADING_ORDER as f64;
    let alpha = 1.0 / alpha_inv;
    let epsilon = 3.0 * alpha / (6.0 + 3.0 * alpha);
    let spectral_gap = epsilon / (1.0 - epsilon);

    Json(DoeblinMassGapResponse {
        alpha_inv,
        doeblin_epsilon: epsilon,
        doeblin_formula: "\u{03b5} = 3\u{03b1}/(6 + 3\u{03b1}) where \u{03b1} = 1/137",
        bridge_theorems: vec![
            BridgeTheorem {
                index: 1,
                name: "Wilson to Z_3",
                statement: "SU(3) Wilson lattice gauge theory in strong coupling maps to Z_3 spin model with bounded error",
            },
            BridgeTheorem {
                index: 2,
                name: "Z_3 Markov chain",
                statement: "Z_3 lattice model transfer matrix defines an irreducible aperiodic Markov chain satisfying Doeblin condition with epsilon > 0",
            },
            BridgeTheorem {
                index: 3,
                name: "spectral gap inheritance",
                statement: "Doeblin coefficient epsilon > 0 implies spectral gap Delta >= epsilon/(1-epsilon) > 0 for the transfer matrix, yielding exponential decay of correlations and mass gap",
            },
        ],
        mass_gap_lower_bound: "m_gap >= (epsilon/(1-epsilon)) * Lambda_QCD",
        spectral_gap: "\u{0394} >= \u{03b5}/(1-\u{03b5})",
        spectral_gap_numeric: spectral_gap,
        status: "structural pathway, not yet full proof",
        plane_references: vec!["GRAND-308", "GRAND-309", "GRAND-310", "GRAND-311"],
        proof_status: "structural",
        derivation: "Doeblin coupling argument: Z_3 transfer matrix from Cl(1,3) orbit structure has epsilon = 3*alpha/(6+3*alpha); spectral gap Delta >= epsilon/(1-epsilon) implies mass gap m >= Delta * Lambda_QCD",
    })
}

// ── Static file serving ─────────────────────────────────────────────────────

async fn serve_index() -> axum::response::Html<String> {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let path = format!("{manifest_dir}/../frontend/index.html");
    let html = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("Error reading {path}: {e}");
        "<html><body><h1>Void</h1></body></html>".to_string()
    });
    axum::response::Html(html)
}

// ── Constant builders ────────────────────────────────────────────────────────

fn alpha_constant() -> Constant {
    Constant {
        name: "alpha_inverse",
        value: gutoe_physics::constants::ALPHA_INV_LEADING_ORDER as f64,
        derivation_method: "T(16) + 1 = 16*17/2 + 1 = 136 + 1 = 137, where T(n) = n(n+1)/2 and 16 = dim Cl(1,3)",
        experimental_value: 137.035_999_084,
        error_percent: ((137.0_f64 - 137.035_999_084) / 137.035_999_084).abs() * 100.0,
        lean_theorem: "FineStructure.lean: alpha_inv_eq_137",
        proof_status: "proven",
    }
}

fn weinberg_constant() -> Constant {
    let predicted = gutoe_em::sin2_weinberg();
    let experimental = 0.23122_f64;
    Constant {
        name: "sin2_weinberg",
        value: predicted,
        derivation_method: "3/13: from Cl(1,3) Z3 orbits: |magnetic-triplet| / (2^4 - |magnetic-triplet|) = 3/(16-3)",
        experimental_value: experimental,
        error_percent: ((predicted - experimental) / experimental).abs() * 100.0,
        lean_theorem: "Weinberg.lean: weinberg_from_z3_orbits",
        proof_status: "proven",
    }
}

fn mass_ratio_constant() -> Constant {
    let predicted = 1836.0_f64;
    let experimental = 1836.152_673_43_f64;
    Constant {
        name: "mp_me_ratio",
        value: predicted,
        derivation_method: "12 * T(17) = 12 * 153 = 1836, where 12 = dim(SU(3))+dim(SU(2))+dim(U(1)), T(17) = 17*18/2",
        experimental_value: experimental,
        error_percent: ((predicted - experimental) / experimental).abs() * 100.0,
        lean_theorem: "MassRatio.lean: mp_me_eq_1836",
        proof_status: "proven",
    }
}

fn w_z_ratio_constant() -> Constant {
    let predicted = gutoe_em::w_z_mass_ratio();
    let experimental = 80.377 / 91.1876;
    Constant {
        name: "mw_mz_ratio",
        value: predicted,
        derivation_method: "sqrt(10/13): from cos(theta_W) = sqrt(1 - 3/13), where 13 = 16-3, 10 = 13-3",
        experimental_value: experimental,
        error_percent: ((predicted - experimental) / experimental).abs() * 100.0,
        lean_theorem: "Weinberg.lean: w_z_mass_ratio_from_weinberg",
        proof_status: "proven",
    }
}

fn koide_constant() -> Constant {
    let predicted = 2.0 / 3.0;
    let experimental = 0.666_61_f64;
    Constant {
        name: "koide_ratio",
        value: predicted,
        derivation_method: "grade-1/grade-2 = 4/6 = 2/3: from Cl(1,3) grade counts",
        experimental_value: experimental,
        error_percent: ((predicted - experimental) / experimental).abs() * 100.0,
        lean_theorem: "Koide.lean: koide_from_grade_counts",
        proof_status: "proven",
    }
}

fn higgs_quartic_constant() -> Constant {
    let predicted = gutoe_physics::constants::HIGGS_QUARTIC_STRUCTURAL;
    let experimental = 0.129_f64;
    Constant {
        name: "higgs_quartic_lambda",
        value: predicted,
        derivation_method: "(16-3)/(4+6)^2 = 13/100: from Cl(1,3) dimension minus SU(2) over (grade-1+grade-2)^2",
        experimental_value: experimental,
        error_percent: ((predicted - experimental) / experimental).abs() * 100.0,
        lean_theorem: "HiggsQuartic.lean: higgs_quartic_structural",
        proof_status: "proven",
    }
}

// ── Main ─────────────────────────────────────────────────────────────────────

#[tokio::main]
async fn main() {
    let cors = CorsLayer::permissive();

    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let frontend_dir = format!("{manifest_dir}/../frontend");

    let app = Router::new()
        .route("/", get(serve_index))
        .route("/api/health", get(health))
        .route("/api/constants", get(all_constants))
        .route("/api/constants/alpha", get(alpha_endpoint))
        .route("/api/constants/weinberg", get(weinberg_endpoint))
        .route("/api/constants/mass-ratio", get(mass_ratio_endpoint))
        .route("/api/flavor/ckm", get(ckm_endpoint))
        .route("/api/flavor/pmns", get(pmns_endpoint))
        .route("/api/cosmology/lambda", get(lambda_endpoint))
        .route("/api/cosmology/hubble", get(hubble_endpoint))
        .route("/api/clifford/z3-orbits", get(z3_orbits_endpoint))
        .route("/api/gauge/generators", get(gauge_generators_endpoint))
        .route("/api/fermions/content", get(fermion_content_endpoint))
        .route("/api/gravity/regge", get(regge_gravity_endpoint))
        .route("/api/stellar/pp-chain", get(pp_chain_endpoint))
        .route("/api/massgap/doeblin", get(doeblin_massgap_endpoint))
        .nest_service("/static", ServeDir::new(&frontend_dir))
        .layer(cors);

    let listener = tokio::net::TcpListener::bind("0.0.0.0:9000")
        .await
        .unwrap();
    println!("Debemus Server running at http://localhost:9000");
    println!("  Omnia cognoscibilia. Omnia computabilia. Possumus. Debemus.");
    println!("  API: /api/constants, /api/flavor/ckm, /api/cosmology/lambda, ...");
    println!("  Additional: /api/clifford/z3-orbits, /api/gauge/generators, /api/fermions/content");
    println!("  Additional: /api/gravity/regge, /api/stellar/pp-chain, /api/massgap/doeblin");

    axum::serve(listener, app).await.unwrap();
}
