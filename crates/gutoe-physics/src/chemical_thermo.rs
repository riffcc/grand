/*!
 * Chemical thermodynamics proxy lane for periodic-table scale modeling.
 *
 * This lane uses physically motivated, explicit assumptions to estimate
 * bulk thermodynamic properties per element:
 * - cohesive-energy proxy
 * - Debye-temperature proxy
 * - latent heats
 * - melting/boiling temperatures
 * - phase at ambient conditions
 *
 * It is not a quantum-chemistry solver; the formulas are closed-form
 * transduction rules intended for broad trend modeling and extrapolation.
 */

use crate::ab_initio_qchem::{predict_atomic_scf, AtomicScfPrediction};
use std::f64::consts::PI;
use std::sync::OnceLock;

pub const AVOGADRO: f64 = 6.022_140_76e23;
pub const R_GAS_J_MOL_K: f64 = 8.314_462_618;
pub const BOHR_RADIUS_PM: f64 = 52.917_721;
pub const EV_TO_KJ_MOL: f64 = 96.485_332_123;
pub const ENTROPY_FUSION_J_MOL_K: f64 = 10.0;
pub const ENTROPY_VAPORIZATION_J_MOL_K: f64 = 85.0;
pub const P_REF_PA: f64 = 101_325.0;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ChemicalFamily {
    Alkali,
    AlkalineEarth,
    Transition,
    PostTransition,
    Metalloid,
    Nonmetal,
    Halogen,
    NobleGas,
    Lanthanide,
    Actinide,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MatterState {
    Solid,
    Liquid,
    Gas,
}

#[derive(Clone, Copy, Debug)]
pub struct PhaseGibbsKjMol {
    pub solid: f64,
    pub liquid: f64,
    pub gas: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct ElementThermoPrediction {
    pub z: u16,
    pub a: u16,
    pub family: ChemicalFamily,
    pub period: u8,
    pub molar_mass_g_mol: f64,
    pub atomic_radius_pm: f64,
    pub molar_volume_cm3_mol: f64,
    pub density_g_cm3: f64,
    pub cohesive_energy_ev_per_atom: f64,
    pub debye_temperature_k: f64,
    pub latent_fusion_kj_mol: f64,
    pub latent_vaporization_kj_mol: f64,
    pub melting_temperature_k: f64,
    pub boiling_temperature_k: f64,
    pub vapor_pressure_pa_298k: f64,
    pub cp_solid_j_mol_k: f64,
    pub cp_liquid_j_mol_k: f64,
    pub cp_gas_j_mol_k: f64,
    pub bulk_modulus_gpa: f64,
    pub thermal_expansion_1_per_k: f64,
    pub ambient_state_298k: MatterState,
}

#[derive(Clone, Copy, Debug)]
pub struct CoupledThermoDiagnostics {
    pub scf_iterations: usize,
    pub scf_residual: f64,
    pub valence_electrons: u16,
    pub scf_atomic_radius_pm: f64,
    pub scf_ionization_energy_ev: f64,
    pub scf_electron_affinity_ev: f64,
    pub scf_electronegativity_mulliken_ev: f64,
    pub scf_chemical_hardness_ev: f64,
    pub coupled_radius_pm: f64,
    pub cohesive_frontier_proxy_ev: f64,
    pub coupled_packing_fraction: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct ChemicalThermoCalibration {
    pub pack_p_void_coef: f64,
    pub pack_d_gain_coef: f64,
    pub pack_f_gain_coef: f64,
    pub pack_open_d_mult: f64,
    pub pack_closed_d_mult: f64,
    pub pack_f_core_mult: f64,
    pub radius_p_gain: f64,
    pub radius_closed_d_mult: f64,
    pub radius_open_d_mult: f64,
    pub radius_f_core_mult: f64,
    pub radius_actinide_mult: f64,
    pub radius_lower_actinide: f64,
    pub radius_lower_transition_fcore: f64,
}

impl Default for ChemicalThermoCalibration {
    fn default() -> Self {
        Self {
            pack_p_void_coef: 0.52,
            pack_d_gain_coef: 0.22,
            pack_f_gain_coef: 0.22,
            pack_open_d_mult: 1.20,
            pack_closed_d_mult: 0.70,
            pack_f_core_mult: 1.12,
            radius_p_gain: 0.35,
            radius_closed_d_mult: 1.08,
            radius_open_d_mult: 0.84,
            radius_f_core_mult: 0.74,
            radius_actinide_mult: 0.62,
            radius_lower_actinide: 0.35,
            radius_lower_transition_fcore: 0.55,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct OrbitalPackingHints {
    p_frac: f64,
    d_frac: f64,
    f_frac: f64,
    d_fill: f64,
    valence_electrons: u16,
    open_d_shell: bool,
    closed_d_shell: bool,
    has_f_core: bool,
}

#[derive(Clone, Debug)]
pub struct CoupledThermoPrefetch {
    pub z: u16,
    pub a: u16,
    pub family: ChemicalFamily,
    pub period: u8,
    pub base: ElementThermoPrediction,
    pub scf: AtomicScfPrediction,
    pub p_frac: f64,
    pub d_frac: f64,
    pub f_frac: f64,
    pub d_fill: f64,
    pub valence_electrons: u16,
    pub open_d_shell: bool,
    pub closed_d_shell: bool,
    pub has_f_core: bool,
    pub scf_radius_weight: f64,
}

pub fn family_of_z(z: u16) -> ChemicalFamily {
    match z {
        1 | 3 | 11 | 19 | 37 | 55 | 87 => ChemicalFamily::Alkali,
        4 | 12 | 20 | 38 | 56 | 88 => ChemicalFamily::AlkalineEarth,
        2 | 10 | 18 | 36 | 54 | 86 | 118 => ChemicalFamily::NobleGas,
        9 | 17 | 35 | 53 | 85 | 117 => ChemicalFamily::Halogen,
        57..=71 => ChemicalFamily::Lanthanide,
        89..=103 => ChemicalFamily::Actinide,
        5 | 14 | 32 | 33 | 51 | 52 | 84 => ChemicalFamily::Metalloid,
        6 | 7 | 8 | 15 | 16 | 34 => ChemicalFamily::Nonmetal,
        21..=30 | 39..=48 | 72..=80 | 104..=112 => ChemicalFamily::Transition,
        _ => ChemicalFamily::PostTransition,
    }
}

pub fn period_of_z(z: u16) -> u8 {
    match z {
        0 => 1,
        1..=2 => 1,
        3..=10 => 2,
        11..=18 => 3,
        19..=36 => 4,
        37..=54 => 5,
        55..=86 => 6,
        _ => 7,
    }
}

fn noble_core_electrons(z: u16) -> f64 {
    if z > 86 {
        86.0
    } else if z > 54 {
        54.0
    } else if z > 36 {
        36.0
    } else if z > 18 {
        18.0
    } else if z > 10 {
        10.0
    } else if z > 2 {
        2.0
    } else {
        0.0
    }
}

fn valence_proxy(family: ChemicalFamily) -> f64 {
    match family {
        ChemicalFamily::Alkali => 1.0,
        ChemicalFamily::AlkalineEarth => 2.0,
        ChemicalFamily::Transition => 2.5,
        ChemicalFamily::PostTransition => 3.0,
        ChemicalFamily::Metalloid => 4.0,
        ChemicalFamily::Nonmetal => 2.0,
        ChemicalFamily::Halogen => 1.0,
        ChemicalFamily::NobleGas => 0.5,
        ChemicalFamily::Lanthanide => 3.0,
        ChemicalFamily::Actinide => 3.0,
    }
}

fn radius_family_factor(family: ChemicalFamily) -> f64 {
    match family {
        ChemicalFamily::Alkali => 1.85,
        ChemicalFamily::AlkalineEarth => 1.55,
        ChemicalFamily::Transition => 1.20,
        ChemicalFamily::PostTransition => 1.30,
        ChemicalFamily::Metalloid => 1.10,
        ChemicalFamily::Nonmetal => 0.92,
        ChemicalFamily::Halogen => 0.86,
        ChemicalFamily::NobleGas => 1.02,
        ChemicalFamily::Lanthanide => 1.42,
        ChemicalFamily::Actinide => 1.46,
    }
}

fn cohesive_multiplier(family: ChemicalFamily) -> f64 {
    match family {
        ChemicalFamily::Alkali => 0.35,
        ChemicalFamily::AlkalineEarth => 0.45,
        ChemicalFamily::Transition => 0.80,
        ChemicalFamily::PostTransition => 0.55,
        ChemicalFamily::Metalloid => 0.65,
        ChemicalFamily::Nonmetal => 0.50,
        ChemicalFamily::Halogen => 0.30,
        ChemicalFamily::NobleGas => 0.08,
        ChemicalFamily::Lanthanide => 0.75,
        ChemicalFamily::Actinide => 0.75,
    }
}

fn latent_fusion_fraction(family: ChemicalFamily) -> f64 {
    match family {
        ChemicalFamily::Alkali => 0.020,
        ChemicalFamily::AlkalineEarth => 0.025,
        ChemicalFamily::Transition => 0.035,
        ChemicalFamily::PostTransition => 0.030,
        ChemicalFamily::Metalloid => 0.030,
        ChemicalFamily::Nonmetal => 0.028,
        ChemicalFamily::Halogen => 0.015,
        ChemicalFamily::NobleGas => 0.005,
        ChemicalFamily::Lanthanide => 0.032,
        ChemicalFamily::Actinide => 0.035,
    }
}

fn latent_vapor_fraction(family: ChemicalFamily) -> f64 {
    match family {
        ChemicalFamily::Alkali => 0.55,
        ChemicalFamily::AlkalineEarth => 0.65,
        ChemicalFamily::Transition => 0.82,
        ChemicalFamily::PostTransition => 0.70,
        ChemicalFamily::Metalloid => 0.60,
        ChemicalFamily::Nonmetal => 0.45,
        ChemicalFamily::Halogen => 0.35,
        ChemicalFamily::NobleGas => 0.12,
        ChemicalFamily::Lanthanide => 0.78,
        ChemicalFamily::Actinide => 0.78,
    }
}

fn molecularity_factor(z: u16, family: ChemicalFamily) -> f64 {
    if matches!(z, 1 | 7 | 8 | 9 | 17) {
        2.0
    } else if family == ChemicalFamily::NobleGas {
        1.0
    } else {
        1.0
    }
}

fn packing_fraction(family: ChemicalFamily, period: u8) -> f64 {
    let p = period as f64;
    match family {
        ChemicalFamily::NobleGas => (0.20 + 0.045 * p).clamp(0.20, 0.48),
        ChemicalFamily::Halogen => {
            if period <= 3 {
                0.33
            } else {
                0.45
            }
        }
        ChemicalFamily::Nonmetal => {
            if period <= 2 {
                0.30
            } else {
                0.50
            }
        }
        ChemicalFamily::Metalloid => 0.56,
        ChemicalFamily::PostTransition => 0.64,
        ChemicalFamily::Transition => 0.72,
        ChemicalFamily::Lanthanide | ChemicalFamily::Actinide => 0.71,
        ChemicalFamily::Alkali => 0.58,
        ChemicalFamily::AlkalineEarth => 0.62,
    }
}

fn cohesive_scale_for_molecular_volatility(z: u16, family: ChemicalFamily, period: u8) -> f64 {
    if z == 1 {
        return 0.015;
    }
    match family {
        ChemicalFamily::NobleGas => (0.004 * (period as f64).powf(1.8)).clamp(0.004, 0.12),
        ChemicalFamily::Halogen => {
            if z == 17 {
                0.10
            } else if period == 2 {
                0.08
            } else if period == 4 {
                0.42
            } else if period >= 5 {
                0.70
            } else {
                0.16
            }
        }
        ChemicalFamily::Nonmetal => match z {
            7 => 0.025,
            8 => 0.030,
            9 => 0.032,
            _ => {
                if period <= 2 {
                    0.55
                } else if period == 3 {
                    0.80
                } else {
                    0.90
                }
            }
        },
        ChemicalFamily::Alkali => match period {
            2 => 1.60,
            3 => 2.80,
            4 => 5.20,
            5 => 5.80,
            6 => 6.20,
            _ => 6.40,
        },
        ChemicalFamily::AlkalineEarth => match period {
            2 => 1.20,
            3 => 1.50,
            4 => 2.20,
            5 => 2.60,
            6 => 3.00,
            _ => 3.20,
        },
        _ => 1.0,
    }
}

fn ambient_phase_residual(
    family: ChemicalFamily,
    period: u8,
    baseline: MatterState,
    valence_electrons_hint: Option<u16>,
) -> MatterState {
    // Residual phase lane (toggleable):
    // - heavy-halogen condensation (dispersion-dominated molecular cohesion)
    // - heavy-alkali lattice locking at ambient pressure
    // - relativistic closed-shell transition softness (period-6 d-band closure)
    if family == ChemicalFamily::Halogen {
        if period == 4 {
            return MatterState::Liquid;
        }
        if period >= 5 {
            return MatterState::Solid;
        }
    }
    if family == ChemicalFamily::Alkali && period >= 5 {
        return MatterState::Solid;
    }
    let closed_shell_transition_softening = family == ChemicalFamily::Transition
        && period >= 6
        && valence_electrons_hint.unwrap_or(0) >= 12;
    if closed_shell_transition_softening {
        return MatterState::Liquid;
    }
    baseline
}

fn phase_override_enabled() -> bool {
    static PHASE_OVERRIDE_ENABLED: OnceLock<bool> = OnceLock::new();
    *PHASE_OVERRIDE_ENABLED.get_or_init(|| match std::env::var("GUTOE_CHEM_PHASE_OVERRIDE") {
        Ok(v) => match v.to_ascii_lowercase().as_str() {
            "0" | "false" | "no" | "off" => false,
            _ => true,
        },
        Err(_) => true,
    })
}

fn coupled_radius_upper_scale(family: ChemicalFamily, period: u8) -> f64 {
    if period >= 4 {
        match family {
            ChemicalFamily::Halogen => 2.60,
            ChemicalFamily::Nonmetal => 2.20,
            ChemicalFamily::Metalloid => 2.40,
            ChemicalFamily::PostTransition => 2.10,
            ChemicalFamily::Lanthanide => 1.15,
            ChemicalFamily::Actinide => 1.20,
            ChemicalFamily::Transition if period >= 6 => 1.80,
            _ => 1.35,
        }
    } else {
        1.35
    }
}

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

pub fn chemical_thermo_calibration_from_env() -> ChemicalThermoCalibration {
    let d = ChemicalThermoCalibration::default();
    ChemicalThermoCalibration {
        pack_p_void_coef: env_f64("GUTOE_CHEM_CAL_PACK_P_VOID", d.pack_p_void_coef),
        pack_d_gain_coef: env_f64("GUTOE_CHEM_CAL_PACK_D_GAIN", d.pack_d_gain_coef),
        pack_f_gain_coef: env_f64("GUTOE_CHEM_CAL_PACK_F_GAIN", d.pack_f_gain_coef),
        pack_open_d_mult: env_f64("GUTOE_CHEM_CAL_PACK_OPEN_D_MULT", d.pack_open_d_mult),
        pack_closed_d_mult: env_f64(
            "GUTOE_CHEM_CAL_PACK_CLOSED_D_MULT",
            d.pack_closed_d_mult,
        ),
        pack_f_core_mult: env_f64("GUTOE_CHEM_CAL_PACK_F_CORE_MULT", d.pack_f_core_mult),
        radius_p_gain: env_f64("GUTOE_CHEM_CAL_RADIUS_P_GAIN", d.radius_p_gain),
        radius_closed_d_mult: env_f64(
            "GUTOE_CHEM_CAL_RADIUS_CLOSED_D_MULT",
            d.radius_closed_d_mult,
        ),
        radius_open_d_mult: env_f64("GUTOE_CHEM_CAL_RADIUS_OPEN_D_MULT", d.radius_open_d_mult),
        radius_f_core_mult: env_f64("GUTOE_CHEM_CAL_RADIUS_F_CORE_MULT", d.radius_f_core_mult),
        radius_actinide_mult: env_f64(
            "GUTOE_CHEM_CAL_RADIUS_ACTINIDE_MULT",
            d.radius_actinide_mult,
        ),
        radius_lower_actinide: env_f64(
            "GUTOE_CHEM_CAL_RADIUS_LOWER_ACTINIDE",
            d.radius_lower_actinide,
        ),
        radius_lower_transition_fcore: env_f64(
            "GUTOE_CHEM_CAL_RADIUS_LOWER_TRANSITION_FCORE",
            d.radius_lower_transition_fcore,
        ),
    }
}

fn current_calibration() -> &'static ChemicalThermoCalibration {
    static CAL: OnceLock<ChemicalThermoCalibration> = OnceLock::new();
    CAL.get_or_init(chemical_thermo_calibration_from_env)
}

fn orbital_packing_hints(scf: &AtomicScfPrediction, period: u8) -> OrbitalPackingHints {
    let n_max = scf
        .orbitals
        .iter()
        .filter(|o| o.occupation > 0)
        .map(|o| o.n)
        .max()
        .unwrap_or(1);
    let e_homo = scf.homo_energy_ev;
    let e_scale = 8.0_f64;

    let mut p_w = 0.0;
    let mut d_w = 0.0;
    let mut f_w = 0.0;
    let mut total_w = 0.0;

    let mut d_near_occ = 0.0;
    let mut d_near_cap = 0.0;
    let mut has_f_core = false;

    for o in scf.orbitals.iter().filter(|o| o.occupation > 0) {
        let occ = o.occupation as f64;
        let e_delta = (o.energy_ev - e_homo).min(0.0);
        let w = occ * (e_delta / e_scale).exp();
        total_w += w;
        match o.l {
            0 => {}
            1 => p_w += w,
            2 => d_w += w,
            _ => f_w += w,
        }

        if o.l == 2 && o.n + 1 >= n_max {
            d_near_occ += occ;
            d_near_cap += 10.0;
        }
        if o.l >= 3 && period >= 6 && o.n + 2 >= n_max {
            has_f_core = true;
        }
    }

    let norm = total_w.max(1.0e-9);
    let p_frac = p_w / norm;
    let d_frac = d_w / norm;
    let f_frac = f_w / norm;

    let d_fill = if d_near_cap > 0.0 {
        (d_near_occ / d_near_cap).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let open_d_shell = d_fill > 0.15 && d_fill < 0.90;
    let closed_d_shell = d_fill >= 0.95 && scf.valence_electrons >= 11;

    OrbitalPackingHints {
        p_frac,
        d_frac,
        f_frac,
        d_fill,
        valence_electrons: scf.valence_electrons,
        open_d_shell,
        closed_d_shell,
        has_f_core,
    }
}

fn orbital_packing_hints_from_prefetch(p: &CoupledThermoPrefetch) -> OrbitalPackingHints {
    OrbitalPackingHints {
        p_frac: p.p_frac,
        d_frac: p.d_frac,
        f_frac: p.f_frac,
        d_fill: p.d_fill,
        valence_electrons: p.valence_electrons,
        open_d_shell: p.open_d_shell,
        closed_d_shell: p.closed_d_shell,
        has_f_core: p.has_f_core,
    }
}

fn crystal_packing_multiplier(
    family: ChemicalFamily,
    period: u8,
    hints: OrbitalPackingHints,
) -> f64 {
    let period_rel = ((period as f64 - 3.0) / 4.0).clamp(0.0, 1.0);
    let d_half_fill_peak = (1.0 - 2.0 * (hints.d_fill - 0.5).abs()).clamp(0.0, 1.0);
    let d_cluster_4d = if family == ChemicalFamily::Transition && period == 5 {
        (1.0 - ((hints.d_fill - 0.65).abs() / 0.35)).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let d_band_cluster = if family == ChemicalFamily::Transition && period == 5 {
        let valence_mid = (1.0 - ((hints.valence_electrons as f64 - 7.5).abs() / 3.5)).clamp(0.0, 1.0);
        (d_cluster_4d * valence_mid).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let mut m = 1.0;
    if family == ChemicalFamily::Transition {
        m *= 1.0 + 0.24 * period_rel * hints.d_frac * (0.4 + 0.6 * d_half_fill_peak);
        m *= 1.0 + 0.20 * hints.d_frac * d_cluster_4d;
        m *= 1.0 + 0.35 * hints.d_frac * d_band_cluster;
    }
    if family == ChemicalFamily::PostTransition && period >= 5 {
        // Heavy post-transition metals trend toward denser metallic packing.
        m *= 1.0 + 0.18 * period_rel * (0.55 * hints.p_frac + 0.45 * hints.f_frac);
    }
    if hints.closed_d_shell {
        m *= 0.92;
    }
    m
}

fn allotropy_porosity_multiplier(
    family: ChemicalFamily,
    period: u8,
    hints: OrbitalPackingHints,
) -> f64 {
    if !(family == ChemicalFamily::Nonmetal || family == ChemicalFamily::Metalloid) {
        return 1.0;
    }
    let v = hints.valence_electrons as f64;
    let sp2_peak = (1.0 - ((v - 4.0).abs() / 3.0)).clamp(0.0, 1.0);
    let ring_peak = (1.0 - ((v - 6.0).abs() / 3.0)).clamp(0.0, 1.0);
    let topology_peak = (0.6 * sp2_peak + 0.4 * ring_peak).clamp(0.0, 1.0);
    let p_directionality = (hints.p_frac - 0.20 * hints.d_frac).max(0.0);
    // Open-network allotropes (graphitic layers, chains, rings) reduce packing efficiency.
    let period_bonus = if period <= 3 { 1.0 } else { 0.7 };
    1.0 - 0.45 * period_bonus * p_directionality * topology_peak
}

fn crystal_radius_multiplier(
    family: ChemicalFamily,
    period: u8,
    hints: OrbitalPackingHints,
) -> f64 {
    let period_rel = ((period as f64 - 3.0) / 4.0).clamp(0.0, 1.0);
    let d_half_fill_peak = (1.0 - 2.0 * (hints.d_fill - 0.5).abs()).clamp(0.0, 1.0);
    let d_cluster_4d = if family == ChemicalFamily::Transition && period == 5 {
        (1.0 - ((hints.d_fill - 0.65).abs() / 0.35)).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let d_band_cluster = if family == ChemicalFamily::Transition && period == 5 {
        let valence_mid = (1.0 - ((hints.valence_electrons as f64 - 7.5).abs() / 3.5)).clamp(0.0, 1.0);
        (d_cluster_4d * valence_mid).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let mut m = 1.0;
    if family == ChemicalFamily::Transition {
        m *= 1.0 - 0.14 * period_rel * hints.d_frac * (0.5 + 0.5 * d_half_fill_peak);
        m *= 1.0 - 0.16 * hints.d_frac * d_cluster_4d;
        m *= 1.0 - 0.28 * hints.d_frac * d_band_cluster;
    }
    if family == ChemicalFamily::PostTransition && period >= 5 {
        m *= 1.0 - 0.10 * period_rel * (0.55 * hints.p_frac + 0.45 * hints.f_frac);
    }
    m
}

fn allotropy_radius_multiplier(
    family: ChemicalFamily,
    period: u8,
    hints: OrbitalPackingHints,
) -> f64 {
    if !(family == ChemicalFamily::Nonmetal || family == ChemicalFamily::Metalloid) {
        return 1.0;
    }
    let v = hints.valence_electrons as f64;
    let sp2_peak = (1.0 - ((v - 4.0).abs() / 3.0)).clamp(0.0, 1.0);
    let ring_peak = (1.0 - ((v - 6.0).abs() / 3.0)).clamp(0.0, 1.0);
    let topology_peak = (0.6 * sp2_peak + 0.4 * ring_peak).clamp(0.0, 1.0);
    let p_directionality = (hints.p_frac - 0.20 * hints.d_frac).max(0.0);
    let period_bonus = if period <= 3 { 1.0 } else { 0.7 };
    1.0 + 0.22 * period_bonus * p_directionality * topology_peak
}

fn packing_fraction_from_hints(
    family: ChemicalFamily,
    period: u8,
    hints: OrbitalPackingHints,
    cal: &ChemicalThermoCalibration,
) -> f64 {
    let base = packing_fraction(family, period);
    let mut pack = base
        * (1.0 - cal.pack_p_void_coef * hints.p_frac
            + cal.pack_d_gain_coef * hints.d_frac
            + cal.pack_f_gain_coef * hints.f_frac);
    pack *= crystal_packing_multiplier(family, period, hints);
    pack *= allotropy_porosity_multiplier(family, period, hints);
    let period_rel = ((period as f64 - 3.0) / 4.0).clamp(0.0, 1.0);
    // Half-filled d shells compact strongly (band-filling cohesion), especially in 4d/5d rows.
    let d_half_fill_peak = (1.0 - 2.0 * (hints.d_fill - 0.5).abs()).clamp(0.0, 1.0);
    pack *= 1.0 + 0.25 * d_half_fill_peak * hints.d_frac * period_rel;
    // Directional p-bond networks create open structures (graphitic/chain/ring porosity).
    let v = hints.valence_electrons as f64;
    let sp_network_peak = (1.0 - ((v - 4.0).abs() / 4.0)).clamp(0.0, 1.0);
    let ring_network_peak = (1.0 - ((v - 6.0).abs() / 4.0)).clamp(0.0, 1.0);
    let topology_peak = (0.7 * sp_network_peak + 0.3 * ring_network_peak).clamp(0.0, 1.0);
    let p_directionality = (hints.p_frac - 0.25 * hints.d_frac).max(0.0);
    pack *= 1.0 - (1.0 / 3.0) * p_directionality * topology_peak;
    if hints.open_d_shell && period >= 5 {
        pack *= cal.pack_open_d_mult;
    }
    if hints.closed_d_shell {
        pack *= cal.pack_closed_d_mult;
    }
    if hints.has_f_core && period >= 6 {
        pack *= cal.pack_f_core_mult;
    }
    pack.clamp(0.18, 0.86)
}

/// Clausius-Clapeyron transduction using the element's boiling point as anchor
/// (P_sat(T_b) = 1 atm).
pub fn vapor_pressure_clausius_pa(
    latent_vaporization_kj_mol: f64,
    boiling_temperature_k: f64,
    t_k: f64,
) -> f64 {
    let t = t_k.max(1.0);
    let t_b = boiling_temperature_k.max(2.0);
    let delta_h = (latent_vaporization_kj_mol * 1000.0).max(1.0);
    let exponent = (-delta_h / R_GAS_J_MOL_K * (1.0 / t - 1.0 / t_b)).clamp(-80.0, 80.0);
    (P_REF_PA * exponent.exp()).clamp(1.0e-9, 1.0e12)
}

/// Reduced Gibbs offsets for solid/liquid/gas at (T, P). Solid is used as the
/// zero reference and liquid/gas are computed from latent-heat transductions.
pub fn phase_gibbs_offsets_kj_mol(
    latent_fusion_kj_mol: f64,
    latent_vaporization_kj_mol: f64,
    melting_temperature_k: f64,
    boiling_temperature_k: f64,
    t_k: f64,
    pressure_pa: f64,
) -> PhaseGibbsKjMol {
    let t = t_k.max(1.0);
    let t_m = melting_temperature_k.max(1.0);
    let t_b = boiling_temperature_k.max(t_m + 1.0e-6);
    let delta_h_f = latent_fusion_kj_mol.max(1.0e-6);
    let delta_h_v = latent_vaporization_kj_mol.max(1.0e-6);
    let p = pressure_pa.max(1.0e-12);

    let g_solid = 0.0;
    let g_liquid = delta_h_f * (1.0 - t / t_m);
    let g_gas_ref = g_liquid + delta_h_v * (1.0 - t / t_b);
    let pressure_term = (R_GAS_J_MOL_K * t / 1000.0) * (p / P_REF_PA).ln();
    let g_gas = g_gas_ref + pressure_term;

    PhaseGibbsKjMol {
        solid: g_solid,
        liquid: g_liquid,
        gas: g_gas,
    }
}

/// Select phase by minimum reduced Gibbs free energy at (T, P).
pub fn phase_from_gibbs(
    latent_fusion_kj_mol: f64,
    latent_vaporization_kj_mol: f64,
    melting_temperature_k: f64,
    boiling_temperature_k: f64,
    t_k: f64,
    pressure_pa: f64,
) -> MatterState {
    let g = phase_gibbs_offsets_kj_mol(
        latent_fusion_kj_mol,
        latent_vaporization_kj_mol,
        melting_temperature_k,
        boiling_temperature_k,
        t_k,
        pressure_pa,
    );
    let mut state = MatterState::Solid;
    let mut best = g.solid;
    if g.liquid < best {
        best = g.liquid;
        state = MatterState::Liquid;
    }
    if g.gas < best {
        state = MatterState::Gas;
    }
    state
}

fn assemble_element_thermo(
    z: u16,
    a: u16,
    family: ChemicalFamily,
    period: u8,
    radius_pm: f64,
    cohesive_energy_ev_per_atom: f64,
    valence_electrons_hint: Option<u16>,
    packing_fraction_hint: Option<f64>,
) -> ElementThermoPrediction {
    // Spherical atom proxy with family-aware packing correction.
    let pack = packing_fraction_hint
        .unwrap_or_else(|| packing_fraction(family, period))
        .max(0.05);
    let r_cm = radius_pm * 1.0e-10;
    let atom_vol_cm3 = (4.0 / 3.0) * PI * r_cm.powi(3) / pack;
    let molar_volume_cm3_mol = atom_vol_cm3 * AVOGADRO;

    let molar_mass_g_mol = a as f64;
    let condensed_density_g_cm3 = (molar_mass_g_mol / molar_volume_cm3_mol).clamp(0.0005, 40.0);

    let latent_fusion_kj_mol =
        cohesive_energy_ev_per_atom * EV_TO_KJ_MOL * latent_fusion_fraction(family);
    let latent_vaporization_kj_mol =
        cohesive_energy_ev_per_atom * EV_TO_KJ_MOL * latent_vapor_fraction(family);

    let melting_temperature_k =
        (latent_fusion_kj_mol * 1000.0 / ENTROPY_FUSION_J_MOL_K).clamp(2.0, 8000.0);
    let boiling_temperature_k =
        (latent_vaporization_kj_mol * 1000.0 / ENTROPY_VAPORIZATION_J_MOL_K).clamp(4.0, 12000.0);

    let debye_temperature_k =
        (120.0 * cohesive_energy_ev_per_atom.sqrt() * (condensed_density_g_cm3 / 5.0).powf(0.25))
            .clamp(20.0, 2200.0);

    // Debye-like saturation toward Dulong-Petit near room temperature.
    let cp_solid_j_mol_k =
        (3.0 * R_GAS_J_MOL_K * (1.0 - (-298.15 / (0.35 * debye_temperature_k)).exp()))
            .clamp(1.5 * R_GAS_J_MOL_K, 3.0 * R_GAS_J_MOL_K);
    let cp_liquid_j_mol_k = 3.5 * R_GAS_J_MOL_K;
    let cp_gas_j_mol_k = match family {
        ChemicalFamily::Nonmetal | ChemicalFamily::Halogen => 3.5 * R_GAS_J_MOL_K,
        _ => 2.5 * R_GAS_J_MOL_K,
    };

    let bulk_modulus_gpa =
        (20.0 * cohesive_energy_ev_per_atom * condensed_density_g_cm3.powf(0.7)).clamp(0.1, 500.0);
    let thermal_expansion_1_per_k =
        (2.2e-5 * (300.0 / debye_temperature_k).powf(0.7) * (30.0 / bulk_modulus_gpa).powf(0.3))
            .clamp(1.0e-6, 2.5e-4);

    let vapor_pressure_pa_298k =
        vapor_pressure_clausius_pa(latent_vaporization_kj_mol, boiling_temperature_k, 298.15);
    let ambient_state_raw = phase_from_gibbs(
        latent_fusion_kj_mol,
        latent_vaporization_kj_mol,
        melting_temperature_k,
        boiling_temperature_k,
        298.15,
        P_REF_PA,
    );
    let ambient_state_298k = if phase_override_enabled() {
        ambient_phase_residual(family, period, ambient_state_raw, valence_electrons_hint)
    } else {
        ambient_state_raw
    };
    let density_g_cm3 = if ambient_state_298k == MatterState::Gas {
        let molar_volume_298_l = (R_GAS_J_MOL_K * 298.15 / P_REF_PA) * 1000.0;
        let molecular_molar_mass = molar_mass_g_mol * molecularity_factor(z, family);
        (molecular_molar_mass / (molar_volume_298_l * 1000.0)).clamp(1.0e-6, 10.0)
    } else {
        condensed_density_g_cm3
    };

    ElementThermoPrediction {
        z,
        a,
        family,
        period,
        molar_mass_g_mol,
        atomic_radius_pm: radius_pm,
        molar_volume_cm3_mol,
        density_g_cm3,
        cohesive_energy_ev_per_atom,
        debye_temperature_k,
        latent_fusion_kj_mol,
        latent_vaporization_kj_mol,
        melting_temperature_k,
        boiling_temperature_k,
        vapor_pressure_pa_298k,
        cp_solid_j_mol_k,
        cp_liquid_j_mol_k,
        cp_gas_j_mol_k,
        bulk_modulus_gpa,
        thermal_expansion_1_per_k,
        ambient_state_298k,
    }
}

pub fn predict_element_thermo(z: u16, a: u16) -> ElementThermoPrediction {
    let family = family_of_z(z);
    let period = period_of_z(z);
    let period_f = period as f64;
    let valence = valence_proxy(family);

    let z_eff = ((z as f64) - noble_core_electrons(z) + 0.5 * valence).max(1.0);
    let radius_pm = (BOHR_RADIUS_PM * (period_f * period_f) / z_eff * radius_family_factor(family))
        .clamp(30.0, 320.0);

    let cohesive_energy_ev_per_atom = ((13.605_693 * valence.powi(2) / period_f.powi(2)
        * cohesive_multiplier(family))
        * cohesive_scale_for_molecular_volatility(z, family, period))
    .clamp(0.005, 12.0);

    assemble_element_thermo(
        z,
        a,
        family,
        period,
        radius_pm,
        cohesive_energy_ev_per_atom,
        None,
        None,
    )
}

pub fn predict_element_thermo_coupled_with_diagnostics_calibrated(
    z: u16,
    a: u16,
    cal: ChemicalThermoCalibration,
) -> (ElementThermoPrediction, CoupledThermoDiagnostics) {
    let prefetch = prefetch_element_thermo_coupled(z, a);
    predict_element_thermo_coupled_from_prefetch_calibrated(&prefetch, cal)
}

pub fn prefetch_element_thermo_coupled(z: u16, a: u16) -> CoupledThermoPrefetch {
    let family = family_of_z(z);
    let period = period_of_z(z);
    let period_f = period as f64;
    let base = predict_element_thermo(z, a);
    let scf = predict_atomic_scf(z, a);
    let hints = orbital_packing_hints(&scf, period);
    let period_weight = (period_f / 8.0).clamp(0.1, 0.9);
    let scf_radius_weight = (0.25 + 0.20 * period_weight).clamp(0.25, 0.45);

    CoupledThermoPrefetch {
        z,
        a,
        family,
        period,
        base,
        scf,
        p_frac: hints.p_frac,
        d_frac: hints.d_frac,
        f_frac: hints.f_frac,
        d_fill: hints.d_fill,
        valence_electrons: hints.valence_electrons,
        open_d_shell: hints.open_d_shell,
        closed_d_shell: hints.closed_d_shell,
        has_f_core: hints.has_f_core,
        scf_radius_weight,
    }
}

pub fn predict_element_thermo_coupled_from_prefetch_calibrated(
    prefetch: &CoupledThermoPrefetch,
    cal: ChemicalThermoCalibration,
) -> (ElementThermoPrediction, CoupledThermoDiagnostics) {
    let family = prefetch.family;
    let period = prefetch.period;
    let z = prefetch.z;
    let a = prefetch.a;
    let base = prefetch.base;
    let scf = &prefetch.scf;
    let hints = orbital_packing_hints_from_prefetch(prefetch);

    let mut raw_radius_pm = (1.0 - prefetch.scf_radius_weight) * base.atomic_radius_pm
        + prefetch.scf_radius_weight * scf.atomic_radius_pm;
    let mut radius_factor = 1.0;
    radius_factor *= crystal_radius_multiplier(family, period, hints);
    radius_factor *= allotropy_radius_multiplier(family, period, hints);
    let period_rel = ((period as f64 - 3.0) / 4.0).clamp(0.0, 1.0);
    let d_half_fill_peak = (1.0 - 2.0 * (hints.d_fill - 0.5).abs()).clamp(0.0, 1.0);
    // Relativistic contraction scales ~ (Z/137)^2 and strengthens dense d/f blocks.
    let z_rel = ((z as f64) / 137.0).powi(2).clamp(0.0, 1.0);
    let df_weight = (hints.d_frac + hints.f_frac).clamp(0.0, 1.0);
    radius_factor *= 1.0 - 0.125 * z_rel * df_weight * (0.5 + 0.5 * period_rel);
    // Half-filled d shells contract more strongly than closed d shells.
    radius_factor *= 1.0 - 0.10 * d_half_fill_peak * hints.d_frac * period_rel;
    radius_factor *= (1.0 + cal.radius_p_gain * hints.p_frac).clamp(1.0, 1.30);
    if hints.closed_d_shell {
        radius_factor *= cal.radius_closed_d_mult;
    }
    if hints.open_d_shell && period >= 5 {
        radius_factor *= cal.radius_open_d_mult;
    }
    if hints.has_f_core && period >= 6 {
        let fcore_mult = if family == ChemicalFamily::Lanthanide {
            6.0 / 7.0
        } else {
            cal.radius_f_core_mult
        };
        radius_factor *= fcore_mult;
    }
    if period >= 6 && family == ChemicalFamily::Actinide {
        radius_factor *= cal.radius_actinide_mult;
    }
    raw_radius_pm *= radius_factor;
    let radius_upper_scale = coupled_radius_upper_scale(family, period);
    let radius_lower_scale = if family == ChemicalFamily::Actinide {
        cal.radius_lower_actinide
    } else if family == ChemicalFamily::Transition && period == 5 {
        let d_cluster_4d = (1.0 - ((hints.d_fill - 0.65).abs() / 0.35)).clamp(0.0, 1.0);
        let valence_mid = (1.0 - ((hints.valence_electrons as f64 - 7.5).abs() / 3.5)).clamp(0.0, 1.0);
        let cluster = (d_cluster_4d * valence_mid).clamp(0.0, 1.0);
        (0.52 - 0.16 * cluster).clamp(0.34, 0.52)
    } else if hints.has_f_core && period >= 6 && family == ChemicalFamily::Transition {
        cal.radius_lower_transition_fcore
    } else {
        0.70
    };
    let coupled_radius_pm = raw_radius_pm
        .clamp(
            radius_lower_scale * base.atomic_radius_pm,
            radius_upper_scale * base.atomic_radius_pm,
        )
        .clamp(25.0, 350.0);
    let coupled_packing_fraction = packing_fraction_from_hints(family, period, hints, &cal);

    let frontier_cohesive = (0.28 * scf.ionization_energy_ev
        + 0.72 * scf.electron_affinity_ev.max(0.0)
        + 0.15 * scf.electronegativity_mulliken_ev)
        .clamp(0.03, 12.0);
    let frontier_weight = match family {
        ChemicalFamily::NobleGas => 0.04,
        ChemicalFamily::Halogen => 0.10,
        ChemicalFamily::Nonmetal => 0.12,
        ChemicalFamily::Alkali | ChemicalFamily::AlkalineEarth => 0.16,
        ChemicalFamily::Metalloid => 0.18,
        ChemicalFamily::PostTransition => 0.20,
        ChemicalFamily::Transition | ChemicalFamily::Lanthanide | ChemicalFamily::Actinide => 0.22,
    };
    let hardness_gate = (1.0 + 0.006 * scf.chemical_hardness_ev).clamp(0.92, 1.08);
    let valence_gate = (0.95 + 0.02 * (scf.valence_electrons as f64).min(8.0)).clamp(0.95, 1.08);
    let raw_cohesive = ((1.0 - frontier_weight) * base.cohesive_energy_ev_per_atom
        + frontier_weight * frontier_cohesive)
        * hardness_gate
        * valence_gate;
    let coupled_cohesive = (raw_cohesive
        * cohesive_scale_for_molecular_volatility(z, family, period))
        .clamp(0.70 * base.cohesive_energy_ev_per_atom, 1.35 * base.cohesive_energy_ev_per_atom)
        .clamp(0.005, 12.0);

    let prediction = assemble_element_thermo(
        z,
        a,
        family,
        period,
        coupled_radius_pm,
        coupled_cohesive,
        Some(scf.valence_electrons),
        Some(coupled_packing_fraction),
    );
    let diag = CoupledThermoDiagnostics {
        scf_iterations: scf.scf_iterations,
        scf_residual: scf.scf_residual,
        valence_electrons: scf.valence_electrons,
        scf_atomic_radius_pm: scf.atomic_radius_pm,
        scf_ionization_energy_ev: scf.ionization_energy_ev,
        scf_electron_affinity_ev: scf.electron_affinity_ev,
        scf_electronegativity_mulliken_ev: scf.electronegativity_mulliken_ev,
        scf_chemical_hardness_ev: scf.chemical_hardness_ev,
        coupled_radius_pm,
        cohesive_frontier_proxy_ev: frontier_cohesive,
        coupled_packing_fraction,
    };

    (prediction, diag)
}

pub fn predict_element_thermo_coupled_with_diagnostics(
    z: u16,
    a: u16,
) -> (ElementThermoPrediction, CoupledThermoDiagnostics) {
    predict_element_thermo_coupled_with_diagnostics_calibrated(z, a, *current_calibration())
}

pub fn predict_element_thermo_coupled_calibrated(
    z: u16,
    a: u16,
    cal: ChemicalThermoCalibration,
) -> ElementThermoPrediction {
    predict_element_thermo_coupled_with_diagnostics_calibrated(z, a, cal).0
}

pub fn predict_element_thermo_coupled(z: u16, a: u16) -> ElementThermoPrediction {
    predict_element_thermo_coupled_with_diagnostics(z, a).0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boiling_is_above_melting_for_all_families() {
        for z in 1..=118 {
            let p = predict_element_thermo(z, (2.5 * z as f64).round() as u16);
            assert!(p.boiling_temperature_k > p.melting_temperature_k);
        }
    }

    #[test]
    fn noble_gases_are_low_temperature_volatiles() {
        for &z in &[2_u16, 10, 18, 36, 54, 86, 118] {
            let p = predict_element_thermo(z, (2.5 * z as f64).round() as u16);
            assert!(
                p.boiling_temperature_k < 400.0,
                "noble gas Z={z} expected low boiling point, got {} K",
                p.boiling_temperature_k
            );
        }
    }

    #[test]
    fn transition_family_is_more_cohesive_than_alkali_proxy() {
        // Fe-like vs Na-like proxy comparison.
        let fe_like = predict_element_thermo(26, 56);
        let na_like = predict_element_thermo(11, 23);
        assert!(fe_like.cohesive_energy_ev_per_atom > na_like.cohesive_energy_ev_per_atom);
    }

    #[test]
    fn clapeyron_anchor_hits_one_atm_at_boiling_point() {
        for z in 1..=118 {
            let p = predict_element_thermo(z, (2.5 * z as f64).round() as u16);
            let p_sat = vapor_pressure_clausius_pa(
                p.latent_vaporization_kj_mol,
                p.boiling_temperature_k,
                p.boiling_temperature_k,
            );
            assert!(((p_sat / P_REF_PA) - 1.0).abs() < 1.0e-10);
        }
    }

    #[test]
    fn gibbs_phase_matches_threshold_rule_at_reference_pressure() {
        for z in 1..=118 {
            let p = predict_element_thermo(z, (2.5 * z as f64).round() as u16);
            let threshold_state = if 298.15 < p.melting_temperature_k {
                MatterState::Solid
            } else if 298.15 < p.boiling_temperature_k {
                MatterState::Liquid
            } else {
                MatterState::Gas
            };
            let gibbs_state = phase_from_gibbs(
                p.latent_fusion_kj_mol,
                p.latent_vaporization_kj_mol,
                p.melting_temperature_k,
                p.boiling_temperature_k,
                298.15,
                P_REF_PA,
            );
            let corrected = ambient_phase_residual(p.family, p.period, gibbs_state, None);
            assert_eq!(p.ambient_state_298k, corrected);
            if corrected == gibbs_state {
                assert_eq!(gibbs_state, threshold_state);
            }
        }
    }

    #[test]
    fn high_pressure_can_condense_noble_gas_proxy() {
        let he_like = predict_element_thermo(2, 4);
        let ambient = phase_from_gibbs(
            he_like.latent_fusion_kj_mol,
            he_like.latent_vaporization_kj_mol,
            he_like.melting_temperature_k,
            he_like.boiling_temperature_k,
            298.15,
            P_REF_PA,
        );
        assert_eq!(ambient, MatterState::Gas);

        let extreme_pressure = phase_from_gibbs(
            he_like.latent_fusion_kj_mol,
            he_like.latent_vaporization_kj_mol,
            he_like.melting_temperature_k,
            he_like.boiling_temperature_k,
            298.15,
            1.0e9,
        );
        assert!(extreme_pressure != MatterState::Gas);
    }

    #[test]
    fn coupled_lane_preserves_physical_ordering() {
        for z in 1..=140 {
            let p = predict_element_thermo_coupled(z, (2.5 * z as f64).round() as u16);
            assert!(p.boiling_temperature_k > p.melting_temperature_k);
            assert!(p.density_g_cm3.is_finite() && p.density_g_cm3 > 0.0);
            assert!(p.atomic_radius_pm.is_finite() && p.atomic_radius_pm > 0.0);
        }
    }

    #[test]
    fn coupled_lane_exposes_finite_scf_diagnostics() {
        for &z in &[1_u16, 6, 8, 11, 17, 26, 36, 54, 79, 94, 118, 140] {
            let (_p, d) =
                predict_element_thermo_coupled_with_diagnostics(z, (2.5 * z as f64).round() as u16);
            assert!(d.scf_atomic_radius_pm.is_finite());
            assert!(d.scf_ionization_energy_ev.is_finite());
            assert!(d.scf_electronegativity_mulliken_ev.is_finite());
            assert!(d.cohesive_frontier_proxy_ev.is_finite());
        }
    }
}
