use gutoe_physics::{gamow_factor, ALPHA_LEADING_ORDER};
use serde_json::Value;
use std::env;
use std::f64::consts::PI;
use std::fs;
use std::path::Path;

// ─── A. Nuclear Mass Table (AME2020, MeV) ───────────────────────────────────

const PROTON: f64 = 938.272_088_16;
const NEUTRON: f64 = 939.565_420_52;
const DEUTERON: f64 = 1875.612_928;
const TRITON: f64 = 2808.921_12;
const HE3: f64 = 2808.391_586;
const HE4: f64 = 3727.379_40;
const LI6: f64 = 5601.518_6;
const B10: f64 = 9324.436_2;
const B11: f64 = 10252.547_5;
const O16: f64 = 14895.079_8;

// ─── Constants ───────────────────────────────────────────────────────────────

const MEV_TO_J: f64 = 1.602_176_634e-13;
const KB: f64 = 1.380_649e-23;
const KEV_TO_K: f64 = 1.160_451_812e7;
const MU0: f64 = 4.0 * PI * 1.0e-7;
const NUCLEON_MASS_MEV: f64 = 938.272_088_16;

// RTSC embedded fallback (Cr toroidal_honeycomb best point)
const RTSC_R_MAJOR_M: f64 = 1.54;
const RTSC_A_MINOR_M: f64 = 0.64;
const RTSC_T_KEV: f64 = 150.0;
const RTSC_B_OPERATING_T: f64 = 26.8;
const RTSC_CONFINEMENT: f64 = 0.548;
const RTSC_P_NET_MW: f64 = 29.37;

const CALIBRATED_SIGMAV_DT_100KEV: f64 = 5.0e-22; // m^3/s, D+T reference
const CALIBRATED_SIGMAV_DHE3_100KEV: f64 = 1.0e-22; // m^3/s, D+He3 reference
const MAGNET_LOAD_COEFF: f64 = 60.0;
const PARTICLE_PRESSURE_FACTOR: f64 = 2.5;
const ELECTRIC_EFF: f64 = 0.62;
const SECONDS_PER_YEAR: f64 = 365.25 * 86400.0;

// ─── B. Reaction Catalog ─────────────────────────────────────────────────────

#[derive(Clone, Debug)]
struct CascadeReaction {
    id: &'static str,
    z1: u16,
    a1: u16,
    z2: u16,
    a2: u16,
    reactant1_mass: f64,
    reactant2_mass: f64,
    q_fusion_mev: f64,
    product_masses: Vec<f64>, // masses of antimatter products (gamma = 0.0)
    product_labels: &'static str,
    calibrated_sigmav_100kev: f64,
}

fn reaction_catalog() -> Vec<CascadeReaction> {
    vec![
        CascadeReaction {
            id: "anti-D + anti-T -> anti-He4 + anti-n",
            z1: 1,
            a1: 2,
            z2: 1,
            a2: 3,
            reactant1_mass: DEUTERON,
            reactant2_mass: TRITON,
            q_fusion_mev: DEUTERON + TRITON - HE4 - NEUTRON,
            product_masses: vec![HE4, NEUTRON],
            product_labels: "anti-He4, anti-n",
            calibrated_sigmav_100kev: CALIBRATED_SIGMAV_DT_100KEV,
        },
        CascadeReaction {
            id: "anti-D + anti-He3 -> anti-He4 + anti-p",
            z1: 1,
            a1: 2,
            z2: 2,
            a2: 3,
            reactant1_mass: DEUTERON,
            reactant2_mass: HE3,
            q_fusion_mev: DEUTERON + HE3 - HE4 - PROTON,
            product_masses: vec![HE4, PROTON],
            product_labels: "anti-He4, anti-p",
            calibrated_sigmav_100kev: CALIBRATED_SIGMAV_DHE3_100KEV,
        },
        CascadeReaction {
            id: "anti-D + anti-D -> anti-T + anti-p",
            z1: 1,
            a1: 2,
            z2: 1,
            a2: 2,
            reactant1_mass: DEUTERON,
            reactant2_mass: DEUTERON,
            q_fusion_mev: 2.0 * DEUTERON - TRITON - PROTON,
            product_masses: vec![TRITON, PROTON],
            product_labels: "anti-T, anti-p",
            calibrated_sigmav_100kev: CALIBRATED_SIGMAV_DHE3_100KEV * 0.5,
        },
        CascadeReaction {
            id: "anti-D + anti-D -> anti-He3 + anti-n",
            z1: 1,
            a1: 2,
            z2: 1,
            a2: 2,
            reactant1_mass: DEUTERON,
            reactant2_mass: DEUTERON,
            q_fusion_mev: 2.0 * DEUTERON - HE3 - NEUTRON,
            product_masses: vec![HE3, NEUTRON],
            product_labels: "anti-He3, anti-n",
            calibrated_sigmav_100kev: CALIBRATED_SIGMAV_DHE3_100KEV * 0.5,
        },
        CascadeReaction {
            id: "anti-D + anti-D -> anti-He4 + gamma",
            z1: 1,
            a1: 2,
            z2: 1,
            a2: 2,
            reactant1_mass: DEUTERON,
            reactant2_mass: DEUTERON,
            q_fusion_mev: 2.0 * DEUTERON - HE4,
            product_masses: vec![HE4, 0.0], // gamma has zero mass
            product_labels: "anti-He4, gamma",
            // EM suppression for gamma branch
            calibrated_sigmav_100kev: CALIBRATED_SIGMAV_DHE3_100KEV
                * dd_gamma_branch_suppression(2.0 * DEUTERON - HE4),
        },
        CascadeReaction {
            id: "anti-p + anti-B11 -> 3x anti-He4",
            z1: 1,
            a1: 1,
            z2: 5,
            a2: 11,
            reactant1_mass: PROTON,
            reactant2_mass: B11,
            q_fusion_mev: PROTON + B11 - 3.0 * HE4,
            product_masses: vec![HE4, HE4, HE4],
            product_labels: "3x anti-He4",
            calibrated_sigmav_100kev: CALIBRATED_SIGMAV_DHE3_100KEV * 0.1,
        },
        CascadeReaction {
            id: "anti-B10 + anti-B10 -> anti-O16 + anti-He4",
            z1: 5,
            a1: 10,
            z2: 5,
            a2: 10,
            reactant1_mass: B10,
            reactant2_mass: B10,
            q_fusion_mev: 2.0 * B10 - O16 - HE4,
            product_masses: vec![O16, HE4],
            product_labels: "anti-O16, anti-He4",
            calibrated_sigmav_100kev: CALIBRATED_SIGMAV_DHE3_100KEV * 0.001,
        },
        CascadeReaction {
            id: "anti-D + anti-Li6 -> 2x anti-He4",
            z1: 1,
            a1: 2,
            z2: 3,
            a2: 6,
            reactant1_mass: DEUTERON,
            reactant2_mass: LI6,
            q_fusion_mev: DEUTERON + LI6 - 2.0 * HE4,
            product_masses: vec![HE4, HE4],
            product_labels: "2x anti-He4",
            calibrated_sigmav_100kev: CALIBRATED_SIGMAV_DHE3_100KEV * 0.3,
        },
    ]
}

fn dd_gamma_branch_suppression(q_gamma_mev: f64) -> f64 {
    (ALPHA_LEADING_ORDER * (q_gamma_mev / NUCLEON_MASS_MEV).powi(3)).clamp(0.0, 1.0)
}

// ─── C. Per-Reaction Cascade Computation ─────────────────────────────────────

#[derive(Clone, Debug)]
struct CascadeResult {
    q_fusion: f64,
    q_from_masses: f64,
    annihilation_per_product: Vec<f64>,
    total_annihilation: f64,
    total_cascade: f64,
    multiplier: f64,
    fuel_rest_mass: f64,
    direct_annihilation: f64,
    cascade_penalty: f64,
    break_even_eta: f64,
}

fn compute_cascade(r: &CascadeReaction) -> CascadeResult {
    let q_from_masses = r.reactant1_mass + r.reactant2_mass
        - r.product_masses.iter().sum::<f64>();
    let q_fusion = r.q_fusion_mev;

    // Each antimatter product annihilates with its matter counterpart: E = 2 * mass
    let annihilation_per_product: Vec<f64> =
        r.product_masses.iter().map(|m| 2.0 * m).collect();
    let total_annihilation: f64 = annihilation_per_product.iter().sum();
    let total_cascade = q_fusion + total_annihilation;

    let multiplier = total_cascade / q_fusion;

    let fuel_rest_mass = r.reactant1_mass + r.reactant2_mass;
    let direct_annihilation = 2.0 * fuel_rest_mass;

    // Mass conservation identity: cascade = direct - Q_fusion
    // Because: products mass = reactants mass - Q  =>  2*products = 2*(reactants - Q) = direct - 2Q
    // total_cascade = Q + 2*products = Q + direct - 2Q = direct - Q
    let cascade_penalty = direct_annihilation - total_cascade;

    let break_even_eta = fuel_rest_mass / total_cascade;

    CascadeResult {
        q_fusion,
        q_from_masses,
        annihilation_per_product,
        total_annihilation,
        total_cascade,
        multiplier,
        fuel_rest_mass,
        direct_annihilation,
        cascade_penalty,
        break_even_eta,
    }
}

// ─── D. Reactor-Level Integration ────────────────────────────────────────────

fn reduced_mass_mev(a1: f64, a2: f64) -> f64 {
    (a1 * a2 / (a1 + a2)) * 931.494
}

fn sigma_v_general(z1: u16, z2: u16, a1: u16, a2: u16, t_kev: f64, ref_sv: f64) -> f64 {
    if t_kev <= 0.0 {
        return 0.0;
    }
    let alpha_eff = ALPHA_LEADING_ORDER * (z1 as f64) * (z2 as f64);
    let m_reduced = reduced_mass_mev(a1 as f64, a2 as f64);
    let e_cm_mev = t_kev / 1000.0;
    let g = gamow_factor(alpha_eff, m_reduced, e_cm_mev).unwrap_or(0.0);

    // Reference Gamow at 100 keV for this pair
    let g_ref = gamow_factor(alpha_eff, m_reduced, 0.1).unwrap_or(1.0e-30).max(1.0e-30);
    let thermal = (t_kev / 100.0).max(0.05).sqrt();
    let v = ref_sv * (g / g_ref) * thermal;
    v.clamp(1.0e-30, 1.0e-18)
}

struct RtscConfig {
    r_major_m: f64,
    a_minor_m: f64,
    t_kev: f64,
    b_operating_t: f64,
    confinement: f64,
    p_net_regular_w: f64,
    volume_m3: f64,
    n_fuel_m3: f64,
    p_recirc_w: f64,
}

fn load_rtsc_config() -> RtscConfig {
    // Try loading from JSON artifacts, fall back to embedded defaults
    let json_path = env::var("GUTOE_CASCADE_RTSC_JSON")
        .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_honeycomb_fusion/rtsc_honeycomb_fusion.json".to_string());

    if let Some(cfg) = load_rtsc_from_json(Path::new(&json_path)) {
        return cfg;
    }

    // Embedded fallback from RTSC best point
    let volume_m3 = 2.0 * PI * PI * RTSC_R_MAJOR_M * RTSC_A_MINOR_M.powi(2);
    let t_k = RTSC_T_KEV * KEV_TO_K;
    let p_magnetic = RTSC_B_OPERATING_T * RTSC_B_OPERATING_T / (2.0 * MU0);
    let beta_target = 0.045 * 1.15; // toroidal_honeycomb beta factor
    let p_plasma = beta_target * p_magnetic;
    let n_fuel_m3 = p_plasma / (PARTICLE_PRESSURE_FACTOR * KB * t_k);

    // Surface of torus
    let surface_m2 = 4.0 * PI * PI * RTSC_R_MAJOR_M * RTSC_A_MINOR_M;
    let p_magnet = surface_m2 * RTSC_B_OPERATING_T.powi(2) * MAGNET_LOAD_COEFF * 1.06;
    let p_fusion_w = RTSC_P_NET_MW * 1.0e6 / (ELECTRIC_EFF - 0.16); // back-compute
    let p_heating = p_fusion_w * (0.12 + 0.28 * (1.0 - RTSC_CONFINEMENT));
    let p_aux = p_fusion_w * 0.06;
    let p_recirc_w = p_magnet + p_heating + p_aux;

    RtscConfig {
        r_major_m: RTSC_R_MAJOR_M,
        a_minor_m: RTSC_A_MINOR_M,
        t_kev: RTSC_T_KEV,
        b_operating_t: RTSC_B_OPERATING_T,
        confinement: RTSC_CONFINEMENT,
        p_net_regular_w: RTSC_P_NET_MW * 1.0e6,
        volume_m3,
        n_fuel_m3,
        p_recirc_w,
    }
}

fn load_rtsc_from_json(json_path: &Path) -> Option<RtscConfig> {
    let txt = fs::read_to_string(json_path).ok()?;
    let v: Value = serde_json::from_str(&txt).ok()?;
    let b = v.get("best_overall")?.as_object()?;

    let r_major_m = b.get("r_major_m").and_then(Value::as_f64)?;
    let a_minor_m = b.get("a_minor_m").and_then(Value::as_f64)?;
    let t_kev = b.get("t_kev").and_then(Value::as_f64)?;
    let p_net_w = b.get("p_net_w").and_then(Value::as_f64)?;
    let b_t = b.get("b_operating_t").and_then(Value::as_f64).unwrap_or(RTSC_B_OPERATING_T);
    let confinement = b.get("confinement").and_then(Value::as_f64).unwrap_or(RTSC_CONFINEMENT);

    let volume_m3 = 2.0 * PI * PI * r_major_m * a_minor_m.powi(2);
    let t_k = t_kev * KEV_TO_K;
    let p_magnetic = b_t * b_t / (2.0 * MU0);
    let beta_target = 0.045 * 1.15;
    let p_plasma = beta_target * p_magnetic;
    let n_fuel_m3 = p_plasma / (PARTICLE_PRESSURE_FACTOR * KB * t_k);

    let surface_m2 = 4.0 * PI * PI * r_major_m * a_minor_m;
    let p_magnet = surface_m2 * b_t.powi(2) * MAGNET_LOAD_COEFF * 1.06;
    let p_fusion_w = p_net_w / (ELECTRIC_EFF - 0.16);
    let p_heating = p_fusion_w * (0.12 + 0.28 * (1.0 - confinement));
    let p_aux = p_fusion_w * 0.06;
    let p_recirc_w = p_magnet + p_heating + p_aux;

    Some(RtscConfig {
        r_major_m,
        a_minor_m,
        t_kev,
        b_operating_t: b_t,
        confinement,
        p_net_regular_w: p_net_w,
        volume_m3,
        n_fuel_m3,
        p_recirc_w,
    })
}

struct ReactorResult {
    sigma_v: f64,
    reaction_rate_m3s: f64,
    p_fusion_w: f64,
    p_cascade_w: f64,
    p_cascade_electric_w: f64,
    reactions_per_s: f64,
    fuel_kg_per_year: f64,
    fuel_ng_per_year: f64,
}

fn reactor_level(r: &CascadeReaction, cascade: &CascadeResult, cfg: &RtscConfig) -> ReactorResult {
    let sv = sigma_v_general(
        r.z1,
        r.z2,
        r.a1,
        r.a2,
        cfg.t_kev,
        r.calibrated_sigmav_100kev,
    );

    let q_j = r.q_fusion_mev * MEV_TO_J;
    let reaction_rate_m3s = 0.25 * cfg.n_fuel_m3.powi(2) * sv;
    let p_fusion_w = reaction_rate_m3s * q_j * cfg.volume_m3 * cfg.confinement;
    let reactions_per_s = if q_j > 0.0 { p_fusion_w / q_j } else { 0.0 };

    let p_cascade_w = p_fusion_w * cascade.multiplier;
    let p_cascade_electric_w = p_cascade_w * ELECTRIC_EFF;

    // Fuel consumption: each reaction consumes one pair of antimatter nuclei
    let fuel_mass_per_reaction_kg = cascade.fuel_rest_mass * MEV_TO_J / (3.0e8_f64.powi(2));
    let fuel_kg_per_year = reactions_per_s * fuel_mass_per_reaction_kg * SECONDS_PER_YEAR;
    let fuel_ng_per_year = fuel_kg_per_year * 1.0e12;

    ReactorResult {
        sigma_v: sv,
        reaction_rate_m3s,
        p_fusion_w,
        p_cascade_w,
        p_cascade_electric_w,
        reactions_per_s,
        fuel_kg_per_year,
        fuel_ng_per_year,
    }
}

// ─── E. Production Cost eta Sweep ────────────────────────────────────────────

struct EtaSweepPoint {
    eta: f64,
    p_production_w: f64,
    p_net_w: f64,
    system_q: f64,
}

fn eta_sweep(
    reactor: &ReactorResult,
    cascade: &CascadeResult,
    cfg: &RtscConfig,
) -> Vec<EtaSweepPoint> {
    let fuel_rest_mass_j = cascade.fuel_rest_mass * MEV_TO_J;
    let fuel_power_baseline = fuel_rest_mass_j * reactor.reactions_per_s;

    let mut points = Vec::new();
    // Logarithmic sweep from 1e-4 to 0.5
    let n_points = 30;
    for i in 0..n_points {
        let log_min = (1.0e-4_f64).ln();
        let log_max = (0.5_f64).ln();
        let log_eta = log_min + (i as f64 / (n_points - 1) as f64) * (log_max - log_min);
        let eta = log_eta.exp();

        let p_production_w = fuel_power_baseline / eta;
        let p_net_w = reactor.p_cascade_electric_w - p_production_w - cfg.p_recirc_w;
        let system_q = if p_production_w + cfg.p_recirc_w > 0.0 {
            reactor.p_cascade_electric_w / (p_production_w + cfg.p_recirc_w)
        } else {
            f64::INFINITY
        };

        points.push(EtaSweepPoint {
            eta,
            p_production_w,
            p_net_w,
            system_q,
        });
    }
    points
}

fn find_break_even_eta(
    reactor: &ReactorResult,
    cascade: &CascadeResult,
    cfg: &RtscConfig,
) -> f64 {
    // P_net = 0 => P_cascade_electric = P_production + P_recirc
    // P_cascade_electric - P_recirc = fuel_power / eta
    // eta = fuel_power / (P_cascade_electric - P_recirc)
    let fuel_rest_mass_j = cascade.fuel_rest_mass * MEV_TO_J;
    let fuel_power = fuel_rest_mass_j * reactor.reactions_per_s;
    let available = reactor.p_cascade_electric_w - cfg.p_recirc_w;
    if available <= 0.0 {
        return f64::INFINITY;
    }
    fuel_power / available
}

fn find_crossover_eta(
    reactor: &ReactorResult,
    cascade: &CascadeResult,
    cfg: &RtscConfig,
) -> f64 {
    // P_net_cascade > P_net_regular
    // P_cascade_electric - fuel_power/eta - P_recirc > P_net_regular
    // fuel_power/eta < P_cascade_electric - P_recirc - P_net_regular
    let fuel_rest_mass_j = cascade.fuel_rest_mass * MEV_TO_J;
    let fuel_power = fuel_rest_mass_j * reactor.reactions_per_s;
    let headroom = reactor.p_cascade_electric_w - cfg.p_recirc_w - cfg.p_net_regular_w;
    if headroom <= 0.0 {
        return f64::INFINITY;
    }
    fuel_power / headroom
}

// ─── F & G. Output Generation ────────────────────────────────────────────────

fn main() {
    let out_dir = env::var("GUTOE_CASCADE_SIM_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/antimatter_fusion_cascade".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");

    let reactions = reaction_catalog();
    let cfg = load_rtsc_config();

    // ─── H. Embedded Assertions on Cascade Physics ───────────────────────────

    let mut txt = String::new();
    let mut csv_reactions = String::from(
        "reaction_id,q_fusion_mev,q_from_masses_mev,q_delta_mev,\
         total_annihilation_mev,total_cascade_mev,multiplier,\
         fuel_rest_mass_mev,direct_annihilation_mev,cascade_penalty_mev,\
         break_even_eta,products,\
         sigma_v_m3s,p_fusion_mw,p_cascade_mw,p_cascade_electric_mw,\
         reactions_per_s,fuel_ng_per_year\n",
    );
    let mut csv_eta = String::from(
        "reaction_id,eta,p_production_mw,p_net_mw,system_q\n",
    );

    txt.push_str("╔══════════════════════════════════════════════════════════════════════╗\n");
    txt.push_str("║          ANTIMATTER FUSION CASCADE SIMULATOR                        ║\n");
    txt.push_str("║          GUTOE-Physics — Thermodynamic Analysis                     ║\n");
    txt.push_str("╚══════════════════════════════════════════════════════════════════════╝\n\n");

    txt.push_str("[rtsc_reference_config]\n");
    txt.push_str(&format!("R_major       = {:.3} m\n", cfg.r_major_m));
    txt.push_str(&format!("a_minor       = {:.3} m\n", cfg.a_minor_m));
    txt.push_str(&format!("T             = {:.1} keV\n", cfg.t_kev));
    txt.push_str(&format!("B_operating   = {:.1} T\n", cfg.b_operating_t));
    txt.push_str(&format!("confinement   = {:.4}\n", cfg.confinement));
    txt.push_str(&format!("volume        = {:.4} m^3\n", cfg.volume_m3));
    txt.push_str(&format!("n_fuel        = {:.4e} m^-3\n", cfg.n_fuel_m3));
    txt.push_str(&format!("P_net_regular = {:.3} MW\n", cfg.p_net_regular_w / 1.0e6));
    txt.push_str(&format!("P_recirc      = {:.3} MW\n\n", cfg.p_recirc_w / 1.0e6));

    // ─── Per-reaction cascade analysis ───────────────────────────────────────

    txt.push_str("[cascade_reactions]\n");
    txt.push_str(&format!(
        "{:<45} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n",
        "Reaction", "Q_fus", "Annihil", "Cascade", "Mult", "Direct", "Penalty"
    ));
    txt.push_str(&format!(
        "{:<45} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n",
        "", "(MeV)", "(MeV)", "(MeV)", "(x)", "(MeV)", "(MeV)"
    ));
    txt.push_str(&"-".repeat(105));
    txt.push('\n');

    let mut all_cascades: Vec<(CascadeReaction, CascadeResult, ReactorResult)> = Vec::new();
    let mut assertion_failures = Vec::new();

    for r in &reactions {
        let cascade = compute_cascade(r);
        let reactor = reactor_level(r, &cascade, &cfg);

        // Assertion 1: Q-value cross-check
        let q_delta = (cascade.q_fusion - cascade.q_from_masses).abs();
        if q_delta > 0.1 {
            assertion_failures.push(format!(
                "FAIL: {} Q cross-check: catalog={:.6}, from_masses={:.6}, delta={:.6} > 0.1 MeV",
                r.id, cascade.q_fusion, cascade.q_from_masses, q_delta
            ));
        }

        // Assertion 2: Mass conservation identity: penalty = Q_fusion
        let penalty_delta = (cascade.cascade_penalty - cascade.q_fusion).abs();
        if penalty_delta > 0.01 {
            assertion_failures.push(format!(
                "FAIL: {} mass conservation: penalty={:.6}, Q={:.6}, delta={:.6}",
                r.id, cascade.cascade_penalty, cascade.q_fusion, penalty_delta
            ));
        }

        // Assertion 3: break-even eta > 0.5 (cascade worse than direct annihilation)
        if cascade.break_even_eta <= 0.5 {
            assertion_failures.push(format!(
                "FAIL: {} break_even_eta={:.6} <= 0.5 (cascade should be worse than direct)",
                r.id, cascade.break_even_eta
            ));
        }

        let per_prod: Vec<String> = cascade
            .annihilation_per_product
            .iter()
            .map(|v| format!("{:.1}", v))
            .collect();

        txt.push_str(&format!(
            "{:<45} {:>10.3} {:>10.1} {:>10.1} {:>10.1} {:>10.1} {:>10.3}\n",
            r.id,
            cascade.q_fusion,
            cascade.total_annihilation,
            cascade.total_cascade,
            cascade.multiplier,
            cascade.direct_annihilation,
            cascade.cascade_penalty
        ));
        txt.push_str(&format!(
            "  products: {} | per-product annihilation: [{}] MeV\n",
            r.product_labels,
            per_prod.join(", ")
        ));

        csv_reactions.push_str(&format!(
            "{},{:.6},{:.6},{:.6},{:.3},{:.3},{:.3},{:.3},{:.3},{:.6},{:.8},{},{:.6e},{:.6},{:.6},{:.6},{:.6e},{:.4}\n",
            r.id,
            cascade.q_fusion,
            cascade.q_from_masses,
            q_delta,
            cascade.total_annihilation,
            cascade.total_cascade,
            cascade.multiplier,
            cascade.fuel_rest_mass,
            cascade.direct_annihilation,
            cascade.cascade_penalty,
            cascade.break_even_eta,
            r.product_labels,
            reactor.sigma_v,
            reactor.p_fusion_w / 1.0e6,
            reactor.p_cascade_w / 1.0e6,
            reactor.p_cascade_electric_w / 1.0e6,
            reactor.reactions_per_s,
            reactor.fuel_ng_per_year
        ));

        all_cascades.push((r.clone(), cascade, reactor));
    }

    // ─── Assertion report ────────────────────────────────────────────────────

    txt.push('\n');
    txt.push_str("[assertions]\n");
    if assertion_failures.is_empty() {
        txt.push_str("ALL ASSERTIONS PASSED\n");
        txt.push_str("  - Q cross-check: |Q_catalog - Q_from_masses| < 0.1 MeV for all reactions\n");
        txt.push_str("  - Mass conservation: cascade_penalty = Q_fusion (within 0.01 MeV)\n");
        txt.push_str("  - Break-even: eta > 0.5 for all reactions (cascade strictly worse than direct annihilation)\n");
    } else {
        for f in &assertion_failures {
            txt.push_str(f);
            txt.push('\n');
        }
    }
    txt.push('\n');

    // ─── Thermodynamic verdict ───────────────────────────────────────────────

    txt.push_str("[thermodynamic_verdict]\n");
    txt.push_str("The antimatter fusion cascade is thermodynamically INFERIOR to direct annihilation.\n");
    txt.push_str("Mass conservation forces: cascade_total = direct_annihilation - Q_fusion.\n");
    txt.push_str("The fusion step burns rest mass into kinetic energy, losing the 2x amplification\n");
    txt.push_str("on that portion. The cascade IS a 531x multiplier vs regular D+T fusion,\n");
    txt.push_str("but offers no thermodynamic edge over simply annihilating the raw fuel.\n\n");

    txt.push_str(&format!(
        "{:<45} {:>12} {:>12} {:>12}\n",
        "Reaction", "eta_break", "eta_direct", "delta_eta"
    ));
    txt.push_str(&"-".repeat(81));
    txt.push('\n');

    for (r, cascade, _) in &all_cascades {
        let eta_direct = 0.5; // direct annihilation break-even
        let delta = cascade.break_even_eta - eta_direct;
        txt.push_str(&format!(
            "{:<45} {:>12.6} {:>12.6} {:>12.6}\n",
            r.id, cascade.break_even_eta, eta_direct, delta
        ));
    }
    txt.push('\n');

    // ─── Reactor-level comparison table ──────────────────────────────────────

    txt.push_str("[reactor_level_comparison]\n");
    txt.push_str(&format!(
        "{:<45} {:>12} {:>12} {:>12} {:>12} {:>12}\n",
        "Reaction", "sigma_v", "P_fus(MW)", "P_cas(MW)", "P_cas_e(MW)", "fuel(ng/yr)"
    ));
    txt.push_str(&"-".repeat(105));
    txt.push('\n');

    for (r, _, reactor) in &all_cascades {
        txt.push_str(&format!(
            "{:<45} {:>12.4e} {:>12.4} {:>12.1} {:>12.1} {:>12.4}\n",
            r.id,
            reactor.sigma_v,
            reactor.p_fusion_w / 1.0e6,
            reactor.p_cascade_w / 1.0e6,
            reactor.p_cascade_electric_w / 1.0e6,
            reactor.fuel_ng_per_year
        ));
    }
    txt.push('\n');

    // ─── Power density comparison ────────────────────────────────────────────

    let rtsc_power_density = cfg.p_net_regular_w / cfg.volume_m3 / 1.0e6;
    txt.push_str("[power_density_comparison]\n");
    txt.push_str(&format!("RTSC D+He3 fusion power density: {:.4} MW/m^3\n", rtsc_power_density));
    txt.push_str(&format!(
        "{:<45} {:>12} {:>12} {:>12}\n",
        "Reaction", "eta=0.3%", "eta=16.7%", "eta=50%"
    ));
    txt.push_str(&"-".repeat(81));
    txt.push('\n');

    let eta_scenarios = [0.003, 0.167, 0.5];
    for (r, cascade, reactor) in &all_cascades {
        let fuel_rest_mass_j = cascade.fuel_rest_mass * MEV_TO_J;
        let fuel_power_base = fuel_rest_mass_j * reactor.reactions_per_s;
        let mut densities = Vec::new();
        for &eta in &eta_scenarios {
            let p_prod = fuel_power_base / eta;
            let p_net = reactor.p_cascade_electric_w - p_prod - cfg.p_recirc_w;
            let density = p_net / cfg.volume_m3 / 1.0e6;
            densities.push(density);
        }
        txt.push_str(&format!(
            "{:<45} {:>12.4} {:>12.4} {:>12.4}\n",
            r.id, densities[0], densities[1], densities[2]
        ));
    }
    txt.push('\n');

    // ─── Eta sweep per reaction ──────────────────────────────────────────────

    txt.push_str("[eta_sweep_summary]\n");
    txt.push_str(&format!(
        "{:<45} {:>12} {:>12} {:>12} {:>12} {:>12}\n",
        "Reaction", "eta_break", "eta_cross", "Q@0.3%", "Q@16.7%", "Q@50%"
    ));
    txt.push_str(&"-".repeat(105));
    txt.push('\n');

    for (r, cascade, reactor) in &all_cascades {
        let be = find_break_even_eta(reactor, cascade, &cfg);
        let co = find_crossover_eta(reactor, cascade, &cfg);

        let sweep = eta_sweep(reactor, cascade, &cfg);
        for pt in &sweep {
            csv_eta.push_str(&format!(
                "{},{:.8e},{:.6},{:.6},{:.6}\n",
                r.id,
                pt.eta,
                pt.p_production_w / 1.0e6,
                pt.p_net_w / 1.0e6,
                pt.system_q
            ));
        }

        // Compute system Q at the three scenarios
        let fuel_rest_mass_j = cascade.fuel_rest_mass * MEV_TO_J;
        let fuel_power_base = fuel_rest_mass_j * reactor.reactions_per_s;
        let mut qs = Vec::new();
        for &eta in &eta_scenarios {
            let p_prod = fuel_power_base / eta;
            let sys_q = if p_prod + cfg.p_recirc_w > 0.0 {
                reactor.p_cascade_electric_w / (p_prod + cfg.p_recirc_w)
            } else {
                f64::INFINITY
            };
            qs.push(sys_q);
        }

        txt.push_str(&format!(
            "{:<45} {:>12.6} {:>12.6} {:>12.4} {:>12.4} {:>12.4}\n",
            r.id,
            be,
            co,
            qs[0],
            qs[1],
            qs[2]
        ));
    }
    txt.push('\n');

    // ─── Consumption rates ───────────────────────────────────────────────────

    txt.push_str("[consumption_rates]\n");
    txt.push_str(&format!(
        "{:<45} {:>15} {:>15} {:>15} {:>15}\n",
        "Reaction", "rate(m^-3 s^-1)", "kg/year", "ng/year", "mg/year"
    ));
    txt.push_str(&"-".repeat(105));
    txt.push('\n');

    for (r, _, reactor) in &all_cascades {
        let ng = reactor.fuel_ng_per_year;
        let mg = ng / 1.0e6;
        txt.push_str(&format!(
            "{:<45} {:>15.4e} {:>15.4e} {:>15.4} {:>15.9}\n",
            r.id, reactor.reaction_rate_m3s, reactor.fuel_kg_per_year, ng, mg
        ));
    }
    txt.push('\n');

    // ─── JSON output ─────────────────────────────────────────────────────────

    let json_reactions: Vec<serde_json::Value> = all_cascades
        .iter()
        .map(|(r, cascade, reactor)| {
            serde_json::json!({
                "reaction_id": r.id,
                "q_fusion_mev": cascade.q_fusion,
                "total_annihilation_mev": cascade.total_annihilation,
                "total_cascade_mev": cascade.total_cascade,
                "multiplier": cascade.multiplier,
                "fuel_rest_mass_mev": cascade.fuel_rest_mass,
                "direct_annihilation_mev": cascade.direct_annihilation,
                "cascade_penalty_mev": cascade.cascade_penalty,
                "break_even_eta": cascade.break_even_eta,
                "sigma_v_m3s": reactor.sigma_v,
                "p_fusion_mw": reactor.p_fusion_w / 1.0e6,
                "p_cascade_mw": reactor.p_cascade_w / 1.0e6,
                "p_cascade_electric_mw": reactor.p_cascade_electric_w / 1.0e6,
                "reactions_per_s": reactor.reactions_per_s,
                "fuel_ng_per_year": reactor.fuel_ng_per_year,
                "break_even_eta_cascade": find_break_even_eta(reactor, cascade, &cfg),
                "crossover_eta": find_crossover_eta(reactor, cascade, &cfg),
            })
        })
        .collect();

    let json_out = serde_json::json!({
        "simulator": "antimatter_fusion_cascade",
        "rtsc_config": {
            "r_major_m": cfg.r_major_m,
            "a_minor_m": cfg.a_minor_m,
            "t_kev": cfg.t_kev,
            "b_operating_t": cfg.b_operating_t,
            "confinement": cfg.confinement,
            "volume_m3": cfg.volume_m3,
            "p_net_regular_mw": cfg.p_net_regular_w / 1.0e6,
            "p_recirc_mw": cfg.p_recirc_w / 1.0e6,
        },
        "reactions": json_reactions,
        "assertions_passed": assertion_failures.is_empty(),
        "thermodynamic_verdict": "cascade is INFERIOR to direct annihilation by exactly Q_fusion per reaction",
    });

    // ─── Write all outputs ───────────────────────────────────────────────────

    let txt_path = format!("{out_dir}/antimatter_fusion_cascade.txt");
    let csv_r_path = format!("{out_dir}/antimatter_fusion_cascade_reactions.csv");
    let csv_e_path = format!("{out_dir}/antimatter_fusion_cascade_eta_sweep.csv");
    let json_path = format!("{out_dir}/antimatter_fusion_cascade.json");

    fs::write(&txt_path, &txt).expect("write txt");
    fs::write(&csv_r_path, &csv_reactions).expect("write reactions csv");
    fs::write(&csv_e_path, &csv_eta).expect("write eta sweep csv");

    let json_str = serde_json::to_string_pretty(&json_out).expect("json serialize");
    fs::write(&json_path, &json_str).expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {csv_r_path}");
    println!("wrote {csv_e_path}");
    println!("wrote {json_path}");

    // Hard-fail on assertion violations
    if !assertion_failures.is_empty() {
        eprintln!("\nASSERTION FAILURES:");
        for f in &assertion_failures {
            eprintln!("  {f}");
        }
        std::process::exit(1);
    }
}
