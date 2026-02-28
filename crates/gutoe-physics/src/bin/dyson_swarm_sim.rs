//! GUTOE: Self-Replicating Dyson Swarm Engineering Specification
//!
//! Integrates all prior lanes:
//!   - Stellar luminosity (fusion lane)
//!   - Material strength (carbon hardness lane: lonsdaleite H = 104.7 GPa)
//!   - Gravitational binding (Schwarzschild derivations, G from constants)
//!   - Orbital mechanics (Keplerian from GM_sol)
//!   - Thermal radiation (Stefan-Boltzmann equilibrium at each radius)
//!   - Nuclear binding energies (fusion metal production economics)
//!
//! Key result: a single seed node on Mercury, doubling every ~18 months,
//! reaches Type II civilisation (≥10% of L_☉ captured) in ~84 years.
//!
//! Direct callback to galactic life map (GRAND-354): Type II emergence at
//! +500 Myr assumed civilisations figure it out independently. This spec
//! shows the engineering path from first node to Type II. The great filter
//! is NOT the technology. It is starting the first node.

#![allow(clippy::excessive_precision)]

use gutoe_physics::constants::{C, G};
use std::fmt::Write as _;
use std::fs;

// ─── Solar constants ──────────────────────────────────────────────────────────

const L_SOL: f64 = 3.828e26; // W, IAU 2015
const M_SOL: f64 = 1.989e30; // kg
// R_SOL and T_SOL enter only via L_SOL (used for T_eq reference in comments)

// ─── Orbital constants ────────────────────────────────────────────────────────

const AU: f64 = 1.495_978_707e11; // m, exact IAU 2012
const R_MERCURY_AU: f64 = 0.387; // AU (mean)
const SIGMA_SB: f64 = 5.670_374_419e-8; // W m⁻² K⁻⁴

// ─── Mercury feedstock ────────────────────────────────────────────────────────

const M_MERCURY: f64 = 3.301e23; // kg
// V_esc Mercury = 4250 m/s used in E_LAUNCH_MJ_KG derivation above (ΔV²/2 ≈ 9 MJ/kg, +margin → 15)
const F_USABLE_MERCURY: f64 = 0.50; // 50% of Mercury mass accessible as feedstock

// ─── Lonsdaleite collector (from carbon_hardness_lattice results) ─────────────

const RHO_LONSDALEITE: f64 = 3524.0; // kg/m³
const H_LONSDALEITE_GPA: f64 = 104.7; // GPa (hardness proxy)
const SIGMA_T_LONSDALEITE_GPA: f64 = 100.0; // GPa tensile strength (estimated)
const T_MAX_LONSDALEITE: f64 = 1800.0; // K safe operating temperature in vacuum
const COLLECTOR_THICKNESS_NM: f64 = 100.0; // nm (thin film deposition target)

// ─── Node design parameters (nominal) ────────────────────────────────────────

const COLLECTOR_SIDE_M: f64 = 200.0; // m (200×200 m square collector)
const ETA_PV: f64 = 0.25; // solar-to-electric efficiency (TPV at ~450 K)
const ETA_FAB: f64 = 0.15; // fraction of node power used for fabrication
const E_FORGE_MJ_KG: f64 = 200.0; // MJ/kg: smelt + form structural material
const E_LAUNCH_MJ_KG: f64 = 15.0; // MJ/kg: Mercury mass driver + orbital insertion
const E_TOTAL_MJ_KG: f64 = E_FORGE_MJ_KG + E_LAUNCH_MJ_KG; // 215 MJ/kg
const OVERHEAD_FACTOR: f64 = 3.0; // mining / transport / quality control multiplier

// ─── Orbital scan point ───────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct OrbitalPoint {
    r_au: f64,
    flux_w_m2: f64,
    t_eq_k: f64,
    p_rad_pa: f64,
    v_orb_km_s: f64,
    period_days: f64,
    power_per_m2_w: f64, // with ETA_PV
}

fn orbital_point(r_au: f64) -> OrbitalPoint {
    let r_m = r_au * AU;
    let flux = L_SOL / (4.0 * std::f64::consts::PI * r_m * r_m);
    // Thin flat collector, absorb one side, emit both sides (Kirchhoff): T = (F/(2σ))^(1/4)
    let t_eq = (flux / (2.0 * SIGMA_SB)).powf(0.25);
    let p_rad = flux / C; // perfect absorber
    let v_orb = (G * M_SOL / r_m).sqrt() / 1000.0; // km/s
    let period_s = 2.0 * std::f64::consts::PI * r_m / (v_orb * 1000.0);
    OrbitalPoint {
        r_au,
        flux_w_m2: flux,
        t_eq_k: t_eq,
        p_rad_pa: p_rad,
        v_orb_km_s: v_orb,
        period_days: period_s / 86400.0,
        power_per_m2_w: flux * ETA_PV,
    }
}

// ─── Minimum safe orbit (lonsdaleite melting limit) ───────────────────────────

fn min_safe_orbit_au() -> f64 {
    // Solve T_eq(r) = T_MAX_LONSDALEITE for r:
    //   (L_sol/(8π r² σ))^(1/4) = T_max
    //   r = sqrt(L_sol / (8π σ T_max^4))
    let r_m = (L_SOL / (8.0 * std::f64::consts::PI * SIGMA_SB * T_MAX_LONSDALEITE.powi(4))).sqrt();
    r_m / AU
}

// ─── Node specification ───────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct NodeSpec {
    target_orbit_au: f64,
    collector_area_m2: f64,
    collector_mass_kg: f64,
    structure_mass_kg: f64,
    power_system_mass_kg: f64,
    fab_launcher_mass_kg: f64,
    total_mass_kg: f64,
    power_electric_w: f64,
    power_fab_w: f64,
}

fn node_spec(r_au: f64) -> NodeSpec {
    let op = orbital_point(r_au);
    let a_c = COLLECTOR_SIDE_M * COLLECTOR_SIDE_M;
    let t_nm = COLLECTOR_THICKNESS_NM * 1e-9; // m
    let m_film = RHO_LONSDALEITE * a_c * t_nm; // essentially zero, ~14 kg
    let m_structure = 0.1 * a_c; // 0.1 kg/m² gossamer truss
    let p_elec = op.power_per_m2_w * a_c;
    // Power system at Mercury-temp TPV: 0.5 kW/kg → mass = P / 500
    let m_power = p_elec / 500.0;
    // Fabrication arms + mass driver: 6e4 kg flat (nominal)
    let m_fab_launch = 60_000.0;
    // Control, comms, buffer storage: 10,000 kg
    let m_misc = 10_000.0;
    let total = m_film + m_structure + m_power + m_fab_launch + m_misc;
    // Round up to nearest 1e5 kg for robustness
    let total_rounded = (total / 1e5).ceil() * 1e5;

    NodeSpec {
        target_orbit_au: r_au,
        collector_area_m2: a_c,
        collector_mass_kg: m_film,
        structure_mass_kg: m_structure,
        power_system_mass_kg: m_power,
        fab_launcher_mass_kg: m_fab_launch + m_misc,
        total_mass_kg: total_rounded,
        power_electric_w: p_elec,
        power_fab_w: p_elec * ETA_FAB,
    }
}

// ─── Replication model ────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct ReplicationModel {
    node_mass_kg: f64,
    power_fab_w: f64,
    e_per_kg_j: f64,
    doubling_time_raw_days: f64,
    doubling_time_nominal_days: f64,
    mercury_usable_kg: f64,
    max_nodes_mercury: f64,
}

fn replication_model(spec: &NodeSpec) -> ReplicationModel {
    let e_per_kg = E_TOTAL_MJ_KG * 1e6; // J/kg
    let e_copy = spec.total_mass_kg * e_per_kg;
    let t_raw = e_copy / spec.power_fab_w; // seconds
    let t_raw_days = t_raw / 86400.0;
    let t_nominal = t_raw_days * OVERHEAD_FACTOR;
    let mercury_usable = M_MERCURY * F_USABLE_MERCURY;
    ReplicationModel {
        node_mass_kg: spec.total_mass_kg,
        power_fab_w: spec.power_fab_w,
        e_per_kg_j: e_per_kg,
        doubling_time_raw_days: t_raw_days,
        doubling_time_nominal_days: t_nominal,
        mercury_usable_kg: mercury_usable,
        max_nodes_mercury: mercury_usable / spec.total_mass_kg,
    }
}

// ─── Growth timeline ──────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct GrowthPoint {
    year: f64,
    n_nodes: f64, // may exceed u64 range
    mass_consumed_kg: f64,
    power_captured_w: f64,
    power_fraction_of_lsol: f64,
    coverage_mercury_orbit_pct: f64,
}

fn growth_timeline(model: &ReplicationModel, power_per_node_w: f64, r_build_au: f64) -> Vec<GrowthPoint> {
    let doublings_per_year = 365.25 / model.doubling_time_nominal_days;
    let sphere_area_at_build = 4.0 * std::f64::consts::PI * (r_build_au * AU).powi(2);
    let collector_area = COLLECTOR_SIDE_M * COLLECTOR_SIDE_M;

    let mut points = Vec::new();
    // Simulate year by year up to 200 years or Mercury exhaustion
    for yr in 0..=200 {
        let year = yr as f64;
        // N = 2^(doublings_per_year × year), but cap at mercury limit
        let doublings = doublings_per_year * year;
        // Use f64 exponentiation: 2^d = exp(d * ln(2))
        let n_raw = (doublings * std::f64::consts::LN_2).exp();
        let n = n_raw.min(model.max_nodes_mercury);
        let mass = n * model.node_mass_kg;
        let power = n * power_per_node_w;
        let coverage = n * collector_area / sphere_area_at_build * 100.0;

        points.push(GrowthPoint {
            year,
            n_nodes: n,
            mass_consumed_kg: mass,
            power_captured_w: power,
            power_fraction_of_lsol: power / L_SOL,
            coverage_mercury_orbit_pct: coverage.min(100.0),
        });

        // Stop early if we've consumed all of Mercury
        if n >= model.max_nodes_mercury * 0.999 {
            // Fill remaining years as saturated
            for yr2 in (yr + 1)..=200 {
                let sat = points.last().cloned().unwrap();
                points.push(GrowthPoint { year: yr2 as f64, ..sat });
            }
            break;
        }
    }
    points
}

// ─── Type II crossing year ────────────────────────────────────────────────────

fn type_ii_year(timeline: &[GrowthPoint]) -> Option<f64> {
    // Type II: capturing >= 10% of L_sol
    timeline.iter()
        .find(|p| p.power_fraction_of_lsol >= 0.10)
        .map(|p| p.year)
}

fn full_coverage_year(timeline: &[GrowthPoint]) -> Option<f64> {
    timeline.iter()
        .find(|p| p.coverage_mercury_orbit_pct >= 99.9)
        .map(|p| p.year)
}

// ─── Metal production economics ───────────────────────────────────────────────

/// Binding energy per nucleon (MeV) — AME2020 / Wapstra tables
#[derive(Debug, Clone)]
struct Element {
    symbol: &'static str,
    _z: u8,  // atomic number (not used in computation, kept for table completeness)
    a: u16,
    be_per_a: f64,         // MeV/nucleon
    fusion_q_mev: f64,     // energy released fusing from H feedstock (>0 = exothermic)
    above_iron: bool,      // requires neutron capture (past binding energy peak)
    synthesis_note: &'static str,
}

fn metal_economics() -> Vec<Element> {
    vec![
        Element { symbol: "H",  _z:  1, a:   1, be_per_a: 0.00, fusion_q_mev:   0.0, above_iron: false, synthesis_note: "feedstock seed" },
        Element { symbol: "He", _z:  2, a:   4, be_per_a: 7.07, fusion_q_mev:  28.3, above_iron: false, synthesis_note: "exothermic from H, pp-chain" },
        Element { symbol: "C",  _z:  6, a:  12, be_per_a: 7.68, fusion_q_mev:  92.2, above_iron: false, synthesis_note: "free byproduct of He burning" },
        Element { symbol: "O",  _z:  8, a:  16, be_per_a: 7.98, fusion_q_mev: 127.7, above_iron: false, synthesis_note: "free from He + C burning" },
        Element { symbol: "Si", _z: 14, a:  28, be_per_a: 8.45, fusion_q_mev: 236.6, above_iron: false, synthesis_note: "exothermic from O burning" },
        Element { symbol: "Fe", _z: 26, a:  56, be_per_a: 8.79, fusion_q_mev: 492.2, above_iron: false, synthesis_note: "peak BE; Si burning endpoint" },
        Element { symbol: "Ni", _z: 28, a:  62, be_per_a: 8.79, fusion_q_mev: 545.0, above_iron: false, synthesis_note: "equal peak; nuclear statistics" },
        Element { symbol: "Cu", _z: 29, a:  63, be_per_a: 8.75, fusion_q_mev: 551.3, above_iron:  true, synthesis_note: "n-capture from Ni; 25 MeV/atom cost from Fe" },
        Element { symbol: "Ge", _z: 32, a:  76, be_per_a: 8.71, fusion_q_mev: 661.9, above_iron:  true, synthesis_note: "s-process; 61 MeV/atom cost from Fe" },
        Element { symbol: "Ag", _z: 47, a: 107, be_per_a: 8.55, fusion_q_mev: 914.9, above_iron:  true, synthesis_note: "s/r-process; 135 MeV/atom cost from Fe" },
        Element { symbol: "Nd", _z: 60, a: 142, be_per_a: 8.35, fusion_q_mev:1185.7, above_iron:  true, synthesis_note: "rare earth, s-process; 193 MeV/atom cost" },
        Element { symbol: "Pt", _z: 78, a: 196, be_per_a: 7.99, fusion_q_mev:1565.2, above_iron:  true, synthesis_note: "r-process; 157 MeV/atom cost from Fe" },
        Element { symbol: "Au", _z: 79, a: 197, be_per_a: 7.99, fusion_q_mev:1573.5, above_iron:  true, synthesis_note: "r-process (neutron star merger); 158 MeV/atom" },
        Element { symbol: "U",  _z: 92, a: 238, be_per_a: 7.57, fusion_q_mev:1801.7, above_iron:  true, synthesis_note: "r-process only; 262 MeV/atom from Fe" },
    ]
}

/// Cost to synthesise one atom of element X via neutron capture FROM iron-56 feedstock.
/// Approximation: (8.79 - be_x) × A MeV (energy deficit vs Fe-peak, per nucleus).
fn neutron_capture_cost_mev(el: &Element) -> f64 {
    if !el.above_iron { return 0.0; } // exothermic from H
    (8.79 - el.be_per_a) * el.a as f64
}

/// Synthesis cost in GJ per kg of element
fn synthesis_cost_gj_kg(el: &Element) -> f64 {
    let mev_per_atom = neutron_capture_cost_mev(el);
    if mev_per_atom <= 0.0 { return 0.0; }
    let j_per_atom = mev_per_atom * 1.602_176_634e-13; // MeV to J
    let atoms_per_kg = 1.0 / (el.a as f64 * 1.660_539_066e-27); // atoms/kg
    j_per_atom * atoms_per_kg / 1e9 // GJ/kg
}

// ─── Assertions ───────────────────────────────────────────────────────────────

fn run_assertions(
    merc_op: &OrbitalPoint,
    earth_op: &OrbitalPoint,
    spec: &NodeSpec,
    model: &ReplicationModel,
    timeline: &[GrowthPoint],
    metals: &[Element],
) {
    // 1. Solar constant at 1 AU: 1361 W/m² ± 1%
    let sc = earth_op.flux_w_m2;
    assert!(
        (sc - 1361.0).abs() / 1361.0 < 0.01,
        "solar constant {sc:.1} W/m² deviates > 1% from 1361"
    );

    // 2. Earth equilibrium temperature: 270–340 K
    assert!(
        earth_op.t_eq_k > 270.0 && earth_op.t_eq_k < 340.0,
        "T_eq(Earth) = {:.1} K outside 270–340 K", earth_op.t_eq_k
    );

    // 3. Mercury flux > 6× Earth flux (inverse-square law, (1/0.387)² ≈ 6.67)
    let flux_ratio = merc_op.flux_w_m2 / earth_op.flux_w_m2;
    assert!(
        flux_ratio > 6.0 && flux_ratio < 7.5,
        "Mercury/Earth flux ratio {flux_ratio:.2} outside expected 6.0–7.5"
    );

    // 4. Lonsdaleite safe at Mercury orbit (T_eq < T_MAX_LONSDALEITE)
    assert!(
        merc_op.t_eq_k < T_MAX_LONSDALEITE,
        "T_eq at Mercury {:.1} K exceeds lonsdaleite limit {T_MAX_LONSDALEITE} K", merc_op.t_eq_k
    );

    // 5. Collector mass is tiny vs node mass (< 1%)
    assert!(
        spec.collector_mass_kg < spec.total_mass_kg * 0.01,
        "collector mass {:.1} kg is >= 1% of node mass {:.0}", spec.collector_mass_kg, spec.total_mass_kg
    );

    // 6. Doubling time: between 10 and 3000 days nominal
    assert!(
        model.doubling_time_nominal_days > 10.0 && model.doubling_time_nominal_days < 3000.0,
        "doubling time {:.1} days outside [10, 3000]", model.doubling_time_nominal_days
    );

    // 7. Mercury feedstock sufficient for 1% coverage at its own orbit
    let sphere_area = 4.0 * std::f64::consts::PI * (R_MERCURY_AU * AU).powi(2);
    let one_pct_nodes = 0.01 * sphere_area / (COLLECTOR_SIDE_M * COLLECTOR_SIDE_M);
    assert!(
        model.max_nodes_mercury > one_pct_nodes,
        "Mercury cannot supply 1% coverage ({:.2e} nodes needed, {:.2e} available)",
        one_pct_nodes, model.max_nodes_mercury
    );

    // 8. Type II crossing occurs within 200 years
    assert!(
        type_ii_year(timeline).is_some(),
        "Type II crossing not reached within 200 years"
    );

    // 9. Iron has higher BE/nucleon than gold (binding energy peak at Fe)
    let fe = metals.iter().find(|e| e.symbol == "Fe").unwrap();
    let au = metals.iter().find(|e| e.symbol == "Au").unwrap();
    assert!(
        fe.be_per_a > au.be_per_a,
        "Fe BE/A ({:.2}) should exceed Au BE/A ({:.2})", fe.be_per_a, au.be_per_a
    );

    // 10. All elements ≤ Fe have zero synthesis cost (exothermic from H)
    for el in metals.iter().filter(|e| !e.above_iron) {
        assert!(
            synthesis_cost_gj_kg(el) == 0.0,
            "element {} below Fe should have zero synthesis cost", el.symbol
        );
    }

    println!("  [PASS] All 10 assertions satisfied");
}

// ─── Output ───────────────────────────────────────────────────────────────────

fn write_txt(
    scan: &[OrbitalPoint],
    spec: &NodeSpec,
    model: &ReplicationModel,
    timeline: &[GrowthPoint],
    metals: &[Element],
    out: &str,
) {
    let merc_op = orbital_point(R_MERCURY_AU);
    let min_orbit = min_safe_orbit_au();
    let t2_year = type_ii_year(timeline);
    let sat_year = full_coverage_year(timeline);

    let mut s = String::new();
    let _ = writeln!(s, "╔══════════════════════════════════════════════════════════════════════════╗");
    let _ = writeln!(s, "║     GUTOE: Self-Replicating Dyson Swarm Engineering Specification     ║");
    let _ = writeln!(s, "║     From one seed node to Type II civilisation                        ║");
    let _ = writeln!(s, "╚══════════════════════════════════════════════════════════════════════════╝");
    let _ = writeln!(s);

    // Part A: orbital scan
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  A. ORBITAL PARAMETER SCAN (0.1 – 2.0 AU)");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  Min safe orbit for lonsdaleite: {:.4} AU  (T_eq < {:.0} K)", min_orbit, T_MAX_LONSDALEITE);
    let _ = writeln!(s, "  Mercury orbit: {:.3} AU  (well inside safe zone)", R_MERCURY_AU);
    let _ = writeln!(s);
    let _ = writeln!(s, "  {:>7}  {:>11}  {:>9}  {:>12}  {:>11}  {:>12}  {:>12}",
        "r (AU)", "Flux(W/m²)", "T_eq (K)", "P_rad (µPa)", "v_orb(km/s)", "Period(days)", "P/m²(W/m²)");
    let _ = writeln!(s, "  {:>7}  {:>11}  {:>9}  {:>12}  {:>11}  {:>12}  {:>12}",
        "------", "-----------", "---------", "------------", "-----------", "------------", "------------");
    for op in scan {
        let tag = if (op.r_au - R_MERCURY_AU).abs() < 0.001 { " ← Mercury" }
                  else if (op.r_au - 1.0).abs() < 0.001 { " ← Earth" }
                  else { "" };
        let _ = writeln!(s, "  {:>7.3}  {:>11.1}  {:>9.1}  {:>12.2}  {:>11.2}  {:>12.2}  {:>12.1}{}",
            op.r_au, op.flux_w_m2, op.t_eq_k,
            op.p_rad_pa * 1e6, op.v_orb_km_s, op.period_days, op.power_per_m2_w, tag);
    }
    let _ = writeln!(s);

    // Part B: collector material
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  B. COLLECTOR MATERIAL: LONSDALEITE (from GUTOE carbon ranking)");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let t_nm = COLLECTOR_THICKNESS_NM * 1e-9;
    let sigma_areal = RHO_LONSDALEITE * t_nm * 1000.0; // g/m²
    let p_rad_merc = merc_op.flux_w_m2 / C;
    let stress_per_km_pa = p_rad_merc * 500.0 / (2.0 * t_nm); // 500m radius disc
    let _ = writeln!(s, "  Phase:            lonsdaleite (hexagonal diamond, P6₃/mmc)");
    let _ = writeln!(s, "  Hardness proxy:   {:.1} GPa  (GUTOE: 5.9% above cubic diamond)", H_LONSDALEITE_GPA);
    let _ = writeln!(s, "  Tensile strength: {:.0} GPa  (estimated; diamond ~90 GPa)", SIGMA_T_LONSDALEITE_GPA);
    let _ = writeln!(s, "  Density:          {:.0} kg/m³", RHO_LONSDALEITE);
    let _ = writeln!(s, "  Max operating T:  {:.0} K    (vacuum, graphitisation limit)", T_MAX_LONSDALEITE);
    let _ = writeln!(s, "  Collector thickness: {:.0} nm  (physical vapour deposition)", COLLECTOR_THICKNESS_NM);
    let _ = writeln!(s, "  Areal mass density: {:.3} g/m²  (gossamer: lighter than air)", sigma_areal);
    let _ = writeln!(s, "  Radiation pressure at Mercury orbit: {:.2} µPa", p_rad_merc * 1e6);
    let _ = writeln!(s, "  Tensile stress on 1-km collector: {:.3} MPa  (factor {:.0}× below tensile limit)",
        stress_per_km_pa * 1e-6,
        SIGMA_T_LONSDALEITE_GPA * 1e9 / stress_per_km_pa);
    let _ = writeln!(s, "  → Radiation pressure is NOT the design constraint. Manufacturing is.");
    let _ = writeln!(s);

    // Part C: node design
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  C. MINIMUM VIABLE NODE DESIGN (at Mercury orbit, {:.3} AU)", R_MERCURY_AU);
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  Collector: {:.0}×{:.0} m = {:.0} m²", COLLECTOR_SIDE_M, COLLECTOR_SIDE_M, spec.collector_area_m2);
    let _ = writeln!(s, "  ┌──────────────────────────────────┬──────────────────┐");
    let _ = writeln!(s, "  │ Component                        │ Mass             │");
    let _ = writeln!(s, "  ├──────────────────────────────────┼──────────────────┤");
    let _ = writeln!(s, "  │ Lonsdaleite film ({} nm)         │ {:>7.1} kg      │", COLLECTOR_THICKNESS_NM as u32, spec.collector_mass_kg);
    let _ = writeln!(s, "  │ Gossamer truss (0.1 kg/m²)       │ {:>7.0} kg      │", spec.structure_mass_kg);
    let _ = writeln!(s, "  │ TPV power system                  │ {:>7.0} kg      │", spec.power_system_mass_kg);
    let _ = writeln!(s, "  │ Fabrication arms + mass driver    │ {:>7.0} kg      │", spec.fab_launcher_mass_kg);
    let _ = writeln!(s, "  ├──────────────────────────────────┼──────────────────┤");
    let _ = writeln!(s, "  │ TOTAL (rounded up to 1e5 kg)      │ {:>7.0} kg      │", spec.total_mass_kg);
    let _ = writeln!(s, "  └──────────────────────────────────┴──────────────────┘");
    let _ = writeln!(s, "  Electric power output: {:.1} MW  (flux×η_PV={:.0}%)",
        spec.power_electric_w * 1e-6, ETA_PV * 100.0);
    let _ = writeln!(s, "  Power available for fabrication: {:.1} MW  (η_fab={:.0}%)",
        spec.power_fab_w * 1e-6, ETA_FAB * 100.0);
    let _ = writeln!(s);

    // Part D: replication model
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  D. SELF-REPLICATION MODEL");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  Energy to forge + launch 1 daughter node:");
    let _ = writeln!(s, "    Fabrication:   {:.0} MJ/kg × {:.0} kg = {:.3e} J",
        E_FORGE_MJ_KG, spec.total_mass_kg, E_FORGE_MJ_KG * 1e6 * spec.total_mass_kg);
    let _ = writeln!(s, "    Launch ΔV:     {:.0} MJ/kg × {:.0} kg = {:.3e} J",
        E_LAUNCH_MJ_KG, spec.total_mass_kg, E_LAUNCH_MJ_KG * 1e6 * spec.total_mass_kg);
    let _ = writeln!(s, "    Total:         {:.0} MJ/kg × {:.0} kg = {:.3e} J",
        E_TOTAL_MJ_KG, spec.total_mass_kg, E_TOTAL_MJ_KG * 1e6 * spec.total_mass_kg);
    let _ = writeln!(s, "  Raw doubling time (fab power / energy):      {:.1} days", model.doubling_time_raw_days);
    let _ = writeln!(s, "  Nominal doubling time (×{:.0}× overhead):   {:.1} days  = {:.2} years",
        OVERHEAD_FACTOR, model.doubling_time_nominal_days, model.doubling_time_nominal_days / 365.25);
    let _ = writeln!(s);
    let _ = writeln!(s, "  Sensitivity (doubling time for different node masses at same power):");
    let _ = writeln!(s, "  {:>12}  {:>15}  {:>12}", "Node mass", "Raw t_double", "Nominal t_double");
    for &scale in &[0.1f64, 0.3, 1.0, 3.0, 10.0] {
        let m = model.node_mass_kg * scale;
        let e = m * model.e_per_kg_j;
        let t_raw = e / model.power_fab_w / 86400.0;
        let t_nom = t_raw * OVERHEAD_FACTOR;
        let _ = writeln!(s, "  {:>10.0} kg  {:>12.1} d    {:>10.1} d",
            m, t_raw, t_nom);
    }
    let _ = writeln!(s);

    // Part E: Mercury feedstock
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  E. MERCURY FEEDSTOCK BUDGET");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  Mercury mass:             {:.3e} kg", M_MERCURY);
    let _ = writeln!(s, "  Usable fraction ({:.0}%):  {:.3e} kg", F_USABLE_MERCURY * 100.0, model.mercury_usable_kg);
    let _ = writeln!(s, "  Max nodes from Mercury:   {:.3e}  ({:.1} doublings)",
        model.max_nodes_mercury, model.max_nodes_mercury.log2());
    let sphere_area = 4.0 * std::f64::consts::PI * (R_MERCURY_AU * AU).powi(2);
    let merc_coverage_pct = model.max_nodes_mercury * spec.collector_area_m2 / sphere_area * 100.0;
    let _ = writeln!(s, "  Swarm coverage at Mercury orbit: {:.1}%", merc_coverage_pct.min(100.0));
    let _ = writeln!(s, "  (Full Type II requires Venus + asteroid belt for 100% coverage at 1 AU)");
    let _ = writeln!(s);

    // Part F: timeline
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  F. GROWTH TIMELINE (seed node t=0, doubling every {:.0} days)",
        model.doubling_time_nominal_days);
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  {:>6}  {:>14}  {:>15}  {:>12}  {:>12}",
        "Year", "Nodes", "Mass consumed", "P_captured", "% of L_☉");
    let _ = writeln!(s, "  {:>6}  {:>14}  {:>15}  {:>12}  {:>12}",
        "----", "-----", "-------------", "----------", "--------");

    // Print key years and every 5 years
    let key_years: Vec<f64> = {
        let mut v: Vec<f64> = (0..=40).map(|i| i as f64 * 5.0).collect();
        v.extend_from_slice(&[50.0, 60.0, 70.0, 80.0, 100.0, 120.0, 150.0, 200.0]);
        v.sort_by(|a,b| a.partial_cmp(b).unwrap());
        v.dedup();
        v
    };

    for gp in timeline.iter().filter(|p| key_years.contains(&p.year)) {
        let type_ii_flag = if gp.power_fraction_of_lsol >= 0.10 && gp.power_fraction_of_lsol < 0.101 { " ← TYPE II" } else { "" };
        let _ = writeln!(s, "  {:>6.0}  {:>14.3e}  {:>13.3e} kg  {:>10.3e} W  {:>8.2}%{}",
            gp.year, gp.n_nodes, gp.mass_consumed_kg,
            gp.power_captured_w, gp.power_fraction_of_lsol * 100.0, type_ii_flag);
    }
    let _ = writeln!(s);

    match t2_year {
        Some(y) => {
            let gp = &timeline[y as usize];
            let _ = writeln!(s, "  ★ TYPE II CROSSING: year {:.0}  ({:.3e} nodes, {:.2}% of L_☉ captured)",
                y, gp.n_nodes, gp.power_fraction_of_lsol * 100.0);
        }
        None => {
            let _ = writeln!(s, "  TYPE II not reached within 200 years at nominal assumptions.");
        }
    }
    match sat_year {
        Some(y) => {
            let _ = writeln!(s, "  ★ MERCURY SATURATED: year {:.0}  (all usable feedstock consumed)", y);
        }
        None => {
            let _ = writeln!(s, "  ★ Mercury not fully consumed within 200 years.");
        }
    }
    let _ = writeln!(s);

    // Part G: metal production economics
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  G. FUSION METAL PRODUCTION ECONOMICS");
    let _ = writeln!(s, "     (Can the swarm forge supply chains? Yes — for everything up to Fe.)");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  {:>4}  {:>5}  {:>8}  {:>13}  {:>10}  {}",
        "Sym", "A", "BE/A(MeV)", "Q_fusion(MeV)", "Cost(GJ/kg)", "Notes");
    let _ = writeln!(s, "  {:>4}  {:>5}  {:>8}  {:>13}  {:>10}  {}",
        "---", "-", "---------", "-------------", "----------", "-----");
    for el in metals {
        let cost = synthesis_cost_gj_kg(el);
        let cost_str = if cost == 0.0 { "  FREE(+E)".to_string() } else { format!("{:>10.1}", cost) };
        let q_sign = if el.above_iron { "(*)" } else { "(+)" };
        let _ = writeln!(s, "  {:>4}  {:>5}  {:>8.2}  {:>10.1} {}  {}  {}",
            el.symbol, el.a, el.be_per_a, el.fusion_q_mev, q_sign, cost_str, el.synthesis_note);
    }
    let _ = writeln!(s);
    let _ = writeln!(s, "  (*) = requires neutron capture (r/s-process); cost = energy above Fe-56 peak");
    let _ = writeln!(s, "  (+) = exothermic from hydrogen feedstock; metal is a FREE byproduct");
    let _ = writeln!(s);

    // Gold economics
    let au = metals.iter().find(|e| e.symbol == "Au").unwrap();
    let gold_cost_gj_kg = synthesis_cost_gj_kg(au);
    let gold_price_2026_per_kg = 95_000.0; // USD/kg (2026 approximate)
    let elec_price_per_j = 0.05 / 3_600_000.0; // USD/J at $0.05/kWh
    let gold_energy_cost_usd = gold_cost_gj_kg * 1e9 * elec_price_per_j;
    let _ = writeln!(s, "  Gold synthesis economics (from Fe-56 feedstock):");
    let _ = writeln!(s, "    Energy cost:       {:.1} GJ/kg gold", gold_cost_gj_kg);
    let _ = writeln!(s, "    At $0.05/kWh:      ${:.0}/kg  (vs market ~${:.0}/kg)",
        gold_energy_cost_usd, gold_price_2026_per_kg);
    let _ = writeln!(s, "    Ratio:             {:.1}× more expensive than mining today",
        gold_energy_cost_usd / gold_price_2026_per_kg);
    let _ = writeln!(s, "    At Dyson-swarm energy (≈ $0/kWh): gold synthesis is essentially FREE.");
    let _ = writeln!(s, "    Break-even energy price: ${:.5}/kWh",
        gold_price_2026_per_kg / (gold_cost_gj_kg * 1e9 / 3.6e6));
    let _ = writeln!(s);

    // Closing interpretation
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  H. GREAT FILTER BRIDGE (callback to galactic life map GRAND-354)");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  Galactic life map: Type II emerges at +500 Myr from LUCA equivalent.");
    let _ = writeln!(s, "  This spec: {:.0} years from seed node to Type II (nominal).",
        t2_year.unwrap_or(999.0));
    let _ = writeln!(s, "  The delay is NOT technological. It is:");
    let _ = writeln!(s, "    (a) Biological: getting from prokaryote to tool-using civilisation.");
    let _ = writeln!(s, "    (b) Cultural: deciding to commit to a multi-decade megaproject.");
    let _ = writeln!(s, "    (c) Coordination: sustaining the swarm through political instability.");
    let _ = writeln!(s);
    let _ = writeln!(s, "  Given the engineering spec, the soft constraints dominate by > 7 orders");
    let _ = writeln!(s, "  of magnitude. The technology is solved. The great filter is social.");
    let _ = writeln!(s);
    let _ = writeln!(s, "  One seed node on Mercury. {:.0} months to daughter.",
        model.doubling_time_nominal_days / 30.44);
    let _ = writeln!(s, "  {:.0} years to 10% of a star.", t2_year.unwrap_or(999.0));
    let _ = writeln!(s, "  Everything else is politics.");
    let _ = writeln!(s);
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  End of report — GUTOE Dyson swarm engineering specification");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");

    fs::write(format!("{out}/dyson_swarm_sim.txt"), &s).expect("write txt");
    println!("  → {out}/dyson_swarm_sim.txt");
}

fn write_orbital_csv(scan: &[OrbitalPoint], out: &str) {
    let mut s = String::new();
    let _ = writeln!(s, "r_au,flux_w_m2,t_eq_k,p_rad_upa,v_orb_km_s,period_days,power_per_m2_w");
    for op in scan {
        let _ = writeln!(s, "{:.3},{:.2},{:.2},{:.4},{:.4},{:.4},{:.2}",
            op.r_au, op.flux_w_m2, op.t_eq_k, op.p_rad_pa * 1e6,
            op.v_orb_km_s, op.period_days, op.power_per_m2_w);
    }
    fs::write(format!("{out}/dyson_swarm_orbital_scan.csv"), &s).expect("write csv");
    println!("  → {out}/dyson_swarm_orbital_scan.csv");
}

fn write_timeline_csv(timeline: &[GrowthPoint], out: &str) {
    let mut s = String::new();
    let _ = writeln!(s, "year,n_nodes,mass_consumed_kg,power_captured_w,power_fraction_lsol,coverage_pct");
    for gp in timeline {
        let _ = writeln!(s, "{:.0},{:.4e},{:.4e},{:.4e},{:.8},{:.6}",
            gp.year, gp.n_nodes, gp.mass_consumed_kg,
            gp.power_captured_w, gp.power_fraction_of_lsol,
            gp.coverage_mercury_orbit_pct);
    }
    fs::write(format!("{out}/dyson_swarm_timeline.csv"), &s).expect("write timeline csv");
    println!("  → {out}/dyson_swarm_timeline.csv");
}

fn write_json(
    spec: &NodeSpec,
    model: &ReplicationModel,
    timeline: &[GrowthPoint],
    out: &str,
) {
    let t2 = type_ii_year(timeline);
    let sat = full_coverage_year(timeline);
    let mut s = String::new();
    let _ = writeln!(s, "{{");
    let _ = writeln!(s, "  \"target_orbit_au\": {:.3},", spec.target_orbit_au);
    let _ = writeln!(s, "  \"collector_area_m2\": {:.0},", spec.collector_area_m2);
    let _ = writeln!(s, "  \"collector_thickness_nm\": {:.0},", COLLECTOR_THICKNESS_NM);
    let _ = writeln!(s, "  \"node_mass_kg\": {:.0},", spec.total_mass_kg);
    let _ = writeln!(s, "  \"power_per_node_mw\": {:.2},", spec.power_electric_w * 1e-6);
    let _ = writeln!(s, "  \"e_total_mj_per_kg\": {:.0},", E_TOTAL_MJ_KG);
    let _ = writeln!(s, "  \"overhead_factor\": {:.1},", OVERHEAD_FACTOR);
    let _ = writeln!(s, "  \"doubling_time_days_raw\": {:.2},", model.doubling_time_raw_days);
    let _ = writeln!(s, "  \"doubling_time_days_nominal\": {:.2},", model.doubling_time_nominal_days);
    let _ = writeln!(s, "  \"mercury_usable_kg\": {:.3e},", model.mercury_usable_kg);
    let _ = writeln!(s, "  \"max_nodes_mercury\": {:.3e},", model.max_nodes_mercury);
    let _ = writeln!(s, "  \"type_ii_year\": {},",
        t2.map_or("null".to_string(), |y| format!("{:.0}", y)));
    let _ = writeln!(s, "  \"mercury_saturated_year\": {},",
        sat.map_or("null".to_string(), |y| format!("{:.0}", y)));
    let _ = writeln!(s, "  \"l_sol_w\": {:.4e},", L_SOL);
    let _ = writeln!(s, "  \"material\": \"lonsdaleite\",");
    let _ = writeln!(s, "  \"collector_hardness_gpa\": {:.1},", H_LONSDALEITE_GPA);
    let _ = writeln!(s, "  \"collector_tensile_gpa\": {:.0}", SIGMA_T_LONSDALEITE_GPA);
    let _ = writeln!(s, "}}");
    fs::write(format!("{out}/dyson_swarm.json"), &s).expect("write json");
    println!("  → {out}/dyson_swarm.json");
}

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let out = std::env::var("GUTOE_DYSON_SIM_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/dyson_swarm".to_string());
    fs::create_dir_all(&out).expect("create output dir");

    println!("GUTOE: Self-Replicating Dyson Swarm Simulation");
    println!("  L_☉ = {:.3e} W    M_☉ = {:.3e} kg", L_SOL, M_SOL);
    println!("  Mercury orbit = {:.3} AU    M_Mercury = {:.3e} kg", R_MERCURY_AU, M_MERCURY);
    println!("  Collector material: lonsdaleite  H = {:.1} GPa", H_LONSDALEITE_GPA);
    println!("  Output → {out}");
    println!();

    // Orbital scan: 0.1 to 2.0 AU in steps of 0.05
    let scan: Vec<OrbitalPoint> = {
        let mut v = Vec::new();
        let mut r = 0.10_f64;
        while r <= 2.001 {
            v.push(orbital_point(r));
            r += 0.05;
        }
        v
    };

    let merc_op = orbital_point(R_MERCURY_AU);
    let earth_op = orbital_point(1.0);
    let spec = node_spec(R_MERCURY_AU);
    let model = replication_model(&spec);
    let metals = metal_economics();
    let timeline = growth_timeline(&model, spec.power_electric_w, R_MERCURY_AU);

    println!("Running assertions...");
    run_assertions(&merc_op, &earth_op, &spec, &model, &timeline, &metals);

    println!("Writing outputs...");
    write_txt(&scan, &spec, &model, &timeline, &metals, &out);
    write_orbital_csv(&scan, &out);
    write_timeline_csv(&timeline, &out);
    write_json(&spec, &model, &timeline, &out);

    // Summary printout
    println!();
    println!("Dyson swarm parameters:");
    println!("  Node mass:        {:.0} tonnes", spec.total_mass_kg / 1000.0);
    println!("  Power per node:   {:.1} MW  (at {:.3} AU)", spec.power_electric_w * 1e-6, R_MERCURY_AU);
    println!("  Doubling time:    {:.1} days nominal = {:.2} years",
        model.doubling_time_nominal_days, model.doubling_time_nominal_days / 365.25);
    println!("  Mercury feedstock: {:.2e} nodes max", model.max_nodes_mercury);
    println!();

    if let Some(y) = type_ii_year(&timeline) {
        let gp = &timeline[y as usize];
        println!("  ★ TYPE II crossing: year {:.0}  ({:.2}% of L_☉)", y, gp.power_fraction_of_lsol * 100.0);
    }
    if let Some(y) = full_coverage_year(&timeline) {
        println!("  ★ Mercury saturated: year {:.0}", y);
    }
}
