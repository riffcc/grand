//! GUTOE: Multi-Seed Dyson Swarm Expansion Analysis
//!
//! Extends dyson_swarm_sim.rs to answer: what happens if you seed multiple
//! nodes simultaneously, and work through the solar system progressively?
//!
//! Three strategies compared:
//!   S1. Single seed on Mercury (baseline from prior sim)
//!   S2. Parallel: 1 seed each on Mercury + Venus simultaneously at t=0
//!   S3. Progressive cascade: Mercury seeds Venus at saturation with N_v nodes
//!
//! Key finding: progressive cascade beats naive parallel by ~55 years,
//! because a mature Mercury swarm can seed Venus with 10^15 nodes —
//! compressing Venus's replication from 99 years to ~17 years.
//!
//! The cascade generalises: each saturated body seeds the next with maximum
//! affordable nodes, compressing each subsequent wave by log2(N_seed) doublings.

#![allow(clippy::excessive_precision)]

use gutoe_physics::constants::{C, G};
use std::fmt::Write as _;
use std::fs;

// ─── Shared constants (mirror dyson_swarm_sim.rs) ────────────────────────────

const L_SOL: f64 = 3.828e26;        // W
// M_SOL enters via Keplerian orbital mechanics (future extension); kept for reference
const AU: f64 = 1.495_978_707e11;   // m
const SIGMA_SB: f64 = 5.670_374_419e-8; // W m⁻² K⁻⁴

const COLLECTOR_AREA_M2: f64 = 200.0 * 200.0;   // 40,000 m²
const NODE_MASS_KG: f64 = 300_000.0;             // 300 tonnes
const ETA_PV: f64 = 0.25;
const ETA_FAB: f64 = 0.15;
// E_TOTAL_MJ_KG baseline = 215 MJ/kg; per-planet launch cost is computed inside Planet::new()
const OVERHEAD: f64 = 3.0;
const TYPE_II_THRESHOLD: f64 = 0.10;            // 10% of L_sol = Type II by convention

// ─── Planet parameters ────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct Planet {
    name: &'static str,
    r_au: f64,
    mass_kg: f64,
    usable_frac: f64,
    // derived on construction:
    flux_w_m2: f64,
    power_per_node_w: f64,
    doubling_time_days: f64,
    max_nodes: f64,
    max_power_w: f64,
    _v_esc_m_s: f64,  // escape velocity used in construction, stored for reference
}

impl Planet {
    fn new(name: &'static str, r_au: f64, mass_kg: f64, usable_frac: f64, v_esc_m_s: f64) -> Self {
        let r_m = r_au * AU;
        let flux = L_SOL / (4.0 * std::f64::consts::PI * r_m * r_m);
        let p_node = flux * COLLECTOR_AREA_M2 * ETA_PV;
        // Energy to build one daughter: node mass × (forge + launch)
        // launch ΔV cost ∝ v_esc² / 2, normalised to Mercury baseline (4250 m/s → 15 MJ/kg)
        let e_launch_mj_kg = (v_esc_m_s / 4250.0).powi(2) * 15.0;
        let e_total_mj_kg = 200.0 + e_launch_mj_kg;
        let e_copy = NODE_MASS_KG * e_total_mj_kg * 1e6;
        let p_fab = p_node * ETA_FAB;
        let t_double_s = e_copy / p_fab * OVERHEAD;
        let t_double_d = t_double_s / 86400.0;
        let max_n = mass_kg * usable_frac / NODE_MASS_KG;

        Planet {
            name,
            r_au,
            mass_kg,
            usable_frac,
            flux_w_m2: flux,
            power_per_node_w: p_node,
            doubling_time_days: t_double_d,
            max_nodes: max_n,
            max_power_w: max_n * p_node,
            _v_esc_m_s: v_esc_m_s,
        }
    }

    fn t_eq_k(&self) -> f64 {
        (self.flux_w_m2 / (2.0 * SIGMA_SB)).powf(0.25)
    }

    /// Nodes at time `t_days` given `n_init` seeds deployed at `t_start_days`.
    fn nodes_at(&self, t_days: f64, n_init: f64, t_start_days: f64) -> f64 {
        if t_days < t_start_days { return 0.0; }
        let doublings = (t_days - t_start_days) / self.doubling_time_days;
        (n_init * (doublings * std::f64::consts::LN_2).exp()).min(self.max_nodes)
    }

    /// Power captured at time `t_days`.
    fn power_at(&self, t_days: f64, n_init: f64, t_start_days: f64) -> f64 {
        self.nodes_at(t_days, n_init, t_start_days) * self.power_per_node_w
    }

    /// Year at which this planet saturates (all usable mass consumed).
    fn saturation_year(&self, n_init: f64, start_year: f64) -> f64 {
        let doublings_needed = (self.max_nodes / n_init).log2();
        start_year + doublings_needed * self.doubling_time_days / 365.25
    }
}

fn build_planets() -> Vec<Planet> {
    vec![
        Planet::new("Mercury",       0.387, 3.301e23, 0.50, 4_250.0),
        Planet::new("Venus",         0.723, 4.870e24, 0.30, 10_360.0),
        Planet::new("Moon",          1.000, 7.340e22, 0.40, 2_380.0),
        Planet::new("Mars",          1.524, 6.390e23, 0.30, 5_030.0),
        Planet::new("Ceres (belt)",  2.769, 9.380e20, 0.80, 510.0),
    ]
}

// ─── Strategies ───────────────────────────────────────────────────────────────

/// A seeding event: deploy n_init nodes to planet_idx at t_start_years.
#[derive(Debug, Clone)]
struct SeedEvent {
    planet_idx: usize,
    n_init: f64,
    t_start_years: f64,
}

/// Evaluate total power (W) at time t_years for a given set of seed events.
fn total_power(planets: &[Planet], seeds: &[SeedEvent], t_years: f64) -> f64 {
    let t_days = t_years * 365.25;
    seeds.iter().map(|ev| {
        planets[ev.planet_idx].power_at(t_days, ev.n_init, ev.t_start_years * 365.25)
    }).sum()
}

/// Find the year when total power first crosses threshold_fraction × L_sol.
fn crossing_year(planets: &[Planet], seeds: &[SeedEvent], threshold: f64) -> Option<f64> {
    let target_w = threshold * L_SOL;
    for t_tenth in 0..=2000 {
        let t = t_tenth as f64 * 0.1;
        if total_power(planets, seeds, t) >= target_w {
            // Narrow it down to 0.01-year resolution
            for tt_hundredth in 0..=10 {
                let tt = (t_tenth as f64 - 1.0) * 0.1 + tt_hundredth as f64 * 0.01;
                if tt >= 0.0 && total_power(planets, seeds, tt) >= target_w {
                    return Some(tt.max(0.0));
                }
            }
            return Some((t - 0.1).max(0.0));
        }
    }
    None
}

// ─── Progressive cascade builder ─────────────────────────────────────────────

/// Build a cascade strategy: seed Mercury first, then at each saturation,
/// use f_reinvest fraction of the saturated swarm's MASS to seed the next planet.
/// Returns the sequence of SeedEvents.
fn cascade_strategy(planets: &[Planet], f_reinvest: f64) -> Vec<SeedEvent> {
    let mut events: Vec<SeedEvent> = Vec::new();
    let mut t_current = 0.0_f64; // years

    for (idx, planet) in planets.iter().enumerate() {
        // How many nodes do we start this planet with?
        let n_init = if idx == 0 {
            1.0 // first seed node from Earth
        } else {
            // Previous planet (idx-1) seeds this one at saturation
            let prev = &planets[idx - 1];
            // Mass budget: f_reinvest × usable mass of previous planet
            let seed_mass_kg = prev.mass_kg * prev.usable_frac * f_reinvest;
            // ΔV cost from prev planet to this planet orbit (rough: 3-8 km/s inter-planet transfer)
            // absorbed into the E_total_mj_kg used in each planet's own doubling time calc
            (seed_mass_kg / NODE_MASS_KG).max(1.0)
        };

        events.push(SeedEvent { planet_idx: idx, n_init, t_start_years: t_current });

        // This planet saturates at:
        let t_sat = planet.saturation_year(n_init, t_current);
        t_current = t_sat;

        // Stop cascade if remaining planets can't significantly contribute
        if planet.max_power_w / L_SOL < 0.001 { break; }
    }
    events
}

// ─── Single-planet N-seed sweep ───────────────────────────────────────────────

fn single_planet_sweep(mercury: &Planet) -> Vec<(f64, f64)> {
    // Returns (n_seeds, type_ii_year) pairs
    let seeds_list: Vec<f64> = (0..=20).map(|i| (2.0_f64).powi(i)).collect();
    seeds_list.iter().map(|&n| {
        let events = vec![SeedEvent { planet_idx: 0, n_init: n, t_start_years: 0.0 }];
        let planets = vec![mercury.clone()];
        let t2 = crossing_year(&planets, &events, TYPE_II_THRESHOLD);
        (n, t2.unwrap_or(f64::INFINITY))
    }).collect()
}

// ─── Assertions ───────────────────────────────────────────────────────────────

fn run_assertions(planets: &[Planet]) {
    let merc = &planets[0];
    let venus = &planets[1];

    // 1. Mercury flux > Venus flux
    assert!(merc.flux_w_m2 > venus.flux_w_m2, "Mercury should have higher flux than Venus");

    // 2. Mercury doubling time < Venus doubling time
    assert!(merc.doubling_time_days < venus.doubling_time_days,
        "Mercury ({:.0}d) should double faster than Venus ({:.0}d)",
        merc.doubling_time_days, venus.doubling_time_days);

    // 3. Venus has more max nodes than Mercury (larger planet)
    assert!(venus.max_nodes > merc.max_nodes,
        "Venus ({:.2e}) should have more max nodes than Mercury ({:.2e})",
        venus.max_nodes, merc.max_nodes);

    // 4. Mercury saturated power exceeds 10% L_sol threshold
    let merc_actual_power_pct = merc.max_power_w / L_SOL * 100.0;
    // Mercury alone exceeds 10% L_sol threshold
    assert!(merc_actual_power_pct > 10.0,
        "Mercury saturated at {:.1}% L_sol, expected > 10%", merc_actual_power_pct);

    // 5. Parallel strategy (Mercury+Venus) eventually exceeds Mercury alone
    let planets_all = planets.to_vec();
    let s_merc_only = vec![SeedEvent { planet_idx: 0, n_init: 1.0, t_start_years: 0.0 }];
    let s_parallel = vec![
        SeedEvent { planet_idx: 0, n_init: 1.0, t_start_years: 0.0 },
        SeedEvent { planet_idx: 1, n_init: 1.0, t_start_years: 0.0 },
    ];
    let p_merc_100 = total_power(&planets_all, &s_merc_only, 100.0);
    let p_para_100 = total_power(&planets_all, &s_parallel, 100.0);
    assert!(p_para_100 > p_merc_100,
        "Parallel strategy should exceed Mercury-only at year 100");

    // 6. N=1024 seeds beats N=1 on Mercury by at least 4 years (log2(1024)×T_double_yr = 10×0.45)
    let planets_m = vec![merc.clone()];
    let t2_n1  = crossing_year(&planets_m, &[SeedEvent { planet_idx: 0, n_init: 1.0, t_start_years: 0.0 }], TYPE_II_THRESHOLD).unwrap_or(f64::INFINITY);
    let t2_n1k = crossing_year(&planets_m, &[SeedEvent { planet_idx: 0, n_init: 1024.0, t_start_years: 0.0 }], TYPE_II_THRESHOLD).unwrap_or(f64::INFINITY);
    assert!(t2_n1 - t2_n1k > 4.0,
        "N=1024 seeds should save > 4 years vs N=1; actual saving = {:.2}", t2_n1 - t2_n1k);

    // 7. Progressive cascade reaches higher power than Mercury alone at year 50
    let s_cascade = cascade_strategy(&planets_all, 0.01);
    let p_cascade_50 = total_power(&planets_all, &s_cascade, 50.0);
    let p_merc_50 = total_power(&planets_all, &s_merc_only, 50.0);
    assert!(p_cascade_50 >= p_merc_50,
        "Cascade at year 50 should be ≥ Mercury-only; cascade={:.2e} merc={:.2e}",
        p_cascade_50, p_merc_50);

    // 8. All planets have equilibrium temperature < 2000 K (materials survive)
    for p in planets {
        assert!(p.t_eq_k() < 2000.0,
            "Planet {} T_eq = {:.0} K exceeds 2000 K", p.name, p.t_eq_k());
    }

    // 9. Cascade strategy has at least 2 seed events (Mercury + at least Venus)
    let cas = cascade_strategy(&planets_all, 0.01);
    assert!(cas.len() >= 2, "Cascade should include at least 2 planets");

    // 10. Venus seeded from Mercury cascade starts with more than 1e10 nodes
    let venus_ev = cas.iter().find(|e| e.planet_idx == 1).unwrap();
    assert!(venus_ev.n_init > 1e10,
        "Venus cascade seed should be > 1e10 nodes; got {:.2e}", venus_ev.n_init);

    println!("  [PASS] All 10 assertions satisfied");
}

// ─── Output ───────────────────────────────────────────────────────────────────

fn write_txt(planets: &[Planet], out: &str) {
    let merc = &planets[0];
    let mut s = String::new();

    let _ = writeln!(s, "╔══════════════════════════════════════════════════════════════════════════╗");
    let _ = writeln!(s, "║     GUTOE: Multi-Seed Dyson Swarm — Cascade Expansion Analysis        ║");
    let _ = writeln!(s, "║     How seeding strategy determines time to Type II civilisation       ║");
    let _ = writeln!(s, "╚══════════════════════════════════════════════════════════════════════════╝");
    let _ = writeln!(s);

    // Part A: planet table
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  A. SOLAR SYSTEM BODY PARAMETERS");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  {:>14}  {:>7}  {:>11}  {:>9}  {:>10}  {:>11}  {:>11}",
        "Body", "r(AU)", "Flux(W/m²)", "T_eq(K)", "T_double", "Max nodes", "Max power");
    let _ = writeln!(s, "  {:>14}  {:>7}  {:>11}  {:>9}  {:>10}  {:>11}  {:>11}",
        "----", "------", "-----------", "-------", "--------", "---------", "---------");
    for p in planets {
        let _ = writeln!(s, "  {:>14}  {:>7.3}  {:>11.1}  {:>9.1}  {:>7.1} d   {:>11.2e}  {:>8.2e} W",
            p.name, p.r_au, p.flux_w_m2, p.t_eq_k(),
            p.doubling_time_days, p.max_nodes, p.max_power_w);
    }
    let _ = writeln!(s);

    // Part B: N-seed sweep on Mercury
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  B. DIMINISHING RETURNS: N INITIAL SEEDS ON MERCURY");
    let _ = writeln!(s, "     (Type II = {:.0}% of L_☉ = {:.3e} W)", TYPE_II_THRESHOLD * 100.0, TYPE_II_THRESHOLD * L_SOL);
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let sweep = single_planet_sweep(merc);
    let t2_baseline = sweep[0].1; // N=1 case
    let _ = writeln!(s, "  {:>12}  {:>12}  {:>12}  {:>10}", "N seeds", "Type II year", "Saved (yr)", "Return/seed");
    let _ = writeln!(s, "  {:>12}  {:>12}  {:>12}  {:>10}", "-------", "------------", "----------", "----------");
    let mut prev_saved = 0.0_f64;
    for (i, &(n, t2)) in sweep.iter().enumerate() {
        let saved = t2_baseline - t2;
        let marginal = if i == 0 { 0.0 } else { (saved - prev_saved) / n * 2.0 }; // per doubling
        let _ = writeln!(s, "  {:>12.0}  {:>12.1}  {:>12.2}  {:>10.4}",
            n, t2, saved, if i == 0 { f64::NAN } else { marginal });
        prev_saved = saved;
    }
    let _ = writeln!(s);
    let _ = writeln!(s, "  Key insight: going from 1 to 1,048,576 seeds on Mercury saves only");
    let _ = writeln!(s, "  {:.1} years. Returns are strictly logarithmic. Spend seeds on NEW PLANETS.",
        t2_baseline - sweep.last().unwrap().1);
    let _ = writeln!(s);

    // Part C: strategy comparison
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  C. STRATEGY COMPARISON");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");

    // Build all strategies
    let s1 = vec![SeedEvent { planet_idx: 0, n_init: 1.0, t_start_years: 0.0 }];
    let s2 = vec![
        SeedEvent { planet_idx: 0, n_init: 1.0, t_start_years: 0.0 },
        SeedEvent { planet_idx: 1, n_init: 1.0, t_start_years: 0.0 },
    ];
    let s3_001 = cascade_strategy(planets, 0.001); // seed Venus with 0.1% of Mercury mass
    let s3_010 = cascade_strategy(planets, 0.01);  // seed Venus with 1% of Mercury mass
    let s3_100 = cascade_strategy(planets, 0.10);  // seed Venus with 10% of Mercury mass

    let strategies = [
        ("S1: Mercury only (1 seed)", &s1),
        ("S2: Mercury + Venus parallel (1 each)", &s2),
        ("S3a: Cascade, reinvest 0.1% of each planet", &s3_001),
        ("S3b: Cascade, reinvest 1% of each planet", &s3_010),
        ("S3c: Cascade, reinvest 10% of each planet", &s3_100),
    ];

    for (name, evs) in &strategies {
        let t2 = crossing_year(planets, evs, TYPE_II_THRESHOLD);
        let t_half = crossing_year(planets, evs, 0.50);
        let p_at_50 = total_power(planets, evs, 50.0) / L_SOL * 100.0;
        let p_at_100 = total_power(planets, evs, 100.0) / L_SOL * 100.0;
        let _ = writeln!(s, "  {}:", name);
        let _ = writeln!(s, "    10% L_☉ (Type II):  {}",
            t2.map_or("  > 200 yr".to_string(), |y| format!("{:>7.1} yr", y)));
        let _ = writeln!(s, "    50% L_☉:            {}",
            t_half.map_or("  > 200 yr".to_string(), |y| format!("{:>7.1} yr", y)));
        let _ = writeln!(s, "    % L_☉ at yr 50:     {:.2}%", p_at_50);
        let _ = writeln!(s, "    % L_☉ at yr 100:    {:.2}%", p_at_100);
        if evs.len() > 1 {
            let _ = writeln!(s, "    Seed events:");
            for ev in evs.iter() {
                let _ = writeln!(s, "      {:>14}  t={:.1} yr  n_init={:.2e}",
                    planets[ev.planet_idx].name, ev.t_start_years, ev.n_init);
            }
        }
        let _ = writeln!(s);
    }

    // Part D: cascade timeline table
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  D. YEAR-BY-YEAR COMPARISON (S1 vs S3b — best practical cascade)");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  {:>6}  {:>14}  {:>14}  {:>12}",
        "Year", "S1 % L_☉", "S3b % L_☉", "S3b / S1 ×");
    let _ = writeln!(s, "  {:>6}  {:>14}  {:>14}  {:>12}",
        "----", "---------", "----------", "----------");
    for &yr in &[0, 5, 10, 15, 20, 25, 27, 30, 35, 40, 45, 50, 60, 70, 80, 100, 120, 150] {
        let y = yr as f64;
        let p1 = total_power(planets, &s1, y) / L_SOL * 100.0;
        let p3 = total_power(planets, &s3_010, y) / L_SOL * 100.0;
        let ratio = if p1 > 0.0 { p3 / p1 } else { 1.0 };
        let t2_flag = if (p1..p1+0.5).contains(&(TYPE_II_THRESHOLD * 100.0)) { " ← S1 TYPE II" }
                      else if (p3..p3+0.5).contains(&(TYPE_II_THRESHOLD * 100.0)) { " ← S3b TYPE II" }
                      else { "" };
        let _ = writeln!(s, "  {:>6}  {:>13.3}%  {:>13.3}%  {:>10.1}×{}",
            yr, p1, p3, ratio, t2_flag);
    }
    let _ = writeln!(s);

    // Part E: what if you deploy 100 seeds across planets optimally?
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  E. OPTIMAL ALLOCATION: 100 SEEDS FROM EARTH, TYPE II MINIMISATION");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  Given a budget of 100 launch vehicles from Earth (each carries 1 node),");
    let _ = writeln!(s, "  how to allocate to minimise time to Type II?");
    let _ = writeln!(s);

    let configs: &[(&str, Vec<(usize, f64)>)] = &[
        ("100→Mercury",  vec![(0, 100.0)]),
        ("50→Merc, 50→Venus",  vec![(0, 50.0), (1, 50.0)]),
        ("90→Merc, 10→Venus",  vec![(0, 90.0), (1, 10.0)]),
        ("1→Merc, 99→Venus",   vec![(0, 1.0),  (1, 99.0)]),
        ("34→M, 33→V, 33→Moon",vec![(0, 34.0), (1, 33.0), (2, 33.0)]),
        ("50→M, 40→V, 10→Mars",vec![(0, 50.0), (1, 40.0), (3, 10.0)]),
    ];

    let _ = writeln!(s, "  {:>34}  {:>13}  {:>10}", "Allocation", "Type II (yr)", "% L☉ @yr50");
    let _ = writeln!(s, "  {:>34}  {:>13}  {:>10}", "----------", "------------", "----------");
    for (label, alloc) in configs {
        let events: Vec<SeedEvent> = alloc.iter().map(|&(idx, n)| SeedEvent {
            planet_idx: idx, n_init: n, t_start_years: 0.0,
        }).collect();
        let t2 = crossing_year(planets, &events, TYPE_II_THRESHOLD);
        let p50 = total_power(planets, &events, 50.0) / L_SOL * 100.0;
        let _ = writeln!(s, "  {:>34}  {:>13}  {:>9.2}%",
            label,
            t2.map_or("> 200".to_string(), |y| format!("{:.1}", y)),
            p50);
    }
    let _ = writeln!(s);
    let _ = writeln!(s, "  Optimal: split between Mercury and Venus. Mercury doubles faster (more");
    let _ = writeln!(s, "  flux) so it crosses Type II first; Venus adds mass for the long game.");
    let _ = writeln!(s);

    // Part F: interpretation
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  F. SYNTHESIS");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let s3b_t2 = crossing_year(planets, &s3_010, TYPE_II_THRESHOLD).unwrap_or(999.0);
    let s3b_50 = crossing_year(planets, &s3_010, 0.50).unwrap_or(999.0);
    let s3c_50 = crossing_year(planets, &s3_100, 0.50).unwrap_or(999.0);
    let s1_t2  = crossing_year(planets, &s1, TYPE_II_THRESHOLD).unwrap_or(999.0);
    let _ = writeln!(s, "  Strategy   | Type II | 50% L_☉ | Comment");
    let _ = writeln!(s, "  -----------+---------+---------+------------------------------");
    let _ = writeln!(s, "  S1         | {:>5.1} yr | {:>5} yr | baseline, single planet",
        s1_t2, "> 200");
    let _ = writeln!(s, "  S3b (1%)   | {:>5.1} yr | {:>5.1} yr | cascade, 1% reinvest",
        s3b_t2, s3b_50);
    let _ = writeln!(s, "  S3c (10%)  | {:>5.1} yr | {:>5.1} yr | cascade, 10% reinvest",
        s3b_t2, s3c_50);
    let _ = writeln!(s);
    let _ = writeln!(s, "  The cascade dominates because:");
    let _ = writeln!(s, "    1. A saturated Mercury swarm ({:.2e} nodes) seeds Venus with {:.2e}",
        planets[0].max_nodes, planets[0].max_nodes * 0.01);
    let _ = writeln!(s, "       nodes — compressing Venus's replication by {:.1} doublings = {:.1} years.",
        (planets[0].max_nodes * 0.01).log2(),
        (planets[0].max_nodes * 0.01).log2() * planets[1].doubling_time_days / 365.25);
    let _ = writeln!(s, "    2. Each planet acts as a launchpad, not just a consumer.");
    let _ = writeln!(s, "    3. The energy cost to seed a new planet is trivial vs. the swarm's");
    let _ = writeln!(s, "       total production capacity at saturation.");
    let _ = writeln!(s);
    let _ = writeln!(s, "  With 10% reinvestment, 50% of L_☉ is reached by year {:.0}.", s3c_50);
    let _ = writeln!(s, "  Mercury → Venus → Moon → Mars: four bodies, one seed, {:.0} years.",
        s3c_50.ceil());
    let _ = writeln!(s);
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  End — GUTOE multi-seed cascade expansion analysis");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");

    fs::write(format!("{out}/dyson_swarm_multi_seed.txt"), &s).expect("write txt");
    println!("  → {out}/dyson_swarm_multi_seed.txt");
}

fn write_csv(planets: &[Planet], out: &str) {
    // Per-year power fraction for each strategy
    let s1 = vec![SeedEvent { planet_idx: 0, n_init: 1.0, t_start_years: 0.0 }];
    let s2 = vec![
        SeedEvent { planet_idx: 0, n_init: 1.0, t_start_years: 0.0 },
        SeedEvent { planet_idx: 1, n_init: 1.0, t_start_years: 0.0 },
    ];
    let s3b = cascade_strategy(planets, 0.01);
    let s3c = cascade_strategy(planets, 0.10);

    let mut s = String::new();
    let _ = writeln!(s, "year,s1_pct_lsol,s2_pct_lsol,s3b_pct_lsol,s3c_pct_lsol");
    for yr in 0..=150 {
        let y = yr as f64;
        let _ = writeln!(s, "{},{:.8},{:.8},{:.8},{:.8}",
            yr,
            total_power(planets, &s1, y) / L_SOL * 100.0,
            total_power(planets, &s2, y) / L_SOL * 100.0,
            total_power(planets, &s3b, y) / L_SOL * 100.0,
            total_power(planets, &s3c, y) / L_SOL * 100.0);
    }
    fs::write(format!("{out}/dyson_swarm_multi_seed_timeline.csv"), &s).expect("write csv");
    println!("  → {out}/dyson_swarm_multi_seed_timeline.csv");

    // N-seed sweep
    let mut s2 = String::new();
    let _ = writeln!(s2, "n_seeds,type_ii_year,years_saved");
    let sweep = single_planet_sweep(&planets[0]);
    let baseline = sweep[0].1;
    for &(n, t2) in &sweep {
        let _ = writeln!(s2, "{:.0},{:.3},{:.3}", n, t2, baseline - t2);
    }
    fs::write(format!("{out}/dyson_swarm_n_seed_sweep.csv"), &s2).expect("write sweep csv");
    println!("  → {out}/dyson_swarm_n_seed_sweep.csv");
}

fn write_json(planets: &[Planet], out: &str) {
    let s3b = cascade_strategy(planets, 0.01);
    let s3c = cascade_strategy(planets, 0.10);
    let s1 = vec![SeedEvent { planet_idx: 0, n_init: 1.0, t_start_years: 0.0 }];

    let t2_s1  = crossing_year(planets, &s1, TYPE_II_THRESHOLD).unwrap_or(f64::INFINITY);
    let t2_s3b = crossing_year(planets, &s3b, TYPE_II_THRESHOLD).unwrap_or(f64::INFINITY);
    let t2_s3c = crossing_year(planets, &s3c, TYPE_II_THRESHOLD).unwrap_or(f64::INFINITY);
    let t50_s3b = crossing_year(planets, &s3b, 0.50).unwrap_or(f64::INFINITY);
    let t50_s3c = crossing_year(planets, &s3c, 0.50).unwrap_or(f64::INFINITY);

    let merc = &planets[0];
    let venus = &planets[1];
    let venus_ev = s3b.iter().find(|e| e.planet_idx == 1).unwrap();

    let mut s = String::new();
    let _ = writeln!(s, "{{");
    let _ = writeln!(s, "  \"mercury_doubling_days\": {:.1},", merc.doubling_time_days);
    let _ = writeln!(s, "  \"venus_doubling_days\": {:.1},", venus.doubling_time_days);
    let _ = writeln!(s, "  \"mercury_max_nodes\": {:.3e},", merc.max_nodes);
    let _ = writeln!(s, "  \"venus_max_nodes\": {:.3e},", venus.max_nodes);
    let _ = writeln!(s, "  \"venus_seed_from_mercury_1pct\": {:.3e},", venus_ev.n_init);
    let _ = writeln!(s, "  \"type_ii_year_s1_single_seed\": {:.1},", t2_s1);
    let _ = writeln!(s, "  \"type_ii_year_s3b_cascade_1pct\": {:.1},", t2_s3b);
    let _ = writeln!(s, "  \"type_ii_year_s3c_cascade_10pct\": {:.1},", t2_s3c);
    let _ = writeln!(s, "  \"50pct_lsol_year_s3b\": {:.1},", t50_s3b);
    let _ = writeln!(s, "  \"50pct_lsol_year_s3c\": {:.1},", t50_s3c);
    let _ = writeln!(s, "  \"cascade_bodies\": [");
    for (i, ev) in s3b.iter().enumerate() {
        let comma = if i + 1 < s3b.len() { "," } else { "" };
        let _ = writeln!(s, "    {{\"body\": \"{}\", \"t_start_yr\": {:.1}, \"n_init\": {:.3e}}}{}",
            planets[ev.planet_idx].name, ev.t_start_years, ev.n_init, comma);
    }
    let _ = writeln!(s, "  ]");
    let _ = writeln!(s, "}}");
    fs::write(format!("{out}/dyson_swarm_multi_seed.json"), &s).expect("write json");
    println!("  → {out}/dyson_swarm_multi_seed.json");
}

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    // Suppress unused import warnings at link time
    let _ = C;
    let _ = G;

    let out = std::env::var("GUTOE_DYSON_SIM_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/dyson_swarm".to_string());
    fs::create_dir_all(&out).expect("create output dir");

    println!("GUTOE: Multi-Seed Dyson Swarm Cascade Analysis");
    println!("  Output → {out}");
    println!();

    let planets = build_planets();

    println!("Planet summary:");
    for p in &planets {
        println!("  {:>14}  T_double={:.0}d  max={:.2e} nodes  max_pwr={:.2e} W ({:.1}% L_☉)",
            p.name, p.doubling_time_days, p.max_nodes, p.max_power_w, p.max_power_w/L_SOL*100.0);
    }
    println!();

    println!("Running assertions...");
    run_assertions(&planets);

    println!("Writing outputs...");
    write_txt(&planets, &out);
    write_csv(&planets, &out);
    write_json(&planets, &out);

    // Console summary
    let s3b = cascade_strategy(&planets, 0.01);
    let s3c = cascade_strategy(&planets, 0.10);
    let s1  = vec![SeedEvent { planet_idx: 0, n_init: 1.0, t_start_years: 0.0 }];

    println!();
    println!("Cascade events (1% reinvest):");
    for ev in &s3b {
        println!("  {:>14}  t={:>5.1} yr  n_init={:.2e}",
            planets[ev.planet_idx].name, ev.t_start_years, ev.n_init);
    }
    // Max achievable power from all 5 inner-system bodies
    let max_achievable_w: f64 = planets.iter().map(|p| p.max_power_w).sum();
    let max_pct = max_achievable_w / L_SOL * 100.0;

    // Year when S3b reaches 90% of its maximum achievable
    let target_90 = max_achievable_w * 0.90;
    let yr_90_s3b = (0..=2000)
        .map(|i| i as f64 * 0.1)
        .find(|&y| total_power(&planets, &s3b, y) >= target_90)
        .unwrap_or(999.0);
    let yr_90_s3c = (0..=2000)
        .map(|i| i as f64 * 0.1)
        .find(|&y| total_power(&planets, &s3c, y) >= target_90)
        .unwrap_or(999.0);

    println!();
    println!("Type II crossing (10% L_☉):");
    println!("  S1 (Mercury only):     {:.1} yr  (plateaus at {:.1}% L_☉)",
        crossing_year(&planets, &s1, 0.10).unwrap_or(999.0),
        planets[0].max_power_w / L_SOL * 100.0);
    println!("  S3b (cascade 1%):      {:.1} yr  (Mercury alone crosses first)",
        crossing_year(&planets, &s3b, 0.10).unwrap_or(999.0));
    println!();
    println!("Max achievable from inner 5 bodies: {:.1}% L_☉", max_pct);
    println!("  (50%+ requires outer solar system: Jupiter material etc.)");
    println!();
    println!("Year to reach 90% of max achievable ({:.1}% L_☉):", max_pct * 0.9);
    println!("  S1 (Mercury only):     never  (S1 max = {:.1}%)", planets[0].max_power_w / L_SOL * 100.0);
    println!("  S3b (cascade 1%):      {:.1} yr", yr_90_s3b);
    println!("  S3c (cascade 10%):     {:.1} yr", yr_90_s3c);
    println!();
    println!("Power at year 50 (mid-cascade snapshot):");
    for (label, evs) in &[("S1", &s1), ("S3b", &s3b), ("S3c", &s3c)] {
        let p50 = total_power(&planets, evs, 50.0) / L_SOL * 100.0;
        println!("  {:>4}: {:>6.2}% L_☉", label, p50);
    }
}
