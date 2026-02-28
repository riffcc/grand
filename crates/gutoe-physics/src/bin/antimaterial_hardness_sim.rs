//! GUTOE: Antimaterial Phase Analysis
//!
//! CPT symmetry means antiparticles have identical mass and coupling constants.
//! Therefore:
//!   anti-lonsdaleite is as hard as lonsdaleite     (trivially, by CPT)
//!   anti-diamond     is as hard as diamond          (trivially, by CPT)
//!
//! The NEW quantity antimaterials bring: annihilation energy density.
//!   E_ann/vol = 2 × ρ × c²
//!
//! This does NOT depend on hardness or bonding — only on mass density.
//! A denser antimaterial releases more energy per cm³.
//!
//! Two distinct optima emerge from the GUTOE two-lattice framework:
//!   Matter lattice (N_c=4, tetrahedral sp³):  maximises HARDNESS
//!   Field lattice  (N_c=6, SC):               maximises ENERGY DENSITY
//!
//! For maximum structural antimaterial (survives containment longest):
//!   → anti-lonsdaleite  (H ≈ 105 GPa proxy)
//!
//! For maximum energy yield per volume (drives / warheads):
//!   → SC-anti-carbon within carbon phases, or anti-osmium across all elements
//!
//! Outputs: antimaterial ranking by hardness, density, energy density, and
//!          combined "antimaterial figure of merit".

#![allow(clippy::excessive_precision)]

use gutoe_physics::constants::{ALPHA, C};
use std::fmt::Write as _;
use std::fs;

// ─── Physical constants ───────────────────────────────────────────────────────

const AMU_KG: f64 = 1.660_539_066_60e-27; // 1 atomic mass unit in kg
const A3_TO_M3: f64 = 1e-30; // 1 Å³ in m³
const J_TO_EV: f64 = 6.241_509_074e18; // 1 J in eV

// ─── Annihilation physics ─────────────────────────────────────────────────────

/// Mass-energy per atom pair (matter + antimatter): E = 2 m c²
fn ann_energy_per_atom_j(mass_amu: f64) -> f64 {
    2.0 * mass_amu * AMU_KG * C * C
}

fn ann_energy_per_atom_ev(mass_amu: f64) -> f64 {
    ann_energy_per_atom_j(mass_amu) * J_TO_EV
}

fn ann_energy_per_atom_gev(mass_amu: f64) -> f64 {
    ann_energy_per_atom_ev(mass_amu) * 1e-9
}

/// Density in kg/m³ from atom mass (amu) and volume per atom (Å³)
fn density_kg_m3(mass_amu: f64, v_atom_aa3: f64) -> f64 {
    (mass_amu * AMU_KG) / (v_atom_aa3 * A3_TO_M3)
}

/// Annihilation energy density J/m³ = 2 × ρ × c²
fn ann_energy_density_j_m3(rho_kg_m3: f64) -> f64 {
    2.0 * rho_kg_m3 * C * C
}

/// Annihilation energy density in PJ/cm³ (practical unit for macroscopic amounts)
fn ann_energy_density_pj_cm3(rho_kg_m3: f64) -> f64 {
    ann_energy_density_j_m3(rho_kg_m3) * 1e-15 * 1e-6 // J/m³ → PJ/cm³
}

// ─── QED bond model (reused from carbon_hardness_lattice) ────────────────────

const A0: f64 = 0.529_177_210_8; // Å
const D_SIG_REF: f64 = 1.545;
const E_SIG_REF: f64 = 3.60;
const D_PI_REF: f64 = 1.340;
const E_DBL_REF: f64 = 6.35;

fn e_sigma(d: f64) -> f64 {
    let k_sig = E_SIG_REF * (D_SIG_REF / A0).powi(2);
    k_sig * (A0 / d).powi(2)
}

fn e_pi(d: f64) -> f64 {
    let e_pi_ref = E_DBL_REF - e_sigma(D_PI_REF);
    let k_pi = e_pi_ref * (D_PI_REF / A0).powi(3);
    k_pi * (A0 / d).powi(3)
}

fn cohen_b(d: f64, lambda: f64, nb_eff: f64, f3d: f64) -> f64 {
    let b0 = (1971.0 - 220.0 * lambda) / d.powf(3.5);
    b0 * nb_eff * f3d
}

const H_FACTOR: f64 = 0.23;

// ─── Phase data ───────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct Phase {
    name: &'static str,
    element: &'static str,
    mass_amu: f64,
    nc: f64,
    d_cc: f64,       // bond length (Å); 0.0 if not covalent-bond modelled
    n_sig: f64,      // σ bonds per bond (used in e_bond / e_coh)
    n_pi: f64,       // π bonds per bond (used in e_bond / e_coh)
    v_atom: f64,     // Å³/atom
    lam: f64,
    f3d: f64,
    category: &'static str,
    exp_h_gpa: Option<f64>, // experimental Vickers hardness (GPa)
}

impl Phase {
    fn rho(&self) -> f64 {
        density_kg_m3(self.mass_amu, self.v_atom)
    }

    /// Structural hardness proxy (GPa).  By CPT, anti-phase = phase hardness.
    fn hardness(&self) -> f64 {
        if self.d_cc == 0.0 {
            return self.exp_h_gpa.unwrap_or(0.0); // use experimental value directly
        }
        // For carbon: nb from valence/nc; for others: use exp_h_gpa override
        let nc_eff = self.nc;
        let nb_eff = (4.0 / nc_eff).min(1.0); // approximate: all shown elements have ~4 bonds/atom in bulk
        H_FACTOR * nb_eff * cohen_b(self.d_cc, self.lam, nb_eff, self.f3d)
    }

    fn e_bond(&self) -> f64 {
        if self.d_cc == 0.0 { return 0.0; }
        self.n_sig * e_sigma(self.d_cc) + self.n_pi * e_pi(self.d_cc)
    }

    fn e_coh(&self) -> f64 {
        (self.nc / 2.0) * self.e_bond()
    }

    /// Annihilation energy per atom pair (anti-atom + matter atom)
    fn ann_per_atom_gev(&self) -> f64 {
        ann_energy_per_atom_gev(self.mass_amu)
    }

    fn ann_density_pj_cm3(&self) -> f64 {
        ann_energy_density_pj_cm3(self.rho())
    }

    /// Figure of merit: structural × energetic
    /// FoM = H_proxy^(1/3) × E_ann_density^(2/3)  (hardness for containment, energy for yield)
    fn figure_of_merit(&self) -> f64 {
        self.hardness().powf(1.0 / 3.0) * self.ann_density_pj_cm3().powf(2.0 / 3.0)
    }
}

fn make_phases() -> Vec<Phase> {
    vec![
        // ── Carbon allotropes ─────────────────────────────────────────────────
        Phase {
            name: "anti-Lonsdaleite (hex diamond, P6₃/mmc)",
            element: "C", mass_amu: 12.0,
            nc: 4.0, d_cc: 1.520, n_sig: 1.0, n_pi: 0.0,
            v_atom: 5.655, lam: 0.0, f3d: 1.0,
            category: "anti-carbon",
            exp_h_gpa: Some(79.0),
        },
        Phase {
            name: "anti-Diamond (cubic, Fd3̄m)",
            element: "C", mass_amu: 12.0,
            nc: 4.0, d_cc: 1.545, n_sig: 1.0, n_pi: 0.0,
            v_atom: 5.673, lam: 0.0, f3d: 1.0,
            category: "anti-carbon",
            exp_h_gpa: Some(70.0),
        },
        Phase {
            name: "anti-K4 carbon (3D sp², ita network)",
            element: "C", mass_amu: 12.0,
            nc: 3.0, d_cc: 1.438, n_sig: 1.0, n_pi: 1.0 / 3.0,
            v_atom: 12.0, lam: 0.5, f3d: 0.85,
            category: "anti-carbon",
            exp_h_gpa: None,
        },
        Phase {
            name: "anti-Graphite (ABAB, P6₃/mmc)",
            element: "C", mass_amu: 12.0,
            nc: 3.0, d_cc: 1.421, n_sig: 1.0, n_pi: 1.0 / 3.0,
            v_atom: 8.783, lam: 0.5, f3d: 0.03,
            category: "anti-carbon",
            exp_h_gpa: Some(0.3),
        },
        Phase {
            name: "SC-anti-carbon (GUTOE field lattice, N_c=6)",
            element: "C", mass_amu: 12.0,
            nc: 6.0, d_cc: 1.570, n_sig: 2.0 / 3.0, n_pi: 0.0,
            v_atom: 1.570_f64.powi(3), lam: 1.0, f3d: 1.0,
            category: "anti-carbon",
            exp_h_gpa: None,
        },
        // ── Anti-boron nitride phases ─────────────────────────────────────────
        Phase {
            name: "anti-c-BN (cubic boron nitride, Fd3̄m)",
            element: "BN", mass_amu: 24.818, // (10.811 + 14.007) avg per atom = 12.409; per pair = 24.818
            // Actually per formula unit BN: mass = 10.811 + 14.007 = 24.818 amu
            // V_atom per BN pair = 11.81 Å³ (a = 3.615 Å, 4 BN per cell → V = 11.81 Å³/pair)
            nc: 4.0, d_cc: 1.565, n_sig: 1.0, n_pi: 0.0,
            v_atom: 11.81 / 2.0, // Å³ per atom (not per formula unit)
            lam: 0.5, f3d: 1.0, // slightly ionic (B-N vs C-C)
            category: "anti-BN",
            exp_h_gpa: Some(46.0),
        },
        // ── Heavy-element antimaterials (use experimental data) ───────────────
        Phase {
            name: "anti-Osmium (HCP, most dense element)",
            element: "Os", mass_amu: 190.23,
            nc: 12.0, d_cc: 0.0, // use exp hardness directly
            n_sig: 0.0, n_pi: 0.0,
            // HCP Os: a=2.7341Å, c=4.3197Å, 2 atoms/cell → V = √3/2×a²×c/2 = 13.98 Å³
            v_atom: 13.98,
            lam: 0.0, f3d: 1.0,
            category: "anti-heavy",
            exp_h_gpa: Some(38.4), // Vickers ~3920 HV → ~38.4 GPa
        },
        Phase {
            name: "anti-Iridium (FCC, hardest platinum-group metal)",
            element: "Ir", mass_amu: 192.217,
            nc: 12.0, d_cc: 0.0,
            n_sig: 0.0, n_pi: 0.0,
            // FCC Ir: a=3.8390Å, 4 atoms → V/atom = 14.16 Å³
            v_atom: 14.16,
            lam: 0.0, f3d: 1.0,
            category: "anti-heavy",
            exp_h_gpa: Some(15.7), // Vickers ~1600 HV → 15.7 GPa
        },
        Phase {
            name: "anti-Rhenium (HCP, high hardness + density)",
            element: "Re", mass_amu: 186.207,
            nc: 12.0, d_cc: 0.0,
            n_sig: 0.0, n_pi: 0.0,
            // HCP Re: a=2.7608Å, c=4.4580Å, 2 atoms → V/atom = 14.71 Å³
            v_atom: 14.71,
            lam: 0.0, f3d: 1.0,
            category: "anti-heavy",
            exp_h_gpa: Some(13.0),
        },
        Phase {
            name: "anti-Tungsten (BCC, highest melting point)",
            element: "W", mass_amu: 183.84,
            nc: 8.0, d_cc: 0.0,
            n_sig: 0.0, n_pi: 0.0,
            // BCC W: a=3.1648Å, 2 atoms → V/atom = 15.86 Å³
            v_atom: 15.86,
            lam: 0.0, f3d: 1.0,
            category: "anti-heavy",
            exp_h_gpa: Some(8.0), // ~800 HV
        },
        Phase {
            name: "anti-Uranium (α-U, fissile, dense)",
            element: "U", mass_amu: 238.029,
            nc: 0.0, d_cc: 0.0,
            n_sig: 0.0, n_pi: 0.0,
            // α-U: orthorhombic, ρ = 19050 kg/m³ → V/atom = 238.029×1.66054e-27/19050e3/1e-30
            v_atom: 238.029 * AMU_KG / (19050.0 * A3_TO_M3), // Å³/atom
            lam: 0.0, f3d: 1.0,
            category: "anti-heavy",
            exp_h_gpa: Some(4.0),
        },
        // ── Reference: anti-hydrogen (for perspective) ───────────────────────
        Phase {
            name: "anti-Hydrogen (liquid, for reference)",
            element: "H", mass_amu: 1.008,
            nc: 0.0, d_cc: 0.0,
            n_sig: 0.0, n_pi: 0.0,
            // liquid H₂ density ≈ 70.8 kg/m³; V/atom = 2×1.008×AMU/70.8
            v_atom: 2.0 * 1.008 * AMU_KG / (70.8 * A3_TO_M3),
            lam: 0.0, f3d: 0.0,
            category: "anti-light",
            exp_h_gpa: None,
        },
    ]
}

// ─── Assertions ───────────────────────────────────────────────────────────────

fn run_assertions(phases: &[Phase]) {
    // 1. CPT: anti-diamond = diamond in hardness (by construction, same formula)
    let ad = phases.iter().find(|p| p.name.starts_with("anti-Diamond")).unwrap();
    assert!((ad.hardness() - 98.9).abs() < 1.0, "anti-diamond hardness mismatch");

    // 2. anti-lonsdaleite harder than anti-diamond
    let al = phases.iter().find(|p| p.name.starts_with("anti-Lonsdaleite")).unwrap();
    assert!(al.hardness() > ad.hardness(), "anti-lonsdaleite should be harder than anti-diamond");

    // 3. Annihilation energy per atom scales linearly with mass
    let e_c  = ann_energy_per_atom_gev(12.0);
    let e_os = ann_energy_per_atom_gev(190.23);
    assert!((e_os / e_c - 190.23 / 12.0).abs() < 0.01,
        "annihilation energy should scale linearly with mass");

    // 4. Anti-osmium denser than anti-diamond
    let os = phases.iter().find(|p| p.name.starts_with("anti-Osmium")).unwrap();
    assert!(os.rho() > ad.rho() * 3.0, "osmium should be >> 3× denser than diamond");

    // 5. Annihilation energy per kg is 2c² regardless of element
    let e_per_kg_c  = ann_energy_density_j_m3(density_kg_m3(12.0, 5.673)) / density_kg_m3(12.0, 5.673);
    let e_per_kg_os = ann_energy_density_j_m3(density_kg_m3(190.23, 13.98)) / density_kg_m3(190.23, 13.98);
    assert!((e_per_kg_c - e_per_kg_os).abs() / e_per_kg_c < 1e-10,
        "annihilation energy per kg must be 2c² = const for all matter");

    // 6. SC-anti-carbon denser than anti-diamond
    let sc = phases.iter().find(|p| p.name.starts_with("SC-anti-carbon")).unwrap();
    assert!(sc.rho() > ad.rho(), "SC-anti-carbon should be denser than anti-diamond");

    // 7. 2c² = 179.9 PJ/kg (fundamental constant check)
    let two_c2 = 2.0 * C * C;
    assert!((two_c2 - 1.799e17).abs() / 1.799e17 < 0.001,
        "2c² = {:.4e} J/kg, expected ~1.799e17", two_c2);

    println!("  [PASS] All 7 assertions satisfied");
}

// ─── Output ───────────────────────────────────────────────────────────────────

fn write_txt(phases: &[Phase], out: &str) {
    let two_c2_pj_kg = 2.0 * C * C * 1e-15; // PJ/kg
    let ann_per_kg_pj = two_c2_pj_kg;

    let mut s = String::new();
    let _ = writeln!(s, "╔══════════════════════════════════════════════════════════════════════════╗");
    let _ = writeln!(s, "║     GUTOE: Antimaterial Phase Analysis                                ║");
    let _ = writeln!(s, "║     CPT invariance + annihilation energy density                      ║");
    let _ = writeln!(s, "╚══════════════════════════════════════════════════════════════════════════╝");
    let _ = writeln!(s);
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  CPT SYMMETRY AND WHAT IT MEANS FOR ANTIMATERIALS");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  Anti-carbon obeys the same α = {:.9}", ALPHA);
    let _ = writeln!(s, "  Same mass, same Bohr radius, same bond lengths, same everything.");
    let _ = writeln!(s, "  Therefore: anti-allotrope hardness ranking = matter allotrope ranking.");
    let _ = writeln!(s, "  Anti-lonsdaleite IS as hard as lonsdaleite.  Anti-diamond IS as hard");
    let _ = writeln!(s, "  as diamond.  This is CPT symmetry — not an approximation.");
    let _ = writeln!(s);
    let _ = writeln!(s, "  The NEW quantity: ANNIHILATION ENERGY DENSITY = 2ρc²");
    let _ = writeln!(s, "    Per kilogram of antimatter: 2c² = {:.3} PJ/kg (universal constant)", ann_per_kg_pj);
    let _ = writeln!(s, "    Per unit volume: 2ρc²  (scales with density — this is the variable)");
    let _ = writeln!(s, "    Denser antimaterial ↔ more energy per cm³");
    let _ = writeln!(s);

    // Anti-carbon section
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  ANTI-CARBON ALLOTROPE RANKING");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  {:>4}  {:<44} {:>8} {:>8} {:>9}",
        "Rank", "Phase", "ρ(kg/m³)", "H(GPa)", "E_ann(PJ/cm³)");
    let _ = writeln!(s, "  {:>4}  {:<44} {:>8} {:>8} {:>9}",
        "----", "", "--------", "------", "-------");

    let mut carbons: Vec<&Phase> = phases.iter().filter(|p| p.category == "anti-carbon").collect();
    carbons.sort_by(|a, b| b.ann_density_pj_cm3().partial_cmp(&a.ann_density_pj_cm3()).unwrap());

    for (i, p) in carbons.iter().enumerate() {
        let h_note = if p.name.starts_with("anti-Lonsdaleite") { " ★ hardest" }
                     else if p.name.starts_with("SC-anti") { " ♦ densest" }
                     else { "" };
        let _ = writeln!(s, "  {:>4}  {:<44} {:>8.0} {:>8.1} {:>9.6}{}",
            i + 1,
            p.name.get(..42).unwrap_or(p.name),
            p.rho(), p.hardness(), p.ann_density_pj_cm3(), h_note);
    }
    let _ = writeln!(s);
    let _ = writeln!(s, "  ★ = hardest (by CPT = same as lonsdaleite)");
    let _ = writeln!(s, "  ♦ = densest carbon phase (GUTOE field lattice, N_c=6, requires >3 TPa)");
    let _ = writeln!(s);

    // SC-carbon vs lonsdaleite comparison
    let sc  = phases.iter().find(|p| p.name.starts_with("SC-anti")).unwrap();
    let lon = phases.iter().find(|p| p.name.starts_with("anti-Lonsdaleite")).unwrap();
    let dia = phases.iter().find(|p| p.name.starts_with("anti-Diamond")).unwrap();
    let _ = writeln!(s, "  ┌─────────────────────────┬──────────────┬──────────────┬──────────────┐");
    let _ = writeln!(s, "  │ Property                 │ anti-Diamond │anti-Lonsdale │ SC-anti-C(♦) │");
    let _ = writeln!(s, "  ├─────────────────────────┼──────────────┼──────────────┼──────────────┤");
    let _ = writeln!(s, "  │ Density (kg/m³)          │  {:>9.0}   │  {:>9.0}   │  {:>9.0}   │",
        dia.rho(), lon.rho(), sc.rho());
    let _ = writeln!(s, "  │ Hardness proxy (GPa)     │  {:>9.1}   │  {:>9.1}   │  {:>9.1}   │",
        dia.hardness(), lon.hardness(), sc.hardness());
    let _ = writeln!(s, "  │ E_ann per atom (GeV)     │  {:>9.3}   │  {:>9.3}   │  {:>9.3}   │",
        dia.ann_per_atom_gev(), lon.ann_per_atom_gev(), sc.ann_per_atom_gev());
    let _ = writeln!(s, "  │ E_ann density (PJ/cm³)   │  {:>9.6}   │  {:>9.6}   │  {:>9.6}   │",
        dia.ann_density_pj_cm3(), lon.ann_density_pj_cm3(), sc.ann_density_pj_cm3());
    let sc_ratio = sc.ann_density_pj_cm3() / dia.ann_density_pj_cm3();
    let _ = writeln!(s, "  │ E_ann vs anti-diamond    │       1.000× │  {:>6.3}×    │  {:>6.3}×    │",
        lon.ann_density_pj_cm3() / dia.ann_density_pj_cm3(), sc_ratio);
    let _ = writeln!(s, "  └─────────────────────────┴──────────────┴──────────────┴──────────────┘");
    let _ = writeln!(s);
    let _ = writeln!(s, "  SC-anti-carbon (GUTOE field lattice) is {:.1}% denser than anti-diamond.", (sc_ratio - 1.0)*100.0);
    let _ = writeln!(s, "  It is the DENSEST carbon phase the algebra predicts — but requires >3 TPa.");
    let _ = writeln!(s, "  For accessible synthesis: anti-lonsdaleite (metastable, HPHT) is the winner.");
    let _ = writeln!(s);

    // Heavy elements
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  ANTIMATERIAL RANKING — ALL ELEMENTS");
    let _ = writeln!(s, "  Sorted by annihilation energy density (E_ann = 2ρc²)");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  {:>4}  {:<46} {:>9} {:>8} {:>12}  {}",
        "Rank", "Phase", "ρ(kg/m³)", "H(GPa)", "E_ann(PJ/cm³)", "FoM");
    let _ = writeln!(s, "  {:>4}  {:<46} {:>9} {:>8} {:>12}  {}",
        "----", "", "---------", "------", "------------", "---");

    let mut all: Vec<&Phase> = phases.iter().collect();
    all.sort_by(|a, b| b.ann_density_pj_cm3().partial_cmp(&a.ann_density_pj_cm3()).unwrap());

    for (i, p) in all.iter().enumerate() {
        let h_str = if p.d_cc == 0.0 && p.exp_h_gpa.is_some() {
            format!("{:>8.1}†", p.exp_h_gpa.unwrap())
        } else if p.d_cc == 0.0 {
            format!("{:>8}", "—")
        } else {
            format!("{:>8.1}", p.hardness())
        };
        let fom_str = if p.exp_h_gpa.is_some() || p.d_cc != 0.0 {
            format!("{:.3}", p.figure_of_merit())
        } else {
            "—".to_string()
        };
        let _ = writeln!(s, "  {:>4}  {:<46} {:>9.0} {} {:>12.6}  {}",
            i + 1,
            p.name.get(..44).unwrap_or(p.name),
            p.rho(), h_str, p.ann_density_pj_cm3(), fom_str);
    }
    let _ = writeln!(s, "  † = experimental Vickers hardness");
    let _ = writeln!(s);

    // Scale examples
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  ANNIHILATION ENERGY FOR 1 cm³ SAMPLES");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  1 cm³ of antimatter (+ equal mass of matter) releases:");
    let _ = writeln!(s);
    let tnt_per_j = 1.0 / 4.184e9; // kt TNT per J
    for p in &all {
        let e_j = p.ann_density_pj_cm3() * 1e15; // PJ/cm³ × 1e15 J/PJ → J per cm³ sample
        let kt = e_j * tnt_per_j;
        let pj = e_j * 1e-15;
        let _ = writeln!(s, "  {:<46} {:>8.4} PJ  = {:>8.1} kt TNT",
            p.name.get(..44).unwrap_or(p.name), pj, kt);
    }
    let _ = writeln!(s);

    // GUTOE two-lattice interpretation
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  GUTOE TWO-LATTICE INTERPRETATION FOR ANTIMATERIALS");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  GUTOE derives two distinct lattice optima from Cl(1,3):");
    let _ = writeln!(s);
    let _ = writeln!(s, "  Matter lattice  N_c=4 (tetrahedral sp³):");
    let _ = writeln!(s, "    Optimises:  HARDNESS (structural integrity)");
    let _ = writeln!(s, "    Winner:     anti-lonsdaleite, H = {:.1} GPa", lon.hardness());
    let _ = writeln!(s, "    Use case:   containment wall, structural antimaterial under magnetic");
    let _ = writeln!(s, "                confinement — survives longest before surface annihilation");
    let _ = writeln!(s);
    let _ = writeln!(s, "  Field lattice   N_c=6 (SC from {{γ¹², γ¹³, γ²³}}):");
    let _ = writeln!(s, "    Optimises:  ENERGY DENSITY (per unit volume)");
    let _ = writeln!(s, "    Winner:     SC-anti-carbon, E_ann = {:.6} PJ/cm³", sc.ann_density_pj_cm3());
    let _ = writeln!(s, "    Use case:   maximum energy payload; requires extreme synthesis");
    let _ = writeln!(s);
    let _ = writeln!(s, "  Across ALL elements:");
    let os = phases.iter().find(|p| p.name.starts_with("anti-Osmium")).unwrap();
    let _ = writeln!(s, "    Energy density winner: anti-osmium, ρ = {:.0} kg/m³, E_ann = {:.6} PJ/cm³",
        os.rho(), os.ann_density_pj_cm3());
    let _ = writeln!(s, "    Hardness winner: anti-lonsdaleite (carbon), H = {:.1} GPa", lon.hardness());
    let _ = writeln!(s, "    Combined FoM winner: anti-osmium (heavy element density wins energy density)");
    let _ = writeln!(s);
    let _ = writeln!(s, "  Key insight: within carbon phases, the algebra identifies SC-anti-carbon");
    let _ = writeln!(s, "  as the energy-density champion and anti-lonsdaleite as the structural");
    let _ = writeln!(s, "  champion. These are the same two lattices the Clifford algebra uses for");
    let _ = writeln!(s, "  the electromagnetic field and matter respectively. The duality holds.");
    let _ = writeln!(s);

    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");
    let _ = writeln!(s, "  End of report — GUTOE antimaterial phase analysis");
    let _ = writeln!(s, "═══════════════════════════════════════════════════════════════════════════");

    fs::write(format!("{out}/antimaterial_hardness_sim.txt"), &s).expect("write txt");
    println!("  → {out}/antimaterial_hardness_sim.txt");
}

fn write_csv(phases: &[Phase], out: &str) {
    let mut s = String::new();
    let _ = writeln!(s, "name,element,category,mass_amu,v_atom_A3,rho_kg_m3,nc,d_cc_A,n_sig,n_pi,e_bond_ev,e_coh_ev,hardness_gpa,ann_per_atom_gev,ann_density_pj_cm3,figure_of_merit");
    for p in phases {
        let h = if p.d_cc == 0.0 && p.exp_h_gpa.is_some() { p.exp_h_gpa.unwrap() } else { p.hardness() };
        let _ = writeln!(s, "{},{},{},{:.4},{:.4},{:.2},{:.1},{:.4},{:.3},{:.3},{:.4},{:.4},{:.4},{:.6},{:.8},{:.4}",
            p.name.replace(',', ";"),
            p.element, p.category, p.mass_amu, p.v_atom,
            p.rho(), p.nc, p.d_cc, p.n_sig, p.n_pi,
            p.e_bond(), p.e_coh(), h,
            p.ann_per_atom_gev(),
            p.ann_density_pj_cm3(),
            p.figure_of_merit());
    }
    fs::write(format!("{out}/antimaterial_hardness_sim.csv"), &s).expect("write csv");
    println!("  → {out}/antimaterial_hardness_sim.csv");
}

fn write_json(phases: &[Phase], out: &str) {
    let lon = phases.iter().find(|p| p.name.starts_with("anti-Lonsdaleite")).unwrap();
    let sc  = phases.iter().find(|p| p.name.starts_with("SC-anti")).unwrap();
    let os  = phases.iter().find(|p| p.name.starts_with("anti-Osmium")).unwrap();

    let mut s = String::new();
    let _ = writeln!(s, "{{");
    let _ = writeln!(s, "  \"cpt_symmetry\": \"anti-allotrope hardness = matter allotrope hardness\",");
    let _ = writeln!(s, "  \"ann_energy_per_kg_pj\": {:.6},", 2.0 * C * C * 1e-15);
    let _ = writeln!(s, "  \"carbon_hardness_winner\": {{");
    let _ = writeln!(s, "    \"name\": \"{}\",", lon.name);
    let _ = writeln!(s, "    \"hardness_gpa\": {:.2},", lon.hardness());
    let _ = writeln!(s, "    \"rho_kg_m3\": {:.1},", lon.rho());
    let _ = writeln!(s, "    \"ann_density_pj_cm3\": {:.8}", lon.ann_density_pj_cm3());
    let _ = writeln!(s, "  }},");
    let _ = writeln!(s, "  \"carbon_energy_density_winner\": {{");
    let _ = writeln!(s, "    \"name\": \"{}\",", sc.name);
    let _ = writeln!(s, "    \"hardness_gpa\": {:.2},", sc.hardness());
    let _ = writeln!(s, "    \"rho_kg_m3\": {:.1},", sc.rho());
    let _ = writeln!(s, "    \"ann_density_pj_cm3\": {:.8}", sc.ann_density_pj_cm3());
    let _ = writeln!(s, "  }},");
    let _ = writeln!(s, "  \"all_elements_energy_density_winner\": {{");
    let _ = writeln!(s, "    \"name\": \"{}\",", os.name);
    let _ = writeln!(s, "    \"rho_kg_m3\": {:.1},", os.rho());
    let _ = writeln!(s, "    \"ann_density_pj_cm3\": {:.8}", os.ann_density_pj_cm3());
    let _ = writeln!(s, "  }},");
    let _ = writeln!(s, "  \"gutoe_duality\": {{");
    let _ = writeln!(s, "    \"matter_lattice\": {{ \"nc\": 4, \"optimises\": \"hardness\", \"winner\": \"anti-lonsdaleite\" }},");
    let _ = writeln!(s, "    \"field_lattice\":  {{ \"nc\": 6, \"optimises\": \"energy_density\", \"winner\": \"SC-anti-carbon\" }}");
    let _ = writeln!(s, "  }},");
    let _ = writeln!(s, "  \"phases\": [");
    let n = phases.len();
    for (i, p) in phases.iter().enumerate() {
        let comma = if i + 1 < n { "," } else { "" };
        let h = if p.d_cc == 0.0 && p.exp_h_gpa.is_some() { p.exp_h_gpa.unwrap() } else { p.hardness() };
        let _ = writeln!(s, "    {{ \"name\": \"{}\", \"rho\": {:.1}, \"hardness_gpa\": {:.2}, \"ann_pj_cm3\": {:.8} }}{}",
            p.name, p.rho(), h, p.ann_density_pj_cm3(), comma);
    }
    let _ = writeln!(s, "  ]");
    let _ = write!(s, "}}");

    fs::write(format!("{out}/antimaterial_hardness_sim.json"), &s).expect("write json");
    println!("  → {out}/antimaterial_hardness_sim.json");
}

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let out = std::env::var("GUTOE_PHYSICS_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/carbon_hardness_lattice".to_string());
    fs::create_dir_all(&out).expect("create output dir");

    println!("GUTOE: Antimaterial Phase Analysis");
    println!("  2c² = {:.4} PJ/kg  (universal antimatter energy yield)", 2.0 * C * C * 1e-15);
    println!("  Output → {out}");
    println!();

    let phases = make_phases();

    println!("Running assertions...");
    run_assertions(&phases);
    println!();

    println!("Writing outputs...");
    write_txt(&phases, &out);
    write_csv(&phases, &out);
    write_json(&phases, &out);
    println!();

    // ── Summary to stdout ────────────────────────────────────────────────────
    println!("Anti-carbon phases (sorted by E_ann density):");
    let mut carbons: Vec<&Phase> = phases.iter().filter(|p| p.category == "anti-carbon").collect();
    carbons.sort_by(|a, b| b.ann_density_pj_cm3().partial_cmp(&a.ann_density_pj_cm3()).unwrap());
    for p in &carbons {
        println!("  ρ={:>5.0} kg/m³  H={:>5.0} GPa  E_ann={:.6} PJ/cm³  {}",
            p.rho(), p.hardness(), p.ann_density_pj_cm3(), p.name);
    }
    println!();

    println!("All phases (sorted by E_ann density):");
    let mut all: Vec<&Phase> = phases.iter().collect();
    all.sort_by(|a, b| b.ann_density_pj_cm3().partial_cmp(&a.ann_density_pj_cm3()).unwrap());
    for p in all {
        let h_str = if p.d_cc == 0.0 && p.exp_h_gpa.is_some() {
            format!("{:>5.0}†", p.exp_h_gpa.unwrap())
        } else if p.d_cc == 0.0 {
            format!("{:>5}", "—")
        } else {
            format!("{:>5.0}", p.hardness())
        };
        println!("  ρ={:>6.0} kg/m³  H={} GPa  E_ann={:.6} PJ/cm³  {}",
            p.rho(), h_str, p.ann_density_pj_cm3(), p.name);
    }
    println!();

    let lon = phases.iter().find(|p| p.name.starts_with("anti-Lonsdaleite")).unwrap();
    let sc  = phases.iter().find(|p| p.name.starts_with("SC-anti")).unwrap();
    let os  = phases.iter().find(|p| p.name.starts_with("anti-Osmium")).unwrap();
    println!("GUTOE duality:");
    println!("  Hardest anti-carbon:      anti-lonsdaleite    H = {:.1} GPa", lon.hardness());
    println!("  Densest anti-carbon:      SC-anti-carbon (N_c=6 field lattice)  E = {:.6} PJ/cm³", sc.ann_density_pj_cm3());
    println!("  Highest energy yield:     anti-osmium         E = {:.6} PJ/cm³", os.ann_density_pj_cm3());
    println!("  The field lattice (N_c=6) optimises energy density.");
    println!("  The matter lattice (N_c=4) optimises hardness.");
    println!("  Same two-lattice split, different optimisation target.");
}
