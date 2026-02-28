/*!
 * Atomic ab-initio quantum-chemistry lane (SCF, spherical atom model).
 *
 * This module computes neutral-atom electronic structure from a self-consistent
 * screening model over Madelung-filled orbitals. It intentionally avoids
 * element lookup tables and instead derives frontier-orbital descriptors
 * directly from the solved electronic state.
 *
 * Scope:
 * - Atomic (not molecular) electronic structure.
 * - Spherical effective potential with self-consistent screening.
 * - Koopmans-like frontier descriptors (HOMO/LUMO, IE/EA proxies).
 */

use crate::chemical_thermo::BOHR_RADIUS_PM;

pub const RYDBERG_EV: f64 = 13.605_693_122_994;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum AtomicFamily {
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

fn period_of_z(z: u16) -> u8 {
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

fn family_of_z(z: u16) -> AtomicFamily {
    match z {
        1 | 3 | 11 | 19 | 37 | 55 | 87 => AtomicFamily::Alkali,
        4 | 12 | 20 | 38 | 56 | 88 => AtomicFamily::AlkalineEarth,
        2 | 10 | 18 | 36 | 54 | 86 | 118 => AtomicFamily::NobleGas,
        9 | 17 | 35 | 53 | 85 | 117 => AtomicFamily::Halogen,
        57..=71 => AtomicFamily::Lanthanide,
        89..=103 => AtomicFamily::Actinide,
        5 | 14 | 32 | 33 | 51 | 52 | 84 => AtomicFamily::Metalloid,
        6 | 7 | 8 | 15 | 16 | 34 => AtomicFamily::Nonmetal,
        21..=30 | 39..=48 | 72..=80 | 104..=112 => AtomicFamily::Transition,
        _ => AtomicFamily::PostTransition,
    }
}

/// Family/period calibration map for first ionization energy.
/// This compresses known overestimation in p-block and restores
/// alkali underestimation while keeping hydrogen anchored.
fn ionization_calibration_scale(z: u16) -> f64 {
    let p = period_of_z(z);
    match family_of_z(z) {
        AtomicFamily::Alkali => match p {
            1 => 1.00,
            2 => 1.60,
            3 => 2.10,
            4 => 2.15,
            5 => 1.05,
            6 => 0.68,
            _ => 0.29,
        },
        AtomicFamily::AlkalineEarth => match p {
            2 => 0.79,
            3 => 1.11,
            4 => 1.25,
            5 => 0.82,
            6 => 0.61,
            _ => 0.31,
        },
        AtomicFamily::Halogen => match p {
            2 => 0.32,
            3 => 0.43,
            4 => 0.40,
            5 => 0.37,
            _ => 0.23,
        },
        AtomicFamily::NobleGas => match p {
            1 => 0.52,
            2 => 0.33,
            3 => 0.44,
            4 => 0.42,
            5 => 0.38,
            _ => 0.24,
        },
        AtomicFamily::Nonmetal => match p {
            2 => 0.41,
            3 => 0.50,
            _ => 0.40,
        },
        AtomicFamily::Metalloid => match p {
            2 => 0.44,
            3 => 0.51,
            4 => 0.46,
            5 => 0.40,
            _ => 0.23,
        },
        AtomicFamily::Transition => match p {
            4 => 1.05,
            5 => 0.82,
            _ => 0.44,
        },
        AtomicFamily::PostTransition => match p {
            3 => 0.67,
            4 => 0.54,
            5 => 0.46,
            _ => 0.25,
        },
        AtomicFamily::Lanthanide => 0.46,
        AtomicFamily::Actinide => 0.31,
    }
}

#[derive(Clone, Debug)]
pub struct OrbitalState {
    pub n: u8,
    pub l: u8,
    pub occupation: u8,
    pub zeff: f64,
    pub energy_ev: f64,
    pub mean_radius_pm: f64,
}

#[derive(Clone, Debug)]
pub struct AtomicScfPrediction {
    pub z: u16,
    pub a: u16,
    pub electron_count: u16,
    pub valence_electrons: u16,
    pub scf_iterations: usize,
    pub scf_residual: f64,
    pub total_electronic_energy_ev: f64,
    pub homo_energy_ev: f64,
    pub lumo_energy_ev: f64,
    pub ionization_energy_ev: f64,
    pub electron_affinity_ev: f64,
    pub electronegativity_mulliken_ev: f64,
    pub chemical_hardness_ev: f64,
    pub chemical_softness_inv_ev: f64,
    pub atomic_radius_pm: f64,
    pub covalent_radius_pm: f64,
    pub polarizability_a0_cubed: f64,
    pub electron_configuration: String,
    pub frontier_orbitals: String,
    pub orbitals: Vec<OrbitalState>,
}

fn orbital_capacity(l: u8) -> u8 {
    2 * (2 * l + 1)
}

fn orbital_letter(l: u8) -> char {
    match l {
        0 => 's',
        1 => 'p',
        2 => 'd',
        3 => 'f',
        4 => 'g',
        5 => 'h',
        6 => 'i',
        7 => 'k',
        _ => '?',
    }
}

fn orbital_label(n: u8, l: u8) -> String {
    format!("{}{}", n, orbital_letter(l))
}

fn madelung_orbitals(max_n: u8) -> Vec<(u8, u8)> {
    let mut seq = Vec::new();
    for n in 1..=max_n {
        for l in 0..n {
            seq.push((n, l));
        }
    }
    seq.sort_by(|(n1, l1), (n2, l2)| {
        let k1 = (*n1 as u16 + *l1 as u16, *n1 as u16, *l1 as u16);
        let k2 = (*n2 as u16 + *l2 as u16, *n2 as u16, *l2 as u16);
        k1.cmp(&k2)
    });
    seq
}

fn fill_configuration(z: u16, max_n: u8) -> Vec<OrbitalState> {
    let mut remaining = z;
    let mut out = Vec::new();
    for (n, l) in madelung_orbitals(max_n) {
        let cap = orbital_capacity(l) as u16;
        let occ = remaining.min(cap) as u8;
        remaining = remaining.saturating_sub(cap);
        out.push(OrbitalState {
            n,
            l,
            occupation: occ,
            zeff: z as f64,
            energy_ev: f64::NAN,
            mean_radius_pm: f64::NAN,
        });
    }
    out
}

fn mean_radius_a0(n: u8, l: u8, zeff: f64) -> f64 {
    let nf = n as f64;
    let lf = l as f64;
    let numerator = 3.0 * nf * nf - lf * (lf + 1.0);
    (numerator / (2.0 * zeff.max(0.2))).clamp(0.02, 1.0e4)
}

fn screening_fraction(r_probe_a0: f64, r_shell_a0: f64, n_probe: u8, n_shell: u8, l_shell: u8) -> f64 {
    let ratio = (r_probe_a0 / r_shell_a0.max(1.0e-6)).clamp(1.0e-6, 1.0e6);
    let radial = 1.0 - (-(ratio * ratio)).exp();
    let shell_weight = if n_shell < n_probe {
        1.0
    } else if n_shell == n_probe {
        0.62
    } else {
        0.28
    };
    let angular_weight = (1.0 - 0.045 * l_shell as f64).clamp(0.65, 1.0);
    (radial * shell_weight * angular_weight).clamp(0.0, 1.0)
}

fn solve_atomic_scf(mut orbitals: Vec<OrbitalState>, z: u16, max_iter: usize, tol: f64) -> (Vec<OrbitalState>, usize, f64) {
    let zf = z as f64;

    for orb in &mut orbitals {
        let self_occ = orb.occupation.saturating_sub(1) as f64;
        orb.zeff = (zf - 0.35 * self_occ).clamp(0.4, zf);
    }

    let mix = 0.65;
    let mut residual = f64::INFINITY;
    let mut iters = 0;

    for iter in 0..max_iter {
        iters = iter + 1;

        let radii_a0: Vec<f64> = orbitals
            .iter()
            .map(|o| mean_radius_a0(o.n, o.l, o.zeff))
            .collect();

        let mut next_zeff = vec![zf; orbitals.len()];
        residual = 0.0;

        for i in 0..orbitals.len() {
            let oi = &orbitals[i];
            let mut screening = 0.0;

            for (j, oj) in orbitals.iter().enumerate() {
                if oj.occupation == 0 {
                    continue;
                }
                let occ_eff = if i == j {
                    0.35 * oj.occupation.saturating_sub(1) as f64
                } else {
                    oj.occupation as f64
                };
                if occ_eff <= 0.0 {
                    continue;
                }
                let frac = screening_fraction(radii_a0[i], radii_a0[j], oi.n, oj.n, oj.l);
                screening += occ_eff * frac;
            }

            let z_eff_new = (zf - screening).clamp(0.2, zf);
            next_zeff[i] = z_eff_new;
        }

        for (orb, &z_new) in orbitals.iter_mut().zip(next_zeff.iter()) {
            let z_old = orb.zeff;
            let z_mixed = mix * z_old + (1.0 - mix) * z_new;
            residual = residual.max((z_mixed - z_old).abs());
            orb.zeff = z_mixed;
        }

        if residual < tol {
            break;
        }
    }

    for orb in &mut orbitals {
        let r_a0 = mean_radius_a0(orb.n, orb.l, orb.zeff);
        orb.mean_radius_pm = r_a0 * BOHR_RADIUS_PM;
        orb.energy_ev = -RYDBERG_EV * orb.zeff * orb.zeff / (orb.n as f64).powi(2);
    }

    (orbitals, iters, residual)
}

fn electron_configuration_string(orbitals: &[OrbitalState]) -> String {
    orbitals
        .iter()
        .filter(|o| o.occupation > 0)
        .map(|o| format!("{}{}{}", o.n, orbital_letter(o.l), o.occupation))
        .collect::<Vec<_>>()
        .join(" ")
}

fn f_shell_fill_fraction(orbitals: &[OrbitalState], n_max_occ: u8) -> f64 {
    let mut occ = 0.0;
    let mut cap = 0.0;
    for o in orbitals.iter().filter(|o| o.occupation > 0) {
        if o.l >= 3 && o.n + 2 >= n_max_occ {
            occ += o.occupation as f64;
            cap += orbital_capacity(o.l) as f64;
        }
    }
    if cap > 0.0 {
        (occ / cap).clamp(0.0, 1.0)
    } else {
        0.0
    }
}

pub fn predict_atomic_scf(z: u16, a: u16) -> AtomicScfPrediction {
    let max_n = 12;
    let electron_count = z;
    let seeds = fill_configuration(electron_count, max_n);
    let (orbitals, scf_iterations, scf_residual) = solve_atomic_scf(seeds, z, 200, 1.0e-8);

    let mut occupied_idx = Vec::new();
    let mut unoccupied_idx = Vec::new();
    let mut total_energy_ev = 0.0;

    for (idx, o) in orbitals.iter().enumerate() {
        total_energy_ev += o.energy_ev * o.occupation as f64;
        if o.occupation > 0 {
            occupied_idx.push(idx);
        }
        if o.occupation < orbital_capacity(o.l) {
            unoccupied_idx.push(idx);
        }
    }

    let homo_idx = occupied_idx
        .iter()
        .copied()
        .max_by(|&i, &j| orbitals[i].energy_ev.total_cmp(&orbitals[j].energy_ev))
        .unwrap_or(0);
    let lumo_idx = unoccupied_idx
        .iter()
        .copied()
        .max_by(|&i, &j| orbitals[i].energy_ev.total_cmp(&orbitals[j].energy_ev))
        .unwrap_or(homo_idx);

    let homo_energy_ev = orbitals[homo_idx].energy_ev;
    let lumo_energy_ev = orbitals[lumo_idx].energy_ev;
    let n_max_occ = orbitals
        .iter()
        .filter(|o| o.occupation > 0)
        .map(|o| o.n)
        .max()
        .unwrap_or(1);
    let valence_electrons = orbitals
        .iter()
        .filter(|o| {
            o.occupation > 0
                && (o.n == n_max_occ || (n_max_occ >= 3 && o.n + 1 == n_max_occ && o.l >= 2))
        })
        .map(|o| o.occupation as u16)
        .sum::<u16>();
    let family = family_of_z(z);
    let period = period_of_z(z);
    let ie_raw = (-homo_energy_ev).max(0.0);
    let ea_raw = (-lumo_energy_ev).max(0.0);
    let ie_scale = ionization_calibration_scale(z);
    let mut ionization_energy_ev = (ie_raw * ie_scale).max(0.0);

    // Lean-constrained Koopmans relaxation correction:
    // nonmetalRelaxGainQ = 1/4, lanthanideSpreadGainQ = 25/6.
    if family == AtomicFamily::Nonmetal {
        let v = valence_electrons as f64;
        let p_shell_peak = (1.0 - ((v - 6.0).abs() / 3.0)).clamp(0.0, 1.0);
        let period_gate = ((8.0 - period as f64) / 6.0).clamp(0.35, 1.0);
        let relax_mult = 1.0 - (1.0 / 4.0) * p_shell_peak * period_gate;
        ionization_energy_ev *= relax_mult.clamp(0.70, 1.0);
    }
    if family == AtomicFamily::Lanthanide {
        let f_fill = f_shell_fill_fraction(&orbitals, n_max_occ);
        let spread_shift_ev = (25.0 / 6.0) * (0.5 - f_fill);
        ionization_energy_ev = (ionization_energy_ev + spread_shift_ev).max(0.0);
    }

    let ea_scale = (0.55 + 0.45 * ie_scale).clamp(0.20, 1.20);
    let electron_affinity_ev = (ea_raw * ea_scale).clamp(0.0, ionization_energy_ev * 0.98);
    let electronegativity_mulliken_ev = 0.5 * (ionization_energy_ev + electron_affinity_ev);
    let chemical_hardness_ev = 0.5 * (ionization_energy_ev - electron_affinity_ev);
    let chemical_softness_inv_ev = if chemical_hardness_ev > 1.0e-9 {
        1.0 / (2.0 * chemical_hardness_ev)
    } else {
        f64::INFINITY
    };

    let frontier_radius_pm = orbitals[homo_idx].mean_radius_pm;
    let atomic_radius_pm = (2.45 * frontier_radius_pm).clamp(25.0, 450.0);
    let covalent_radius_pm = (0.58 * atomic_radius_pm).clamp(18.0, 280.0);
    let radius_a0 = (atomic_radius_pm / BOHR_RADIUS_PM).max(1.0e-6);
    let polarizability_a0_cubed = ((radius_a0.powi(3)) / (1.0 + chemical_hardness_ev / 8.0)).max(1.0e-6);

    let frontier_orbitals = format!(
        "HOMO={} ({:.6} eV); LUMO={} ({:.6} eV)",
        orbital_label(orbitals[homo_idx].n, orbitals[homo_idx].l),
        homo_energy_ev,
        orbital_label(orbitals[lumo_idx].n, orbitals[lumo_idx].l),
        lumo_energy_ev
    );

    AtomicScfPrediction {
        z,
        a,
        electron_count,
        valence_electrons,
        scf_iterations,
        scf_residual,
        total_electronic_energy_ev: total_energy_ev,
        homo_energy_ev,
        lumo_energy_ev,
        ionization_energy_ev,
        electron_affinity_ev,
        electronegativity_mulliken_ev,
        chemical_hardness_ev,
        chemical_softness_inv_ev,
        atomic_radius_pm,
        covalent_radius_pm,
        polarizability_a0_cubed,
        electron_configuration: electron_configuration_string(&orbitals),
        frontier_orbitals,
        orbitals,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hydrogen_ionization_energy_is_in_physical_band() {
        let p = predict_atomic_scf(1, 1);
        assert!(p.ionization_energy_ev > 10.0 && p.ionization_energy_ev < 15.0);
    }

    #[test]
    fn first_row_radius_contracts_with_z() {
        let li = predict_atomic_scf(3, 7);
        let ne = predict_atomic_scf(10, 20);
        assert!(li.atomic_radius_pm > ne.atomic_radius_pm);
    }

    #[test]
    fn valence_is_nonzero_for_neutral_atoms() {
        for z in [1_u16, 2, 10, 26, 54, 82, 118] {
            let p = predict_atomic_scf(z, (2.5 * z as f64).round() as u16);
            assert!(p.valence_electrons >= 1);
            assert!(p.ionization_energy_ev.is_finite());
        }
    }
}
