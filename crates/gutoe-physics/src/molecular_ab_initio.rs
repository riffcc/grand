/*!
 * Molecular ab-initio lane: restricted Hartree-Fock (RHF) + MP2 correlation.
 *
 * This lane uses a compact Gaussian atomic-orbital basis (s-type primitives)
 * with explicit one- and two-electron integrals, generalized-eigen SCF, and
 * a post-HF MP2 correlation correction on converged molecular orbitals.
 */

use std::f64::consts::PI;

pub const HARTREE_TO_EV: f64 = 27.211_386_245_988;
pub const BOHR_PER_ANGSTROM: f64 = 1.889_726_124_625_770_2;
pub const DEBYE_PER_EA0: f64 = 2.541_746;

#[derive(Clone, Copy, Debug)]
pub struct Atom3D {
    pub z: u16,
    pub x_ang: f64,
    pub y_ang: f64,
    pub z_ang: f64,
}

#[derive(Clone, Debug)]
pub struct MoleculeInput {
    pub name: String,
    pub atoms: Vec<Atom3D>,
    pub charge: i32,
    /// Spin multiplicity 2S+1. Use 1 for closed-shell singlets.
    pub multiplicity: u8,
}

#[derive(Clone, Debug)]
struct BasisFn {
    atom_index: usize,
    alpha: f64,
    center_bohr: [f64; 3],
}

#[derive(Clone, Debug)]
pub struct MolecularAbInitioResult {
    pub name: String,
    pub method: String,
    pub spin_multiplicity: u8,
    pub basis_functions: usize,
    pub electron_count: usize,
    pub alpha_electrons: usize,
    pub beta_electrons: usize,
    pub electron_pairs: usize,
    pub scf_iterations: usize,
    pub scf_residual: f64,
    pub s2_expectation: f64,
    pub nuclear_repulsion_hartree: f64,
    pub electronic_energy_hartree: f64,
    pub total_energy_hartree: f64,
    pub mp2_correlation_hartree: f64,
    pub total_energy_mp2_hartree: f64,
    pub homo_energy_ev: f64,
    pub lumo_energy_ev: f64,
    pub homo_lumo_gap_ev: f64,
    pub dipole_debye: f64,
    pub mulliken_charges: Vec<f64>,
    pub orbital_energies_ev: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct GeometryOptimizationResult {
    pub optimized_molecule: MoleculeInput,
    pub final_result: MolecularAbInitioResult,
    pub iterations: usize,
    pub converged: bool,
    pub final_gradient_norm_hartree_per_angstrom: f64,
}

fn dist2(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn erf_approx(x: f64) -> f64 {
    // Abramowitz-Stegun 7.1.26
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let t = 1.0 / (1.0 + 0.327_591_1 * x.abs());
    let a1 = 0.254_829_592;
    let a2 = -0.284_496_736;
    let a3 = 1.421_413_741;
    let a4 = -1.453_152_027;
    let a5 = 1.061_405_429;
    let poly = (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t;
    sign * (1.0 - poly * (-(x * x)).exp())
}

fn boys0(t: f64) -> f64 {
    if t < 1.0e-8 {
        1.0 - t / 3.0
    } else {
        0.5 * (PI / t).sqrt() * erf_approx(t.sqrt())
    }
}

fn gaussian_norm(alpha: f64) -> f64 {
    (2.0 * alpha / PI).powf(0.75)
}

fn overlap_ss(ai: f64, aj: f64, ra: [f64; 3], rb: [f64; 3]) -> f64 {
    let p = ai + aj;
    let mu = ai * aj / p;
    let rab2 = dist2(ra, rb);
    let ni = gaussian_norm(ai);
    let nj = gaussian_norm(aj);
    ni * nj * (PI / p).powf(1.5) * (-mu * rab2).exp()
}

fn kinetic_ss(ai: f64, aj: f64, ra: [f64; 3], rb: [f64; 3]) -> f64 {
    let p = ai + aj;
    let mu = ai * aj / p;
    let rab2 = dist2(ra, rb);
    let s = overlap_ss(ai, aj, ra, rb);
    mu * (3.0 - 2.0 * mu * rab2) * s
}

fn nuclear_attraction_ss(ai: f64, aj: f64, ra: [f64; 3], rb: [f64; 3], rc: [f64; 3], zc: f64) -> f64 {
    let p = ai + aj;
    let mu = ai * aj / p;
    let rab2 = dist2(ra, rb);
    let pcenter = scale3(add3(scale3(ra, ai), scale3(rb, aj)), 1.0 / p);
    let rpc2 = dist2(pcenter, rc);
    let ni = gaussian_norm(ai);
    let nj = gaussian_norm(aj);
    let pref = -zc * ni * nj * (2.0 * PI / p) * (-mu * rab2).exp();
    pref * boys0(p * rpc2)
}

fn eri_ssss(ai: f64, aj: f64, ak: f64, al: f64, ra: [f64; 3], rb: [f64; 3], rc: [f64; 3], rd: [f64; 3]) -> f64 {
    let p = ai + aj;
    let q = ak + al;
    let mu = ai * aj / p;
    let nu = ak * al / q;

    let rab2 = dist2(ra, rb);
    let rcd2 = dist2(rc, rd);

    let pcenter = scale3(add3(scale3(ra, ai), scale3(rb, aj)), 1.0 / p);
    let qcenter = scale3(add3(scale3(rc, ak), scale3(rd, al)), 1.0 / q);
    let rpq2 = dist2(pcenter, qcenter);

    let ni = gaussian_norm(ai);
    let nj = gaussian_norm(aj);
    let nk = gaussian_norm(ak);
    let nl = gaussian_norm(al);

    let pref = ni
        * nj
        * nk
        * nl
        * (2.0 * PI.powf(2.5))
        / (p * q * (p + q).sqrt())
        * (-(mu * rab2 + nu * rcd2)).exp();

    let t = p * q / (p + q) * rpq2;
    pref * boys0(t)
}

fn idx(i: usize, j: usize, n: usize) -> usize {
    i * n + j
}

fn idx4(i: usize, j: usize, k: usize, l: usize, n: usize) -> usize {
    (((i * n + j) * n + k) * n) + l
}

fn mat_mul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[idx(i, k, n)];
            if aik.abs() < 1.0e-15 {
                continue;
            }
            for j in 0..n {
                c[idx(i, j, n)] += aik * b[idx(k, j, n)];
            }
        }
    }
    c
}

fn mat_transpose(a: &[f64], n: usize) -> Vec<f64> {
    let mut t = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            t[idx(j, i, n)] = a[idx(i, j, n)];
        }
    }
    t
}

fn jacobi_eigen_symmetric(a_in: &[f64], n: usize, max_iter: usize, tol: f64) -> (Vec<f64>, Vec<f64>) {
    let mut a = a_in.to_vec();
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[idx(i, i, n)] = 1.0;
    }

    for _ in 0..max_iter {
        let mut p = 0;
        let mut q = 1.min(n.saturating_sub(1));
        let mut max_off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let val = a[idx(i, j, n)].abs();
                if val > max_off {
                    max_off = val;
                    p = i;
                    q = j;
                }
            }
        }
        if max_off < tol || n <= 1 {
            break;
        }

        let app = a[idx(p, p, n)];
        let aqq = a[idx(q, q, n)];
        let apq = a[idx(p, q, n)];

        let tau = (aqq - app) / (2.0 * apq);
        let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
        } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;

        for k in 0..n {
            if k != p && k != q {
                let akp = a[idx(k, p, n)];
                let akq = a[idx(k, q, n)];
                a[idx(k, p, n)] = c * akp - s * akq;
                a[idx(p, k, n)] = a[idx(k, p, n)];
                a[idx(k, q, n)] = s * akp + c * akq;
                a[idx(q, k, n)] = a[idx(k, q, n)];
            }
        }

        let app_new = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        let aqq_new = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[idx(p, p, n)] = app_new;
        a[idx(q, q, n)] = aqq_new;
        a[idx(p, q, n)] = 0.0;
        a[idx(q, p, n)] = 0.0;

        for k in 0..n {
            let vkp = v[idx(k, p, n)];
            let vkq = v[idx(k, q, n)];
            v[idx(k, p, n)] = c * vkp - s * vkq;
            v[idx(k, q, n)] = s * vkp + c * vkq;
        }
    }

    let mut evals = vec![0.0; n];
    for i in 0..n {
        evals[i] = a[idx(i, i, n)];
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| evals[i].total_cmp(&evals[j]));

    let evals_sorted: Vec<f64> = order.iter().map(|&i| evals[i]).collect();
    let mut evecs_sorted = vec![0.0; n * n];
    for (col_new, &col_old) in order.iter().enumerate() {
        for row in 0..n {
            evecs_sorted[idx(row, col_new, n)] = v[idx(row, col_old, n)];
        }
    }

    (evals_sorted, evecs_sorted)
}

fn symmetric_orthogonalizer(s: &[f64], n: usize) -> Vec<f64> {
    let (evals, u) = jacobi_eigen_symmetric(s, n, 200, 1.0e-12);
    let mut d_inv_sqrt = vec![0.0; n * n];
    for i in 0..n {
        let ev = evals[i].max(1.0e-8);
        d_inv_sqrt[idx(i, i, n)] = 1.0 / ev.sqrt();
    }
    let ut = mat_transpose(&u, n);
    let ud = mat_mul(&u, &d_inv_sqrt, n);
    mat_mul(&ud, &ut, n)
}

fn generalized_eigensolve(f: &[f64], s: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let x = symmetric_orthogonalizer(s, n);
    let xt = mat_transpose(&x, n);
    let fprime = mat_mul(&mat_mul(&xt, f, n), &x, n);
    let (eps, cprime) = jacobi_eigen_symmetric(&fprime, n, 400, 1.0e-11);
    let c = mat_mul(&x, &cprime, n);
    (eps, c)
}

fn electron_partition(electrons: usize, multiplicity: u8) -> anyhow::Result<(usize, usize)> {
    if multiplicity == 0 {
        return Err(anyhow::anyhow!("spin multiplicity must be >= 1"));
    }
    let delta = (multiplicity - 1) as usize; // n_alpha - n_beta
    if delta > electrons {
        return Err(anyhow::anyhow!(
            "invalid multiplicity {} for electron count {}",
            multiplicity,
            electrons
        ));
    }
    if (electrons + delta) % 2 != 0 {
        return Err(anyhow::anyhow!(
            "parity mismatch for multiplicity {} and electron count {}",
            multiplicity,
            electrons
        ));
    }
    let n_alpha = (electrons + delta) / 2;
    let n_beta = electrons - n_alpha;
    Ok((n_alpha, n_beta))
}

fn mulliken_and_dipole(
    mol: &MoleculeInput,
    basis: &[BasisFn],
    s: &[f64],
    p_tot: &[f64],
) -> (Vec<f64>, f64) {
    let n = basis.len();
    let mut atom_pop = vec![0.0; mol.atoms.len()];
    for mu in 0..n {
        let a = basis[mu].atom_index;
        for nu in 0..n {
            atom_pop[a] += 0.5 * p_tot[idx(mu, nu, n)] * s[idx(mu, nu, n)];
        }
    }
    let mut mulliken = vec![0.0; mol.atoms.len()];
    for (aidx, a) in mol.atoms.iter().enumerate() {
        mulliken[aidx] = a.z as f64 - atom_pop[aidx];
    }

    let mut mu_bohr = [0.0_f64; 3];
    for (aidx, a) in mol.atoms.iter().enumerate() {
        let r = [
            a.x_ang * BOHR_PER_ANGSTROM,
            a.y_ang * BOHR_PER_ANGSTROM,
            a.z_ang * BOHR_PER_ANGSTROM,
        ];
        mu_bohr[0] += mulliken[aidx] * r[0];
        mu_bohr[1] += mulliken[aidx] * r[1];
        mu_bohr[2] += mulliken[aidx] * r[2];
    }
    let dipole =
        (mu_bohr[0].powi(2) + mu_bohr[1].powi(2) + mu_bohr[2].powi(2)).sqrt() * DEBYE_PER_EA0;
    (mulliken, dipole)
}

fn valence_electrons(z: u16) -> usize {
    if z == 1 {
        return 1;
    }
    if z == 2 {
        return 2;
    }
    let core = if z <= 10 {
        2
    } else if z <= 18 {
        10
    } else if z <= 36 {
        18
    } else if z <= 54 {
        36
    } else if z <= 86 {
        54
    } else {
        86
    };
    (z.saturating_sub(core)).clamp(1, 8) as usize
}

fn basis_for_molecule(mol: &MoleculeInput) -> (Vec<BasisFn>, usize) {
    let mut basis = Vec::new();
    let mut electrons = 0usize;

    for (ai, a) in mol.atoms.iter().enumerate() {
        let val = valence_electrons(a.z);
        electrons += val;

        let nbf = if a.z <= 2 {
            1
        } else if val <= 2 {
            2
        } else {
            4
        };

        let z_eff = (a.z as f64).min(20.0);
        let alpha_base = (0.24 * z_eff.powf(1.18) + 0.10).clamp(0.18, 12.0);
        for k in 0..nbf {
            let alpha = (alpha_base * 0.42_f64.powi(k as i32)).max(0.06);
            basis.push(BasisFn {
                atom_index: ai,
                alpha,
                center_bohr: [
                    a.x_ang * BOHR_PER_ANGSTROM,
                    a.y_ang * BOHR_PER_ANGSTROM,
                    a.z_ang * BOHR_PER_ANGSTROM,
                ],
            });
        }
    }

    let e = (electrons as i32 - mol.charge).max(0) as usize;
    (basis, e)
}

fn build_integrals(mol: &MoleculeInput, basis: &[BasisFn]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, f64) {
    let n = basis.len();
    let mut s = vec![0.0; n * n];
    let mut t = vec![0.0; n * n];
    let mut v = vec![0.0; n * n];
    let mut eri = vec![0.0; n * n * n * n];

    let atom_centers: Vec<[f64; 3]> = mol
        .atoms
        .iter()
        .map(|a| [
            a.x_ang * BOHR_PER_ANGSTROM,
            a.y_ang * BOHR_PER_ANGSTROM,
            a.z_ang * BOHR_PER_ANGSTROM,
        ])
        .collect();

    for i in 0..n {
        for j in 0..=i {
            let bi = &basis[i];
            let bj = &basis[j];
            let sij = overlap_ss(bi.alpha, bj.alpha, bi.center_bohr, bj.center_bohr);
            let tij = kinetic_ss(bi.alpha, bj.alpha, bi.center_bohr, bj.center_bohr);
            let mut vij = 0.0;
            for (aidx, a) in mol.atoms.iter().enumerate() {
                vij += nuclear_attraction_ss(
                    bi.alpha,
                    bj.alpha,
                    bi.center_bohr,
                    bj.center_bohr,
                    atom_centers[aidx],
                    a.z as f64,
                );
            }
            s[idx(i, j, n)] = sij;
            s[idx(j, i, n)] = sij;
            t[idx(i, j, n)] = tij;
            t[idx(j, i, n)] = tij;
            v[idx(i, j, n)] = vij;
            v[idx(j, i, n)] = vij;
        }
    }

    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                for l in 0..n {
                    let vijkl = eri_ssss(
                        basis[i].alpha,
                        basis[j].alpha,
                        basis[k].alpha,
                        basis[l].alpha,
                        basis[i].center_bohr,
                        basis[j].center_bohr,
                        basis[k].center_bohr,
                        basis[l].center_bohr,
                    );
                    eri[idx4(i, j, k, l, n)] = vijkl;
                }
            }
        }
    }

    let mut e_nuc = 0.0;
    for a in 0..mol.atoms.len() {
        for b in (a + 1)..mol.atoms.len() {
            let ra = atom_centers[a];
            let rb = atom_centers[b];
            let r = dist2(ra, rb).sqrt().max(1.0e-8);
            e_nuc += (mol.atoms[a].z as f64) * (mol.atoms[b].z as f64) / r;
        }
    }

    (s, t, v, eri, e_nuc)
}

fn rhf_scf(mol: &MoleculeInput, basis: &[BasisFn], electrons: usize, s: &[f64], h: &[f64], eri: &[f64], e_nuc: f64) -> MolecularAbInitioResult {
    let n = basis.len();
    let nocc = electrons / 2;

    let mut p = vec![0.0; n * n];
    let mut f = h.to_vec();
    let mut c = vec![0.0; n * n];
    let mut eps = vec![0.0; n];

    let mut e_total = 0.0;
    let mut e_elec = 0.0;
    let mut residual = f64::INFINITY;
    let mut iter_count = 0;

    for iter in 0..200 {
        iter_count = iter + 1;

        for i in 0..n {
            for j in 0..n {
                let mut gij = 0.0;
                for k in 0..n {
                    for l in 0..n {
                        let pkl = p[idx(k, l, n)];
                        let coul = eri[idx4(i, j, k, l, n)];
                        let exch = eri[idx4(i, k, j, l, n)];
                        gij += pkl * (coul - 0.5 * exch);
                    }
                }
                f[idx(i, j, n)] = h[idx(i, j, n)] + gij;
            }
        }

        let (eps_new, c_new) = generalized_eigensolve(&f, s, n);
        eps = eps_new;
        c = c_new;

        let mut p_new = vec![0.0; n * n];
        for mu in 0..n {
            for nu in 0..n {
                let mut sum = 0.0;
                for m in 0..nocc {
                    sum += c[idx(mu, m, n)] * c[idx(nu, m, n)];
                }
                p_new[idx(mu, nu, n)] = 2.0 * sum;
            }
        }

        let mut de = 0.0;
        e_elec = 0.0;
        for i in 0..n {
            for j in 0..n {
                let hij = h[idx(i, j, n)];
                let fij = f[idx(i, j, n)];
                e_elec += 0.5 * p_new[idx(i, j, n)] * (hij + fij);
                let dp = p_new[idx(i, j, n)] - p[idx(i, j, n)];
                de += dp * dp;
            }
        }
        residual = de.sqrt();
        e_total = e_elec + e_nuc;

        p = p_new;
        if residual < 1.0e-8 {
            break;
        }
    }

    // MP2 correlation energy (closed shell)
    let mut emp2 = 0.0;
    if nocc < n {
        for i in 0..nocc {
            for j in 0..nocc {
                for a in nocc..n {
                    for b in nocc..n {
                        let mut iajb = 0.0;
                        let mut ibja = 0.0;
                        for mu in 0..n {
                            let cmi = c[idx(mu, i, n)];
                            for nu in 0..n {
                                let cna = c[idx(nu, a, n)];
                                let cnb = c[idx(nu, b, n)];
                                for lam in 0..n {
                                    let clj = c[idx(lam, j, n)];
                                    for sig in 0..n {
                                        let csb = c[idx(sig, b, n)];
                                        let csa = c[idx(sig, a, n)];
                                        let eri_mnls = eri[idx4(mu, nu, lam, sig, n)];
                                        iajb += cmi * cna * clj * csb * eri_mnls;
                                        ibja += cmi * cnb * clj * csa * eri_mnls;
                                    }
                                }
                            }
                        }
                        let denom = eps[i] + eps[j] - eps[a] - eps[b];
                        if denom.abs() > 1.0e-10 {
                            emp2 += iajb * (2.0 * iajb - ibja) / denom;
                        }
                    }
                }
            }
        }
    }

    let (mulliken, dipole) = mulliken_and_dipole(mol, basis, s, &p);

    let homo = if nocc > 0 { eps[nocc - 1] * HARTREE_TO_EV } else { f64::NAN };
    let lumo = if nocc < n { eps[nocc] * HARTREE_TO_EV } else { f64::NAN };

    MolecularAbInitioResult {
        name: mol.name.clone(),
        method: "RHF+MP2".to_string(),
        spin_multiplicity: mol.multiplicity,
        basis_functions: n,
        electron_count: electrons,
        alpha_electrons: nocc,
        beta_electrons: nocc,
        electron_pairs: nocc,
        scf_iterations: iter_count,
        scf_residual: residual,
        s2_expectation: 0.0,
        nuclear_repulsion_hartree: e_nuc,
        electronic_energy_hartree: e_elec,
        total_energy_hartree: e_total,
        mp2_correlation_hartree: emp2,
        total_energy_mp2_hartree: e_total + emp2,
        homo_energy_ev: homo,
        lumo_energy_ev: lumo,
        homo_lumo_gap_ev: lumo - homo,
        dipole_debye: dipole,
        mulliken_charges: mulliken,
        orbital_energies_ev: eps.iter().map(|e| e * HARTREE_TO_EV).collect(),
    }
}

fn uhf_scf(
    mol: &MoleculeInput,
    basis: &[BasisFn],
    n_alpha: usize,
    n_beta: usize,
    s: &[f64],
    h: &[f64],
    eri: &[f64],
    e_nuc: f64,
) -> MolecularAbInitioResult {
    let n = basis.len();
    let mut p_alpha = vec![0.0; n * n];
    let mut p_beta = vec![0.0; n * n];
    let mut f_alpha = h.to_vec();
    let mut f_beta = h.to_vec();

    let mut c_alpha = vec![0.0; n * n];
    let mut c_beta = vec![0.0; n * n];
    let mut eps_alpha = vec![0.0; n];
    let mut eps_beta = vec![0.0; n];

    let mut e_elec = 0.0;
    let mut e_total = 0.0;
    let mut residual = f64::INFINITY;
    let mut iter_count = 0;

    for iter in 0..240 {
        iter_count = iter + 1;
        let mut p_tot = vec![0.0; n * n];
        for i in 0..n * n {
            p_tot[i] = p_alpha[i] + p_beta[i];
        }

        for i in 0..n {
            for j in 0..n {
                let mut g_alpha = 0.0;
                let mut g_beta = 0.0;
                for k in 0..n {
                    for l in 0..n {
                        let pkl = p_tot[idx(k, l, n)];
                        let coul = eri[idx4(i, j, k, l, n)];
                        let exch = eri[idx4(i, k, j, l, n)];
                        g_alpha += pkl * coul - p_alpha[idx(k, l, n)] * exch;
                        g_beta += pkl * coul - p_beta[idx(k, l, n)] * exch;
                    }
                }
                f_alpha[idx(i, j, n)] = h[idx(i, j, n)] + g_alpha;
                f_beta[idx(i, j, n)] = h[idx(i, j, n)] + g_beta;
            }
        }

        let (eps_a_new, c_a_new) = generalized_eigensolve(&f_alpha, s, n);
        let (eps_b_new, c_b_new) = generalized_eigensolve(&f_beta, s, n);
        eps_alpha = eps_a_new;
        eps_beta = eps_b_new;
        c_alpha = c_a_new;
        c_beta = c_b_new;

        let mut p_alpha_new = vec![0.0; n * n];
        let mut p_beta_new = vec![0.0; n * n];
        for mu in 0..n {
            for nu in 0..n {
                let mut sa = 0.0;
                for m in 0..n_alpha {
                    sa += c_alpha[idx(mu, m, n)] * c_alpha[idx(nu, m, n)];
                }
                let mut sb = 0.0;
                for m in 0..n_beta {
                    sb += c_beta[idx(mu, m, n)] * c_beta[idx(nu, m, n)];
                }
                p_alpha_new[idx(mu, nu, n)] = sa;
                p_beta_new[idx(mu, nu, n)] = sb;
            }
        }

        let mut de = 0.0;
        e_elec = 0.0;
        for i in 0..n {
            for j in 0..n {
                let p_tot_new = p_alpha_new[idx(i, j, n)] + p_beta_new[idx(i, j, n)];
                e_elec += 0.5
                    * (p_tot_new * h[idx(i, j, n)]
                        + p_alpha_new[idx(i, j, n)] * f_alpha[idx(i, j, n)]
                        + p_beta_new[idx(i, j, n)] * f_beta[idx(i, j, n)]);

                let da = p_alpha_new[idx(i, j, n)] - p_alpha[idx(i, j, n)];
                let db = p_beta_new[idx(i, j, n)] - p_beta[idx(i, j, n)];
                de += da * da + db * db;
            }
        }
        residual = de.sqrt();
        e_total = e_elec + e_nuc;

        p_alpha = p_alpha_new;
        p_beta = p_beta_new;

        if residual < 1.0e-8 {
            break;
        }
    }

    let mut p_tot = vec![0.0; n * n];
    for i in 0..n * n {
        p_tot[i] = p_alpha[i] + p_beta[i];
    }
    let (mulliken, dipole) = mulliken_and_dipole(mol, basis, s, &p_tot);

    let homo_a = if n_alpha > 0 {
        eps_alpha[n_alpha - 1] * HARTREE_TO_EV
    } else {
        f64::NEG_INFINITY
    };
    let homo_b = if n_beta > 0 {
        eps_beta[n_beta - 1] * HARTREE_TO_EV
    } else {
        f64::NEG_INFINITY
    };
    let homo = homo_a.max(homo_b);
    let lumo_a = if n_alpha < n {
        eps_alpha[n_alpha] * HARTREE_TO_EV
    } else {
        f64::INFINITY
    };
    let lumo_b = if n_beta < n {
        eps_beta[n_beta] * HARTREE_TO_EV
    } else {
        f64::INFINITY
    };
    let lumo = lumo_a.min(lumo_b);

    // UHF spin contamination estimate:
    // <S^2> = Sz(Sz+1) + Nβ - Σ_ij | <φα_i | φβ_j> |^2
    let sz = 0.5 * (n_alpha as f64 - n_beta as f64);
    let mut overlap_sum_sq = 0.0;
    for i in 0..n_alpha {
        for j in 0..n_beta {
            let mut ov = 0.0;
            for mu in 0..n {
                for nu in 0..n {
                    ov += c_alpha[idx(mu, i, n)] * s[idx(mu, nu, n)] * c_beta[idx(nu, j, n)];
                }
            }
            overlap_sum_sq += ov * ov;
        }
    }
    let s2 = sz * (sz + 1.0) + n_beta as f64 - overlap_sum_sq;

    let mut orb_energies = eps_alpha
        .iter()
        .map(|e| e * HARTREE_TO_EV)
        .collect::<Vec<_>>();
    orb_energies.extend(eps_beta.iter().map(|e| e * HARTREE_TO_EV));

    MolecularAbInitioResult {
        name: mol.name.clone(),
        method: "UHF".to_string(),
        spin_multiplicity: mol.multiplicity,
        basis_functions: n,
        electron_count: n_alpha + n_beta,
        alpha_electrons: n_alpha,
        beta_electrons: n_beta,
        electron_pairs: n_beta,
        scf_iterations: iter_count,
        scf_residual: residual,
        s2_expectation: s2,
        nuclear_repulsion_hartree: e_nuc,
        electronic_energy_hartree: e_elec,
        total_energy_hartree: e_total,
        mp2_correlation_hartree: f64::NAN,
        total_energy_mp2_hartree: f64::NAN,
        homo_energy_ev: homo,
        lumo_energy_ev: lumo,
        homo_lumo_gap_ev: lumo - homo,
        dipole_debye: dipole,
        mulliken_charges: mulliken,
        orbital_energies_ev: orb_energies,
    }
}

pub fn run_molecular_ab_initio(mol: MoleculeInput) -> anyhow::Result<MolecularAbInitioResult> {
    let (basis, electrons) = basis_for_molecule(&mol);
    if basis.is_empty() {
        return Err(anyhow::anyhow!("molecule has no basis functions"));
    }
    if electrons == 0 {
        return Err(anyhow::anyhow!("molecule has zero electrons"));
    }
    let (n_alpha, n_beta) = electron_partition(electrons, mol.multiplicity)?;
    if n_alpha > basis.len() || n_beta > basis.len() {
        return Err(anyhow::anyhow!(
            "insufficient basis for spin occupation: n_alpha={}, n_beta={}, basis={}",
            n_alpha,
            n_beta,
            basis.len(),
        ));
    }

    let n = basis.len();
    let (s, t, v, eri, e_nuc) = build_integrals(&mol, &basis);
    let mut h = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            h[idx(i, j, n)] = t[idx(i, j, n)] + v[idx(i, j, n)];
        }
    }

    if mol.multiplicity == 1 && electrons % 2 == 0 {
        Ok(rhf_scf(&mol, &basis, electrons, &s, &h, &eri, e_nuc))
    } else {
        Ok(uhf_scf(
            &mol, &basis, n_alpha, n_beta, &s, &h, &eri, e_nuc,
        ))
    }
}

pub fn optimize_molecule_geometry(
    mol: MoleculeInput,
    max_iter: usize,
    step_size: f64,
) -> anyhow::Result<GeometryOptimizationResult> {
    if mol.atoms.len() <= 1 {
        let r = run_molecular_ab_initio(mol.clone())?;
        return Ok(GeometryOptimizationResult {
            optimized_molecule: mol,
            final_result: r,
            iterations: 0,
            converged: true,
            final_gradient_norm_hartree_per_angstrom: 0.0,
        });
    }

    let mut current = mol.clone();
    let mut result = run_molecular_ab_initio(current.clone())?;
    let mut converged = false;
    let mut final_grad_norm = f64::INFINITY;
    let h = 1.0e-3;
    let max_step = 0.05;
    let alpha0 = step_size.max(1.0e-4);

    for it in 0..max_iter {
        let mut grad = vec![[0.0_f64; 3]; current.atoms.len()];
        for a in 1..current.atoms.len() {
            for c in 0..3 {
                let mut plus = current.clone();
                let mut minus = current.clone();
                match c {
                    0 => {
                        plus.atoms[a].x_ang += h;
                        minus.atoms[a].x_ang -= h;
                    }
                    1 => {
                        plus.atoms[a].y_ang += h;
                        minus.atoms[a].y_ang -= h;
                    }
                    _ => {
                        plus.atoms[a].z_ang += h;
                        minus.atoms[a].z_ang -= h;
                    }
                }
                let e_plus = run_molecular_ab_initio(plus)?.total_energy_hartree;
                let e_minus = run_molecular_ab_initio(minus)?.total_energy_hartree;
                grad[a][c] = (e_plus - e_minus) / (2.0 * h);
            }
        }

        let mut gn2 = 0.0;
        for g in &grad {
            gn2 += g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
        }
        final_grad_norm = gn2.sqrt();
        if final_grad_norm < 1.0e-3 {
            converged = true;
            return Ok(GeometryOptimizationResult {
                optimized_molecule: current,
                final_result: result,
                iterations: it + 1,
                converged,
                final_gradient_norm_hartree_per_angstrom: final_grad_norm,
            });
        }

        let mut trial = current.clone();
        let mut alpha = alpha0;
        let e_current = result.total_energy_hartree;
        let mut accepted = false;
        for _ in 0..6 {
            for a in 1..trial.atoms.len() {
                let dx = (-alpha * grad[a][0]).clamp(-max_step, max_step);
                let dy = (-alpha * grad[a][1]).clamp(-max_step, max_step);
                let dz = (-alpha * grad[a][2]).clamp(-max_step, max_step);
                trial.atoms[a].x_ang = current.atoms[a].x_ang + dx;
                trial.atoms[a].y_ang = current.atoms[a].y_ang + dy;
                trial.atoms[a].z_ang = current.atoms[a].z_ang + dz;
            }
            let trial_result = run_molecular_ab_initio(trial.clone())?;
            if trial_result.total_energy_hartree <= e_current {
                current = trial;
                result = trial_result;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if !accepted {
            break;
        }
    }

    Ok(GeometryOptimizationResult {
        optimized_molecule: current,
        final_result: result,
        iterations: max_iter,
        converged,
        final_gradient_norm_hartree_per_angstrom: final_grad_norm,
    })
}

pub fn benchmark_molecules() -> Vec<MoleculeInput> {
    vec![
        MoleculeInput {
            name: "H2".to_string(),
            atoms: vec![
                Atom3D {
                    z: 1,
                    x_ang: -0.37,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
                Atom3D {
                    z: 1,
                    x_ang: 0.37,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
            ],
            charge: 0,
            multiplicity: 1,
        },
        MoleculeInput {
            name: "LiH".to_string(),
            atoms: vec![
                Atom3D {
                    z: 3,
                    x_ang: -0.80,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
                Atom3D {
                    z: 1,
                    x_ang: 0.80,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
            ],
            charge: 0,
            multiplicity: 1,
        },
        MoleculeInput {
            name: "HF".to_string(),
            atoms: vec![
                Atom3D {
                    z: 1,
                    x_ang: -0.458,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
                Atom3D {
                    z: 9,
                    x_ang: 0.458,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
            ],
            charge: 0,
            multiplicity: 1,
        },
        MoleculeInput {
            name: "H2O".to_string(),
            atoms: vec![
                Atom3D {
                    z: 8,
                    x_ang: 0.0,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
                Atom3D {
                    z: 1,
                    x_ang: 0.757,
                    y_ang: 0.586,
                    z_ang: 0.0,
                },
                Atom3D {
                    z: 1,
                    x_ang: -0.757,
                    y_ang: 0.586,
                    z_ang: 0.0,
                },
            ],
            charge: 0,
            multiplicity: 1,
        },
        MoleculeInput {
            name: "NH3".to_string(),
            atoms: vec![
                Atom3D {
                    z: 7,
                    x_ang: 0.0,
                    y_ang: 0.0,
                    z_ang: 0.10,
                },
                Atom3D {
                    z: 1,
                    x_ang: 0.0,
                    y_ang: 0.94,
                    z_ang: -0.32,
                },
                Atom3D {
                    z: 1,
                    x_ang: 0.814,
                    y_ang: -0.47,
                    z_ang: -0.32,
                },
                Atom3D {
                    z: 1,
                    x_ang: -0.814,
                    y_ang: -0.47,
                    z_ang: -0.32,
                },
            ],
            charge: 0,
            multiplicity: 1,
        },
        MoleculeInput {
            name: "CH4".to_string(),
            atoms: vec![
                Atom3D {
                    z: 6,
                    x_ang: 0.0,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
                Atom3D {
                    z: 1,
                    x_ang: 0.629,
                    y_ang: 0.629,
                    z_ang: 0.629,
                },
                Atom3D {
                    z: 1,
                    x_ang: -0.629,
                    y_ang: -0.629,
                    z_ang: 0.629,
                },
                Atom3D {
                    z: 1,
                    x_ang: -0.629,
                    y_ang: 0.629,
                    z_ang: -0.629,
                },
                Atom3D {
                    z: 1,
                    x_ang: 0.629,
                    y_ang: -0.629,
                    z_ang: -0.629,
                },
            ],
            charge: 0,
            multiplicity: 1,
        },
        MoleculeInput {
            name: "N2".to_string(),
            atoms: vec![
                Atom3D {
                    z: 7,
                    x_ang: -0.55,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
                Atom3D {
                    z: 7,
                    x_ang: 0.55,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
            ],
            charge: 0,
            multiplicity: 1,
        },
        MoleculeInput {
            name: "CO2".to_string(),
            atoms: vec![
                Atom3D {
                    z: 8,
                    x_ang: -1.16,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
                Atom3D {
                    z: 6,
                    x_ang: 0.0,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
                Atom3D {
                    z: 8,
                    x_ang: 1.16,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
            ],
            charge: 0,
            multiplicity: 1,
        },
        MoleculeInput {
            name: "OH_radical".to_string(),
            atoms: vec![
                Atom3D {
                    z: 8,
                    x_ang: -0.485,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
                Atom3D {
                    z: 1,
                    x_ang: 0.485,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
            ],
            charge: 0,
            multiplicity: 2,
        },
        MoleculeInput {
            name: "NO_radical".to_string(),
            atoms: vec![
                Atom3D {
                    z: 7,
                    x_ang: -0.575,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
                Atom3D {
                    z: 8,
                    x_ang: 0.575,
                    y_ang: 0.0,
                    z_ang: 0.0,
                },
            ],
            charge: 0,
            multiplicity: 2,
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn h2_runs_and_binds() {
        let h2 = benchmark_molecules()
            .into_iter()
            .find(|m| m.name == "H2")
            .unwrap();
        let r = run_molecular_ab_initio(h2).unwrap();
        assert!(r.total_energy_hartree.is_finite());
        assert!(r.homo_lumo_gap_ev > 0.0);
    }

    #[test]
    fn water_has_nonzero_dipole_proxy() {
        let h2o = benchmark_molecules()
            .into_iter()
            .find(|m| m.name == "H2O")
            .unwrap();
        let r = run_molecular_ab_initio(h2o).unwrap();
        assert!(r.dipole_debye > 0.05);
    }

    #[test]
    fn open_shell_uhf_runs() {
        let oh = benchmark_molecules()
            .into_iter()
            .find(|m| m.name == "OH_radical")
            .unwrap();
        let r = run_molecular_ab_initio(oh).unwrap();
        assert_eq!(r.method, "UHF");
        assert!(r.s2_expectation > 0.1);
    }

    #[test]
    fn geometry_optimization_lowers_h2_energy() {
        let mut h2 = benchmark_molecules()
            .into_iter()
            .find(|m| m.name == "H2")
            .unwrap();
        h2.atoms[0].x_ang = -0.60;
        h2.atoms[1].x_ang = 0.60;
        let e0 = run_molecular_ab_initio(h2.clone()).unwrap().total_energy_hartree;
        let opt = optimize_molecule_geometry(h2, 8, 0.08).unwrap();
        let e1 = opt.final_result.total_energy_hartree;
        assert!(e1 <= e0 + 1.0e-8);
    }
}
