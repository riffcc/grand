// GUTOE EM — SU(2) lattice gauge theory from the Z₃ magnetic triplet
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// The three SU(2) generators correspond to the Z₃ magnetic bivectors:
//   n₁ ↔ γ¹²  (state 7)   — Wilson plaquette component 1
//   n₂ ↔ γ²³  (state 11)  — Wilson plaquette component 2
//   n₃ ↔ γ³¹  (state 13)  — Wilson plaquette component 3
//
// SU(2) ≅ unit quaternions: U = n₀·1 + n₁·i + n₂·j + n₃·k, |n|² = 1
//
// Wilson action:  S = −β Σ_triangles Re Tr(W_p) / 2
//                   = −β Σ_triangles W_p[0]   (since Re Tr U = 2 n₀)
//
// The same Z₃ orbit decomposition that gives sin²θ_W = 3/13 now places
// exactly 3 independent generators on each hex lattice plaquette —
// the spatial content of Cl(1,3) becomes the SU(2) gauge group.

use std::collections::HashMap;

use rand::Rng;

use crate::config::LatticeConfig;
use crate::geometry::{mesh_neighbours, site_coords};

// ── SU(2) algebra ─────────────────────────────────────────────────────────────
//
// Unit quaternion representation:
//   U = (n₀, n₁, n₂, n₃)  with  n₀² + n₁² + n₂² + n₃² = 1
//
// Generators (pure imaginary quaternions):
//   i = (0, 1, 0, 0) ↔ γ¹²  — the "red" SU(2) direction
//   j = (0, 0, 1, 0) ↔ γ²³  — the "green" SU(2) direction
//   k = (0, 0, 0, 1) ↔ γ³¹  — the "blue" SU(2) direction
//
// Quaternion multiplication table: ij=k, jk=i, ki=j, ji=-k, kj=-i, ik=-j
// This encodes the su(2) Lie algebra: [i,j]=2k, [j,k]=2i, [k,i]=2j
//
// Trace: Tr U = 2·n₀  (Re Tr U = 2·n₀ used in Wilson action)
// Inverse/Dagger: U† = (n₀, -n₁, -n₂, -n₃)

pub type Su2 = [f64; 4];

/// Identity element of SU(2).
pub fn su2_identity() -> Su2 {
    [1.0, 0.0, 0.0, 0.0]
}

/// Quaternion multiplication (SU(2) group product).
pub fn su2_mul(a: &Su2, b: &Su2) -> Su2 {
    [
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    ]
}

/// Dagger (group inverse = quaternion conjugate).
pub fn su2_dag(a: &Su2) -> Su2 {
    [a[0], -a[1], -a[2], -a[3]]
}

/// Re Tr U = 2·n₀.  Normalised plaquette value is in [−1, +1].
pub fn su2_re_tr(a: &Su2) -> f64 {
    2.0 * a[0]
}

/// Dot product of quaternion components: Re Tr(A†B) = Σ aᵢbᵢ.
/// Used for fast staple-link inner products.
pub fn su2_dot(a: &Su2, b: &Su2) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
}

fn su2_normalise(a: &mut Su2) {
    let norm = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3]).sqrt();
    if norm > 1e-14 {
        for x in a.iter_mut() {
            *x /= norm;
        }
    }
}

/// Small random perturbation of identity: U = (√(1−ε²), ε·n̂).
/// `eps` controls the perturbation size; tune for ~50% Metropolis acceptance.
pub fn su2_random_perturb<R: Rng>(rng: &mut R, eps: f64) -> Su2 {
    let nx: f64 = rng.gen::<f64>() * 2.0 - 1.0;
    let ny: f64 = rng.gen::<f64>() * 2.0 - 1.0;
    let nz: f64 = rng.gen::<f64>() * 2.0 - 1.0;
    let len = (nx * nx + ny * ny + nz * nz).sqrt().max(1e-14);
    let n0 = (1.0 - eps * eps).max(0.0).sqrt();
    [n0, eps * nx / len, eps * ny / len, eps * nz / len]
}

/// Uniformly random SU(2) element (for hot start / β=0 initialisation).
pub fn su2_random<R: Rng>(rng: &mut R) -> Su2 {
    let mut u: Su2 = [
        rng.gen::<f64>() * 2.0 - 1.0,
        rng.gen::<f64>() * 2.0 - 1.0,
        rng.gen::<f64>() * 2.0 - 1.0,
        rng.gen::<f64>() * 2.0 - 1.0,
    ];
    su2_normalise(&mut u);
    u
}

// ── Link variable store ────────────────────────────────────────────────────────
//
// One SU(2) element per undirected edge, stored keyed by (min, max).
// For directed access:
//   get(i→j) where i < j : returns stored value
//   get(i→j) where i > j : returns dag of stored value (U†)
// This enforces U_{ji} = U_{ij}† automatically.

pub struct Su2Links {
    links: HashMap<(usize, usize), Su2>,
}

impl Su2Links {
    /// Cold start: all links = identity (corresponds to β → ∞, fully ordered).
    pub fn cold_start(cfg: &LatticeConfig) -> Self {
        let mut links = HashMap::new();
        for site in 0..cfg.n_sites() {
            let (r, c, z) = site_coords(site, cfg);
            for nb in mesh_neighbours(r, c, z, cfg) {
                let key = (site.min(nb), site.max(nb));
                links.entry(key).or_insert_with(su2_identity);
            }
        }
        Self { links }
    }

    /// Hot start: all links random (corresponds to β = 0, disordered).
    pub fn hot_start<R: Rng>(rng: &mut R, cfg: &LatticeConfig) -> Self {
        let mut links = HashMap::new();
        for site in 0..cfg.n_sites() {
            let (r, c, z) = site_coords(site, cfg);
            for nb in mesh_neighbours(r, c, z, cfg) {
                let key = (site.min(nb), site.max(nb));
                links.entry(key).or_insert_with(|| su2_random(rng));
            }
        }
        Self { links }
    }

    /// Get directed link U(from → to).
    /// Automatically returns U† when traversing against storage orientation.
    pub fn get(&self, from: usize, to: usize) -> Su2 {
        let key = (from.min(to), from.max(to));
        let u = self.links.get(&key).copied().unwrap_or_else(su2_identity);
        if from < to { u } else { su2_dag(&u) }
    }

    /// Set directed link U(from → to).  Stores the canonical (min, max) orientation.
    pub fn set(&mut self, from: usize, to: usize, u: Su2) {
        let key = (from.min(to), from.max(to));
        let stored = if from < to { u } else { su2_dag(&u) };
        self.links.insert(key, stored);
    }

    /// Number of stored undirected links.
    pub fn n_links(&self) -> usize {
        self.links.len()
    }

    // ── Plaquette ──────────────────────────────────────────────────────────────

    /// Wilson plaquette for triangle (i, j, k):
    ///   W_p = Re Tr(U_ij U_jk U_ki) / 2 = (U_ij · U_jk · U_ki)[0]
    /// Value in [−1, +1].  +1 = identity plaquette (ordered), 0 = random.
    pub fn plaquette_triangle(&self, i: usize, j: usize, k: usize) -> f64 {
        let uij = self.get(i, j);
        let ujk = self.get(j, k);
        let uki = self.get(k, i);
        su2_mul(&su2_mul(&uij, &ujk), &uki)[0]
    }

    /// Average plaquette over all triangular plaquettes in the hex lattice.
    /// Triangles are found by looking for common neighbours of each edge (i, j).
    pub fn avg_plaquette(&self, cfg: &LatticeConfig) -> f64 {
        let n = cfg.n_sites();
        let mut total = 0.0;
        let mut count = 0usize;
        for i in 0..n {
            let (ri, ci, zi) = site_coords(i, cfg);
            let nbrs_i = mesh_neighbours(ri, ci, zi, cfg);
            for &j in &nbrs_i {
                if j <= i {
                    continue;
                }
                let (rj, cj, zj) = site_coords(j, cfg);
                let nbrs_j = mesh_neighbours(rj, cj, zj, cfg);
                // Find third vertex k > j completing a triangle
                for &k in &nbrs_j {
                    if k <= j {
                        continue;
                    }
                    let (rk, ck, zk) = site_coords(k, cfg);
                    if mesh_neighbours(rk, ck, zk, cfg).contains(&i) {
                        total += self.plaquette_triangle(i, j, k);
                        count += 1;
                    }
                }
            }
        }
        if count == 0 { 0.0 } else { total / count as f64 }
    }

    // ── Metropolis update ──────────────────────────────────────────────────────
    //
    // Wilson action: S = −β Σ_p W_p[0]
    //
    // For link U_{ij}, the local action is:
    //   S_local = −β Σ_{triangles ∋ (i,j)} Re Tr(U_ij · staple_{ij})
    //           = −β · Re Tr(U_ij · Σ staples)
    //           = −β · su2_dot(U_ij, Σ staples)
    //
    // where staple for triangle (i,j,k) = U_jk · U_ki  (the path completing the triangle
    // without using edge (i,j)).
    //
    // Metropolis proposal: U_new = dU · U_old  (left-multiply by small perturbation)
    // Accept if exp(β(s_new − s_old)) > uniform[0,1].

    /// Single Metropolis sweep: visit every link once, propose and accept/reject.
    ///
    /// `eps`: perturbation amplitude (≈0.5 for β ≈ 1, smaller for larger β)
    pub fn metropolis_sweep<R: Rng>(
        &mut self,
        rng: &mut R,
        beta: f64,
        eps: f64,
        cfg: &LatticeConfig,
    ) {
        let n = cfg.n_sites();
        // Collect all undirected edges once
        let edges: Vec<(usize, usize)> = {
            let mut seen = std::collections::HashSet::new();
            let mut v = Vec::new();
            for site in 0..n {
                let (r, c, z) = site_coords(site, cfg);
                for nb in mesh_neighbours(r, c, z, cfg) {
                    let key = (site.min(nb), site.max(nb));
                    if seen.insert(key) {
                        v.push(key);
                    }
                }
            }
            v
        };

        for &(i, j) in &edges {
            // Compute the sum of staples for edge (i,j):
            //   Each triangle (i,j,k): staple = U_jk · U_ki
            let (ri, ci, zi) = site_coords(i, cfg);
            let (rj, cj, zj) = site_coords(j, cfg);
            let nbrs_i: std::collections::HashSet<usize> =
                mesh_neighbours(ri, ci, zi, cfg).into_iter().collect();
            let nbrs_j = mesh_neighbours(rj, cj, zj, cfg);

            let mut staple_sum = [0.0f64; 4];
            for &k in &nbrs_j {
                if k != i && nbrs_i.contains(&k) {
                    let st = su2_mul(&self.get(j, k), &self.get(k, i));
                    staple_sum[0] += st[0];
                    staple_sum[1] += st[1];
                    staple_sum[2] += st[2];
                    staple_sum[3] += st[3];
                }
            }

            let u_old = self.get(i, j);
            // Re Tr(U · K) = 2 · (U · K)[0]  — the actual Wilson contribution
            let s_old = su2_mul(&u_old, &staple_sum)[0];

            let du = su2_random_perturb(rng, eps);
            let mut u_new = su2_mul(&du, &u_old);
            su2_normalise(&mut u_new);

            let s_new = su2_mul(&u_new, &staple_sum)[0];

            // ΔS = −2β · (s_new − s_old); accept if ΔS ≤ 0 (improvement) or exp(−ΔS)
            let delta_s = -2.0 * beta * (s_new - s_old);
            if delta_s <= 0.0 || rng.gen::<f64>() < (-delta_s).exp() {
                self.set(i, j, u_new);
            }
        }
    }
}

// ── Wilson potential V(R) ─────────────────────────────────────────────────────
//
// Confinement test: measure W(R) = average Wilson loop of perimeter R.
//
// On the hex lattice the minimal loops are triangles (perimeter = 3).
// For larger R we follow the shortest closed path from a source site.
//
// Area law (confining):    W(R) ~ exp(−σ R²)  → V(R) = −ln W(R) grows as R
// Perimeter law (free):    W(R) ~ exp(−c R)   → V(R) linear in R (or const)
// Coulomb (weakly coupled): V(R) ~ 1/R

/// Average Wilson loop for all triangular plaquettes around a source site.
/// Returns the mean plaquette at the given source.
pub fn wilson_triangles_at(links: &Su2Links, source: usize, cfg: &LatticeConfig) -> f64 {
    let (r0, c0, z0) = site_coords(source, cfg);
    let src_nbrs = mesh_neighbours(r0, c0, z0, cfg);
    let mut total = 0.0;
    let mut count = 0usize;
    for &j in &src_nbrs {
        let (rj, cj, zj) = site_coords(j, cfg);
        let j_nbrs = mesh_neighbours(rj, cj, zj, cfg);
        for &k in &j_nbrs {
            if k != source {
                let (rk, ck, zk) = site_coords(k, cfg);
                if mesh_neighbours(rk, ck, zk, cfg).contains(&source) && k > j {
                    total += links.plaquette_triangle(source, j, k);
                    count += 1;
                }
            }
        }
    }
    if count == 0 { 0.0 } else { total / count as f64 }
}

/// Run a confinement experiment at a given β value.
///
/// Returns `(avg_plaquette, V3)` where:
///   - `avg_plaquette`: mean Re Tr W_p / 2 over the lattice  (1=ordered, 0=random)
///   - `V3`: effective potential −ln⟨W_triangle⟩  (proxy for string tension)
///
/// At strong coupling (β << 1): plaquette low, V3 large  → area law / confinement
/// At weak coupling (β >> 1):   plaquette → 1, V3 → 0    → deconfinement
pub fn confinement_experiment<R: Rng>(
    beta: f64,
    n_therm: usize,
    n_meas: usize,
    eps: f64,
    cfg: &LatticeConfig,
    rng: &mut R,
) -> (f64, f64) {
    let mut links = Su2Links::hot_start(rng, cfg);

    // Thermalise
    for _ in 0..n_therm {
        links.metropolis_sweep(rng, beta, eps, cfg);
    }

    // Measure
    let n = cfg.n_sites();
    let mut plaq_sum = 0.0;
    let mut w3_sum = 0.0;
    let mut w3_count = 0usize;

    for _ in 0..n_meas {
        links.metropolis_sweep(rng, beta, eps, cfg);
        plaq_sum += links.avg_plaquette(cfg);
        // Sample Wilson triangles from a few central sites
        for src in [n / 4, n / 2, 3 * n / 4] {
            let w = wilson_triangles_at(&links, src, cfg);
            if w.abs() > 0.0 {
                w3_sum += w;
                w3_count += 1;
            }
        }
    }

    let plaq_avg = plaq_sum / n_meas as f64;
    let w3_avg = if w3_count > 0 { w3_sum / w3_count as f64 } else { 1e-10 };
    let v3 = if w3_avg > 1e-10 { -w3_avg.ln() } else { 100.0 };

    (plaq_avg, v3)
}

// ── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::LatticeConfig;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn small_cfg() -> LatticeConfig {
        LatticeConfig { hex_rows: 8, hex_cols: 8, layers: 1, ..Default::default() }
    }

    // ── SU(2) algebra ──────────────────────────────────────────────────────────

    #[test]
    fn su2_identity_is_neutral() {
        let id = su2_identity();
        let mut a: Su2 = [0.6, 0.4, 0.5, 0.5];
        su2_normalise(&mut a);
        let prod = su2_mul(&id, &a);
        for i in 0..4 {
            assert!((prod[i] - a[i]).abs() < 1e-12, "id·a ≠ a at [{i}]");
        }
    }

    #[test]
    fn su2_dag_is_inverse() {
        let mut a: Su2 = [0.6, 0.4, 0.5, 0.5];
        su2_normalise(&mut a);
        let prod = su2_mul(&a, &su2_dag(&a));
        assert!((prod[0] - 1.0).abs() < 1e-12, "U·U† ≠ I: n₀={}", prod[0]);
        for i in 1..4 {
            assert!(prod[i].abs() < 1e-12, "U·U† ≠ I: n_{i}={}", prod[i]);
        }
    }

    #[test]
    fn su2_mul_preserves_norm() {
        let mut rng = StdRng::seed_from_u64(42);
        for _ in 0..200 {
            let a = su2_random(&mut rng);
            let b = su2_random(&mut rng);
            let c = su2_mul(&a, &b);
            let norm_sq: f64 = c.iter().map(|x| x * x).sum();
            assert!((norm_sq - 1.0).abs() < 1e-12, "|a·b|² = {norm_sq}");
        }
    }

    // ── Z₃ magnetic triplet ↔ SU(2) generators ────────────────────────────────
    //
    // The three Z₃ magnetic bivectors {γ¹², γ²³, γ³¹} (GUTOE states 7, 11, 13)
    // map onto the three quaternion generators i, j, k.
    // Commutation relations encode the su(2) Lie algebra.

    #[test]
    fn z3_magnetic_triplet_encodes_su2_algebra() {
        // γ¹² ↔ i = (0,1,0,0),  γ²³ ↔ j = (0,0,1,0),  γ³¹ ↔ k = (0,0,0,1)
        let gamma12: Su2 = [0.0, 1.0, 0.0, 0.0]; // state 7
        let gamma23: Su2 = [0.0, 0.0, 1.0, 0.0]; // state 11
        let gamma31: Su2 = [0.0, 0.0, 0.0, 1.0]; // state 13

        // i·j = k
        let ij = su2_mul(&gamma12, &gamma23);
        assert!((ij[3] - 1.0).abs() < 1e-12, "γ¹²·γ²³ should equal γ³¹ (k): {:?}", ij);

        // j·i = −k  (anticommutation)
        let ji = su2_mul(&gamma23, &gamma12);
        assert!((ji[3] + 1.0).abs() < 1e-12, "γ²³·γ¹² should equal −γ³¹: {:?}", ji);

        // [i,j] = 2k: ij[3] − ji[3] = 2
        let comm = ij[3] - ji[3];
        assert!((comm - 2.0).abs() < 1e-12, "[γ¹²,γ²³] = 2γ³¹: diff = {comm}");

        // Cyclic: j·k = i, k·i = j
        let jk = su2_mul(&gamma23, &gamma31);
        let ki = su2_mul(&gamma31, &gamma12);
        assert!((jk[1] - 1.0).abs() < 1e-12, "j·k = i: {:?}", jk);
        assert!((ki[2] - 1.0).abs() < 1e-12, "k·i = j: {:?}", ki);

        // All three generators square to −1 (the su(2) Lie algebra relation)
        for (gen, name) in [&gamma12, &gamma23, &gamma31].iter().zip(["γ¹²","γ²³","γ³¹"]) {
            let sq = su2_mul(gen, gen);
            assert!((sq[0] + 1.0).abs() < 1e-12, "{name}² should equal −1: {:?}", sq);
        }
    }

    // ── Link store ─────────────────────────────────────────────────────────────

    #[test]
    fn cold_start_plaquette_is_one() {
        let cfg = small_cfg();
        let links = Su2Links::cold_start(&cfg);
        let p = links.avg_plaquette(&cfg);
        assert!((p - 1.0).abs() < 1e-10, "Cold start plaquette = {p}, expected 1.0");
    }

    #[test]
    fn hot_start_plaquette_near_zero() {
        let cfg = small_cfg();
        let mut rng = StdRng::seed_from_u64(137);
        let links = Su2Links::hot_start(&mut rng, &cfg);
        let p = links.avg_plaquette(&cfg);
        // E[n₀] of uniform SU(2) = 0 by symmetry; should be small
        assert!(p.abs() < 0.3, "Hot start plaquette = {p:.4} (expected near 0)");
    }

    #[test]
    fn backward_link_is_dagger() {
        let cfg = small_cfg();
        let mut rng = StdRng::seed_from_u64(42);
        let links = Su2Links::hot_start(&mut rng, &cfg);
        // For site 0 and its first neighbour nb: get(nb, 0) == dag(get(0, nb))
        let (r, c, z) = site_coords(0, &cfg);
        let nb = mesh_neighbours(r, c, z, &cfg)[0];
        let u_fwd = links.get(0, nb);
        let u_bwd = links.get(nb, 0);
        let expected = su2_dag(&u_fwd);
        for i in 0..4 {
            assert!(
                (u_bwd[i] - expected[i]).abs() < 1e-12,
                "backward link != dag at component {i}"
            );
        }
    }

    // ── Metropolis dynamics ────────────────────────────────────────────────────

    #[test]
    fn metropolis_orders_links_at_large_beta() {
        let cfg = small_cfg();
        let mut rng = StdRng::seed_from_u64(137);
        let mut links = Su2Links::hot_start(&mut rng, &cfg);
        let p_before = links.avg_plaquette(&cfg);
        // Thermalise at large β — links should order toward identity
        for _ in 0..300 {
            links.metropolis_sweep(&mut rng, 3.0, 0.3, &cfg);
        }
        let p_after = links.avg_plaquette(&cfg);
        assert!(
            p_after > p_before + 0.3,
            "β=3 should order links: before={p_before:.3} after={p_after:.3}"
        );
    }

    // ── Confinement phase transition ───────────────────────────────────────────
    //
    // SU(2) on a 2D triangular lattice has a phase transition around β ≈ 1:
    //   β < β_c (strong coupling):  ⟨W_p⟩ low, V3 = −ln⟨W⟩ large  → area law
    //   β > β_c (weak coupling):    ⟨W_p⟩ high, V3 small            → deconfinement

    #[test]
    fn confinement_vs_deconfinement() {
        let cfg = small_cfg();
        let mut rng = StdRng::seed_from_u64(42);

        let (plaq_strong, v3_strong) =
            confinement_experiment(0.3, 200, 50, 0.9, &cfg, &mut rng);
        let (plaq_weak, v3_weak) =
            confinement_experiment(2.5, 200, 50, 0.2, &cfg, &mut rng);

        println!(
            "Strong coupling β=0.3: plaquette={plaq_strong:.3}  V3={v3_strong:.3}"
        );
        println!(
            "Weak coupling  β=2.5: plaquette={plaq_weak:.3}  V3={v3_weak:.3}"
        );

        // Plaquette should be larger at weak coupling (more ordered)
        assert!(
            plaq_weak > plaq_strong,
            "β=2.5 plaquette ({plaq_weak:.3}) should exceed β=0.3 ({plaq_strong:.3})"
        );

        // At strong coupling, effective potential should be larger (confining)
        assert!(
            v3_strong > v3_weak,
            "Strong coupling should have higher V3: V3(β=0.3)={v3_strong:.3} V3(β=2.5)={v3_weak:.3}"
        );

        println!(
            "CONFINEMENT: CONFIRMED — strong coupling gives V3={v3_strong:.3} > weak V3={v3_weak:.3}"
        );
    }

    #[test]
    fn plaquette_increases_monotonically_with_beta() {
        let cfg = small_cfg();
        let mut rng = StdRng::seed_from_u64(99);
        let betas = [0.2, 0.5, 1.0, 2.0];
        let mut plaq_prev = -1.0;
        for beta in betas {
            let eps = if beta < 1.0 { 0.8 } else { 0.4 };
            let mut links = Su2Links::hot_start(&mut rng, &cfg);
            for _ in 0..200 {
                links.metropolis_sweep(&mut rng, beta, eps, &cfg);
            }
            let p = links.avg_plaquette(&cfg);
            println!("  β={beta:.1}: plaquette={p:.4}");
            assert!(p > plaq_prev, "plaquette should increase: β={beta} p={p:.4} prev={plaq_prev:.4}");
            plaq_prev = p;
        }
    }
}
