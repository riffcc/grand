/*!
 * Coarse global nuclear chart scanner for Islands-of-Stability work.
 *
 * Phase-1 scope:
 * - Semi-empirical mass formula (SEMF) baseline.
 * - Lightweight shell-correction layer (magic-number attractors).
 * - Derived observables for falsification gates:
 *   - valley of stability by A
 *   - S2n / S2p separation energies
 *   - superheavy island ranking
 */

use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::Path;

/// Canonical magic numbers used by the shell correction layer.
pub const MAGIC_NUMBERS: [u16; 8] = [2, 8, 20, 28, 50, 82, 126, 184];

#[derive(Clone, Copy, Debug)]
pub struct SemfParams {
    pub a_v: f64,
    pub a_s: f64,
    pub a_c: f64,
    pub a_a: f64,
    pub a_p: f64,
}

impl Default for SemfParams {
    fn default() -> Self {
        Self {
            a_v: 15.8,
            a_s: 18.3,
            a_c: 0.714,
            a_a: 23.2,
            a_p: 12.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ShellParams {
    pub amplitude_z: f64,
    pub amplitude_n: f64,
    pub sigma_z: f64,
    pub sigma_n: f64,
}

impl Default for ShellParams {
    fn default() -> Self {
        Self {
            amplitude_z: 2.2,
            amplitude_n: 2.8,
            sigma_z: 4.0,
            sigma_n: 5.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ScanConfig {
    pub z_min: u16,
    pub z_max: u16,
    pub n_min: u16,
    pub n_max: u16,
    pub semf: SemfParams,
    pub shell: ShellParams,
}

impl Default for ScanConfig {
    fn default() -> Self {
        Self {
            z_min: 1,
            z_max: 140,
            n_min: 1,
            n_max: 260,
            semf: SemfParams::default(),
            shell: ShellParams::default(),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct NucleusRecord {
    pub z: u16,
    pub n: u16,
    pub a: u16,
    pub binding_mev: f64,
    pub binding_per_nucleon_mev: f64,
    pub shell_bonus_mev: f64,
    pub pairing_mev: f64,
    pub s2n_mev: Option<f64>,
    pub s2p_mev: Option<f64>,
    pub beta_optimal_for_a: bool,
    pub fissility: f64,
    pub stability_score: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct IslandRankingConfig {
    pub min_z: u16,
    pub max_fissility: f64,
    pub target_z: u16,
    pub target_n: u16,
    pub sigma_z: f64,
    pub sigma_n: f64,
    pub proximity_weight: f64,
    pub score_threshold: f64,
}

impl Default for IslandRankingConfig {
    fn default() -> Self {
        Self {
            min_z: 104,
            max_fissility: 1.1,
            target_z: 114,
            target_n: 184,
            sigma_z: 10.0,
            sigma_n: 18.0,
            proximity_weight: 0.35,
            score_threshold: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct MagicDiscontinuity {
    pub magic_n: u16,
    pub z: u16,
    pub delta_s2n_mev: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct ShellGateMetrics {
    pub top_delta_s2n_mev: f64,
    pub avg_top5_delta_s2n_mev: f64,
    pub strongest_n184_delta_s2n_mev: f64,
}

fn pairing_term(z: u16, n: u16, a: f64, a_p: f64) -> f64 {
    let z_even = z % 2 == 0;
    let n_even = n % 2 == 0;
    if z_even && n_even {
        a_p / a.sqrt()
    } else if !z_even && !n_even {
        -a_p / a.sqrt()
    } else {
        0.0
    }
}

fn shell_bonus(x: u16, magics: &[u16], amplitude: f64, sigma: f64) -> f64 {
    let xf = x as f64;
    magics
        .iter()
        .map(|&m| {
            let dx = xf - m as f64;
            amplitude * (-(dx * dx) / (2.0 * sigma * sigma)).exp()
        })
        .sum()
}

fn semf_binding_mev(z: u16, n: u16, semf: SemfParams, shell: ShellParams) -> (f64, f64, f64) {
    let a_u16 = z + n;
    let a = a_u16 as f64;
    let zf = z as f64;
    let n_asym = n as f64 - zf;

    let volume = semf.a_v * a;
    let surface = semf.a_s * a.powf(2.0 / 3.0);
    let coulomb = semf.a_c * zf * (zf - 1.0) / a.powf(1.0 / 3.0);
    let asymmetry = semf.a_a * n_asym * n_asym / a;
    let pairing = pairing_term(z, n, a, semf.a_p);

    let shell_z = shell_bonus(z, &MAGIC_NUMBERS, shell.amplitude_z, shell.sigma_z);
    let shell_n = shell_bonus(n, &MAGIC_NUMBERS, shell.amplitude_n, shell.sigma_n);
    let shell_total = shell_z + shell_n;

    let binding = volume - surface - coulomb - asymmetry + pairing + shell_total;
    (binding, shell_total, pairing)
}

fn fissility(z: u16, a: u16) -> f64 {
    let zf = z as f64;
    let af = a as f64;
    zf * zf / af / 50.0
}

fn stability_score(binding_per_nucleon: f64, s2n: Option<f64>, s2p: Option<f64>, fissility: f64) -> f64 {
    let s2n_term = s2n.unwrap_or(-10.0).clamp(-10.0, 20.0);
    let s2p_term = s2p.unwrap_or(-10.0).clamp(-10.0, 20.0);
    let fissility_penalty = if fissility > 1.0 { (fissility - 1.0) * 2.5 } else { 0.0 };
    binding_per_nucleon + 0.02 * s2n_term + 0.02 * s2p_term - fissility_penalty
}

/// Build a full nuclide table with SEMF+shell observables.
pub fn scan_nuclear_chart(cfg: ScanConfig) -> Vec<NucleusRecord> {
    let mut binding_map: BTreeMap<(u16, u16), (f64, f64, f64)> = BTreeMap::new();
    for z in cfg.z_min..=cfg.z_max {
        for n in cfg.n_min..=cfg.n_max {
            let (binding, shell_bonus_mev, pairing_mev) = semf_binding_mev(z, n, cfg.semf, cfg.shell);
            binding_map.insert((z, n), (binding, shell_bonus_mev, pairing_mev));
        }
    }

    let mut best_by_a: BTreeMap<u16, (u16, f64)> = BTreeMap::new();
    for (&(z, n), &(b, _, _)) in &binding_map {
        let a = z + n;
        match best_by_a.get(&a) {
            Some((_, best_b)) if *best_b >= b => {}
            _ => {
                best_by_a.insert(a, (z, b));
            }
        }
    }

    let mut out = Vec::with_capacity(binding_map.len());
    for (&(z, n), &(binding, shell_bonus_mev, pairing_mev)) in &binding_map {
        let a = z + n;
        let s2n = if n >= cfg.n_min + 2 {
            binding_map.get(&(z, n - 2)).map(|(b_prev, _, _)| binding - *b_prev)
        } else {
            None
        };
        let s2p = if z >= cfg.z_min + 2 {
            binding_map.get(&(z - 2, n)).map(|(b_prev, _, _)| binding - *b_prev)
        } else {
            None
        };
        let f = fissility(z, a);
        let bpa = binding / a as f64;
        let score = stability_score(bpa, s2n, s2p, f);
        let beta_optimal_for_a = best_by_a
            .get(&a)
            .map(|(best_z, _)| *best_z == z)
            .unwrap_or(false);

        out.push(NucleusRecord {
            z,
            n,
            a,
            binding_mev: binding,
            binding_per_nucleon_mev: bpa,
            shell_bonus_mev,
            pairing_mev,
            s2n_mev: s2n,
            s2p_mev: s2p,
            beta_optimal_for_a,
            fissility: f,
            stability_score: score,
        });
    }

    out.sort_by_key(|r| (r.z, r.n));
    out
}

/// Return strongest neutron shell closure signatures around magic N.
pub fn magic_s2n_discontinuities(records: &[NucleusRecord], top_k: usize) -> Vec<MagicDiscontinuity> {
    let mut by_zn: BTreeMap<(u16, u16), &NucleusRecord> = BTreeMap::new();
    for r in records {
        by_zn.insert((r.z, r.n), r);
    }

    let mut out = Vec::new();
    for &magic_n in &MAGIC_NUMBERS {
        for r in records {
            if r.n == magic_n {
                let Some(s2n_here) = r.s2n_mev else {
                    continue;
                };
                let Some(next) = by_zn.get(&(r.z, magic_n + 2)) else {
                    continue;
                };
                let Some(s2n_next) = next.s2n_mev else {
                    continue;
                };
                out.push(MagicDiscontinuity {
                    magic_n,
                    z: r.z,
                    delta_s2n_mev: s2n_here - s2n_next,
                });
            }
        }
    }

    out.sort_by(|a, b| b.delta_s2n_mev.total_cmp(&a.delta_s2n_mev));
    out.truncate(top_k);
    out
}

fn gaussian_proximity(x: f64, target: f64, sigma: f64) -> f64 {
    let d = x - target;
    (-(d * d) / (2.0 * sigma * sigma)).exp()
}

/// Rank superheavy candidates with an optional proximity pull toward a target `(Z,N)` island.
pub fn rank_island_candidates_with_config(
    records: &[NucleusRecord],
    cfg: IslandRankingConfig,
    top_k: usize,
) -> Vec<NucleusRecord> {
    let mut candidates: Vec<(f64, NucleusRecord)> = records
        .iter()
        .copied()
        .filter(|r| {
            r.z >= cfg.min_z
                && r.s2n_mev.unwrap_or(-1.0) > 0.0
                && r.s2p_mev.unwrap_or(-1.0) > 0.0
                && r.fissility < cfg.max_fissility
        })
        .filter_map(|r| {
            let proximity = gaussian_proximity(r.z as f64, cfg.target_z as f64, cfg.sigma_z)
                * gaussian_proximity(r.n as f64, cfg.target_n as f64, cfg.sigma_n);
            let score = r.stability_score + cfg.proximity_weight * proximity;
            if score >= cfg.score_threshold {
                Some((score, r))
            } else {
                None
            }
        })
        .collect();

    candidates.sort_by(|a, b| b.0.total_cmp(&a.0));
    candidates.truncate(top_k);
    candidates.into_iter().map(|(_, r)| r).collect()
}

/// Backwards-compatible ranker with default target pull toward the superheavy region.
pub fn rank_island_candidates(records: &[NucleusRecord], min_z: u16, top_k: usize) -> Vec<NucleusRecord> {
    let cfg = IslandRankingConfig {
        min_z,
        ..IslandRankingConfig::default()
    };
    rank_island_candidates_with_config(records, cfg, top_k)
}

/// Score shell/discontinuity quality for quick calibration loops.
pub fn shell_gate_metrics(records: &[NucleusRecord]) -> ShellGateMetrics {
    let rows = magic_s2n_discontinuities(records, records.len().min(200));
    if rows.is_empty() {
        return ShellGateMetrics {
            top_delta_s2n_mev: 0.0,
            avg_top5_delta_s2n_mev: 0.0,
            strongest_n184_delta_s2n_mev: 0.0,
        };
    }
    let top_delta = rows[0].delta_s2n_mev;
    let top5_len = rows.len().min(5);
    let avg_top5 = rows.iter().take(top5_len).map(|r| r.delta_s2n_mev).sum::<f64>() / top5_len as f64;
    let n184 = rows
        .iter()
        .filter(|r| r.magic_n == 184)
        .map(|r| r.delta_s2n_mev)
        .fold(0.0_f64, f64::max);
    ShellGateMetrics {
        top_delta_s2n_mev: top_delta,
        avg_top5_delta_s2n_mev: avg_top5,
        strongest_n184_delta_s2n_mev: n184,
    }
}

pub fn valley_of_stability(records: &[NucleusRecord]) -> Vec<NucleusRecord> {
    records.iter().copied().filter(|r| r.beta_optimal_for_a).collect()
}

pub fn closest_to_target_island(records: &[NucleusRecord], target_z: u16, target_n: u16) -> Option<NucleusRecord> {
    records
        .iter()
        .copied()
        .filter(|r| r.s2n_mev.unwrap_or(-1.0) > 0.0 && r.s2p_mev.unwrap_or(-1.0) > 0.0)
        .min_by(|a, b| {
            let da = ((a.z as i32 - target_z as i32).abs() + (a.n as i32 - target_n as i32).abs()) as i64;
            let db = ((b.z as i32 - target_z as i32).abs() + (b.n as i32 - target_n as i32).abs()) as i64;
            da.cmp(&db)
                .then_with(|| b.stability_score.total_cmp(&a.stability_score))
        })
}

/// Rank superheavy candidates by coarse stability score.
pub fn rank_island_candidates_legacy(records: &[NucleusRecord], min_z: u16, top_k: usize) -> Vec<NucleusRecord> {
    let mut candidates: Vec<NucleusRecord> = records
        .iter()
        .copied()
        .filter(|r| {
            r.z >= min_z
                && r.s2n_mev.unwrap_or(-1.0) > 0.0
                && r.s2p_mev.unwrap_or(-1.0) > 0.0
                && r.fissility < 1.1
        })
        .collect();

    candidates.sort_by(|a, b| b.stability_score.total_cmp(&a.stability_score));
    candidates.truncate(top_k);
    candidates
}

pub fn write_records_csv(path: impl AsRef<Path>, records: &[NucleusRecord]) -> std::io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(
        file,
        "Z,N,A,binding_mev,binding_per_nucleon_mev,shell_bonus_mev,pairing_mev,s2n_mev,s2p_mev,beta_optimal_for_a,fissility,stability_score"
    )?;
    for r in records {
        let s2n = r.s2n_mev.map(|v| format!("{v:.6}")).unwrap_or_default();
        let s2p = r.s2p_mev.map(|v| format!("{v:.6}")).unwrap_or_default();
        writeln!(
            file,
            "{},{},{},{:.6},{:.6},{:.6},{:.6},{},{},{},{:.6},{:.6}",
            r.z,
            r.n,
            r.a,
            r.binding_mev,
            r.binding_per_nucleon_mev,
            r.shell_bonus_mev,
            r.pairing_mev,
            s2n,
            s2p,
            r.beta_optimal_for_a,
            r.fissility,
            r.stability_score
        )?;
    }
    Ok(())
}

pub fn write_magic_discontinuities_csv(
    path: impl AsRef<Path>,
    rows: &[MagicDiscontinuity],
) -> std::io::Result<()> {
    let mut file = fs::File::create(path)?;
    writeln!(file, "magic_n,Z,delta_s2n_mev")?;
    for row in rows {
        writeln!(file, "{},{},{:.6}", row.magic_n, row.z, row.delta_s2n_mev)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iron_region_binding_is_reasonable() {
        let cfg = ScanConfig::default();
        let records = scan_nuclear_chart(cfg);
        let iron56 = records
            .iter()
            .find(|r| r.z == 26 && r.n == 30)
            .expect("Fe-56 must be in scan range");
        assert!(iron56.binding_per_nucleon_mev > 8.2);
        assert!(iron56.binding_per_nucleon_mev < 9.2);
    }

    #[test]
    fn shell_discontinuities_exist() {
        let cfg = ScanConfig::default();
        let records = scan_nuclear_chart(cfg);
        let rows = magic_s2n_discontinuities(&records, 20);
        assert!(!rows.is_empty());
        assert!(rows.iter().any(|r| r.delta_s2n_mev > 1.0));
    }

    #[test]
    fn ranking_with_config_returns_nonempty_superheavy_list() {
        let cfg = ScanConfig::default();
        let records = scan_nuclear_chart(cfg);
        let ranked = rank_island_candidates_with_config(&records, IslandRankingConfig::default(), 20);
        assert!(!ranked.is_empty());
        assert!(ranked.iter().all(|r| r.z >= 104));
    }
}
