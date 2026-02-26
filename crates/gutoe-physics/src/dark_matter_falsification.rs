/*!
 * GUTOE Physics - Dark Matter Dataset Scoring + Falsification Gates
 * Copyright (C) 2026  Riff Labs
 *
 * Dataset source:
 *   SPARC Mass Models (Lelli+2016c), mirrored snapshot in
 *   `data/sparc_massmodels_2016c_baryon.csv`.
 */

use crate::constants::{C, DARK_TO_VISIBLE_COUNT_RATIO, DARK_TO_VISIBLE_GEOMETRIC_RATIO};
use crate::dark_sector::DarkSectorBranch;

/// Observed baryon and dark-matter density fractions (Planck-era baseline used
/// across the existing GRAND-346 harness).
pub const OMEGA_BARYON_OBS: f64 = 0.0493;
pub const OMEGA_DM_OBS: f64 = 0.264;
pub const OMEGA_MATTER_OBS: f64 = OMEGA_BARYON_OBS + OMEGA_DM_OBS;
pub const DM_FRACTION_OBS: f64 = OMEGA_DM_OBS / OMEGA_MATTER_OBS;

/// One SPARC mass-model data row (rotation curve + baryonic decomposition).
#[derive(Debug, Clone, PartialEq)]
pub struct SparcMassRow {
    pub galaxy: String,
    pub radius_kpc: f64,
    pub v_obs_kms: f64,
    pub e_vobs_kms: f64,
    pub v_gas_kms: f64,
    pub v_disk_kms: f64,
    pub v_bulge_kms: f64,
    pub v_baryon_kms: f64,
}

/// Aggregate fit metrics for one dark-sector branch over the SPARC snapshot.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DarkMatterBranchFitMetrics {
    pub n_points: usize,
    pub rotation_rmse_kms: f64,
    pub rotation_mape: f64,
    pub rotation_chi2_ndof: f64,
    pub lensing_proxy_rmse_rad: f64,
    pub lensing_proxy_mape: f64,
    pub predicted_dm_fraction: f64,
    pub observed_dm_fraction: f64,
    pub dm_fraction_delta: f64,
}

/// Explicit thresholds for the dark-matter falsification gate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DarkMatterFalsificationWindows {
    pub rotation_mape_max: f64,
    pub lensing_proxy_mape_max: f64,
    pub dm_fraction_delta_abs_max: f64,
}

impl Default for DarkMatterFalsificationWindows {
    fn default() -> Self {
        Self {
            // Keep thresholds explicit and conservative: this is a first-pass
            // dataset gate, not a tuned halo-fit pipeline.
            rotation_mape_max: 0.35,
            lensing_proxy_mape_max: 0.80,
            dm_fraction_delta_abs_max: 0.01,
        }
    }
}

/// Branch scorecard against dataset + CMB constraints.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DarkMatterBranchScorecard {
    pub branch: DarkSectorBranch,
    pub metrics: DarkMatterBranchFitMetrics,
    pub rotation_ok: bool,
    pub lensing_ok: bool,
    pub cmb_fraction_ok: bool,
}

impl DarkMatterBranchScorecard {
    pub const fn passes_all(self) -> bool {
        self.rotation_ok && self.lensing_ok && self.cmb_fraction_ok
    }
}

/// Parse the vendored SPARC CSV snapshot.
pub fn parse_sparc_massmodels_csv(csv_data: &str) -> Vec<SparcMassRow> {
    let mut out = Vec::new();
    for (idx, line) in csv_data.lines().enumerate() {
        if idx == 0 {
            continue;
        }
        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() != 8 {
            continue;
        }
        let parse = |s: &str| s.trim().parse::<f64>().ok();
        let row = match (
            parse(cols[1]),
            parse(cols[2]),
            parse(cols[3]),
            parse(cols[4]),
            parse(cols[5]),
            parse(cols[6]),
            parse(cols[7]),
        ) {
            (Some(radius_kpc), Some(v_obs_kms), Some(e_vobs_kms), Some(v_gas_kms), Some(v_disk_kms), Some(v_bulge_kms), Some(v_baryon_kms))
                if radius_kpc > 0.0 && v_obs_kms > 0.0 && v_baryon_kms > 0.0 =>
            {
                SparcMassRow {
                    galaxy: cols[0].trim().to_string(),
                    radius_kpc,
                    v_obs_kms,
                    e_vobs_kms,
                    v_gas_kms,
                    v_disk_kms,
                    v_bulge_kms,
                    v_baryon_kms,
                }
            }
            _ => continue,
        };
        out.push(row);
    }
    out
}

/// Load the vendored SPARC mass-model snapshot.
pub fn sparc_massmodels_dataset() -> Vec<SparcMassRow> {
    parse_sparc_massmodels_csv(include_str!("../data/sparc_massmodels_2016c_baryon.csv"))
}

/// Structural dark/visible mass ratio by branch.
pub fn branch_dark_to_visible_ratio(branch: DarkSectorBranch) -> f64 {
    match branch {
        DarkSectorBranch::Particle => DARK_TO_VISIBLE_COUNT_RATIO,
        DarkSectorBranch::Geometric => DARK_TO_VISIBLE_GEOMETRIC_RATIO,
    }
}

/// Structural total speed prediction from baryonic speed and branch ratio.
///
/// With `v^2 ~ GM/r`, scaling enclosed mass by `(1 + r_dm)` scales speed by
/// `sqrt(1 + r_dm)`.
pub fn predicted_speed_kms(branch: DarkSectorBranch, v_baryon_kms: f64) -> f64 {
    let ratio = branch_dark_to_visible_ratio(branch);
    (1.0 + ratio).sqrt() * v_baryon_kms
}

/// Lensing-deflection proxy from circular speed (weak-field SIS-like relation):
/// α ≈ 4(v/c)^2.
pub fn lensing_proxy_from_speed(v_kms: f64) -> f64 {
    let v_m_s = v_kms * 1.0e3;
    4.0 * (v_m_s / C).powi(2)
}

/// Branch-level fit metrics over the SPARC snapshot.
pub fn evaluate_branch_fit(branch: DarkSectorBranch, data: &[SparcMassRow]) -> DarkMatterBranchFitMetrics {
    if data.is_empty() {
        return DarkMatterBranchFitMetrics {
            n_points: 0,
            rotation_rmse_kms: f64::NAN,
            rotation_mape: f64::NAN,
            rotation_chi2_ndof: f64::NAN,
            lensing_proxy_rmse_rad: f64::NAN,
            lensing_proxy_mape: f64::NAN,
            predicted_dm_fraction: f64::NAN,
            observed_dm_fraction: DM_FRACTION_OBS,
            dm_fraction_delta: f64::NAN,
        };
    }

    let mut sq_v = 0.0;
    let mut ape_v = 0.0;
    let mut chi2 = 0.0;
    let mut sq_alpha = 0.0;
    let mut ape_alpha = 0.0;

    for row in data {
        let v_pred = predicted_speed_kms(branch, row.v_baryon_kms);
        let dv = v_pred - row.v_obs_kms;
        sq_v += dv * dv;
        ape_v += (dv.abs() / row.v_obs_kms).min(1.0e6);

        if row.e_vobs_kms > 0.0 {
            let z = dv / row.e_vobs_kms;
            chi2 += z * z;
        }

        let alpha_obs = lensing_proxy_from_speed(row.v_obs_kms);
        let alpha_pred = lensing_proxy_from_speed(v_pred);
        let da = alpha_pred - alpha_obs;
        sq_alpha += da * da;
        ape_alpha += (da.abs() / alpha_obs).min(1.0e6);
    }

    let n = data.len() as f64;
    let ratio = branch_dark_to_visible_ratio(branch);
    let predicted_dm_fraction = ratio / (1.0 + ratio);
    let dm_fraction_delta = predicted_dm_fraction - DM_FRACTION_OBS;

    DarkMatterBranchFitMetrics {
        n_points: data.len(),
        rotation_rmse_kms: (sq_v / n).sqrt(),
        rotation_mape: ape_v / n,
        rotation_chi2_ndof: chi2 / (n - 1.0).max(1.0),
        lensing_proxy_rmse_rad: (sq_alpha / n).sqrt(),
        lensing_proxy_mape: ape_alpha / n,
        predicted_dm_fraction,
        observed_dm_fraction: DM_FRACTION_OBS,
        dm_fraction_delta,
    }
}

/// Evaluate one branch against explicit falsification windows.
pub fn evaluate_branch_scorecard(
    branch: DarkSectorBranch,
    data: &[SparcMassRow],
    windows: DarkMatterFalsificationWindows,
) -> DarkMatterBranchScorecard {
    let metrics = evaluate_branch_fit(branch, data);
    DarkMatterBranchScorecard {
        branch,
        rotation_ok: metrics.rotation_mape <= windows.rotation_mape_max,
        lensing_ok: metrics.lensing_proxy_mape <= windows.lensing_proxy_mape_max,
        cmb_fraction_ok: metrics.dm_fraction_delta.abs() <= windows.dm_fraction_delta_abs_max,
        metrics,
    }
}

/// Evaluate both currently-supported dark-sector branches.
pub fn evaluate_dark_matter_gate(
    windows: DarkMatterFalsificationWindows,
) -> [DarkMatterBranchScorecard; 2] {
    let data = sparc_massmodels_dataset();
    [
        evaluate_branch_scorecard(DarkSectorBranch::Particle, &data, windows),
        evaluate_branch_scorecard(DarkSectorBranch::Geometric, &data, windows),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{DARK_FRACTION_GEOMETRIC_STRUCTURAL, DARK_FRACTION_TOTAL_STATE_SPLIT};

    #[test]
    fn dataset_parse_has_expected_size_scale() {
        let data = sparc_massmodels_dataset();
        assert!(
            data.len() > 3000,
            "SPARC snapshot unexpectedly small: {}",
            data.len()
        );
    }

    #[test]
    fn structural_dm_fractions_match_lean_parity_values() {
        let p = branch_dark_to_visible_ratio(DarkSectorBranch::Particle);
        let g = branch_dark_to_visible_ratio(DarkSectorBranch::Geometric);
        assert!((p - 5.0 / 11.0).abs() < 1e-15);
        assert!((g - 60.0 / 11.0).abs() < 1e-15);
        assert!((DARK_FRACTION_TOTAL_STATE_SPLIT - 5.0 / 16.0).abs() < 1e-15);
        assert!((DARK_FRACTION_GEOMETRIC_STRUCTURAL - 60.0 / 71.0).abs() < 1e-15);
    }

    #[test]
    fn particle_branch_fits_rotation_better_than_geometric_branch_on_sparc_snapshot() {
        let data = sparc_massmodels_dataset();
        let m_p = evaluate_branch_fit(DarkSectorBranch::Particle, &data);
        let m_g = evaluate_branch_fit(DarkSectorBranch::Geometric, &data);
        assert!(m_p.rotation_mape < m_g.rotation_mape);
        assert!(m_p.lensing_proxy_mape < m_g.lensing_proxy_mape);
    }

    #[test]
    fn geometric_branch_hits_cmb_fraction_better_than_particle_branch() {
        let data = sparc_massmodels_dataset();
        let m_p = evaluate_branch_fit(DarkSectorBranch::Particle, &data);
        let m_g = evaluate_branch_fit(DarkSectorBranch::Geometric, &data);
        assert!(m_g.dm_fraction_delta.abs() < m_p.dm_fraction_delta.abs());
    }

    #[test]
    fn gate_exposes_particle_vs_geometric_tension() {
        let [particle, geometric] = evaluate_dark_matter_gate(DarkMatterFalsificationWindows::default());
        assert!(particle.rotation_ok && particle.lensing_ok);
        assert!(!particle.cmb_fraction_ok);
        assert!(!geometric.rotation_ok && !geometric.lensing_ok);
        assert!(geometric.cmb_fraction_ok);
    }
}
