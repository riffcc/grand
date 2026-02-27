/*!
 * GUTOE Physics - CLASS Full CMB Shape Harness
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-355:
 *   Run a full Boltzmann lane (CLASS) using derived GUTOE cosmology inputs
 *   and compare full TT spectrum shape against Planck binned data.
 */

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::Command;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClassRunInputs {
    pub h: f64,
    pub omega_b: f64,
    pub omega_cdm: f64,
    pub omega_k: f64,
    pub omega_lambda: f64,
    pub n_s: f64,
    pub a_s: f64,
    /// Reionization optical depth.
    /// NOTE: currently explicit assumption until derived in-framework.
    pub tau_reio: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlanckTtPoint {
    pub ell: u32,
    pub d_ell_tt_uk2: f64,
    pub sigma_uk2: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClassTtPoint {
    pub ell: u32,
    pub d_ell_tt_uk2: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TtResidual {
    pub ell: u32,
    pub observed_uk2: f64,
    pub predicted_uk2: f64,
    pub sigma_uk2: f64,
    pub pull: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ClassPlanckFit {
    pub n_points: usize,
    pub chi2: f64,
    pub reduced_chi2: f64,
    pub mean_abs_pull: f64,
    pub max_abs_pull: f64,
    pub rms_residual_uk2: f64,
    pub residuals: Vec<TtResidual>,
}

pub fn write_class_ini(
    ini_path: &Path,
    output_root: &str,
    lmax: u32,
    i: ClassRunInputs,
) -> Result<(), String> {
    let mut f =
        File::create(ini_path).map_err(|e| format!("create CLASS ini {:?}: {e}", ini_path))?;

    // CAMB format writes D_ell = l(l+1)C_ell/(2pi) in uK^2, matching the
    // common Planck binned TT products.
    writeln!(f, "h = {:.12}", i.h).map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "omega_b = {:.12}", i.omega_b).map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "omega_cdm = {:.12}", i.omega_cdm).map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "Omega_k = {:.12}", i.omega_k).map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "Omega_Lambda = {:.12}", i.omega_lambda).map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "A_s = {:.12e}", i.a_s).map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "n_s = {:.12}", i.n_s).map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "tau_reio = {:.12}", i.tau_reio).map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "output = tCl,lCl,pCl").map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "lensing = yes").map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "l_max_scalars = {}", lmax).map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "format = camb").map_err(|e| format!("write ini: {e}"))?;
    writeln!(f, "root = {}", output_root).map_err(|e| format!("write ini: {e}"))?;
    Ok(())
}

pub fn run_class(class_bin: &str, ini_path: &Path) -> Result<(), String> {
    let status = Command::new(class_bin)
        .arg(ini_path)
        .status()
        .map_err(|e| format!("failed to run CLASS binary '{class_bin}': {e}"))?;
    if !status.success() {
        return Err(format!(
            "CLASS failed for {:?} with exit status {}",
            ini_path, status
        ));
    }
    Ok(())
}

pub fn run_classy_fallback(
    python_bin: &str,
    out_path: &Path,
    lmax: u32,
    i: ClassRunInputs,
) -> Result<(), String> {
    // Keep this script tiny and pure-stdlib so only `classy` is required.
    let script = r##"
import math
import pathlib
import sys

try:
    from classy import Class
except Exception as e:
    raise RuntimeError(f"classy import failed: {e}")

out_path = pathlib.Path(sys.argv[1])
lmax = int(sys.argv[2])
h = float(sys.argv[3])
omega_b = float(sys.argv[4])
omega_cdm = float(sys.argv[5])
omega_k = float(sys.argv[6])
omega_lambda = float(sys.argv[7])
n_s = float(sys.argv[8])
a_s = float(sys.argv[9])
tau_reio = float(sys.argv[10])

params = {
    "h": h,
    "omega_b": omega_b,
    "omega_cdm": omega_cdm,
    "Omega_k": omega_k,
    "Omega_Lambda": omega_lambda,
    "n_s": n_s,
    "A_s": a_s,
    "tau_reio": tau_reio,
    "output": "tCl,lCl,pCl",
    "lensing": "yes",
    "l_max_scalars": lmax,
}

cosmo = Class()
cosmo.set(params)
cosmo.compute()

try:
    cls = cosmo.lensed_cl(lmax, CMB_unit="muK")
except TypeError:
    cls = cosmo.lensed_cl(lmax)

ells = cls["ell"]
tt = cls["tt"]

out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w", encoding="utf-8") as f:
    f.write("# ell D_ell_tt_uk2 (classy fallback)\n")
    for ell, c_tt in zip(ells, tt):
        ell_i = int(round(float(ell)))
        if ell_i < 2:
            continue
        # Convert C_ell -> D_ell = ell(ell+1)C_ell/(2pi)
        d_ell = (ell_i * (ell_i + 1) / (2.0 * math.pi)) * float(c_tt)
        f.write(f"{ell_i} {d_ell:.16e}\n")

cosmo.struct_cleanup()
cosmo.empty()
"##;

    let status = Command::new(python_bin)
        .arg("-c")
        .arg(script)
        .arg(out_path)
        .arg(lmax.to_string())
        .arg(i.h.to_string())
        .arg(i.omega_b.to_string())
        .arg(i.omega_cdm.to_string())
        .arg(i.omega_k.to_string())
        .arg(i.omega_lambda.to_string())
        .arg(i.n_s.to_string())
        .arg(i.a_s.to_string())
        .arg(i.tau_reio.to_string())
        .status()
        .map_err(|e| format!("failed to run python fallback '{python_bin}': {e}"))?;
    if !status.success() {
        return Err(format!(
            "classy fallback failed with exit status {} (python='{}')",
            status, python_bin
        ));
    }
    Ok(())
}

pub fn read_planck_tt_csv(path: &Path) -> Result<Vec<PlanckTtPoint>, String> {
    let f = File::open(path).map_err(|e| format!("open Planck CSV {:?}: {e}", path))?;
    let mut out = Vec::new();
    for (idx, line) in BufReader::new(f).lines().enumerate() {
        let line = line.map_err(|e| format!("read {:?} line {}: {e}", path, idx + 1))?;
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') {
            continue;
        }

        // Support both:
        //   CSV: ell,d_ell_tt_uk2,sigma_uk2
        // and Planck R3 whitespace text:
        //   l Dl -dDl +dDl BestFit
        let parsed = if s.contains(',') {
            let fields: Vec<&str> = s.split(',').map(|x| x.trim()).collect();
            if fields.len() < 3 {
                if idx == 0 {
                    continue;
                }
                return Err(format!(
                    "expected at least 3 CSV columns in {:?} line {}",
                    path,
                    idx + 1
                ));
            }
            let ell: u32 = match fields[0].parse() {
                Ok(v) => v,
                Err(_) if idx == 0 => continue,
                Err(e) => {
                    return Err(format!(
                        "parse ell in {:?} line {} ('{}'): {e}",
                        path,
                        idx + 1,
                        fields[0]
                    ));
                }
            };
            let d_ell_tt_uk2: f64 = fields[1].parse().map_err(|e| {
                format!(
                    "parse D_ell in {:?} line {} ('{}'): {e}",
                    path,
                    idx + 1,
                    fields[1]
                )
            })?;
            let sigma_uk2: f64 = fields[2].parse().map_err(|e| {
                format!(
                    "parse sigma in {:?} line {} ('{}'): {e}",
                    path,
                    idx + 1,
                    fields[2]
                )
            })?;
            (ell, d_ell_tt_uk2, sigma_uk2)
        } else {
            let fields: Vec<&str> = s.split_whitespace().collect();
            if fields.len() < 4 {
                if idx == 0 {
                    continue;
                }
                return Err(format!(
                    "expected >=4 whitespace columns in {:?} line {}",
                    path,
                    idx + 1
                ));
            }
            let ell_f: f64 = fields[0].parse().map_err(|e| {
                format!(
                    "parse l in {:?} line {} ('{}'): {e}",
                    path,
                    idx + 1,
                    fields[0]
                )
            })?;
            let ell = ell_f.round() as u32;
            let d_ell_tt_uk2: f64 = fields[1].parse().map_err(|e| {
                format!(
                    "parse Dl in {:?} line {} ('{}'): {e}",
                    path,
                    idx + 1,
                    fields[1]
                )
            })?;
            let sigma_minus: f64 = fields[2].parse().map_err(|e| {
                format!(
                    "parse -dDl in {:?} line {} ('{}'): {e}",
                    path,
                    idx + 1,
                    fields[2]
                )
            })?;
            let sigma_plus: f64 = fields[3].parse().map_err(|e| {
                format!(
                    "parse +dDl in {:?} line {} ('{}'): {e}",
                    path,
                    idx + 1,
                    fields[3]
                )
            })?;
            let sigma_uk2 = 0.5 * (sigma_minus.abs() + sigma_plus.abs());
            (ell, d_ell_tt_uk2, sigma_uk2)
        };

        let (ell, d_ell_tt_uk2, sigma_uk2) = parsed;
        if sigma_uk2 <= 0.0 {
            return Err(format!(
                "non-positive sigma in {:?} line {}: {}",
                path,
                idx + 1,
                sigma_uk2
            ));
        }
        out.push(PlanckTtPoint {
            ell,
            d_ell_tt_uk2,
            sigma_uk2,
        });
    }
    if out.is_empty() {
        return Err(format!("no usable points found in Planck CSV {:?}", path));
    }
    out.sort_by_key(|p| p.ell);
    Ok(out)
}

/// Generic Planck `D_ell` parser; supports both CSV and Planck whitespace files.
/// Kept separate from `read_planck_tt_csv` so TE/EE channels can use the same code path.
pub fn read_planck_dl_csv(path: &Path) -> Result<Vec<PlanckTtPoint>, String> {
    read_planck_tt_csv(path)
}

pub fn read_class_tt_camb(
    path: &Path,
    ell_min: u32,
    ell_max: u32,
) -> Result<Vec<ClassTtPoint>, String> {
    read_class_dl_camb_column(path, ell_min, ell_max, 2)
}

/// Read a CLASS CAMB-format CMB spectrum column.
///
/// The `column_1_based` index follows the CLASS header convention:
/// - 1: ell
/// - 2: TT
/// - 3: EE
/// - 4: BB
/// - 5: TE
pub fn read_class_dl_camb_column(
    path: &Path,
    ell_min: u32,
    ell_max: u32,
    column_1_based: usize,
) -> Result<Vec<ClassTtPoint>, String> {
    if column_1_based < 2 {
        return Err(format!(
            "invalid CLASS column index {} (must be >=2: 2=TT,3=EE,5=TE)",
            column_1_based
        ));
    }
    let f = File::open(path).map_err(|e| format!("open CLASS output {:?}: {e}", path))?;
    let mut out = Vec::new();
    for (idx, line) in BufReader::new(f).lines().enumerate() {
        let line = line.map_err(|e| format!("read {:?} line {}: {e}", path, idx + 1))?;
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = s.split_whitespace().collect();
        if fields.len() < column_1_based {
            continue;
        }
        let ell: u32 = match fields[0].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if ell < ell_min || ell > ell_max {
            continue;
        }
        let d_ell_tt_uk2: f64 = fields[column_1_based - 1].parse().map_err(|e| {
            format!(
                "parse CLASS D_ell column {} in {:?} line {} ('{}'): {e}",
                column_1_based,
                path,
                idx + 1,
                fields[column_1_based - 1]
            )
        })?;
        out.push(ClassTtPoint { ell, d_ell_tt_uk2 });
    }
    if out.is_empty() {
        return Err(format!(
            "no TT points in {:?} for ell range [{ell_min}, {ell_max}]",
            path
        ));
    }
    out.sort_by_key(|p| p.ell);
    Ok(out)
}

fn interp_class_tt(points: &[ClassTtPoint], ell: u32) -> Option<f64> {
    if points.is_empty() {
        return None;
    }
    let x = ell as f64;
    let first = points[0];
    let last = points[points.len() - 1];
    if x < first.ell as f64 || x > last.ell as f64 {
        return None;
    }
    for w in points.windows(2) {
        let a = w[0];
        let b = w[1];
        let xa = a.ell as f64;
        let xb = b.ell as f64;
        if x >= xa && x <= xb {
            let t = if xb > xa { (x - xa) / (xb - xa) } else { 0.0 };
            return Some(a.d_ell_tt_uk2 * (1.0 - t) + b.d_ell_tt_uk2 * t);
        }
    }
    Some(last.d_ell_tt_uk2)
}

pub fn compare_class_to_planck(
    class_tt: &[ClassTtPoint],
    planck_tt: &[PlanckTtPoint],
) -> Result<ClassPlanckFit, String> {
    let mut residuals = Vec::new();
    for p in planck_tt {
        let pred = interp_class_tt(class_tt, p.ell).ok_or_else(|| {
            format!(
                "CLASS spectrum does not bracket Planck multipole ell={}",
                p.ell
            )
        })?;
        let pull = (pred - p.d_ell_tt_uk2) / p.sigma_uk2;
        residuals.push(TtResidual {
            ell: p.ell,
            observed_uk2: p.d_ell_tt_uk2,
            predicted_uk2: pred,
            sigma_uk2: p.sigma_uk2,
            pull,
        });
    }
    if residuals.is_empty() {
        return Err("no overlap between CLASS and Planck points".to_string());
    }
    let n = residuals.len();
    let chi2 = residuals.iter().map(|r| r.pull * r.pull).sum::<f64>();
    let mean_abs_pull = residuals.iter().map(|r| r.pull.abs()).sum::<f64>() / n as f64;
    let max_abs_pull = residuals
        .iter()
        .map(|r| r.pull.abs())
        .fold(0.0_f64, f64::max);
    let rms_residual_uk2 = (residuals
        .iter()
        .map(|r| {
            let d = r.predicted_uk2 - r.observed_uk2;
            d * d
        })
        .sum::<f64>()
        / n as f64)
        .sqrt();
    let ndof = (n as i64 - 1).max(1) as f64;
    Ok(ClassPlanckFit {
        n_points: n,
        chi2,
        reduced_chi2: chi2 / ndof,
        mean_abs_pull,
        max_abs_pull,
        rms_residual_uk2,
        residuals,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn interpolation_and_fit_are_stable() {
        let class_tt = vec![
            ClassTtPoint {
                ell: 2,
                d_ell_tt_uk2: 100.0,
            },
            ClassTtPoint {
                ell: 10,
                d_ell_tt_uk2: 200.0,
            },
            ClassTtPoint {
                ell: 20,
                d_ell_tt_uk2: 300.0,
            },
        ];
        let planck = vec![
            PlanckTtPoint {
                ell: 2,
                d_ell_tt_uk2: 100.0,
                sigma_uk2: 10.0,
            },
            PlanckTtPoint {
                ell: 6,
                d_ell_tt_uk2: 150.0,
                sigma_uk2: 10.0,
            },
            PlanckTtPoint {
                ell: 20,
                d_ell_tt_uk2: 300.0,
                sigma_uk2: 10.0,
            },
        ];
        let fit = compare_class_to_planck(&class_tt, &planck).expect("fit");
        assert_eq!(fit.n_points, 3);
        assert!(
            fit.chi2 < 1e-9,
            "expected exact interpolation fit, got {fit:?}"
        );
    }

    #[test]
    fn read_planck_parser_supports_whitespace_and_csv() {
        let dir = std::env::temp_dir();
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos();

        let txt_path = dir.join(format!("planck_tt_test_{stamp}.txt"));
        let csv_path = dir.join(format!("planck_tt_test_{stamp}.csv"));

        fs::write(
            &txt_path,
            "# l Dl -dDl +dDl BestFit\n1.00e+02 3.00e+03 5.0e+01 5.0e+01 2.95e+03\n",
        )
        .expect("write txt");
        fs::write(&csv_path, "ell,d_ell_tt_uk2,sigma_uk2\n100,3000,50\n").expect("write csv");

        let txt = read_planck_tt_csv(&txt_path).expect("parse txt");
        let csv = read_planck_tt_csv(&csv_path).expect("parse csv");
        assert_eq!(txt.len(), 1);
        assert_eq!(csv.len(), 1);
        assert_eq!(txt[0].ell, 100);
        assert_eq!(csv[0].ell, 100);
        assert!((txt[0].d_ell_tt_uk2 - 3000.0).abs() < 1e-9);
        assert!((csv[0].d_ell_tt_uk2 - 3000.0).abs() < 1e-9);
        assert!((txt[0].sigma_uk2 - 50.0).abs() < 1e-9);
        assert!((csv[0].sigma_uk2 - 50.0).abs() < 1e-9);

        let _ = fs::remove_file(&txt_path);
        let _ = fs::remove_file(&csv_path);
    }
}
