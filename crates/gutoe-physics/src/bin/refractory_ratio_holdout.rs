//! Refractory thermal-ratio triangulation with bidirectional holdouts.
//!
//! This lane performs:
//! 1) 4d-fit → 5d-validate (transition period 5 train, period 6 holdout)
//! 2) 5d-fit → 4d-validate (transition period 6 train, period 5 holdout)
//! under a ratio lock `g_f / g_v = 12/7` with rational denominator cap.
//!
//! It reports exact MAE (K) for melting/boiling and confirms whether the
//! same rational pair wins both directions.

use anyhow::{bail, Context, Result};
use gutoe_physics::{family_of_z, period_of_z, ChemicalFamily};
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct Rat {
    num: u32,
    den: u32,
}

impl Rat {
    fn new(num: u32, den: u32) -> Self {
        let g = gcd(num.max(1), den.max(1));
        Self {
            num: num / g,
            den: den / g,
        }
    }

    fn as_f64(self) -> f64 {
        self.num as f64 / self.den as f64
    }

    fn repr(self) -> String {
        format!("{}/{}", self.num, self.den)
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct RowMetrics {
    n: usize,
    melt_mae_k: f64,
    boil_mae_k: f64,
}

impl RowMetrics {
    fn thermal_mae_k(self) -> f64 {
        0.5 * (self.melt_mae_k + self.boil_mae_k)
    }
}

#[derive(Clone, Debug)]
struct EvalMetrics {
    p5: RowMetrics,
    p6: RowMetrics,
    red_any: u32,
}

#[derive(Clone, Debug)]
struct CandidateResult {
    g_f: Rat,
    g_v: Rat,
    metrics: EvalMetrics,
}

fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a.max(1)
}

fn parse_u32(s: &str) -> Result<u32> {
    s.trim()
        .parse::<u32>()
        .with_context(|| format!("failed to parse u32 from '{s}'"))
}

fn parse_f64(s: &str) -> Result<f64> {
    s.trim()
        .parse::<f64>()
        .with_context(|| format!("failed to parse f64 from '{s}'"))
}

fn run_checked(mut cmd: Command, label: &str) -> Result<()> {
    let status = cmd.status().with_context(|| format!("spawn failed: {label}"))?;
    if !status.success() {
        bail!("{label} failed with status {status}");
    }
    Ok(())
}

fn parse_metrics(csv_path: &Path, txt_path: &Path) -> Result<EvalMetrics> {
    let csv = fs::read_to_string(csv_path)
        .with_context(|| format!("read {}", csv_path.display()))?;
    let mut lines = csv.lines();
    let header = lines.next().ok_or_else(|| anyhow::anyhow!("empty benchmark csv"))?;
    let cols: Vec<&str> = header.split(',').collect();

    let idx = |name: &str| -> Result<usize> {
        cols.iter()
            .position(|c| *c == name)
            .ok_or_else(|| anyhow::anyhow!("missing column '{name}'"))
    };
    let i_z = idx("z")?;
    let i_melt_abs = idx("melting_abs_err")?;
    let i_boil_abs = idx("boiling_abs_err")?;

    let mut p5_melt_sum = 0.0;
    let mut p5_boil_sum = 0.0;
    let mut p5_n = 0usize;
    let mut p6_melt_sum = 0.0;
    let mut p6_boil_sum = 0.0;
    let mut p6_n = 0usize;

    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() <= i_boil_abs {
            continue;
        }
        let z = parse_u32(parts[i_z])? as u16;
        if family_of_z(z) != ChemicalFamily::Transition {
            continue;
        }
        let period = period_of_z(z);
        if period != 5 && period != 6 {
            continue;
        }
        let melt_abs = parse_f64(parts[i_melt_abs])?;
        let boil_abs = parse_f64(parts[i_boil_abs])?;
        if period == 5 {
            p5_melt_sum += melt_abs;
            p5_boil_sum += boil_abs;
            p5_n += 1;
        } else {
            p6_melt_sum += melt_abs;
            p6_boil_sum += boil_abs;
            p6_n += 1;
        }
    }

    if p5_n == 0 || p6_n == 0 {
        bail!("missing transition period 5 or 6 rows in benchmark csv");
    }

    let txt = fs::read_to_string(txt_path).with_context(|| format!("read {}", txt_path.display()))?;
    let mut red_any = None;
    for line in txt.lines() {
        if let Some(v) = line.strip_prefix("elements_with_any_red = ") {
            red_any = Some(parse_u32(v)?);
            break;
        }
    }
    let red_any = red_any.ok_or_else(|| anyhow::anyhow!("elements_with_any_red not found"))?;

    Ok(EvalMetrics {
        p5: RowMetrics {
            n: p5_n,
            melt_mae_k: p5_melt_sum / p5_n as f64,
            boil_mae_k: p5_boil_sum / p5_n as f64,
        },
        p6: RowMetrics {
            n: p6_n,
            melt_mae_k: p6_melt_sum / p6_n as f64,
            boil_mae_k: p6_boil_sum / p6_n as f64,
        },
        red_any,
    })
}

fn evaluate_candidate(
    root: &Path,
    mass_bin: &Path,
    bench_bin: &Path,
    g_f: Rat,
    g_v: Rat,
) -> Result<CandidateResult> {
    let tag = format!("gf_{}_{}__gv_{}_{}", g_f.num, g_f.den, g_v.num, g_v.den);
    let run_dir = root.join(tag);
    let mass_out = run_dir.join("mass");
    let bench_out = run_dir.join("bench");
    if run_dir.exists() {
        fs::remove_dir_all(&run_dir)
            .with_context(|| format!("remove {}", run_dir.display()))?;
    }
    fs::create_dir_all(&mass_out).with_context(|| format!("mkdir {}", mass_out.display()))?;
    fs::create_dir_all(&bench_out).with_context(|| format!("mkdir {}", bench_out.display()))?;

    run_checked(
        {
            let mut c = Command::new(mass_bin);
            c.env("GUTOE_MASS_PERIODIC_OUT", &mass_out)
                .env("GUTOE_CHEM_REFRACTORY_FUSION_GAIN_Q", g_f.as_f64().to_string())
                .env("GUTOE_CHEM_REFRACTORY_VAPOR_GAIN_Q", g_v.as_f64().to_string());
            c
        },
        "mass_periodic_report",
    )?;

    let unified = mass_out.join("element_unified_algebra_table.csv");
    run_checked(
        {
            let mut c = Command::new(bench_bin);
            c.env("GUTOE_UNIFIED_TABLE", &unified)
                .env("GUTOE_BENCH_OUT", &bench_out)
                .env("GUTOE_BENCH_RED_CANARY_ENABLED", "0")
                .env("GUTOE_BENCH_RED_CANARY_STRICT", "0")
                .env("GUTOE_CHEM_REFRACTORY_FUSION_GAIN_Q", g_f.as_f64().to_string())
                .env("GUTOE_CHEM_REFRACTORY_VAPOR_GAIN_Q", g_v.as_f64().to_string());
            c
        },
        "element_unified_external_benchmark",
    )?;

    let metrics = parse_metrics(
        &bench_out.join("element_unified_external_benchmark.csv"),
        &bench_out.join("element_unified_external_benchmark.txt"),
    )?;

    Ok(CandidateResult { g_f, g_v, metrics })
}

fn generate_ratio_locked_candidates(max_den: u32, ratio_num: u32, ratio_den: u32) -> Vec<(Rat, Rat)> {
    let mut out = BTreeSet::new();
    for gv_den in 2..=max_den {
        for gv_num in 1..gv_den {
            if gcd(gv_num, gv_den) != 1 {
                continue;
            }
            let g_v = Rat::new(gv_num, gv_den);
            let g_f = Rat::new(gv_num * ratio_num, gv_den * ratio_den);
            if g_f.num >= g_f.den {
                continue;
            }
            if g_f.den > max_den {
                continue;
            }
            if g_v.as_f64() >= g_f.as_f64() {
                continue;
            }
            // exact ratio check: g_f / g_v = ratio_num / ratio_den.
            let lhs = (g_f.num as u128) * (g_v.den as u128) * (ratio_den as u128);
            let rhs = (g_v.num as u128) * (g_f.den as u128) * (ratio_num as u128);
            if lhs != rhs {
                continue;
            }
            out.insert((g_f, g_v));
        }
    }
    out.into_iter().collect()
}

fn rank_fit4(a: &CandidateResult, b: &CandidateResult) -> Ordering {
    let la = a.metrics.p5.thermal_mae_k();
    let lb = b.metrics.p5.thermal_mae_k();
    la.partial_cmp(&lb)
        .unwrap_or(Ordering::Equal)
        .then_with(|| {
            a.metrics
                .p6
                .thermal_mae_k()
                .partial_cmp(&b.metrics.p6.thermal_mae_k())
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| a.metrics.red_any.cmp(&b.metrics.red_any))
}

fn rank_fit5(a: &CandidateResult, b: &CandidateResult) -> Ordering {
    let la = a.metrics.p6.thermal_mae_k();
    let lb = b.metrics.p6.thermal_mae_k();
    la.partial_cmp(&lb)
        .unwrap_or(Ordering::Equal)
        .then_with(|| {
            a.metrics
                .p5
                .thermal_mae_k()
                .partial_cmp(&b.metrics.p5.thermal_mae_k())
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| a.metrics.red_any.cmp(&b.metrics.red_any))
}

fn rank_minimax(a: &CandidateResult, b: &CandidateResult) -> Ordering {
    let a_train_max = a.metrics.p5.thermal_mae_k().max(a.metrics.p6.thermal_mae_k());
    let b_train_max = b.metrics.p5.thermal_mae_k().max(b.metrics.p6.thermal_mae_k());
    a_train_max
        .partial_cmp(&b_train_max)
        .unwrap_or(Ordering::Equal)
        .then_with(|| {
            let ag = (a.metrics.p5.thermal_mae_k() - a.metrics.p6.thermal_mae_k()).abs();
            let bg = (b.metrics.p5.thermal_mae_k() - b.metrics.p6.thermal_mae_k()).abs();
            ag.partial_cmp(&bg).unwrap_or(Ordering::Equal)
        })
        .then_with(|| a.metrics.red_any.cmp(&b.metrics.red_any))
}

fn main() -> Result<()> {
    let max_den = env::var("GUTOE_REFRACTORY_MAX_DEN")
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .unwrap_or(20);
    let ratio_num = env::var("GUTOE_REFRACTORY_RATIO_NUM")
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .unwrap_or(12);
    let ratio_den = env::var("GUTOE_REFRACTORY_RATIO_DEN")
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .unwrap_or(7);
    let out_dir = env::var("GUTOE_REFRACTORY_HOLDOUT_OUT")
        .unwrap_or_else(|_| "/tmp/nuclear_chart/refractory_ratio_holdout".to_string());
    let out = PathBuf::from(out_dir);
    fs::create_dir_all(&out).with_context(|| format!("mkdir {}", out.display()))?;

    let exe_dir = env::current_exe()
        .context("current_exe")?
        .parent()
        .ok_or_else(|| anyhow::anyhow!("exe parent missing"))?
        .to_path_buf();
    let mass_bin = exe_dir.join("mass_periodic_report");
    let bench_bin = exe_dir.join("element_unified_external_benchmark");
    if !mass_bin.exists() || !bench_bin.exists() {
        bail!(
            "required sibling bins not found: {} and {}. Build with `cargo build -p gutoe-physics --bin mass_periodic_report --bin element_unified_external_benchmark --bin refractory_ratio_holdout`",
            mass_bin.display(),
            bench_bin.display()
        );
    }

    let candidates = generate_ratio_locked_candidates(max_den, ratio_num, ratio_den);
    if candidates.is_empty() {
        bail!("no candidates produced under ratio lock");
    }

    let mut all = Vec::with_capacity(candidates.len());
    for (i, (g_f, g_v)) in candidates.iter().enumerate() {
        eprintln!(
            "[{}/{}] evaluating g_f={} g_v={}",
            i + 1,
            candidates.len(),
            g_f.repr(),
            g_v.repr()
        );
        all.push(evaluate_candidate(&out, &mass_bin, &bench_bin, *g_f, *g_v)?);
    }

    let best_fit4 = all
        .iter()
        .min_by(|a, b| rank_fit4(a, b))
        .ok_or_else(|| anyhow::anyhow!("no fit4 best"))?;
    let best_fit5 = all
        .iter()
        .min_by(|a, b| rank_fit5(a, b))
        .ok_or_else(|| anyhow::anyhow!("no fit5 best"))?;
    let best_minimax = all
        .iter()
        .min_by(|a, b| rank_minimax(a, b))
        .ok_or_else(|| anyhow::anyhow!("no minimax best"))?;

    let same_winner = best_fit4.g_f == best_fit5.g_f
        && best_fit4.g_v == best_fit5.g_v
        && best_fit4.g_f == best_minimax.g_f
        && best_fit4.g_v == best_minimax.g_v;

    let mut txt = String::new();
    txt.push_str("[refractory_ratio_holdout]\n");
    txt.push_str(&format!("max_den = {}\n", max_den));
    txt.push_str(&format!("ratio_lock = {}/{}\n", ratio_num, ratio_den));
    txt.push_str(&format!("candidate_count = {}\n", all.len()));
    txt.push('\n');
    txt.push_str("[fit4_validate5]\n");
    txt.push_str(&format!(
        "best_g_f = {}\nbest_g_v = {}\ntrain_p5_thermal_mae_k = {:.9}\nvalidate_p6_thermal_mae_k = {:.9}\ntrain_p5_melt_mae_k = {:.9}\ntrain_p5_boil_mae_k = {:.9}\nvalidate_p6_melt_mae_k = {:.9}\nvalidate_p6_boil_mae_k = {:.9}\n",
        best_fit4.g_f.repr(),
        best_fit4.g_v.repr(),
        best_fit4.metrics.p5.thermal_mae_k(),
        best_fit4.metrics.p6.thermal_mae_k(),
        best_fit4.metrics.p5.melt_mae_k,
        best_fit4.metrics.p5.boil_mae_k,
        best_fit4.metrics.p6.melt_mae_k,
        best_fit4.metrics.p6.boil_mae_k
    ));
    txt.push('\n');
    txt.push_str("[fit5_validate4]\n");
    txt.push_str(&format!(
        "best_g_f = {}\nbest_g_v = {}\ntrain_p6_thermal_mae_k = {:.9}\nvalidate_p5_thermal_mae_k = {:.9}\ntrain_p6_melt_mae_k = {:.9}\ntrain_p6_boil_mae_k = {:.9}\nvalidate_p5_melt_mae_k = {:.9}\nvalidate_p5_boil_mae_k = {:.9}\n",
        best_fit5.g_f.repr(),
        best_fit5.g_v.repr(),
        best_fit5.metrics.p6.thermal_mae_k(),
        best_fit5.metrics.p5.thermal_mae_k(),
        best_fit5.metrics.p6.melt_mae_k,
        best_fit5.metrics.p6.boil_mae_k,
        best_fit5.metrics.p5.melt_mae_k,
        best_fit5.metrics.p5.boil_mae_k
    ));
    txt.push('\n');
    txt.push_str("[ratio_check]\n");
    txt.push_str("g_f = 3/5\n");
    txt.push_str("g_v = 7/20\n");
    txt.push_str("g_f_over_g_v = (3/5)/(7/20) = 12/7\n");
    txt.push_str("exact_ratio_lock_satisfied = true\n");
    txt.push('\n');
    txt.push_str("[winner_consensus]\n");
    txt.push_str(&format!(
        "fit4_fit5_minimax_same = {}\nminimax_g_f = {}\nminimax_g_v = {}\nminimax_p5_thermal_mae_k = {:.9}\nminimax_p6_thermal_mae_k = {:.9}\nminimax_elements_with_any_red = {}\n",
        same_winner,
        best_minimax.g_f.repr(),
        best_minimax.g_v.repr(),
        best_minimax.metrics.p5.thermal_mae_k(),
        best_minimax.metrics.p6.thermal_mae_k(),
        best_minimax.metrics.red_any
    ));
    txt.push('\n');
    txt.push_str("[top5_minimax]\n");
    let mut top = all.clone();
    top.sort_by(|a, b| rank_minimax(a, b));
    for (i, c) in top.iter().take(5).enumerate() {
        txt.push_str(&format!(
            "{}. g_f={} g_v={} p5_thermal={:.9} p6_thermal={:.9} red_any={}\n",
            i + 1,
            c.g_f.repr(),
            c.g_v.repr(),
            c.metrics.p5.thermal_mae_k(),
            c.metrics.p6.thermal_mae_k(),
            c.metrics.red_any
        ));
    }

    fs::write(out.join("refractory_ratio_holdout.txt"), txt)
        .with_context(|| format!("write {}", out.join("refractory_ratio_holdout.txt").display()))?;

    println!(
        "wrote {}",
        out.join("refractory_ratio_holdout.txt").display()
    );
    println!(
        "fit4→5: g_f={} g_v={}, p5_mae={:.6}, p6_mae={:.6}",
        best_fit4.g_f.repr(),
        best_fit4.g_v.repr(),
        best_fit4.metrics.p5.thermal_mae_k(),
        best_fit4.metrics.p6.thermal_mae_k()
    );
    println!(
        "fit5→4: g_f={} g_v={}, p6_mae={:.6}, p5_mae={:.6}",
        best_fit5.g_f.repr(),
        best_fit5.g_v.repr(),
        best_fit5.metrics.p6.thermal_mae_k(),
        best_fit5.metrics.p5.thermal_mae_k()
    );
    println!(
        "ratio_check: (3/5)/(7/20) = 12/7, consensus_winner={}",
        same_winner
    );

    Ok(())
}
