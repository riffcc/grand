//! Propagate quark-mass uncertainties into s²_down.
//!
//! Computes
//!   s²_down = 6K(d,s,b) - 2
//! where
//!   K = Σm / (Σ√m)^2.
//!
//! Defaults use broad PDG-style light-quark uncertainties (as discussed in-lane),
//! and can be overridden via env vars:
//!   GUTOE_MD_MEV, GUTOE_MS_MEV, GUTOE_MB_MEV
//!   GUTOE_SIG_MD_MEV, GUTOE_SIG_MS_MEV, GUTOE_SIG_MB_MEV
//!   GUTOE_MC_SAMPLES

use serde::Serialize;
use std::env;
use std::fs;

#[derive(Debug, Clone, Serialize)]
struct Config {
    md_mev: f64,
    ms_mev: f64,
    mb_mev: f64,
    sig_md_mev: f64,
    sig_ms_mev: f64,
    sig_mb_mev: f64,
    mc_samples: usize,
}

#[derive(Debug, Clone, Serialize)]
struct MonteCarloSummary {
    mean: f64,
    sd: f64,
    q16: f64,
    q84: f64,
    q025: f64,
    q975: f64,
    z_target: f64,
    p_s2_ge_target: f64,
    target_in_1sigma: bool,
    target_in_95: bool,
}

#[derive(Debug, Clone, Serialize)]
struct BoxSummary {
    min_s2: f64,
    max_s2: f64,
    target_in_box: bool,
}

#[derive(Debug, Clone, Serialize)]
struct Report {
    config: Config,
    s2_center: f64,
    s2_target: f64,
    delta_center_minus_target: f64,
    rel_percent_center_minus_target: f64,
    sensitivity_ds2_dmd: f64,
    sensitivity_ds2_dms: f64,
    sensitivity_ds2_dmb: f64,
    monte_carlo: MonteCarloSummary,
    box_scan: BoxSummary,
    summary: String,
}

fn env_f64(name: &str, default: f64) -> f64 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_usize(name: &str, default: usize) -> usize {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
}

fn koide3(m1: f64, m2: f64, m3: f64) -> f64 {
    let sum = m1 + m2 + m3;
    let sum_sqrt = m1.sqrt() + m2.sqrt() + m3.sqrt();
    sum / (sum_sqrt * sum_sqrt)
}

fn s2_down(md: f64, ms: f64, mb: f64) -> f64 {
    6.0 * koide3(md, ms, mb) - 2.0
}

fn normal_sample(seed: &mut u64) -> f64 {
    // Box-Muller using deterministic LCG source (no external deps).
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    let u1 = (((*seed >> 11) as f64) + 1.0) / ((1u64 << 53) as f64 + 1.0);
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    let u2 = (((*seed >> 11) as f64) + 1.0) / ((1u64 << 53) as f64 + 1.0);
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    let n = sorted.len();
    let idx = ((q.clamp(0.0, 1.0)) * ((n - 1) as f64)).round() as usize;
    sorted[idx]
}

fn main() {
    let out_dir =
        env::var("GUTOE_DOWN_UNCERT_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);

    // Broad PDG-style defaults discussed in lane.
    let cfg = Config {
        md_mev: env_f64("GUTOE_MD_MEV", 4.67),
        ms_mev: env_f64("GUTOE_MS_MEV", 93.0),
        mb_mev: env_f64("GUTOE_MB_MEV", 4180.0),
        sig_md_mev: env_f64("GUTOE_SIG_MD_MEV", 0.48),
        sig_ms_mev: env_f64("GUTOE_SIG_MS_MEV", 11.0),
        sig_mb_mev: env_f64("GUTOE_SIG_MB_MEV", 30.0),
        mc_samples: env_usize("GUTOE_MC_SAMPLES", 600_000),
    };

    let s2_target = 2.0 + 4.0 / 9.0;
    let s2_center = s2_down(cfg.md_mev, cfg.ms_mev, cfg.mb_mev);
    let delta = s2_center - s2_target;
    let rel_pct = 100.0 * delta / s2_target;

    // Local sensitivities at center.
    let h = 1.0e-5;
    let dmd = cfg.md_mev * h;
    let dms = cfg.ms_mev * h;
    let dmb = cfg.mb_mev * h;
    let ds2_dmd = (s2_down(cfg.md_mev + dmd, cfg.ms_mev, cfg.mb_mev)
        - s2_down(cfg.md_mev - dmd, cfg.ms_mev, cfg.mb_mev))
        / (2.0 * dmd);
    let ds2_dms = (s2_down(cfg.md_mev, cfg.ms_mev + dms, cfg.mb_mev)
        - s2_down(cfg.md_mev, cfg.ms_mev - dms, cfg.mb_mev))
        / (2.0 * dms);
    let ds2_dmb = (s2_down(cfg.md_mev, cfg.ms_mev, cfg.mb_mev + dmb)
        - s2_down(cfg.md_mev, cfg.ms_mev, cfg.mb_mev - dmb))
        / (2.0 * dmb);

    // Gaussian MC.
    let mut seed = 0xC0FFEEu64;
    let mut vals = Vec::with_capacity(cfg.mc_samples);
    for _ in 0..cfg.mc_samples {
        let md = (cfg.md_mev + cfg.sig_md_mev * normal_sample(&mut seed)).max(1.0e-12);
        let ms = (cfg.ms_mev + cfg.sig_ms_mev * normal_sample(&mut seed)).max(1.0e-12);
        let mb = (cfg.mb_mev + cfg.sig_mb_mev * normal_sample(&mut seed)).max(1.0e-12);
        vals.push(s2_down(md, ms, mb));
    }
    vals.sort_by(|a, b| a.total_cmp(b));
    let mean = vals.iter().sum::<f64>() / vals.len() as f64;
    let sd = (vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / vals.len() as f64).sqrt();
    let q16 = quantile(&vals, 0.158_655);
    let q84 = quantile(&vals, 0.841_345);
    let q025 = quantile(&vals, 0.025);
    let q975 = quantile(&vals, 0.975);
    let z_target = if sd > 0.0 {
        (s2_target - mean) / sd
    } else {
        0.0
    };
    let n_ge = vals.iter().filter(|&&v| v >= s2_target).count();
    let p_s2_ge_target = n_ge as f64 / vals.len() as f64;
    let mc = MonteCarloSummary {
        mean,
        sd,
        q16,
        q84,
        q025,
        q975,
        z_target,
        p_s2_ge_target,
        target_in_1sigma: q16 <= s2_target && s2_target <= q84,
        target_in_95: q025 <= s2_target && s2_target <= q975,
    };

    // Uniform ±1σ box scan (corners + dense random).
    let mut min_s2 = f64::INFINITY;
    let mut max_s2 = f64::NEG_INFINITY;
    for &sd_sign in &[-1.0, 1.0] {
        for &ss_sign in &[-1.0, 1.0] {
            for &sb_sign in &[-1.0, 1.0] {
                let md = cfg.md_mev + sd_sign * cfg.sig_md_mev;
                let ms = cfg.ms_mev + ss_sign * cfg.sig_ms_mev;
                let mb = cfg.mb_mev + sb_sign * cfg.sig_mb_mev;
                let v = s2_down(md.max(1.0e-12), ms.max(1.0e-12), mb.max(1.0e-12));
                min_s2 = min_s2.min(v);
                max_s2 = max_s2.max(v);
            }
        }
    }
    for _ in 0..250_000usize {
        let u1 = 0.5 * (normal_sample(&mut seed).tanh() + 1.0);
        let u2 = 0.5 * (normal_sample(&mut seed).tanh() + 1.0);
        let u3 = 0.5 * (normal_sample(&mut seed).tanh() + 1.0);
        let md = cfg.md_mev + (2.0 * u1 - 1.0) * cfg.sig_md_mev;
        let ms = cfg.ms_mev + (2.0 * u2 - 1.0) * cfg.sig_ms_mev;
        let mb = cfg.mb_mev + (2.0 * u3 - 1.0) * cfg.sig_mb_mev;
        let v = s2_down(md.max(1.0e-12), ms.max(1.0e-12), mb.max(1.0e-12));
        min_s2 = min_s2.min(v);
        max_s2 = max_s2.max(v);
    }
    let box_scan = BoxSummary {
        min_s2,
        max_s2,
        target_in_box: min_s2 <= s2_target && s2_target <= max_s2,
    };

    let summary = format!(
        "down-uncertainty propagation: center={:.9}, target(2+4/9)={:.9}, z_target={:.3}, in_1sigma={}, in_95={}",
        s2_center, s2_target, mc.z_target, mc.target_in_1sigma, mc.target_in_95
    );

    let report = Report {
        config: cfg,
        s2_center,
        s2_target,
        delta_center_minus_target: delta,
        rel_percent_center_minus_target: rel_pct,
        sensitivity_ds2_dmd: ds2_dmd,
        sensitivity_ds2_dms: ds2_dms,
        sensitivity_ds2_dmb: ds2_dmb,
        monte_carlo: mc,
        box_scan,
        summary,
    };

    let txt_path = format!("{out_dir}/yukawa_down_uncertainty_propagation.txt");
    let json_path = format!("{out_dir}/yukawa_down_uncertainty_propagation.json");
    fs::write(
        &txt_path,
        format!(
            "[yukawa_down_uncertainty_propagation]\n\
             md={:.6}±{:.6} MeV\n\
             ms={:.6}±{:.6} MeV\n\
             mb={:.6}±{:.6} MeV\n\
             s2_center={:.12}\n\
             s2_target_2plus4over9={:.12}\n\
             delta_center_minus_target={:.12} ({:+.6}%)\n\
             ds2_dmd={:.12e} ds2_dms={:.12e} ds2_dmb={:.12e}\n\
             mc_mean={:.12}\n\
             mc_sd={:.12}\n\
             mc_1sigma=[{:.12}, {:.12}]\n\
             mc_95=[{:.12}, {:.12}]\n\
             mc_z_target={:.12}\n\
             mc_p_s2_ge_target={:.12}\n\
             mc_target_in_1sigma={}\n\
             mc_target_in_95={}\n\
             box_range=[{:.12}, {:.12}] target_in_box={}\n\
             summary={}\n",
            report.config.md_mev,
            report.config.sig_md_mev,
            report.config.ms_mev,
            report.config.sig_ms_mev,
            report.config.mb_mev,
            report.config.sig_mb_mev,
            report.s2_center,
            report.s2_target,
            report.delta_center_minus_target,
            report.rel_percent_center_minus_target,
            report.sensitivity_ds2_dmd,
            report.sensitivity_ds2_dms,
            report.sensitivity_ds2_dmb,
            report.monte_carlo.mean,
            report.monte_carlo.sd,
            report.monte_carlo.q16,
            report.monte_carlo.q84,
            report.monte_carlo.q025,
            report.monte_carlo.q975,
            report.monte_carlo.z_target,
            report.monte_carlo.p_s2_ge_target,
            report.monte_carlo.target_in_1sigma,
            report.monte_carlo.target_in_95,
            report.box_scan.min_s2,
            report.box_scan.max_s2,
            report.box_scan.target_in_box,
            report.summary
        ),
    )
    .expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&report).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!("{}", report.summary);
}
