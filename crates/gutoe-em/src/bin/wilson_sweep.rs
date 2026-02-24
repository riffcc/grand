use std::fs::File;
use std::io::Write;

use gutoe_em::{confinement_experiment, LatticeConfig};
use rand::rngs::StdRng;
use rand::SeedableRng;

fn mean(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        0.0
    } else {
        xs.iter().sum::<f64>() / xs.len() as f64
    }
}

fn stddev(xs: &[f64]) -> f64 {
    if xs.len() < 2 {
        0.0
    } else {
        let m = mean(xs);
        let var = xs.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / xs.len() as f64;
        var.sqrt()
    }
}

fn main() {
    let cfg = LatticeConfig {
        hex_rows: std::env::var("WILSON_ROWS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(24),
        hex_cols: std::env::var("WILSON_COLS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(24),
        layers: std::env::var("WILSON_LAYERS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(1),
        ..Default::default()
    };
    let seeds: usize = std::env::var("WILSON_SEEDS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    let n_therm: usize = std::env::var("WILSON_THERM")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(300);
    let n_meas: usize = std::env::var("WILSON_MEAS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(140);

    let betas = [0.2, 0.3, 0.5, 0.8, 1.0, 1.4, 2.0, 2.5, 3.0];
    let out = std::env::var("WILSON_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/confinement_wilson_sweep.csv".to_string());
    let mut f = File::create(&out).expect("create sweep csv");
    writeln!(
        f,
        "beta,seeds,therm,meas,eps,plaquette_mean,plaquette_std,v3_mean,v3_std"
    )
    .expect("csv header");

    println!(
        "Wilson confinement sweep: cfg={}x{}x{} seeds={} therm={} meas={}",
        cfg.hex_rows, cfg.hex_cols, cfg.layers, seeds, n_therm, n_meas
    );
    for beta in betas {
        let eps = if beta < 0.5 {
            0.9
        } else if beta < 1.0 {
            0.75
        } else if beta < 2.0 {
            0.45
        } else {
            0.25
        };
        let mut p = Vec::with_capacity(seeds);
        let mut v = Vec::with_capacity(seeds);
        for s in 0..seeds {
            let mut rng = StdRng::seed_from_u64(0xC0FFEE + s as u64 * 7919);
            let (plaq, v3) = confinement_experiment(beta, n_therm, n_meas, eps, &cfg, &mut rng);
            p.push(plaq);
            v.push(v3);
        }
        let pm = mean(&p);
        let ps = stddev(&p);
        let vm = mean(&v);
        let vs = stddev(&v);
        println!(
            "  beta={beta:.2} eps={eps:.2}  plaquette={pm:.5}±{ps:.5}  V3={vm:.5}±{vs:.5}"
        );
        writeln!(
            f,
            "{beta:.6},{seeds},{n_therm},{n_meas},{eps:.6},{pm:.9},{ps:.9},{vm:.9},{vs:.9}"
        )
        .expect("csv row");
    }
    println!("wrote {out}");
}
