//! RH Phase-3 all-lanes report:
//! determinant lock, de Branges proxy, Weil-positivity proxy,
//! prime-trace probe, inverse spectral reconstruction,
//! functional-equation symmetry embedding, and train/holdout/freeze checks.

use gutoe_physics::riemann_lane::{
    build_operator, fit_against_reference, map_eigenvalues, zeta_zero_reference, RiemannMapParams,
    RiemannOperatorParams,
};
use nalgebra::{DMatrix, SymmetricEigen};
use num_complex::Complex64;
use serde_json::json;
use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::Write;

const TRAIN_N: usize = 40;
const HOLD_N: usize = 40;
const FREEZE_N: usize = 40;

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
}

fn sorted_eigs(h: DMatrix<f64>) -> Vec<f64> {
    let eig = SymmetricEigen::new(h);
    let mut evals: Vec<f64> = eig.eigenvalues.iter().copied().collect();
    evals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    evals
}

fn nearest_dist_stats(target: &[f64], pred: &[f64]) -> (f64, f64) {
    let mut mean = 0.0;
    let mut maxd: f64 = 0.0;
    for &t in target {
        let mut best = f64::INFINITY;
        for &p in pred {
            let d = (t - p).abs();
            if d < best {
                best = d;
            }
        }
        mean += best;
        maxd = maxd.max(best);
    }
    mean /= target.len().max(1) as f64;
    (mean, maxd)
}

fn fit_affine(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len().min(y.len()).max(1) as f64;
    let mx = x.iter().take(n as usize).sum::<f64>() / n;
    let my = y.iter().take(n as usize).sum::<f64>() / n;
    let mut cov = 0.0;
    let mut var = 0.0;
    for i in 0..(n as usize) {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        cov += dx * dy;
        var += dx * dx;
    }
    let slope = if var > 1e-20 { cov / var } else { 1.0 };
    let shift = my - slope * mx;
    (slope, shift)
}

fn complex_even_product(z: Complex64, zeros: &[f64]) -> Complex64 {
    let mut a = Complex64::new(1.0, 0.0);
    for &g in zeros {
        let gg = g * g;
        a *= Complex64::new(1.0, 0.0) - (z * z) / Complex64::new(gg, 0.0);
    }
    a
}

fn complex_even_product_derivative(z: Complex64, zeros: &[f64]) -> Complex64 {
    let a = complex_even_product(z, zeros);
    let mut sum = Complex64::new(0.0, 0.0);
    for &g in zeros {
        let gg = g * g;
        let den = Complex64::new(gg, 0.0) - z * z;
        sum += (-2.0 * z) / den;
    }
    a * sum
}

fn debranges_hb_pass_fraction(zeros: &[f64]) -> f64 {
    let mut pass = 0usize;
    let mut total = 0usize;
    let ys = [0.3, 0.7, 1.3];
    let mut x = 10.0;
    while x <= 260.0 {
        for &y in &ys {
            let z = Complex64::new(x, y);
            let a = complex_even_product(z, zeros);
            let ap = complex_even_product_derivative(z, zeros);
            let e = a - Complex64::new(0.0, 1.0) * ap;
            let estar = a + Complex64::new(0.0, 1.0) * ap;
            if e.norm() > estar.norm() {
                pass += 1;
            }
            total += 1;
        }
        x += 5.0;
    }
    pass as f64 / total.max(1) as f64
}

fn min_separation(zeros: &[f64]) -> f64 {
    zeros
        .windows(2)
        .map(|w| w[1] - w[0])
        .fold(f64::INFINITY, |a, b| a.min(b))
}

fn weil_proxy_min_eig(zeros: &[f64], u: f64, beta: f64) -> f64 {
    let n = zeros.len();
    let mut m = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let d = zeros[i] - zeros[j];
            let k = (-u * d * d).exp() * (1.0 - beta * d * d);
            m[(i, j)] = k;
        }
    }
    let eig = SymmetricEigen::new(m);
    eig.eigenvalues
        .iter()
        .copied()
        .fold(f64::INFINITY, |a, b| a.min(b))
}

fn primes_upto_300() -> Vec<u64> {
    let mut p = Vec::new();
    'outer: for n in 2..=300u64 {
        let r = (n as f64).sqrt() as u64;
        for d in 2..=r {
            if n % d == 0 {
                continue 'outer;
            }
        }
        p.push(n);
    }
    p
}

fn trace_signal(zeros: &[f64], omega: f64) -> f64 {
    zeros.iter().map(|g| (omega * *g).cos()).sum::<f64>()
}

fn prime_trace_zscore_mean(zeros: &[f64], primes: &[u64]) -> f64 {
    let mut zsum = 0.0;
    let mut n = 0usize;
    for &p in primes {
        let w0 = (p as f64).ln();
        let peak = trace_signal(zeros, w0);
        let mut bg = Vec::new();
        let mut t: f64 = -0.30;
        while t <= 0.30 {
            if t.abs() > 0.03 {
                bg.push(trace_signal(zeros, w0 + t));
            }
            t += 0.03;
        }
        let m = bg.iter().sum::<f64>() / bg.len().max(1) as f64;
        let v = bg
            .iter()
            .map(|x| {
                let d = *x - m;
                d * d
            })
            .sum::<f64>()
            / bg.len().max(1) as f64;
        let s = v.sqrt().max(1e-9);
        zsum += (peak - m) / s;
        n += 1;
    }
    zsum / n.max(1) as f64
}

fn build_symmetry_embedded(op: RiemannOperatorParams) -> DMatrix<f64> {
    let a = build_operator(op);
    let n = a.nrows();
    let mut h = DMatrix::<f64>::zeros(2 * n, 2 * n);
    for i in 0..n {
        for j in 0..n {
            h[(i, n + j)] = a[(i, j)];
            h[(n + i, j)] = a[(j, i)];
        }
    }
    h
}

fn symmetry_residual(evals: &[f64]) -> f64 {
    let n = evals.len();
    let mut maxr: f64 = 0.0;
    for i in 0..(n / 2) {
        maxr = maxr.max((evals[i] + evals[n - 1 - i]).abs());
    }
    maxr
}

fn main() {
    let n = env_usize("GUTOE_RIEMANN_DIM", 512);
    let refs = zeta_zero_reference();
    let needed = TRAIN_N + HOLD_N + FREEZE_N;
    assert!(refs.len() >= needed, "need at least {needed} reference zeros");

    let train = &refs[..TRAIN_N];
    let hold = &refs[TRAIN_N..TRAIN_N + HOLD_N];
    let freeze = &refs[TRAIN_N + HOLD_N..needed];

    let op = RiemannOperatorParams {
        n,
        hop1_scale: 0.5,
        hop2_scale: 0.0,
        potential_scale: 0.0,
    };
    let map_struct = RiemannMapParams {
        slope: 11.0 / 18.0,
        shift: 13.0 * 24.0 + 8.0 / 17.0,
    };

    // Baseline spectrum and structural map.
    let h_base = build_operator(op);
    let sym_base = (&h_base - h_base.transpose()).norm();
    let raw_base = sorted_eigs(h_base);
    let raw120 = &raw_base[..needed.min(raw_base.len())];
    let pred_struct = map_eigenvalues(raw120, map_struct);

    let fit_train_struct = fit_against_reference(&pred_struct[..TRAIN_N], train);
    let fit_hold_struct = fit_against_reference(&pred_struct[TRAIN_N..TRAIN_N + HOLD_N], hold);
    let fit_freeze_struct = fit_against_reference(&pred_struct[TRAIN_N + HOLD_N..needed], freeze);

    // 1) Determinant-lock proxy via nearest root distance.
    let (det_train_mean, det_train_max) = nearest_dist_stats(train, &pred_struct);
    let (det_hold_mean, det_hold_max) = nearest_dist_stats(hold, &pred_struct);
    let (det_freeze_mean, det_freeze_max) = nearest_dist_stats(freeze, &pred_struct);

    // 2) de Branges proxy on predicted positive ordinates.
    let debranges_sep_min = min_separation(&pred_struct);
    let debranges_hb_fraction = debranges_hb_pass_fraction(&pred_struct);

    // 3) Weil-positivity proxy.
    let u = 0.0007;
    let beta = 0.00015;
    let min_eig_target = weil_proxy_min_eig(train, u, beta);
    let min_eig_pred = weil_proxy_min_eig(&pred_struct[..TRAIN_N], u, beta);

    // 4) Prime trace-formula probe.
    let primes = primes_upto_300();
    let prime_z_target = prime_trace_zscore_mean(train, &primes);
    let prime_z_pred = prime_trace_zscore_mean(&pred_struct[..TRAIN_N], &primes);

    // 5) Inverse spectral reconstruction (fit affine map on train only).
    let (fit_slope, fit_shift) = fit_affine(&raw120[..TRAIN_N], train);
    let map_fit = RiemannMapParams {
        slope: fit_slope,
        shift: fit_shift,
    };
    let pred_fit = map_eigenvalues(raw120, map_fit);
    let fit_hold_inv = fit_against_reference(&pred_fit[TRAIN_N..TRAIN_N + HOLD_N], hold);
    let fit_freeze_inv = fit_against_reference(&pred_fit[TRAIN_N + HOLD_N..needed], freeze);

    // 6) Functional-equation symmetry embedding.
    let h_sym = build_symmetry_embedded(op);
    let evals_sym = sorted_eigs(h_sym);
    let sym_embed_resid = symmetry_residual(&evals_sym);
    let mut pos_branch: Vec<f64> = evals_sym.into_iter().filter(|x| *x > 0.0).collect();
    pos_branch.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let pos120 = &pos_branch[..needed.min(pos_branch.len())];
    let (sym_slope, sym_shift) = fit_affine(&pos120[..TRAIN_N], train);
    let pred_sym = map_eigenvalues(
        pos120,
        RiemannMapParams {
            slope: sym_slope,
            shift: sym_shift,
        },
    );
    let fit_hold_sym = fit_against_reference(&pred_sym[TRAIN_N..TRAIN_N + HOLD_N], hold);
    let fit_freeze_sym = fit_against_reference(&pred_sym[TRAIN_N + HOLD_N..needed], freeze);

    // 7) Hilbert–Pólya formal bridge status (report lane only).
    let hilbert_polya_bridge_status = "exploratory_report_lane";

    // 8) Train/holdout/freeze summary is embedded across all lanes above.

    let report = json!({
        "lane": "riemann_phase3_all_lanes",
        "splits": {
            "train": [1, TRAIN_N],
            "holdout": [TRAIN_N + 1, TRAIN_N + HOLD_N],
            "freeze": [TRAIN_N + HOLD_N + 1, needed],
            "total_points": needed
        },
        "baseline": {
            "operator": {
                "n": n,
                "hop1": op.hop1_scale,
                "hop2": op.hop2_scale,
                "pot": op.potential_scale,
                "self_adjoint_residual": sym_base
            },
            "structural_map": {
                "slope": map_struct.slope,
                "shift": map_struct.shift
            },
            "fit_train_struct": fit_train_struct,
            "fit_hold_struct": fit_hold_struct,
            "fit_freeze_struct": fit_freeze_struct
        },
        "lane1_determinant_lock": {
            "train_mean_nearest_root_distance": det_train_mean,
            "train_max_nearest_root_distance": det_train_max,
            "hold_mean_nearest_root_distance": det_hold_mean,
            "hold_max_nearest_root_distance": det_hold_max,
            "freeze_mean_nearest_root_distance": det_freeze_mean,
            "freeze_max_nearest_root_distance": det_freeze_max
        },
        "lane2_debranges_proxy": {
            "min_separation": debranges_sep_min,
            "hb_pass_fraction": debranges_hb_fraction
        },
        "lane3_weil_positivity_proxy": {
            "kernel_u": u,
            "kernel_beta": beta,
            "min_eig_target_train": min_eig_target,
            "min_eig_pred_train": min_eig_pred,
            "min_eig_gap_pred_minus_target": min_eig_pred - min_eig_target
        },
        "lane4_prime_trace_probe": {
            "num_primes": primes.len(),
            "mean_peak_z_target": prime_z_target,
            "mean_peak_z_pred": prime_z_pred,
            "peak_z_gap": prime_z_pred - prime_z_target
        },
        "lane5_inverse_spectral_reconstruction": {
            "train_fitted_map": {
                "slope": fit_slope,
                "shift": fit_shift
            },
            "fit_hold": fit_hold_inv,
            "fit_freeze": fit_freeze_inv
        },
        "lane6_symmetry_embedding": {
            "spectrum_pair_symmetry_residual": sym_embed_resid,
            "train_fitted_map": {
                "slope": sym_slope,
                "shift": sym_shift
            },
            "fit_hold": fit_hold_sym,
            "fit_freeze": fit_freeze_sym
        },
        "lane7_hilbert_polya_bridge": {
            "status": hilbert_polya_bridge_status,
            "note": "formal proof lane remains open; this report tracks numerical/operator-side prerequisites"
        }
    });

    let out_dir = std::env::var("GUTOE_RIEMANN_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    fs::create_dir_all(&out_dir).expect("create output directory");
    let txt_path = format!("{out_dir}/riemann_phase3_all_lanes_report.txt");
    let json_path = format!("{out_dir}/riemann_phase3_all_lanes_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "Riemann Phase-3 All Lanes Report").expect("write");
    writeln!(txt, "===============================").expect("write");
    writeln!(txt, "n                            = {}", n).expect("write");
    writeln!(txt, "train/hold/freeze            = {}/{}/{}", TRAIN_N, HOLD_N, FREEZE_N).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "Baseline structural map:").expect("write");
    writeln!(txt, "  hold MAPE                 = {:.12e}", fit_hold_struct.mape).expect("write");
    writeln!(txt, "  freeze MAPE               = {:.12e}", fit_freeze_struct.mape).expect("write");

    writeln!(txt, "\nLane1 determinant-lock proxy:").expect("write");
    writeln!(txt, "  hold mean nearest root    = {:.12e}", det_hold_mean).expect("write");
    writeln!(txt, "  freeze mean nearest root  = {:.12e}", det_freeze_mean).expect("write");

    writeln!(txt, "\nLane2 de Branges proxy:").expect("write");
    writeln!(txt, "  min separation            = {:.12e}", debranges_sep_min).expect("write");
    writeln!(txt, "  HB pass fraction          = {:.12e}", debranges_hb_fraction).expect("write");

    writeln!(txt, "\nLane3 Weil positivity proxy:").expect("write");
    writeln!(txt, "  min eig target            = {:.12e}", min_eig_target).expect("write");
    writeln!(txt, "  min eig pred              = {:.12e}", min_eig_pred).expect("write");

    writeln!(txt, "\nLane4 prime trace probe:").expect("write");
    writeln!(txt, "  mean peak z target        = {:.12e}", prime_z_target).expect("write");
    writeln!(txt, "  mean peak z pred          = {:.12e}", prime_z_pred).expect("write");

    writeln!(txt, "\nLane5 inverse spectral:").expect("write");
    writeln!(txt, "  fitted slope              = {:.12e}", fit_slope).expect("write");
    writeln!(txt, "  fitted shift              = {:.12e}", fit_shift).expect("write");
    writeln!(txt, "  hold MAPE                 = {:.12e}", fit_hold_inv.mape).expect("write");
    writeln!(txt, "  freeze MAPE               = {:.12e}", fit_freeze_inv.mape).expect("write");

    writeln!(txt, "\nLane6 symmetry embedding:").expect("write");
    writeln!(txt, "  pair symmetry residual    = {:.12e}", sym_embed_resid).expect("write");
    writeln!(txt, "  hold MAPE                 = {:.12e}", fit_hold_sym.mape).expect("write");
    writeln!(txt, "  freeze MAPE               = {:.12e}", fit_freeze_sym.mape).expect("write");

    writeln!(txt, "\nLane7 Hilbert-Polya bridge:").expect("write");
    writeln!(txt, "  status                    = {}", hilbert_polya_bridge_status).expect("write");

    fs::write(
        &json_path,
        serde_json::to_string_pretty(&report).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
}
