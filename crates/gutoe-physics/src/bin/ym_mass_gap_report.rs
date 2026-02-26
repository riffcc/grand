//! GRAND-297/298 scaffold report:
//! finite-volume transfer-matrix spectral gap diagnostics.

use gutoe_physics::mass_gap::{
    continuum_stability_band, monotone_nonincreasing_in_volume, transfer_matrix_diagnostics,
    DenseSymmetricMatrix, VolumeGapPoint,
};
use std::fs::{self, File};
use std::io::Write;

fn toy_transfer_matrix_for_volume(l: usize) -> DenseSymmetricMatrix {
    // Synthetic finite-volume spectrum:
    // m1(L) = m_inf + c / L^3, m2(L) heavier channel.
    let l3 = (l as f64).powi(3);
    let m_inf = 0.68;
    let c = 180.0;
    let m1 = m_inf + c / l3;
    let m2 = m1 + 0.55;
    let d0 = 1.0;
    let d1 = (-m1).exp();
    let d2 = (-m2).exp();
    let eps = 0.01 * d2;

    DenseSymmetricMatrix::from_rows(&[
        vec![d0, eps, eps],
        vec![eps, d1, eps],
        vec![eps, eps, d2],
    ])
    .expect("valid 3x3 matrix")
}

fn main() {
    let a_t = 1.0;
    let tol = 1e-12;
    let max_iters = 10_000;
    let volumes = [8usize, 10, 12, 14, 16];

    let mut points = Vec::<VolumeGapPoint>::new();
    let mut rows = Vec::new();

    for &l in &volumes {
        let t = toy_transfer_matrix_for_volume(l);
        let d = transfer_matrix_diagnostics(&t, a_t, max_iters, tol).expect("diagnostics");
        let g = d.gap.clone().expect("gap estimate");
        let err = (g.gap_est - g.gap_lower_bound.unwrap_or(g.gap_est)).abs();
        points.push(VolumeGapPoint {
            volume_l3: l * l * l,
            gap_est: g.gap_est,
            gap_err: err,
        });
        rows.push((l, d, g, err));
    }

    let monotone_ok = monotone_nonincreasing_in_volume(&points, 1e-6);
    let band = continuum_stability_band(&points);

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let txt_path = format!("{out_dir}/ym_mass_gap_report.txt");
    let json_path = format!("{out_dir}/ym_mass_gap_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "GRAND-297/298 mass-gap scaffold report").expect("write");
    writeln!(txt, "a_t={:.6}, tol={:.2e}", a_t, tol).expect("write");
    writeln!(txt).expect("write");
    for (l, d, g, err) in &rows {
        writeln!(
            txt,
            "L={:>2}, dim={}, symmetric={}, nonnegative={}, gersh_lb={:.6e}, λ0={:.9}, λ1={:.9}, gap={:.9}, gap_lb={:.9}, err={:.3e}",
            l,
            d.dim,
            d.symmetric,
            d.entrywise_nonnegative,
            d.gershgorin_lower_bound,
            g.lambda0_est,
            g.lambda1_est,
            g.gap_est,
            g.gap_lower_bound.unwrap_or(f64::NAN),
            err
        )
        .expect("write");
    }
    writeln!(txt).expect("write");
    writeln!(txt, "monotone_nonincreasing_in_volume={}", monotone_ok).expect("write");
    if let Some((lo, hi)) = band {
        writeln!(txt, "continuum_stability_band=[{:.9}, {:.9}]", lo, hi).expect("write");
    } else {
        writeln!(txt, "continuum_stability_band=null").expect("write");
    }

    let mut json = File::create(&json_path).expect("create json");
    writeln!(json, "{{").expect("write");
    writeln!(json, "  \"a_t\": {:.9},", a_t).expect("write");
    writeln!(json, "  \"tol\": {:.3e},", tol).expect("write");
    writeln!(json, "  \"rows\": [").expect("write");
    for (idx, (l, d, g, err)) in rows.iter().enumerate() {
        writeln!(
            json,
            "    {{\"L\":{},\"volume_l3\":{},\"dim\":{},\"symmetric\":{},\"entrywise_nonnegative\":{},\"gershgorin_lb\":{:.12e},\"lambda0\":{:.12e},\"lambda1\":{:.12e},\"lambda0_residual\":{:.12e},\"lambda1_residual\":{:.12e},\"gap_est\":{:.12e},\"gap_lower_bound\":{:.12e},\"gap_err\":{:.12e}}}{}",
            l,
            l * l * l,
            d.dim,
            d.symmetric,
            d.entrywise_nonnegative,
            d.gershgorin_lower_bound,
            g.lambda0_est,
            g.lambda1_est,
            g.lambda0_residual,
            g.lambda1_residual,
            g.gap_est,
            g.gap_lower_bound.unwrap_or(f64::NAN),
            err,
            if idx + 1 == rows.len() { "" } else { "," }
        )
        .expect("write");
    }
    writeln!(json, "  ],").expect("write");
    writeln!(json, "  \"monotone_nonincreasing_in_volume\": {},", monotone_ok).expect("write");
    match band {
        Some((lo, hi)) => {
            writeln!(
                json,
                "  \"continuum_stability_band\": {{\"lo\": {:.12e}, \"hi\": {:.12e}}}",
                lo, hi
            )
            .expect("write");
        }
        None => {
            writeln!(json, "  \"continuum_stability_band\": null").expect("write");
        }
    }
    writeln!(json, "}}").expect("write");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "YM scaffold: monotone={} continuum_band={:?}",
        monotone_ok, band
    );
}
