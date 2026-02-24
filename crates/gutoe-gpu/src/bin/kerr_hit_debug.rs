use gutoe_gpu::kerr::KerrMetric;
use gutoe_gpu::tracer::{trace_photon_kerr, TraceResult};

fn main() {
    let width = 1280usize;
    let height = 720usize;
    let inc_deg = 17.0f64;
    let az_deg = 0.0f64;
    let fov_rs = 7.0f64;
    let dphi = 0.005f64;
    let disk_inner = 3.0f64;
    let disk_outer = 16.0f64;
    let k = KerrMetric::new(1.0, 0.94).expect("kerr");
    let max_phi = 60.0 * std::f64::consts::PI;
    let ca = az_deg.to_radians().cos();
    let sa = az_deg.to_radians().sin();

    let scale = 2.0 * fov_rs / width as f64;
    let mut hit = 0usize;
    let mut cap = 0usize;
    let mut esc = 0usize;
    let mut hit_top = 0usize;
    let mut hit_bottom = 0usize;
    for iy in 0..height {
        for ix in 0..width {
            let sx = (ix as f64 - 0.5 * (width as f64 - 1.0)) * scale;
            let sy = (0.5 * (height as f64 - 1.0) - iy as f64) * scale;
            let bx_raw = sx;
            let by_raw = sy;
            let bx = ca * bx_raw - sa * by_raw;
            let by = sa * bx_raw + ca * by_raw;
            let tr = trace_photon_kerr(
                &k,
                disk_inner,
                disk_outer,
                bx,
                by,
                inc_deg,
                max_phi,
                dphi,
            );
            match tr {
                TraceResult::DiskHit { .. } => {
                    hit += 1;
                    if by_raw >= 0.0 {
                        hit_top += 1;
                    } else {
                        hit_bottom += 1;
                    }
                }
                TraceResult::Captured => cap += 1,
                TraceResult::Escaped { .. } => esc += 1,
            }
        }
    }
    let n = (width * height) as f64;
    println!(
        "hits={} ({:.2}%) cap={} ({:.2}%) esc={} ({:.2}%)",
        hit,
        100.0 * hit as f64 / n,
        cap,
        100.0 * cap as f64 / n,
        esc,
        100.0 * esc as f64 / n
    );
    println!("hit_top={} hit_bottom={}", hit_top, hit_bottom);
}
