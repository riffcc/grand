//! GUTOE Black Hole Screenshot Gallery
//!
//! Renders 6 views of the GUTOE black hole from different angles using the
//! CPU ray tracer, saves as PNG files, and serves an HTML gallery via HTTP
//! on 0.0.0.0:52345.
//!
//! Views:
//!   1. 85° inclination, az=0°  — classic edge-on (photon ring visible)
//!   2. 70° inclination, az=0°  — slightly tilted
//!   3. 50° inclination, az=0°  — medium tilt
//!   4. 30° inclination, az=0°  — more face-on
//!   5. 85° inclination, az=60° — edge-on, disk rotated 60°
//!   6. 10° inclination, az=0°  — nearly face-on (shadow becomes circular)

use std::{
    fs,
    io::{Read, Write},
    net::{TcpListener, TcpStream},
    path::{Path, PathBuf},
    process::Command,
    sync::Arc,
};

use gutoe_gpu::{
    metric::GutoeMetric,
    tracer::{render, write_ppm, RenderConfig},
};

// ── View definitions ─────────────────────────────────────────────────────────

struct View {
    label: &'static str,
    slug:  &'static str,
    inc:   f64,   // inclination in degrees (90 = edge-on, 0 = face-on)
    az:    f64,   // azimuth in degrees (disk rotation on screen)
    fov:   f64,   // half-width in r_s units
}

static VIEWS: &[View] = &[
    View { label: "Edge-on (85°, classic photon ring)",   slug: "v1_edge85",   inc: 85.0, az:  0.0, fov: 12.0 },
    View { label: "Slightly tilted (70°)",                slug: "v2_tilt70",   inc: 70.0, az:  0.0, fov: 12.0 },
    View { label: "Medium tilt (50°)",                    slug: "v3_tilt50",   inc: 50.0, az:  0.0, fov: 12.0 },
    View { label: "Tilted toward face-on (30°)",          slug: "v4_tilt30",   inc: 30.0, az:  0.0, fov: 12.0 },
    View { label: "Edge-on rotated 60° (85°, az=60°)",    slug: "v5_edge_rot", inc: 85.0, az: 60.0, fov: 12.0 },
    View { label: "Nearly face-on (10°, circular shadow)", slug: "v6_face10",  inc: 10.0, az:  0.0, fov: 12.0 },
];

// ── Render & write ────────────────────────────────────────────────────────────

fn render_view(out_dir: &Path, view: &View) -> PathBuf {
    eprintln!("  rendering {}  (inc={:.0}°, az={:.0}°) …", view.label, view.inc, view.az);

    let metric = GutoeMetric::planck_units(1.0);

    let cfg = RenderConfig {
        width:          1200,
        height:         1200,
        fov_rs:         view.fov,
        inclination_deg: view.inc,
        max_phi:        30.0 * std::f64::consts::PI,
        dphi:           0.005,
    };

    // The CPU tracer doesn't take azimuth directly — the azimuth in bh_viewer
    // is applied in the shader. For the CPU path we rotate the impact parameters
    // ourselves: (bx, by) → (bx·cos(az) − by·sin(az), bx·sin(az) + by·cos(az)).
    // We achieve this by rendering with az=0 and a rotated pixel grid, OR by
    // post-rotating the image. The cleanest approach: rotate bx/by here by
    // embedding an az-aware wrapper around the pixel loop.
    //
    // Since RenderConfig has no az field, we render the view twice with a
    // thin wrapper that rotates screen coordinates before passing to trace_photon.
    let pixels = render_with_az(&metric, 3.0, 10.0, &cfg, view.az);

    let ppm_path = out_dir.join(format!("{}.ppm", view.slug));
    let png_path = out_dir.join(format!("{}.png", view.slug));

    fs::write(&ppm_path, write_ppm(&pixels, cfg.width, cfg.height))
        .expect("write ppm");

    let status = Command::new("convert")
        .arg(&ppm_path)
        .arg(&png_path)
        .status()
        .expect("ImageMagick convert not found");
    assert!(status.success(), "convert failed for {}", view.slug);

    fs::remove_file(&ppm_path).ok();
    eprintln!("  → saved {}", png_path.display());
    png_path
}

// ── Az-aware render ───────────────────────────────────────────────────────────
//
// Identical to tracer::render() but rotates each pixel's impact parameters
// by `az_deg` before tracing.  This reproduces the WGSL shader's azimuth
// rotation: (bx,by) ← R(-az)·(bx,by).

fn render_with_az(
    metric: &GutoeMetric,
    disk_inner_rs: f64,
    disk_outer_rs: f64,
    cfg: &RenderConfig,
    az_deg: f64,
) -> Vec<[u8; 3]> {
    use gutoe_gpu::tracer::{trace_photon, TraceResult};

    if az_deg.abs() < 1e-9 {
        return render(metric, disk_inner_rs, disk_outer_rs, cfg);
    }

    let r_s          = metric.r_s;
    let disk_inner   = disk_inner_rs * r_s;
    let disk_outer   = disk_outer_rs * r_s;
    let r_isco       = 3.0 * r_s;
    let sin_inc      = cfg.inclination_deg.to_radians().sin();
    let scale        = 2.0 * cfg.fov_rs * r_s / cfg.width as f64;
    let az_rad       = az_deg.to_radians();
    let (ca, sa)     = (az_rad.cos(), az_rad.sin());

    let mut pixels = vec![[0u8; 3]; cfg.width * cfg.height];
    for iy in 0..cfg.height {
        for ix in 0..cfg.width {
            let sx = (ix as f64 - 0.5 * (cfg.width  as f64 - 1.0)) * scale;
            let sy = (0.5 * (cfg.height as f64 - 1.0) - iy as f64) * scale;
            let bx_raw = sx;
            let by_raw = sy * sin_inc;
            // rotate by az in screen plane
            let bx = ca * bx_raw - sa * by_raw;
            let by = sa * bx_raw + ca * by_raw;

            let result = trace_photon(
                metric, disk_inner, disk_outer, bx, by,
                cfg.max_phi, cfg.dphi,
            );

            pixels[iy * cfg.width + ix] = match result {
                TraceResult::Captured          => [0, 0, 0],
                TraceResult::Escaped { .. }    => [5, 5, 20],
                TraceResult::DiskHit { r_eff, n_cross, .. } => disk_temp_color(r_eff, r_isco, n_cross),
            };
        }
    }
    pixels
}

fn disk_temp_color(r_eff: f64, r_isco: f64, n_cross: u32) -> [u8; 3] {
    let t_rel = (r_isco / r_eff).powf(0.75).clamp(0.005, 1.0);
    let fade  = 0.7_f64.powi(n_cross as i32 - 1);
    let b     = (t_rel * fade).clamp(0.0, 1.0);
    let r  = (255.0 * b.powf(0.4)).clamp(0.0, 255.0) as u8;
    let g  = (200.0 * b.powf(0.7)).clamp(0.0, 255.0) as u8;
    let bl = (120.0 * b.powf(1.8)).clamp(0.0, 255.0) as u8;
    [r, g, bl]
}

// ── HTML gallery ──────────────────────────────────────────────────────────────

fn build_html() -> String {
    let mut imgs = String::new();
    for v in VIEWS {
        imgs.push_str(&format!(
            r#"<figure>
  <img src="/img/{slug}.png" alt="{label}" loading="lazy">
  <figcaption>{label}</figcaption>
</figure>
"#,
            slug  = v.slug,
            label = v.label,
        ));
    }

    format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>GUTOE Black Hole — Six Views</title>
<style>
  body {{ background:#050508; color:#ddd; font-family:monospace; padding:2rem; }}
  h1   {{ color:#ffd080; margin-bottom:.25rem; }}
  p    {{ color:#aaa; margin-top:0; }}
  .gallery {{ display:flex; flex-wrap:wrap; gap:1rem; margin-top:1.5rem; }}
  figure {{ margin:0; }}
  img  {{ display:block; width:400px; height:400px; border:1px solid #333; }}
  figcaption {{ font-size:.8rem; color:#888; margin-top:.3rem; text-align:center; }}
</style>
</head>
<body>
<h1>GUTOE Black Hole — Six Views</h1>
<p>Cl(1,3) Schwarzschild metric with lattice core r_c = √C_∞ · l_P &nbsp;|&nbsp;
   CPU geodesic ray tracer &nbsp;|&nbsp; 600×600 &nbsp;|&nbsp; dphi=0.01</p>
<div class="gallery">
{imgs}
</div>
</body>
</html>
"#,
        imgs = imgs,
    )
}

// ── Minimal HTTP server ───────────────────────────────────────────────────────

fn serve_http(out_dir: Arc<PathBuf>, html: Arc<String>) {
    let addr = "0.0.0.0:52345";
    let listener = TcpListener::bind(addr).expect("bind 0.0.0.0:52345");
    eprintln!("\nGallery live at  http://10.7.1.200:52345/");
    eprintln!("Press Ctrl-C to stop.\n");

    for stream in listener.incoming() {
        match stream {
            Ok(s)  => handle_connection(s, &out_dir, &html),
            Err(e) => eprintln!("accept error: {e}"),
        }
    }
}

fn handle_connection(mut stream: TcpStream, out_dir: &Path, html: &str) {
    // Read the request line
    let mut buf = [0u8; 4096];
    let n = match stream.read(&mut buf) {
        Ok(n) => n,
        Err(_) => return,
    };
    let req = std::str::from_utf8(&buf[..n]).unwrap_or("");
    let path = req.lines().next()
        .and_then(|l| l.split_whitespace().nth(1))
        .unwrap_or("/");

    if path == "/" || path == "/index.html" {
        let body = html.as_bytes();
        let _ = write!(
            stream,
            "HTTP/1.1 200 OK\r\nContent-Type: text/html; charset=utf-8\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
            body.len()
        );
        let _ = stream.write_all(body);
        return;
    }

    if let Some(slug) = path.strip_prefix("/img/").and_then(|s| s.strip_suffix(".png")) {
        let file_path = out_dir.join(format!("{slug}.png"));
        match fs::read(&file_path) {
            Ok(data) => {
                let _ = write!(
                    stream,
                    "HTTP/1.1 200 OK\r\nContent-Type: image/png\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                    data.len()
                );
                let _ = stream.write_all(&data);
            }
            Err(_) => {
                let _ = stream.write_all(b"HTTP/1.1 404 Not Found\r\nContent-Length: 9\r\n\r\nNot Found");
            }
        }
        return;
    }

    let _ = stream.write_all(b"HTTP/1.1 404 Not Found\r\nContent-Length: 9\r\n\r\nNot Found");
}

// ── Entry point ───────────────────────────────────────────────────────────────

fn main() {
    let out_dir = PathBuf::from("/tmp/bh_renders");
    fs::create_dir_all(&out_dir).expect("create /tmp/bh_renders");

    eprintln!("GUTOE Black Hole — rendering 6 views …\n");
    for view in VIEWS {
        render_view(&out_dir, view);
    }
    eprintln!("\nAll 6 renders complete.");

    let html = Arc::new(build_html());
    let dir  = Arc::new(out_dir);
    serve_http(dir, html);
}
