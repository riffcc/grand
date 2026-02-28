use std::env;
use std::f64::consts::PI;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::{SystemTime, UNIX_EPOCH};

use gutoe_physics::{
    evaluate_baryogenesis_gate, evaluate_universe_gate_with_depth, BaryogenesisWindows,
    UniverseAssumptions, UniverseSimulationDepth, UniverseWindows,
};

const PLANCK_TIME: f64 = 5.391_247e-44;
const T_INFLATION: f64 = 1e-36;
const T_EW_BREAK: f64 = 1e-12;
const T_QCD: f64 = 1e-6;
const T_NEUTRINO: f64 = 1e-2;

#[derive(Debug, Clone)]
struct Config {
    out_dir: PathBuf,
    work_dir: PathBuf,
    clip_seconds: f64,
    first_seconds: f64,
    fps: usize,
    size: usize,
    width: usize,
    height: usize,
    grid: GridKind,
    skip_gif: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GridKind {
    Square,
    Hex,
}

impl GridKind {
    fn as_str(self) -> &'static str {
        match self {
            Self::Square => "square",
            Self::Hex => "hex",
        }
    }
}

#[derive(Clone, Copy)]
struct FrameState {
    temperature_k: f64,
    omega_r: f64,
    omega_m: f64,
    omega_lambda: f64,
    z: f64,
}

fn smoothstep(x: f64, c: f64, w: f64) -> f64 {
    1.0 / (1.0 + (-(x - c) / w.max(1e-9)).exp())
}

fn temp_radiation_era_k(t: f64) -> f64 {
    1.160_45e10 / t.max(1e-40).sqrt()
}

fn phase_fractions(t: f64) -> (f64, f64, f64, f64) {
    let logt = t.max(PLANCK_TIME).log10();
    let s_infl = smoothstep(logt, T_INFLATION.log10(), 0.35);
    let s_ew = smoothstep(logt, T_EW_BREAK.log10(), 0.25);
    let s_qcd = smoothstep(logt, T_QCD.log10(), 0.30);
    let foam = 1.0 - s_infl;
    let inflation = s_infl * (1.0 - s_ew);
    let plasma = s_ew * (1.0 - s_qcd);
    let hadronic = s_qcd;
    let norm = foam + inflation + plasma + hadronic;
    if norm <= 0.0 {
        return (1.0, 0.0, 0.0, 0.0);
    }
    (foam / norm, inflation / norm, plasma / norm, hadronic / norm)
}

fn stage_name(t: f64) -> &'static str {
    if t < T_INFLATION {
        "Quantum foam"
    } else if t < T_EW_BREAK {
        "Inflation / expansion"
    } else if t < T_QCD {
        "Electroweak broken phase"
    } else if t < T_NEUTRINO {
        "QCD confinement onset"
    } else {
        "Approaching neutrino decoupling"
    }
}

fn clamp01(x: f64) -> f64 {
    x.clamp(0.0, 1.0)
}

fn colormap(name: &str, v: f64) -> [u8; 3] {
    let x = clamp01(v);
    match name {
        "coolwarm" => {
            let r = clamp01(0.25 + 0.95 * x);
            let g = clamp01(0.3 + 0.6 * (1.0 - (2.0 * x - 1.0).abs()));
            let b = clamp01(0.25 + 0.95 * (1.0 - x));
            [(255.0 * r) as u8, (255.0 * g) as u8, (255.0 * b) as u8]
        }
        "magma" => {
            let r = clamp01(1.2 * x.powf(0.7));
            let g = clamp01(0.8 * x.powf(1.4));
            let b = clamp01(0.5 * x.powf(2.0));
            [(255.0 * r) as u8, (255.0 * g) as u8, (255.0 * b) as u8]
        }
        "inferno" => {
            let r = clamp01(1.3 * x.powf(0.65));
            let g = clamp01(0.9 * x.powf(1.2));
            let b = clamp01(0.45 * x.powf(2.4));
            [(255.0 * r) as u8, (255.0 * g) as u8, (255.0 * b) as u8]
        }
        "plasma" => {
            let r = clamp01(0.75 + 0.35 * x);
            let g = clamp01(0.2 + 0.8 * (1.0 - (x - 0.55).abs()));
            let b = clamp01(0.4 + 0.6 * (1.0 - x));
            [(255.0 * r) as u8, (255.0 * g) as u8, (255.0 * b) as u8]
        }
        "viridis" => {
            let r = clamp01(0.2 + 0.8 * x.powf(1.7));
            let g = clamp01(0.1 + 1.0 * x.powf(0.8));
            let b = clamp01(0.35 + 0.55 * (1.0 - x).powf(0.9));
            [(255.0 * r) as u8, (255.0 * g) as u8, (255.0 * b) as u8]
        }
        "turbo" => {
            let r = clamp01(1.15 * (1.0 - (x - 0.8).abs()));
            let g = clamp01(1.1 * (1.0 - (x - 0.5).abs()));
            let b = clamp01(1.1 * (1.0 - (x - 0.2).abs()));
            [(255.0 * r) as u8, (255.0 * g) as u8, (255.0 * b) as u8]
        }
        "cubehelix" => {
            let a = 2.0 * PI * (0.5 + 1.3 * x);
            let amp = 0.5 * x * (1.0 - x);
            let r = clamp01(x + amp * (-0.14861 * a.cos() + 1.78277 * a.sin()));
            let g = clamp01(x + amp * (-0.29227 * a.cos() - 0.90649 * a.sin()));
            let b = clamp01(x + amp * (1.97294 * a.cos()));
            [(255.0 * r) as u8, (255.0 * g) as u8, (255.0 * b) as u8]
        }
        "Spectral_r" => {
            let u = 1.0 - x;
            let r = clamp01(0.3 + 1.1 * (1.0 - (u - 0.8).abs()));
            let g = clamp01(0.2 + 1.0 * (1.0 - (u - 0.5).abs()));
            let b = clamp01(0.2 + 1.0 * (1.0 - (u - 0.2).abs()));
            [(255.0 * r) as u8, (255.0 * g) as u8, (255.0 * b) as u8]
        }
        _ => {
            let g = clamp01(0.2 + 0.9 * x);
            [(255.0 * g) as u8, (255.0 * g) as u8, (255.0 * g) as u8]
        }
    }
}

fn parse_args() -> Result<Config, String> {
    let mut out_dir = PathBuf::from("/tmp/bh_renders/first_second_multispectral");
    let mut work_dir = env::temp_dir();
    let mut clip_seconds = 120.0;
    let mut first_seconds = 1.0;
    let mut fps = 24usize;
    let mut size = 420usize;
    let mut width = 1920usize;
    let mut height = 1080usize;
    let mut grid = GridKind::Square;
    let mut skip_gif = false;
    let mut it = env::args().skip(1);
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--help" | "-h" => {
                println!(
                    "universe_first_second_multispectral_movie\n\
                     \n\
                     Flags:\n\
                     \n\
                     --out-dir <PATH>        Output directory (default: /tmp/bh_renders/first_second_multispectral)\n\
                     --work-dir <PATH>       Working directory for intermediate frames (default: system temp)\n\
                     --clip-seconds <FLOAT>  Full movie wall-clock duration (default: 120)\n\
                     --first-seconds <FLOAT> Universe-time horizon for the full clip (default: 1.0)\n\
                     --fps <INT>             Frames per second (default: 24)\n\
                     --size <INT>            Internal field grid size (default: 420)\n\
                     --width <INT>           Output width (default: 1920)\n\
                     --height <INT>          Output height (default: 1080)\n\
                     --grid <square|hex>    Sampling geometry (default: square)\n\
                     --skip-gif              Skip GIF artifact generation\n\
                     --help, -h              Show this help"
                );
                std::process::exit(0);
            }
            "--out-dir" => out_dir = PathBuf::from(it.next().ok_or("missing --out-dir value")?),
            "--work-dir" => work_dir = PathBuf::from(it.next().ok_or("missing --work-dir value")?),
            "--clip-seconds" => {
                clip_seconds = it
                    .next()
                    .ok_or("missing --clip-seconds value")?
                    .parse::<f64>()
                    .map_err(|e| format!("invalid --clip-seconds: {e}"))?;
            }
            "--first-seconds" => {
                first_seconds = it
                    .next()
                    .ok_or("missing --first-seconds value")?
                    .parse::<f64>()
                    .map_err(|e| format!("invalid --first-seconds: {e}"))?;
            }
            "--fps" => {
                fps = it
                    .next()
                    .ok_or("missing --fps value")?
                    .parse::<usize>()
                    .map_err(|e| format!("invalid --fps: {e}"))?;
            }
            "--size" => {
                size = it
                    .next()
                    .ok_or("missing --size value")?
                    .parse::<usize>()
                    .map_err(|e| format!("invalid --size: {e}"))?;
            }
            "--width" => {
                width = it
                    .next()
                    .ok_or("missing --width value")?
                    .parse::<usize>()
                    .map_err(|e| format!("invalid --width: {e}"))?;
            }
            "--height" => {
                height = it
                    .next()
                    .ok_or("missing --height value")?
                    .parse::<usize>()
                    .map_err(|e| format!("invalid --height: {e}"))?;
            }
            "--grid" => {
                let v = it.next().ok_or("missing --grid value")?;
                grid = match v.as_str() {
                    "square" => GridKind::Square,
                    "hex" => GridKind::Hex,
                    _ => return Err(format!("invalid --grid: {v} (expected square|hex)")),
                };
            }
            "--skip-gif" => skip_gif = true,
            _ => return Err(format!("unknown arg: {arg}")),
        }
    }
    Ok(Config {
        out_dir,
        work_dir,
        clip_seconds,
        first_seconds,
        fps,
        size,
        width,
        height,
        grid,
        skip_gif,
    })
}

fn build_frame_states(nframes: usize, sim_end_s: f64) -> Result<Vec<FrameState>, String> {
    let depth = UniverseSimulationDepth {
        history_points: 4096,
        history_z_max: 1.0e12,
        integral_z_max: 1.0e12,
    };
    let score = evaluate_universe_gate_with_depth(
        UniverseAssumptions::default(),
        UniverseWindows::default(),
        depth,
    );
    let mut hist = score.history.clone();
    hist.sort_by(|a, b| a.age_seconds.total_cmp(&b.age_seconds));
    if hist.is_empty() {
        return Err("universe history is empty".to_string());
    }

    let mut out = Vec::with_capacity(nframes);
    for i in 0..nframes {
        let u = if nframes <= 1 {
            0.0
        } else {
            i as f64 / (nframes - 1) as f64
        };
        let logt = PLANCK_TIME.log10() + u * (sim_end_s.log10() - PLANCK_TIME.log10());
        let t = 10f64.powf(logt);

        // Below the simulated history floor, clamp to earliest solved row but keep
        // radiation-era T(t) scaling for honest early-time temperature trend.
        let idx = hist.partition_point(|r| r.age_seconds < t);
        let (temperature_k, omega_r, omega_m, omega_lambda) = if idx == 0 {
            let r = &hist[0];
            (temp_radiation_era_k(t), r.omega_r, r.omega_m, r.omega_lambda)
        } else if idx >= hist.len() {
            let r = &hist[hist.len() - 1];
            (r.temperature_k, r.omega_r, r.omega_m, r.omega_lambda)
        } else {
            let a = &hist[idx - 1];
            let b = &hist[idx];
            let denom = (b.age_seconds - a.age_seconds).max(1e-30);
            let w = ((t - a.age_seconds) / denom).clamp(0.0, 1.0);
            (
                a.temperature_k + w * (b.temperature_k - a.temperature_k),
                a.omega_r + w * (b.omega_r - a.omega_r),
                a.omega_m + w * (b.omega_m - a.omega_m),
                a.omega_lambda + w * (b.omega_lambda - a.omega_lambda),
            )
        };
        out.push(FrameState {
            temperature_k,
            omega_r,
            omega_m,
            omega_lambda,
            z: if idx == 0 {
                hist[0].z
            } else if idx >= hist.len() {
                hist[hist.len() - 1].z
            } else {
                let a = &hist[idx - 1];
                let b = &hist[idx];
                let denom = (b.age_seconds - a.age_seconds).max(1e-30);
                let w = ((t - a.age_seconds) / denom).clamp(0.0, 1.0);
                a.z + w * (b.z - a.z)
            },
        });
    }
    Ok(out)
}

fn ffmpeg_path() -> Result<String, String> {
    let out = Command::new("sh")
        .arg("-lc")
        .arg("command -v ffmpeg")
        .output()
        .map_err(|e| format!("failed to locate ffmpeg: {e}"))?;
    if !out.status.success() {
        return Err("ffmpeg not found in PATH".to_string());
    }
    Ok(String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn smooth_box(grid: &[f64], n: usize, radius: usize, x: usize, y: usize) -> f64 {
    let mut acc = 0.0;
    let mut cnt = 0usize;
    let y0 = y.saturating_sub(radius);
    let y1 = (y + radius).min(n - 1);
    let x0 = x.saturating_sub(radius);
    let x1 = (x + radius).min(n - 1);
    for yy in y0..=y1 {
        for xx in x0..=x1 {
            acc += grid[yy * n + xx];
            cnt += 1;
        }
    }
    acc / cnt as f64
}

fn render_frame(
    path: &Path,
    i: usize,
    nframes: usize,
    size: usize,
    sim_end_s: f64,
    state: FrameState,
    cp_asymmetry: f64,
    width: usize,
    height: usize,
    grid: GridKind,
) -> Result<(), String> {
    let u = if nframes <= 1 { 0.0 } else { i as f64 / (nframes - 1) as f64 };
    let logt = PLANCK_TIME.log10() + u * (sim_end_s.log10() - PLANCK_TIME.log10());
    let t = 10f64.powf(logt);
    let temp_k = state.temperature_k.max(1.0);
    let temp_norm = clamp01((temp_k.max(1.0).log10() - 10.0) / 18.0);
    let (mut foam_w, mut infl_w, mut plasma_w, mut had_w) = phase_fractions(t);
    // Physics-lane coupling: use solved energy fractions to modulate phase blend.
    let rad_boost = state.omega_r.clamp(0.0, 1.0);
    let matter_boost = state.omega_m.clamp(0.0, 1.0);
    let lambda_boost = state.omega_lambda.clamp(0.0, 1.0);
    foam_w *= 0.7 + 0.6 * rad_boost;
    infl_w *= 0.7 + 0.4 * lambda_boost;
    plasma_w *= 0.7 + 0.6 * rad_boost;
    had_w *= 0.7 + 0.6 * matter_boost;
    let norm = (foam_w + infl_w + plasma_w + had_w).max(1e-12);
    foam_w /= norm;
    infl_w /= norm;
    plasma_w /= norm;
    had_w /= norm;

    // Force structural dependence on lane state so timeline horizon changes are visible.
    let z_norm = ((state.z + 1.0).log10() / 12.0).clamp(0.0, 1.0);
    let late_time = clamp01((t.log10() + 2.0) / 4.0); // turns on near t~1e-2..1e2
    let morph = (0.55 * state.omega_m + 0.35 * state.omega_lambda + 0.10 * late_time).clamp(0.0, 1.0);
    let infl_prog = clamp01((logt - PLANCK_TIME.log10()) / (T_EW_BREAK.log10() - PLANCK_TIME.log10()));
    let swirl_freq = 7.0 + 9.0 * morph + 6.0 * (1.0 - z_norm);
    let ring_center = 0.08 + 0.58 * infl_prog * (1.0 - 0.35 * morph) + 0.20 * morph;
    let turb_gain = 0.10 + 0.55 * (state.omega_r + 0.35 * state.omega_m).clamp(0.0, 1.0);
    let post_qcd_shell = clamp01((t.log10() + 6.0) / 4.0);
    let phase = (logt - PLANCK_TIME.log10()) * (1.8 + 1.2 * morph);
    let handed = if cp_asymmetry >= 0.0 { 1.0 } else { -1.0 };
    let cp_mag = cp_asymmetry.abs().clamp(0.0, 1.0);

    let mut base = vec![0.0; size * size];
    for y in 0..size {
        let y0 = -1.0 + 2.0 * (y as f64 / (size - 1) as f64);
        for x in 0..size {
            let x0 = -1.0 + 2.0 * (x as f64 / (size - 1) as f64);
            let (xf, yf) = match grid {
                GridKind::Square => (x0, y0),
                GridKind::Hex => {
                    // Axial-like remap in pixel space (pointy-top hex embedding).
                    let q = x0 - 0.5 * y0;
                    let r_hex = 0.866_025_403_784_438_6 * y0;
                    (q, r_hex)
                }
            };
            let r = (xf * xf + yf * yf).sqrt();
            let a = yf.atan2(xf);
            let foam = ((22.0 * xf + 17.0 * yf + phase).sin()
                + (31.0 * xf - 12.0 * yf + 1.7 * phase).sin()
                + (19.0 * (xf + yf) - 0.9 * phase).sin()
                + 3.0)
                / 6.0;
            let ring = (-24.0 * (r - ring_center).powi(2)).exp();
            let rays = 0.5 + 0.5 * (swirl_freq * a * handed + 3.0 * phase).sin();
            let core = (-20.0 * r * r).exp();
            let infl = clamp01(0.7 * ring * rays + 0.6 * core);
            let plasma_core = (-4.0 * r * r).exp();
            let turb = 0.45 * (8.0 * xf + 6.5 * yf + handed * phase).sin()
                + 0.35 * (14.0 * r - 1.4 * handed * phase).sin()
                + 0.20 * (11.0 * xf - 13.0 * yf + 0.7 * handed * phase).sin();
            let shell = (-16.0 * (r - 0.58).powi(2)).exp();
            let mut plasma = plasma_core * (1.0 + 1.6 * temp_norm + 0.35 * state.omega_r)
                + (0.12 + 0.35 * post_qcd_shell) * shell
                + turb_gain * (1.0 + 0.5 * cp_mag) * turb;
            if r > 1.0 {
                plasma = 0.0;
            }
            let idx = y * size + x;
            base[idx] = (0.95 * foam_w * foam
                + 0.95 * infl_w * infl
                + 1.05 * plasma_w * plasma
                + 0.30 * had_w * smooth_box(&base, size, 1, x, y))
                .max(0.0);
        }
    }
    let bmin = base.iter().copied().fold(f64::INFINITY, f64::min);
    let bmax = base.iter().copied().fold(f64::NEG_INFINITY, f64::max).max(bmin + 1e-9);
    for v in &mut base {
        *v = (*v - bmin) / (bmax - bmin);
    }

    let mut radio = vec![0.0; size * size];
    let mut bmag = vec![0.0; size * size];
    let mut entropy = vec![0.0; size * size];
    for y in 0..size {
        for x in 0..size {
            let idx = y * size + x;
            radio[idx] = smooth_box(&base, size, 3, x, y);
            let xl = x.saturating_sub(1);
            let xr = (x + 1).min(size - 1);
            let yu = y.saturating_sub(1);
            let yd = (y + 1).min(size - 1);
            let gx = base[y * size + xr] - base[y * size + xl];
            let gy = base[yd * size + x] - base[yu * size + x];
            bmag[idx] = (gx * gx + gy * gy).sqrt();
            let p = base[idx].max(1e-9);
            entropy[idx] = -p * p.ln();
        }
    }
    let bmag_max = bmag.iter().copied().fold(0.0, f64::max).max(1e-9);
    let ent_max = entropy.iter().copied().fold(0.0, f64::max).max(1e-9);
    for v in &mut bmag {
        *v /= bmag_max;
    }
    for v in &mut entropy {
        *v /= ent_max;
    }

    let channels: [(&str, &str, Box<dyn Fn(usize) -> f64>); 9] = [
        ("Radio", "cividis", Box::new(|idx| radio[idx])),
        ("Microwave", "coolwarm", Box::new(|idx| base[idx].sqrt())),
        ("Infrared", "magma", Box::new(|idx| base[idx].powf(1.35))),
        ("Visible", "inferno", Box::new(|idx| base[idx])),
        ("Ultraviolet", "plasma", Box::new(|idx| clamp01((base[idx] - 0.45) / 0.55))),
        ("X-ray", "viridis", Box::new(|idx| clamp01((base[idx] - 0.62) / 0.38).powf(1.15))),
        ("Gamma", "turbo", Box::new(|idx| clamp01((base[idx] - 0.8) / 0.2).powf(1.35))),
        ("Magnetic |B|", "cubehelix", Box::new(|idx| bmag[idx])),
        ("Entropy Proxy", "Spectral_r", Box::new(|idx| entropy[idx])),
    ];

    let w = width.max(3);
    let h = height.max(3);
    let tw = w / 3;
    let th = h / 3;
    let mut rgb = vec![0u8; w * h * 3];
    for ty in 0..3 {
        for tx in 0..3 {
            let cidx = ty * 3 + tx;
            for yy in 0..th {
                let sy = yy * size / th;
                for xx in 0..tw {
                    let sx = xx * size / tw;
                    let src_idx = sy * size + sx;
                    let v = channels[cidx].2(src_idx);
                    let color = colormap(channels[cidx].1, v);
                    let px = tx * tw + xx;
                    let py = ty * th + yy;
                    let di = (py * w + px) * 3;
                    rgb[di] = color[0];
                    rgb[di + 1] = color[1];
                    rgb[di + 2] = color[2];
                }
            }
        }
    }

    let mut ppm = format!("P6\n{} {}\n255\n", w, h).into_bytes();
    ppm.extend_from_slice(&rgb);
    fs::write(path, ppm).map_err(|e| format!("write frame failed: {e}"))?;
    let _ = stage_name(t);
    Ok(())
}

fn run_ffmpeg(ffmpeg: &str, args: &[&str]) -> Result<(), String> {
    let st = Command::new(ffmpeg)
        .args(args)
        .status()
        .map_err(|e| format!("ffmpeg launch failed: {e}"))?;
    if !st.success() {
        return Err(format!("ffmpeg failed with status {st}"));
    }
    Ok(())
}

fn main() -> Result<(), String> {
    let cfg = parse_args()?;
    let baryo = evaluate_baryogenesis_gate(BaryogenesisWindows::default());
    let cp_asymmetry = ((baryo.eta_predicted - baryo.eta_observed) / baryo.eta_observed)
        .clamp(-1.0, 1.0);
    let sim_end_s = cfg.first_seconds.max(PLANCK_TIME);
    fs::create_dir_all(&cfg.out_dir).map_err(|e| format!("create out dir failed: {e}"))?;
    let ffmpeg = ffmpeg_path()?;
    let nframes = (cfg.clip_seconds * cfg.fps as f64).round().max(2.0) as usize;
    let frame_states = Arc::new(build_frame_states(nframes, sim_end_s)?);

    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| format!("clock error: {e}"))?
        .as_nanos();
    fs::create_dir_all(&cfg.work_dir).map_err(|e| format!("create work dir failed: {e}"))?;
    let frame_dir = cfg.work_dir.join(format!("first_second_multi_frames_{nonce}"));
    fs::create_dir_all(&frame_dir).map_err(|e| format!("create temp dir failed: {e}"))?;

    let worker_count = thread::available_parallelism().map(|n| n.get()).unwrap_or(4).max(1);
    let next = Arc::new(AtomicUsize::new(0));
    let mut handles = Vec::with_capacity(worker_count);
    for _ in 0..worker_count {
        let next = Arc::clone(&next);
        let frame_dir = frame_dir.clone();
        let frame_states = Arc::clone(&frame_states);
        let size = cfg.size;
        handles.push(thread::spawn(move || -> Result<(), String> {
            loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= nframes {
                    break;
                }
                let frame = frame_dir.join(format!("frame_{i:05}.ppm"));
                render_frame(
                    &frame,
                    i,
                    nframes,
                    size,
                    sim_end_s,
                    frame_states[i],
                    cp_asymmetry,
                    cfg.width,
                    cfg.height,
                    cfg.grid,
                )?;
            }
            Ok(())
        }));
    }
    for h in handles {
        let res = h.join().map_err(|_| "worker thread panicked".to_string())?;
        res?;
    }

    let mp4 = cfg.out_dir.join("universe_first_second_multispectral_1080p.mp4");
    run_ffmpeg(
        &ffmpeg,
        &[
            "-y",
            "-framerate",
            &cfg.fps.to_string(),
            "-i",
            &frame_dir.join("frame_%05d.ppm").to_string_lossy(),
            "-c:v",
            "libx264",
            "-preset",
            "slow",
            "-crf",
            "18",
            "-pix_fmt",
            "yuv420p",
            &mp4.to_string_lossy(),
        ],
    )?;

    let gif_path = cfg.out_dir.join("universe_first_second_multispectral_1080p.gif");
    if !cfg.skip_gif {
        let palette = frame_dir.join("palette.png");
        run_ffmpeg(
            &ffmpeg,
            &[
                "-y",
                "-i",
                &frame_dir.join("frame_%05d.ppm").to_string_lossy(),
                "-vf",
                "palettegen",
                &palette.to_string_lossy(),
            ],
        )?;
        run_ffmpeg(
            &ffmpeg,
            &[
                "-y",
                "-framerate",
                &cfg.fps.to_string(),
                "-i",
                &frame_dir.join("frame_%05d.ppm").to_string_lossy(),
                "-i",
                &palette.to_string_lossy(),
                "-lavfi",
                "paletteuse",
                &gif_path.to_string_lossy(),
            ],
        )?;
    }

    let summary_path = cfg.out_dir.join("universe_first_second_multispectral_summary.json");
    let summary = format!(
        "{{\n  \"sim_start_s\": {:.9e},\n  \"sim_end_s\": {:.9e},\n  \"clip_seconds\": {:.3},\n  \"first_seconds\": {:.9e},\n  \"fps\": {},\n  \"frames\": {},\n  \"grid\": \"{}\",\n  \"baryogenesis\": {{\"eta_predicted\": {:.9e}, \"eta_observed\": {:.9e}, \"cp_asymmetry\": {:.9e}}},\n  \"artifacts\": {{\n    \"mp4\": \"{}\",\n    \"gif\": {}\n  }}\n}}\n",
        PLANCK_TIME,
        sim_end_s,
        cfg.clip_seconds,
        cfg.first_seconds,
        cfg.fps,
        nframes,
        cfg.grid.as_str(),
        baryo.eta_predicted,
        baryo.eta_observed,
        cp_asymmetry,
        mp4.display(),
        if cfg.skip_gif {
            "null".to_string()
        } else {
            format!("\"{}\"", gif_path.display())
        },
    );
    fs::write(&summary_path, summary).map_err(|e| format!("write summary failed: {e}"))?;

    let _ = fs::remove_dir_all(&frame_dir);
    println!("wrote {}", mp4.display());
    if !cfg.skip_gif {
        println!("wrote {}", gif_path.display());
    }
    println!("wrote {}", summary_path.display());
    Ok(())
}
