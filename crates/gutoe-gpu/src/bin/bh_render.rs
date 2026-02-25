//! GUTOE Black Hole Gallery — multi-view ray tracer
//!
//! Renders 10 views of the GUTOE Schwarzschild black hole from the CPU ray
//! tracer, saves as PNG files, and serves an HTML gallery via HTTP on
//! 0.0.0.0:52345.  Each render is 1200×1200, dphi=0.005 (0.003 for sub-ring).
//!
//! New in this version:
//!   • Relativistic Doppler beaming (M87*, Sgr A*) — approaching side boosted by D⁴
//!   • Photon sub-ring coloring — n=1 (orange), n=2 (cyan), n=3 (purple)
//!   • GR comparison — pure Schwarzschild (l_P=0) side-by-side with GUTOE
//!   • Reinhard tone mapping — handles Doppler over-exposure gracefully

use std::{
    collections::BTreeMap,
    fs,
    io::{Read, Write},
    net::{TcpListener, TcpStream},
    path::{Path, PathBuf},
    process::Command,
    sync::{Arc, OnceLock},
};

use rayon::prelude::*;

use gutoe_gpu::{
    geodesic3d::{reduce_3d_to_axisym, trace_photon_3d_schwarzschild, CameraFrame, Vec3},
    kerr::KerrMetric,
    metric::{GutoeMetric, C_INF, LAMBDA_QG, WATSON_SC},
    synchrotron::{band_tint, band_weight_with_exposure, RenderSpectrum as SpectralBand},
    transfer::{covariant_absorption, covariant_emissivity, transfer_step},
    tracer::{
        b_critical, trace_photon, trace_photon_interior, trace_photon_interior_core,
        trace_photon_kerr, trace_photon_kerr_activity, trace_photon_kerr_debug,
        trace_photon_kerr_wave_branches, trace_photon_kerr_wave_samples, write_ppm,
        RenderConfig, TraceResult,
    },
};

#[derive(Debug, Clone, Copy)]
struct KerrImageMetrics {
    shadow_cx: f64,
    shadow_cy: f64,
    shadow_radius: f64,
    shadow_diameter: f64,
    ring_thickness: f64,
    flux_left: f64,
    flux_right: f64,
    flux_asymmetry: f64,
    closure_phase_deg: f64,
}

struct StarMap {
    w: usize,
    h: usize,
    rgb: Vec<[u8; 3]>,
}

impl StarMap {
    fn sample(&self, u: f64, v: f64) -> [u8; 3] {
        let uu = u.rem_euclid(1.0);
        let vv = v.clamp(0.0, 1.0);
        let x = uu * (self.w as f64 - 1.0);
        let y = vv * (self.h as f64 - 1.0);
        let x0 = x.floor() as usize;
        let y0 = y.floor() as usize;
        let x1 = (x0 + 1).min(self.w - 1);
        let y1 = (y0 + 1).min(self.h - 1);
        let tx = x - x0 as f64;
        let ty = y - y0 as f64;
        let p00 = self.rgb[y0 * self.w + x0];
        let p10 = self.rgb[y0 * self.w + x1];
        let p01 = self.rgb[y1 * self.w + x0];
        let p11 = self.rgb[y1 * self.w + x1];
        let lerp = |a: f64, b: f64, t: f64| a * (1.0 - t) + b * t;
        let ch = |c: usize| -> u8 {
            let a = lerp(p00[c] as f64, p10[c] as f64, tx);
            let b = lerp(p01[c] as f64, p11[c] as f64, tx);
            lerp(a, b, ty).round().clamp(0.0, 255.0) as u8
        };
        [ch(0), ch(1), ch(2)]
    }
}

static STARMAP: OnceLock<Option<Arc<StarMap>>> = OnceLock::new();
static TRUTH_MODE: OnceLock<TruthMode> = OnceLock::new();

fn starmap() -> Option<&'static Arc<StarMap>> {
    STARMAP
        .get_or_init(|| {
            let path = std::env::var("BH_STARMAP_PATH").ok()?;
            let img = image::open(&path).ok()?.to_rgb8();
            let (w, h) = img.dimensions();
            let mut rgb = Vec::with_capacity((w * h) as usize);
            for p in img.pixels() {
                rgb.push([p[0], p[1], p[2]]);
            }
            eprintln!("    starmap=loaded {} ({}x{})", path, w, h);
            Some(Arc::new(StarMap {
                w: w as usize,
                h: h as usize,
                rgb,
            }))
        })
        .as_ref()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TruthMode {
    Off,
    Class,
    Impact,
    Crossings,
    Angle,
    Lean,
    Transfer,
    Tau,
}

impl TruthMode {
    fn from_env() -> Self {
        let raw = std::env::var("BH_TRUTH_MODE")
            .unwrap_or_default()
            .trim()
            .to_ascii_lowercase();
        match raw.as_str() {
            "class" | "classes" | "hitclass" => Self::Class,
            "impact" | "b" | "impact_radius" => Self::Impact,
            "cross" | "crossings" | "n_cross" => Self::Crossings,
            "angle" | "phi" | "phi_total" => Self::Angle,
            "lean" | "oracle" | "lean_oracle" => Self::Lean,
            "transfer" | "g4" | "g^4" => Self::Transfer,
            "tau" | "optical_depth" => Self::Tau,
            _ => Self::Off,
        }
    }

    fn as_label(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::Class => "class",
            Self::Impact => "impact",
            Self::Crossings => "crossings",
            Self::Angle => "angle",
            Self::Lean => "lean",
            Self::Transfer => "transfer",
            Self::Tau => "tau",
        }
    }
}

fn truth_mode() -> TruthMode {
    *TRUTH_MODE.get_or_init(TruthMode::from_env)
}

fn false_color_ramp01(x: f64) -> [u8; 3] {
    let u = x.clamp(0.0, 1.0);
    let r = (1.5 - (4.0 * u - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * u - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * u - 1.0).abs()).clamp(0.0, 1.0);
    [(255.0 * r) as u8, (255.0 * g) as u8, (255.0 * b) as u8]
}

fn lean_warm_color(u_raw: f64) -> [u8; 3] {
    let clamp01 = |x: f64| x.clamp(0.0, 1.0);
    // Mirror Lean kerr_ref_frame defaults so Rust can run in Lean-compatible mode.
    let exposure = 1.15_f64;
    let gamma = 1.55_f64;
    let black_level = 0.10_f64;
    let u = clamp01(u_raw).sqrt();
    let lifted = clamp01((u - black_level) / (1.0 - black_level).max(1e-9));
    let t = clamp01((exposure * lifted).powf(gamma)).sqrt();
    let r = clamp01(1.35 * t);
    let g = clamp01(0.92 * t * t + 0.08 * t);
    let b = clamp01(0.22 * t * t * t);
    [
        (255.0 * r).round() as u8,
        (255.0 * g).round() as u8,
        (255.0 * b).round() as u8,
    ]
}

#[allow(clippy::too_many_arguments)]
fn truth_color_from_trace(
    mode: TruthMode,
    result: TraceResult,
    lean_activity: Option<f64>,
    b_mag: f64,
    b_crit: f64,
    r_s: f64,
    r_isco: f64,
    disk_outer: f64,
    bx_raw: f64,
    sin_inc: f64,
    doppler: bool,
    plasma_model: PlasmaModel,
    tau_scale: f64,
    max_phi: f64,
    kerr_astar: f64,
) -> Option<[u8; 3]> {
    if mode == TruthMode::Off {
        return None;
    }

    let impact_norm = (b_mag / b_crit.max(1e-9)).clamp(0.0, 2.0) * 0.5;
    let angle_norm = |phi: f64| -> f64 { (phi / max_phi.max(1e-9)).clamp(0.0, 1.0) };
    let crossing_norm = |n: u32| -> f64 { (n as f64 / 6.0).clamp(0.0, 1.0) };

    let color = match (mode, result) {
        (TruthMode::Class, TraceResult::Captured) => [0, 0, 0],
        (TruthMode::Class, TraceResult::Escaped { .. }) => [50, 130, 255],
        (TruthMode::Class, TraceResult::DiskHit { .. }) => [255, 190, 40],

        (TruthMode::Impact, _) => false_color_ramp01(impact_norm),

        (TruthMode::Crossings, TraceResult::DiskHit { n_cross, .. }) => {
            false_color_ramp01(crossing_norm(n_cross))
        }
        (TruthMode::Crossings, _) => [0, 0, 0],

        (TruthMode::Angle, TraceResult::Escaped { phi_total }) => false_color_ramp01(angle_norm(phi_total)),
        (TruthMode::Angle, TraceResult::DiskHit { phi_orb, .. }) => false_color_ramp01(angle_norm(phi_orb.abs())),
        (TruthMode::Angle, TraceResult::Captured) => [0, 0, 0],

        (TruthMode::Lean, _) => lean_warm_color(lean_activity.unwrap_or(0.0)),

        (TruthMode::Transfer, TraceResult::DiskHit { r_eff, phi_orb, .. }) => {
            let g4 = disk_transfer_factor(r_eff, r_s, bx_raw, phi_orb, sin_inc, doppler, kerr_astar);
            let u = ((g4.max(1e-9).log10() + 6.0) / 8.5).clamp(0.0, 1.0);
            false_color_ramp01(u)
        }
        (TruthMode::Transfer, _) => [0, 0, 0],

        (TruthMode::Tau, TraceResult::DiskHit { r_eff, n_cross, .. }) => {
            let (_j_scale, a_scale) = plasma_profile_scales(r_eff, r_s, n_cross, plasma_model);
            let alpha_base = 0.35
                * tau_scale.max(0.0)
                * (1.0 + (n_cross.saturating_sub(1)) as f64 * 0.15)
                * a_scale;
            let path_scale = (r_eff / r_s.max(1e-9)).max(1e-9);
            let tau_eff = (alpha_base * path_scale).max(0.0);
            false_color_ramp01(1.0 - (-tau_eff).exp())
        }
        (TruthMode::Tau, _) => [0, 0, 0],
        (TruthMode::Off, _) => unreachable!("handled above"),
    };

    // Keep disk bounds in the transfer/tau channels honest.
    if matches!(mode, TruthMode::Transfer | TruthMode::Tau) {
        if let TraceResult::DiskHit { r_eff, .. } = result {
            if !(r_isco..=disk_outer).contains(&r_eff) {
                return Some([0, 0, 0]);
            }
        }
    }

    Some(color)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum DiskModel {
    Thin,
    Riaf,
}

impl DiskModel {
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    fn as_i32(self) -> i32 {
        match self {
            Self::Thin => 0,
            Self::Riaf => 1,
        }
    }

    fn from_env() -> Self {
        match std::env::var("BH_DISK_MODEL")
            .ok()
            .map(|s| s.to_ascii_lowercase())
            .as_deref()
        {
            Some("riaf") | Some("volumetric") => Self::Riaf,
            _ => Self::Thin,
        }
    }
    fn as_label(self) -> &'static str {
        match self {
            Self::Thin => "thin",
            Self::Riaf => "riaf",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PlasmaModel {
    Nt,
    Grmhd,
}

impl PlasmaModel {
    fn from_env() -> Self {
        match std::env::var("BH_PLASMA_MODEL")
            .ok()
            .map(|s| s.to_ascii_lowercase())
            .as_deref()
        {
            Some("grmhd") | Some("mad") | Some("sane") => Self::Grmhd,
            _ => Self::Nt,
        }
    }
    fn as_label(self) -> &'static str {
        match self {
            Self::Nt => "nt",
            Self::Grmhd => "grmhd",
        }
    }
}

// ── GPU FFI (available with --features cuda or --features rocm) ──────────────

#[cfg(any(feature = "cuda", feature = "rocm"))]
extern "C" {
    /// GPU black hole renderer — one CUDA thread per pixel.
    /// Writes RGB triples to `out_pixels` (host buffer, width*height*3 bytes).
    fn gutoe_render_bh(
        width: i32,
        height: i32,
        frame_width: i32,
        frame_height: i32,
        x_off: i32,
        y_off: i32,
        fov_rs: f64,
        inclination_deg: f64,
        r_s: f64,
        r_c: f64,
        disk_inner_rs: f64,
        disk_outer_rs: f64,
        max_phi: f64,
        dphi: f64,
        az_deg: f64,
        doppler: i32,
        ring_mode: i32,
        interior_mode: i32,
        core_look_mode: i32,
        spectral_band: i32, // 0..7, see SpectralBand
        disk_model: i32,    // 0=thin, 1=riaf
        plasma_model: i32,  // 0=nt, 1=grmhd profile proxy
        use_transfer: i32,  // 1 => use transfer_step path for disk intensity
        tau_scale: f64,     // optical-depth scale for transfer path
        fixed_exposure: f64, // <0 => per-band default exposure; >=0 force fixed exposure
        adaptive_dphi: i32, // 1 => reduce dphi near critical impact parameter
        kerr_enable: i32,   // 1 => use Kerr tracer for exterior rays
        kerr_astar: f64,    // dimensionless spin a*
        r_cam_rs: f64, // 0.0 = exterior; >0 = interior camera at r_cam_rs × r_s
        jitter_x: f64, // subpixel jitter in [0,1)
        jitter_y: f64, // subpixel jitter in [0,1)
        out_pixels: *mut u8,
    );
}

/// GPU render path — calls the CUDA/HIP kernel and returns the pixel buffer.
#[cfg(any(feature = "cuda", feature = "rocm"))]
fn render_with_options_gpu(
    metric: &GutoeMetric,
    disk_inner_rs: f64,
    disk_outer_rs: f64,
    cfg: &RenderConfig,
    az_deg: f64,
    doppler: bool,
    ring_mode: bool,
    interior_mode: bool,
    core_look_mode: bool,
    spectral_band: SpectralBand,
    disk_model: DiskModel,
    plasma_model: PlasmaModel,
    use_transfer: bool,
    tau_scale: f64,
    adaptive_dphi: bool,
    kerr: Option<&KerrMetric>,
    r_cam_rs: f64,
    jitter_x: f64,
    jitter_y: f64,
) -> Vec<[u8; 3]> {
    render_with_options_gpu_window(
        metric,
        disk_inner_rs,
        disk_outer_rs,
        cfg.width,
        cfg.height,
        0,
        0,
        cfg,
        az_deg,
        doppler,
        ring_mode,
        interior_mode,
        core_look_mode,
        spectral_band,
        disk_model,
        plasma_model,
        use_transfer,
        tau_scale,
        adaptive_dphi,
        kerr,
        r_cam_rs,
        jitter_x,
        jitter_y,
    )
}

#[cfg(any(feature = "cuda", feature = "rocm"))]
fn render_with_options_gpu_window(
    metric: &GutoeMetric,
    disk_inner_rs: f64,
    disk_outer_rs: f64,
    frame_width: usize,
    frame_height: usize,
    x_off: usize,
    y_off: usize,
    cfg: &RenderConfig,
    az_deg: f64,
    doppler: bool,
    ring_mode: bool,
    interior_mode: bool,
    core_look_mode: bool,
    spectral_band: SpectralBand,
    disk_model: DiskModel,
    plasma_model: PlasmaModel,
    use_transfer: bool,
    tau_scale: f64,
    adaptive_dphi: bool,
    kerr: Option<&KerrMetric>,
    r_cam_rs: f64,
    jitter_x: f64,
    jitter_y: f64,
) -> Vec<[u8; 3]> {
    let n = cfg.width * cfg.height;
    let mut flat = vec![0u8; n * 3];
    unsafe {
        gutoe_render_bh(
            cfg.width as i32,
            cfg.height as i32,
            frame_width as i32,
            frame_height as i32,
            x_off as i32,
            y_off as i32,
            cfg.fov_rs,
            cfg.inclination_deg,
            metric.r_s,
            metric.r_core(),
            disk_inner_rs,
            disk_outer_rs,
            cfg.max_phi,
            cfg.dphi,
            az_deg,
            doppler as i32,
            ring_mode as i32,
            interior_mode as i32,
            core_look_mode as i32,
            spectral_band.as_i32(),
            disk_model.as_i32(),
            match plasma_model {
                PlasmaModel::Nt => 0,
                PlasmaModel::Grmhd => 1,
            },
            use_transfer as i32,
            tau_scale,
            fixed_exposure_override().unwrap_or(-1.0),
            adaptive_dphi as i32,
            kerr.is_some() as i32,
            kerr.map_or(0.0, |k| k.a_star),
            r_cam_rs,
            jitter_x,
            jitter_y,
            flat.as_mut_ptr(),
        );
    }
    flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

// ── View definitions ─────────────────────────────────────────────────────────

struct View {
    /// Short name shown in page heading.
    label: &'static str,
    /// File slug: /tmp/bh_renders/{slug}.png
    slug: &'static str,
    /// Physics description shown under the image.
    caption: &'static str,
    /// Observer inclination from disk normal, degrees.  0 = face-on, 90 = edge-on.
    inc: f64,
    /// Disk azimuth rotation on screen, degrees.
    az: f64,
    /// Image half-width in r_s.
    fov: f64,
    /// Inner disk edge in r_s (areal radius).  3.0 = ISCO.
    disk_inner: f64,
    /// Outer disk edge in r_s (areal radius).
    disk_outer: f64,
    /// max_phi as a multiple of π.
    max_phi_pi: f64,
    /// RK4 step in radians.  0.005 good, 0.003 for sub-ring detail.
    dphi: f64,
    /// Apply relativistic Keplerian Doppler beaming D⁴.
    doppler: bool,
    /// Color pixels by photon ring order n instead of temperature.
    ring_mode: bool,
    /// Use pure Schwarzschild metric (l_P=0) instead of GUTOE.
    gr_mode: bool,
    /// False-colour the shadow interior by b/b_crit (impact-parameter map).
    interior_mode: bool,
    /// Interior camera aimed at core/lattice floor rather than outward sky.
    core_look_mode: bool,
    /// Literal interior camera at r_cam = r_cam_frac × r_horizon.
    /// 0.0 = exterior view (default).  0.5 = halfway between horizon and core.
    r_cam_frac: f64,
}

static VIEWS: &[View] = &[
    View {
        label: "Classic edge-on (85°)",
        slug: "v1_edge85",
        caption: "85° inclination — the classic view. The bright ring is the \
                  Novikov–Thorne accretion disk, T ∝ (r_ISCO/r)^{3/4}. \
                  The dark central region is the shadow: photons with b < b_crit = (3√3/2) r_s \
                  are captured. The glowing photon ring sits at the shadow boundary.",
        inc: 85.0,
        az: 0.0,
        fov: 12.0,
        disk_inner: 3.0,
        disk_outer: 10.0,
        max_phi_pi: 30.0,
        dphi: 0.005,
        doppler: false,
        ring_mode: false,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "Tilted 70°",
        slug: "v2_tilt70",
        caption: "70° inclination. The far side of the disk bends into view \
                  through gravitational lensing — the inner disk becomes visible \
                  above and below the shadow.",
        inc: 70.0,
        az: 0.0,
        fov: 12.0,
        disk_inner: 3.0,
        disk_outer: 10.0,
        max_phi_pi: 30.0,
        dphi: 0.005,
        doppler: false,
        ring_mode: false,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "Medium tilt (50°)",
        slug: "v3_tilt50",
        caption: "50° inclination. Strong lensing from the far side of the disk \
                  creates the characteristic asymmetric ring structure.",
        inc: 50.0,
        az: 0.0,
        fov: 12.0,
        disk_inner: 3.0,
        disk_outer: 10.0,
        max_phi_pi: 30.0,
        dphi: 0.005,
        doppler: false,
        ring_mode: false,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "Tilted toward face-on (30°)",
        slug: "v4_tilt30",
        caption: "30° inclination. The shadow becomes nearly circular and \
                  the disk face is visible as a bright annulus.",
        inc: 30.0,
        az: 0.0,
        fov: 12.0,
        disk_inner: 3.0,
        disk_outer: 10.0,
        max_phi_pi: 30.0,
        dphi: 0.005,
        doppler: false,
        ring_mode: false,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "Edge-on rotated 60° (85°, az=60°)",
        slug: "v5_edge_rot",
        caption: "85° inclination, disk rotated 60° in the image plane. \
                  Confirms the axial symmetry of the GUTOE Schwarzschild metric.",
        inc: 85.0,
        az: 60.0,
        fov: 12.0,
        disk_inner: 3.0,
        disk_outer: 10.0,
        max_phi_pi: 30.0,
        dphi: 0.005,
        doppler: false,
        ring_mode: false,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "Nearly face-on (10°, circular shadow)",
        slug: "v6_face10",
        caption: "10° inclination — almost face-on. The shadow is nearly circular \
                  (b_crit is the same in every azimuthal direction for Schwarzschild). \
                  The disk appears as a wide bright ring.",
        inc: 10.0,
        az: 0.0,
        fov: 12.0,
        disk_inner: 3.0,
        disk_outer: 10.0,
        max_phi_pi: 30.0,
        dphi: 0.005,
        doppler: false,
        ring_mode: false,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "M87★ — 17° inclination (Doppler crescent)",
        slug: "m87star",
        caption: "M87★ geometry: jet axis 17° from the line of sight → disk at 17° \
                  from face-on. Relativistic Doppler beaming (D⁴ factor, prograde orbit) \
                  creates a crescent asymmetry: the approaching side (right) is boosted \
                  by up to 57% at the ISCO (v_K ≈ 0.41c, β_obs ≈ 0.12 at 17°). \
                  Compare with EHT 2019/2024 images of M87★.",
        inc: 17.0,
        az: 0.0,
        fov: 7.0,
        disk_inner: 3.0,
        disk_outer: 16.0,
        max_phi_pi: 40.0,
        dphi: 0.005,
        doppler: true,
        ring_mode: false,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "M87★ EHT-2017 aligned (230 GHz, 17°, pullback)",
        slug: "m87_eht2017",
        caption: "EHT-style observer geometry for M87★: 230 GHz band, 17° inclination, \
                  wider framing for ring-scale comparison and synthetic-observation workflows.",
        inc: 17.0,
        az: 0.0,
        fov: 8.5,
        disk_inner: 3.0,
        disk_outer: 16.0,
        max_phi_pi: 40.0,
        dphi: 0.003,
        doppler: true,
        ring_mode: false,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "Sgr A★ — 50° inclination (strong Doppler)",
        slug: "sgr_astar",
        caption: "Sgr A★ geometry: ~50° inclination. At this angle Doppler \
                  beaming is much stronger (sin 50° ≈ 0.77, β_obs ≈ 0.31 at ISCO) \
                  giving a 4.5× brightness ratio between approaching and receding disk. \
                  The EHT 2022 image of Sgr A★ shows a similar bright crescent.",
        inc: 50.0,
        az: 0.0,
        fov: 7.0,
        disk_inner: 3.0,
        disk_outer: 12.0,
        max_phi_pi: 30.0,
        dphi: 0.005,
        doppler: true,
        ring_mode: false,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "Photon sub-rings n=1,2,3 (85°, zoomed)",
        slug: "rings85",
        caption: "Photon ring structure at 85°, zoomed to ±6 r_s, dphi=0.003. \
                  Direct image (n=1, orange): photon hits disk on first inward pass. \
                  Secondary ring (n=2, cyan): photon orbited the BH once before \
                  hitting the disk — sits just outside the shadow boundary. \
                  Tertiary ring (n=3, purple): orbited twice. Each sub-ring n carries \
                  an independent image of the entire accretion disk, exponentially \
                  demagnified. Future ngEHT observations may resolve n=2.",
        inc: 85.0,
        az: 0.0,
        fov: 6.0,
        disk_inner: 3.0,
        disk_outer: 10.0,
        max_phi_pi: 50.0,
        dphi: 0.003,
        doppler: false,
        ring_mode: true,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "GR comparison — pure Schwarzschild (85°)",
        slug: "gr_compare",
        caption: "Pure Schwarzschild metric (l_P = 0, classical GR). \
                  For macroscopic black holes (M87★: r_s/l_P ≈ 10⁴⁴), \
                  the GUTOE lattice correction (r_core/r_s)² ≈ (l_P/r_s)² is \
                  unmeasurably small. In Planck units (r_s = l_P = 1) the \
                  correction is r_core = √C_∞ ≈ 0.74 — visible only near the \
                  shadow boundary where geodesics orbit for many revolutions.",
        inc: 85.0,
        az: 0.0,
        fov: 12.0,
        disk_inner: 3.0,
        disk_outer: 10.0,
        max_phi_pi: 30.0,
        dphi: 0.005,
        doppler: false,
        ring_mode: false,
        gr_mode: true,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "Shadow interior — fractal impact-parameter map",
        slug: "shadow_interior",
        caption: "False-colour map of the shadow interior (outside-observer view). \
                  Captured photons are coloured by half-orbit count n ≈ −ln(1 − b/b_crit)/π \
                  before falling in. Each colour band is a complete image of the disk, \
                  exponentially demagnified toward the photon ring — the fractal self-similar \
                  structure of the black hole shadow. Orange = n=1, cyan = n=2, violet = n=3.",
        inc: 85.0,
        az: 0.0,
        fov: 12.0,
        disk_inner: 3.0,
        disk_outer: 10.0,
        max_phi_pi: 50.0,
        dphi: 0.003,
        doppler: false,
        ring_mode: false,
        gr_mode: false,
        interior_mode: true,
        core_look_mode: false,
        r_cam_frac: 0.0,
    },
    View {
        label: "Inside the horizon — camera at r = 0.5 r_h",
        slug: "camera_inside",
        caption: "Camera placed at coordinate radius r_cam = 0.5 × r_horizon — deep inside \
                  the GUTOE event horizon.  Photons are fired outward in all directions. \
                  Centre disc (b < b_crit): the entire outside universe compressed overhead — \
                  stars, the accretion disk, and the photon ring form a bright window to the \
                  cosmos above.  Outer ring (b > b_crit): photons turn around before reaching \
                  the photon sphere and fall back — they illuminate the GUTOE lattice floor at \
                  r_core = √C_∞ l_P glowing amber-orange.  The ring between disc and floor is \
                  the photon ring seen from below — identical critical impact parameter as from \
                  outside. This image has never before been computed.",
        inc: 85.0,
        az: 0.0,
        fov: 5.0,
        disk_inner: 3.0,
        disk_outer: 16.0,
        max_phi_pi: 50.0,
        dphi: 0.003,
        doppler: true,
        ring_mode: false,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: false,
        r_cam_frac: 0.5,
    },
    View {
        label: "Inside the horizon — looking down at the lattice core",
        slug: "camera_core",
        caption: "Interior camera aimed toward the regularized GUTOE core instead of outward. \
                  Rays are traced on the plunging branch from r_cam = 0.72 × r_horizon toward \
                  r_core = √C_∞ l_P. Colour is physics-only from traced invariants \
                  (impact parameter, orbital angle, capture radius), with no procedural texture.",
        inc: 85.0,
        az: 0.0,
        fov: 3.6,
        disk_inner: 3.0,
        disk_outer: 16.0,
        max_phi_pi: 36.0,
        dphi: 0.003,
        doppler: false,
        ring_mode: false,
        gr_mode: false,
        interior_mode: false,
        core_look_mode: true,
        r_cam_frac: 0.72,
    },
];

// ── Colour model ─────────────────────────────────────────────────────────────

/// Returns the RGB colour for a disk pixel.
///
/// # Arguments
/// - `r_eff`     — areal radius of disk crossing (same units as r_s)
/// - `r_isco`    — ISCO areal radius (= 3 r_s)
/// - `r_outer`   — soft outer disk edge (emission tapers exponentially beyond this)
/// - `r_s`       — Schwarzschild radius
/// - `bx`        — horizontal impact parameter (signed; positive = approaching side)
/// - `sin_inc`   — sin(observer inclination from disk normal)
/// - `n_cross`   — photon ring order (1 = direct, 2 = secondary, …)
/// - `doppler`   — if true, apply relativistic Doppler D⁴ factor
/// - `ring_mode` — if true, colour by ring order instead of temperature
fn pixel_color(
    r_eff: f64,
    r_isco: f64,
    r_outer: f64,
    r_s: f64,
    bx: f64,
    phi_orb: f64,
    sin_inc: f64,
    n_cross: u32,
    doppler: bool,
    ring_mode: bool,
    spectral_band: SpectralBand,
    plasma_model: PlasmaModel,
    use_transfer: bool,
    tau_scale: f64,
    kerr_astar: f64,
) -> [u8; 3] {
    if ring_mode {
        return ring_order_color(n_cross);
    }

    // Novikov–Thorne temperature profile: T ∝ (r_ISCO / r_eff)^{3/4}
    let t_rel = (r_isco / r_eff).powf(0.75);

    // Smooth outer disk taper: emission falls as exp(-((r-r_outer)/(0.5·r_outer))²)
    // Avoids the unnatural hard edge at disk_outer.  Inside the nominal outer
    // radius the factor is 1.0; it decays gently beyond it.
    let outer_taper = {
        let excess = (r_eff - r_outer).max(0.0) / (0.5 * r_outer);
        (-excess * excess).exp()
    };

    // Higher-order images are dimmer (photon lost energy looping around)
    let fade = 0.65_f64.powi(n_cross as i32 - 1);

    // Relativistic transfer factor g⁴ (gravitational redshift × Doppler beaming).
    let transfer = disk_transfer_factor(r_eff, r_s, bx, phi_orb, sin_inc, doppler, kerr_astar);

    let spectral = band_weight_with_exposure(spectral_band, t_rel, fixed_exposure_override());
    let (j_scale, a_scale) = plasma_profile_scales(r_eff, r_s, n_cross, plasma_model);
    // Local covariant source proxy.
    let source_local = (t_rel * fade * outer_taper * spectral * j_scale).max(0.0);
    let g_cov = transfer.max(1e-9).powf(0.25);
    let alpha_base = 0.35
        * tau_scale.max(0.0)
        * (1.0 + (n_cross.saturating_sub(1)) as f64 * 0.15)
        * a_scale;
    let luminance_raw = if use_transfer {
        // Multi-step covariant transfer integration along an effective path
        // segment through the emitting flow.
        let steps = 8usize;
        let path_scale = (r_eff / r_s.max(1e-9)).max(1e-9);
        let mut intensity = 0.0_f64;
        for si in 0..steps {
            let u = (si as f64 + 0.5) / steps as f64; // 0..1
            let local_mod = 1.0 + 0.20 * (1.0 - u);
            let j_obs = covariant_emissivity((source_local * local_mod).max(0.0), g_cov);
            let alpha_obs = covariant_absorption(alpha_base * (0.7 + 0.6 * u), g_cov);
            let tau_seg = (alpha_obs * path_scale / steps as f64).max(0.0);
            let source_fn = if alpha_obs > 1e-12 {
                j_obs / alpha_obs
            } else {
                j_obs
            };
            intensity = transfer_step(intensity, source_fn, tau_seg);
        }
        intensity.max(0.0)
    } else {
        source_local * transfer
    };
    // Reinhard tone mapping: luminance → luminance / (1 + luminance).
    let b = luminance_raw.max(0.0) / (1.0 + luminance_raw.max(0.0));

    // Orange-white thermal palette (hot inner disk white, outer disk orange-red)
    let [tr, tg, tb] = band_tint(spectral_band);
    let r = (255.0 * b.powf(0.35) * tr).clamp(0.0, 255.0) as u8;
    let g = (210.0 * b.powf(0.60) * tg).clamp(0.0, 255.0) as u8;
    let bl = (130.0 * b.powf(1.60) * tb).clamp(0.0, 255.0) as u8;
    [r, g, bl]
}

fn plasma_profile_scales(
    r_eff: f64,
    r_s: f64,
    n_cross: u32,
    plasma_model: PlasmaModel,
) -> (f64, f64) {
    fn parse_env_f64(key: &str, default: f64) -> f64 {
        std::env::var(key)
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(default)
    }

    // EHT-standard electron heating proxy:
    // Ti/Te = (R_high * beta^2 + R_low) / (1 + beta^2)
    // with beta = P_gas / P_mag.  We use local profile proxies for P_gas/P_mag.
    fn te_over_ti_rbeta(beta: f64, r_low: f64, r_high: f64) -> f64 {
        let b2 = beta * beta;
        let ti_over_te = (r_high * b2 + r_low) / (1.0 + b2).max(1e-9);
        (1.0 / ti_over_te).clamp(1e-3, 1.0)
    }

    match plasma_model {
        PlasmaModel::Nt => (1.0, 1.0),
        PlasmaModel::Grmhd => {
            let x = (r_eff / r_s.max(1e-9)).max(1e-9);
            // Simple GRMHD-inspired profile proxy:
            // density ~ r^-1.1, temperature ~ r^-0.8, magnetic energy ~ r^-1.0
            // emissivity boost ∝ n_e * B * T_e^0.5, absorption tracks n_e * B / T_e.
            let ne = x.powf(-1.1);
            let te = x.powf(-0.8);
            let b = x.powf(-1.0);
            let r_low = parse_env_f64("BH_R_LOW", 1.0).max(1e-3);
            let r_high = parse_env_f64("BH_R_HIGH", 40.0).max(r_low);
            let p_gas = ne * te;
            let p_mag = (b * b).max(1e-9);
            let beta = (p_gas / p_mag).max(1e-6);
            let te_ti = te_over_ti_rbeta(beta, r_low, r_high);
            let electron_boost = te_ti.sqrt();
            let ring = 1.0 + n_cross.saturating_sub(1) as f64 * 0.08;
            let j_scale = (ne * b * te.sqrt() * electron_boost * ring).clamp(0.08, 6.0);
            let a_scale = (ne * b / (te * te_ti).max(1e-6) * ring).clamp(0.05, 8.0);
            (j_scale, a_scale)
        }
    }
}

fn riaf_composite_color(
    r_eff: f64,
    r_isco: f64,
    r_outer: f64,
    r_s: f64,
    bx_raw: f64,
    bx: f64,
    by: f64,
    sin_inc: f64,
    n_cross: u32,
    phi_orb: f64,
    doppler: bool,
    ring_mode: bool,
    spectral_band: SpectralBand,
    plasma_model: PlasmaModel,
    tau_scale: f64,
    kerr_astar: f64,
) -> [u8; 3] {
    let disk = pixel_color(
        r_eff,
        r_isco,
        r_outer,
        r_s,
        bx_raw,
        phi_orb,
        sin_inc,
        n_cross,
        doppler,
        ring_mode,
        spectral_band,
        plasma_model,
        true,
        tau_scale,
        kerr_astar,
    );
    let bg = star_field_color(bx, by, phi_orb + std::f64::consts::PI);
    // Low optical-depth RIAF proxy: hot diffuse flow, not an opaque wall.
    let tau = (0.45 * tau_scale.max(0.0))
        * (r_s / r_eff.max(1e-9)).powf(0.7)
        * (1.0 + n_cross.saturating_sub(1) as f64 * 0.10);
    let trans = (-tau).exp().clamp(0.0, 1.0);
    let gain = 1.6;
    [
        (((disk[0] as f64) * gain).min(255.0) * (1.0 - trans) + (bg[0] as f64) * trans)
            .round()
            .clamp(0.0, 255.0) as u8,
        (((disk[1] as f64) * gain).min(255.0) * (1.0 - trans) + (bg[1] as f64) * trans)
            .round()
            .clamp(0.0, 255.0) as u8,
        (((disk[2] as f64) * gain).min(255.0) * (1.0 - trans) + (bg[2] as f64) * trans)
            .round()
            .clamp(0.0, 255.0) as u8,
    ]
}

/// Disk radiative transfer proxy `g^4` where:
/// - `g_gr = sqrt(1 - r_s/r)` (static gravitational redshift),
/// - `g_dop = 1/(gamma * (1 - beta_obs))` (SR beaming factor),
/// - total transfer is `(g_gr * g_dop)^4`.
///
/// This is still a Schwarzschild disk approximation, but closer to physical
/// transfer than pure `D^4` and keeps GPU/CPU parity.
fn disk_transfer_factor(
    r_eff: f64,
    r_s: f64,
    bx: f64,
    phi_orb: f64,
    sin_inc: f64,
    doppler: bool,
    kerr_astar: f64,
) -> f64 {
    let r_safe = r_eff.max(1e-12);
    let g_gr = (1.0 - r_s / r_safe).max(0.0).sqrt();
    // Kerr frame-drag proxy in emissive flow (strong near horizon, fades outward).
    let mut beta = (r_s / (2.0 * r_safe)).sqrt().min(0.7);
    let mut omega_boost = 0.0_f64;
    if kerr_astar.abs() > 1e-9 {
        omega_boost = 0.35 * kerr_astar.abs() * (r_s / r_safe).clamp(0.0, 1.0).powf(1.5);
        beta = (beta * (1.0 + omega_boost)).min(0.85);
    }
    if !doppler {
        let g4 = g_gr.powi(4).clamp(1e-6, 300.0);
        return g4 * (1.0 + 0.15 * omega_boost);
    }

    // Keplerian orbital speed proxy; capped for numerical robustness.
    let gamma = 1.0 / (1.0 - beta * beta).sqrt();
    // More physical LOS projection: dominant azimuthal emitter angle, with a
    // small screen-space blend to preserve continuity near direct hits.
    let mu_phi = phi_orb.sin().clamp(-1.0, 1.0);
    let mu_screen = (bx / r_safe).clamp(-1.0, 1.0);
    let mu = (0.8 * mu_phi + 0.2 * mu_screen).clamp(-1.0, 1.0);
    let beta_obs = beta * sin_inc * mu;
    let g_dop = 1.0 / (gamma * (1.0 - beta_obs));
    let g = g_gr * g_dop;
    let g4 = g.powi(4).clamp(1e-6, 300.0);
    g4 * (1.0 + 0.15 * omega_boost)
}

/// False-colour the shadow interior by capture depth.
///
/// Photons with b < b_crit are captured. The number of approximate half-orbits
/// before capture: n ≈ −ln(1 − b/b_crit) / π  (logarithmic divergence near b_crit).
///
/// Each band (integer n) is a complete image of the entire disk, exponentially
/// demagnified. The bands nest infinitely toward the photon ring — the fractal
/// self-similar structure of the black hole shadow.
///
/// Colour cycling (4 bands per cycle): orange → cyan → purple → warm-white → repeat.
/// Brightness increases with n (captures with more orbits are highlighted).
fn shadow_interior_color(bx: f64, by: f64, r_s: f64) -> [u8; 3] {
    use std::f64::consts::PI;
    let b = (bx * bx + by * by).sqrt();
    let b_crit = 1.5 * 3.0_f64.sqrt() * r_s; // (3√3/2) r_s
    let ratio = (b / b_crit).clamp(0.0, 1.0 - 1e-9);
    // Approximate half-orbit count before capture (exact at photon sphere limit)
    let n_float = -(1.0 - ratio).ln() / PI;
    let n = n_float as u32; // floor: which band
    let frac = n_float.fract() as f32; // position within band [0..1)

    // Global brightness: peaks at n=1 (the most visible inner ring), fades toward center
    let brightness = ((n as f32 + 1.0) * 0.5).tanh() * (0.3 + 0.7 * frac);

    // Band colour cycles every 4 half-orbits (orange → cyan → purple → warm-white)
    let (rr, gg, bb) = match n % 4 {
        0 => (1.0_f32, 0.55, 0.10), // deep orange — direct inner ring
        1 => (0.15, 0.90, 0.80),    // cyan / teal — secondary inner ring
        2 => (0.75, 0.20, 1.00),    // violet — tertiary
        _ => (1.00, 0.95, 0.60),    // warm yellow-white — quaternary+
    };
    [
        (rr * brightness * 255.0).clamp(0.0, 255.0) as u8,
        (gg * brightness * 255.0).clamp(0.0, 255.0) as u8,
        (bb * brightness * 255.0).clamp(0.0, 255.0) as u8,
    ]
}

/// GUTOE lattice-floor glow for the interior-camera view (b > b_crit photons).
///
/// These are outward-fired photons from inside the horizon that turned around
/// before reaching the photon sphere.  Photons with b ≈ b_crit orbited many
/// times near r_ph before falling back → they are brightest (have gathered the
/// most energy from the photon ring).  Photons with b >> b_crit turned quickly
/// → dimmer.  Palette: hot amber-orange → dark red → black at large excess.
fn gutoe_core_color(b: f64, b_crit: f64) -> [u8; 3] {
    let excess = (b - b_crit).max(0.0) / b_crit;
    let glow = (-excess * excess * 3.0).exp();
    // Reinhard: peak L = 4.0 → bv ≈ 0.8 (near-white hot at the photon ring)
    let l = glow * 4.0;
    let bv = l / (1.0 + l);
    [
        (255.0 * bv.powf(0.35)).clamp(0.0, 255.0) as u8,
        (160.0 * bv.powf(0.65)).clamp(0.0, 255.0) as u8,
        (30.0 * bv.powf(2.00)).clamp(0.0, 255.0) as u8,
    ]
}

/// Core-facing interior palette derived only from traced geodesic invariants.
fn gutoe_core_physics_color(
    b: f64,
    b_crit: f64,
    r_eff_hit: f64,
    phi_orb: f64,
    r_cam: f64,
    r_core: f64,
) -> [u8; 3] {
    let eta = (b / b_crit.max(1e-9)).clamp(0.0, 2.0); // impact class
    let n_half = (phi_orb / std::f64::consts::PI).max(0.0); // winding count proxy
    let re_span = (r_cam - r_core).max(1e-9);
    let depth = ((r_eff_hit - r_core) / re_span).clamp(0.0, 1.0); // 0 core, 1 near camera

    // Near-critical impact parameters linger longer before capture.
    let near_crit = (-(eta - 1.0).powi(2) / 0.08).exp();
    let winding = (n_half / 12.0).min(1.0);
    let plunge = 1.0 - depth;
    let luminance = (0.10 + 0.75 * near_crit + 0.45 * winding + 0.35 * plunge).max(0.0);
    let tone = luminance / (1.0 + luminance);

    [
        (255.0 * tone.powf(0.36)).clamp(0.0, 255.0) as u8,
        (170.0 * tone.powf(0.62)).clamp(0.0, 255.0) as u8,
        (45.0 * tone.powf(1.50)).clamp(0.0, 255.0) as u8,
    ]
}

/// Colour-by-photon-ring-order: n=1 orange, n=2 cyan, n=3 purple, n≥4 grey.
fn ring_order_color(n_cross: u32) -> [u8; 3] {
    match n_cross {
        1 => [255, 165, 40],  // direct image — warm orange
        2 => [40, 220, 200],  // secondary ring — cyan / teal
        3 => [190, 60, 255],  // tertiary ring — purple
        _ => [100, 100, 100], // higher order — grey
    }
}

// ── Procedural star field ─────────────────────────────────────────────────────

/// splitmix64-style spatial hash of two integer sky-grid coordinates.
#[inline]
fn star_hash(x: i64, y: i64) -> u64 {
    let s = (x as u64)
        .wrapping_mul(0x9e3779b97f4a7c15)
        .wrapping_add((y as u64).wrapping_mul(0x6c62272e07bb0142));
    let mut h = s ^ (s >> 30);
    h = h.wrapping_mul(0xbf58476d1ce4e5b9);
    h ^= h >> 27;
    h = h.wrapping_mul(0x94d049bb133111eb);
    h ^= h >> 31;
    h
}

#[inline]
fn parity_legacy_stars_enabled() -> bool {
    static LEGACY: OnceLock<bool> = OnceLock::new();
    *LEGACY.get_or_init(|| {
        std::env::var("BH_PARITY_LEGACY_STARS")
            .ok()
            .is_some_and(|v| v == "1" || v.eq_ignore_ascii_case("true"))
    })
}

#[inline]
fn fixed_exposure_override() -> Option<f64> {
    static EXP: OnceLock<Option<f64>> = OnceLock::new();
    *EXP.get_or_init(|| {
        std::env::var("BH_FIXED_EXPOSURE")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .map(|v| v.max(0.0))
    })
}

/// Returns an RGB pixel for an escaped (background) ray.
///
/// Uses the **sky direction** (computed by rotating the impact parameter vector
/// by the total orbital angle phi_total) as the hash key.  This ensures that:
///
/// - Background stars sit at fixed sky positions (consistent across the image)
/// - Gravitational lensing is visible: near the photon ring, phi_total ≈ 2π, 4π, …
///   → the same stars appear in multiple Einstein-ring copies at different radii
/// - Near the shadow boundary phi_total spins rapidly → photon ring shimmers with
///   rapidly changing star backgrounds (characteristic lensing noise)
///
/// Star density ≈ 1.5 % of sky-grid cells.  Spectral distribution:
///   deep-red M-star · orange K · warm-white F/G · pure-white A · blue-white B/O
fn star_field_color(bx: f64, by: f64, phi_total: f64) -> [u8; 3] {
    // Diagnostics mode: match CUDA legacy stars exactly so parity checks isolate
    // transfer/tone-map differences rather than background generator divergence.
    if parity_legacy_stars_enabled() {
        let sky_x = bx * phi_total.cos() - by * phi_total.sin();
        let sky_y = bx * phi_total.sin() + by * phi_total.cos();
        let hx = (sky_x * 50.0).round() as i64;
        let hy = (sky_y * 50.0).round() as i64;
        let h = star_hash(hx, hy);
        if (h & 0xFFFF) >= 1000 {
            let v = ((h >> 40) & 7) as u8;
            return [v >> 2, v >> 2, (v >> 1).saturating_add(8)];
        }
        let bright_raw = ((h >> 16) & 0xFF) as u8;
        let bright = 80u8.saturating_add(((bright_raw as u16 * 175) / 255) as u8);
        let spec = ((h >> 24) & 0xF) as u8;
        return if spec <= 1 {
            [bright / 2, bright / 5, bright / 10]
        } else if spec <= 4 {
            [
                bright,
                ((bright as u16 * 68) / 100) as u8,
                ((bright as u16 * 38) / 100) as u8,
            ]
        } else if spec <= 9 {
            [bright, bright, ((bright as u16 * 88) / 100) as u8]
        } else if spec <= 12 {
            [bright, bright, bright]
        } else {
            let rg = ((bright as u16 * 82) / 100) as u8;
            [rg, rg, bright]
        };
    }

    // Rotate impact-parameter vector (bx, by) by phi_total to get sky direction.
    // In the Schwarzschild geometry a photon starting from angle 0 exits at
    // phi_total → source is at angle phi_total from the observer direction.
    // phi_total ≈ π : direct (weakly deflected) ray — sky ≈ "behind" BH
    // phi_total ≈ 2π: once-orbiting ray   — Einstein ring; same stars as direct
    // phi_total ≫ π : photon ring regime  — rapidly changing star backgrounds
    let sky_x = bx * phi_total.cos() - by * phi_total.sin();
    let sky_y = bx * phi_total.sin() + by * phi_total.cos();
    let map_rgb = starmap().map(|m| {
        // Map asymptotic ray direction to equirectangular (lon/lat) UV.
        let f = 1.0_f64;
        let inv = 1.0 / (sky_x * sky_x + sky_y * sky_y + f * f).sqrt();
        let dx = sky_x * inv;
        let dy = sky_y * inv;
        let dz = f * inv;
        let lon = dx.atan2(dz); // [-pi,pi]
        let lat = dy.asin(); // [-pi/2,pi/2]
        let u = lon / (2.0 * std::f64::consts::PI) + 0.5;
        let v = 0.5 - lat / std::f64::consts::PI;
        m.sample(u, v)
    });

    let mut r = 0.0015_f64;
    let mut g = 0.0020_f64;
    let mut b = 0.0060_f64;

    // Smooth Milky-Way-style band (no blocky per-cell "clouds").
    let gx = sky_x * 0.035;
    let gy = sky_y * 0.035;
    let lat = (gy + 0.28 * (0.7 * gx).sin() + 0.14 * (1.3 * gx + 0.5).sin()).abs();
    let band = (-(lat * lat) / 0.05).exp();
    // Continuous pseudo-noise dust modulation.
    let dust = ((sky_x * 0.11).sin() * (sky_y * 0.07 + 1.7).cos() * 0.5 + 0.5).powf(1.4);
    let gal = band * (0.25 + 0.75 * dust);
    r += 0.16 * gal;
    g += 0.14 * gal;
    b += 0.20 * gal;

    // Sparse point stars only (avoid unresolved noise blanket).
    let layers = [
        (42.0_f64, 0.008_f64, 0.55_f64),
        (86.0_f64, 0.020_f64, 0.20_f64),
        (150.0_f64, 0.045_f64, 0.08_f64),
    ];
    for (scale, density, amp) in layers {
        let hx = (sky_x * scale).round() as i64;
        let hy = (sky_y * scale).round() as i64;
        let h = star_hash(hx, hy);
        let gate = ((h & 0xFFFF) as f64) / 65535.0;
        if gate > density {
            continue;
        }
        let bright_u = ((h >> 16) & 0xFF) as f64 / 255.0;
        let bright = amp * (0.25 + 0.75 * bright_u.powf(1.8));
        let (sr, sg, sb) = match (h >> 24) & 0xF {
            0..=2 => (1.00, 0.62, 0.35),
            3..=6 => (1.00, 0.92, 0.78),
            7..=11 => (1.00, 1.00, 1.00),
            _ => (0.78, 0.86, 1.00),
        };
        let boost = 1.0 + 0.6 * band;
        r += bright * sr * boost;
        g += bright * sg * boost;
        b += bright * sb * boost;
    }

    // Tone map + tiny dither to suppress visible banding.
    let rr = (r / (1.0 + r)).powf(1.0 / 2.0);
    let gg = (g / (1.0 + g)).powf(1.0 / 2.0);
    let bb = (b / (1.0 + b)).powf(1.0 / 2.0);
    let dh = star_hash((sky_x * 121.0) as i64, (sky_y * 157.0) as i64);
    let d = (((dh >> 8) & 0xFF) as f64 / 255.0 - 0.5) * (0.7 / 255.0);
    let proc = [
        ((rr + d) * 255.0).clamp(0.0, 255.0) as u8,
        ((gg + d) * 255.0).clamp(0.0, 255.0) as u8,
        ((bb + d) * 255.0).clamp(0.0, 255.0) as u8,
    ];
    if let Some(map) = map_rgb {
        // Keep a bit of procedural sparkle over real-sky texture.
        return [
            ((map[0] as f64 * 0.88) + (proc[0] as f64 * 0.12)).round().clamp(0.0, 255.0) as u8,
            ((map[1] as f64 * 0.88) + (proc[1] as f64 * 0.12)).round().clamp(0.0, 255.0) as u8,
            ((map[2] as f64 * 0.88) + (proc[2] as f64 * 0.12)).round().clamp(0.0, 255.0) as u8,
        ];
    }
    proc
}

// ── Gaussian PSF blur ─────────────────────────────────────────────────────────

/// Separable 2D Gaussian blur.  sigma < 0.5 → returns a copy unchanged.
///
/// Used to simulate the EHT telescope beam (FWHM ≈ 20 μas ≈ 40% shadow size).
/// σ = 6  → light anti-aliasing glow at 4K.
/// σ = 30 → realistic EHT M87★ beam smear (~40% shadow, fuzzy donut).
///
/// Runtime: O(W × H × 6σ) per channel.  Two passes (H then V).
fn gaussian_blur(pixels: &[[u8; 3]], w: usize, h: usize, sigma: f64) -> Vec<[u8; 3]> {
    if sigma < 0.5 {
        return pixels.to_vec();
    }
    let radius = (3.0 * sigma).ceil() as usize;
    let ksize = 2 * radius + 1;

    // Build normalised 1D kernel
    let kernel: Vec<f64> = (0..ksize)
        .map(|i| {
            let x = i as f64 - radius as f64;
            (-x * x / (2.0 * sigma * sigma)).exp()
        })
        .collect();
    let ksum: f64 = kernel.iter().sum();
    let kernel: Vec<f64> = kernel.iter().map(|&k| k / ksum).collect();

    // Horizontal pass: pixels → temp (f64 accumulator)
    let mut temp = vec![[0f64; 3]; w * h];
    for y in 0..h {
        for x in 0..w {
            let mut acc = [0f64; 3];
            for (ki, &kw) in kernel.iter().enumerate() {
                let cx = {
                    let raw = x + ki;
                    if raw < radius {
                        0
                    } else {
                        (raw - radius).min(w - 1)
                    }
                };
                let p = pixels[y * w + cx];
                acc[0] += kw * p[0] as f64;
                acc[1] += kw * p[1] as f64;
                acc[2] += kw * p[2] as f64;
            }
            temp[y * w + x] = acc;
        }
    }

    // Vertical pass: temp → out (u8)
    let mut out = vec![[0u8; 3]; w * h];
    for y in 0..h {
        for x in 0..w {
            let mut acc = [0f64; 3];
            for (ki, &kw) in kernel.iter().enumerate() {
                let cy = {
                    let raw = y + ki;
                    if raw < radius {
                        0
                    } else {
                        (raw - radius).min(h - 1)
                    }
                };
                let p = temp[cy * w + x];
                acc[0] += kw * p[0];
                acc[1] += kw * p[1];
                acc[2] += kw * p[2];
            }
            out[y * w + x] = [
                acc[0].clamp(0.0, 255.0) as u8,
                acc[1].clamp(0.0, 255.0) as u8,
                acc[2].clamp(0.0, 255.0) as u8,
            ];
        }
    }
    out
}

/// Simple unsharp mask to recover ring detail after beam blur.
#[allow(dead_code)]
fn unsharp_mask(
    pixels: &[[u8; 3]],
    w: usize,
    h: usize,
    blur_sigma: f64,
    amount: f64,
) -> Vec<[u8; 3]> {
    if blur_sigma < 0.5 || amount <= 0.0 {
        return pixels.to_vec();
    }
    let low = gaussian_blur(pixels, w, h, blur_sigma);
    let mut out = vec![[0u8; 3]; w * h];
    for i in 0..(w * h) {
        let mut p = [0u8; 3];
        for c in 0..3 {
            let hi = pixels[i][c] as f64;
            let lo = low[i][c] as f64;
            let v = (hi + amount * (hi - lo)).clamp(0.0, 255.0);
            p[c] = v as u8;
        }
        out[i] = p;
    }
    out
}

fn save_png_rgb(path: &Path, pixels: &[[u8; 3]], w: usize, h: usize) {
    let mut img = image::RgbImage::new(w as u32, h as u32);
    for (i, p) in pixels.iter().enumerate() {
        let x = (i % w) as u32;
        let y = (i / w) as u32;
        img.put_pixel(x, y, image::Rgb([p[0], p[1], p[2]]));
    }
    img.save(path).expect("save png");
}

fn save_png_gray_u8(path: &Path, pixels: &[u8], w: usize, h: usize) {
    let mut img = image::GrayImage::new(w as u32, h as u32);
    for (i, p) in pixels.iter().enumerate() {
        let x = (i % w) as u32;
        let y = (i / w) as u32;
        img.put_pixel(x, y, image::Luma([*p]));
    }
    img.save(path).expect("save gray png");
}

fn save_png_gray_f64_norm(path: &Path, vals: &[f64], w: usize, h: usize) {
    if vals.is_empty() {
        return;
    }
    let mut vmin = f64::INFINITY;
    let mut vmax = f64::NEG_INFINITY;
    for &v in vals {
        if v.is_finite() {
            vmin = vmin.min(v);
            vmax = vmax.max(v);
        }
    }
    if !vmin.is_finite() || !vmax.is_finite() {
        vmin = 0.0;
        vmax = 1.0;
    }
    let scale = if (vmax - vmin).abs() < 1e-12 {
        0.0
    } else {
        255.0 / (vmax - vmin)
    };
    let mut img = image::GrayImage::new(w as u32, h as u32);
    for (i, &v) in vals.iter().enumerate() {
        let u = if scale == 0.0 || !v.is_finite() {
            0u8
        } else {
            ((v - vmin) * scale).round().clamp(0.0, 255.0) as u8
        };
        let x = (i % w) as u32;
        let y = (i / w) as u32;
        img.put_pixel(x, y, image::Luma([u]));
    }
    img.save(path).expect("save norm gray png");
}

fn moffat_kernel(alpha: f64, beta: f64, radius: usize) -> Vec<f64> {
    let n = 2 * radius + 1;
    let mut k = vec![0.0; n * n];
    let mut sum = 0.0;
    for y in 0..n {
        for x in 0..n {
            let dx = x as f64 - radius as f64;
            let dy = y as f64 - radius as f64;
            let r2 = dx * dx + dy * dy;
            let v = (1.0 + r2 / (alpha * alpha).max(1e-6)).powf(-beta);
            k[y * n + x] = v;
            sum += v;
        }
    }
    if sum > 0.0 {
        for v in &mut k {
            *v /= sum;
        }
    }
    k
}

fn convolve_rgb_f64(input: &[[f64; 3]], w: usize, h: usize, kernel: &[f64], radius: usize) -> Vec<[f64; 3]> {
    let n = 2 * radius + 1;
    let mut out = vec![[0.0; 3]; w * h];
    for y in 0..h {
        for x in 0..w {
            let mut acc = [0.0_f64; 3];
            for ky in 0..n {
                for kx in 0..n {
                    let sx = (x + kx).saturating_sub(radius).min(w - 1);
                    let sy = (y + ky).saturating_sub(radius).min(h - 1);
                    let kw = kernel[ky * n + kx];
                    let p = input[sy * w + sx];
                    acc[0] += kw * p[0];
                    acc[1] += kw * p[1];
                    acc[2] += kw * p[2];
                }
            }
            out[y * w + x] = acc;
        }
    }
    out
}

fn richardson_lucy_rgb(observed: &[[u8; 3]], w: usize, h: usize, kernel: &[f64], radius: usize, iters: usize) -> Vec<[u8; 3]> {
    let obs: Vec<[f64; 3]> = observed
        .iter()
        .map(|p| [p[0] as f64 / 255.0, p[1] as f64 / 255.0, p[2] as f64 / 255.0])
        .collect();
    let mut est = obs.clone();
    for _ in 0..iters {
        let conv = convolve_rgb_f64(&est, w, h, kernel, radius);
        let mut ratio = vec![[0.0_f64; 3]; w * h];
        for i in 0..(w * h) {
            for c in 0..3 {
                ratio[i][c] = obs[i][c] / conv[i][c].max(1e-6);
            }
        }
        let corr = convolve_rgb_f64(&ratio, w, h, kernel, radius);
        for i in 0..(w * h) {
            for c in 0..3 {
                est[i][c] = (est[i][c] * corr[i][c]).clamp(0.0, 2.0);
            }
        }
    }
    est.into_iter()
        .map(|p| {
            [
                (p[0] * 255.0).clamp(0.0, 255.0) as u8,
                (p[1] * 255.0).clamp(0.0, 255.0) as u8,
                (p[2] * 255.0).clamp(0.0, 255.0) as u8,
            ]
        })
        .collect()
}

/// Box-downsample an integer superscaled image back to base resolution.
fn downsample_box(
    pixels_hi: &[[u8; 3]],
    width_hi: usize,
    height_hi: usize,
    superscale: usize,
) -> Vec<[u8; 3]> {
    assert!(superscale >= 1);
    if superscale == 1 {
        return pixels_hi.to_vec();
    }
    assert_eq!(width_hi % superscale, 0);
    assert_eq!(height_hi % superscale, 0);

    let width = width_hi / superscale;
    let height = height_hi / superscale;
    let norm = (superscale * superscale) as f64;
    let mut out = vec![[0u8; 3]; width * height];

    for y in 0..height {
        for x in 0..width {
            let mut acc = [0.0_f64; 3];
            for sy in 0..superscale {
                let yy = y * superscale + sy;
                for sx in 0..superscale {
                    let xx = x * superscale + sx;
                    let p = pixels_hi[yy * width_hi + xx];
                    acc[0] += p[0] as f64;
                    acc[1] += p[1] as f64;
                    acc[2] += p[2] as f64;
                }
            }
            out[y * width + x] = [
                (acc[0] / norm).round().clamp(0.0, 255.0) as u8,
                (acc[1] / norm).round().clamp(0.0, 255.0) as u8,
                (acc[2] / norm).round().clamp(0.0, 255.0) as u8,
            ];
        }
    }
    out
}

/// Deterministic subpixel jitter pattern.
///
/// Returns `(jx, jy)` in `[0,1)×[0,1)` for sample `i`.
fn sample_jitter(i: usize, spp: usize) -> (f64, f64) {
    if spp <= 1 {
        return (0.5, 0.5);
    }
    let n = (spp as f64).sqrt().ceil() as usize;
    let sx = i % n;
    let sy = i / n;
    (((sx as f64) + 0.5) / n as f64, ((sy as f64) + 0.5) / n as f64)
}

/// Reduce integration step near the critical impact parameter to sharpen rings.
fn adaptive_dphi_for_b(base_dphi: f64, b: f64, b_crit: f64) -> f64 {
    let rel = ((b / b_crit.max(1e-12)) - 1.0).abs();
    let scale = if rel < 0.01 {
        0.20
    } else if rel < 0.03 {
        0.35
    } else if rel < 0.08 {
        0.60
    } else {
        1.0
    };
    (base_dphi * scale).max(8e-4)
}

// ── Render & write ────────────────────────────────────────────────────────────

/// Render one view and save as PNG.
///
/// - `width_override`  — pixel width;  0 → default 1200
/// - `height_override` — pixel height; 0 → same as width (square)
/// - `blur_sigma`      — Gaussian PSF σ in pixels; 0.0 → no blur
fn render_view(
    out_dir: &Path,
    view: &View,
    width_override: usize,
    height_override: usize,
    blur_sigma: f64,
) -> PathBuf {
    use std::f64::consts::PI;
    let parse_env_f64 = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<f64>().ok());
    let parse_env_usize = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<usize>().ok());
    let parse_env_bool = |k: &str| {
        std::env::var(k).ok().is_some_and(|s| {
            matches!(s.as_str(), "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON")
        })
    };
    let width = if width_override > 0 {
        width_override
    } else {
        1200
    };
    let height = if height_override > 0 {
        height_override
    } else {
        width
    };
    let r_cam_override = parse_env_f64("BH_R_CAM_FRAC_OVERRIDE");
    let fov_override = parse_env_f64("BH_FOV_OVERRIDE");
    let max_phi_pi_override = parse_env_f64("BH_MAX_PHI_PI_OVERRIDE");
    let dphi_override = parse_env_f64("BH_DPHI_OVERRIDE");
    let az_override = parse_env_f64("BH_AZ_OVERRIDE");
    let inc_override = parse_env_f64("BH_INC_OVERRIDE");
    let detail_preset = std::env::var("BH_DETAIL_PRESET")
        .ok()
        .map(|s| s.to_ascii_lowercase());
    let mut superscale = parse_env_usize("BH_SUPERSCALE").unwrap_or(1).clamp(1, 8);
    let mut spp = parse_env_usize("BH_SPP").unwrap_or(1).clamp(1, 64);
    let mut adaptive_dphi = parse_env_bool("BH_ADAPTIVE_DPHI");
    let spectral_band = SpectralBand::from_env();
    let disk_model_env = std::env::var("BH_DISK_MODEL").ok();
    let mut disk_model = DiskModel::from_env();
    let plasma_model = PlasmaModel::from_env();
    let use_transfer_env = std::env::var("BH_USE_TRANSFER").ok();
    let mut use_transfer = use_transfer_env
        .as_ref()
        .is_some_and(|s| matches!(s.as_str(), "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON"));
    let tau_scale_env = std::env::var("BH_TAU_SCALE").ok();
    let mut tau_scale = tau_scale_env
        .as_ref()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.0)
        .max(0.0);
    let force_gr = parse_env_bool("BH_FORCE_GR");
    let effective_r_cam_frac = r_cam_override.unwrap_or(view.r_cam_frac);
    let effective_fov = fov_override.unwrap_or(view.fov);
    let mut effective_max_phi_pi = max_phi_pi_override.unwrap_or(view.max_phi_pi);
    let effective_az = az_override.unwrap_or(view.az);
    let effective_inc = inc_override.unwrap_or(view.inc);
    // Use finer integration step at high resolution so the photon ring stays crisp
    let base_dphi = if width >= 2048 {
        view.dphi.min(0.003)
    } else {
        view.dphi
    };
    let mut dphi = dphi_override.unwrap_or(base_dphi);
    if detail_preset.as_deref() == Some("imax") {
        superscale = superscale.max(2);
        spp = spp.max(4);
        adaptive_dphi = true;
        dphi = dphi.min(0.0025);
        if max_phi_pi_override.is_none() {
            effective_max_phi_pi *= 1.5;
        }
    }
    let blur_str = if blur_sigma >= 0.5 {
        format!("blur={blur_sigma:.1}")
    } else {
        "no blur".into()
    };
    eprintln!(
        "  rendering {}  (inc={:.0}°, az={:.0}°, doppler={}, rings={}, gr={}, fov={:.2}, r_cam_frac={:.3}, max_phi_pi={:.1}, {}×{}, dphi={:.4}, adaptive_dphi={}, superscale={}x, spp={}, {}) …",
        view.label, effective_inc, effective_az, view.doppler, view.ring_mode, view.gr_mode,
        effective_fov, effective_r_cam_frac, effective_max_phi_pi, width, height, dphi,
        adaptive_dphi, superscale, spp, blur_str,
    );
    let metric = if view.gr_mode || force_gr {
        GutoeMetric::schwarzschild(1.0)
    } else {
        GutoeMetric::planck_units(1.0)
    };

    let kerr_metric = if let Some(a_star) = parse_env_f64("BH_KERR_ASTAR") {
        let kerr = KerrMetric::new(metric.r_s, a_star)
            .unwrap_or_else(|| panic!("invalid BH_KERR_ASTAR={a_star}: expected |a*| <= 1"));
        let (r_plus, r_minus) = kerr.horizons();
        let r_erg_eq = kerr.ergosphere_radius(PI / 2.0);
        let r_ph_pro = kerr.equatorial_photon_orbit_radius(true);
        let r_ph_ret = kerr.equatorial_photon_orbit_radius(false);
        eprintln!(
            "    [KERR MODE] a*={:.3} r+={:.4} r-={:.4} r_erg,eq={:.4} r_ph,pro={:.4} r_ph,ret={:.4} Ω_H={:.5}",
            a_star, r_plus, r_minus, r_erg_eq, r_ph_pro, r_ph_ret, kerr.horizon_angular_velocity()
        );
        Some(kerr)
    } else {
        None
    };
    if kerr_metric.is_some() {
        // In Kerr mode, prefer volumetric flow defaults unless user explicitly overrides.
        if disk_model_env.is_none() {
            disk_model = DiskModel::Riaf;
        }
        if use_transfer_env.is_none() {
            use_transfer = true;
        }
        if tau_scale_env.is_none() {
            tau_scale = 0.35;
        }
    }
    eprintln!("    spectrum={} (BH_SPECTRUM)", spectral_band.as_label());
    eprintln!("    disk_model={} (BH_DISK_MODEL)", disk_model.as_label());
    eprintln!("    plasma_model={} (BH_PLASMA_MODEL)", plasma_model.as_label());
    eprintln!(
        "    transfer={} tau_scale={:.3} (BH_USE_TRANSFER/BH_TAU_SCALE)",
        use_transfer, tau_scale
    );
    eprintln!("    truth_mode={} (BH_TRUTH_MODE)", truth_mode().as_label());

    // Interior camera: compute r_cam in units of r_s from the horizon fraction.
    // r_cam_frac = 0 → exterior view.
    // r_cam_frac > 0 → camera at r_cam = r_cam_frac × r_horizon (always inside).
    let r_cam_rs = if effective_r_cam_frac > 0.0 {
        let r_h = metric
            .r_horizon()
            .expect("horizon must exist for interior view");
        (effective_r_cam_frac * r_h / metric.r_s).min(0.99 * r_h / metric.r_s)
    } else {
        0.0
    };
    if kerr_metric.is_some() && r_cam_rs > 0.0 {
        eprintln!("    [KERR MODE] interior camera currently uses Schwarzschild interior tracer.");
    }

    let cfg_base = RenderConfig {
        width,
        height,
        fov_rs: effective_fov,
        inclination_deg: effective_inc,
        max_phi: effective_max_phi_pi * PI,
        dphi,
    };

    let cfg_render = if superscale > 1 {
        let width_hi = width
            .checked_mul(superscale)
            .expect("BH_SUPERSCALE causes width overflow");
        let height_hi = height
            .checked_mul(superscale)
            .expect("BH_SUPERSCALE causes height overflow");
        RenderConfig {
            width: width_hi,
            height: height_hi,
            ..cfg_base.clone()
        }
    } else {
        cfg_base.clone()
    };

    let raw_hi = if spp == 1 {
        render_with_options(
            &metric,
            kerr_metric.as_ref(),
            view.disk_inner,
            view.disk_outer,
            &cfg_render,
            effective_az,
            view.doppler,
            view.ring_mode,
            view.interior_mode,
            view.core_look_mode,
            spectral_band,
            disk_model,
            plasma_model,
            use_transfer,
            tau_scale,
            adaptive_dphi,
            r_cam_rs,
            0.5,
            0.5,
        )
    } else {
        let npx = cfg_render.width * cfg_render.height;
        let mut acc = vec![[0_u32; 3]; npx];
        for si in 0..spp {
            let (jx, jy) = sample_jitter(si, spp);
            eprintln!(
                "    sample {}/{} (jitter={:.3},{:.3}) …",
                si + 1,
                spp,
                jx,
                jy
            );
            let sample = render_with_options(
                &metric,
                kerr_metric.as_ref(),
                view.disk_inner,
                view.disk_outer,
                &cfg_render,
                effective_az,
                view.doppler,
                view.ring_mode,
                view.interior_mode,
                view.core_look_mode,
                spectral_band,
                disk_model,
                plasma_model,
                use_transfer,
                tau_scale,
                adaptive_dphi,
                r_cam_rs,
                jx,
                jy,
            );
            for (a, p) in acc.iter_mut().zip(sample.iter()) {
                a[0] += p[0] as u32;
                a[1] += p[1] as u32;
                a[2] += p[2] as u32;
            }
        }
        acc.into_iter()
            .map(|a| {
                [
                    ((a[0] as f64 / spp as f64).round().clamp(0.0, 255.0)) as u8,
                    ((a[1] as f64 / spp as f64).round().clamp(0.0, 255.0)) as u8,
                    ((a[2] as f64 / spp as f64).round().clamp(0.0, 255.0)) as u8,
                ]
            })
            .collect()
    };

    let raw = if superscale > 1 {
        eprintln!(
            "    downsampling {}×{} -> {}×{} ({}x box)",
            cfg_render.width, cfg_render.height, cfg_base.width, cfg_base.height, superscale
        );
        downsample_box(&raw_hi, cfg_render.width, cfg_render.height, superscale)
    } else {
        raw_hi
    };

    // Apply Gaussian PSF blur (EHT beam simulation) if requested
    let pixels = if blur_sigma >= 0.5 {
        eprintln!("    applying Gaussian blur σ={blur_sigma:.1} …");
        gaussian_blur(&raw, cfg_base.width, cfg_base.height, blur_sigma)
    } else {
        raw
    };

    let slug_suffix = std::env::var("BH_SLUG_SUFFIX")
        .ok()
        .filter(|s| !s.trim().is_empty());
    let output_slug = match slug_suffix {
        Some(sfx) => format!("{}__{}", view.slug, sfx),
        None => view.slug.to_string(),
    };

    let ppm_path = out_dir.join(format!("{}.ppm", output_slug));
    let png_path = out_dir.join(format!("{}.png", output_slug));

    fs::write(
        &ppm_path,
        write_ppm(&pixels, cfg_base.width, cfg_base.height),
    )
    .expect("write ppm");

    let status = Command::new("convert")
        .arg(&ppm_path)
        .arg(&png_path)
        .status()
        .expect("ImageMagick convert not found");
    assert!(status.success(), "convert failed for {}", output_slug);

    // Per-render provenance sidecar for reproducibility/audit trails.
    let git_sha = Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "unknown".to_string());
    let sidecar_path = out_dir.join(format!("{}.json", output_slug));
    let mut sidecar = String::new();
    sidecar.push_str("{\n");
    sidecar.push_str(&format!("  \"slug\": \"{}\",\n", output_slug.replace('"', "\\\"")));
    sidecar.push_str(&format!("  \"label\": \"{}\",\n", view.label.replace('"', "\\\"")));
    sidecar.push_str(&format!("  \"git_sha\": \"{}\",\n", git_sha));
    sidecar.push_str(&format!("  \"png\": \"{}\",\n", png_path.display()));
    sidecar.push_str(&format!("  \"width\": {},\n", cfg_base.width));
    sidecar.push_str(&format!("  \"height\": {},\n", cfg_base.height));
    sidecar.push_str(&format!("  \"inclination_deg\": {:.6},\n", effective_inc));
    sidecar.push_str(&format!("  \"azimuth_deg\": {:.6},\n", effective_az));
    sidecar.push_str(&format!("  \"fov_rs\": {:.6},\n", effective_fov));
    sidecar.push_str(&format!("  \"r_cam_frac\": {:.6},\n", effective_r_cam_frac));
    sidecar.push_str(&format!("  \"max_phi_pi\": {:.6},\n", effective_max_phi_pi));
    sidecar.push_str(&format!("  \"dphi\": {:.8},\n", dphi));
    sidecar.push_str(&format!("  \"adaptive_dphi\": {},\n", adaptive_dphi));
    sidecar.push_str(&format!("  \"superscale\": {},\n", superscale));
    sidecar.push_str(&format!("  \"spp\": {},\n", spp));
    sidecar.push_str(&format!("  \"blur_sigma\": {:.6},\n", blur_sigma));
    sidecar.push_str(&format!("  \"disk_model\": \"{}\",\n", disk_model.as_label()));
    sidecar.push_str(&format!("  \"plasma_model\": \"{}\",\n", plasma_model.as_label()));
    sidecar.push_str(&format!("  \"spectrum\": \"{}\",\n", spectral_band.as_label()));
    sidecar.push_str(&format!("  \"use_transfer\": {},\n", use_transfer));
    sidecar.push_str(&format!("  \"tau_scale\": {:.6},\n", tau_scale));
    sidecar.push_str(&format!("  \"truth_mode\": \"{}\",\n", truth_mode().as_label()));
    sidecar.push_str(&format!("  \"gutoe_mode\": {},\n", !(view.gr_mode || force_gr)));
    sidecar.push_str(&format!("  \"kerr_astar\": {},\n", parse_env_f64("BH_KERR_ASTAR").unwrap_or(0.0)));
    sidecar.push_str(&format!("  \"starmap_path\": \"{}\",\n", std::env::var("BH_STARMAP_PATH").unwrap_or_default().replace('"', "\\\"")));
    let parity_check = parse_env_bool("BH_KERR_PARITY") || parse_env_bool("BH_VALIDATE_GPU");
    sidecar.push_str(&format!("  \"cpu_gpu_parity_check\": {}\n", parity_check));
    sidecar.push_str("}\n");
    fs::write(&sidecar_path, sidecar).expect("write provenance sidecar json");

    if parse_env_bool("BH_DEBUG_DUMP") {
        dump_debug_channels(
            out_dir,
            &output_slug,
            &metric,
            kerr_metric.as_ref(),
            view.disk_inner,
            view.disk_outer,
            &cfg_base,
            effective_az,
            view.interior_mode,
            view.core_look_mode,
            adaptive_dphi,
            r_cam_rs,
        );
    }

    fs::remove_file(&ppm_path).ok();
    eprintln!(
        "  → saved {} (+ {})",
        png_path.display(),
        sidecar_path.display()
    );
    png_path
}

fn render_view_tiled(
    out_dir: &Path,
    view: &View,
    width: usize,
    height: usize,
    tile_px: usize,
    blur_sigma: f64,
) -> PathBuf {
    use std::f64::consts::PI;
    use std::io::{Seek, SeekFrom, Write};
    let parse_env_f64 = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<f64>().ok());
    let parse_env_usize = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<usize>().ok());
    let parse_env_bool = |k: &str| {
        std::env::var(k).ok().is_some_and(|s| {
            matches!(s.as_str(), "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON")
        })
    };
    let fov_override = parse_env_f64("BH_FOV_OVERRIDE");
    let inc_override = parse_env_f64("BH_INC_OVERRIDE");
    let az_override = parse_env_f64("BH_AZ_OVERRIDE");
    let max_phi_pi_override = parse_env_f64("BH_MAX_PHI_PI_OVERRIDE");
    let effective_fov = fov_override.unwrap_or(view.fov);
    let effective_inc = inc_override.unwrap_or(view.inc);
    let effective_az = az_override.unwrap_or(view.az);
    let effective_max_phi_pi = max_phi_pi_override.unwrap_or(view.max_phi_pi);
    let spectral_band = SpectralBand::from_env();
    let disk_model = DiskModel::from_env();
    let plasma_model = PlasmaModel::from_env();
    let use_transfer = std::env::var("BH_USE_TRANSFER")
        .ok()
        .is_some_and(|s| matches!(s.as_str(), "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON"));
    let tau_scale = std::env::var("BH_TAU_SCALE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.0)
        .max(0.0);
    let force_gr = parse_env_bool("BH_FORCE_GR");
    let metric = if view.gr_mode || force_gr {
        GutoeMetric::schwarzschild(1.0)
    } else {
        GutoeMetric::planck_units(1.0)
    };
    let kerr_metric = parse_env_f64("BH_KERR_ASTAR")
        .and_then(|a| KerrMetric::new(metric.r_s, a));
    let r_cam_rs = if view.r_cam_frac > 0.0 {
        metric
            .r_horizon()
            .map(|rh| (view.r_cam_frac * rh / metric.r_s).min(0.99 * rh / metric.r_s))
            .unwrap_or(0.0)
    } else {
        0.0
    };
    let cfg = RenderConfig {
        width,
        height,
        fov_rs: effective_fov,
        inclination_deg: effective_inc,
        max_phi: effective_max_phi_pi * PI,
        dphi: view.dphi,
    };
    let mut full = vec![[0_u8; 3]; width * height];
    let output_slug = format!("{}__tiled", view.slug);
    let progress_ppm_path = out_dir.join(format!("{}.progress.ppm", output_slug));
    let mut progress_ppm = std::fs::OpenOptions::new()
        .create(true)
        .truncate(true)
        .read(true)
        .write(true)
        .open(&progress_ppm_path)
        .expect("create tiled progress ppm");
    let ppm_header = format!("P6\n{} {}\n255\n", width, height);
    progress_ppm
        .write_all(ppm_header.as_bytes())
        .expect("write tiled progress ppm header");
    let ppm_header_len = ppm_header.len() as u64;
    let zero_row = vec![0u8; width * 3];
    for _ in 0..height {
        progress_ppm
            .write_all(&zero_row)
            .expect("initialize tiled progress ppm body");
    }
    progress_ppm.flush().ok();

    let tile = tile_px.clamp(128, 4096);
    let tx_count = width.div_ceil(tile);
    let ty_count = height.div_ceil(tile);
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    let tiled_gpu = std::env::var("BH_TILED_GPU")
        .ok()
        .map(|s| !matches!(s.as_str(), "0" | "false" | "FALSE" | "off" | "OFF"))
        .unwrap_or(true);
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    let tiled_gpu = if truth_mode() != TruthMode::Off {
        eprintln!(
            "  [GPU] truth mode '{}' active; tiled render forced to CPU",
            truth_mode().as_label()
        );
        false
    } else {
        tiled_gpu
    };
    #[cfg(not(any(feature = "cuda", feature = "rocm")))]
    let tiled_gpu = false;
    let faithful_guard = parse_env_bool("BH_TILED_FAITHFUL");
    let probe_px = parse_env_usize("BH_TILED_PROBE_PX").unwrap_or(48).clamp(16, 128);
    let probe_mad_max = parse_env_f64("BH_TILED_PROBE_MAD_MAX")
        .unwrap_or(8.0)
        .clamp(0.5, 64.0);
    eprintln!(
        "  tiled render {}: {}x{} in {}x{} tiles (tile={}px, backend={})",
        view.slug,
        width,
        height,
        tx_count,
        ty_count,
        tile,
        if tiled_gpu { "gpu" } else { "cpu" }
    );
    for ty in 0..ty_count {
        for tx in 0..tx_count {
            let x0 = tx * tile;
            let y0 = ty * tile;
            let tw = (width - x0).min(tile);
            let th = (height - y0).min(tile);
            eprintln!("    tile ({},{})  x={} y={}  {}x{}", tx + 1, ty + 1, x0, y0, tw, th);
            #[cfg(any(feature = "cuda", feature = "rocm"))]
            let tile_cfg = RenderConfig {
                width: tw,
                height: th,
                fov_rs: cfg.fov_rs,
                inclination_deg: cfg.inclination_deg,
                max_phi: cfg.max_phi,
                dphi: cfg.dphi,
            };
            let mut tile_pixels = if tiled_gpu {
                #[cfg(any(feature = "cuda", feature = "rocm"))]
                {
                    render_with_options_gpu_window(
                        &metric,
                        view.disk_inner,
                        view.disk_outer,
                        width,
                        height,
                        x0,
                        y0,
                        &tile_cfg,
                        effective_az,
                        view.doppler,
                        view.ring_mode,
                        view.interior_mode,
                        view.core_look_mode,
                        spectral_band,
                        disk_model,
                        plasma_model,
                        use_transfer,
                        tau_scale,
                        true,
                        kerr_metric.as_ref(),
                        r_cam_rs,
                        0.5,
                        0.5,
                    )
                }
                #[cfg(not(any(feature = "cuda", feature = "rocm")))]
                {
                    unreachable!()
                }
            } else {
                render_with_options_cpu_window(
                    &metric,
                    kerr_metric.as_ref(),
                    view.disk_inner,
                    view.disk_outer,
                    width,
                    height,
                    x0,
                    y0,
                    tw,
                    th,
                    &cfg,
                    effective_az,
                    view.doppler,
                    view.ring_mode,
                    view.interior_mode,
                    view.core_look_mode,
                    spectral_band,
                    disk_model,
                    plasma_model,
                    use_transfer,
                    tau_scale,
                    true,
                    r_cam_rs,
                    0.5,
                    0.5,
                )
            };
            if tiled_gpu && faithful_guard {
                let pw = probe_px.min(tw);
                let ph = probe_px.min(th);
                let px0 = x0 + (tw - pw) / 2;
                let py0 = y0 + (th - ph) / 2;
                let cpu_probe = render_with_options_cpu_window(
                    &metric,
                    kerr_metric.as_ref(),
                    view.disk_inner,
                    view.disk_outer,
                    width,
                    height,
                    px0,
                    py0,
                    pw,
                    ph,
                    &cfg,
                    effective_az,
                    view.doppler,
                    view.ring_mode,
                    view.interior_mode,
                    view.core_look_mode,
                    spectral_band,
                    disk_model,
                    plasma_model,
                    use_transfer,
                    tau_scale,
                    true,
                    r_cam_rs,
                    0.5,
                    0.5,
                );
                let mut mad = 0.0_f64;
                for py in 0..ph {
                    let gy = py0 + py;
                    let ty_local = gy - y0;
                    let gpu_row = ty_local * tw;
                    let cpu_row = py * pw;
                    let tx_local0 = px0 - x0;
                    for px in 0..pw {
                        let g = tile_pixels[gpu_row + tx_local0 + px];
                        let c = cpu_probe[cpu_row + px];
                        mad += ((g[0].abs_diff(c[0]) as f64)
                            + (g[1].abs_diff(c[1]) as f64)
                            + (g[2].abs_diff(c[2]) as f64))
                            / 3.0;
                    }
                }
                mad /= (pw * ph).max(1) as f64;
                if mad > probe_mad_max {
                    eprintln!(
                        "      [faithful-guard] probe MAD {:.2} > {:.2} on tile ({},{}); rerendering tile on CPU",
                        mad,
                        probe_mad_max,
                        tx + 1,
                        ty + 1
                    );
                    tile_pixels = render_with_options_cpu_window(
                        &metric,
                        kerr_metric.as_ref(),
                        view.disk_inner,
                        view.disk_outer,
                        width,
                        height,
                        x0,
                        y0,
                        tw,
                        th,
                        &cfg,
                        effective_az,
                        view.doppler,
                        view.ring_mode,
                        view.interior_mode,
                        view.core_look_mode,
                        spectral_band,
                        disk_model,
                        plasma_model,
                        use_transfer,
                        tau_scale,
                        true,
                        r_cam_rs,
                        0.5,
                        0.5,
                    );
                }
            }
            for row in 0..th {
                let dst_off = (y0 + row) * width + x0;
                let src_off = row * tw;
                full[dst_off..dst_off + tw].copy_from_slice(&tile_pixels[src_off..src_off + tw]);
                let mut row_bytes = vec![0u8; tw * 3];
                for (i, px) in tile_pixels[src_off..src_off + tw].iter().enumerate() {
                    let o = i * 3;
                    row_bytes[o] = px[0];
                    row_bytes[o + 1] = px[1];
                    row_bytes[o + 2] = px[2];
                }
                let file_off = ppm_header_len + (((y0 + row) * width + x0) * 3) as u64;
                progress_ppm
                    .seek(SeekFrom::Start(file_off))
                    .expect("seek tiled progress ppm row");
                progress_ppm
                    .write_all(&row_bytes)
                    .expect("write tiled progress ppm row");
            }
            progress_ppm.flush().ok();
        }
    }
    let pixels = if blur_sigma >= 0.5 {
        gaussian_blur(&full, width, height, blur_sigma)
    } else {
        full
    };
    let ppm_path = out_dir.join(format!("{}.ppm", output_slug));
    let png_path = out_dir.join(format!("{}.png", output_slug));
    fs::write(&ppm_path, write_ppm(&pixels, width, height)).expect("write tiled ppm");
    let status = Command::new("convert")
        .arg(&ppm_path)
        .arg(&png_path)
        .status()
        .expect("ImageMagick convert not found");
    assert!(status.success(), "convert failed for {}", output_slug);
    fs::remove_file(&ppm_path).ok();
    eprintln!("  → saved {}", png_path.display());
    png_path
}

/// Controlled interior-core sweep for rapid parameter scans.
/// Writes files as `camera_core__r062_f014.png` where:
///   r062 = r_cam = 0.62 * r_h, f014 = fov = 1.4 * r_s
fn render_camera_core_sweep(
    out_dir: &Path,
    width_override: usize,
    height_override: usize,
    blur_sigma: f64,
) {
    fn parse_csv_f64(var: &str) -> Option<Vec<f64>> {
        let raw = std::env::var(var).ok()?;
        let vals: Vec<f64> = raw
            .split(',')
            .filter_map(|s| s.trim().parse::<f64>().ok())
            .filter(|v| *v > 0.0)
            .collect();
        if vals.is_empty() {
            None
        } else {
            Some(vals)
        }
    }

    let Some(core_view) = VIEWS.iter().find(|v| v.slug == "camera_core") else {
        eprintln!("camera_core view not found");
        return;
    };

    // Bracket "see core + context" regime.
    // Optional overrides:
    //   BH_SWEEP_R_CAMS="0.60,0.68,0.76"
    //   BH_SWEEP_FOVS="1.4,2.2,3.6"
    let r_cam_fracs = parse_csv_f64("BH_SWEEP_R_CAMS")
        .unwrap_or_else(|| vec![0.58, 0.62, 0.66, 0.70, 0.74, 0.78]);
    let fov_vals = parse_csv_f64("BH_SWEEP_FOVS").unwrap_or_else(|| vec![1.4, 1.8, 2.3, 3.0, 4.0]);
    let max_phi_pi = std::env::var("BH_SWEEP_MAX_PHI_PI")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(140.0);
    let dphi = std::env::var("BH_SWEEP_DPHI")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0015);

    std::env::set_var("BH_MAX_PHI_PI_OVERRIDE", format!("{max_phi_pi:.4}"));
    std::env::set_var("BH_DPHI_OVERRIDE", format!("{dphi:.6}"));

    let mut total = 0usize;
    for &r_cam in &r_cam_fracs {
        for &fov in &fov_vals {
            std::env::set_var("BH_R_CAM_FRAC_OVERRIDE", format!("{r_cam:.6}"));
            std::env::set_var("BH_FOV_OVERRIDE", format!("{fov:.6}"));
            std::env::set_var(
                "BH_SLUG_SUFFIX",
                format!(
                    "r{:03}_f{:03}",
                    (r_cam * 100.0).round() as i32,
                    (fov * 10.0).round() as i32
                ),
            );
            render_view(
                out_dir,
                core_view,
                width_override,
                height_override,
                blur_sigma,
            );
            total += 1;
        }
    }

    std::env::remove_var("BH_R_CAM_FRAC_OVERRIDE");
    std::env::remove_var("BH_FOV_OVERRIDE");
    std::env::remove_var("BH_MAX_PHI_PI_OVERRIDE");
    std::env::remove_var("BH_DPHI_OVERRIDE");
    std::env::remove_var("BH_SLUG_SUFFIX");

    eprintln!(
        "Core sweep complete: {total} frames in {}",
        out_dir.display()
    );
}

/// Cinematic sequence with physically traced rays and rotating accretion-flow view.
///
/// Note: this is still the Schwarzschild/GUTOE metric (non-spinning BH). The
/// apparent rotation is from disk orientation + camera choreography.
fn render_interstellar_spin(
    out_dir: &Path,
    width_override: usize,
    height_override: usize,
    blur_sigma: f64,
    frames: usize,
) {
    let Some(base_view) = VIEWS.iter().find(|v| v.slug == "sgr_astar") else {
        eprintln!("sgr_astar view not found");
        return;
    };
    let n = frames.max(1);

    let prev_spectrum = std::env::var("BH_SPECTRUM").ok();
    let prev_fixed_exposure = std::env::var("BH_FIXED_EXPOSURE").ok();
    std::env::set_var("BH_SPECTRUM", "optical");
    std::env::set_var("BH_FIXED_EXPOSURE", "1.4");
    std::env::set_var("BH_MAX_PHI_PI_OVERRIDE", "60.0");
    std::env::set_var("BH_DPHI_OVERRIDE", "0.0030");

    for i in 0..n {
        let t = i as f64 / n as f64;
        let phase = std::f64::consts::TAU * t;
        let az = 360.0 * t;
        let inc = 52.0 + 8.0 * (phase * 0.7).sin();
        let fov = 6.8 + 0.8 * (phase * 0.5).sin();

        std::env::set_var("BH_AZ_OVERRIDE", format!("{az:.6}"));
        std::env::set_var("BH_INC_OVERRIDE", format!("{inc:.6}"));
        std::env::set_var("BH_FOV_OVERRIDE", format!("{fov:.6}"));
        std::env::set_var("BH_SLUG_SUFFIX", format!("spin_{i:04}"));
        render_view(
            out_dir,
            base_view,
            width_override,
            height_override,
            blur_sigma,
        );
    }

    std::env::remove_var("BH_AZ_OVERRIDE");
    std::env::remove_var("BH_INC_OVERRIDE");
    std::env::remove_var("BH_FOV_OVERRIDE");
    std::env::remove_var("BH_MAX_PHI_PI_OVERRIDE");
    std::env::remove_var("BH_DPHI_OVERRIDE");
    std::env::remove_var("BH_SLUG_SUFFIX");
    if let Some(v) = prev_spectrum {
        std::env::set_var("BH_SPECTRUM", v);
    } else {
        std::env::remove_var("BH_SPECTRUM");
    }
    if let Some(v) = prev_fixed_exposure {
        std::env::set_var("BH_FIXED_EXPOSURE", v);
    } else {
        std::env::remove_var("BH_FIXED_EXPOSURE");
    }

    eprintln!(
        "Interstellar spin complete: {n} frames in {}",
        out_dir.display()
    );
}

/// Render a fixed camera view across spectrum bands.
///
/// Output slugs: `{view.slug}__spec_{band}`.
fn render_spectrum_sweep(
    out_dir: &Path,
    view: &View,
    width_override: usize,
    height_override: usize,
    blur_sigma: f64,
) {
    let bands: &[(SpectralBand, &str)] = &[
        (SpectralBand::Bolometric, "bolo"),
        (SpectralBand::Radio, "radio"),
        (SpectralBand::Millimeter, "mm"),
        (SpectralBand::Infrared, "ir"),
        (SpectralBand::Optical, "optical"),
        (SpectralBand::Ultraviolet, "uv"),
        (SpectralBand::Xray, "xray"),
        (SpectralBand::Gamma, "gamma"),
    ];

    let prev_spec = std::env::var("BH_SPECTRUM").ok();
    let prev_suffix = std::env::var("BH_SLUG_SUFFIX").ok();

    for &(band, slug) in bands {
        std::env::set_var("BH_SPECTRUM", band.as_label());
        std::env::set_var("BH_SLUG_SUFFIX", format!("spec_{slug}"));
        render_view(out_dir, view, width_override, height_override, blur_sigma);
    }

    if let Some(v) = prev_spec {
        std::env::set_var("BH_SPECTRUM", v);
    } else {
        std::env::remove_var("BH_SPECTRUM");
    }
    if let Some(v) = prev_suffix {
        std::env::set_var("BH_SLUG_SUFFIX", v);
    } else {
        std::env::remove_var("BH_SLUG_SUFFIX");
    }

    eprintln!("Spectrum sweep complete for {}.", view.slug);
}

fn make_diff_map_png(base: &Path, other: &Path, out: &Path) {
    let status = Command::new("convert")
        .arg(base)
        .arg(other)
        .arg("-compose")
        .arg("difference")
        .arg("-composite")
        .arg(out)
        .status()
        .expect("ImageMagick convert not found");
    assert!(status.success(), "difference-map convert failed");
}

/// For one view, render GUTOE + forced-GR per band and write residual maps.
fn render_spectrum_diff_sweep(
    out_dir: &Path,
    view: &View,
    width_override: usize,
    height_override: usize,
    blur_sigma: f64,
) {
    let bands: &[(SpectralBand, &str)] = &[
        (SpectralBand::Bolometric, "bolo"),
        (SpectralBand::Radio, "radio"),
        (SpectralBand::Millimeter, "mm"),
        (SpectralBand::Infrared, "ir"),
        (SpectralBand::Optical, "optical"),
        (SpectralBand::Ultraviolet, "uv"),
        (SpectralBand::Xray, "xray"),
        (SpectralBand::Gamma, "gamma"),
    ];

    let prev_spec = std::env::var("BH_SPECTRUM").ok();
    let prev_suffix = std::env::var("BH_SLUG_SUFFIX").ok();
    let prev_force_gr = std::env::var("BH_FORCE_GR").ok();

    for &(band, slug) in bands {
        std::env::set_var("BH_SPECTRUM", band.as_label());

        std::env::set_var("BH_FORCE_GR", "0");
        std::env::set_var("BH_SLUG_SUFFIX", format!("spec_{slug}_gutoe"));
        let gutoe_png = render_view(out_dir, view, width_override, height_override, blur_sigma);

        std::env::set_var("BH_FORCE_GR", "1");
        std::env::set_var("BH_SLUG_SUFFIX", format!("spec_{slug}_gr"));
        let gr_png = render_view(out_dir, view, width_override, height_override, blur_sigma);

        let diff_png = out_dir.join(format!("{}__spec_{}_diff.png", view.slug, slug));
        make_diff_map_png(&gutoe_png, &gr_png, &diff_png);
        eprintln!("  → saved {}", diff_png.display());
    }

    if let Some(v) = prev_spec {
        std::env::set_var("BH_SPECTRUM", v);
    } else {
        std::env::remove_var("BH_SPECTRUM");
    }
    if let Some(v) = prev_suffix {
        std::env::set_var("BH_SLUG_SUFFIX", v);
    } else {
        std::env::remove_var("BH_SLUG_SUFFIX");
    }
    if let Some(v) = prev_force_gr {
        std::env::set_var("BH_FORCE_GR", v);
    } else {
        std::env::remove_var("BH_FORCE_GR");
    }

    eprintln!("Spectrum+diff sweep complete for {}.", view.slug);
}

/// Run the current BH science campaign for target views.
/// Produces `{slug}__spec_*_{gutoe,gr,diff}.png`.
fn render_bh_campaign(
    out_dir: &Path,
    width_override: usize,
    height_override: usize,
    blur_sigma: f64,
) {
    for slug in ["m87star", "sgr_astar"] {
        let Some(view) = VIEWS.iter().find(|v| v.slug == slug) else {
            eprintln!("campaign view missing: {slug}");
            continue;
        };
        eprintln!("\n== Campaign target: {} ==", view.slug);
        render_spectrum_diff_sweep(out_dir, view, width_override, height_override, blur_sigma);
    }
}

// ── Full render with options ──────────────────────────────────────────────────
//
// Pixel loop: azimuth rotation + Doppler + ring-mode colouring.
// When az=0 and doppler=false and ring_mode=false, delegates to the fast
// tracer::render() which skips the per-pixel boxing overhead.

fn render_with_options_cpu(
    metric: &GutoeMetric,
    kerr: Option<&KerrMetric>,
    disk_inner_rs: f64,
    disk_outer_rs: f64,
    cfg: &RenderConfig,
    az_deg: f64,
    doppler: bool,
    ring_mode: bool,
    interior_mode: bool,
    core_look_mode: bool,
    spectral_band: SpectralBand,
    disk_model: DiskModel,
    plasma_model: PlasmaModel,
    use_transfer: bool,
    tau_scale: f64,
    adaptive_dphi: bool,
    r_cam_rs: f64, // 0.0 = exterior; >0 = camera inside at r_cam = r_cam_rs × r_s
    jitter_x: f64, // subpixel jitter in [0,1)
    jitter_y: f64, // subpixel jitter in [0,1)
) -> Vec<[u8; 3]> {
    let r_s = metric.r_s;
    let disk_inner = disk_inner_rs * r_s;
    let disk_outer = disk_outer_rs * r_s;
    let r_isco = 3.0 * r_s;
    let b_crit = b_critical(r_s);
    let sin_inc = cfg.inclination_deg.to_radians().sin();
    let scale = 2.0 * cfg.fov_rs * r_s / (cfg.width.min(cfg.height) as f64);
    let az_rad = az_deg.to_radians();
    let (ca, sa) = (az_rad.cos(), az_rad.sin());
    let width = cfg.width;
    let height = cfg.height;
    let r_cam = r_cam_rs * r_s; // 0.0 = exterior
    let kerr_astar = kerr.map_or(0.0, |km| km.a_star);
    let true3d = std::env::var("BH_TRUE3D_TRACE")
        .ok()
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        // Kerr should default to the full camera projection unless explicitly disabled.
        .unwrap_or(kerr.is_some());
    let z_obs = 60.0 * r_s;
    let true3d_cam = if true3d && r_cam <= 0.0 {
        let inc = cfg.inclination_deg.to_radians();
        let obs = Vec3::new(
            z_obs * inc.sin() * az_rad.cos(),
            z_obs * inc.sin() * az_rad.sin(),
            z_obs * inc.cos(),
        );
        let fwd = Vec3::new(-obs.x, -obs.y, -obs.z);
        let world_up = if inc.sin().abs() > 0.995 {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        let half_view_y = cfg.fov_rs * r_s * (height as f64 / width.min(height) as f64);
        let fov_y = 2.0 * (half_view_y / z_obs).atan();
        CameraFrame::new(
            obs,
            fwd,
            world_up,
            fov_y.clamp(1e-4, std::f64::consts::PI - 1e-4),
            width as f64 / height as f64,
        )
    } else {
        None
    };

    (0..height * width)
        .into_par_iter()
        .map(|idx| {
            let iy = idx / width;
            let ix = idx % width;
            let dx = jitter_x - 0.5;
            let dy = 0.5 - jitter_y;
            let sx = (ix as f64 - 0.5 * (width as f64 - 1.0) + dx) * scale;
            let sy = (0.5 * (height as f64 - 1.0) - iy as f64 + dy) * scale;
            let bx_raw = sx;
            // Schwarzschild tracer uses by scaled by sin(i) to encode orbit-plane
            // inclination in its reduced 2D formulation. Kerr tracer expects raw
            // screen beta with inclination supplied separately.
            let by_raw = if kerr.is_some() { sy } else { sy * sin_inc };
            let mut bx = ca * bx_raw - sa * by_raw;
            let mut by = sa * bx_raw + ca * by_raw;
            let mut b_mag = (bx * bx + by * by).sqrt();
            let dphi_ray = if adaptive_dphi {
                adaptive_dphi_for_b(cfg.dphi, b_mag, b_crit)
            } else {
                cfg.dphi
            };

            let result = if let Some(cam3d) = true3d_cam.as_ref() {
                let x = ((ix as f64 + jitter_x) / width as f64).clamp(0.0, 1.0);
                let y = ((iy as f64 + jitter_y) / height as f64).clamp(0.0, 1.0);
                let ndc_x = 2.0 * x - 1.0;
                let ndc_y = 1.0 - 2.0 * y;
                let ray_dir = cam3d.ray_dir_from_ndc(ndc_x, ndc_y);
                if let Some(red) = reduce_3d_to_axisym(cam3d.position, ray_dir) {
                    bx = red.bx;
                    by = red.by;
                    b_mag = red.b;
                }
                trace_photon_3d_schwarzschild(
                    metric,
                    disk_inner,
                    disk_outer,
                    cam3d.position,
                    ray_dir,
                    cfg.max_phi,
                    dphi_ray,
                )
            } else if r_cam > 0.0 {
                if core_look_mode {
                    trace_photon_interior_core(metric, r_cam, bx, by, cfg.max_phi, dphi_ray)
                } else {
                    trace_photon_interior(
                        metric,
                        disk_inner,
                        disk_outer,
                        r_cam,
                        bx,
                        by,
                        cfg.max_phi,
                        dphi_ray,
                    )
                }
            } else if let Some(km) = kerr {
                trace_photon_kerr(
                    km,
                    disk_inner,
                    disk_outer,
                    bx,
                    by,
                    cfg.inclination_deg,
                    cfg.max_phi * 1.2,
                    dphi_ray,
                )
            } else {
                trace_photon(metric, disk_inner, disk_outer, bx, by, cfg.max_phi, dphi_ray)
            };

            let lean_activity = if truth_mode() == TruthMode::Lean {
                if let Some(km) = kerr {
                    // Keep Rust oracle-mode numerics aligned with Lean reference defaults.
                    let lean_max_lambda = 30.0 * std::f64::consts::PI;
                    let lean_dlambda = 0.02;
                    Some(trace_photon_kerr_activity(
                        km,
                        disk_inner,
                        disk_outer,
                        bx,
                        by,
                        cfg.inclination_deg,
                        lean_max_lambda,
                        lean_dlambda,
                    ))
                } else {
                    None
                }
            } else {
                None
            };

            if let Some(c) = truth_color_from_trace(
                truth_mode(),
                result,
                lean_activity,
                b_mag,
                b_crit,
                r_s,
                r_isco,
                disk_outer,
                bx_raw,
                sin_inc,
                doppler,
                plasma_model,
                tau_scale,
                cfg.max_phi,
                kerr_astar,
            ) {
                c
            } else {
                match result {
                    TraceResult::Captured => {
                        if r_cam > 0.0 {
                            gutoe_core_color(b_mag, b_crit)
                        } else if interior_mode {
                            shadow_interior_color(bx, by, r_s)
                        } else {
                            [0, 0, 0]
                        }
                    }
                    TraceResult::Escaped { phi_total } => star_field_color(bx, by, phi_total),
                    TraceResult::DiskHit {
                        r_eff,
                        n_cross,
                        phi_orb,
                    } => {
                        if r_cam > 0.0 && core_look_mode {
                            gutoe_core_physics_color(
                                b_mag,
                                b_crit,
                                r_eff,
                                phi_orb,
                                r_cam,
                                metric.r_core(),
                            )
                        } else {
                            match disk_model {
                                DiskModel::Thin => pixel_color(
                                    r_eff,
                                    r_isco,
                                    disk_outer,
                                    r_s,
                                    bx_raw,
                                    phi_orb,
                                    sin_inc,
                                    n_cross,
                                    doppler,
                                    ring_mode,
                                    spectral_band,
                                    plasma_model,
                                    use_transfer,
                                    tau_scale,
                                    kerr_astar,
                                ),
                                DiskModel::Riaf => riaf_composite_color(
                                    r_eff,
                                    r_isco,
                                    disk_outer,
                                    r_s,
                                    bx_raw,
                                    bx,
                                    by,
                                    sin_inc,
                                    n_cross,
                                    phi_orb,
                                    doppler,
                                    ring_mode,
                                    spectral_band,
                                    plasma_model,
                                    tau_scale,
                                    kerr_astar,
                                ),
                            }
                        }
                    }
                }
            }
        })
        .collect()
}

fn render_with_options_cpu_window(
    metric: &GutoeMetric,
    kerr: Option<&KerrMetric>,
    disk_inner_rs: f64,
    disk_outer_rs: f64,
    total_width: usize,
    total_height: usize,
    x0: usize,
    y0: usize,
    tile_width: usize,
    tile_height: usize,
    cfg: &RenderConfig,
    az_deg: f64,
    doppler: bool,
    ring_mode: bool,
    interior_mode: bool,
    core_look_mode: bool,
    spectral_band: SpectralBand,
    disk_model: DiskModel,
    plasma_model: PlasmaModel,
    use_transfer: bool,
    tau_scale: f64,
    adaptive_dphi: bool,
    r_cam_rs: f64,
    jitter_x: f64,
    jitter_y: f64,
) -> Vec<[u8; 3]> {
    let r_s = metric.r_s;
    let disk_inner = disk_inner_rs * r_s;
    let disk_outer = disk_outer_rs * r_s;
    let r_isco = 3.0 * r_s;
    let b_crit = b_critical(r_s);
    let sin_inc = cfg.inclination_deg.to_radians().sin();
    let scale = 2.0 * cfg.fov_rs * r_s / (total_width.min(total_height) as f64);
    let az_rad = az_deg.to_radians();
    let (ca, sa) = (az_rad.cos(), az_rad.sin());
    let r_cam = r_cam_rs * r_s;
    let kerr_astar = kerr.map_or(0.0, |km| km.a_star);
    let true3d = std::env::var("BH_TRUE3D_TRACE")
        .ok()
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        // Kerr should default to the full camera projection unless explicitly disabled.
        .unwrap_or(kerr.is_some());
    let z_obs = 60.0 * r_s;
    let true3d_cam = if true3d && r_cam <= 0.0 {
        let inc = cfg.inclination_deg.to_radians();
        let obs = Vec3::new(
            z_obs * inc.sin() * az_rad.cos(),
            z_obs * inc.sin() * az_rad.sin(),
            z_obs * inc.cos(),
        );
        let fwd = Vec3::new(-obs.x, -obs.y, -obs.z);
        let world_up = if inc.sin().abs() > 0.995 {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        let half_view_y = cfg.fov_rs * r_s * (total_height as f64 / total_width.min(total_height) as f64);
        let fov_y = 2.0 * (half_view_y / z_obs).atan();
        CameraFrame::new(
            obs,
            fwd,
            world_up,
            fov_y.clamp(1e-4, std::f64::consts::PI - 1e-4),
            total_width as f64 / total_height as f64,
        )
    } else {
        None
    };

    (0..tile_height * tile_width)
        .into_par_iter()
        .map(|idx| {
            let iy = idx / tile_width;
            let ix = idx % tile_width;
            let gx = x0 + ix;
            let gy = y0 + iy;
            let dx = jitter_x - 0.5;
            let dy = 0.5 - jitter_y;
            let sx = (gx as f64 - 0.5 * (total_width as f64 - 1.0) + dx) * scale;
            let sy = (0.5 * (total_height as f64 - 1.0) - gy as f64 + dy) * scale;
            let bx_raw = sx;
            let by_raw = if kerr.is_some() { sy } else { sy * sin_inc };
            let mut bx = ca * bx_raw - sa * by_raw;
            let mut by = sa * bx_raw + ca * by_raw;
            let mut b_mag = (bx * bx + by * by).sqrt();
            let dphi_ray = if adaptive_dphi {
                adaptive_dphi_for_b(cfg.dphi, b_mag, b_crit)
            } else {
                cfg.dphi
            };

            let result = if let Some(cam3d) = true3d_cam.as_ref() {
                let x = ((gx as f64 + jitter_x) / total_width as f64).clamp(0.0, 1.0);
                let y = ((gy as f64 + jitter_y) / total_height as f64).clamp(0.0, 1.0);
                let ndc_x = 2.0 * x - 1.0;
                let ndc_y = 1.0 - 2.0 * y;
                let ray_dir = cam3d.ray_dir_from_ndc(ndc_x, ndc_y);
                if let Some(red) = reduce_3d_to_axisym(cam3d.position, ray_dir) {
                    bx = red.bx;
                    by = red.by;
                    b_mag = red.b;
                }
                trace_photon_3d_schwarzschild(
                    metric,
                    disk_inner,
                    disk_outer,
                    cam3d.position,
                    ray_dir,
                    cfg.max_phi,
                    dphi_ray,
                )
            } else if r_cam > 0.0 {
                if core_look_mode {
                    trace_photon_interior_core(metric, r_cam, bx, by, cfg.max_phi, dphi_ray)
                } else {
                    trace_photon_interior(
                        metric,
                        disk_inner,
                        disk_outer,
                        r_cam,
                        bx,
                        by,
                        cfg.max_phi,
                        dphi_ray,
                    )
                }
            } else if let Some(km) = kerr {
                trace_photon_kerr(
                    km,
                    disk_inner,
                    disk_outer,
                    bx,
                    by,
                    cfg.inclination_deg,
                    cfg.max_phi * 1.2,
                    dphi_ray,
                )
            } else {
                trace_photon(metric, disk_inner, disk_outer, bx, by, cfg.max_phi, dphi_ray)
            };

            let lean_activity = if truth_mode() == TruthMode::Lean {
                if let Some(km) = kerr {
                    // Keep Rust oracle-mode numerics aligned with Lean reference defaults.
                    let lean_max_lambda = 30.0 * std::f64::consts::PI;
                    let lean_dlambda = 0.02;
                    Some(trace_photon_kerr_activity(
                        km,
                        disk_inner,
                        disk_outer,
                        bx,
                        by,
                        cfg.inclination_deg,
                        lean_max_lambda,
                        lean_dlambda,
                    ))
                } else {
                    None
                }
            } else {
                None
            };

            if let Some(c) = truth_color_from_trace(
                truth_mode(),
                result,
                lean_activity,
                b_mag,
                b_crit,
                r_s,
                r_isco,
                disk_outer,
                bx_raw,
                sin_inc,
                doppler,
                plasma_model,
                tau_scale,
                cfg.max_phi,
                kerr_astar,
            ) {
                c
            } else {
                match result {
                    TraceResult::Captured => {
                        if r_cam > 0.0 {
                            gutoe_core_color(b_mag, b_crit)
                        } else if interior_mode {
                            shadow_interior_color(bx, by, r_s)
                        } else {
                            [0, 0, 0]
                        }
                    }
                    TraceResult::Escaped { phi_total } => star_field_color(bx, by, phi_total),
                    TraceResult::DiskHit {
                        r_eff,
                        n_cross,
                        phi_orb,
                    } => {
                        if r_cam > 0.0 && core_look_mode {
                            gutoe_core_physics_color(
                                b_mag,
                                b_crit,
                                r_eff,
                                phi_orb,
                                r_cam,
                                metric.r_core(),
                            )
                        } else {
                            match disk_model {
                                DiskModel::Thin => pixel_color(
                                    r_eff,
                                    r_isco,
                                    disk_outer,
                                    r_s,
                                    bx_raw,
                                    phi_orb,
                                    sin_inc,
                                    n_cross,
                                    doppler,
                                    ring_mode,
                                    spectral_band,
                                    plasma_model,
                                    use_transfer,
                                    tau_scale,
                                    kerr_astar,
                                ),
                                DiskModel::Riaf => riaf_composite_color(
                                    r_eff,
                                    r_isco,
                                    disk_outer,
                                    r_s,
                                    bx_raw,
                                    bx,
                                    by,
                                    sin_inc,
                                    n_cross,
                                    phi_orb,
                                    doppler,
                                    ring_mode,
                                    spectral_band,
                                    plasma_model,
                                    tau_scale,
                                    kerr_astar,
                                ),
                            }
                        }
                    }
                }
            }
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn dump_debug_channels(
    out_dir: &Path,
    output_slug: &str,
    metric: &GutoeMetric,
    kerr: Option<&KerrMetric>,
    disk_inner_rs: f64,
    disk_outer_rs: f64,
    cfg: &RenderConfig,
    az_deg: f64,
    interior_mode: bool,
    core_look_mode: bool,
    adaptive_dphi: bool,
    r_cam_rs: f64,
) {
    let w = cfg.width;
    let h = cfg.height;
    let r_s = metric.r_s;
    let disk_inner = disk_inner_rs * r_s;
    let disk_outer = disk_outer_rs * r_s;
    let b_crit = b_critical(r_s);
    let sin_inc = cfg.inclination_deg.to_radians().sin();
    let scale = 2.0 * cfg.fov_rs * r_s / (w.min(h) as f64);
    let az_rad = az_deg.to_radians();
    let (ca, sa) = (az_rad.cos(), az_rad.sin());
    let r_cam = r_cam_rs * r_s;

    let mut hit_class = vec![0u8; w * h];
    let mut n_cross = vec![0u8; w * h];
    let mut branch_sign = vec![0u8; w * h]; // 0=none, 64=-1, 192=+1
    let mut phi_vals = vec![0.0_f64; w * h];
    let mut r_eff_vals = vec![0.0_f64; w * h];

    for iy in 0..h {
        for ix in 0..w {
            let idx = iy * w + ix;
            let sx = (ix as f64 - 0.5 * (w as f64 - 1.0)) * scale;
            let sy = (0.5 * (h as f64 - 1.0) - iy as f64) * scale;
            let bx_raw = sx;
            let by_raw = if kerr.is_some() { sy } else { sy * sin_inc };
            let bx = ca * bx_raw - sa * by_raw;
            let by = sa * bx_raw + ca * by_raw;
            let b_mag = (bx * bx + by * by).sqrt();
            let dphi_ray = if adaptive_dphi {
                adaptive_dphi_for_b(cfg.dphi, b_mag, b_crit)
            } else {
                cfg.dphi
            };

            let (result, branch) = if r_cam > 0.0 {
                let res = if core_look_mode {
                    trace_photon_interior_core(metric, r_cam, bx, by, cfg.max_phi, dphi_ray)
                } else {
                    trace_photon_interior(
                        metric,
                        disk_inner,
                        disk_outer,
                        r_cam,
                        bx,
                        by,
                        cfg.max_phi,
                        dphi_ray,
                    )
                };
                (res, 0)
            } else if let Some(km) = kerr {
                let kd = trace_photon_kerr_debug(
                    km,
                    disk_inner,
                    disk_outer,
                    bx,
                    by,
                    cfg.inclination_deg,
                    cfg.max_phi * 1.2,
                    dphi_ray,
                );
                (kd.result, kd.branch_sign)
            } else {
                (
                    trace_photon(metric, disk_inner, disk_outer, bx, by, cfg.max_phi, dphi_ray),
                    0,
                )
            };

            match result {
                TraceResult::Captured => {
                    hit_class[idx] = if interior_mode { 1 } else { 0 };
                }
                TraceResult::Escaped { phi_total } => {
                    hit_class[idx] = 1;
                    phi_vals[idx] = phi_total;
                }
                TraceResult::DiskHit {
                    r_eff,
                    n_cross: nc,
                    phi_orb,
                } => {
                    hit_class[idx] = 2;
                    n_cross[idx] = nc.min(255) as u8;
                    phi_vals[idx] = phi_orb;
                    r_eff_vals[idx] = r_eff;
                }
            }
            branch_sign[idx] = match branch {
                -1 => 64,
                1 => 192,
                _ => 0,
            };
        }
    }

    let prefix = out_dir.join(format!("{output_slug}__debug"));
    save_png_gray_u8(&prefix.with_file_name(format!("{output_slug}__debug_hit_class.png")), &hit_class, w, h);
    save_png_gray_u8(&prefix.with_file_name(format!("{output_slug}__debug_n_cross.png")), &n_cross, w, h);
    save_png_gray_u8(
        &prefix.with_file_name(format!("{output_slug}__debug_branch_sign.png")),
        &branch_sign,
        w,
        h,
    );
    save_png_gray_f64_norm(&prefix.with_file_name(format!("{output_slug}__debug_phi.png")), &phi_vals, w, h);
    save_png_gray_f64_norm(
        &prefix.with_file_name(format!("{output_slug}__debug_r_eff.png")),
        &r_eff_vals,
        w,
        h,
    );

    let csv_path = prefix.with_file_name(format!("{output_slug}__debug_channels.csv"));
    let mut csv = String::from("x,y,hit_class,n_cross,branch_sign,phi,r_eff\n");
    for iy in 0..h {
        for ix in 0..w {
            let idx = iy * w + ix;
            let branch = match branch_sign[idx] {
                64 => -1,
                192 => 1,
                _ => 0,
            };
            csv.push_str(&format!(
                "{ix},{iy},{},{},{},{:.9},{:.9}\n",
                hit_class[idx], n_cross[idx], branch, phi_vals[idx], r_eff_vals[idx]
            ));
        }
    }
    std::fs::write(&csv_path, csv).expect("write debug csv");
    eprintln!("    [debug] wrote {}", csv_path.display());
}

fn render_with_options(
    metric: &GutoeMetric,
    kerr: Option<&KerrMetric>,
    disk_inner_rs: f64,
    disk_outer_rs: f64,
    cfg: &RenderConfig,
    az_deg: f64,
    doppler: bool,
    ring_mode: bool,
    interior_mode: bool,
    core_look_mode: bool,
    spectral_band: SpectralBand,
    disk_model: DiskModel,
    plasma_model: PlasmaModel,
    use_transfer: bool,
    tau_scale: f64,
    adaptive_dphi: bool,
    r_cam_rs: f64,
    jitter_x: f64,
    jitter_y: f64,
) -> Vec<[u8; 3]> {
    let force_cpu = std::env::var("BH_FORCE_CPU")
        .ok()
        .is_some_and(|v| v == "1" || v.eq_ignore_ascii_case("true"));
    // GPU path with optional CPU parity check.
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    {
        if force_cpu {
            eprintln!("    [GPU] BH_FORCE_CPU=1 -> using CPU path");
        } else {
        let tmode = truth_mode();
        if tmode == TruthMode::Off {
            if kerr.is_some() {
                eprintln!("    [GPU] Kerr mode enabled (experimental).");
            }
            eprintln!(
                "    [GPU] launching GPU kernel ({} × {} = {} pixels) …",
                cfg.width,
                cfg.height,
                cfg.width * cfg.height
            );
            let t0 = std::time::Instant::now();
            let gpu = render_with_options_gpu(
                metric,
                disk_inner_rs,
                disk_outer_rs,
                cfg,
                az_deg,
                doppler,
                ring_mode,
                interior_mode,
                core_look_mode,
                spectral_band,
                disk_model,
                plasma_model,
                use_transfer,
                tau_scale,
                adaptive_dphi,
                kerr,
                r_cam_rs,
                jitter_x,
                jitter_y,
            );
            eprintln!("    [GPU] done in {:.2}s", t0.elapsed().as_secs_f64());

            let parity = std::env::var("BH_KERR_PARITY")
                .ok()
                .or_else(|| std::env::var("BH_VALIDATE_GPU").ok())
                .is_some_and(|v| v == "1" || v.eq_ignore_ascii_case("true"));
            if parity {
                eprintln!("    [PARITY] running CPU reference …");
                let cpu = render_with_options_cpu(
                    metric,
                    kerr,
                    disk_inner_rs,
                    disk_outer_rs,
                    cfg,
                    az_deg,
                    doppler,
                    ring_mode,
                    interior_mode,
                    core_look_mode,
                    spectral_band,
                    disk_model,
                    plasma_model,
                    use_transfer,
                    tau_scale,
                    adaptive_dphi,
                    r_cam_rs,
                    jitter_x,
                    jitter_y,
                );
                let mut mad = 0.0_f64;
                let mut maxd = 0_u8;
                for (g, c) in gpu.iter().zip(cpu.iter()) {
                    let d0 = g[0].abs_diff(c[0]);
                    let d1 = g[1].abs_diff(c[1]);
                    let d2 = g[2].abs_diff(c[2]);
                    mad += (d0 as f64 + d1 as f64 + d2 as f64) / 3.0;
                    maxd = maxd.max(d0.max(d1.max(d2)));
                }
                mad /= gpu.len().max(1) as f64;
                eprintln!("    [PARITY] GPU↔CPU mean|Δ|={mad:.3} max|Δ|={maxd}");
            }
            return gpu;
        }
        eprintln!(
            "    [GPU] truth mode '{}' active; forcing CPU path for faithful channels",
            tmode.as_label()
        );
        }
    }

    #[allow(unreachable_code)]
    render_with_options_cpu(
        metric,
        kerr,
        disk_inner_rs,
        disk_outer_rs,
        cfg,
        az_deg,
        doppler,
        ring_mode,
        interior_mode,
        core_look_mode,
        spectral_band,
        disk_model,
        plasma_model,
        use_transfer,
        tau_scale,
        adaptive_dphi,
        r_cam_rs,
        jitter_x,
        jitter_y,
    )
}

fn pixel_diff_stats(a: &[[u8; 3]], b: &[[u8; 3]]) -> (f64, u8) {
    assert_eq!(a.len(), b.len());
    let mut mad = 0.0_f64;
    let mut maxd = 0_u8;
    for (pa, pb) in a.iter().zip(b.iter()) {
        let d0 = pa[0].abs_diff(pb[0]);
        let d1 = pa[1].abs_diff(pb[1]);
        let d2 = pa[2].abs_diff(pb[2]);
        mad += (d0 as f64 + d1 as f64 + d2 as f64) / 3.0;
        maxd = maxd.max(d0.max(d1.max(d2)));
    }
    (mad / a.len().max(1) as f64, maxd)
}

#[cfg(any(feature = "cuda", feature = "rocm"))]
fn mean_luma_img(img: &[[u8; 3]]) -> f64 {
    if img.is_empty() {
        return 0.0;
    }
    let mut sum = 0.0_f64;
    for &px in img {
        sum += luminance(px);
    }
    sum / img.len() as f64
}

#[cfg(any(feature = "cuda", feature = "rocm"))]
fn centered_luma_mad(a: &[[u8; 3]], b: &[[u8; 3]]) -> f64 {
    if a.is_empty() || b.is_empty() || a.len() != b.len() {
        return 0.0;
    }
    let ma = mean_luma_img(a);
    let mb = mean_luma_img(b);
    let mut acc = 0.0_f64;
    for (&pa, &pb) in a.iter().zip(b.iter()) {
        let da = luminance(pa) - ma;
        let db = luminance(pb) - mb;
        acc += (da - db).abs();
    }
    acc / a.len() as f64
}

#[cfg(any(feature = "cuda", feature = "rocm"))]
fn luma_mean_std(img: &[[u8; 3]]) -> (f64, f64) {
    if img.is_empty() {
        return (0.0, 1.0);
    }
    let m = mean_luma_img(img);
    let mut var = 0.0_f64;
    for &px in img {
        let d = luminance(px) - m;
        var += d * d;
    }
    var /= img.len() as f64;
    (m, var.sqrt().max(1e-12))
}

#[cfg(any(feature = "cuda", feature = "rocm"))]
fn affine_luma_mad(a: &[[u8; 3]], b: &[[u8; 3]]) -> f64 {
    if a.is_empty() || b.is_empty() || a.len() != b.len() {
        return 0.0;
    }
    let (ma, sa) = luma_mean_std(a);
    let (mb, sb) = luma_mean_std(b);
    // Affine-normalize a to b's global luminance moments:
    // a' = ((a - ma) * (sb/sa)) + mb
    let scale = sb / sa;
    let mut acc = 0.0_f64;
    for (&pa, &pb) in a.iter().zip(b.iter()) {
        let la = luminance(pa);
        let lb = luminance(pb);
        let la_aff = (la - ma) * scale + mb;
        acc += (la_aff - lb).abs();
    }
    acc / a.len() as f64
}

#[cfg(any(feature = "cuda", feature = "rocm"))]
fn bright_mask_luma_mad(a: &[[u8; 3]], b: &[[u8; 3]]) -> (f64, f64) {
    if a.is_empty() || b.is_empty() || a.len() != b.len() {
        return (0.0, 0.0);
    }
    let (ma, sa) = luma_mean_std(a);
    let (mb, sb) = luma_mean_std(b);
    let ta = ma + 2.0 * sa;
    let tb = mb + 2.0 * sb;
    let mut acc = 0.0_f64;
    let mut n = 0usize;
    for (&pa, &pb) in a.iter().zip(b.iter()) {
        let la = luminance(pa);
        let lb = luminance(pb);
        if la >= ta && lb >= tb {
            acc += (la - lb).abs();
            n += 1;
        }
    }
    if n == 0 {
        (0.0, 0.0)
    } else {
        (acc / n as f64, n as f64 / a.len() as f64)
    }
}

fn luminance(px: [u8; 3]) -> f64 {
    let r = px[0] as f64 / 255.0;
    let g = px[1] as f64 / 255.0;
    let b = px[2] as f64 / 255.0;
    0.2126 * r + 0.7152 * g + 0.0722 * b
}

fn visibility_dft(img: &[[u8; 3]], width: usize, height: usize, u: f64, v: f64) -> (f64, f64) {
    let mut re = 0.0_f64;
    let mut im = 0.0_f64;
    let mut norm = 0.0_f64;
    let w = width as f64;
    let h = height as f64;
    for y in 0..height {
        let yn = (y as f64 + 0.5) / h - 0.5;
        for x in 0..width {
            let idx = y * width + x;
            let l = luminance(img[idx]);
            let xn = (x as f64 + 0.5) / w - 0.5;
            let ph = -2.0 * std::f64::consts::PI * (u * xn + v * yn);
            re += l * ph.cos();
            im += l * ph.sin();
            norm += l;
        }
    }
    if norm > 0.0 {
        (re / norm, im / norm)
    } else {
        (0.0, 0.0)
    }
}

fn closure_phase_proxy_deg(img: &[[u8; 3]], width: usize, height: usize) -> f64 {
    // Simple EHT-style closure proxy using three synthetic baselines in uv space.
    // u,v units are cycles across the image extent.
    let (v1r, v1i) = visibility_dft(img, width, height, 8.0, 0.0);
    let (v2r, v2i) = visibility_dft(img, width, height, 0.0, 8.0);
    let (v3r, v3i) = visibility_dft(img, width, height, -8.0, -8.0);
    // Product V12 * V23 * V31
    let p12r = v1r * v2r - v1i * v2i;
    let p12i = v1r * v2i + v1i * v2r;
    let pr = p12r * v3r - p12i * v3i;
    let pi = p12r * v3i + p12i * v3r;
    pi.atan2(pr).to_degrees()
}

fn compute_kerr_image_metrics(img: &[[u8; 3]], width: usize, height: usize) -> KerrImageMetrics {
    let mut flux_left = 0.0_f64;
    let mut flux_right = 0.0_f64;
    let mut dark_mass = 0.0_f64;
    let mut dark_count = 0usize;
    let mut dark_x = 0.0_f64;
    let mut dark_y = 0.0_f64;
    let mut radial_lum = vec![0.0_f64; width.max(height)];
    let mut radial_count = vec![0_u32; width.max(height)];

    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let l = luminance(img[idx]);
            if x < width / 2 {
                flux_left += l;
            } else {
                flux_right += l;
            }
            if l < 0.045 {
                dark_count += 1;
                let w = 1.0 - (l / 0.045).clamp(0.0, 1.0);
                dark_mass += w;
                dark_x += w * x as f64;
                dark_y += w * y as f64;
            }
        }
    }

    let cx = if dark_mass > 0.0 {
        dark_x / dark_mass
    } else {
        (width as f64 - 1.0) * 0.5
    };
    let cy = if dark_mass > 0.0 {
        dark_y / dark_mass
    } else {
        (height as f64 - 1.0) * 0.5
    };

    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let l = luminance(img[idx]);
            let dx = x as f64 - cx;
            let dy = y as f64 - cy;
            let r = (dx * dx + dy * dy).sqrt();
            let ri = r.floor() as usize;
            if ri < radial_lum.len() {
                radial_lum[ri] += l;
                radial_count[ri] += 1;
            }
        }
    }

    let mut radial = Vec::new();
    for i in 0..radial_lum.len() {
        if radial_count[i] > 0 {
            radial.push(radial_lum[i] / radial_count[i] as f64);
        }
    }
    let mut shadow_r = radial
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i as f64)
        .unwrap_or((width.min(height) as f64) * 0.2);
    // Robust fallback: if radial-min degenerates near 0, use equivalent radius
    // of the dark-pixel area.
    if shadow_r < 1.0 && dark_count > 0 {
        shadow_r = ((dark_count as f64) / std::f64::consts::PI).sqrt();
    }
    let peak = radial
        .iter()
        .copied()
        .fold(0.0_f64, |m, v| if v > m { v } else { m });
    let half = 0.5 * peak;
    let mut rin = shadow_r;
    let mut rout = shadow_r;
    for (i, v) in radial.iter().enumerate() {
        let r = i as f64;
        if r <= shadow_r && *v >= half {
            rin = r;
        }
        if r >= shadow_r && *v >= half {
            rout = r;
            break;
        }
    }
    let ring_thickness = (rout - rin).max(0.0);
    let flux_asymmetry =
        (flux_left - flux_right) / (flux_left + flux_right).max(1e-12);
    let shadow_diameter = 2.0 * shadow_r;
    let closure_phase_deg = closure_phase_proxy_deg(img, width, height);

    KerrImageMetrics {
        shadow_cx: cx,
        shadow_cy: cy,
        shadow_radius: shadow_r,
        shadow_diameter,
        ring_thickness,
        flux_left,
        flux_right,
        flux_asymmetry,
        closure_phase_deg,
    }
}

fn estimate_bcrit_numeric(metric: &GutoeMetric, max_phi: f64, dphi: f64) -> f64 {
    let disk_inner = 10.0 * metric.r_s;
    let disk_outer = 3.0 * metric.r_s; // disable disk hits (inner > outer)
    let bcrit = b_critical(metric.r_s);
    let classify_capture = |b: f64| -> bool {
        matches!(
            trace_photon(metric, disk_inner, disk_outer, b, 0.0, max_phi, dphi),
            TraceResult::Captured
        )
    };

    let mut lo = 0.7 * bcrit;
    let mut hi = 1.3 * bcrit;
    // Ensure brackets: lo captured, hi escaped
    for _ in 0..24 {
        if !classify_capture(lo) {
            lo *= 0.95;
        }
        if classify_capture(hi) {
            hi *= 1.05;
        }
    }
    for _ in 0..48 {
        let mid = 0.5 * (lo + hi);
        if classify_capture(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn run_validation_report(view: &View, width: usize, height: usize) {
    use std::f64::consts::PI;
    let metric = if view.gr_mode {
        GutoeMetric::schwarzschild(1.0)
    } else {
        GutoeMetric::planck_units(1.0)
    };
    let kerr = std::env::var("BH_KERR_ASTAR")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .and_then(|a| KerrMetric::new(metric.r_s, a));
    let spectral_band = SpectralBand::from_env();
    let use_transfer = std::env::var("BH_USE_TRANSFER")
        .ok()
        .is_some_and(|s| matches!(s.as_str(), "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON"));
    let tau_scale = std::env::var("BH_TAU_SCALE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.0)
        .max(0.0);
    let disk_model = DiskModel::from_env();
    let plasma_model = PlasmaModel::from_env();
    let cfg_lo = RenderConfig {
        width,
        height,
        fov_rs: view.fov,
        inclination_deg: view.inc,
        max_phi: view.max_phi_pi * PI,
        dphi: view.dphi,
    };
    let cfg_hi = RenderConfig {
        width: width * 2,
        height: height * 2,
        ..cfg_lo.clone()
    };
    let r_cam_rs = if view.r_cam_frac > 0.0 {
        metric
            .r_horizon()
            .map(|rh| (view.r_cam_frac * rh / metric.r_s).min(0.99 * rh / metric.r_s))
            .unwrap_or(0.0)
    } else {
        0.0
    };

    let lo = render_with_options(
        &metric,
        kerr.as_ref(),
        view.disk_inner,
        view.disk_outer,
        &cfg_lo,
        view.az,
        view.doppler,
        view.ring_mode,
        view.interior_mode,
        view.core_look_mode,
        spectral_band,
        disk_model,
        plasma_model,
        use_transfer,
        tau_scale,
        true,
        r_cam_rs,
        0.5,
        0.5,
    );
    let hi = render_with_options(
        &metric,
        kerr.as_ref(),
        view.disk_inner,
        view.disk_outer,
        &cfg_hi,
        view.az,
        view.doppler,
        view.ring_mode,
        view.interior_mode,
        view.core_look_mode,
        spectral_band,
        disk_model,
        plasma_model,
        use_transfer,
        tau_scale,
        true,
        r_cam_rs,
        0.5,
        0.5,
    );
    let hi_ds = downsample_box(&hi, cfg_hi.width, cfg_hi.height, 2);
    let (mad, maxd) = pixel_diff_stats(&lo, &hi_ds);

    let b_analytic = b_critical(metric.r_s);
    let b_numeric = estimate_bcrit_numeric(&metric, 80.0 * PI, 0.0015);
    let rel = ((b_numeric - b_analytic) / b_analytic).abs();

    eprintln!("\nValidation report: {}", view.slug);
    eprintln!("  convergence MAD={mad:.4} max|Δ|={maxd}");
    eprintln!(
        "  b_crit analytic={:.8} numeric={:.8} rel_err={:.3e}",
        b_analytic, b_numeric, rel
    );
}

fn run_transfer_parity_report(view: &View, width: usize, height: usize, out_dir: &Path) {
    #[cfg(not(any(feature = "cuda", feature = "rocm")))]
    {
        let _ = (view, width, height, out_dir);
        eprintln!(
            "transfer_parity requires `--features cuda` or `--features rocm` so CPU/GPU covariant transfer can be compared."
        );
        return;
    }

    #[cfg(any(feature = "cuda", feature = "rocm"))]
    {
        use std::f64::consts::PI;
        use std::io::Write as _;
        use std::time::Instant;

        let metric = if view.gr_mode {
            GutoeMetric::schwarzschild(1.0)
        } else {
            GutoeMetric::planck_units(1.0)
        };
        let kerr = std::env::var("BH_KERR_ASTAR")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .and_then(|a| KerrMetric::new(metric.r_s, a));
        let cfg = RenderConfig {
            width,
            height,
            fov_rs: view.fov,
            inclination_deg: view.inc,
            max_phi: view.max_phi_pi * PI,
            dphi: view.dphi,
        };
        let r_cam_rs = if view.r_cam_frac > 0.0 {
            metric
                .r_horizon()
                .map(|rh| (view.r_cam_frac * rh / metric.r_s).min(0.99 * rh / metric.r_s))
                .unwrap_or(0.0)
        } else {
            0.0
        };

        let backend_tag = std::env::var("BH_BACKEND_TAG").ok().unwrap_or_else(|| {
            if cfg!(all(feature = "cuda", not(feature = "rocm"))) {
                "cuda".to_string()
            } else if cfg!(all(feature = "rocm", not(feature = "cuda"))) {
                "rocm".to_string()
            } else if cfg!(all(feature = "cuda", feature = "rocm")) {
                "multi".to_string()
            } else {
                "gpu".to_string()
            }
        });
        let csv_path = out_dir.join(format!("transfer_parity_{}_{}.csv", view.slug, backend_tag));
        let tmp_path =
            out_dir.join(format!("transfer_parity_{}_{}.csv.tmp", view.slug, backend_tag));
        let mut f = std::fs::File::create(&tmp_path).expect("create transfer parity tmp csv");
        writeln!(
            f,
            "backend,disk_model,use_transfer,tau_scale,mad,max_delta,centered_luma_mad,affine_luma_mad,bright_mask_luma_mad,bright_mask_coverage,gpu_mean_luma,cpu_mean_luma,gpu_delta_from_base,cpu_delta_from_base,transfer_delta_parity_abs,gpu_ms,cpu_ms,width,height"
        )
        .expect("write transfer parity header");

        let combos = [
            (DiskModel::Thin, false, 0.0),
            (DiskModel::Thin, true, 1.0),
            (DiskModel::Thin, true, 1.5),
            (DiskModel::Riaf, false, 0.0),
            (DiskModel::Riaf, true, 1.0),
            (DiskModel::Riaf, true, 1.5),
        ];

        eprintln!("\nTransfer parity report: {} @ {}x{}", view.slug, width, height);
        let prev_legacy = std::env::var("BH_PARITY_LEGACY_STARS").ok();
        std::env::set_var("BH_PARITY_LEGACY_STARS", "1");
        let mut thin_gpu_base_luma: Option<f64> = None;
        let mut thin_cpu_base_luma: Option<f64> = None;
        let mut riaf_gpu_base_luma: Option<f64> = None;
        let mut riaf_cpu_base_luma: Option<f64> = None;
        for (disk_model, use_transfer, tau_scale) in combos {
            let t_gpu = Instant::now();
            let gpu = render_with_options_gpu(
                &metric,
                view.disk_inner,
                view.disk_outer,
                &cfg,
                view.az,
                view.doppler,
                view.ring_mode,
                view.interior_mode,
                view.core_look_mode,
                SpectralBand::Bolometric,
                disk_model,
                PlasmaModel::Grmhd,
                use_transfer,
                tau_scale,
                false,
                kerr.as_ref(),
                r_cam_rs,
                0.5,
                0.5,
            );
            let gpu_ms = t_gpu.elapsed().as_secs_f64() * 1000.0;
            let t_cpu = Instant::now();
            let cpu = render_with_options_cpu(
                &metric,
                kerr.as_ref(),
                view.disk_inner,
                view.disk_outer,
                &cfg,
                view.az,
                view.doppler,
                view.ring_mode,
                view.interior_mode,
                view.core_look_mode,
                SpectralBand::Bolometric,
                disk_model,
                PlasmaModel::Grmhd,
                use_transfer,
                tau_scale,
                false,
                r_cam_rs,
                0.5,
                0.5,
            );
            let cpu_ms = t_cpu.elapsed().as_secs_f64() * 1000.0;
            let (mad, maxd) = pixel_diff_stats(&gpu, &cpu);
            let centered = centered_luma_mad(&gpu, &cpu);
            let affine = affine_luma_mad(&gpu, &cpu);
            let (bright_mask_mad, bright_mask_cov) = bright_mask_luma_mad(&gpu, &cpu);
            let gpu_mean_luma = mean_luma_img(&gpu);
            let cpu_mean_luma = mean_luma_img(&cpu);
            let (gpu_base, cpu_base) = match disk_model {
                DiskModel::Thin => (&mut thin_gpu_base_luma, &mut thin_cpu_base_luma),
                DiskModel::Riaf => (&mut riaf_gpu_base_luma, &mut riaf_cpu_base_luma),
            };
            if !use_transfer {
                *gpu_base = Some(gpu_mean_luma);
                *cpu_base = Some(cpu_mean_luma);
            }
            let gpu_delta_from_base = gpu_mean_luma - gpu_base.unwrap_or(gpu_mean_luma);
            let cpu_delta_from_base = cpu_mean_luma - cpu_base.unwrap_or(cpu_mean_luma);
            let transfer_delta_parity_abs = (gpu_delta_from_base - cpu_delta_from_base).abs();
            eprintln!(
                "  {} transfer={} tau={:.2} -> MAD={:.4} max|Δ|={} centered={:.5} affine={:.5} bright={:.5} cov={:.4} | luma gpu={:.5} cpu={:.5} Δbase gpu={:+.5} cpu={:+.5} | |Δtransfer|={:.5} | gpu_ms={:.3} cpu_ms={:.3}",
                disk_model.as_label(),
                use_transfer,
                tau_scale,
                mad,
                maxd,
                centered,
                affine,
                bright_mask_mad,
                bright_mask_cov,
                gpu_mean_luma,
                cpu_mean_luma,
                gpu_delta_from_base,
                cpu_delta_from_base,
                transfer_delta_parity_abs,
                gpu_ms,
                cpu_ms
            );
            writeln!(
                f,
                "{},{},{},{:.6},{:.9},{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.3},{:.3},{},{}",
                backend_tag,
                disk_model.as_label(),
                use_transfer,
                tau_scale,
                mad,
                maxd,
                centered,
                affine,
                bright_mask_mad,
                bright_mask_cov,
                gpu_mean_luma,
                cpu_mean_luma,
                gpu_delta_from_base,
                cpu_delta_from_base,
                transfer_delta_parity_abs,
                gpu_ms,
                cpu_ms,
                width,
                height
            )
            .expect("write transfer parity row");
        }
        if let Some(v) = prev_legacy {
            std::env::set_var("BH_PARITY_LEGACY_STARS", v);
        } else {
            std::env::remove_var("BH_PARITY_LEGACY_STARS");
        }

        std::fs::rename(&tmp_path, &csv_path).expect("commit transfer parity csv");
        eprintln!("  wrote {}", csv_path.display());
    }
}

fn parse_csv_f64_env(var: &str, default: &[f64]) -> Vec<f64> {
    std::env::var(var)
        .ok()
        .map(|s| {
            s.split(',')
                .filter_map(|v| v.trim().parse::<f64>().ok())
                .collect::<Vec<_>>()
        })
        .filter(|v| !v.is_empty())
        .unwrap_or_else(|| default.to_vec())
}

fn run_kerr_metrics_report(view: &View, width: usize, height: usize, out_dir: &Path) {
    use std::f64::consts::PI;
    use std::io::Write as _;

    let metric = if view.gr_mode {
        GutoeMetric::schwarzschild(1.0)
    } else {
        GutoeMetric::planck_units(1.0)
    };
    let spectral_band = SpectralBand::from_env();
    let use_transfer = std::env::var("BH_USE_TRANSFER")
        .ok()
        .is_some_and(|s| matches!(s.as_str(), "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON"));
    let tau_scale = std::env::var("BH_TAU_SCALE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.0)
        .max(0.0);
    let disk_model = DiskModel::from_env();
    let plasma_model = PlasmaModel::from_env();
    let spins = parse_csv_f64_env("BH_KERR_SWEEP", &[0.0, 0.3, 0.6, 0.9]);

    let cfg = RenderConfig {
        width,
        height,
        fov_rs: view.fov,
        inclination_deg: view.inc,
        max_phi: view.max_phi_pi * PI,
        dphi: view.dphi,
    };
    let r_cam_rs = if view.r_cam_frac > 0.0 {
        metric
            .r_horizon()
            .map(|rh| (view.r_cam_frac * rh / metric.r_s).min(0.99 * rh / metric.r_s))
            .unwrap_or(0.0)
    } else {
        0.0
    };
    let csv_path = out_dir.join(format!("kerr_metrics_{}.csv", view.slug));
    let mut f = std::fs::File::create(&csv_path).expect("create kerr metrics csv");
    writeln!(
        f,
        "a_star,shadow_cx_px,shadow_cy_px,shadow_radius_px,shadow_diameter_px,ring_thickness_px,flux_left,flux_right,flux_asymmetry,closure_phase_deg"
    )
    .expect("write csv header");

    eprintln!("\nKerr metrics sweep: {} @ {}x{}", view.slug, width, height);
    for a in spins {
        let kerr = KerrMetric::new(metric.r_s, a).expect("invalid sweep a*");
        let img = render_with_options(
            &metric,
            Some(&kerr),
            view.disk_inner,
            view.disk_outer,
            &cfg,
            view.az,
            view.doppler,
            view.ring_mode,
            view.interior_mode,
            view.core_look_mode,
            spectral_band,
            disk_model,
            plasma_model,
            use_transfer,
            tau_scale,
            true,
            r_cam_rs,
            0.5,
            0.5,
        );
        let m = compute_kerr_image_metrics(&img, width, height);
        writeln!(
            f,
            "{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.9},{:.9},{:.9},{:.6}",
            a,
            m.shadow_cx,
            m.shadow_cy,
            m.shadow_radius,
            m.shadow_diameter,
            m.ring_thickness,
            m.flux_left,
            m.flux_right,
            m.flux_asymmetry,
            m.closure_phase_deg,
        )
        .expect("write csv row");
        eprintln!(
            "  a*={:.2}  cx={:.2} cy={:.2}  D_sh={:.2}  ring_w={:.2}  A_LR={:+.4}  CP={:+.2}°",
            a,
            m.shadow_cx,
            m.shadow_cy,
            m.shadow_diameter,
            m.ring_thickness,
            m.flux_asymmetry,
            m.closure_phase_deg
        );
    }
    eprintln!("  wrote {}", csv_path.display());
}

fn run_wave_samples_report(view: &View, width: usize, height: usize, out_dir: &Path) {
    use std::f64::consts::PI;
    use std::io::Write as _;

    let force_gr = std::env::var("BH_FORCE_GR")
        .ok()
        .is_some_and(|s| matches!(s.as_str(), "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON"));
    let metric = if view.gr_mode || force_gr {
        GutoeMetric::schwarzschild(1.0)
    } else {
        GutoeMetric::planck_units(1.0)
    };
    let a_star = std::env::var("BH_KERR_ASTAR")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.9);
    let kerr = KerrMetric::new(metric.r_s, a_star)
        .unwrap_or_else(|| panic!("invalid BH_KERR_ASTAR={a_star}: expected |a*| <= 1"));

    let fov_rs = std::env::var("BH_FOV_OVERRIDE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(view.fov);
    let max_phi_pi = std::env::var("BH_MAX_PHI_PI_OVERRIDE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(view.max_phi_pi);
    let dphi = std::env::var("BH_DPHI_OVERRIDE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(view.dphi);
    let inc_deg = std::env::var("BH_INC_OVERRIDE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(view.inc);
    let az_deg = std::env::var("BH_AZ_OVERRIDE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(view.az);
    let disk_inner = view.disk_inner * metric.r_s;
    let disk_outer = view.disk_outer * metric.r_s;
    let max_lambda = max_phi_pi * PI * 1.2;

    let scale = 2.0 * fov_rs * metric.r_s / (width.min(height) as f64);
    let az = az_deg.to_radians();
    let (ca, sa) = (az.cos(), az.sin());
    let r_isco = 3.0 * metric.r_s;

    #[derive(Default, Clone, Copy)]
    struct RingAgg {
        count: u64,
        re_sum: f64,
        im_sum: f64,
        incoh_sum: f64,
    }

    let mut ring_agg: BTreeMap<u32, RingAgg> = BTreeMap::new();
    let csv_path = out_dir.join(format!("wave_samples_{}.csv", view.slug));
    let mut csv = std::fs::File::create(&csv_path).expect("create wave csv");
    writeln!(
        csv,
        "x,y,sample_idx,ring,hit_class,r_eff,phi_orb,phase,amp,re,im,branch_sign"
    )
    .expect("write wave csv header");

    for iy in 0..height {
        for ix in 0..width {
            let sx = (ix as f64 - 0.5 * (width as f64 - 1.0)) * scale;
            let sy = (0.5 * (height as f64 - 1.0) - iy as f64) * scale;
            let bx = ca * sx - sa * sy;
            let by = sa * sx + ca * sy;
            let kd = trace_photon_kerr_debug(
                &kerr,
                disk_inner,
                disk_outer,
                bx,
                by,
                inc_deg,
                max_lambda,
                dphi,
            );
            let samples = trace_photon_kerr_wave_samples(
                &kerr,
                disk_inner,
                disk_outer,
                bx,
                by,
                inc_deg,
                max_lambda,
                dphi,
            );

            // Preserve compatibility with selected-branch diagnostics by writing
            // both admissible branch samples and then coherent bucket sums.
            for (si, s) in samples.iter().enumerate() {
                let mut ring = 0_u32;
                let mut r_eff = 0.0_f64;
                let mut phi_orb = 0.0_f64;
                let mut hit_class = "escaped";
                let mut amp = 0.0_f64;

                if let TraceResult::DiskHit {
                    r_eff: re,
                    phi_orb: po,
                    n_cross,
                } = s.result
                {
                    ring = n_cross;
                    r_eff = re;
                    phi_orb = po;
                    hit_class = "disk";
                    let radial = (r_isco / re.max(1e-9)).powf(0.75);
                    let order_fade = 0.7_f64.powi(n_cross.saturating_sub(1) as i32);
                    amp = (radial * order_fade).max(0.0);
                } else if matches!(s.result, TraceResult::Captured) {
                    hit_class = "captured";
                }

                let phase = s.phase_accum;
                let re = amp * phase.cos();
                let im = amp * phase.sin();
                if ring > 0 {
                    let e = ring_agg.entry(ring).or_default();
                    e.count += 1;
                    e.re_sum += re;
                    e.im_sum += im;
                    e.incoh_sum += amp * amp;
                }
                writeln!(
                    csv,
                    "{ix},{iy},{si},{ring},{hit_class},{r_eff:.9},{phi_orb:.9},{phase:.9},{amp:.9},{re:.9},{im:.9},{}",
                    s.branch_sign
                )
                .expect("write wave csv row");
            }

            // Record selected-branch marker so downstream tools can compare
            // coherent two-branch sum vs legacy single-branch pick.
            writeln!(
                csv,
                "{ix},{iy},selected,0,selected,0.000000000,0.000000000,{:.9},0.000000000,0.000000000,0.000000000,{}",
                kd.phase_accum, kd.branch_sign
            )
            .expect("write selected marker row");
        }
    }

    let json_path = out_dir.join(format!("wave_samples_{}.json", view.slug));
    let mut js = String::new();
    js.push_str("{\n");
    js.push_str(&format!("  \"slug\": \"{}\",\n", view.slug));
    js.push_str(&format!("  \"width\": {},\n", width));
    js.push_str(&format!("  \"height\": {},\n", height));
    js.push_str(&format!("  \"a_star\": {:.9},\n", a_star));
    js.push_str(&format!("  \"inclination_deg\": {:.9},\n", inc_deg));
    js.push_str(&format!("  \"azimuth_deg\": {:.9},\n", az_deg));
    js.push_str("  \"rings\": [\n");
    let mut first = true;
    for (ring, agg) in ring_agg {
        if !first {
            js.push_str(",\n");
        }
        first = false;
        let coh_i = agg.re_sum * agg.re_sum + agg.im_sum * agg.im_sum;
        js.push_str(&format!(
            "    {{\"ring\": {}, \"count\": {}, \"re_sum\": {:.9}, \"im_sum\": {:.9}, \"coherent_intensity\": {:.9}, \"incoherent_intensity\": {:.9}}}",
            ring, agg.count, agg.re_sum, agg.im_sum, coh_i, agg.incoh_sum
        ));
    }
    js.push_str("\n  ]\n}\n");
    std::fs::write(&json_path, js).expect("write wave json");

    eprintln!("wave_samples: wrote {}", csv_path.display());
    eprintln!("wave_samples: wrote {}", json_path.display());
}

fn run_wave_coherent_report(view: &View, width: usize, height: usize, out_dir: &Path) {
    use std::f64::consts::PI;
    use std::io::Write as _;

    let force_gr = std::env::var("BH_FORCE_GR")
        .ok()
        .is_some_and(|s| matches!(s.as_str(), "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON"));
    let metric = if view.gr_mode || force_gr {
        GutoeMetric::schwarzschild(1.0)
    } else {
        GutoeMetric::planck_units(1.0)
    };
    let a_star = std::env::var("BH_KERR_ASTAR")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.9);
    let kerr = KerrMetric::new(metric.r_s, a_star)
        .unwrap_or_else(|| panic!("invalid BH_KERR_ASTAR={a_star}: expected |a*| <= 1"));
    let fov_rs = std::env::var("BH_FOV_OVERRIDE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(view.fov);
    let max_phi_pi = std::env::var("BH_MAX_PHI_PI_OVERRIDE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(view.max_phi_pi);
    let dphi = std::env::var("BH_DPHI_OVERRIDE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(view.dphi);
    let inc_deg = std::env::var("BH_INC_OVERRIDE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(view.inc);
    let az_deg = std::env::var("BH_AZ_OVERRIDE")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(view.az);
    let disk_inner = view.disk_inner * metric.r_s;
    let disk_outer = view.disk_outer * metric.r_s;
    let max_lambda = max_phi_pi * PI * 1.2;
    let scale = 2.0 * fov_rs * metric.r_s / (width.min(height) as f64);
    let az = az_deg.to_radians();
    let (ca, sa) = (az.cos(), az.sin());
    let r_isco = 3.0 * metric.r_s;
    let n_lambda = std::env::var("BH_WAVE_BAND_SAMPLES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(9)
        .clamp(1, 64);
    let band_frac = std::env::var("BH_WAVE_BANDWIDTH_FRAC")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.08)
        .clamp(0.0, 0.5);
    let mix = std::env::var("BH_WAVE_MIX")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.35)
        .clamp(0.0, 1.0);
    let overlay_gain = std::env::var("BH_WAVE_OVERLAY_GAIN")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.30)
        .clamp(0.0, 2.0);
    let psf_sigma = std::env::var("BH_WAVE_PSF_SIGMA")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.2)
        .max(0.0);
    let min_ring_order = std::env::var("BH_WAVE_MIN_RING_ORDER")
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(2);
    let max_ring_order = std::env::var("BH_WAVE_MAX_RING_ORDER")
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(8)
        .max(min_ring_order);
    let ring_softness = std::env::var("BH_WAVE_RING_SOFTNESS")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.5)
        .max(0.1);
    let spectral_band = SpectralBand::from_env();
    let max_crossings = std::env::var("BH_WAVE_MAX_CROSSINGS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(6)
        .clamp(1, 8);

    let mut coh = vec![0.0_f64; width * height];
    let mut incoh = vec![0.0_f64; width * height];
    let mut contrast = vec![0.0_f64; width * height];
    let mut observed = vec![0.0_f64; width * height];
    let mut coh_max = 0.0_f64;
    let mut incoh_max = 0.0_f64;
    let mut contrast_max = 0.0_f64;
    let mut observed_max = 0.0_f64;

    for iy in 0..height {
        for ix in 0..width {
            let idx = iy * width + ix;
            let sx = (ix as f64 - 0.5 * (width as f64 - 1.0)) * scale;
            let sy = (0.5 * (height as f64 - 1.0) - iy as f64) * scale;
            let bx = ca * sx - sa * sy;
            let by = sa * sx + ca * sy;
            let branches = trace_photon_kerr_wave_branches(
                &kerr,
                disk_inner,
                disk_outer,
                bx,
                by,
                inc_deg,
                max_lambda,
                dphi,
                max_crossings,
            );
            let mut amps = Vec::with_capacity(2 * max_crossings);
            let mut phases = Vec::with_capacity(2 * max_crossings);
            let mut incoh_sum = 0.0_f64;
            for b in &branches {
                for ci in 0..b.count {
                    let c = b.crossings[ci];
                    let n = c.n_cross as f64;
                    let nmin = min_ring_order as f64;
                    let nmax = max_ring_order as f64;
                    // Soft ring window avoids hard cutoffs in the wave field.
                    let ring_weight = if n < nmin {
                        (-(nmin - n) / ring_softness).exp()
                    } else if n > nmax {
                        (-(n - nmax) / ring_softness).exp()
                    } else {
                        1.0
                    };
                    let radial = (r_isco / c.r_eff.max(1e-9)).powf(0.75);
                    let order_fade = 0.7_f64.powi(c.n_cross.saturating_sub(1) as i32);
                    let amp = (radial * order_fade * ring_weight).max(0.0);
                    amps.push(amp);
                    phases.push(c.phase_accum);
                    incoh_sum += amp * amp;
                }
            }
            // Finite detector bandwidth: average coherent intensity over nearby
            // wavelengths instead of a single monochromatic phase.
            let mut coh_avg = 0.0_f64;
            for li in 0..n_lambda {
                let t = if n_lambda <= 1 {
                    0.5
                } else {
                    li as f64 / (n_lambda as f64 - 1.0)
                };
                let f = 1.0 - band_frac + 2.0 * band_frac * t;
                let mut re_sum = 0.0_f64;
                let mut im_sum = 0.0_f64;
                for j in 0..amps.len() {
                    let ph = phases[j] * f;
                    re_sum += amps[j] * ph.cos();
                    im_sum += amps[j] * ph.sin();
                }
                coh_avg += re_sum * re_sum + im_sum * im_sum;
            }
            let coh_i = coh_avg / (n_lambda as f64).max(1.0);
            let c = if incoh_sum > 1e-12 {
                (coh_i / incoh_sum).clamp(0.0, 4.0)
            } else {
                0.0
            };
            // Observer-facing output: mostly incoherent intensity with optional
            // coherent modulation mixed in.
            let obs = incoh_sum * (1.0 + mix * (c - 1.0)).max(0.0);
            coh[idx] = coh_i;
            incoh[idx] = incoh_sum;
            contrast[idx] = c;
            observed[idx] = obs;
            coh_max = coh_max.max(coh_i);
            incoh_max = incoh_max.max(incoh_sum);
            contrast_max = contrast_max.max(c);
            observed_max = observed_max.max(obs);
        }
    }

    let coh_png = out_dir.join(format!("wave_coherent_{}__intensity.png", view.slug));
    let contrast_png = out_dir.join(format!("wave_coherent_{}__contrast.png", view.slug));
    let observed_png = out_dir.join(format!("wave_coherent_{}__observed.png", view.slug));
    let observed_color_png = out_dir.join(format!("wave_coherent_{}__observed_color.png", view.slug));
    let instrument_png = out_dir.join(format!("wave_coherent_{}__instrument.png", view.slug));
    save_png_gray_f64_norm(&coh_png, &coh, width, height);
    save_png_gray_f64_norm(&contrast_png, &contrast, width, height);
    save_png_gray_f64_norm(&observed_png, &observed, width, height);
    let [tr, tg, tb] = band_tint(spectral_band);
    let obs_scale = 1.0 / observed_max.max(1e-9);
    let observed_rgb: Vec<[u8; 3]> = observed
        .iter()
        .map(|&v| {
            let b = (v * obs_scale).clamp(0.0, 1.0);
            [
                (255.0 * b.powf(0.35) * tr).clamp(0.0, 255.0) as u8,
                (210.0 * b.powf(0.60) * tg).clamp(0.0, 255.0) as u8,
                (130.0 * b.powf(1.60) * tb).clamp(0.0, 255.0) as u8,
            ]
        })
        .collect();
    save_png_rgb(&observed_color_png, &observed_rgb, width, height);

    // Synthetic observation: start from physical RIAF base image, then apply
    // wave modulation from ring-restricted coherent field.
    let cfg = RenderConfig {
        width,
        height,
        fov_rs,
        inclination_deg: inc_deg,
        max_phi: max_phi_pi * PI,
        dphi,
    };
    let base_rgb = render_with_options(
        &metric,
        Some(&kerr),
        view.disk_inner,
        view.disk_outer,
        &cfg,
        az_deg,
        view.doppler,
        false,
        view.interior_mode,
        view.core_look_mode,
        spectral_band,
        DiskModel::Riaf,
        PlasmaModel::from_env(),
        true,
        0.35,
        true,
        0.0,
        0.5,
        0.5,
    );
    let obs_gray_u8: Vec<[u8; 3]> = observed
        .iter()
        .map(|&v| {
            let g = (255.0 * (v * obs_scale).clamp(0.0, 1.0)).round().clamp(0.0, 255.0) as u8;
            [g, g, g]
        })
        .collect();
    let obs_psf = if psf_sigma > 0.0 {
        gaussian_blur(&obs_gray_u8, width, height, psf_sigma)
    } else {
        obs_gray_u8
    };
    let instrument_rgb: Vec<[u8; 3]> = base_rgb
        .iter()
        .zip(obs_psf.iter())
        .map(|(bpx, wpx)| {
            let w = (wpx[0] as f64 / 255.0).clamp(0.0, 1.0);
            let gain = (1.0 + overlay_gain * (w - 0.5)).clamp(0.5, 1.5);
            [
                ((bpx[0] as f64) * gain).round().clamp(0.0, 255.0) as u8,
                ((bpx[1] as f64) * gain).round().clamp(0.0, 255.0) as u8,
                ((bpx[2] as f64) * gain).round().clamp(0.0, 255.0) as u8,
            ]
        })
        .collect();
    save_png_rgb(&instrument_png, &instrument_rgb, width, height);

    let cx = 0.5 * (width as f64 - 1.0);
    let cy = 0.5 * (height as f64 - 1.0);
    let rmax = ((cx * cx + cy * cy).sqrt().ceil() as usize).max(1);
    let mut sum_obs = vec![0.0_f64; rmax + 1];
    let mut sum_con = vec![0.0_f64; rmax + 1];
    let mut cnt = vec![0_u64; rmax + 1];
    for iy in 0..height {
        for ix in 0..width {
            let dx = ix as f64 - cx;
            let dy = iy as f64 - cy;
            let rb = ((dx * dx + dy * dy).sqrt().round() as usize).min(rmax);
            let idx = iy * width + ix;
            sum_obs[rb] += observed[idx];
            sum_con[rb] += contrast[idx];
            cnt[rb] += 1;
        }
    }
    let radial_csv = out_dir.join(format!("wave_coherent_{}__radial_profile.csv", view.slug));
    let mut rf = std::fs::File::create(&radial_csv).expect("create wave radial csv");
    writeln!(rf, "radius_px,observed_mean,contrast_mean,count").expect("write radial header");
    let mut peak_r = 0usize;
    let mut peak_val = f64::NEG_INFINITY;
    for r in 0..=rmax {
        let c = cnt[r] as f64;
        let obs_m = if c > 0.0 { sum_obs[r] / c } else { 0.0 };
        let con_m = if c > 0.0 { sum_con[r] / c } else { 0.0 };
        if obs_m > peak_val {
            peak_val = obs_m;
            peak_r = r;
        }
        writeln!(rf, "{r},{obs_m:.9},{con_m:.9},{}", cnt[r]).expect("write radial row");
    }

    let mean = |v: &[f64]| -> f64 {
        if v.is_empty() {
            0.0
        } else {
            v.iter().sum::<f64>() / v.len() as f64
        }
    };
    let summary_path = out_dir.join(format!("wave_coherent_{}.json", view.slug));
    let summary = format!(
        "{{\n  \"slug\": \"{}\",\n  \"width\": {},\n  \"height\": {},\n  \"a_star\": {:.9},\n  \"inclination_deg\": {:.9},\n  \"azimuth_deg\": {:.9},\n  \"band_samples\": {},\n  \"bandwidth_frac\": {:.9},\n  \"wave_mix\": {:.9},\n  \"wave_overlay_gain\": {:.9},\n  \"wave_psf_sigma\": {:.9},\n  \"wave_min_ring_order\": {},\n  \"wave_max_ring_order\": {},\n  \"wave_ring_softness\": {:.9},\n  \"wave_max_crossings\": {},\n  \"spectrum\": \"{}\",\n  \"coherent_mean\": {:.9},\n  \"coherent_max\": {:.9},\n  \"incoherent_mean\": {:.9},\n  \"incoherent_max\": {:.9},\n  \"contrast_mean\": {:.9},\n  \"contrast_max\": {:.9},\n  \"observed_mean\": {:.9},\n  \"observed_max\": {:.9},\n  \"radial_peak_radius_px\": {},\n  \"radial_peak_observed\": {:.9},\n  \"radial_profile_csv\": \"{}\"\n}}\n",
        view.slug,
        width,
        height,
        a_star,
        inc_deg,
        az_deg,
        n_lambda,
        band_frac,
        mix,
        overlay_gain,
        psf_sigma,
        min_ring_order,
        max_ring_order,
        ring_softness,
        max_crossings,
        spectral_band.as_label(),
        mean(&coh),
        coh_max,
        mean(&incoh),
        incoh_max,
        mean(&contrast),
        contrast_max,
        mean(&observed),
        observed_max,
        peak_r,
        peak_val,
        radial_csv.display()
    );
    std::fs::write(&summary_path, summary).expect("write wave coherent summary");
    eprintln!("wave_coherent: wrote {}", coh_png.display());
    eprintln!("wave_coherent: wrote {}", contrast_png.display());
    eprintln!("wave_coherent: wrote {}", observed_png.display());
    eprintln!("wave_coherent: wrote {}", observed_color_png.display());
    eprintln!("wave_coherent: wrote {}", instrument_png.display());
    eprintln!("wave_coherent: wrote {}", radial_csv.display());
    eprintln!("wave_coherent: wrote {}", summary_path.display());
}

fn run_sgr_astar_eht_report(out_dir: &Path, width: usize, height: usize, blur_sigma: f64) {
    use std::f64::consts::PI;
    use std::io::Write as _;

    let Some(view) = VIEWS.iter().find(|v| v.slug == "sgr_astar") else {
        eprintln!("sgr_astar view not found");
        return;
    };

    let prev_spec = std::env::var("BH_SPECTRUM").ok();
    let prev_disk = std::env::var("BH_DISK_MODEL").ok();
    let prev_tx = std::env::var("BH_USE_TRANSFER").ok();
    let prev_tau = std::env::var("BH_TAU_SCALE").ok();
    let prev_suffix = std::env::var("BH_SLUG_SUFFIX").ok();
    let plasma_model = PlasmaModel::from_env();

    std::env::set_var("BH_SPECTRUM", "millimeter");
    std::env::set_var("BH_DISK_MODEL", "riaf");
    std::env::set_var("BH_USE_TRANSFER", "1");
    std::env::set_var("BH_TAU_SCALE", "1.0");

    let psf_alpha = if blur_sigma > 0.0 { blur_sigma } else { 1.2 };

    std::env::set_var("BH_SLUG_SUFFIX", "intrinsic_230ghz");
    let intrinsic_png = render_view(out_dir, view, width, height, 0.0);
    std::env::set_var("BH_SLUG_SUFFIX", "eht_230ghz");
    let observed_png = render_view(out_dir, view, width, height, 0.0);

    let metric = GutoeMetric::planck_units(1.0);
    let cfg = RenderConfig {
        width,
        height,
        fov_rs: view.fov,
        inclination_deg: view.inc,
        max_phi: view.max_phi_pi * PI,
        dphi: view.dphi,
    };
    let raw = render_with_options(
        &metric,
        None,
        view.disk_inner,
        view.disk_outer,
        &cfg,
        view.az,
        view.doppler,
        view.ring_mode,
        view.interior_mode,
        view.core_look_mode,
        SpectralBand::Millimeter,
        DiskModel::Riaf,
        plasma_model,
        true,
        1.0,
        true,
        0.0,
        0.5,
        0.5,
    );
    let m_intr = compute_kerr_image_metrics(&raw, width, height);
    let m_obs = compute_kerr_image_metrics(&raw, width, height);

    // Non-Gaussian observational proxy: Moffat PSF + Richardson-Lucy reconstruction.
    let radius = ((3.0 * psf_alpha).ceil() as usize).clamp(2, 8);
    let k = moffat_kernel(psf_alpha, 2.6, radius);
    let raw_f: Vec<[f64; 3]> = raw
        .iter()
        .map(|p| [p[0] as f64 / 255.0, p[1] as f64 / 255.0, p[2] as f64 / 255.0])
        .collect();
    let dirty_f = convolve_rgb_f64(&raw_f, width, height, &k, radius);
    let dirty: Vec<[u8; 3]> = dirty_f
        .into_iter()
        .map(|p| {
            [
                (p[0] * 255.0).clamp(0.0, 255.0) as u8,
                (p[1] * 255.0).clamp(0.0, 255.0) as u8,
                (p[2] * 255.0).clamp(0.0, 255.0) as u8,
            ]
        })
        .collect();
    let recon = richardson_lucy_rgb(&dirty, width, height, &k, radius, 8);
    let dirty_png = out_dir.join("sgr_astar__eht_dirty_230ghz.png");
    let recon_png = out_dir.join("sgr_astar__eht_recon_230ghz.png");
    save_png_rgb(&dirty_png, &dirty, width, height);
    save_png_rgb(&recon_png, &recon, width, height);
    let m_dirty = compute_kerr_image_metrics(&dirty, width, height);
    let m_recon = compute_kerr_image_metrics(&recon, width, height);

    let csv_path = out_dir.join("sgr_astar_eht_metrics.csv");
    let mut f = std::fs::File::create(&csv_path).expect("create sgr_astar_eht_metrics.csv");
    writeln!(
        f,
        "mode,width,height,beam_sigma_px,shadow_diameter_px,ring_thickness_px,flux_asymmetry,closure_phase_deg"
    )
    .expect("write csv header");
    writeln!(
        f,
        "intrinsic_230ghz,{width},{height},0.0,{:.6},{:.6},{:.9},{:.6}",
        m_intr.shadow_diameter, m_intr.ring_thickness, m_intr.flux_asymmetry, m_intr.closure_phase_deg
    )
    .expect("write intrinsic row");
    writeln!(
        f,
        "eht_observed_230ghz_nobeam,{width},{height},0.000,{:.6},{:.6},{:.9},{:.6}",
        m_obs.shadow_diameter, m_obs.ring_thickness, m_obs.flux_asymmetry, m_obs.closure_phase_deg
    )
    .expect("write observed row");
    writeln!(
        f,
        "eht_dirty_moffat,{width},{height},{psf_alpha:.3},{:.6},{:.6},{:.9},{:.6}",
        m_dirty.shadow_diameter, m_dirty.ring_thickness, m_dirty.flux_asymmetry, m_dirty.closure_phase_deg
    )
    .expect("write dirty row");
    writeln!(
        f,
        "eht_recon_rl,{width},{height},{psf_alpha:.3},{:.6},{:.6},{:.9},{:.6}",
        m_recon.shadow_diameter, m_recon.ring_thickness, m_recon.flux_asymmetry, m_recon.closure_phase_deg
    )
    .expect("write recon row");

    eprintln!("\nSgr A* EHT report:");
    eprintln!("  intrinsic: {}", intrinsic_png.display());
    eprintln!("  observed : {}", observed_png.display());
    eprintln!("  dirty    : {}", dirty_png.display());
    eprintln!("  recon    : {}", recon_png.display());
    eprintln!("  metrics  : {}", csv_path.display());

    if let Some(v) = prev_spec {
        std::env::set_var("BH_SPECTRUM", v);
    } else {
        std::env::remove_var("BH_SPECTRUM");
    }
    if let Some(v) = prev_disk {
        std::env::set_var("BH_DISK_MODEL", v);
    } else {
        std::env::remove_var("BH_DISK_MODEL");
    }
    if let Some(v) = prev_tx {
        std::env::set_var("BH_USE_TRANSFER", v);
    } else {
        std::env::remove_var("BH_USE_TRANSFER");
    }
    if let Some(v) = prev_tau {
        std::env::set_var("BH_TAU_SCALE", v);
    } else {
        std::env::remove_var("BH_TAU_SCALE");
    }
    if let Some(v) = prev_suffix {
        std::env::set_var("BH_SLUG_SUFFIX", v);
    } else {
        std::env::remove_var("BH_SLUG_SUFFIX");
    }
}

fn run_eht_uv_export(view: &View, width: usize, height: usize, out_dir: &Path) {
    use std::collections::HashMap;
    use std::f64::consts::PI;
    use std::io::Write as _;

    let parse_env_f64 = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<f64>().ok());
    let parse_env_bool = |k: &str| {
        std::env::var(k).ok().is_some_and(|s| {
            matches!(s.as_str(), "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON")
        })
    };
    fn hash_u64(s: &str) -> u64 {
        let mut h: u64 = 0xcbf29ce484222325;
        for &b in s.as_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
        h
    }
    fn u01(seed: u64) -> f64 {
        let x = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let v = ((x >> 11) as f64) / ((1u64 << 53) as f64);
        v.clamp(1e-12, 1.0 - 1e-12)
    }
    fn gaussian(seed: u64) -> f64 {
        let u1 = u01(seed);
        let u2 = u01(seed ^ 0x9E3779B97F4A7C15);
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    fn cmul(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
        (ar * br - ai * bi, ar * bi + ai * br)
    }

    let force_gr = parse_env_bool("BH_FORCE_GR");
    let metric = if view.gr_mode || force_gr {
        GutoeMetric::schwarzschild(1.0)
    } else {
        GutoeMetric::planck_units(1.0)
    };
    let kerr_metric = parse_env_f64("BH_KERR_ASTAR").and_then(|a| KerrMetric::new(metric.r_s, a));
    let fov_rs = parse_env_f64("BH_FOV_OVERRIDE").unwrap_or(view.fov);
    let inc_deg = parse_env_f64("BH_INC_OVERRIDE").unwrap_or(view.inc);
    let az_deg = parse_env_f64("BH_AZ_OVERRIDE").unwrap_or(view.az);
    let max_phi_pi = parse_env_f64("BH_MAX_PHI_PI_OVERRIDE").unwrap_or(view.max_phi_pi);
    let dphi = parse_env_f64("BH_DPHI_OVERRIDE").unwrap_or(view.dphi);
    let cfg = RenderConfig {
        width,
        height,
        fov_rs,
        inclination_deg: inc_deg,
        max_phi: max_phi_pi * PI,
        dphi,
    };
    let spectral_band = SpectralBand::from_env();
    let disk_model = DiskModel::from_env();
    let plasma_model = PlasmaModel::from_env();
    let use_transfer = parse_env_bool("BH_USE_TRANSFER");
    let tau_scale = parse_env_f64("BH_TAU_SCALE").unwrap_or(1.0).max(0.0);

    let raw = render_with_options(
        &metric,
        kerr_metric.as_ref(),
        view.disk_inner,
        view.disk_outer,
        &cfg,
        az_deg,
        view.doppler,
        view.ring_mode,
        view.interior_mode,
        view.core_look_mode,
        spectral_band,
        disk_model,
        plasma_model,
        use_transfer,
        tau_scale,
        true,
        0.0,
        0.5,
        0.5,
    );

    // Baselines in normalized uv units (cycles per FOV).
    // Optional input file for EHT-17:
    //   BH_UV_TRACK_CSV=/path/file.csv
    //   CSV forms supported (header optional):
    //     baseline,u,v
    //     station_i,station_j,u,v
    //     baseline,station_i,station_j,u,v
    let mut baselines: Vec<(String, f64, f64)> = Vec::new();
    if let Ok(path) = std::env::var("BH_UV_TRACK_CSV") {
        if let Ok(text) = std::fs::read_to_string(&path) {
            for line in text.lines() {
                let l = line.trim();
                if l.is_empty() || l.starts_with('#') {
                    continue;
                }
                let parts: Vec<&str> = l.split(',').map(|x| x.trim()).collect();
                if parts.len() < 3 {
                    continue;
                }
                if parts[0].eq_ignore_ascii_case("baseline")
                    || parts[0].eq_ignore_ascii_case("station_i")
                    || parts[0].eq_ignore_ascii_case("station1")
                    || parts[1].eq_ignore_ascii_case("u")
                    || parts[2].eq_ignore_ascii_case("v")
                {
                    continue;
                }
                let parsed = if parts.len() >= 5 {
                    let name = if parts[0].is_empty() {
                        format!("{}-{}", parts[1], parts[2])
                    } else {
                        parts[0].to_string()
                    };
                    let Ok(u) = parts[3].parse::<f64>() else {
                        continue;
                    };
                    let Ok(v) = parts[4].parse::<f64>() else {
                        continue;
                    };
                    Some((name, u, v))
                } else if parts.len() >= 4 && parts[1].parse::<f64>().is_err() {
                    let name = format!("{}-{}", parts[0], parts[1]);
                    let Ok(u) = parts[2].parse::<f64>() else {
                        continue;
                    };
                    let Ok(v) = parts[3].parse::<f64>() else {
                        continue;
                    };
                    Some((name, u, v))
                } else {
                    let Ok(u) = parts[1].parse::<f64>() else {
                        continue;
                    };
                    let Ok(v) = parts[2].parse::<f64>() else {
                        continue;
                    };
                    Some((parts[0].to_string(), u, v))
                };
                if let Some(x) = parsed {
                    baselines.push(x);
                }
            }
        }
    }
    if baselines.is_empty() {
        baselines = vec![
            ("ALMA-LMT".into(), 8.0, 2.0),
            ("ALMA-SMT".into(), 7.2, -2.2),
            ("ALMA-PV".into(), 6.1, 4.4),
            ("ALMA-SPT".into(), 0.8, -9.5),
            ("LMT-SMT".into(), 2.7, -1.6),
            ("LMT-PV".into(), 1.1, 5.4),
            ("LMT-JCMT".into(), -5.3, 3.6),
            ("SMT-PV".into(), -1.9, 6.1),
            ("SMT-JCMT".into(), -7.1, 1.1),
            ("SMT-SMA".into(), -6.8, 0.7),
            ("PV-JCMT".into(), -8.4, -2.2),
            ("PV-SMA".into(), -8.1, -2.5),
            ("JCMT-SMA".into(), 0.4, -0.3),
            ("APEX-LMT".into(), 7.8, 1.8),
            ("APEX-PV".into(), 5.8, 4.1),
            ("APEX-SPT".into(), 0.5, -9.1),
        ];
    }

    let uv_csv = out_dir.join(format!("{}_eht_uv.csv", view.slug));
    let mut f = std::fs::File::create(&uv_csv).expect("create eht_uv csv");
    writeln!(f, "baseline,u_norm,v_norm,re,im,amp,phase_deg").expect("write eht_uv header");
    let mut vis = Vec::with_capacity(baselines.len());
    for (name, u, v) in baselines {
        let (re, im) = visibility_dft(&raw, width, height, u, v);
        let amp = (re * re + im * im).sqrt();
        let phase = im.atan2(re).to_degrees();
        writeln!(f, "{name},{u:.6},{v:.6},{re:.9},{im:.9},{amp:.9},{phase:.6}")
            .expect("write eht_uv row");
        vis.push((name, re, im));
    }

    // Optional station gain + thermal noise corruption model (EHT-18/EHT-19 scaffold).
    let noise_sigma = parse_env_f64("BH_EHT_NOISE_SIGMA").unwrap_or(0.0).max(0.0);
    let gain_amp_sigma = parse_env_f64("BH_EHT_GAIN_AMP_SIGMA").unwrap_or(0.0).max(0.0);
    let gain_phase_sigma_deg = parse_env_f64("BH_EHT_GAIN_PHASE_SIGMA_DEG")
        .unwrap_or(0.0)
        .max(0.0);
    let bw_hz = parse_env_f64("BH_EHT_BW_HZ").unwrap_or(2.0e9).max(1.0);
    let tint_s = parse_env_f64("BH_EHT_TINT_SEC").unwrap_or(10.0).max(1e-6);
    let sefd_csv_path = std::env::var("BH_EHT_SEFD_CSV").ok();
    let gain_csv_path = std::env::var("BH_EHT_GAIN_CSV").ok();
    let mut sefd_map: HashMap<String, f64> = HashMap::new();
    if let Some(path) = &sefd_csv_path {
        if let Ok(text) = std::fs::read_to_string(path) {
            for line in text.lines() {
                let l = line.trim();
                if l.is_empty() || l.starts_with('#') {
                    continue;
                }
                let parts: Vec<&str> = l.split(',').map(|x| x.trim()).collect();
                if parts.len() < 2 {
                    continue;
                }
                if parts[0].eq_ignore_ascii_case("station") || parts[1].eq_ignore_ascii_case("sefd_jy") {
                    continue;
                }
                if let Ok(sefd) = parts[1].parse::<f64>() {
                    sefd_map.insert(parts[0].to_string(), sefd.max(0.0));
                }
            }
        }
    }
    let mut gain_map: HashMap<String, (f64, f64)> = HashMap::new();
    if let Some(path) = &gain_csv_path {
        if let Ok(text) = std::fs::read_to_string(path) {
            for line in text.lines() {
                let l = line.trim();
                if l.is_empty() || l.starts_with('#') {
                    continue;
                }
                let parts: Vec<&str> = l.split(',').map(|x| x.trim()).collect();
                if parts.len() < 3 {
                    continue;
                }
                if parts[0].eq_ignore_ascii_case("station")
                    || parts[1].eq_ignore_ascii_case("amp_gain")
                    || parts[2].eq_ignore_ascii_case("phase_deg")
                {
                    continue;
                }
                if let (Ok(amp), Ok(phase_deg)) = (parts[1].parse::<f64>(), parts[2].parse::<f64>()) {
                    gain_map.insert(parts[0].to_string(), (amp, phase_deg.to_radians()));
                }
            }
        }
    }
    let apply_corruption = noise_sigma > 0.0
        || gain_amp_sigma > 0.0
        || gain_phase_sigma_deg > 0.0
        || !sefd_map.is_empty()
        || !gain_map.is_empty();
    let mut vis_corrupt: Vec<(String, f64, f64)> = Vec::with_capacity(vis.len());
    let uv_corrupt_csv = out_dir.join(format!("{}_eht_uv_corrupted.csv", view.slug));
    if apply_corruption {
        let mut fc = std::fs::File::create(&uv_corrupt_csv).expect("create eht_uv_corrupted csv");
        writeln!(
            fc,
            "baseline,re,im,amp,phase_deg,noise_sigma,sefd_i,sefd_j,bw_hz,tint_s,gain_amp_sigma,gain_phase_sigma_deg"
        )
        .expect("write eht_uv_corrupted header");
        for (name, re, im) in &vis {
            let mut cr = *re;
            let mut ci = *im;
            let mut parts = name.split('-');
            let s1 = parts.next().unwrap_or("A");
            let s2 = parts.next().unwrap_or("B");
            let (a1, p1) = if let Some((amp, ph)) = gain_map.get(s1) {
                (*amp, *ph)
            } else {
                (
                    1.0 + gain_amp_sigma * gaussian(hash_u64(&format!("amp1:{s1}"))),
                    gain_phase_sigma_deg.to_radians() * gaussian(hash_u64(&format!("ph1:{s1}"))),
                )
            };
            let (a2, p2) = if let Some((amp, ph)) = gain_map.get(s2) {
                (*amp, *ph)
            } else {
                (
                    1.0 + gain_amp_sigma * gaussian(hash_u64(&format!("amp2:{s2}"))),
                    gain_phase_sigma_deg.to_radians() * gaussian(hash_u64(&format!("ph2:{s2}"))),
                )
            };
            let gamp = (a1 * a2).max(0.01);
            let gph = p1 - p2;
            let (gr, gi) = (gamp * gph.cos(), gamp * gph.sin());
            (cr, ci) = cmul(cr, ci, gr, gi);
            let sefd_i = sefd_map.get(s1).copied().unwrap_or(0.0);
            let sefd_j = sefd_map.get(s2).copied().unwrap_or(0.0);
            let baseline_noise_sigma = if sefd_i > 0.0 && sefd_j > 0.0 {
                ((sefd_i * sefd_j) / (2.0 * bw_hz * tint_s)).sqrt()
            } else {
                noise_sigma
            };
            cr += baseline_noise_sigma * gaussian(hash_u64(&format!("nr:{name}")));
            ci += baseline_noise_sigma * gaussian(hash_u64(&format!("ni:{name}")));
            let amp = (cr * cr + ci * ci).sqrt();
            let phase = ci.atan2(cr).to_degrees();
            writeln!(
                fc,
                "{name},{cr:.9},{ci:.9},{amp:.9},{phase:.6},{baseline_noise_sigma:.6},{sefd_i:.3},{sefd_j:.3},{bw_hz:.3},{tint_s:.3},{gain_amp_sigma:.6},{gain_phase_sigma_deg:.6}"
            )
            .expect("write eht_uv_corrupted row");
            vis_corrupt.push((name.clone(), cr, ci));
        }
    } else {
        vis_corrupt = vis.iter().map(|(n, r, i)| (n.clone(), *r, *i)).collect();
    }

    let closure_csv = out_dir.join(format!("{}_eht_closure.csv", view.slug));
    let mut c = std::fs::File::create(&closure_csv).expect("create eht_closure csv");
    writeln!(c, "triangle,closure_phase_deg,closure_phase_corrupted_deg")
        .expect("write closure header");
    let vis_map: HashMap<String, (f64, f64)> = vis
        .iter()
        .map(|(name, re, im)| (name.clone(), (*re, *im)))
        .collect();
    let vis_corrupt_map: HashMap<String, (f64, f64)> = vis_corrupt
        .iter()
        .map(|(name, re, im)| (name.clone(), (*re, *im)))
        .collect();
    let get_pair = |a: &str, b: &str, map: &HashMap<String, (f64, f64)>| -> Option<(f64, f64)> {
        let k1 = format!("{a}-{b}");
        if let Some((re, im)) = map.get(&k1) {
            return Some((*re, *im));
        }
        let k2 = format!("{b}-{a}");
        map.get(&k2).map(|(re, im)| (*re, -*im))
    };
    let mut stations: Vec<String> = vis
        .iter()
        .flat_map(|(name, _, _)| {
            let mut parts = name.split('-');
            let a = parts.next().unwrap_or("").to_string();
            let b = parts.next().unwrap_or("").to_string();
            [a, b]
        })
        .filter(|s| !s.is_empty())
        .collect();
    stations.sort();
    stations.dedup();
    for i in 0..stations.len() {
        for j in (i + 1)..stations.len() {
            for k in (j + 1)..stations.len() {
                let sa = &stations[i];
                let sb = &stations[j];
                let sc = &stations[k];
                let Some((ar, ai)) = get_pair(sa, sb, &vis_map) else {
                    continue;
                };
                let Some((br, bi)) = get_pair(sb, sc, &vis_map) else {
                    continue;
                };
                let Some((cr, ci)) = get_pair(sc, sa, &vis_map) else {
                    continue;
                };
                let Some((ar2, ai2)) = get_pair(sa, sb, &vis_corrupt_map) else {
                    continue;
                };
                let Some((br2, bi2)) = get_pair(sb, sc, &vis_corrupt_map) else {
                    continue;
                };
                let Some((cr2, ci2)) = get_pair(sc, sa, &vis_corrupt_map) else {
                    continue;
                };
                let name = format!("{sa}-{sb}-{sc}");
                let p1r = ar * br - ai * bi;
                let p1i = ar * bi + ai * br;
                let pr = p1r * cr - p1i * ci;
                let pi = p1r * ci + p1i * cr;
                let closure = pi.atan2(pr).to_degrees();

                let p2r = ar2 * br2 - ai2 * bi2;
                let p2i = ar2 * bi2 + ai2 * br2;
                let pr2 = p2r * cr2 - p2i * ci2;
                let pi2 = p2r * ci2 + p2i * cr2;
                let closure_corrupt = pi2.atan2(pr2).to_degrees();
                writeln!(c, "{name},{closure:.6},{closure_corrupt:.6}").expect("write closure row");
            }
        }
    }

    let closure_amp_csv = out_dir.join(format!("{}_eht_closure_amp.csv", view.slug));
    let mut ca = std::fs::File::create(&closure_amp_csv).expect("create eht_closure_amp csv");
    writeln!(
        ca,
        "quad,closure_amplitude,closure_amplitude_corrupted"
    )
    .expect("write closure amp header");
    let amp = |z: (f64, f64)| -> f64 { (z.0 * z.0 + z.1 * z.1).sqrt().max(1e-12) };
    for i in 0..stations.len() {
        for j in (i + 1)..stations.len() {
            for k in (j + 1)..stations.len() {
                for l in (k + 1)..stations.len() {
                    let sa = &stations[i];
                    let sb = &stations[j];
                    let sc = &stations[k];
                    let sd = &stations[l];
                    // A = |Vab||Vcd| / (|Vac||Vbd|)
                    let Some(vab) = get_pair(sa, sb, &vis_map) else { continue };
                    let Some(vcd) = get_pair(sc, sd, &vis_map) else { continue };
                    let Some(vac) = get_pair(sa, sc, &vis_map) else { continue };
                    let Some(vbd) = get_pair(sb, sd, &vis_map) else { continue };
                    let Some(vab2) = get_pair(sa, sb, &vis_corrupt_map) else { continue };
                    let Some(vcd2) = get_pair(sc, sd, &vis_corrupt_map) else { continue };
                    let Some(vac2) = get_pair(sa, sc, &vis_corrupt_map) else { continue };
                    let Some(vbd2) = get_pair(sb, sd, &vis_corrupt_map) else { continue };
                    let a = (amp(vab) * amp(vcd)) / (amp(vac) * amp(vbd));
                    let a2 = (amp(vab2) * amp(vcd2)) / (amp(vac2) * amp(vbd2));
                    writeln!(ca, "{sa}-{sb}-{sc}-{sd},{a:.9},{a2:.9}")
                        .expect("write closure amp row");
                }
            }
        }
    }

    let summary = out_dir.join(format!("{}_eht_uv.json", view.slug));
    let mut s = std::fs::File::create(&summary).expect("create eht_uv json");
    let m_solar = parse_env_f64("BH_MASS_SOLAR").unwrap_or(6.5e9).max(1.0);
    let d_mpc = parse_env_f64("BH_DISTANCE_MPC").unwrap_or(16.8).max(1e-9);
    let g_si = 6.67430e-11_f64;
    let c_si = 299_792_458.0_f64;
    let m_sun_kg = 1.98847e30_f64;
    let mpc_m = 3.0856775814913673e22_f64;
    let rad_to_uas = 206_264_806_247.09636_f64;
    let rg_m = g_si * m_solar * m_sun_kg / (c_si * c_si);
    let rs_rad = (2.0 * rg_m) / (d_mpc * mpc_m);
    let rs_uas = rs_rad * rad_to_uas;
    let fov_uas = fov_rs * rs_uas;
    let pixel_uas = fov_uas / (width as f64).max(1.0);
    let shadow_diam_gr_uas = 2.0 * (27.0_f64).sqrt() * rg_m / (d_mpc * mpc_m) * rad_to_uas;
    writeln!(
        s,
        "{{\n  \"slug\": \"{}\",\n  \"png\": \"{}\",\n  \"width\": {},\n  \"height\": {},\n  \"fov_rs\": {:.6},\n  \"uv_units\": \"normalized cycles per image FOV\",\n  \"note\": \"Synthetic uv set; replace with 2017 measured tracks in EHT-17\",\n  \"uv_csv\": \"{}\",\n  \"uv_corrupted_csv\": \"{}\",\n  \"closure_csv\": \"{}\",\n  \"closure_amp_csv\": \"{}\",\n  \"noise_sigma\": {:.6},\n  \"gain_amp_sigma\": {:.6},\n  \"gain_phase_sigma_deg\": {:.6},\n  \"bw_hz\": {:.3},\n  \"tint_s\": {:.3},\n  \"sefd_csv\": \"{}\",\n  \"gain_csv\": \"{}\",\n  \"m_solar\": {:.3},\n  \"distance_mpc\": {:.6},\n  \"rs_uas\": {:.6},\n  \"fov_uas\": {:.6},\n  \"pixel_uas\": {:.6},\n  \"shadow_diameter_gr_uas\": {:.6}\n}}",
        view.slug,
        out_dir.join(format!("{}.png", view.slug)).display(),
        width,
        height,
        fov_rs,
        uv_csv.display(),
        uv_corrupt_csv.display(),
        closure_csv.display(),
        closure_amp_csv.display(),
        noise_sigma,
        gain_amp_sigma,
        gain_phase_sigma_deg,
        bw_hz,
        tint_s,
        sefd_csv_path.as_deref().unwrap_or(""),
        gain_csv_path.as_deref().unwrap_or(""),
        m_solar,
        d_mpc,
        rs_uas,
        fov_uas,
        pixel_uas,
        shadow_diam_gr_uas
    )
    .expect("write eht_uv json");

    eprintln!("eht_uv: wrote {}", uv_csv.display());
    if apply_corruption {
        eprintln!("eht_uv: wrote {}", uv_corrupt_csv.display());
    }
    eprintln!("eht_uv: wrote {}", closure_csv.display());
    eprintln!("eht_uv: wrote {}", closure_amp_csv.display());
    eprintln!("eht_uv: wrote {}", summary.display());
}

fn run_metric_report(out_dir: &Path) {
    use std::io::Write as _;

    let r_s = std::env::var("BH_RS")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.0)
        .max(1e-9);
    let l_p = std::env::var("BH_LP")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.0)
        .max(0.0);
    let r_max = std::env::var("BH_METRIC_R_MAX")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(8.0)
        .max(0.2);
    let samples = std::env::var("BH_METRIC_SAMPLES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(64)
        .clamp(8, 4096);

    let metric = if l_p == 0.0 {
        GutoeMetric::schwarzschild(r_s)
    } else {
        GutoeMetric::new(r_s, l_p)
    };
    let rc = metric.r_core();
    let rh = metric.r_horizon();
    let rph = metric.r_photon_sphere();
    let risco = metric.r_isco();
    let th = metric.hawking_temperature();
    let th_gr = metric.gr_hawking_temperature();
    let frac = metric.hawking_correction_fraction();

    let csv_path = out_dir.join("metric_report.csv");
    let mut f = std::fs::File::create(&csv_path).expect("create metric_report.csv");
    writeln!(f, "# GUTOE metric report").expect("write metric report header");
    writeln!(f, "# C_INF,{C_INF:.10}").expect("write metric report header");
    writeln!(f, "# LAMBDA_QG,{LAMBDA_QG:.10}").expect("write metric report header");
    writeln!(f, "# WATSON_SC,{WATSON_SC:.10}").expect("write metric report header");
    writeln!(f, "# r_s,{r_s:.10}").expect("write metric report header");
    writeln!(f, "# l_p,{l_p:.10}").expect("write metric report header");
    writeln!(f, "# r_core,{rc:.10}").expect("write metric report header");
    writeln!(f, "# r_horizon,{}", rh.map(|v| format!("{v:.10}")).unwrap_or_else(|| "none".into()))
        .expect("write metric report header");
    writeln!(
        f,
        "# r_photon_sphere,{}",
        rph.map(|v| format!("{v:.10}")).unwrap_or_else(|| "none".into())
    )
    .expect("write metric report header");
    writeln!(f, "# r_isco,{}", risco.map(|v| format!("{v:.10}")).unwrap_or_else(|| "none".into()))
        .expect("write metric report header");
    writeln!(f, "# T_hawking,{th:.12}").expect("write metric report header");
    writeln!(f, "# T_gr,{th_gr:.12}").expect("write metric report header");
    writeln!(f, "# dT_over_T,{frac:.12}").expect("write metric report header");
    writeln!(f, "r,reff,g_tt,g_rr,g_theta,g_phi_eq").expect("write metric report columns");

    for i in 0..samples {
        let t = i as f64 / (samples.saturating_sub(1).max(1) as f64);
        let r = 1e-6 + t * r_max;
        let re = metric.r_eff(r);
        let gtt = metric.g_tt(r);
        let grr = metric.g_rr(r);
        let gth = metric.g_theta(r);
        let gph = metric.g_phi(r, std::f64::consts::FRAC_PI_2);
        writeln!(f, "{r:.10},{re:.10},{gtt:.12},{grr:.12},{gth:.12},{gph:.12}")
            .expect("write metric report row");
    }

    eprintln!("\nMetric report:");
    eprintln!("  r_s={r_s:.6} l_p={l_p:.6} r_core={rc:.6}");
    eprintln!(
        "  horizon={} photon_sphere={} isco={}",
        rh.map(|v| format!("{v:.6}")).unwrap_or_else(|| "none".into()),
        rph.map(|v| format!("{v:.6}")).unwrap_or_else(|| "none".into()),
        risco.map(|v| format!("{v:.6}")).unwrap_or_else(|| "none".into())
    );
    eprintln!("  T_h={th:.8e} T_gr={th_gr:.8e} dT/T={frac:.8e}");
    eprintln!("  csv: {}", csv_path.display());
}

// ── HTML gallery ──────────────────────────────────────────────────────────────

fn build_html() -> String {
    let mut sections = String::new();

    // Group views into labelled sections
    let groups: &[(&str, &[usize])] = &[
        ("Classic Views — GUTOE Schwarzschild", &[0, 1, 2, 3, 4, 5]),
        ("EHT Comparisons — Relativistic Doppler", &[6, 7]),
        ("Photon Sub-Ring Structure", &[8]),
        ("GR Comparison — Schwarzschild (l_P = 0)", &[9]),
        ("Inside the Shadow — Fractal Impact Map", &[10]),
        ("Inside the Horizon — GUTOE Lattice Floor", &[11, 12]),
    ];

    for &(section_title, indices) in groups {
        let mut imgs = String::new();
        for &i in indices {
            let v = &VIEWS[i];
            imgs.push_str(&format!(
                r#"<figure>
  <img src="/img/{slug}.png" alt="{label}" loading="lazy">
  <figcaption><strong>{label}</strong><br><span class="cap">{caption}</span></figcaption>
</figure>
"#,
                slug = v.slug,
                label = v.label,
                caption = v.caption,
            ));
        }
        sections.push_str(&format!(
            r#"<section>
<h2>{title}</h2>
<div class="gallery">{imgs}</div>
</section>
"#,
            title = section_title,
            imgs = imgs,
        ));
    }

    format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>GUTOE Black Hole Gallery</title>
<style>
  :root {{ --bg:#050508; --fg:#ddd; --accent:#ffd080; --dim:#888; --rim:#333; }}
  * {{ box-sizing:border-box; }}
  body {{ background:var(--bg); color:var(--fg); font-family:monospace;
          padding:2rem; max-width:1800px; margin:0 auto; }}
  h1   {{ color:var(--accent); margin-bottom:.2rem; font-size:1.6rem; }}
  h2   {{ color:var(--accent); border-bottom:1px solid var(--rim);
          padding-bottom:.4rem; margin-top:2rem; font-size:1.1rem; }}
  .sub {{ color:#aaa; margin-top:0; font-size:.85rem; }}
  .gallery {{ display:flex; flex-wrap:wrap; gap:1.2rem; margin-top:1rem; }}
  figure {{ margin:0; }}
  img  {{ display:block; width:400px; height:400px; border:1px solid var(--rim);
           image-rendering:auto; }}
  figcaption {{ width:400px; font-size:.72rem; color:var(--dim); margin-top:.4rem; line-height:1.4; }}
  figcaption strong {{ color:#ccc; display:block; margin-bottom:.2rem; }}
  .cap {{ display:block; }}
  section {{ margin-bottom:2.5rem; }}
</style>
</head>
<body>
<h1>GUTOE Black Hole Gallery</h1>
<p class="sub">
  Cl(1,3) Schwarzschild metric — singularity regularised by SC lattice core r_c = √C_∞ · l_P
  &nbsp;|&nbsp; C_∞ = 0.5466 (GPU Richardson L=161–961)
  &nbsp;|&nbsp; CUDA geodesic ray tracer, RK4
  &nbsp;|&nbsp; Views 1–10: exterior camera at infinity &nbsp;|&nbsp; Views 11–13: interior-horizon perspectives
</p>
{sections}
</body>
</html>
"#,
        sections = sections,
    )
}

// ── Minimal HTTP server ───────────────────────────────────────────────────────

fn serve_http(out_dir: Arc<PathBuf>, html: Arc<String>) {
    let addr = "0.0.0.0:52345";
    let listener = match TcpListener::bind(addr) {
        Ok(l) => l,
        Err(e) => {
            eprintln!("\nGallery server skipped ({} busy): {e}", addr);
            eprintln!("Images are still saved under /tmp/bh_renders\n");
            return;
        }
    };
    eprintln!("\nGallery live at  http://10.7.1.200:52345/");
    eprintln!("Press Ctrl-C to stop.\n");

    for stream in listener.incoming() {
        match stream {
            Ok(s) => handle_connection(s, &out_dir, &html),
            Err(e) => eprintln!("accept error: {e}"),
        }
    }
}

fn handle_connection(mut stream: TcpStream, out_dir: &Path, html: &str) {
    let mut buf = [0u8; 4096];
    let n = match stream.read(&mut buf) {
        Ok(n) => n,
        Err(_) => return,
    };
    let req = std::str::from_utf8(&buf[..n]).unwrap_or("");
    let path = req
        .lines()
        .next()
        .and_then(|l| l.split_whitespace().nth(1))
        .unwrap_or("/");

    if path == "/" || path == "/index.html" {
        let body = html.as_bytes();
        let _ = write!(
            stream,
            "HTTP/1.1 200 OK\r\nContent-Type: text/html; charset=utf-8\r\n\
             Content-Length: {}\r\nConnection: close\r\n\r\n",
            body.len()
        );
        let _ = stream.write_all(body);
        return;
    }

    if let Some(slug) = path
        .strip_prefix("/img/")
        .and_then(|s| s.strip_suffix(".png"))
    {
        let file_path = out_dir.join(format!("{slug}.png"));
        match fs::read(&file_path) {
            Ok(data) => {
                let _ = write!(
                    stream,
                    "HTTP/1.1 200 OK\r\nContent-Type: image/png\r\n\
                     Content-Length: {}\r\nConnection: close\r\n\r\n",
                    data.len()
                );
                let _ = stream.write_all(&data);
            }
            Err(_) => {
                let _ = stream
                    .write_all(b"HTTP/1.1 404 Not Found\r\nContent-Length: 9\r\n\r\nNot Found");
            }
        }
        return;
    }

    let _ = stream.write_all(b"HTTP/1.1 404 Not Found\r\nContent-Length: 9\r\n\r\nNot Found");
}

// ── Entry point ───────────────────────────────────────────────────────────────

/// Parse "WxH" or "W" (or "WxH") into (width, height).  Returns (0,0) on failure.
fn parse_resolution(s: &str) -> (usize, usize) {
    if let Some((w, h)) = s.split_once('x') {
        let w = w.parse::<usize>().unwrap_or(0);
        let h = h.parse::<usize>().unwrap_or(0);
        (w, h)
    } else {
        let w = s.parse::<usize>().unwrap_or(0);
        (w, 0) // 0 height → square
    }
}

/// Returns true if a string looks like a resolution spec ("WxH" or plain integer).
fn looks_like_resolution(s: &str) -> bool {
    if let Some((w, h)) = s.split_once('x') {
        w.parse::<usize>().is_ok() && h.parse::<usize>().is_ok()
    } else {
        s.parse::<usize>().is_ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn downsample_box_2x_preserves_constant_color() {
        let hi = vec![[12, 34, 56]; 4 * 4];
        let lo = downsample_box(&hi, 4, 4, 2);
        assert_eq!(lo.len(), 4);
        assert!(lo.iter().all(|&px| px == [12, 34, 56]));
    }

    #[test]
    fn parse_resolution_and_detector() {
        assert_eq!(parse_resolution("3840x2160"), (3840, 2160));
        assert_eq!(parse_resolution("1200"), (1200, 0));
        assert!(looks_like_resolution("512x512"));
        assert!(looks_like_resolution("768"));
        assert!(!looks_like_resolution("camera_core"));
    }

    #[test]
    fn sample_jitter_bounds() {
        for i in 0..16 {
            let (jx, jy) = sample_jitter(i, 16);
            assert!((0.0..1.0).contains(&jx));
            assert!((0.0..1.0).contains(&jy));
        }
        assert_eq!(sample_jitter(0, 1), (0.5, 0.5));
    }

    #[test]
    fn transfer_factor_has_doppler_asymmetry() {
        let r = 3.0;
        let r_s = 1.0;
        let sin_inc = 0.9;
        let approaching = disk_transfer_factor(
            r,
            r_s,
            0.8,
            std::f64::consts::FRAC_PI_2,
            sin_inc,
            true,
            0.0,
        );
        let receding = disk_transfer_factor(
            r,
            r_s,
            -0.8,
            -std::f64::consts::FRAC_PI_2,
            sin_inc,
            true,
            0.0,
        );
        assert!(approaching > receding);
    }

    #[test]
    fn adaptive_dphi_tightens_near_critical() {
        let base = 0.005;
        let bcrit = 2.598_076_211;
        let near = adaptive_dphi_for_b(base, 1.01 * bcrit, bcrit);
        let far = adaptive_dphi_for_b(base, 2.0 * bcrit, bcrit);
        assert!(near < far);
        assert!(far <= base);
    }

    #[test]
    fn spectral_band_parsing_aliases() {
        assert_eq!(SpectralBand::parse("mm"), Some(SpectralBand::Millimeter));
        assert_eq!(SpectralBand::parse("x-ray"), Some(SpectralBand::Xray));
        assert_eq!(SpectralBand::parse("all"), Some(SpectralBand::Bolometric));
        assert_eq!(SpectralBand::parse("unknown"), None);
    }

    #[test]
    fn cpu_cuda_core_palette_coefficients_stay_in_lockstep() {
        // Guard GRAND-215 parity: CUDA tracer palette constants must match CPU path.
        let cuda = include_str!("../../kernels/tracer.cu");

        // gutoe_core_color / bh_gutoe_core_color
        assert!(cuda.contains("255.0 * pow(bv, 0.35)"));
        assert!(cuda.contains("160.0 * pow(bv, 0.65)"));
        assert!(cuda.contains("30.0 * pow(bv, 2.0)"));

        // gutoe_core_physics_color / bh_gutoe_core_physics_color
        assert!(cuda.contains("255.0 * pow(tone, 0.36)"));
        assert!(cuda.contains("170.0 * pow(tone, 0.62)"));
        assert!(cuda.contains("45.0 * pow(tone, 1.50)"));
    }
}

fn main() {
    // CLI:
    //   bh_render                        → all 10 views at default res (1200)
    //   bh_render 3840x2160              → all 10 views at 3840×2160
    //   bh_render 3840x2160 6            → all 10 views at 3840×2160 + blur σ=6
    //   bh_render m87star                → only M87★ at default res
    //   bh_render m87star 3840x2160      → M87★ at 3840×2160
    //   bh_render m87star 3840x2160 6    → M87★ at 4K + blur σ=6
    //   bh_render camera_core_sweep      → controlled core-camera sweep grid
    //   bh_render camera_core_sweep 1024 → sweep at 1024×1024
    //   bh_render interstellar_spin      → rotating camera/disk sequence (Schwarzschild geodesics)
    //   bh_render interstellar_spin 400  → sequence with 400 frames
    //   bh_render spectrum_sweep m87star [WxH] [blur] → render all bands for one view
    //   bh_render campaign_bh [WxH] [blur] → m87star+sgr_astar, GUTOE+GR+difference per band
    //   bh_render validate_bh [slug] [WxH] → convergence + b_crit validation report
    //   bh_render kerr_metrics [slug] [WxH] → sweep Kerr spins and write asymmetry CSV
    //   bh_render transfer_parity [slug] [WxH] → CPU/CUDA covariant-transfer parity CSV
    //   bh_render wave_samples [slug] [WxH] → per-pixel ring/phase/amplitude complex samples
    //   bh_render wave_coherent [slug] [WxH] → coherent intensity + fringe contrast maps
    //   bh_render sgr_astar_eht [WxH] [beam_sigma] → intrinsic + EHT-like 230 GHz pair + metrics
    //   bh_render eht_uv [slug] [WxH] → synthetic visibility + closure CSV/JSON export
    //   bh_render metric_report → dump GUTOE metric table + invariants CSV
    //   bh_render tiled <slug> <WxH> [tile_px] [blur] → deterministic tiled render
    // Sweep env overrides:
    //   BH_SWEEP_R_CAMS="0.60,0.68,0.76"
    //   BH_SWEEP_FOVS="1.4,2.2,3.6"
    //   BH_SWEEP_MAX_PHI_PI=160
    //   BH_SWEEP_DPHI=0.0012
    // Detail env controls:
    //   BH_SUPERSCALE=2        (render at 2x linear resolution, box downsample)
    //   BH_SPP=4               (4 subpixel jittered samples, averaged)
    //   BH_FIXED_EXPOSURE=1.4  (override per-band exposure for non-bolometric renders)
    //   BH_ADAPTIVE_DPHI=1     (smaller per-ray step near b≈b_crit)
    //   BH_KERR_PARITY=1       (with CUDA: run GPU and CPU, print diff stats)
    //   BH_DETAIL_PRESET=imax  (auto: superscale>=2, spp>=4, tighter dphi)
    //   BH_SPECTRUM=radio|mm|infrared|optical|uv|xray|gamma|bolometric
    let args: Vec<String> = std::env::args().collect();
    let run_core_sweep = args.get(1).is_some_and(|s| s == "camera_core_sweep");
    let run_interstellar_spin = args.get(1).is_some_and(|s| s == "interstellar_spin");
    let run_spectrum_sweep = args.get(1).is_some_and(|s| s == "spectrum_sweep");
    let run_campaign_bh = args.get(1).is_some_and(|s| s == "campaign_bh");
    let run_validate_bh = args.get(1).is_some_and(|s| s == "validate_bh");
    let run_kerr_metrics = args.get(1).is_some_and(|s| s == "kerr_metrics");
    let run_transfer_parity = args.get(1).is_some_and(|s| s == "transfer_parity");
    let run_wave_samples = args.get(1).is_some_and(|s| s == "wave_samples");
    let run_wave_coherent = args.get(1).is_some_and(|s| s == "wave_coherent");
    let run_sgr_astar_eht = args.get(1).is_some_and(|s| s == "sgr_astar_eht");
    let run_eht_uv = args.get(1).is_some_and(|s| s == "eht_uv");
    let run_metric_report_cmd = args.get(1).is_some_and(|s| s == "metric_report");
    let run_tiled = args.get(1).is_some_and(|s| s == "tiled");

    // Arg 1 is a slug if it matches a known view slug.
    // Otherwise, if it looks like a resolution, treat it as global res for all views.
    let (filter_slug, width_override, height_override, blur_sigma) = {
        let a1 = args.get(1).map(|s| s.as_str());
        match a1 {
            None => (None, 0usize, 0usize, 0.0f64),
            Some("camera_core_sweep") => {
                // bh_render camera_core_sweep [WxH] [blur]
                let (w, h) = args.get(2).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                let blur = args
                    .get(3)
                    .and_then(|x| x.parse::<f64>().ok())
                    .unwrap_or(0.0);
                (None, w, h, blur)
            }
            Some("interstellar_spin") => {
                // bh_render interstellar_spin [frames] [WxH] [blur]
                // `frames` is consumed separately below.
                let (w, h) = args.get(3).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                let blur = args
                    .get(4)
                    .and_then(|x| x.parse::<f64>().ok())
                    .unwrap_or(0.0);
                (None, w, h, blur)
            }
            Some("spectrum_sweep") => {
                // bh_render spectrum_sweep <slug> [WxH] [blur]
                let (w, h) = args.get(3).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                let blur = args
                    .get(4)
                    .and_then(|x| x.parse::<f64>().ok())
                    .unwrap_or(0.0);
                (None, w, h, blur)
            }
            Some("campaign_bh") => {
                // bh_render campaign_bh [WxH] [blur]
                let (w, h) = args.get(2).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                let blur = args
                    .get(3)
                    .and_then(|x| x.parse::<f64>().ok())
                    .unwrap_or(0.0);
                (None, w, h, blur)
            }
            Some("validate_bh") => {
                // bh_render validate_bh [slug] [WxH]
                let (w, h) = args.get(3).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                (None, w, h, 0.0)
            }
            Some("kerr_metrics") => {
                // bh_render kerr_metrics [slug] [WxH]
                let (w, h) = args.get(3).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                (None, w, h, 0.0)
            }
            Some("transfer_parity") => {
                // bh_render transfer_parity [slug] [WxH]
                let (w, h) = args.get(3).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                (None, w, h, 0.0)
            }
            Some("wave_samples") => {
                // bh_render wave_samples [slug] [WxH]
                let (w, h) = args.get(3).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                (None, w, h, 0.0)
            }
            Some("wave_coherent") => {
                // bh_render wave_coherent [slug] [WxH]
                let (w, h) = args.get(3).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                (None, w, h, 0.0)
            }
            Some("sgr_astar_eht") => {
                // bh_render sgr_astar_eht [WxH] [beam_sigma]
                let (w, h) = args.get(2).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                let blur = args
                    .get(3)
                    .and_then(|x| x.parse::<f64>().ok())
                    .unwrap_or(0.0);
                (None, w, h, blur)
            }
            Some("eht_uv") => {
                // bh_render eht_uv [slug] [WxH]
                let (w, h) = args.get(3).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                (None, w, h, 0.0)
            }
            Some("metric_report") => (None, 0, 0, 0.0),
            Some("tiled") => {
                // bh_render tiled <slug> <WxH> [tile_px] [blur]
                let (w, h) = args.get(3).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                let blur = args
                    .get(5)
                    .and_then(|x| x.parse::<f64>().ok())
                    .unwrap_or(0.0);
                (None, w, h, blur)
            }
            Some(s) if looks_like_resolution(s) => {
                // bh_render <WxH> [blur]
                let (w, h) = parse_resolution(s);
                let blur = args
                    .get(2)
                    .and_then(|x| x.parse::<f64>().ok())
                    .unwrap_or(0.0);
                (None, w, h, blur)
            }
            Some(s) => {
                // bh_render <slug> [WxH] [blur]
                let (w, h) = args.get(2).map(|x| parse_resolution(x)).unwrap_or((0, 0));
                let blur = args
                    .get(3)
                    .and_then(|x| x.parse::<f64>().ok())
                    .unwrap_or(0.0);
                (Some(s), w, h, blur)
            }
        }
    };

    let views_to_render: Vec<&View> = VIEWS
        .iter()
        .filter(|v| filter_slug.map_or(true, |slug| v.slug == slug))
        .collect();

    if views_to_render.is_empty() {
        eprintln!(
            "No view matching slug {:?}.  Known slugs:",
            filter_slug.unwrap_or("")
        );
        for v in VIEWS {
            eprintln!("  {}", v.slug);
        }
        eprintln!("Special mode:");
        eprintln!("  camera_core_sweep [WxH] [blur]");
        eprintln!("  spectrum_sweep <slug> [WxH] [blur]");
        eprintln!("  kerr_metrics [slug] [WxH]");
        eprintln!("  transfer_parity [slug] [WxH]");
        eprintln!("  wave_samples [slug] [WxH]");
        eprintln!("  wave_coherent [slug] [WxH]");
        eprintln!("  sgr_astar_eht [WxH] [beam_sigma]");
        eprintln!("  eht_uv [slug] [WxH]");
        eprintln!("  metric_report");
        eprintln!("  tiled <slug> <WxH> [tile_px] [blur]");
        std::process::exit(1);
    }

    let out_dir = PathBuf::from("/tmp/bh_renders");
    fs::create_dir_all(&out_dir).expect("create /tmp/bh_renders");

    if run_core_sweep {
        eprintln!("GUTOE Black Hole Gallery — controlled camera_core sweep …\n");
        render_camera_core_sweep(&out_dir, width_override, height_override, blur_sigma);
        return;
    }

    if run_interstellar_spin {
        let frames = args
            .get(2)
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(240);
        eprintln!("GUTOE Black Hole Gallery — interstellar spin sequence …\n");
        render_interstellar_spin(&out_dir, width_override, height_override, blur_sigma, frames);
        return;
    }

    if run_spectrum_sweep {
        let sweep_slug = args.get(2).map(String::as_str).unwrap_or("m87star");
        let Some(view) = VIEWS.iter().find(|v| v.slug == sweep_slug) else {
            eprintln!("spectrum_sweep view slug not found: {sweep_slug}");
            std::process::exit(1);
        };
        eprintln!(
            "GUTOE Black Hole Gallery — spectrum sweep for {} …\n",
            view.slug
        );
        render_spectrum_sweep(&out_dir, view, width_override, height_override, blur_sigma);
        return;
    }

    if run_campaign_bh {
        eprintln!("GUTOE Black Hole Gallery — BH campaign (spectrum + GR diff) …\n");
        render_bh_campaign(&out_dir, width_override, height_override, blur_sigma);
        return;
    }
    if run_validate_bh {
        let slug = args.get(2).map(String::as_str).unwrap_or("m87star");
        let Some(view) = VIEWS.iter().find(|v| v.slug == slug) else {
            eprintln!("validate_bh view slug not found: {slug}");
            std::process::exit(1);
        };
        let w = if width_override > 0 { width_override } else { 960 };
        let h = if height_override > 0 { height_override } else { w };
        run_validation_report(view, w, h);
        return;
    }
    if run_kerr_metrics {
        let slug = args.get(2).map(String::as_str).unwrap_or("m87star");
        let Some(view) = VIEWS.iter().find(|v| v.slug == slug) else {
            eprintln!("kerr_metrics view slug not found: {slug}");
            std::process::exit(1);
        };
        let w = if width_override > 0 { width_override } else { 960 };
        let h = if height_override > 0 { height_override } else { w };
        run_kerr_metrics_report(view, w, h, &out_dir);
        return;
    }
    if run_transfer_parity {
        let slug = args.get(2).map(String::as_str).unwrap_or("m87star");
        let Some(view) = VIEWS.iter().find(|v| v.slug == slug) else {
            eprintln!("transfer_parity view slug not found: {slug}");
            std::process::exit(1);
        };
        let w = if width_override > 0 { width_override } else { 960 };
        let h = if height_override > 0 { height_override } else { w };
        run_transfer_parity_report(view, w, h, &out_dir);
        return;
    }
    if run_wave_samples {
        let slug = args.get(2).map(String::as_str).unwrap_or("m87star");
        let Some(view) = VIEWS.iter().find(|v| v.slug == slug) else {
            eprintln!("wave_samples view slug not found: {slug}");
            std::process::exit(1);
        };
        let w = if width_override > 0 { width_override } else { 512 };
        let h = if height_override > 0 { height_override } else { w };
        run_wave_samples_report(view, w, h, &out_dir);
        return;
    }
    if run_wave_coherent {
        let slug = args.get(2).map(String::as_str).unwrap_or("m87star");
        let Some(view) = VIEWS.iter().find(|v| v.slug == slug) else {
            eprintln!("wave_coherent view slug not found: {slug}");
            std::process::exit(1);
        };
        let w = if width_override > 0 { width_override } else { 512 };
        let h = if height_override > 0 { height_override } else { w };
        run_wave_coherent_report(view, w, h, &out_dir);
        return;
    }
    if run_sgr_astar_eht {
        let w = if width_override > 0 { width_override } else { 1920 };
        let h = if height_override > 0 { height_override } else { 1080 };
        run_sgr_astar_eht_report(&out_dir, w, h, blur_sigma);
        return;
    }
    if run_eht_uv {
        let slug = args.get(2).map(String::as_str).unwrap_or("m87_eht2017");
        let Some(view) = VIEWS.iter().find(|v| v.slug == slug) else {
            eprintln!("eht_uv view slug not found: {slug}");
            std::process::exit(1);
        };
        let w = if width_override > 0 { width_override } else { 1024 };
        let h = if height_override > 0 { height_override } else { w };
        run_eht_uv_export(view, w, h, &out_dir);
        return;
    }
    if run_metric_report_cmd {
        run_metric_report(&out_dir);
        return;
    }
    if run_tiled {
        let slug = args.get(2).map(String::as_str).unwrap_or("m87star");
        let Some(view) = VIEWS.iter().find(|v| v.slug == slug) else {
            eprintln!("tiled view slug not found: {slug}");
            std::process::exit(1);
        };
        let (rw, rh) = args.get(3).map(|x| parse_resolution(x)).unwrap_or((0, 0));
        let w = if rw > 0 { rw } else { 7680 };
        let h = if rh > 0 { rh } else { w };
        let tile_px = args
            .get(4)
            .and_then(|x| x.parse::<usize>().ok())
            .unwrap_or(1024);
        render_view_tiled(&out_dir, view, w, h, tile_px, blur_sigma);
        return;
    }

    eprintln!(
        "GUTOE Black Hole Gallery — rendering {} view(s) …\n",
        views_to_render.len()
    );
    for view in &views_to_render {
        render_view(&out_dir, view, width_override, height_override, blur_sigma);
    }
    eprintln!("\nAll {} render(s) complete.", views_to_render.len());

    let html = Arc::new(build_html());
    let dir = Arc::new(out_dir);
    serve_http(dir, html);
}
