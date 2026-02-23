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
    fs,
    io::{Read, Write},
    net::{TcpListener, TcpStream},
    path::{Path, PathBuf},
    process::Command,
    sync::Arc,
};

use rayon::prelude::*;

use gutoe_gpu::{
    kerr::KerrMetric,
    metric::GutoeMetric,
    tracer::{
        b_critical, trace_photon, trace_photon_interior, trace_photon_interior_core, write_ppm,
        RenderConfig, TraceResult,
    },
};

// ── CUDA FFI (only available with --features cuda) ───────────────────────────

#[cfg(feature = "cuda")]
extern "C" {
    /// GPU black hole renderer — one CUDA thread per pixel.
    /// Writes RGB triples to `out_pixels` (host buffer, width*height*3 bytes).
    fn gutoe_render_bh(
        width: i32,
        height: i32,
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
        r_cam_rs: f64, // 0.0 = exterior; >0 = interior camera at r_cam_rs × r_s
        out_pixels: *mut u8,
    );
}

/// GPU render path — calls the CUDA kernel and returns the pixel buffer.
#[cfg(feature = "cuda")]
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
    r_cam_rs: f64,
) -> Vec<[u8; 3]> {
    let n = cfg.width * cfg.height;
    let mut flat = vec![0u8; n * 3];
    unsafe {
        gutoe_render_bh(
            cfg.width as i32,
            cfg.height as i32,
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
            r_cam_rs,
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
    sin_inc: f64,
    n_cross: u32,
    doppler: bool,
    ring_mode: bool,
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

    // Relativistic Keplerian Doppler: v_K = c · √(r_s / (2 r_eff))
    // Observer-velocity alignment: β_obs ≈ v_K · sin_inc · (bx / r_eff)
    // bx > 0 → disk material moving toward observer → blueshift → D > 1.
    let doppler_d4 = if doppler {
        let beta = (r_s / (2.0 * r_eff)).sqrt().min(0.5);
        let beta_obs = beta * sin_inc * (bx / (r_eff.max(1e-12))).clamp(-1.0, 1.0);
        // D = 1 / (1 − β_obs);  D⁴ for thermal emission (Stefan–Boltzmann)
        (1.0 / (1.0 - beta_obs)).powi(4).clamp(0.01, 200.0)
    } else {
        1.0
    };

    // Reinhard tone mapping: luminance → luminance / (1 + luminance)
    // Never clips; handles both dim and Doppler-boosted pixels uniformly.
    let luminance = (t_rel * fade * doppler_d4 * outer_taper).max(0.0);
    let b = luminance / (1.0 + luminance);

    // Orange-white thermal palette (hot inner disk white, outer disk orange-red)
    let r = (255.0 * b.powf(0.35)).clamp(0.0, 255.0) as u8;
    let g = (210.0 * b.powf(0.60)).clamp(0.0, 255.0) as u8;
    let bl = (130.0 * b.powf(1.60)).clamp(0.0, 255.0) as u8;
    [r, g, bl]
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
    // Rotate impact-parameter vector (bx, by) by phi_total to get sky direction.
    // In the Schwarzschild geometry a photon starting from angle 0 exits at
    // phi_total → source is at angle phi_total from the observer direction.
    // phi_total ≈ π : direct (weakly deflected) ray — sky ≈ "behind" BH
    // phi_total ≈ 2π: once-orbiting ray   — Einstein ring; same stars as direct
    // phi_total ≫ π : photon ring regime  — rapidly changing star backgrounds
    let sky_x = bx * phi_total.cos() - by * phi_total.sin();
    let sky_y = bx * phi_total.sin() + by * phi_total.cos();

    // Map to integer sky-grid (one cell ≈ 0.02 r_s, fine enough for 4K)
    const SKY_SCALE: f64 = 50.0;
    let hx = (sky_x * SKY_SCALE).round() as i64;
    let hy = (sky_y * SKY_SCALE).round() as i64;
    let h = star_hash(hx, hy);

    // low 16 bits → density gate: 1000/65536 ≈ 1.5 %
    if (h & 0xFFFF) as u32 >= 1000 {
        // Background: very dark blue-black with a touch of variation
        let v = ((h >> 40) & 7) as u8; // 0..7
        return [v >> 2, v >> 2, (v >> 1) + 8]; // [0-1, 0-1, 8-11]
    }

    // Star brightness: bits 16-23 → 80..255
    let bright_raw = ((h >> 16) & 0xFF) as u32;
    let bright = (80 + bright_raw * 175 / 255) as u8;

    // Spectral type from bits 24-27 (4 bits → 16 classes)
    match (h >> 24) & 0xF {
        0..=1 => {
            // Deep red — cool M-dwarf
            [bright / 2, bright / 5, bright / 10]
        }
        2..=4 => {
            // Orange — K-star
            let g = (bright as u32 * 68 / 100) as u8;
            let b = (bright as u32 * 38 / 100) as u8;
            [bright, g, b]
        }
        5..=9 => {
            // Warm white — F/G (sun-like)
            let b = (bright as u32 * 88 / 100) as u8;
            [bright, bright, b]
        }
        10..=12 => {
            // Pure white — A-star
            [bright, bright, bright]
        }
        _ => {
            // Blue-white — hot B/O star
            let rg = (bright as u32 * 82 / 100) as u8;
            [rg, rg, bright]
        }
    }
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
    let effective_r_cam_frac = r_cam_override.unwrap_or(view.r_cam_frac);
    let effective_fov = fov_override.unwrap_or(view.fov);
    let effective_max_phi_pi = max_phi_pi_override.unwrap_or(view.max_phi_pi);
    let effective_az = az_override.unwrap_or(view.az);
    let effective_inc = inc_override.unwrap_or(view.inc);
    // Use finer integration step at high resolution so the photon ring stays crisp
    let base_dphi = if width >= 2048 {
        view.dphi.min(0.003)
    } else {
        view.dphi
    };
    let dphi = dphi_override.unwrap_or(base_dphi);
    let blur_str = if blur_sigma >= 0.5 {
        format!("blur={blur_sigma:.1}")
    } else {
        "no blur".into()
    };
    eprintln!(
        "  rendering {}  (inc={:.0}°, az={:.0}°, doppler={}, rings={}, gr={}, fov={:.2}, r_cam_frac={:.3}, max_phi_pi={:.1}, {}×{}, dphi={:.4}, {}) …",
        view.label, effective_inc, effective_az, view.doppler, view.ring_mode, view.gr_mode,
        effective_fov, effective_r_cam_frac, effective_max_phi_pi, width, height, dphi, blur_str,
    );

    let metric = if view.gr_mode {
        GutoeMetric::schwarzschild(1.0)
    } else {
        GutoeMetric::planck_units(1.0)
    };

    // Kerr is not wired into the geodesic integrator yet. If requested, report the
    // exact Kerr invariants and require explicit placeholder opt-in.
    if let Some(a_star) = parse_env_f64("BH_KERR_ASTAR") {
        if let Some(kerr) = KerrMetric::new(metric.r_s, a_star) {
            let (r_plus, r_minus) = kerr.horizons();
            let r_erg_eq = kerr.ergosphere_radius(PI / 2.0);
            let r_ph_pro = kerr.equatorial_photon_orbit_radius(true);
            let r_ph_ret = kerr.equatorial_photon_orbit_radius(false);
            eprintln!(
                "    [KERR BASELINE] a*={:.3} r+={:.4} r-={:.4} r_erg,eq={:.4} r_ph,pro={:.4} r_ph,ret={:.4} Ω_H={:.5}",
                a_star, r_plus, r_minus, r_erg_eq, r_ph_pro, r_ph_ret, kerr.horizon_angular_velocity()
            );
            let allow_placeholder = std::env::var("BH_ALLOW_KERR_PLACEHOLDER")
                .ok()
                .is_some_and(|v| v == "1" || v.eq_ignore_ascii_case("true"));
            assert!(
                allow_placeholder,
                "BH_KERR_ASTAR is set but Kerr geodesic rendering is not implemented yet. \
Set BH_ALLOW_KERR_PLACEHOLDER=1 only if you explicitly want Schwarzschild fallback while prototyping."
            );
        } else {
            panic!("invalid BH_KERR_ASTAR={a_star}: expected |a*| <= 1");
        }
    }

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

    let cfg = RenderConfig {
        width,
        height,
        fov_rs: effective_fov,
        inclination_deg: effective_inc,
        max_phi: effective_max_phi_pi * PI,
        dphi,
    };

    let raw = render_with_options(
        &metric,
        view.disk_inner,
        view.disk_outer,
        &cfg,
        effective_az,
        view.doppler,
        view.ring_mode,
        view.interior_mode,
        view.core_look_mode,
        r_cam_rs,
    );

    // Apply Gaussian PSF blur (EHT beam simulation) if requested
    let pixels = if blur_sigma >= 0.5 {
        eprintln!("    applying Gaussian blur σ={blur_sigma:.1} …");
        gaussian_blur(&raw, cfg.width, cfg.height, blur_sigma)
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

    fs::write(&ppm_path, write_ppm(&pixels, cfg.width, cfg.height)).expect("write ppm");

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

    eprintln!(
        "Interstellar spin complete: {n} frames in {}",
        out_dir.display()
    );
}

// ── Full render with options ──────────────────────────────────────────────────
//
// Pixel loop: azimuth rotation + Doppler + ring-mode colouring.
// When az=0 and doppler=false and ring_mode=false, delegates to the fast
// tracer::render() which skips the per-pixel boxing overhead.

fn render_with_options(
    metric: &GutoeMetric,
    disk_inner_rs: f64,
    disk_outer_rs: f64,
    cfg: &RenderConfig,
    az_deg: f64,
    doppler: bool,
    ring_mode: bool,
    interior_mode: bool,
    core_look_mode: bool,
    r_cam_rs: f64, // 0.0 = exterior; >0 = camera inside at r_cam = r_cam_rs × r_s
) -> Vec<[u8; 3]> {
    // GPU path: delegate to CUDA kernel (all pixels in one launch, ~1 second)
    #[cfg(feature = "cuda")]
    {
        eprintln!(
            "    [GPU] launching CUDA kernel ({} × {} = {} pixels) …",
            cfg.width,
            cfg.height,
            cfg.width * cfg.height
        );
        let t0 = std::time::Instant::now();
        let result = render_with_options_gpu(
            metric,
            disk_inner_rs,
            disk_outer_rs,
            cfg,
            az_deg,
            doppler,
            ring_mode,
            interior_mode,
            core_look_mode,
            r_cam_rs,
        );
        eprintln!("    [GPU] done in {:.2}s", t0.elapsed().as_secs_f64());
        return result;
    }

    // CPU path: rayon parallel pixel loop (16 threads on local machine)
    #[allow(unreachable_code)]
    {
        let r_s = metric.r_s;
        let disk_inner = disk_inner_rs * r_s;
        let disk_outer = disk_outer_rs * r_s;
        let r_isco = 3.0 * r_s;
        let b_crit = b_critical(r_s);
        let sin_inc = cfg.inclination_deg.to_radians().sin();
        let scale = 2.0 * cfg.fov_rs * r_s / cfg.width as f64;
        let az_rad = az_deg.to_radians();
        let (ca, sa) = (az_rad.cos(), az_rad.sin());
        let width = cfg.width;
        let height = cfg.height;
        let r_cam = r_cam_rs * r_s; // 0.0 = exterior

        (0..height * width)
            .into_par_iter()
            .map(|idx| {
                let iy = idx / width;
                let ix = idx % width;
                let sx = (ix as f64 - 0.5 * (width as f64 - 1.0)) * scale;
                let sy = (0.5 * (height as f64 - 1.0) - iy as f64) * scale;
                let bx_raw = sx;
                let by_raw = sy * sin_inc;
                let bx = ca * bx_raw - sa * by_raw;
                let by = sa * bx_raw + ca * by_raw;

                let result = if r_cam > 0.0 {
                    if core_look_mode {
                        trace_photon_interior_core(metric, r_cam, bx, by, cfg.max_phi, cfg.dphi)
                    } else {
                        trace_photon_interior(
                            metric,
                            disk_inner,
                            disk_outer,
                            r_cam,
                            bx,
                            by,
                            cfg.max_phi,
                            cfg.dphi,
                        )
                    }
                } else {
                    trace_photon(
                        metric,
                        disk_inner,
                        disk_outer,
                        bx,
                        by,
                        cfg.max_phi,
                        cfg.dphi,
                    )
                };

                let b_mag = (bx * bx + by * by).sqrt();
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
                            pixel_color(
                                r_eff, r_isco, disk_outer, r_s, bx_raw, sin_inc, n_cross, doppler,
                                ring_mode,
                            )
                        }
                    }
                }
            })
            .collect()
    }
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
    let listener = TcpListener::bind(addr).expect("bind 0.0.0.0:52345");
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
    // Sweep env overrides:
    //   BH_SWEEP_R_CAMS="0.60,0.68,0.76"
    //   BH_SWEEP_FOVS="1.4,2.2,3.6"
    //   BH_SWEEP_MAX_PHI_PI=160
    //   BH_SWEEP_DPHI=0.0012
    let args: Vec<String> = std::env::args().collect();
    let run_core_sweep = args.get(1).is_some_and(|s| s == "camera_core_sweep");
    let run_interstellar_spin = args.get(1).is_some_and(|s| s == "interstellar_spin");

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
