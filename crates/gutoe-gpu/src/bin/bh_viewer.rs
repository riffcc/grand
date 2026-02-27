//! GUTOE Black Hole Viewer — real-time interactive geodesic ray tracer
//!
//! Runs the GUTOE-corrected Schwarzschild metric null-geodesic integrator
//! entirely on the GPU in a WGSL fragment shader.
//!
//! Controls:
//!   Left drag (vertical)   — inclination: 0° face-on ↔ 90° edge-on
//!   Left drag (horizontal) — disk plane rotation (azimuth)
//!   Scroll                 — zoom (field of view in units of r_s)
//!   +/-                    — disk outer radius
//!   G                      — toggle GUTOE lattice core (r_c on/off)
//!   R                      — reset camera
//!   Q / Escape             — quit

use std::{
    panic::{self, AssertUnwindSafe},
    sync::Arc,
    time::{Duration, Instant},
};

use gilrs::{Axis, Button, Gilrs};
use winit::{
    application::ApplicationHandler,
    dpi::PhysicalPosition,
    event::{ElementState, KeyEvent, MouseButton, MouseScrollDelta, WindowEvent},
    event_loop::{ActiveEventLoop, EventLoop},
    keyboard::{Key, NamedKey},
    window::{Window, WindowId},
};

#[derive(Default, Clone, Copy)]
struct PadPrev {
    south: bool,
    east: bool,
    west: bool,
    north: bool,
    left_thumb: bool,
    start: bool,
    select: bool,
    mode: bool,
}

#[derive(Default, Clone, Copy)]
struct HeldKeys {
    fwd: bool,
    back: bool,
    left: bool,
    right: bool,
    up: bool,
    down: bool,
    roll_pos: bool,
    roll_neg: bool,
    zoom_in: bool,
    zoom_out: bool,
}

#[derive(Clone, Copy)]
struct PadTuning {
    deadzone: f32,
    look_sens: f32,
    move_sens: f32,
    trigger_sens: f32,
    invert_look_x: bool,
    invert_look_y: bool,
}

impl Default for PadTuning {
    fn default() -> Self {
        Self {
            deadzone: 0.12,
            look_sens: 1.0,
            move_sens: 1.0,
            trigger_sens: 1.0,
            invert_look_x: false,
            invert_look_y: false,
        }
    }
}

// ── WGSL shader ───────────────────────────────────────────────────────────────

const SHADER: &str = r#"
// GUTOE Black Hole Shader — null geodesic tracer with stars, Doppler, gravitational redshift
//
// Physics:
//   d²r/dφ² = r(2r²+r_c²)/b² − r + r_s·r(r²+2r_c²)/(2r_eff³)   [GUTOE Schwarzschild]
//   Disk: Novikov-Thorne T ∝ (r_ISCO/r)^(3/4)
//   Doppler: D³ factor for approaching / receding disk gas (Keplerian orbit)
//   Gravitational redshift: sqrt(1 - r_s/r_eff)
//   Stars: procedural hash field; gravitationally lensed by tracing exit angle

struct Params {
    r_s      : f32,   // Schwarzschild radius (internal units = 1)
    r_c      : f32,   // GUTOE lattice core radius = sqrt(C_inf) × l_P
    disk_in  : f32,   // disk inner areal radius (× r_s)
    disk_out : f32,   // disk outer areal radius (× r_s)
    sin_inc  : f32,   // sin(observer inclination): 1=edge-on, 0=face-on
    fov      : f32,   // half-width of image in r_s
    width    : f32,   // viewport width  (pixels)
    height   : f32,   // viewport height (pixels)
    max_phi  : f32,   // max integration angle (radians)
    dphi     : f32,   // RK4 step size (radians)
    az       : f32,   // azimuth offset (rotates disk orientation on screen)
    cam_x    : f32,   // freecam world-space x offset (r_s units)
    cam_y    : f32,   // freecam world-space y offset (r_s units)
    cam_z    : f32,   // freecam world-space z offset (r_s units)
    cam_roll : f32,   // camera roll (radians)
    interior_mode : f32, // 1 = interior camera
    core_look_mode: f32, // 1 = interior camera looking at core
    r_cam_frac    : f32, // r_cam = r_cam_frac * r_horizon (coordinate)
    quality_tier  : f32, // 0=low, 1=medium, 2=high
    use_transfer  : f32, // 1 => run multi-step covariant transfer integration
    tau_scale     : f32, // optical-depth scale for transfer mode
    disk_model    : f32, // 0=thin, 1=riaf composite
    local_stars   : f32, // 1 => nearby 3D star volume shells
    riaf_volume   : f32, // 1 => escaped-ray volumetric RIAF blend
    kerr_enable   : f32, // 1 => use Kerr tracer for exterior camera
    kerr_astar    : f32, // dimensionless spin a*
    _pad1         : f32,
    _pad2         : f32,
    _pad3         : f32,
}
@group(0) @binding(0) var<uniform> P : Params;
@group(0) @binding(1) var star_tex : texture_2d<f32>;
@group(0) @binding(2) var star_samp : sampler;

// ── Hash / Noise ──────────────────────────────────────────────────────────────

fn hash21(p: vec2<f32>) -> f32 {
    var q = fract(p * vec2(127.1, 311.7));
    q += dot(q, q.yx + 19.19);
    return fract(q.x * q.y);
}

fn hash31(p: vec3<f32>) -> f32 {
    var q = fract(p * vec3(127.1, 311.7, 269.5));
    q += dot(q, q.yzx + 19.19);
    return fract((q.x + q.y) * q.z);
}

fn local_star_volume(cam: vec3<f32>, dir: vec3<f32>, depth: f32, scale: f32, quality: f32) -> vec3<f32> {
    let steps = select(select(6u, 10u, quality >= 1.0), 14u, quality >= 1.5);
    let star_prob = mix(0.0018, 0.0030, clamp(quality / 2.0, 0.0, 1.0));
    var c = vec3(0.0);
    let inv = 1.0 / max(f32(steps), 1.0);
    for (var i = 0u; i < 20u; i++) {
        if (i >= steps) {
            break;
        }
        let t = (f32(i) + 0.5) * depth * inv;
        let p = (cam + dir * t) * scale;
        let cell = floor(p);
        let local = fract(p) - vec3(0.5);
        let gate = hash31(cell + vec3(17.0 + f32(i) * 3.7, 31.0, 47.0));
        if gate > star_prob {
            continue;
        }
        let center = vec3(
            hash31(cell + vec3(1.3, 2.1, 3.7)),
            hash31(cell + vec3(4.1, 5.9, 6.8)),
            hash31(cell + vec3(7.2, 8.6, 9.4))
        ) - vec3(0.5);
        let d = local - center;
        let r2 = dot(d, d);
        let psf = exp(-r2 * 52.0);
        let temp = hash31(cell + vec3(0.9, 1.1, 1.7));
        let star_col = select(
            select(
                select(vec3(1.00, 0.62, 0.35), vec3(1.00, 0.92, 0.78), temp > 0.20),
                vec3(1.00, 1.00, 1.00),
                temp > 0.45
            ),
            vec3(0.78, 0.86, 1.00),
            temp > 0.75
        );
        let bright = (0.06 + 0.34 * pow(hash31(cell + vec3(3.9, 2.7, 1.5)), 2.0)) * (1.0 + 0.35 * quality);
        c += star_col * bright * psf * inv * 2.4;
    }
    return c;
}

// ── Star field ────────────────────────────────────────────────────────────────
// Procedural star field matched to bh_render's hash/band model for parity.

fn starfield_from_dir(dir_in: vec3<f32>) -> vec3<f32> {
    let dir = safe_normalize(dir_in, vec3(0.0, 0.0, 1.0));
    let lon = atan2(dir.x, dir.z);
    let lat = asin(clamp(dir.y, -1.0, 1.0));
    let sky = vec2(lon * 57.2957795, lat * 114.591559);
    var col = vec3(0.0015, 0.0020, 0.0060); // dark sky baseline

    // Milky-Way-like galactic band with smooth dust modulation.
    let gx = sky.x * 0.035;
    let gy = sky.y * 0.035;
    let gal_lat = abs(gy + 0.28 * sin(0.7 * gx) + 0.14 * sin(1.3 * gx + 0.5));
    let band = exp(-(gal_lat * gal_lat) / 0.05);
    let dust = pow(sin(sky.x * 0.11) * cos(sky.y * 0.07 + 1.7) * 0.5 + 0.5, 1.4);
    let gal = band * (0.25 + 0.75 * dust);
    col += vec3(0.16, 0.14, 0.20) * gal;

    // Multi-scale stellar population with sub-pixel PSF, avoiding "cloudy" blocks.
    {
        let scale = 42.0;
        let p = sky * scale;
        let cell = floor(p);
        let local = fract(p) - 0.5;
        let gate = hash21(cell);
        if (gate <= 0.010) {
            let cseed = cell + vec2(hash21(cell + 0.41), hash21(cell + 0.73));
            let center = vec2(hash21(cseed + 1.2), hash21(cseed + 2.6)) - 0.5;
            let d = local - center;
            let r2 = dot(d, d);
            let psf = exp(-r2 * 320.0);
            let temp = hash21(cell + 0.99);
            let star_col = select(
                select(
                    select(vec3(1.00, 0.62, 0.35), vec3(1.00, 0.92, 0.78), temp > 0.20),
                    vec3(1.00, 1.00, 1.00),
                    temp > 0.45
                ),
                vec3(0.78, 0.86, 1.00),
                temp > 0.75
            );
            let bright = 0.95 * (0.25 + 0.75 * pow(hash21(cell + 0.37), 1.8));
            col += star_col * bright * psf * (1.0 + 0.6 * band);
        }
    }
    {
        let scale = 86.0;
        let p = sky * scale;
        let cell = floor(p);
        let local = fract(p) - 0.5;
        let gate = hash21(cell);
        if (P.quality_tier >= 1.0 && gate <= 0.025) {
            let cseed = cell + vec2(hash21(cell + 1.41), hash21(cell + 1.73));
            let center = vec2(hash21(cseed + 3.2), hash21(cseed + 4.6)) - 0.5;
            let d = local - center;
            let r2 = dot(d, d);
            let psf = exp(-r2 * 420.0);
            let temp = hash21(cell + 0.99);
            let star_col = select(
                select(
                    select(vec3(1.00, 0.62, 0.35), vec3(1.00, 0.92, 0.78), temp > 0.20),
                    vec3(1.00, 1.00, 1.00),
                    temp > 0.45
                ),
                vec3(0.78, 0.86, 1.00),
                temp > 0.75
            );
            let bright = 0.30 * (0.25 + 0.75 * pow(hash21(cell + 0.37), 1.8));
            col += star_col * bright * psf * (1.0 + 0.6 * band);
        }
    }
    {
        let scale = 150.0;
        let p = sky * scale;
        let cell = floor(p);
        let local = fract(p) - 0.5;
        let gate = hash21(cell);
        if (P.quality_tier >= 1.5 && gate <= 0.050) {
            let cseed = cell + vec2(hash21(cell + 2.41), hash21(cell + 2.73));
            let center = vec2(hash21(cseed + 5.2), hash21(cseed + 6.6)) - 0.5;
            let d = local - center;
            let r2 = dot(d, d);
            let psf = exp(-r2 * 520.0);
            let temp = hash21(cell + 0.99);
            let star_col = select(
                select(
                    select(vec3(1.00, 0.62, 0.35), vec3(1.00, 0.92, 0.78), temp > 0.20),
                    vec3(1.00, 1.00, 1.00),
                    temp > 0.45
                ),
                vec3(0.78, 0.86, 1.00),
                temp > 0.75
            );
            let bright = 0.12 * (0.25 + 0.75 * pow(hash21(cell + 0.37), 1.8));
            col += star_col * bright * psf * (1.0 + 0.6 * band);
        }
    }

    // Optional real-sky map overlay (same equirect projection used by bh_render).
    let uv = vec2(lon / (2.0 * 3.14159265359) + 0.5, 0.5 - lat / 3.14159265359);
    let map_sample = textureSampleLevel(star_tex, star_samp, uv, 0.0);
    let map_col = map_sample.rgb;
    let map_mix = 0.88 * clamp(map_sample.a, 0.0, 1.0);
    col = mix(col, map_col, map_mix);

    // Local 3D stellar volume shells: introduces coherent camera parallax for
    // nearby stars while preserving the far-field map/procedural sky.
    if (P.local_stars > 0.5) {
        let cam = vec3(P.cam_x, P.cam_y, P.cam_z);
        let rs = max(P.r_s, 1.0);
        col += local_star_volume(cam, dir, 24.0 * rs, 0.20, P.quality_tier) * 0.16;
        if (P.quality_tier >= 1.0) {
            col += local_star_volume(cam, dir, 52.0 * rs, 0.12, P.quality_tier) * 0.10;
        }
        if (P.quality_tier >= 1.5) {
            col += local_star_volume(cam, dir, 96.0 * rs, 0.08, P.quality_tier) * 0.06;
        }
    }
    return col;
}

// ── Orbit integrator ──────────────────────────────────────────────────────────

fn accel(r: f32, b: f32, r_s: f32, r_c: f32) -> f32 {
    let re2 = r * r + r_c * r_c;
    let re3 = re2 * sqrt(re2);
    return r * (2.0 * r * r + r_c * r_c) / (b * b)
         - r
         + r_s * r * (r * r + 2.0 * r_c * r_c) / (2.0 * re3);
}

fn orbit_vr_sq(r: f32, b: f32, r_s: f32, r_c: f32) -> f32 {
    let re2 = r * r + r_c * r_c;
    let re  = sqrt(re2);
    let f   = 1.0 - r_s / re;
    return r * r * re2 / (b * b) - r * r * f;
}

fn rk4(r: f32, p: f32, b: f32, r_s: f32, r_c: f32, h: f32) -> vec2<f32> {
    let k1r = p;                               let k1p = accel(r,                     b, r_s, r_c);
    let k2r = p + 0.5 * h * k1p;              let k2p = accel(r + 0.5 * h * k1r,    b, r_s, r_c);
    let k3r = p + 0.5 * h * k2p;              let k3p = accel(r + 0.5 * h * k2r,    b, r_s, r_c);
    let k4r = p + h * k3p;                    let k4p = accel(r + h * k3r,           b, r_s, r_c);
    return vec2(
        r + h * (k1r + 2.0*k2r + 2.0*k3r + k4r) / 6.0,
        p + h * (k1p + 2.0*k2p + 2.0*k3p + k4p) / 6.0,
    );
}

fn safe_normalize(v: vec3<f32>, fallback: vec3<f32>) -> vec3<f32> {
    let n = length(v);
    if n > 1e-6 {
        return v / n;
    }
    return fallback;
}

struct TraceHit3D {
    hit: vec4<f32>,
    sky_dir: vec3<f32>,
}

// ── Tracer ────────────────────────────────────────────────────────────────────
// Returns vec4(r_eff_hit, phi_total, f32(n_cross), kind)
//   kind 0 = captured, 1 = disk hit, 2 = escaped
// phi_total is the total orbital angle (used for gravitational lensing of stars).

fn trace(bx_in: f32, by_in: f32) -> vec4<f32> {
    // Rotate impact parameters by azimuth (spins the disk orientation on screen)
    let ca = cos(P.az); let sa = sin(P.az);
    let bx = ca * bx_in - sa * by_in;
    let by = sa * bx_in + ca * by_in;

    let b = sqrt(bx * bx + by * by);
    if b < 0.0001 { return vec4(0.0, 0.0, 0.0, 0.0); }

    let r_s   = P.r_s;
    let r_c   = P.r_c;
    // Fast-path: far-field rays are weakly deflected and do not cross the disk.
    // This removes a huge amount of unnecessary integration in live mode.
    let far_cut = select(10.0 * r_s, 14.0 * r_s, P.quality_tier >= 1.5);
    if b > far_cut {
        return vec4(0.0, 3.14159265359, 0.0, 2.0);
    }

    // Deep-shadow shortcut: b < b_crit/2 → definitely captured.
    // b_crit = (3√3/2) r_s ≈ 2.598 r_s; half is (3√3/4) r_s ≈ 1.299 r_s.
    let b_crit_half = 0.75 * sqrt(3.0) * r_s;
    if b < b_crit_half { return vec4(0.0, 0.0, 0.0, 0.0); }

    let sin_i = by / b;
    let is_eq = abs(sin_i) < 1e-5;

    // r_start = 3b: orbit_accel ≈ 2r³/b², p ≈ r²/b → change per step ≈ 6·dphi = 6%.
    // Scale-independent and stable for all b > 0.5·b_crit (deep shadow already handled).
    let r_start = 3.0 * b;

    // Initial radial velocity: ingoing from r_start
    let re0   = sqrt(r_start * r_start + r_c * r_c);
    let v0sq  = r_start * r_start * re0 * re0 / (b * b)
              - r_start * r_start * (1.0 - r_s / re0);
    var p_init = -r_start * r_start / b;
    if v0sq > 0.0 { p_init = -sqrt(v0sq); }

    // Capture threshold: areal radius r_eff = sqrt(r²+r_c²) drops below the horizon.
    // The event horizon is at areal radius r_s (not the coordinate radius sqrt(r_s²-r_c²)).
    let r_cap = r_s * 0.99;

    var r          = r_start;
    var p          = p_init;
    var phi        = 0.0f;
    var n_cross    = 0u;
    var turned     = false;
    var in_disk_eq = false;
    let b_crit     = 1.5 * sqrt(3.0) * r_s;
    let rel        = abs(b / max(b_crit, 1e-6) - 1.0);
    // Live-mode integrator budget: avoid runaway step counts near b_crit.
    let step_scale =
        select(1.0, 0.85, rel < 0.10) *
        select(1.0, 0.75, rel < 0.05) *
        select(1.0, 0.65, rel < 0.02);
    let h = max(P.dphi * step_scale, P.dphi * 0.55);
    let max_steps_nominal = i32(P.max_phi / h) + 1;
    let max_steps_cap = select(select(2400, 3600, P.quality_tier >= 1.0), 5200, P.quality_tier >= 1.5);
    let max_steps = min(max_steps_nominal, max_steps_cap);

    for (var i = 0; i < max_steps; i++) {
        let s      = rk4(r, p, b, r_s, r_c, h);
        let rn     = s.x;
        let pn_rk4 = s.y;
        // Enforce orbital constraint p² = orbit_vr_sq(r). Prevents centrifugal blowup;
        // preserves direction sign from RK4 → correctly triggers turned / capture.
        let vr2n   = max(orbit_vr_sq(rn, b, r_s, r_c), 0.0);
        let pn     = select(-sqrt(vr2n), sqrt(vr2n), pn_rk4 >= 0.0);
        let phin = phi + h;
        let ren  = sqrt(rn * rn + r_c * r_c);

        // Capture: inside horizon, or coordinate r went below core / negative
        if rn < r_c * 0.01 || ren < r_cap { return vec4(0.0, 0.0, 0.0, 0.0); }

        // Turning point (ingoing → outgoing)
        if !turned && p < 0.0 && pn >= 0.0 { turned = true; }

        // Escape: returned to launch radius
        if turned && rn >= r_start * 0.99 { return vec4(0.0, phin, 0.0, 2.0); }

        // Disk detection
        if is_eq {
            let re_cur = sqrt(r * r + r_c * r_c);
            let in_d   = re_cur >= P.disk_in && re_cur <= P.disk_out;
            if !in_disk_eq && in_d && p < 0.0 {
                return vec4(re_cur, phi, 1.0, 1.0);
            }
            in_disk_eq = in_d;
        } else {
            let phi_cross = f32(n_cross + 1u) * 3.14159265359;
            if phi < phi_cross && phin >= phi_cross {
                let t    = (phi_cross - phi) / h;
                let r_x  = r + t * (rn - r);
                let re_x = sqrt(r_x * r_x + r_c * r_c);
                n_cross += 1u;
                if re_x >= P.disk_in && re_x <= P.disk_out {
                    return vec4(re_x, phi_cross, f32(n_cross), 1.0);
                }
            }
        }

        r   = rn;
        p   = pn;
        phi = phin;
    }

    if r >= r_start * 0.5 { return vec4(0.0, phi, 0.0, 2.0); }
    return vec4(0.0, 0.0, 0.0, 0.0);
}

// Full 3D exterior tracer in spherical symmetry:
// integrates geodesic in orbital plane, then lifts each step to world space and
// checks equatorial plane crossings (z=0) for disk hits.
fn trace_true3d(sx: f32, sy: f32) -> TraceHit3D {
    let r_s = P.r_s;
    let r_c = P.r_c;
    let z_obs = 60.0 * r_s;
    let sin_i = clamp(P.sin_inc, -1.0, 1.0);
    let cos_i = sqrt(max(1.0 - sin_i * sin_i, 0.0));
    let caz = cos(P.az);
    let saz = sin(P.az);
    let obs_base = vec3(z_obs * sin_i * caz, z_obs * sin_i * saz, z_obs * cos_i);
    let obs = obs_base + vec3(P.cam_x, P.cam_y, P.cam_z);
    let fwd = safe_normalize(-obs, vec3(0.0, 0.0, -1.0));
    let world_up = select(vec3(0.0, 0.0, 1.0), vec3(0.0, 1.0, 0.0), abs(sin_i) > 0.995);
    let right = safe_normalize(cross(fwd, world_up), vec3(1.0, 0.0, 0.0));
    let up = safe_normalize(cross(right, fwd), vec3(0.0, 1.0, 0.0));
    let ray = safe_normalize(fwd * z_obs + right * sx + up * sy, fwd);

    let l = cross(obs, ray);
    let b = length(l);
    if b < 1e-4 {
        return TraceHit3D(vec4(0.0, 0.0, 0.0, 0.0), ray);
    }
    let n_hat = safe_normalize(l, vec3(0.0, 0.0, 1.0));

    let p_vec = obs - ray * dot(obs, ray);
    let b_line = length(p_vec);
    if b_line < 1e-4 {
        return TraceHit3D(vec4(0.0, 0.0, 0.0, 0.0), ray);
    }

    let r_start = max(3.0 * b_line, 3.0 * b);
    let s = -sqrt(max(r_start * r_start - b_line * b_line, 0.0));
    let start_world = p_vec + ray * s;
    let ex = safe_normalize(start_world, vec3(0.0, 0.0, 1.0));
    let ey = safe_normalize(cross(n_hat, ex), vec3(0.0, 1.0, 0.0));

    let vr0_sq = max(orbit_vr_sq(r_start, b, r_s, r_c), 0.0);
    var p = select(-r_start * r_start / b, -sqrt(vr0_sq), vr0_sq > 0.0);
    var r = r_start;
    var phi = 0.0f;
    var turned = false;
    var n_cross = 0u;
    let r_cap = r_s * 0.99;
    let max_steps_nominal = i32(P.max_phi / P.dphi) + 1;
    let max_steps_cap = select(select(2200, 3200, P.quality_tier >= 1.0), 4600, P.quality_tier >= 1.5);
    let max_steps = min(max_steps_nominal, max_steps_cap);

    var w_prev = ex * r;
    for (var i = 0; i < max_steps; i++) {
        let s_rk = rk4(r, p, b, r_s, r_c, P.dphi);
        let rn = s_rk.x;
        let pn_rk4 = s_rk.y;
        let vr2n = max(orbit_vr_sq(rn, b, r_s, r_c), 0.0);
        let pn = select(-sqrt(vr2n), sqrt(vr2n), pn_rk4 >= 0.0);
        let phin = phi + P.dphi;
        let ren = sqrt(rn * rn + r_c * r_c);

        let ren_bad = (ren != ren) || (abs(ren) > 1e20);
        if ren_bad || ren < r_cap || rn < r_c * 0.01 {
            return TraceHit3D(vec4(0.0, 0.0, 0.0, 0.0), ray);
        }
        if !turned && p < 0.0 && pn >= 0.0 {
            turned = true;
        }
        let w_new = ex * (rn * cos(phin)) + ey * (rn * sin(phin));
        if turned && rn >= r_start * 0.99 {
            let sky_dir = safe_normalize(w_new - w_prev, ray);
            return TraceHit3D(vec4(0.0, phin, 0.0, 2.0), sky_dir);
        }
        if w_prev.z * w_new.z <= 0.0 {
            let dz = w_new.z - w_prev.z;
            let t = select(0.0, clamp(-w_prev.z / dz, 0.0, 1.0), abs(dz) > 1e-6);
            let r_x = r + t * (rn - r);
            let re_x = sqrt(r_x * r_x + r_c * r_c);
            n_cross += 1u;
            if re_x >= P.disk_in && re_x <= P.disk_out {
                return TraceHit3D(vec4(re_x, phin, f32(n_cross), 1.0), ray);
            }
        }

        r = rn;
        p = pn;
        phi = phin;
        w_prev = w_new;
    }

    if r >= r_start * 0.5 {
        let sky_dir = safe_normalize(w_prev, ray);
        return TraceHit3D(vec4(0.0, phi, 0.0, 2.0), sky_dir);
    }
    return TraceHit3D(vec4(0.0, 0.0, 0.0, 0.0), ray);
}

fn kerr_mass(r_s: f32) -> f32 {
    return 0.5 * r_s;
}
fn kerr_a(r_s: f32, a_star: f32) -> f32 {
    return a_star * kerr_mass(r_s);
}
fn kerr_sigma(r: f32, th: f32, a: f32) -> f32 {
    let c = cos(th);
    r * r + a * a * c * c
}
fn kerr_delta(r: f32, r_s: f32, a: f32) -> f32 {
    r * r - r_s * r + a * a
}

fn trace_kerr3d(sx: f32, sy: f32) -> TraceHit3D {
    let r_s = P.r_s;
    let a = kerr_a(r_s, clamp(P.kerr_astar, -0.999, 0.999));
    let m = kerr_mass(r_s);
    let r_plus = m + sqrt(max(m * m - a * a, 0.0));
    let sin_inc = clamp(P.sin_inc, -1.0, 1.0);
    let theta_obs = clamp(asin(sin_inc), 1e-4, 3.14159265359 - 1e-4);

    let z_obs = 60.0 * r_s;
    let cos_i = sqrt(max(1.0 - sin_inc * sin_inc, 0.0));
    let caz = cos(P.az);
    let saz = sin(P.az);
    let obs_base = vec3(z_obs * sin_inc * caz, z_obs * sin_inc * saz, z_obs * cos_i);
    let obs = obs_base + vec3(P.cam_x, P.cam_y, P.cam_z);
    let fwd = safe_normalize(-obs, vec3(0.0, 0.0, -1.0));
    let world_up = select(vec3(0.0, 0.0, 1.0), vec3(0.0, 1.0, 0.0), abs(sin_inc) > 0.995);
    let right = safe_normalize(cross(fwd, world_up), vec3(1.0, 0.0, 0.0));
    let up = safe_normalize(cross(right, fwd), vec3(0.0, 1.0, 0.0));
    let ray = safe_normalize(fwd * z_obs + right * sx + up * sy, fwd);

    // Use Kerr image constants from camera-plane alpha/beta.
    let bx = sx;
    let by = sy;
    let xi = -bx * sin_inc;
    let c2 = max(1.0 - sin_inc * sin_inc, 0.0);
    let eta = by * by + (bx * bx - a * a) * c2;

    let b = sqrt(bx * bx + by * by);
    let r_start = max(max(40.0 * r_s, 12.0 * b), 20.0);

    var r = r_start;
    var th = theta_obs;
    var ph = 0.0;
    var sgn_r = -1.0;
    // Keep Kerr polar branch consistent with CPU/CUDA tracer conventions:
    // +beta (screen-up) maps to decreasing theta in this basis.
    var sgn_th = select(1.0, -1.0, by >= 0.0);
    var n_cross = 0u;
    let max_steps_nominal = i32(P.max_phi / max(P.dphi, 1e-4)) + 1;
    let max_steps_cap = select(select(2200, 3200, P.quality_tier >= 1.0), 4600, P.quality_tier >= 1.5);
    let max_steps = min(max_steps_nominal, max_steps_cap);

    var w_prev = vec3(r * sin(th) * cos(ph), r * sin(th) * sin(ph), r * cos(th));
    for (var i = 0; i < max_steps; i++) {
        let sig = max(kerr_sigma(r, th, a), 1e-12);
        let del = kerr_delta(r, r_s, a);
        let s = sin(th);
        let c = cos(th);
        let s2 = max(s * s, 1e-9);
        let p = (r * r + a * a) - a * xi;
        let rpot = p * p - del * ((xi - a) * (xi - a) + eta);
        let cot2 = (c * c) / s2;
        let tpot = eta + a * a * c * c - xi * xi * cot2;

        if rpot <= 1e-12 { sgn_r = -sgn_r; }
        if tpot <= 1e-12 { sgn_th = -sgn_th; }

        let rr = sqrt(max(rpot, 0.0));
        let thh = sqrt(max(tpot, 0.0));
        let dr = sgn_r * rr / sig;
        let dth = sgn_th * thh / sig;
        let dph = select(
            ((xi / s2) - a) / sig,
            ((xi / s2) - a + a * p / del) / sig,
            abs(del) >= 1e-9
        );

        let r_new = r + P.dphi * dr;
        let th_new = clamp(th + P.dphi * dth, 1e-4, 3.14159265359 - 1e-4);
        let ph_new = ph + P.dphi * dph;
        let bad = (r_new != r_new) || (th_new != th_new) || (ph_new != ph_new)
            || (abs(r_new) > 1e20) || (abs(ph_new) > 1e20);
        if bad || r_new <= 1.001 * r_plus {
            return TraceHit3D(vec4(0.0, 0.0, 0.0, 0.0), ray);
        }

        let w_new = vec3(
            r_new * sin(th_new) * cos(ph_new),
            r_new * sin(th_new) * sin(ph_new),
            r_new * cos(th_new)
        );

        if sgn_r > 0.0 && r_new >= 0.995 * r_start {
            let sky_dir = safe_normalize(w_new - w_prev, ray);
            return TraceHit3D(vec4(0.0, abs(ph_new), 0.0, 2.0), sky_dir);
        }

        let pi2 = 0.5 * 3.14159265359;
        if (th - pi2) * (th_new - pi2) <= 0.0 {
            n_cross += 1u;
            let t = clamp((pi2 - th) / (th_new - th + 1e-12), 0.0, 1.0);
            let r_cross = r + t * (r_new - r);
            if r_cross >= P.disk_in && r_cross <= P.disk_out {
                return TraceHit3D(vec4(r_cross, abs(ph_new), f32(n_cross), 1.0), ray);
            }
        }

        r = r_new;
        th = th_new;
        ph = ph_new;
        w_prev = w_new;
    }

    if (sgn_r > 0.0 && r > 0.5 * r_start) || r > max(3.0 * r_s, 1.2 * r_plus) {
        let sky_dir = safe_normalize(w_prev, ray);
        return TraceHit3D(vec4(0.0, abs(ph), 0.0, 2.0), sky_dir);
    }
    return TraceHit3D(vec4(0.0, 0.0, 0.0, 0.0), ray);
}

// ── Interior tracer ───────────────────────────────────────────────────────────
// Returns vec4(r_eff_hit, phi_total_or_hit, f32(n_cross), kind)
//   kind 0 = captured (fell back to core floor)
//   kind 1 = surface hit (disk or core shell, depending on mode)
//   kind 2 = escaped to outside sky

fn trace_interior(bx_in: f32, by_in: f32, core_look: bool) -> vec4<f32> {
    let ca = cos(P.az); let sa = sin(P.az);
    let bx = ca * bx_in - sa * by_in;
    let by = sa * bx_in + ca * by_in;
    let b = sqrt(bx * bx + by * by);
    let r_s = P.r_s;
    let r_c = P.r_c;
    let r_h = sqrt(max(r_s * r_s - r_c * r_c, 1e-6));
    let r_cam = clamp(P.r_cam_frac, 0.05, 0.99) * r_h;

    if b < 1e-12 {
        if core_look {
            return vec4(r_c, 0.0, 1.0, 1.0);
        }
        return vec4(0.0, P.max_phi, 0.0, 2.0);
    }

    let sin_i = by / b;
    let is_eq = abs(sin_i) < 1e-6;
    let vr0_sq = max(orbit_vr_sq(r_cam, b, r_s, r_c), 0.0);
    var p = select(sqrt(vr0_sq), -sqrt(vr0_sq), core_look); // outward vs inward branch
    var r = r_cam;
    var phi = 0.0f;
    var n_cross = 0u;
    var turned = false;

    let re_cam = sqrt(r_cam * r_cam + r_c * r_c);
    let re_cap = re_cam * 1.05;
    let re_core_cap = max(1.02 * r_c, r_c + 1e-6);
    let r_escape = max(max(3.0 * b, P.disk_out * 1.5), 20.0 * r_s);
    let max_steps_nominal = i32(P.max_phi / P.dphi) + 1;
    let max_steps_cap = select(select(1800, 2800, P.quality_tier >= 1.0), 4200, P.quality_tier >= 1.5);
    let max_steps = min(max_steps_nominal, max_steps_cap);

    for (var i = 0; i < max_steps; i++) {
        let s      = rk4(r, p, b, r_s, r_c, P.dphi);
        let rn     = s.x;
        let pn_rk4 = s.y;
        let vr2n   = max(orbit_vr_sq(rn, b, r_s, r_c), 0.0);
        let pn     = select(-sqrt(vr2n), sqrt(vr2n), pn_rk4 >= 0.0);
        let phin   = phi + P.dphi;
        let ren    = sqrt(rn * rn + r_c * r_c);

        if core_look {
            if ren <= re_core_cap || rn <= r_c * 0.01 {
                return vec4(ren, phin, 1.0, 1.0);
            }
            if p < 0.0 && pn >= 0.0 && rn > r_cam * 1.01 {
                return vec4(0.0, phin, 0.0, 2.0);
            }
        } else {
            if !turned && p > 0.0 && pn <= 0.0 { turned = true; }
            if turned && ren < re_cap { return vec4(0.0, phin, 0.0, 0.0); }
            if !turned && rn >= r_escape { return vec4(0.0, phin, 0.0, 2.0); }

            if is_eq {
                let re_cur = sqrt(r * r + r_c * r_c);
                if !turned && re_cur >= P.disk_in && re_cur <= P.disk_out && p > 0.0 {
                    return vec4(re_cur, phi, 1.0, 1.0);
                }
            } else {
                let phi_target = (f32(n_cross) + 1.0) * 3.14159265359;
                if phi < phi_target && phin >= phi_target {
                    let t = (phi_target - phi) / P.dphi;
                    let r_x = r + t * (rn - r);
                    let re_x = sqrt(r_x * r_x + r_c * r_c);
                    n_cross += 1u;
                    if re_x >= P.disk_in && re_x <= P.disk_out {
                        return vec4(re_x, phi_target, f32(n_cross), 1.0);
                    }
                }
            }
        }

        r = rn;
        p = pn;
        phi = phin;
    }

    if core_look {
        return vec4(re_core_cap, phi, 1.0, 1.0);
    }
    if r >= r_escape * 0.5 && !turned { return vec4(0.0, phi, 0.0, 2.0); }
    return vec4(0.0, phi, 0.0, 0.0);
}

// ── Disk colour: covariant transfer (CPU/CUDA parity target) ─────────────────

fn disk_transfer_factor(r_eff: f32, r_s: f32, bx_raw: f32, phi_orb: f32, sin_inc: f32) -> f32 {
    let r_safe = max(r_eff, 1e-6);
    let g_gr = sqrt(max(1.0 - r_s / r_safe, 0.0));
    let mut beta = min(sqrt(max(r_s / (2.0 * r_safe), 0.0)), 0.7);
    let mut omega_boost = 0.0;
    if (P.kerr_enable > 0.5) {
        // Frame-drag proxy for Kerr emissive flow (strong near horizon, fades outward).
        omega_boost = 0.35 * abs(P.kerr_astar) * pow(clamp(r_s / r_safe, 0.0, 1.0), 1.5);
        beta = min(beta * (1.0 + omega_boost), 0.85);
    }
    let beta_eff = clamp(beta, -0.88, 0.88);
    let gamma = 1.0 / sqrt(max(1.0 - beta_eff * beta_eff, 1e-6));
    // Match bh_render’s LOS projection model: dominant orbital phase + small
    // screen-space blend for continuity.
    let mu_phi = clamp(sin(phi_orb), -1.0, 1.0);
    let mu_screen = clamp(bx_raw / r_safe, -1.0, 1.0);
    let mu = clamp(0.8 * mu_phi + 0.2 * mu_screen, -1.0, 1.0);
    let beta_obs = beta_eff * sin_inc * mu;
    let g_dop = 1.0 / (gamma * (1.0 - beta_obs));
    let g = g_gr * g_dop;
    let g4 = clamp(pow(g, 4.0), 1e-6, 300.0);
    return g4 * (1.0 + 0.15 * omega_boost);
}

fn disk_color(r_eff: f32, phi_orb: f32, n_cross: u32, bx_raw: f32, bg_stars: vec3<f32>) -> vec3<f32> {
    let r_isco = 3.0 * P.r_s;
    let t_rel = pow(max(r_isco / max(r_eff, 1e-6), 1e-6), 0.75);
    let excess = max(r_eff - P.disk_out, 0.0) / max(0.5 * P.disk_out, 1e-6);
    let outer_taper = exp(-excess * excess);
    let fade = pow(0.65, f32(max(i32(n_cross) - 1, 0)));
    let transfer = disk_transfer_factor(r_eff, P.r_s, bx_raw, phi_orb, P.sin_inc);
    let source_local = max(t_rel * fade * outer_taper, 0.0);
    let g_cov = pow(max(transfer, 1e-9), 0.25);
    let alpha_base = 0.35 * max(P.tau_scale, 0.0) * (1.0 + f32(max(i32(n_cross) - 1, 0)) * 0.4);

    var luminance = 0.0;
    if (P.use_transfer > 0.5) {
        let steps = 8;
        let path_scale = max(r_eff / max(P.r_s, 1e-6), 1e-6);
        var intensity = 0.0;
        for (var si = 0; si < steps; si++) {
            let u = (f32(si) + 0.5) / f32(steps);
            let local_mod = 1.0 + 0.20 * (1.0 - u);
            let j_obs = max(source_local * local_mod, 0.0) * g_cov * g_cov * g_cov;
            let alpha_obs = max(alpha_base * (0.7 + 0.6 * u), 0.0) * g_cov;
            let tau_seg = max(alpha_obs * path_scale / f32(steps), 0.0);
            let source_fn = select(j_obs, j_obs / alpha_obs, alpha_obs > 1e-12);
            let e = exp(-tau_seg);
            intensity = intensity * e + source_fn * (1.0 - e);
        }
        luminance = max(intensity, 0.0);
    } else {
        luminance = max(source_local * transfer, 0.0);
    }

    let b = luminance / (1.0 + luminance);
    var disk = vec3(
        clamp(pow(b, 0.35), 0.0, 1.0),
        clamp((210.0 / 255.0) * pow(b, 0.60), 0.0, 1.0),
        clamp((130.0 / 255.0) * pow(b, 1.60), 0.0, 1.0)
    );

    if (P.disk_model > 0.5) {
        let tau = (0.45 * max(P.tau_scale, 0.0))
            * pow(P.r_s / max(r_eff, 1e-6), 0.7)
            * (1.0 + f32(max(i32(n_cross) - 1, 0)) * 0.25);
        let trans = clamp(exp(-tau), 0.0, 1.0);
        disk = min(disk * 2.5, vec3(1.0));
        return mix(disk, bg_stars, trans);
    }

    return disk;
}

// Optically thin volumetric RIAF glow for escaped rays.
// This avoids the "flat half-plane" look by letting diffuse plasma emission
// contribute on both sides of the hole while still preserving background stars.
fn riaf_volume_color(sx_raw: f32, sy_raw: f32, bg_stars: vec3<f32>) -> vec3<f32> {
    let sx = cos(P.az) * sx_raw - sin(P.az) * sy_raw;
    let sy = sin(P.az) * sx_raw + cos(P.az) * sy_raw;
    let b = sqrt(sx * sx + (sy * P.sin_inc) * (sy * P.sin_inc));
    let r_s = max(P.r_s, 1e-6);
    let q = clamp(P.quality_tier, 0.0, 2.0);
    let steps = select(select(6u, 10u, q >= 1.0), 14u, q >= 1.5);
    let u_max = select(select(16.0 * r_s, 24.0 * r_s, q >= 1.0), 32.0 * r_s, q >= 1.5);
    let h = 0.45 + 0.35 * (1.0 - P.sin_inc);
    var i_obs = vec3(0.0);
    var trans = 1.0;
    for (var i = 0u; i < 20u; i++) {
        if (i >= steps) {
            break;
        }
        let fu = (f32(i) + 0.5) / max(f32(steps), 1.0);
        let u = mix(-u_max, u_max, fu);
        let r_eff = sqrt(b * b + u * u + P.r_c * P.r_c);

        let r_peak = 6.0 * r_s;
        let r_width = 7.2 * r_s;
        let radial_env = exp(-pow((r_eff - r_peak) / max(r_width, 1e-6), 2.0));
        let z_los = sy + 0.17 * u * (1.0 - P.sin_inc);
        let vertical_env = exp(-pow(z_los / max(h * r_s, 1e-6), 2.0));
        let turb = 0.78 + 0.22 * hash31(vec3(floor(vec2(sx, sy) * 11.0), f32(i) * 3.7));
        let transfer = pow(max(disk_transfer_factor(max(r_eff, 3.0 * r_s), r_s, sx, P.sin_inc), 1e-9), 0.58);
        let source = max(0.30 * radial_env * vertical_env * transfer * turb, 0.0);
        let bmap = source / (1.0 + source);
        let glow = vec3(
            clamp(pow(bmap, 0.38), 0.0, 1.0),
            clamp((220.0 / 255.0) * pow(bmap, 0.62), 0.0, 1.0),
            clamp((120.0 / 255.0) * pow(bmap, 1.55), 0.0, 1.0)
        );

        let path = u_max / max(f32(steps) * r_s, 1e-6);
        let tau_seg = max(P.tau_scale, 0.0)
            * (0.12 + 0.10 * radial_env)
            * radial_env
            * (0.45 + 0.55 * vertical_env)
            * path;
        let e = clamp(exp(-max(tau_seg, 0.0)), 0.0, 1.0);
        i_obs = i_obs * e + glow * (1.0 - e);
        trans *= e;
        if trans < 0.02 {
            break;
        }
    }
    return mix(i_obs, bg_stars, clamp(trans, 0.0, 1.0));
}

// ── Vertex / Fragment ─────────────────────────────────────────────────────────

struct VOut { @builtin(position) pos: vec4<f32> }

@vertex
fn vs(@builtin(vertex_index) vi: u32) -> VOut {
    var xy = array<vec2<f32>, 6>(
        vec2(-1.0, -1.0), vec2(1.0, -1.0), vec2(-1.0, 1.0),
        vec2(-1.0,  1.0), vec2(1.0, -1.0), vec2(1.0,  1.0),
    );
    return VOut(vec4(xy[vi], 0.0, 1.0));
}

fn shade_sample(sx: f32, sy: f32) -> vec3<f32> {
    // Camera roll rotates the screen-space sampling basis.
    let cr = cos(P.cam_roll);
    let sr = sin(P.cam_roll);
    let sxr = cr * sx - sr * sy;
    let syr = sr * sx + cr * sy;
    let interior = P.interior_mode > 0.5;
    let core_look = P.core_look_mode > 0.5;
    var hit: vec4<f32>;
    var sky_dir = safe_normalize(vec3(sxr, syr, P.r_s), vec3(0.0, 0.0, 1.0));
    if interior {
        let bx = sxr;
        let by = syr * P.sin_inc;
        hit = trace_interior(bx, by, core_look);
    } else {
        let t3 = if P.kerr_enable > 0.5 {
            trace_kerr3d(sxr, syr)
        } else {
            trace_true3d(sxr, syr)
        };
        hit = t3.hit;
        sky_dir = t3.sky_dir;
    }
    let kind = hit.w;

    // Shadow — pure black
    if kind < 0.5 {
        if interior && !core_look {
            // Interior floor glow for returned (non-escaping) outward rays.
            let bx = sxr;
            let by = syr * P.sin_inc;
            let b = sqrt(bx * bx + by * by);
            let b_crit = 1.5 * sqrt(3.0) * P.r_s;
            let x = clamp((b / max(b_crit, 1e-6) - 1.0) * 2.4 + 0.5, 0.0, 1.0);
            return vec3(0.20 + 0.70 * x, 0.10 + 0.45 * x, 0.03 + 0.18 * x);
        }
        return vec3(0.0, 0.0, 0.0);
    }

    // Disk hit — Novikov-Thorne + Doppler + gravitational redshift
    // For Doppler, use the x-component in the disk frame (rotated by azimuth),
    // so the bright crescent stays on the correct side when azimuth is changed.
    if kind < 1.5 {
        if interior && core_look {
            // Pure-physics-style core palette from traced invariants.
            let rr = clamp(hit.x / max(P.r_c * 2.4, 1e-6), 0.0, 1.0);
            let swirl = 0.5 + 0.5 * sin(8.0 * hit.y + 2.0 * sxr / max(P.r_s, 1e-6));
            return vec3(0.22 + 0.63 * (1.0 - rr), 0.12 + 0.42 * (1.0 - rr), 0.04 + 0.16 * swirl * (1.0 - rr));
        }
        let sx_disk = cos(P.az) * sxr - sin(P.az) * syr;
        let stars = starfield_from_dir(sky_dir);
        return disk_color(hit.x, hit.y, u32(hit.z), sx_disk, stars);
    }

    // Escaped photon — background with gravitationally lensed star field
    //
    // The photon swept phi_total radians in its orbital plane.
    // A straight photon (no BH) sweeps exactly π. The excess delta = phi_total − π
    // is the deflection angle. We reverse-rotate the sky direction by delta so
    // stars appear at their true (source) positions, giving lensing arcs near
    // the photon sphere automatically.
    let stars = starfield_from_dir(sky_dir);
    if P.disk_model > 0.5 && P.riaf_volume > 0.5 {
        return riaf_volume_color(sxr, syr, stars);
    }
    return stars;
}

@fragment
fn fs(@builtin(position) frag: vec4<f32>) -> @location(0) vec4<f32> {
    let min_dim = max(min(P.width, P.height), 1.0);
    let scale = 2.0 * P.fov * P.r_s / min_dim;
    var frag_xy = frag.xy;
    if (P.quality_tier < 0.5) {
        // Quarter-rate shading on low tier (2x2 pixel blocks share one sample).
        frag_xy = floor(frag.xy * 0.5) * 2.0 + vec2(1.0, 1.0);
    }

    var col = vec3(0.0, 0.0, 0.0);
    if (P.quality_tier < 0.5) {
        // Low-memory mode: one sample/pixel.
        col = shade_sample(
            (frag_xy.x - P.width * 0.5) * scale,
            (P.height * 0.5 - frag_xy.y) * scale
        );
    } else if (P.quality_tier < 1.5) {
        // Medium mode: single sample (keeps controls snappy at 1080p+).
        col = shade_sample(
            (frag.x - P.width * 0.5) * scale,
            (P.height * 0.5 - frag.y) * scale
        );
    } else {
        // High mode: 2-sample diagonal supersampling (not 2x2) for speed.
        let o = 0.25;
        let c00 = shade_sample((frag.x - o - P.width * 0.5) * scale, (P.height * 0.5 - (frag.y - o)) * scale);
        let c11 = shade_sample((frag.x + o - P.width * 0.5) * scale, (P.height * 0.5 - (frag.y + o)) * scale);
        col = (c00 + c11) * 0.5;
    }

    // Filmic-ish tone map + display gamma.
    col = col / (vec3(1.0) + col);
    col = pow(col, vec3(1.0 / 2.2));
    return vec4(col, 1.0);
}
"#;

// ── CPU-side uniform struct (must match WGSL layout, std140 / 16-byte align) ──

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    r_s: f32,
    r_c: f32,
    disk_in: f32,
    disk_out: f32,
    sin_inc: f32,
    fov: f32,
    width: f32,
    height: f32,
    max_phi: f32,
    dphi: f32,
    az: f32,
    cam_x: f32,
    cam_y: f32,
    cam_z: f32,
    cam_roll: f32,
    interior_mode: f32,
    core_look_mode: f32,
    r_cam_frac: f32,
    quality_tier: f32,
    use_transfer: f32,
    tau_scale: f32,
    disk_model: f32,
    local_stars: f32,
    riaf_volume: f32,
    kerr_enable: f32,
    kerr_astar: f32,
    _pad1: f32,
    _pad2: f32,
    _pad3: f32,
}

// ── Camera state ──────────────────────────────────────────────────────────────

struct Camera {
    inclination: f32,     // degrees, 0 = face-on, 90 = edge-on
    azimuth: f32,         // radians, rotates disk orientation on screen
    fov_rs: f32,          // half-width of image in r_s
    disk_outer: f32,      // disk outer radius in r_s
    gutoe_core: bool,     // toggle GUTOE lattice correction
    interior_mode: bool,  // camera placed inside horizon
    core_look_mode: bool, // interior camera points at lattice core
    r_cam_frac: f32,      // r_cam = r_cam_frac * r_horizon (coordinate)
    cam_x: f32,           // freecam world-space x offset
    cam_y: f32,           // freecam world-space y offset
    cam_z: f32,           // freecam world-space z offset
    cam_roll: f32,        // camera roll in radians
    use_transfer: bool,
    tau_scale: f32,
    riaf_mode: bool,
    local_stars: bool,
    riaf_volume: bool,
    kerr_enable: bool,
    kerr_astar: f32,
}

impl Default for Camera {
    fn default() -> Self {
        let use_transfer = std::env::var("BH_USE_TRANSFER").ok().is_some_and(|s| {
            matches!(
                s.as_str(),
                "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON"
            )
        });
        let tau_scale = std::env::var("BH_TAU_SCALE")
            .ok()
            .and_then(|s| s.parse::<f32>().ok())
            .unwrap_or(1.0)
            .max(0.0);
        let riaf_mode = std::env::var("BH_DISK_MODEL")
            .ok()
            .map(|s| s.eq_ignore_ascii_case("riaf") || s.eq_ignore_ascii_case("volumetric"))
            .unwrap_or(false);
        let local_stars = std::env::var("BH_LOCAL_STARS")
            .ok()
            .map(|s| !matches!(s.as_str(), "0" | "false" | "FALSE" | "off" | "OFF"))
            .unwrap_or(true);
        let riaf_volume = std::env::var("BH_RIAF_VOLUME")
            .ok()
            .map(|s| !matches!(s.as_str(), "0" | "false" | "FALSE" | "off" | "OFF"))
            .unwrap_or(true);
        let kerr_astar = std::env::var("BH_KERR_ASTAR")
            .ok()
            .and_then(|s| s.parse::<f32>().ok())
            .unwrap_or(0.0)
            .clamp(-0.999, 0.999);
        let kerr_enable = kerr_astar.abs() > 1e-6;
        Self {
            // EHT M87 geometry: ~17° from face-on, bright crescent at bottom.
            // azimuth = -π/2 rotates the disk so the approaching (bright) side is
            // at screen-bottom, matching the classic EHT image orientation.
            inclination: 17.0,
            azimuth: -std::f32::consts::FRAC_PI_2,
            fov_rs: 7.0,
            disk_outer: 20.0,
            gutoe_core: true,
            interior_mode: false,
            core_look_mode: false,
            r_cam_frac: 0.72,
            cam_x: 0.0,
            cam_y: 0.0,
            cam_z: 0.0,
            cam_roll: 0.0,
            use_transfer,
            tau_scale,
            riaf_mode,
            local_stars,
            riaf_volume,
            kerr_enable,
            kerr_astar,
        }
    }
}

impl Camera {
    fn move_local(&mut self, forward: f32, right: f32, up: f32, speed: f32) {
        let inc = self.inclination.to_radians();
        let (sinc, cinc) = inc.sin_cos();
        let (saz, caz) = self.azimuth.sin_cos();
        let fwd = [-sinc * caz, -sinc * saz, -cinc];
        let mut right_v = [fwd[1], -fwd[0], 0.0];
        let rn = (right_v[0] * right_v[0] + right_v[1] * right_v[1] + right_v[2] * right_v[2])
            .sqrt()
            .max(1e-6);
        right_v[0] /= rn;
        right_v[1] /= rn;
        right_v[2] /= rn;
        let mut up_v = [
            right_v[1] * fwd[2] - right_v[2] * fwd[1],
            right_v[2] * fwd[0] - right_v[0] * fwd[2],
            right_v[0] * fwd[1] - right_v[1] * fwd[0],
        ];
        let (sr, cr) = self.cam_roll.sin_cos();
        let rr = [
            right_v[0] * cr + up_v[0] * sr,
            right_v[1] * cr + up_v[1] * sr,
            right_v[2] * cr + up_v[2] * sr,
        ];
        up_v = [
            -right_v[0] * sr + up_v[0] * cr,
            -right_v[1] * sr + up_v[1] * cr,
            -right_v[2] * sr + up_v[2] * cr,
        ];
        self.cam_x += (fwd[0] * forward + rr[0] * right + up_v[0] * up) * speed;
        self.cam_y += (fwd[1] * forward + rr[1] * right + up_v[1] * up) * speed;
        self.cam_z += (fwd[2] * forward + rr[2] * right + up_v[2] * up) * speed;
    }

    fn params(&self, width: f32, height: f32, quality_tier: f32) -> Params {
        // r_s = 1 in internal units; r_core = sqrt(C_inf) * l_P
        const C_INF: f32 = 0.5466;
        let r_c = if self.gutoe_core { C_INF.sqrt() } else { 0.0 };
        // Live-viewer budgets: prioritize responsiveness over offline-grade convergence.
        let (max_phi, dphi) = if quality_tier < 0.5 {
            (8.0 * std::f32::consts::PI, 0.030)
        } else if quality_tier < 1.5 {
            (12.0 * std::f32::consts::PI, 0.018)
        } else {
            (18.0 * std::f32::consts::PI, 0.012)
        };
        Params {
            r_s: 1.0,
            r_c,
            disk_in: 3.0, // r_ISCO = 3 r_s
            disk_out: self.disk_outer,
            sin_inc: self.inclination.to_radians().sin(),
            fov: self.fov_rs,
            width,
            height,
            max_phi,
            dphi,
            az: self.azimuth,
            cam_x: self.cam_x,
            cam_y: self.cam_y,
            cam_z: self.cam_z,
            cam_roll: self.cam_roll,
            interior_mode: if self.interior_mode { 1.0 } else { 0.0 },
            core_look_mode: if self.core_look_mode { 1.0 } else { 0.0 },
            r_cam_frac: self.r_cam_frac,
            quality_tier,
            use_transfer: if self.use_transfer { 1.0 } else { 0.0 },
            tau_scale: self.tau_scale,
            disk_model: if self.riaf_mode { 1.0 } else { 0.0 },
            local_stars: if self.local_stars { 1.0 } else { 0.0 },
            riaf_volume: if self.riaf_volume { 1.0 } else { 0.0 },
            kerr_enable: if self.kerr_enable { 1.0 } else { 0.0 },
            kerr_astar: self.kerr_astar,
            _pad1: 0.0,
            _pad2: 0.0,
            _pad3: 0.0,
        }
    }
}

// ── wgpu resources ────────────────────────────────────────────────────────────

struct Gpu {
    _instance: wgpu::Instance,
    surface: wgpu::Surface<'static>,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    pipeline: wgpu::RenderPipeline,
    bind_group: wgpu::BindGroup,
    uni_buf: wgpu::Buffer,
    _star_tex: wgpu::Texture,
    _star_view: wgpu::TextureView,
    _star_samp: wgpu::Sampler,
    quality_tier: f32,
}

fn detect_quality_tier(adapter_info: &wgpu::AdapterInfo) -> f32 {
    if let Ok(force) = std::env::var("BH_VIEWER_QUALITY") {
        let v = force.to_ascii_lowercase();
        return match v.as_str() {
            "low" => 0.0,
            "medium" | "med" => 1.0,
            "high" => 2.0,
            _ => 2.0,
        };
    }
    if adapter_info.backend == wgpu::Backend::Metal
        && adapter_info.device_type == wgpu::DeviceType::IntegratedGpu
    {
        let n = adapter_info.name.to_ascii_lowercase();
        if n.contains("m1") {
            return 0.0; // safest default for 8 GB M1 class machines
        }
        return 1.0;
    }
    2.0
}

impl Gpu {
    fn new(window: Arc<Window>) -> Self {
        pollster::block_on(Self::new_async(window))
    }

    async fn new_async(window: Arc<Window>) -> Self {
        let size = window.inner_size();
        let power_pref = if cfg!(target_os = "macos") {
            wgpu::PowerPreference::LowPower
        } else {
            wgpu::PowerPreference::HighPerformance
        };
        let try_backend = |backends: wgpu::Backends, window: Arc<Window>| async move {
            let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
                backends,
                ..Default::default()
            });
            let surface = instance.create_surface(window).ok()?;
            let adapter = instance
                .request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference: power_pref,
                    compatible_surface: Some(&surface),
                    force_fallback_adapter: false,
                })
                .await?;
            Some((instance, surface, adapter))
        };

        let backend_override = std::env::var("BH_BACKEND")
            .ok()
            .map(|v| v.to_ascii_lowercase());
        let primary = match backend_override.as_deref() {
            Some("gl") => wgpu::Backends::GL,
            Some("vulkan") => wgpu::Backends::VULKAN,
            Some("metal") => wgpu::Backends::METAL,
            Some("dx12") => wgpu::Backends::DX12,
            Some("primary") => wgpu::Backends::PRIMARY,
            _ => {
                if cfg!(target_os = "linux") {
                    wgpu::Backends::VULKAN
                } else if cfg!(target_os = "windows") {
                    wgpu::Backends::DX12
                } else if cfg!(target_os = "macos") {
                    wgpu::Backends::METAL
                } else {
                    wgpu::Backends::PRIMARY
                }
            }
        };

        let primary_init = try_backend(primary, Arc::clone(&window)).await;
        let fallback_init =
            if primary_init.is_none() && cfg!(any(target_os = "linux", target_os = "windows")) {
                log::warn!("Primary backend init failed; retrying with OpenGL backend");
                try_backend(wgpu::Backends::GL, Arc::clone(&window)).await
            } else {
                None
            };
        let (instance, surface, adapter) = primary_init
            .or(fallback_init)
            .expect("no GPU adapter found — backend init failed");

        let adapter_info = adapter.get_info();
        let quality_tier = detect_quality_tier(&adapter_info);
        log::info!(
            "Adapter={} backend={:?} type={:?} quality_tier={}",
            adapter_info.name,
            adapter_info.backend,
            adapter_info.device_type,
            quality_tier
        );

        let caps = surface.get_capabilities(&adapter);
        let format = caps
            .formats
            .iter()
            .copied()
            .find(|f| f.is_srgb())
            .unwrap_or(caps.formats[0]);

        // Stability first: some NVIDIA/Linux stacks lose device on AutoNoVsync.
        // Default to FIFO unless explicitly overridden.
        let present_mode = match std::env::var("BH_PRESENT_MODE")
            .ok()
            .map(|s| s.to_ascii_lowercase())
            .as_deref()
        {
            Some("immediate") if caps.present_modes.contains(&wgpu::PresentMode::Immediate) => {
                wgpu::PresentMode::Immediate
            }
            Some("mailbox") if caps.present_modes.contains(&wgpu::PresentMode::Mailbox) => {
                wgpu::PresentMode::Mailbox
            }
            Some("auto") if caps.present_modes.contains(&wgpu::PresentMode::AutoNoVsync) => {
                wgpu::PresentMode::AutoNoVsync
            }
            _ => wgpu::PresentMode::Fifo,
        };
        let base_config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format,
            width: size.width.max(1),
            height: size.height.max(1),
            present_mode,
            alpha_mode: caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: if quality_tier < 1.0 { 2 } else { 2 },
        };

        // Some drivers report a stale lost device right at startup under load.
        // Retry device creation + surface configure before giving up.
        let mut configured = None;
        for attempt in 1..=4 {
            let (device, queue) = adapter
                .request_device(&wgpu::DeviceDescriptor::default(), None)
                .await
                .expect("failed to create wgpu device");
            let config = base_config.clone();
            let ok = panic::catch_unwind(AssertUnwindSafe(|| {
                surface.configure(&device, &config);
            }))
            .is_ok();
            if ok {
                configured = Some((device, queue, config));
                break;
            }
            log::warn!(
                "surface.configure failed on startup attempt {attempt}/4 (device lost); retrying"
            );
            std::thread::sleep(Duration::from_millis(120));
        }
        let (device, queue, config) =
            configured.expect("wgpu surface/device init failed after retries");

        // Uniform buffer
        let uni_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("params"),
            size: std::mem::size_of::<Params>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // Optional real-sky starmap texture.
        let (star_tex, star_view) = if let Ok(path) = std::env::var("BH_STARMAP_PATH") {
            if let Ok(img) = image::open(&path).map(|i| i.to_rgba8()) {
                let (tw, th) = img.dimensions();
                let tex = device.create_texture(&wgpu::TextureDescriptor {
                    label: Some("starmap"),
                    size: wgpu::Extent3d {
                        width: tw,
                        height: th,
                        depth_or_array_layers: 1,
                    },
                    mip_level_count: 1,
                    sample_count: 1,
                    dimension: wgpu::TextureDimension::D2,
                    format: wgpu::TextureFormat::Rgba8UnormSrgb,
                    usage: wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_DST,
                    view_formats: &[],
                });
                queue.write_texture(
                    wgpu::ImageCopyTexture {
                        texture: &tex,
                        mip_level: 0,
                        origin: wgpu::Origin3d::ZERO,
                        aspect: wgpu::TextureAspect::All,
                    },
                    &img,
                    wgpu::ImageDataLayout {
                        offset: 0,
                        bytes_per_row: Some(4 * tw),
                        rows_per_image: Some(th),
                    },
                    wgpu::Extent3d {
                        width: tw,
                        height: th,
                        depth_or_array_layers: 1,
                    },
                );
                log::info!("viewer starmap loaded: {} ({}x{})", path, tw, th);
                let view = tex.create_view(&wgpu::TextureViewDescriptor::default());
                (tex, view)
            } else {
                log::warn!("viewer starmap failed to load: {}", path);
                let tex = device.create_texture(&wgpu::TextureDescriptor {
                    label: Some("starmap_fallback"),
                    size: wgpu::Extent3d {
                        width: 1,
                        height: 1,
                        depth_or_array_layers: 1,
                    },
                    mip_level_count: 1,
                    sample_count: 1,
                    dimension: wgpu::TextureDimension::D2,
                    format: wgpu::TextureFormat::Rgba8UnormSrgb,
                    usage: wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_DST,
                    view_formats: &[],
                });
                queue.write_texture(
                    wgpu::ImageCopyTexture {
                        texture: &tex,
                        mip_level: 0,
                        origin: wgpu::Origin3d::ZERO,
                        aspect: wgpu::TextureAspect::All,
                    },
                    &[0_u8, 0_u8, 0_u8, 0_u8],
                    wgpu::ImageDataLayout {
                        offset: 0,
                        bytes_per_row: Some(4),
                        rows_per_image: Some(1),
                    },
                    wgpu::Extent3d {
                        width: 1,
                        height: 1,
                        depth_or_array_layers: 1,
                    },
                );
                let view = tex.create_view(&wgpu::TextureViewDescriptor::default());
                (tex, view)
            }
        } else {
            let tex = device.create_texture(&wgpu::TextureDescriptor {
                label: Some("starmap_fallback"),
                size: wgpu::Extent3d {
                    width: 1,
                    height: 1,
                    depth_or_array_layers: 1,
                },
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format: wgpu::TextureFormat::Rgba8UnormSrgb,
                usage: wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_DST,
                view_formats: &[],
            });
            queue.write_texture(
                wgpu::ImageCopyTexture {
                    texture: &tex,
                    mip_level: 0,
                    origin: wgpu::Origin3d::ZERO,
                    aspect: wgpu::TextureAspect::All,
                },
                &[0_u8, 0_u8, 0_u8, 0_u8],
                wgpu::ImageDataLayout {
                    offset: 0,
                    bytes_per_row: Some(4),
                    rows_per_image: Some(1),
                },
                wgpu::Extent3d {
                    width: 1,
                    height: 1,
                    depth_or_array_layers: 1,
                },
            );
            let view = tex.create_view(&wgpu::TextureViewDescriptor::default());
            (tex, view)
        };
        let star_samp = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("starmap_sampler"),
            address_mode_u: wgpu::AddressMode::Repeat,
            address_mode_v: wgpu::AddressMode::ClampToEdge,
            address_mode_w: wgpu::AddressMode::ClampToEdge,
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            mipmap_filter: wgpu::FilterMode::Nearest,
            ..Default::default()
        });

        // Bind group layout
        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: true },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                    count: None,
                },
            ],
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("bg"),
            layout: &bgl,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: uni_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&star_view),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::Sampler(&star_samp),
                },
            ],
        });

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("bh"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("pl"),
            bind_group_layouts: &[&bgl],
            push_constant_ranges: &[],
        });
        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("bh_pipe"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: "vs",
                compilation_options: Default::default(),
                buffers: &[],
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: "fs",
                compilation_options: Default::default(),
                targets: &[Some(wgpu::ColorTargetState {
                    format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState::default(),
            depth_stencil: None,
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
            cache: None,
        });

        Self {
            _instance: instance,
            surface,
            device,
            queue,
            config,
            pipeline,
            bind_group,
            uni_buf,
            _star_tex: star_tex,
            _star_view: star_view,
            _star_samp: star_samp,
            quality_tier,
        }
    }

    fn resize(&mut self, w: u32, h: u32) {
        self.config.width = w.max(1);
        self.config.height = h.max(1);
        self.surface.configure(&self.device, &self.config);
    }

    fn upload_params(&self, params: &Params) {
        self.queue
            .write_buffer(&self.uni_buf, 0, bytemuck::bytes_of(params));
    }

    fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        let frame = self.surface.get_current_texture()?;
        let view = frame.texture.create_view(&Default::default());
        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut rp = enc.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("frame"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            rp.set_pipeline(&self.pipeline);
            rp.set_bind_group(0, &self.bind_group, &[]);
            rp.draw(0..6, 0..1); // fullscreen quad (2 triangles = 6 vertices)
        }
        self.queue.submit([enc.finish()]);
        frame.present();
        Ok(())
    }
}

// ── Application ───────────────────────────────────────────────────────────────

struct App {
    window: Option<Arc<Window>>,
    gpu: Option<Gpu>,
    camera: Camera,
    mouse_down: bool,
    mouse_right_down: bool,
    last_mouse: PhysicalPosition<f64>,
    win_size: (f32, f32), // (width, height) in physical pixels
    gilrs: Option<Gilrs>,
    pad_prev: PadPrev,
    held: HeldKeys,
    auto_spin: bool,
    quality_tier: f32,
    last_frame_at: Instant,
    min_frame_dt: Duration,
    dynamic_quality: bool,
    perf_window_start: Instant,
    perf_accum: Duration,
    perf_frames: u32,
    avg_frame_ms: f32,
    over_budget_windows: u32,
    under_budget_windows: u32,
    pad_tuning: PadTuning,
    pad_name: Option<String>,
    pad_dualsense: bool,
}

impl App {
    fn recreate_gpu(&mut self) {
        let Some(win) = self.window.as_ref() else {
            return;
        };
        log::warn!("Recreating GPU device/surface after loss");
        let gpu = Gpu::new(Arc::clone(win));
        self.apply_quality_tier(gpu.quality_tier);
        self.gpu = Some(gpu);
    }

    fn safe_resize(&mut self, w: u32, h: u32) {
        let Some(gpu) = self.gpu.as_mut() else { return };
        let resized = panic::catch_unwind(AssertUnwindSafe(|| gpu.resize(w, h)));
        if resized.is_err() {
            log::warn!("Surface resize panicked (likely device lost); recovering");
            self.recreate_gpu();
        }
    }

    fn new() -> Self {
        Self {
            window: None,
            gpu: None,
            camera: Camera::default(),
            mouse_down: false,
            mouse_right_down: false,
            last_mouse: PhysicalPosition::new(0.0, 0.0),
            win_size: (800.0, 800.0),
            gilrs: Gilrs::new().ok(),
            pad_prev: PadPrev::default(),
            held: HeldKeys::default(),
            auto_spin: false,
            quality_tier: 2.0,
            last_frame_at: Instant::now() - Duration::from_millis(16),
            min_frame_dt: Duration::from_millis(16),
            dynamic_quality: std::env::var("BH_DYNAMIC_QUALITY")
                .ok()
                .is_some_and(|v| matches!(v.as_str(), "1" | "true" | "TRUE" | "on" | "ON")),
            perf_window_start: Instant::now(),
            perf_accum: Duration::ZERO,
            perf_frames: 0,
            avg_frame_ms: 0.0,
            over_budget_windows: 0,
            under_budget_windows: 0,
            pad_tuning: Self::read_pad_tuning(),
            pad_name: None,
            pad_dualsense: false,
        }
    }

    fn read_pad_tuning() -> PadTuning {
        let parse_bool = |k: &str| {
            std::env::var(k).ok().is_some_and(|v| {
                matches!(
                    v.as_str(),
                    "1" | "true" | "TRUE" | "on" | "ON" | "yes" | "YES"
                )
            })
        };
        let parse_f32 = |k: &str| std::env::var(k).ok().and_then(|v| v.parse::<f32>().ok());
        let mut t = PadTuning::default();
        if let Some(v) = parse_f32("BH_PAD_DEADZONE") {
            t.deadzone = v.clamp(0.0, 0.5);
        }
        if let Some(v) = parse_f32("BH_PAD_LOOK_SENS") {
            t.look_sens = v.clamp(0.05, 8.0);
        }
        if let Some(v) = parse_f32("BH_PAD_MOVE_SENS") {
            t.move_sens = v.clamp(0.05, 8.0);
        }
        if let Some(v) = parse_f32("BH_PAD_TRIGGER_SENS") {
            t.trigger_sens = v.clamp(0.05, 8.0);
        }
        t.invert_look_x = parse_bool("BH_PAD_INVERT_X");
        t.invert_look_y = parse_bool("BH_PAD_INVERT_Y");
        t
    }

    #[inline]
    fn smooth_deadzone(v: f32, dead: f32) -> f32 {
        let a = v.abs();
        if a <= dead {
            0.0
        } else {
            let n = ((a - dead) / (1.0 - dead)).clamp(0.0, 1.0);
            v.signum() * n * n
        }
    }

    fn apply_quality_tier(&mut self, tier: f32) {
        self.quality_tier = tier.clamp(0.0, 2.0).round();
        self.min_frame_dt = if self.quality_tier < 0.5 {
            Duration::from_millis(33) // ~30 FPS low tier
        } else if self.quality_tier < 1.5 {
            Duration::from_millis(22) // ~45 FPS medium tier
        } else {
            Duration::from_millis(16) // ~60 FPS high tier
        };
    }

    fn button_pressed(now: bool, prev: &mut bool) -> bool {
        let fired = now && !*prev;
        *prev = now;
        fired
    }

    fn apply_gamepad(&mut self, el: &ActiveEventLoop) {
        let Some(gilrs) = self.gilrs.as_mut() else {
            return;
        };
        while gilrs.next_event().is_some() {}
        let Some((_, gp)) = gilrs.gamepads().find(|(_, g)| g.is_connected()) else {
            self.pad_name = None;
            self.pad_dualsense = false;
            return;
        };

        let name = gp.name().to_string();
        let is_ds = name.to_ascii_lowercase().contains("dualsense");
        let mut title_needs_update = false;
        if self.pad_name.as_deref() != Some(name.as_str()) {
            self.pad_name = Some(name);
            self.pad_dualsense = is_ds;
            title_needs_update = true;
        }

        let dead = self.pad_tuning.deadzone;
        let mut fine = 1.0_f32;
        if gp.is_pressed(Button::LeftThumb) {
            fine *= 0.35;
        }
        if gp.is_pressed(Button::RightThumb) {
            fine *= 2.2;
        }

        let lx = Self::smooth_deadzone(gp.value(Axis::LeftStickX), dead);
        let ly = Self::smooth_deadzone(gp.value(Axis::LeftStickY), dead);
        let mut rx = Self::smooth_deadzone(gp.value(Axis::RightStickX), dead);
        let mut ry = Self::smooth_deadzone(gp.value(Axis::RightStickY), dead);
        if self.pad_tuning.invert_look_x {
            rx = -rx;
        }
        if self.pad_tuning.invert_look_y {
            ry = -ry;
        }
        let lt = gp.value(Axis::LeftZ).max(0.0);
        let rt = gp.value(Axis::RightZ).max(0.0);

        let mut dirty = false;
        if rx != 0.0 || ry != 0.0 {
            self.camera.azimuth += rx * 0.06 * fine * self.pad_tuning.look_sens;
            self.camera.inclination = (self.camera.inclination
                + ry * 1.8 * fine * self.pad_tuning.look_sens)
                .clamp(1.0, 90.0);
            dirty = true;
        }
        if lx != 0.0 || ly != 0.0 {
            self.camera
                .move_local(-ly, lx, 0.0, 0.28 * fine * self.pad_tuning.move_sens);
            dirty = true;
        }
        if lt > 0.02 || rt > 0.02 {
            self.camera.move_local(
                0.0,
                0.0,
                rt - lt,
                0.30 * fine * self.pad_tuning.trigger_sens,
            );
            dirty = true;
        }

        if gp.is_pressed(Button::LeftTrigger) {
            self.camera.disk_outer = (self.camera.disk_outer - 0.10 * fine).clamp(3.5, 30.0);
            dirty = true;
        }
        if gp.is_pressed(Button::RightTrigger) {
            self.camera.disk_outer = (self.camera.disk_outer + 0.10 * fine).clamp(3.5, 30.0);
            dirty = true;
        }
        if gp.is_pressed(Button::DPadLeft) {
            self.camera.azimuth -= std::f32::consts::PI / 240.0 * fine;
            dirty = true;
        }
        if gp.is_pressed(Button::DPadRight) {
            self.camera.azimuth += std::f32::consts::PI / 240.0 * fine;
            dirty = true;
        }
        if gp.is_pressed(Button::DPadUp) {
            self.camera.inclination = (self.camera.inclination + 0.20 * fine).min(90.0);
            dirty = true;
        }
        if gp.is_pressed(Button::DPadDown) {
            self.camera.inclination = (self.camera.inclination - 0.20 * fine).max(1.0);
            dirty = true;
        }

        if Self::button_pressed(gp.is_pressed(Button::South), &mut self.pad_prev.south) {
            self.camera.gutoe_core = !self.camera.gutoe_core;
            dirty = true;
        }
        if Self::button_pressed(gp.is_pressed(Button::East), &mut self.pad_prev.east) {
            self.camera = Camera::default();
            self.auto_spin = false;
            dirty = true;
        }
        if Self::button_pressed(gp.is_pressed(Button::West), &mut self.pad_prev.west) {
            self.camera.disk_outer = if self.camera.disk_outer < 12.0 {
                20.0
            } else {
                10.0
            };
            dirty = true;
        }
        if Self::button_pressed(gp.is_pressed(Button::North), &mut self.pad_prev.north) {
            self.camera.fov_rs = if self.camera.fov_rs < 8.0 { 14.0 } else { 7.0 };
            dirty = true;
        }
        if Self::button_pressed(gp.is_pressed(Button::Start), &mut self.pad_prev.start) {
            self.auto_spin = !self.auto_spin;
        }
        if Self::button_pressed(gp.is_pressed(Button::Select), &mut self.pad_prev.select) {
            self.camera.azimuth = -std::f32::consts::FRAC_PI_2;
            dirty = true;
        }

        // Combos: safe and memorable.
        let lb = gp.is_pressed(Button::LeftTrigger);
        let rb = gp.is_pressed(Button::RightTrigger);
        if lb && rb && Self::button_pressed(gp.is_pressed(Button::Mode), &mut self.pad_prev.mode) {
            el.exit();
            return;
        }
        let lz_pressed = lt > 0.7;
        let rz_pressed = rt > 0.7;
        if lz_pressed
            && rz_pressed
            && Self::button_pressed(gp.is_pressed(Button::West), &mut self.pad_prev.left_thumb)
        {
            self.camera = Camera::default();
            self.camera.gutoe_core = true;
            self.auto_spin = false;
            dirty = true;
        }
        if self.auto_spin {
            self.camera.azimuth += 0.008 * fine;
            dirty = true;
        }

        if dirty {
            self.update_title();
            self.push_frame();
        } else if title_needs_update {
            self.update_title();
        }
    }

    fn update_title(&self) {
        let Some(win) = self.window.as_ref() else {
            return;
        };
        let core = if self.camera.gutoe_core {
            "GUTOE r_c"
        } else {
            "GR"
        };
        let disk_model = if self.camera.riaf_mode {
            "RIAF"
        } else {
            "Thin"
        };
        let stars_mode = if self.camera.local_stars {
            "stars3d:on"
        } else {
            "stars3d:off"
        };
        let riaf_vol = if self.camera.riaf_volume {
            "riafV:on"
        } else {
            "riafV:off"
        };
        let kerr = if self.camera.kerr_enable {
            format!("kerr:{:+.2}", self.camera.kerr_astar)
        } else {
            "kerr:off".to_string()
        };
        let cam_mode = if self.camera.interior_mode {
            if self.camera.core_look_mode {
                "inside→core"
            } else {
                "inside→out"
            }
        } else {
            "outside"
        };
        let spin = if self.auto_spin { " spin" } else { "" };
        let az_deg = self.camera.azimuth.to_degrees().rem_euclid(360.0);
        let roll_deg = self.camera.cam_roll.to_degrees();
        let pad = match (self.pad_name.as_deref(), self.pad_dualsense) {
            (Some(_), true) => "DualSense",
            (Some(_), false) => "Gamepad",
            (None, _) => "None",
        };
        win.set_title(&format!(
            "GUTOE BH  |  q{} {:.1}ms  {} r={:.2}r_h  inc {:.0}° az {:.0}° roll {:+.0}°  fov {:.1} r_s  disk {:.0} r_s ({})  {} {} {}  cam({:+.2},{:+.2},{:+.2})  pad:{}  [3D {}{}]",
            self.quality_tier as i32,
            self.avg_frame_ms,
            cam_mode,
            self.camera.r_cam_frac,
            self.camera.inclination,
            az_deg,
            roll_deg,
            self.camera.fov_rs,
            self.camera.disk_outer,
            disk_model,
            stars_mode,
            riaf_vol,
            kerr,
            self.camera.cam_x,
            self.camera.cam_y,
            self.camera.cam_z,
            pad,
            core,
            spin,
        ));
    }

    fn push_frame(&mut self) {
        // Frame pacing: prevents redraw storms from input/event bursts on Metal.
        let now = Instant::now();
        if now.duration_since(self.last_frame_at) < self.min_frame_dt {
            return;
        }
        self.last_frame_at = now;

        let Some(win) = self.window.as_ref() else {
            return;
        };
        let frame_start = Instant::now();
        let sz = win.inner_size();
        let pixels = (sz.width as u64).saturating_mul(sz.height as u64);
        // Fullscreen often jumps to HiDPI backbuffers (effectively ~4x pixels).
        // Cap effective quality by resolution to avoid nonlinear frametime spikes.
        let mut effective_quality = self.quality_tier;
        if pixels >= 6_000_000 {
            effective_quality = (effective_quality - 2.0).max(0.0);
        } else if pixels >= 2_600_000 {
            effective_quality = (effective_quality - 1.0).max(0.0);
        }
        let params = self
            .camera
            .params(sz.width as f32, sz.height as f32, effective_quality);
        let render_result = {
            let Some(gpu) = self.gpu.as_mut() else { return };
            gpu.upload_params(&params);
            gpu.render()
        };
        match render_result {
            Ok(_) => {}
            Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                self.safe_resize(sz.width, sz.height);
            }
            Err(wgpu::SurfaceError::OutOfMemory) => {
                log::error!("Surface out of memory; recreating GPU at lower quality");
                if self.quality_tier > 0.5 {
                    self.apply_quality_tier(self.quality_tier - 1.0);
                }
                self.recreate_gpu();
            }
            Err(e) => log::error!("render error: {e}"),
        }

        // Rolling perf window + adaptive quality downshift/upshift.
        let frame_dt = frame_start.elapsed();
        self.perf_accum += frame_dt;
        self.perf_frames += 1;
        let window_elapsed = self.perf_window_start.elapsed();
        if self.perf_frames >= 45 || window_elapsed >= Duration::from_secs(2) {
            self.avg_frame_ms =
                (self.perf_accum.as_secs_f64() * 1000.0 / (self.perf_frames.max(1) as f64)) as f32;
            if self.dynamic_quality {
                let budget_ms = self.min_frame_dt.as_secs_f32() * 1000.0;
                if self.avg_frame_ms > budget_ms * 1.35 {
                    self.over_budget_windows += 1;
                    self.under_budget_windows = 0;
                } else if self.avg_frame_ms < budget_ms * 0.70 {
                    self.under_budget_windows += 1;
                    self.over_budget_windows = 0;
                } else {
                    self.over_budget_windows = 0;
                    self.under_budget_windows = 0;
                }
                if self.over_budget_windows >= 2 && self.quality_tier > 0.5 {
                    self.apply_quality_tier(self.quality_tier - 1.0);
                    self.over_budget_windows = 0;
                    log::warn!(
                        "Adaptive quality downshift -> q{} (avg {:.1} ms)",
                        self.quality_tier as i32,
                        self.avg_frame_ms
                    );
                } else if self.under_budget_windows >= 6 && self.quality_tier < 1.5 {
                    self.apply_quality_tier(self.quality_tier + 1.0);
                    self.under_budget_windows = 0;
                    log::info!(
                        "Adaptive quality upshift -> q{} (avg {:.1} ms)",
                        self.quality_tier as i32,
                        self.avg_frame_ms
                    );
                }
            }
            self.perf_window_start = Instant::now();
            self.perf_accum = Duration::ZERO;
            self.perf_frames = 0;
            self.update_title();
        }
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, el: &ActiveEventLoop) {
        let win = Arc::new(
            el.create_window(
                winit::window::WindowAttributes::default()
                    .with_title("GUTOE Black Hole — GRAND-151")
                    .with_inner_size(winit::dpi::LogicalSize::new(900u32, 900u32)),
            )
            .unwrap(),
        );
        let sz = win.inner_size();
        self.win_size = (sz.width.max(1) as f32, sz.height.max(1) as f32);
        let gpu = Gpu::new(Arc::clone(&win));
        self.apply_quality_tier(gpu.quality_tier);
        self.window = Some(win);
        self.gpu = Some(gpu);
        self.update_title();
        self.push_frame();
    }

    fn window_event(&mut self, el: &ActiveEventLoop, _id: WindowId, event: WindowEvent) {
        match event {
            // ── Quit ──────────────────────────────────────────────────────────
            WindowEvent::CloseRequested => el.exit(),

            WindowEvent::KeyboardInput {
                event: KeyEvent {
                    logical_key, state, ..
                },
                ..
            } => {
                let pressed = state == ElementState::Pressed;
                match logical_key {
                    Key::Named(NamedKey::Escape) if pressed => el.exit(),
                    Key::Named(NamedKey::ArrowUp) if pressed => {
                        self.camera.inclination = (self.camera.inclination + 5.0).min(90.0);
                        self.update_title();
                        self.push_frame();
                    }
                    Key::Named(NamedKey::ArrowDown) if pressed => {
                        self.camera.inclination = (self.camera.inclination - 5.0).max(1.0);
                        self.update_title();
                        self.push_frame();
                    }
                    Key::Named(NamedKey::ArrowLeft) if pressed => {
                        self.camera.azimuth -= std::f32::consts::PI / 12.0;
                        self.update_title();
                        self.push_frame();
                    }
                    Key::Named(NamedKey::ArrowRight) if pressed => {
                        self.camera.azimuth += std::f32::consts::PI / 12.0;
                        self.update_title();
                        self.push_frame();
                    }
                    Key::Character(ref s) => match s.as_str() {
                        "q" | "Q" if pressed => el.exit(),
                        "r" | "R" if pressed => {
                            self.camera = Camera::default();
                            self.update_title();
                            self.push_frame();
                        }
                        "w" | "W" => self.held.fwd = pressed,
                        "s" | "S" => self.held.back = pressed,
                        "a" | "A" => self.held.left = pressed,
                        "d" | "D" => self.held.right = pressed,
                        "z" | "Z" => self.held.zoom_in = pressed,
                        "x" | "X" => self.held.zoom_out = pressed,
                        "e" | "E" => self.held.roll_pos = pressed,
                        "c" | "C" => self.held.roll_neg = pressed,
                        "i" | "I" if pressed => {
                            self.camera.interior_mode = !self.camera.interior_mode;
                            if !self.camera.interior_mode {
                                self.camera.core_look_mode = false;
                            }
                            self.update_title();
                            self.push_frame();
                        }
                        "o" | "O" if pressed => {
                            self.camera.core_look_mode = !self.camera.core_look_mode;
                            if self.camera.core_look_mode {
                                self.camera.interior_mode = true;
                            }
                            self.update_title();
                            self.push_frame();
                        }
                        "[" if pressed => {
                            self.camera.r_cam_frac =
                                (self.camera.r_cam_frac - 0.02).clamp(0.05, 0.99);
                            self.update_title();
                            self.push_frame();
                        }
                        "]" if pressed => {
                            self.camera.r_cam_frac =
                                (self.camera.r_cam_frac + 0.02).clamp(0.05, 0.99);
                            self.update_title();
                            self.push_frame();
                        }
                        "g" | "G" if pressed => {
                            self.camera.gutoe_core = !self.camera.gutoe_core;
                            log::info!("GUTOE lattice core: {}", self.camera.gutoe_core);
                            self.update_title();
                            self.push_frame();
                        }
                        "t" | "T" if pressed => {
                            self.camera.riaf_mode = !self.camera.riaf_mode;
                            log::info!(
                                "Disk model: {}",
                                if self.camera.riaf_mode {
                                    "RIAF composite"
                                } else {
                                    "Thin disk"
                                }
                            );
                            self.update_title();
                            self.push_frame();
                        }
                        "v" | "V" if pressed => {
                            self.camera.local_stars = !self.camera.local_stars;
                            log::info!("Local 3D stars: {}", self.camera.local_stars);
                            self.update_title();
                            self.push_frame();
                        }
                        "m" | "M" if pressed => {
                            self.camera.riaf_volume = !self.camera.riaf_volume;
                            log::info!("Volumetric RIAF blend: {}", self.camera.riaf_volume);
                            self.update_title();
                            self.push_frame();
                        }
                        "k" | "K" if pressed => {
                            self.camera.kerr_enable = !self.camera.kerr_enable;
                            log::info!(
                                "Kerr mode: {} (a*={:+.3})",
                                self.camera.kerr_enable,
                                self.camera.kerr_astar
                            );
                            self.update_title();
                            self.push_frame();
                        }
                        "," if pressed => {
                            self.camera.kerr_astar =
                                (self.camera.kerr_astar - 0.05).clamp(-0.999, 0.999);
                            self.camera.kerr_enable = self.camera.kerr_astar.abs() > 1e-6;
                            self.update_title();
                            self.push_frame();
                        }
                        "." if pressed => {
                            self.camera.kerr_astar =
                                (self.camera.kerr_astar + 0.05).clamp(-0.999, 0.999);
                            self.camera.kerr_enable = self.camera.kerr_astar.abs() > 1e-6;
                            self.update_title();
                            self.push_frame();
                        }
                        "=" | "+" if pressed => {
                            self.camera.disk_outer = (self.camera.disk_outer + 1.0).min(30.0);
                            self.update_title();
                            self.push_frame();
                        }
                        "-" if pressed => {
                            self.camera.disk_outer = (self.camera.disk_outer - 1.0).max(3.5);
                            self.update_title();
                            self.push_frame();
                        }
                        _ => {}
                    },
                    _ => {}
                }
            }

            // ── Resize ────────────────────────────────────────────────────────
            WindowEvent::Resized(size) => {
                self.win_size = (size.width.max(1) as f32, size.height.max(1) as f32);
                self.safe_resize(size.width, size.height);
                self.push_frame();
            }

            // ── Mouse ─────────────────────────────────────────────────────────
            WindowEvent::MouseInput { button, state, .. } => match button {
                MouseButton::Left => self.mouse_down = state == ElementState::Pressed,
                MouseButton::Right => self.mouse_right_down = state == ElementState::Pressed,
                _ => {}
            },

            WindowEvent::CursorMoved { position, .. } => {
                if self.mouse_down || self.mouse_right_down {
                    let dx = (position.x - self.last_mouse.x) as f32;
                    let dy = (position.y - self.last_mouse.y) as f32;
                    let (w, h) = self.win_size;

                    if self.mouse_down {
                        // Left drag = look.
                        self.camera.inclination =
                            (self.camera.inclination - dy * 90.0 / h).clamp(1.0, 90.0);
                        self.camera.azimuth -= dx * std::f32::consts::TAU / w;
                    }
                    if self.mouse_right_down {
                        // Right drag = translate in camera plane (HL-style noclip feel).
                        let strafe = -dx / w * self.camera.fov_rs * 1.8;
                        let lift = dy / h * self.camera.fov_rs * 1.8;
                        self.camera.move_local(0.0, strafe, lift, 1.0);
                    }

                    self.update_title();
                    self.push_frame();
                }
                self.last_mouse = position;
            }

            // ── Scroll / zoom ─────────────────────────────────────────────────
            WindowEvent::MouseWheel { delta, .. } => {
                // Positive scroll = zoom in (fov shrinks), negative = zoom out.
                // Line delta: 1 notch ≈ 12% fov change — snappy but not jumpy.
                // Pixel delta (trackpad): scale so ~8px movement = same as one notch.
                let scroll = match delta {
                    MouseScrollDelta::LineDelta(_, y) => y,
                    MouseScrollDelta::PixelDelta(p) => p.y as f32 / 8.0,
                };
                self.camera.fov_rs = (self.camera.fov_rs * (1.0 - scroll * 0.12)).clamp(3.0, 40.0);
                self.update_title();
                self.push_frame();
            }

            WindowEvent::RedrawRequested => self.push_frame(),
            _ => {}
        }
    }

    fn about_to_wait(&mut self, _el: &ActiveEventLoop) {
        self.apply_gamepad(_el);
        // Continuous keyboard motion (WASD-style noclip).
        let mut moved = false;
        let speed = 0.12;
        if self.held.fwd
            || self.held.back
            || self.held.left
            || self.held.right
            || self.held.up
            || self.held.down
        {
            let forward = (self.held.fwd as i32 - self.held.back as i32) as f32;
            let right = (self.held.right as i32 - self.held.left as i32) as f32;
            let up = (self.held.up as i32 - self.held.down as i32) as f32;
            self.camera.move_local(forward, right, up, speed);
            moved = true;
        }
        if self.held.roll_pos {
            self.camera.cam_roll += std::f32::consts::PI / 180.0 * 1.8;
            moved = true;
        }
        if self.held.roll_neg {
            self.camera.cam_roll -= std::f32::consts::PI / 180.0 * 1.8;
            moved = true;
        }
        if self.held.zoom_in {
            self.camera.move_local(1.0, 0.0, 0.0, speed);
            moved = true;
        }
        if self.held.zoom_out {
            self.camera.move_local(-1.0, 0.0, 0.0, speed);
            moved = true;
        }
        if moved {
            self.update_title();
        }
        if let Some(win) = self.window.as_ref() {
            // Always redraw for responsive controls; frame pacing still limits GPU load.
            win.request_redraw();
        }
    }
}

// ── Entry point ───────────────────────────────────────────────────────────────

fn main() {
    env_logger::init();

    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  GUTOE Black Hole Viewer — GRAND-151");
    println!("  GUTOE Schwarzschild metric + SC lattice (Cl(1,3)) + EHT-like rendering");
    println!("  Physics: Novikov-Thorne disk · Doppler beaming · gravitational redshift");
    println!("           Lensed star field · Hawking temperature · singularity-free core");
    println!("  Default: M87-like — 17° inclination, bright crescent at bottom");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  Left-drag          — look (yaw/pitch), HL2-style freecam");
    println!("  Right-drag         — translate camera in view plane");
    println!("  Scroll             — zoom (field of view in units of r_s)");
    println!("  WASD               — freecam noclip strafe/forward in camera-local frame");
    println!("  Z / X              — freecam dolly in/out");
    println!("  E / C              — roll camera clockwise/counterclockwise");
    println!("  I                  — toggle outside/inside horizon camera");
    println!("  O                  — toggle inside-core look mode");
    println!("  [ / ]              — interior camera radius (r_cam/r_h)");
    println!("  + / -              — disk outer radius (grow / shrink accretion disk)");
    println!("  T                  — toggle disk model (RIAF ↔ thin)");
    println!("  V                  — toggle local 3D star parallax shells");
    println!("  M                  — toggle escaped-ray volumetric RIAF blend");
    println!("  K                  — toggle Kerr mode");
    println!("  , / .              — Kerr spin a* down / up");
    println!("  G                  — toggle GUTOE lattice core r_c  (GR ↔ GUTOE)");
    println!("  R                  — reset to M87-like defaults");
    println!("  Q / Escape         — quit");
    println!("  Gamepad (Xbox / DualSense):");
    println!("    LS move (strafe/forward) · RS look (yaw/pitch) · LT/RT up/down");
    println!("    L1/R1 (or LB/RB) disk size");
    println!("    D-pad nudge camera · A toggle GUTOE/GR · B reset · X disk preset · Y fov preset");
    println!("    Start auto-spin · Back recenter azimuth · LB+RB+Guide quit");
    println!("    LT+RT+X hard reset (default M87 + GUTOE core)");
    println!(
        "    Env tuning: BH_PAD_DEADZONE, BH_PAD_LOOK_SENS, BH_PAD_MOVE_SENS, BH_PAD_TRIGGER_SENS"
    );
    println!("    Env invert: BH_PAD_INVERT_X, BH_PAD_INVERT_Y");
    println!("  Title bar shows live: inclination, azimuth, fov, disk size, mode");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let ev = match EventLoop::new() {
        Ok(ev) => ev,
        Err(e) => {
            eprintln!("failed to create event loop: {e}");
            return;
        }
    };
    if let Err(e) = ev.run_app(&mut App::new()) {
        eprintln!("viewer runtime error: {e}");
    }
}
