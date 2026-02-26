//! GUTOE Fusion Viewer — real-time multi-zone stellar burn visualization.

use std::{
    panic::{self, AssertUnwindSafe},
    sync::Arc,
    time::{Duration, Instant},
};

use gilrs::{Axis, Button, Gilrs};
use gutoe_physics::{MultiZoneBurn, Species, ZoneState};
use winit::{
    application::ApplicationHandler,
    event::{ElementState, KeyEvent, WindowEvent},
    event_loop::{ActiveEventLoop, EventLoop},
    keyboard::{Key, NamedKey},
    window::{Window, WindowId},
};

#[derive(Default, Clone, Copy)]
struct PadPrev {
    right_thumb: bool,
    west: bool,
    start: bool,
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
    yaw_left: bool,
    yaw_right: bool,
    pitch_up: bool,
    pitch_down: bool,
    roll_left: bool,
    roll_right: bool,
    turbo: bool,
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    width: f32,
    height: f32,
    time: f32,
    fusion_power: f32,
    h1: f32,
    he4: f32,
    running: f32,
    speed_norm: f32,
    zone_p0: f32,
    zone_p1: f32,
    zone_p2: f32,
    zone_p3: f32,
    zone_h10: f32,
    zone_h11: f32,
    zone_h12: f32,
    zone_h13: f32,
    zone_he40: f32,
    zone_he41: f32,
    zone_he42: f32,
    zone_he43: f32,
    zone_count: f32,
    diffusion_log10: f32,
    t9_scale: f32,
    _zone_pad: f32,
    cam_x: f32,
    cam_y: f32,
    cam_z: f32,
    _cam_pad0: f32,
    cam_rx: f32,
    cam_ry: f32,
    cam_rz: f32,
    _cam_pad1: f32,
    cam_ux: f32,
    cam_uy: f32,
    cam_uz: f32,
    _cam_pad2: f32,
    cam_fx: f32,
    cam_fy: f32,
    cam_fz: f32,
    fov_tan: f32,
    quality_tier: f32,
    draw_distance: f32,
    planet_scale: f32,
    planet_time_scale: f32,
    _render_pad0: f32,
    _render_pad1: f32,
    _render_pad2: f32,
}

const STAR_RADIUS_UNITS: f32 = 5.4;

const SHADER: &str = r#"
struct Params {
    width       : f32,
    height      : f32,
    time        : f32,
    fusion_power: f32,
    h1          : f32,
    he4         : f32,
    running     : f32,
    speed_norm  : f32,
    zone_p0     : f32,
    zone_p1     : f32,
    zone_p2     : f32,
    zone_p3     : f32,
    zone_h10    : f32,
    zone_h11    : f32,
    zone_h12    : f32,
    zone_h13    : f32,
    zone_he40   : f32,
    zone_he41   : f32,
    zone_he42   : f32,
    zone_he43   : f32,
    zone_count  : f32,
    diffusion_log10: f32,
    t9_scale    : f32,
    _zone_pad   : f32,
    cam_x       : f32,
    cam_y       : f32,
    cam_z       : f32,
    _cam_pad0   : f32,
    cam_rx      : f32,
    cam_ry      : f32,
    cam_rz      : f32,
    _cam_pad1   : f32,
    cam_ux      : f32,
    cam_uy      : f32,
    cam_uz      : f32,
    _cam_pad2   : f32,
    cam_fx      : f32,
    cam_fy      : f32,
    cam_fz      : f32,
    fov_tan     : f32,
    quality_tier: f32,
    draw_distance: f32,
    planet_scale: f32,
    planet_time_scale: f32,
    _render_pad0: f32,
    _render_pad1: f32,
    _render_pad2: f32,
}
@group(0) @binding(0) var<uniform> P : Params;

fn hash31(p: vec3<f32>) -> f32 {
    var q = fract(p * vec3(127.1, 311.7, 269.5));
    q += dot(q, q.yzx + 19.19);
    return fract((q.x + q.y) * q.z);
}

fn noise3(p: vec3<f32>) -> f32 {
    let i = floor(p);
    let f = fract(p);
    let u = f * f * (3.0 - 2.0 * f);
    let a000 = hash31(i + vec3(0.0, 0.0, 0.0));
    let a100 = hash31(i + vec3(1.0, 0.0, 0.0));
    let a010 = hash31(i + vec3(0.0, 1.0, 0.0));
    let a110 = hash31(i + vec3(1.0, 1.0, 0.0));
    let a001 = hash31(i + vec3(0.0, 0.0, 1.0));
    let a101 = hash31(i + vec3(1.0, 0.0, 1.0));
    let a011 = hash31(i + vec3(0.0, 1.0, 1.0));
    let a111 = hash31(i + vec3(1.0, 1.0, 1.0));
    let x00 = mix(a000, a100, u.x);
    let x10 = mix(a010, a110, u.x);
    let x01 = mix(a001, a101, u.x);
    let x11 = mix(a011, a111, u.x);
    let y0 = mix(x00, x10, u.y);
    let y1 = mix(x01, x11, u.y);
    return mix(y0, y1, u.z);
}

fn fbm3(p0: vec3<f32>) -> f32 {
    var p = p0;
    var acc = 0.0;
    var amp = 0.5;
    for (var i = 0; i < 5; i++) {
        acc += amp * noise3(p);
        p = p * 2.03 + vec3(17.7, 9.2, 13.4);
        amp *= 0.5;
    }
    return acc;
}

fn starfield(dir: vec3<f32>) -> vec3<f32> {
    let d = normalize(dir);
    let lon = atan2(d.x, d.z);
    let lat = asin(clamp(d.y, -1.0, 1.0));
    let uv = vec2(lon / (2.0 * 3.14159265359) + 0.5, 0.5 - lat / 3.14159265359);
    let cell = floor(uv * vec2(1600.0, 800.0));
    let gate = hash31(vec3(cell, floor(d.z * 913.0)));
    var c = vec3(0.003, 0.004, 0.010);
    let star_w = smoothstep(0.9976, 1.0, gate);
    if (star_w > 1e-5) {
        let local = fract(uv * vec2(1600.0, 800.0)) - 0.5;
        let r2 = dot(local, local);
        let psf = exp(-r2 * 70.0);
        let t = hash31(vec3(cell + vec2(5.0, 11.0), 2.0));
        let scol = mix(vec3(0.75, 0.86, 1.00), vec3(1.00, 0.82, 0.60), t);
        c += scol * psf * star_w * (0.2 + 0.8 * pow(hash31(vec3(cell, 7.0)), 2.0));
    }
    return c;
}

// Position-aware local stellar shell volume. Translation of camera produces
// coherent parallax, unlike pure direction-space background sampling.
fn local_star_volume(cam: vec3<f32>, dir: vec3<f32>, depth: f32, scale: f32, quality: f32) -> vec3<f32> {
    let steps = u32(mix(8.0, 16.0, clamp(quality / 2.0, 0.0, 1.0)));
    var c = vec3(0.0);
    let inv = 1.0 / max(f32(steps), 1.0);
    for (var i = 0u; i < 24u; i++) {
        if (i >= steps) {
            break;
        }
        let t = (f32(i) + 0.5) * depth * inv;
        let p = cam + dir * t;
        let cell = floor(p * scale);
        let gate = hash31(cell + vec3(17.0, 31.0, 47.0));
        let star_w = smoothstep(0.9975, 1.0, gate);
        if (star_w > 1e-6) {
            let q = fract(p * scale) - 0.5;
            let r2 = dot(q, q);
            let psf = exp(-r2 * 55.0);
            let tint = hash31(cell + vec3(3.0, 11.0, 19.0));
            let scol = mix(vec3(0.72, 0.85, 1.00), vec3(1.00, 0.82, 0.60), tint);
            c += scol * psf * star_w * (0.8 + 0.2 * tint) * inv * 2.5;
        }
    }
    return c;
}

fn ray_sphere(ro: vec3<f32>, rd: vec3<f32>, center: vec3<f32>, radius: f32) -> f32 {
    let oc = ro - center;
    let b = dot(oc, rd);
    let c = dot(oc, oc) - radius * radius;
    let h = b * b - c;
    if (h < 0.0) {
        return -1.0;
    }
    let s = sqrt(h);
    let t0 = -b - s;
    if (t0 > 0.0) {
        return t0;
    }
    let t1 = -b + s;
    if (t1 > 0.0) {
        return t1;
    }
    return -1.0;
}

struct VOut { @builtin(position) pos: vec4<f32> }

@vertex
fn vs(@builtin(vertex_index) vi: u32) -> VOut {
    var xy = array<vec2<f32>, 6>(
        vec2(-1.0, -1.0), vec2(1.0, -1.0), vec2(-1.0, 1.0),
        vec2(-1.0,  1.0), vec2(1.0, -1.0), vec2(1.0,  1.0),
    );
    return VOut(vec4(xy[vi], 0.0, 1.0));
}

@fragment
fn fs(@builtin(position) frag: vec4<f32>) -> @location(0) vec4<f32> {
    let res = vec2(max(P.width, 1.0), max(P.height, 1.0));
    let burn = clamp(P.fusion_power, 0.0, 3.0);
    let quality = clamp(P.quality_tier, 0.0, 2.0);
    let comp = clamp(P.he4 / max(P.h1 + P.he4, 1e-6), 0.0, 1.0);
    let ndc = vec2(
        (frag.x / max(P.width, 1.0)) * 2.0 - 1.0,
        1.0 - (frag.y / max(P.height, 1.0)) * 2.0
    );
    let aspect = P.width / max(P.height, 1.0);
    let ray_cam = normalize(vec3(ndc.x * aspect * P.fov_tan, ndc.y * P.fov_tan, -1.0));

    let cam_pos = vec3(P.cam_x, P.cam_y, P.cam_z);
    let cam_right = normalize(vec3(P.cam_rx, P.cam_ry, P.cam_rz));
    let cam_up = normalize(vec3(P.cam_ux, P.cam_uy, P.cam_uz));
    let cam_fwd = normalize(vec3(P.cam_fx, P.cam_fy, P.cam_fz));
    let dir = normalize(cam_right * ray_cam.x + cam_up * ray_cam.y + cam_fwd * (-ray_cam.z));

    let draw_dist = max(P.draw_distance, 6.0);
    let planet_time = P.time * max(P.planet_time_scale, 0.0);
    var planet_t = draw_dist + 1.0;
    var planet_hit = -1;
    var planet_center = vec3(0.0, 0.0, 0.0);
    var planet_radius = 0.0;
    let planets_on = P.planet_scale > 1e-4;
    if (planets_on) {
        let ps = max(P.planet_scale, 0.05);
        let p0 = vec3(
            cos(planet_time * 0.55) * 18.0,
            sin(planet_time * 0.21) * 1.3,
            sin(planet_time * 0.55) * 18.0
        );
        let p1 = vec3(
            cos(planet_time * 0.31 + 1.7) * 31.0,
            sin(planet_time * 0.13 + 0.9) * 2.0,
            sin(planet_time * 0.31 + 1.7) * 31.0
        );
        let p2 = vec3(
            cos(planet_time * 0.17 + 4.2) * 52.0,
            sin(planet_time * 0.09 + 2.4) * 3.1,
            sin(planet_time * 0.17 + 4.2) * 52.0
        );

        let t0 = ray_sphere(cam_pos, dir, p0, 0.58 * ps);
        if (t0 > 0.0 && t0 < planet_t && t0 < draw_dist) {
            planet_t = t0;
            planet_hit = 0;
            planet_center = p0;
            planet_radius = 0.58 * ps;
        }
        let t1 = ray_sphere(cam_pos, dir, p1, 0.84 * ps);
        if (t1 > 0.0 && t1 < planet_t && t1 < draw_dist) {
            planet_t = t1;
            planet_hit = 1;
            planet_center = p1;
            planet_radius = 0.84 * ps;
        }
        let t2 = ray_sphere(cam_pos, dir, p2, 1.25 * ps);
        if (t2 > 0.0 && t2 < planet_t && t2 < draw_dist) {
            planet_t = t2;
            planet_hit = 2;
            planet_center = p2;
            planet_radius = 1.25 * ps;
        }
    }

    var col = starfield(dir);
    let local_depth = min(draw_dist, mix(35.0, 95.0, quality / 2.0));
    col += local_star_volume(cam_pos, dir, local_depth, 0.22, quality) * 0.16;
    if (quality > 0.9) {
        col += local_star_volume(cam_pos, dir, local_depth * 0.7, 0.42, quality) * 0.09;
    }
    if (quality > 1.6) {
        col += local_star_volume(cam_pos, dir, local_depth * 0.55, 0.85, quality) * 0.05;
    }
    let hot = mix(vec3(0.72, 0.87, 1.00), vec3(1.00, 0.73, 0.32), comp);
    let zpow = vec4(
        clamp(P.zone_p0, 0.0, 3.0),
        clamp(P.zone_p1, 0.0, 3.0),
        clamp(P.zone_p2, 0.0, 3.0),
        clamp(P.zone_p3, 0.0, 3.0)
    );
    let zh1 = vec4(
        clamp(P.zone_h10, 0.0, 1.0),
        clamp(P.zone_h11, 0.0, 1.0),
        clamp(P.zone_h12, 0.0, 1.0),
        clamp(P.zone_h13, 0.0, 1.0)
    );
    let zhe = vec4(
        clamp(P.zone_he40, 0.0, 1.0),
        clamp(P.zone_he41, 0.0, 1.0),
        clamp(P.zone_he42, 0.0, 1.0),
        clamp(P.zone_he43, 0.0, 1.0)
    );
    let zc0 = clamp(zhe.x / max(zh1.x + zhe.x, 1e-6), 0.0, 1.0);
    let zc1 = clamp(zhe.y / max(zh1.y + zhe.y, 1e-6), 0.0, 1.0);
    let zc2 = clamp(zhe.z / max(zh1.z + zhe.z, 1e-6), 0.0, 1.0);
    let zc3 = clamp(zhe.w / max(zh1.w + zhe.w, 1e-6), 0.0, 1.0);
    let hot0 = mix(vec3(0.72, 0.87, 1.00), vec3(1.00, 0.73, 0.32), zc0);
    let hot1 = mix(vec3(0.72, 0.87, 1.00), vec3(1.00, 0.73, 0.32), zc1);
    let hot2 = mix(vec3(0.72, 0.87, 1.00), vec3(1.00, 0.73, 0.32), zc2);
    let hot3 = mix(vec3(0.72, 0.87, 1.00), vec3(1.00, 0.73, 0.32), zc3);
    let zone_gain = 0.75 + 0.15 * clamp(P.zone_count, 1.0, 12.0);

    let pix_jitter = hash31(vec3(floor(frag.xy), 17.0));
    var t = 0.20 + 0.14 * pix_jitter;
    var trans = 1.0;
    let march_end = min(draw_dist, select(draw_dist, planet_t, planet_hit >= 0));
    let max_steps = mix(30.0, mix(54.0, 76.0, step(1.5, quality)), step(0.5, quality));
    for (var i = 0; i < 84; i++) {
        if (f32(i) >= max_steps) {
            break;
        }
        if (t > march_end) {
            break;
        }
        let pos = cam_pos + dir * t;
        let r = length(pos);
        // Keep volumetric turbulence slow + low-amplitude to avoid aggressive flicker.
        let flow_base = fbm3(pos * (0.52 + 0.06 * burn) + vec3(11.3, 7.1, 5.7));
        let flow_drift = fbm3(
            pos * (0.22 + 0.03 * burn)
            + vec3(0.0035 * P.time, -0.0025 * P.time, 0.0020 * P.time)
        );
        let flow = mix(flow_base, flow_drift, 0.18);
        let turb = 0.82 + 0.18 * flow;
        let r0 = 1.20 + 0.05 * burn;
        let r1 = 2.35 + 0.10 * burn;
        let r2 = 3.80 + 0.18 * burn;
        let r3 = 5.40 + 0.22 * burn;
        let w0 = 0.55;
        let w1 = 0.80;
        let w2 = 1.20;
        let w3 = 1.65;
        let d0 = exp(-pow((r - r0) / w0, 2.0)) * zpow.x;
        let d1 = exp(-pow((r - r1) / w1, 2.0)) * zpow.y;
        let d2 = exp(-pow((r - r2) / w2, 2.0)) * zpow.z;
        let d3 = exp(-pow((r - r3) / w3, 2.0)) * zpow.w;
        let shell_mix = d0 + d1 + d2 + d3;
        let core_fog = exp(-r * r * (0.10 + 0.035 * burn));
        let density = max((zone_gain * shell_mix + 0.35 * core_fog) * turb - 0.24, 0.0);
        let zone_emis = hot0 * d0 + hot1 * d1 + hot2 * d2 + hot3 * d3;
        let emis = (zone_emis + hot * 0.22 * core_fog) * (0.55 + 1.15 * burn);
        let step_base = mix(0.11, 0.42, clamp(t / max(draw_dist, 1.0), 0.0, 1.0));
        let tau_seg = max(density * (0.08 + 0.035 * burn) * (0.8 + 0.25 * step_base), 0.0);
        let trans_seg = clamp(exp(-tau_seg), 0.0, 1.0);
        let alpha = 1.0 - trans_seg;
        col += trans * emis * alpha;
        trans *= trans_seg;
        if trans < 0.02 {
            break;
        }
        let q_step = mix(0.86, mix(1.0, 1.14, step(1.5, quality)), step(0.5, quality));
        t += (step_base + 0.018 * t) * q_step;
    }

    let b = length(cross(cam_pos, dir));
    let toward = select(0.0, 1.0, dot(-cam_pos, dir) > 0.0);
    let core = exp(-b * b * (22.0 + 8.0 * burn)) * toward;
    let corona = exp(-b * b * 3.0) * toward;
    col += hot * core * (0.5 + 1.9 * burn);
    col += vec3(1.0, 0.88, 0.65) * corona * (0.10 + 0.16 * burn + 0.02 * max(P.t9_scale, 0.0));

    if (planet_hit >= 0) {
        let hit_pos = cam_pos + dir * planet_t;
        let nrm = normalize(hit_pos - planet_center);
        let ldir = normalize(-hit_pos);
        let ndotl = clamp(dot(nrm, ldir), 0.0, 1.0);
        let vdot = clamp(dot(nrm, normalize(cam_pos - hit_pos)), 0.0, 1.0);
        let rim = pow(1.0 - vdot, 2.2);

        let uvw = nrm * (8.0 + 3.0 * f32(planet_hit)) + vec3(0.0, 0.0, 0.05 * planet_time);
        let cloud = fbm3(uvw);
        var base_col = vec3(0.40, 0.42, 0.46);
        var atm_col = vec3(0.50, 0.62, 0.92);
        if (planet_hit == 0) {
            base_col = mix(vec3(0.52, 0.38, 0.26), vec3(0.66, 0.47, 0.28), cloud);
            atm_col = vec3(0.95, 0.68, 0.32);
        } else if (planet_hit == 1) {
            base_col = mix(vec3(0.24, 0.39, 0.60), vec3(0.17, 0.25, 0.16), cloud);
            atm_col = vec3(0.52, 0.70, 1.00);
        } else {
            base_col = mix(vec3(0.74, 0.62, 0.42), vec3(0.90, 0.78, 0.55), cloud);
            atm_col = vec3(0.98, 0.88, 0.66);
        }

        let city = smoothstep(0.20, 0.85, cloud) * smoothstep(-0.2, 0.3, dot(nrm, ldir));
        let night = 0.06 + 0.18 * city;
        let diffuse = night + 0.86 * ndotl;
        let spec = pow(max(dot(reflect(-ldir, nrm), normalize(cam_pos - hit_pos)), 0.0), 24.0);
        let planet_col = base_col * diffuse + 0.08 * spec + atm_col * rim * (0.35 + 0.15 * ndotl);
        col = planet_col;

        let dist_to_star = length(hit_pos);
        let transit = smoothstep(planet_radius * 1.2, planet_radius * 0.35, dist_to_star);
        col *= 1.0 - 0.28 * transit;
    }

    // Pause mode: subtly desaturate/dim so state is obvious.
    if (P.running < 0.5) {
        let gray = dot(col, vec3(0.299, 0.587, 0.114));
        col = mix(col, vec3(gray), 0.30);
        col *= 0.82;
    }

    col = col / (vec3(1.0) + col);
    col = pow(col, vec3(1.0 / 2.2));
    return vec4(col, 1.0);
}
"#;

#[derive(Clone, Copy)]
struct FlightCamera {
    pos: [f32; 3],
    yaw: f32,
    pitch: f32,
    roll: f32,
    fov_deg: f32,
}

impl Default for FlightCamera {
    fn default() -> Self {
        Self {
            pos: [0.0, 0.0, 8.0],
            yaw: std::f32::consts::PI,
            pitch: 0.0,
            roll: 0.0,
            fov_deg: 72.0,
        }
    }
}

fn v_add(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn v_scale(v: [f32; 3], s: f32) -> [f32; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn v_cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn v_norm(v: [f32; 3]) -> [f32; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-6);
    [v[0] / n, v[1] / n, v[2] / n]
}

impl FlightCamera {
    fn basis(&self) -> ([f32; 3], [f32; 3], [f32; 3]) {
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        let forward = v_norm([cp * sy, sp, cp * cy]);
        let world_up = if forward[1].abs() > 0.97 {
            [0.0, 0.0, 1.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let right = v_norm(v_cross(forward, world_up));
        let up0 = v_norm(v_cross(right, forward));
        let (sr, cr) = self.roll.sin_cos();
        let right_r = v_norm(v_add(v_scale(right, cr), v_scale(up0, sr)));
        let up_r = v_norm(v_add(v_scale(up0, cr), v_scale(right, -sr)));
        (right_r, up_r, forward)
    }

    fn move_local(&mut self, forward: f32, right: f32, up: f32, speed: f32, dt: f32) {
        let (rv, uv, fv) = self.basis();
        let v = v_add(
            v_add(v_scale(fv, forward), v_scale(rv, right)),
            v_scale(uv, up),
        );
        self.pos = v_add(self.pos, v_scale(v, speed * dt));
    }

    fn look(&mut self, yaw_delta: f32, pitch_delta: f32) {
        self.yaw += yaw_delta;
        self.pitch = (self.pitch + pitch_delta).clamp(
            -std::f32::consts::FRAC_PI_2 + 0.01,
            std::f32::consts::FRAC_PI_2 - 0.01,
        );
    }
}

struct FusionRuntime {
    burn: MultiZoneBurn,
    zones: Vec<ZoneState>,
    zone_count: usize,
    running: bool,
    speed: f32,
    t9_scale: f64,
    diffusion_coeff: f64,
    wall_to_model_scale: f64,
    sim_time_s: f64,
    prev_h1_mean: Option<f64>,
    prev_he4_mean: Option<f64>,
    diag: FusionDiagnostics,
}

#[derive(Clone, Copy, Default)]
struct FusionDiagnostics {
    sum_x_mean: f64,
    sum_x_max_dev: f64,
    min_x: f64,
    dh_dt: f64,
    dhe4_dt: f64,
    pp_power: f64,
    cno_power: f64,
    triple_alpha_power: f64,
}

impl Default for FusionRuntime {
    fn default() -> Self {
        let parse_bool = |k: &str, dflt: bool| {
            std::env::var(k).ok().map_or(dflt, |s| {
                matches!(
                    s.as_str(),
                    "1" | "true" | "TRUE" | "on" | "ON" | "yes" | "YES"
                )
            })
        };
        let parse_f32 = |k: &str, dflt: f32| {
            std::env::var(k)
                .ok()
                .and_then(|s| s.parse::<f32>().ok())
                .unwrap_or(dflt)
        };
        let parse_f64 = |k: &str, dflt: f64| {
            std::env::var(k)
                .ok()
                .and_then(|s| s.parse::<f64>().ok())
                .unwrap_or(dflt)
        };
        let zone_count = std::env::var("FUSION_VIEWER_ZONES")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(64)
            .clamp(8, 256);
        let t9_scale = (parse_f64("FUSION_VIEWER_T9", 0.02) / 0.02).clamp(0.05, 8.0);
        let diffusion_coeff = parse_f64("FUSION_VIEWER_DIFFUSION", 1.0e-10).max(0.0);
        let mut burn = MultiZoneBurn::baseline();
        burn.cfg.diffusion_coeff = diffusion_coeff;
        burn.cfg.zone_temperatures_t9 = Self::build_temperature_profile(zone_count, t9_scale);
        burn.cfg.zone_rate_scales = Self::build_rate_profile(zone_count);
        let zones = burn.seed_zones(zone_count);
        let mut s = Self {
            burn,
            zones,
            zone_count,
            running: parse_bool("FUSION_VIEWER_RUN", true),
            speed: parse_f32("FUSION_VIEWER_SPEED", 1.0).clamp(0.05, 100.0),
            t9_scale,
            diffusion_coeff,
            wall_to_model_scale: parse_f64("FUSION_VIEWER_TIME_SCALE", 5.0e7).max(1.0),
            sim_time_s: 0.0,
            prev_h1_mean: None,
            prev_he4_mean: None,
            diag: FusionDiagnostics::default(),
        };
        s.refresh_diagnostics(0.0);
        s
    }
}

impl FusionRuntime {
    fn build_temperature_profile(zone_count: usize, t9_scale: f64) -> Vec<f64> {
        let mut temps = Vec::with_capacity(zone_count);
        for i in 0..zone_count {
            let f = if zone_count <= 1 {
                0.0
            } else {
                i as f64 / (zone_count as f64 - 1.0)
            };
            // Solar-like radial profile: hot core with rapid falloff toward envelope.
            // T_core ≈ 1.57e7 K (0.0157 T9), photosphere floor ≈ 5.8e-6 T9.
            let baseline = 0.0157 * (-7.2 * f.powf(1.85)).exp() + 5.8e-6;
            temps.push((baseline * t9_scale).clamp(1.0e-5, 5.0));
        }
        temps
    }

    fn build_rate_profile(zone_count: usize) -> Vec<f64> {
        let mut scales = Vec::with_capacity(zone_count);
        for i in 0..zone_count {
            let f = if zone_count <= 1 {
                0.0
            } else {
                i as f64 / (zone_count as f64 - 1.0)
            };
            // Density-inspired burn scaling: reaction power tracks rho^2.
            let rho_rel = (1.0 - f.powf(1.32)).max(0.0).powf(2.6);
            let scale = (rho_rel * rho_rel).max(1.0e-8);
            scales.push(scale);
        }
        scales
    }

    fn rebuild_profile(&mut self) {
        self.burn.cfg.diffusion_coeff = self.diffusion_coeff.max(0.0);
        self.burn.cfg.zone_temperatures_t9 =
            Self::build_temperature_profile(self.zone_count, self.t9_scale);
        self.burn.cfg.zone_rate_scales = Self::build_rate_profile(self.zone_count);
    }

    fn mean_h1_he4(&self) -> (f64, f64) {
        if self.zones.is_empty() {
            return (0.0, 0.0);
        }
        let inv = 1.0 / self.zones.len() as f64;
        let h1 = self.zones.iter().map(|z| z.get(Species::P1)).sum::<f64>() * inv;
        let he4 = self.zones.iter().map(|z| z.get(Species::He4)).sum::<f64>() * inv;
        (h1, he4)
    }

    fn mean_power(&self) -> f64 {
        if self.zones.is_empty() {
            return 0.0;
        }
        self.zones.iter().map(|z| z.thermal_power).sum::<f64>() / self.zones.len() as f64
    }

    fn core_t9(&self) -> f64 {
        self.burn
            .cfg
            .zone_temperatures_t9
            .first()
            .copied()
            .unwrap_or(0.02)
    }

    fn outer_t9(&self) -> f64 {
        self.burn
            .cfg
            .zone_temperatures_t9
            .last()
            .copied()
            .unwrap_or(0.02)
    }

    fn refresh_diagnostics(&mut self, sim_dt_total: f64) {
        let (h1, he4) = self.mean_h1_he4();
        let mut dh_dt = 0.0;
        let mut dhe4_dt = 0.0;
        if sim_dt_total > 0.0 {
            if let Some(prev) = self.prev_h1_mean {
                dh_dt = (h1 - prev) / sim_dt_total;
            }
            if let Some(prev) = self.prev_he4_mean {
                dhe4_dt = (he4 - prev) / sim_dt_total;
            }
        }
        self.prev_h1_mean = Some(h1);
        self.prev_he4_mean = Some(he4);

        let mut sum_x = 0.0;
        let mut max_dev: f64 = 0.0;
        let mut min_x = f64::INFINITY;
        for z in &self.zones {
            let s = z.abund.values().sum::<f64>();
            sum_x += s;
            max_dev = max_dev.max((s - 1.0).abs());
            for v in z.abund.values() {
                min_x = min_x.min(*v);
            }
        }
        let zone_n = self.zones.len().max(1) as f64;
        if !min_x.is_finite() {
            min_x = 0.0;
        }

        let (pp_power, cno_power, triple_alpha_power) = self.channel_powers();
        self.diag = FusionDiagnostics {
            sum_x_mean: sum_x / zone_n,
            sum_x_max_dev: max_dev,
            min_x,
            dh_dt,
            dhe4_dt,
            pp_power,
            cno_power,
            triple_alpha_power,
        };
    }

    fn channel_powers(&self) -> (f64, f64, f64) {
        let mut pp = 0.0;
        let mut cno = 0.0;
        let mut triple_alpha = 0.0;
        for (zi, z) in self.zones.iter().enumerate() {
            let t9 = self
                .burn
                .cfg
                .zone_temperatures_t9
                .get(zi)
                .copied()
                .or_else(|| self.burn.cfg.zone_temperatures_t9.last().copied())
                .unwrap_or(0.02);
            let rate_scale = self
                .burn
                .cfg
                .zone_rate_scales
                .get(zi)
                .copied()
                .or_else(|| self.burn.cfg.zone_rate_scales.last().copied())
                .unwrap_or(1.0)
                .max(0.0);
            for r in &self.burn.core.net.reactions {
                let base_rate =
                    self.burn.core.rates.rate_for(r.id, t9).unwrap_or(0.0) * r.branching_weight;
                let abund_factor = r.reactants.iter().fold(1.0_f64, |acc, st| {
                    acc * z.get(st.species).powi(st.coeff.max(0))
                });
                let flux = base_rate * abund_factor * rate_scale;
                let neutrino_frac = if r
                    .products
                    .iter()
                    .any(|s| s.species == Species::ElectronNeutrino)
                {
                    0.35
                } else {
                    0.0
                };
                let thermal_power = (flux * r.q_mev * (1.0 - neutrino_frac)).max(0.0);
                match r.channel {
                    "pp" => pp += thermal_power,
                    "cno" => cno += thermal_power,
                    "triple_alpha" => triple_alpha += thermal_power,
                    _ => {}
                }
            }
        }
        (pp, cno, triple_alpha)
    }

    fn reset(&mut self) {
        self.zones = self.burn.seed_zones(self.zone_count);
        self.sim_time_s = 0.0;
        self.prev_h1_mean = None;
        self.prev_he4_mean = None;
        self.diag = FusionDiagnostics::default();
        self.refresh_diagnostics(0.0);
    }

    fn step_for_wall_dt(&mut self, wall_dt: Duration) {
        if !self.running {
            return;
        }
        let sim_dt_total = wall_dt.as_secs_f64() * self.wall_to_model_scale * self.speed as f64;
        if sim_dt_total <= 0.0 {
            return;
        }
        let mut remaining = sim_dt_total.min(3.0e8);
        let dt_chunk = 1.0e6_f64;
        let mut guard = 0u32;
        while remaining > 0.0 && guard < 4096 {
            let dt = remaining.min(dt_chunk);
            self.burn.step(&mut self.zones, dt);
            remaining -= dt;
            guard += 1;
        }
        self.sim_time_s += sim_dt_total;
        self.refresh_diagnostics(sim_dt_total);
    }

    fn params(
        &self,
        width: f32,
        height: f32,
        time: f32,
        camera: &FlightCamera,
        quality_tier: f32,
        draw_distance: f32,
        planet_scale: f32,
        planet_time_scale: f32,
    ) -> Params {
        let n = self.zones.len().max(1);
        let i0 = 0usize;
        let i1 = (((n - 1) as f32) * 0.33).round() as usize;
        let i2 = (((n - 1) as f32) * 0.67).round() as usize;
        let i3 = n - 1;
        let z0 = &self.zones[i0.min(n - 1)];
        let z1 = &self.zones[i1.min(n - 1)];
        let z2 = &self.zones[i2.min(n - 1)];
        let z3 = &self.zones[i3.min(n - 1)];
        let w0 = 0.42_f32;
        let w1 = 0.28_f32;
        let w2 = 0.18_f32;
        let w3 = 0.12_f32;
        let h1 = (w0 * z0.get(Species::P1) as f32
            + w1 * z1.get(Species::P1) as f32
            + w2 * z2.get(Species::P1) as f32
            + w3 * z3.get(Species::P1) as f32)
            .clamp(0.0, 1.0);
        let he4 = (w0 * z0.get(Species::He4) as f32
            + w1 * z1.get(Species::He4) as f32
            + w2 * z2.get(Species::He4) as f32
            + w3 * z3.get(Species::He4) as f32)
            .clamp(0.0, 1.0);
        let raw_power = (w0 * z0.thermal_power as f32
            + w1 * z1.thermal_power as f32
            + w2 * z2.thermal_power as f32
            + w3 * z3.thermal_power as f32)
            .max(0.0);
        let power_norm = (raw_power / 2.0e-5).powf(0.55).clamp(0.0, 3.0);
        let zpow0 = (z0.thermal_power as f32 / 2.0e-5)
            .powf(0.55)
            .clamp(0.0, 3.0);
        let zpow1 = (z1.thermal_power as f32 / 2.0e-5)
            .powf(0.55)
            .clamp(0.0, 3.0);
        let zpow2 = (z2.thermal_power as f32 / 2.0e-5)
            .powf(0.55)
            .clamp(0.0, 3.0);
        let zpow3 = (z3.thermal_power as f32 / 2.0e-5)
            .powf(0.55)
            .clamp(0.0, 3.0);
        let (rv, uv, fv) = camera.basis();
        let fov_tan = (camera.fov_deg.to_radians() * 0.5).tan().clamp(0.1, 4.0);
        Params {
            width,
            height,
            time,
            fusion_power: power_norm,
            h1,
            he4,
            running: if self.running { 1.0 } else { 0.0 },
            speed_norm: self.speed.clamp(0.0, 30.0) / 10.0,
            zone_p0: zpow0,
            zone_p1: zpow1,
            zone_p2: zpow2,
            zone_p3: zpow3,
            zone_h10: z0.get(Species::P1) as f32,
            zone_h11: z1.get(Species::P1) as f32,
            zone_h12: z2.get(Species::P1) as f32,
            zone_h13: z3.get(Species::P1) as f32,
            zone_he40: z0.get(Species::He4) as f32,
            zone_he41: z1.get(Species::He4) as f32,
            zone_he42: z2.get(Species::He4) as f32,
            zone_he43: z3.get(Species::He4) as f32,
            zone_count: self.zone_count as f32,
            diffusion_log10: self.diffusion_coeff.max(1e-20).log10() as f32,
            t9_scale: self.t9_scale as f32,
            _zone_pad: 0.0,
            cam_x: camera.pos[0],
            cam_y: camera.pos[1],
            cam_z: camera.pos[2],
            _cam_pad0: 0.0,
            cam_rx: rv[0],
            cam_ry: rv[1],
            cam_rz: rv[2],
            _cam_pad1: 0.0,
            cam_ux: uv[0],
            cam_uy: uv[1],
            cam_uz: uv[2],
            _cam_pad2: 0.0,
            cam_fx: fv[0],
            cam_fy: fv[1],
            cam_fz: fv[2],
            fov_tan,
            quality_tier: quality_tier.clamp(0.0, 2.0),
            draw_distance: draw_distance.max(1.0),
            planet_scale: planet_scale.max(0.0),
            planet_time_scale: planet_time_scale.max(0.0),
            _render_pad0: 0.0,
            _render_pad1: 0.0,
            _render_pad2: 0.0,
        }
    }
}

struct Gpu {
    _instance: wgpu::Instance,
    surface: wgpu::Surface<'static>,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    pipeline: wgpu::RenderPipeline,
    bind_group: wgpu::BindGroup,
    uni_buf: wgpu::Buffer,
}

impl Gpu {
    fn new(window: Arc<Window>) -> Self {
        pollster::block_on(Self::new_async(window))
    }

    async fn new_async(window: Arc<Window>) -> Self {
        let size = window.inner_size();
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::default());
        let surface = instance
            .create_surface(window)
            .expect("failed to create surface");
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .expect("no GPU adapter found");
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor::default(), None)
            .await
            .expect("failed to create wgpu device");

        let caps = surface.get_capabilities(&adapter);
        let format = caps
            .formats
            .iter()
            .copied()
            .find(|f| f.is_srgb())
            .unwrap_or(caps.formats[0]);
        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format,
            width: size.width.max(1),
            height: size.height.max(1),
            present_mode: wgpu::PresentMode::Fifo,
            alpha_mode: caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        let uni_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("fusion_params"),
            size: std::mem::size_of::<Params>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("fusion_bgl"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("fusion_bg"),
            layout: &bgl,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uni_buf.as_entire_binding(),
            }],
        });

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("fusion_shader"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("fusion_pl"),
            bind_group_layouts: &[&bgl],
            push_constant_ranges: &[],
        });
        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("fusion_pipe"),
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
                label: Some("fusion_frame"),
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
            rp.draw(0..6, 0..1);
        }
        self.queue.submit([enc.finish()]);
        frame.present();
        Ok(())
    }
}

struct App {
    window: Option<Arc<Window>>,
    gpu: Option<Gpu>,
    fusion: FusionRuntime,
    camera: FlightCamera,
    quality_tier: f32,
    draw_distance: f32,
    planet_scale: f32,
    planet_time_scale: f32,
    gilrs: Option<Gilrs>,
    pad_prev: PadPrev,
    held: HeldKeys,
    flight_speed: f32,
    last_frame_at: Instant,
    min_frame_dt: Duration,
    start_at: Instant,
    last_fusion_tick_at: Instant,
    last_input_at: Instant,
    avg_frame_ms: f32,
    perf_window_start: Instant,
    perf_accum: Duration,
    perf_frames: u32,
}

impl App {
    fn new() -> Self {
        let parse_f32 = |k: &str, dflt: f32| {
            std::env::var(k)
                .ok()
                .and_then(|s| s.parse::<f32>().ok())
                .unwrap_or(dflt)
        };
        Self {
            window: None,
            gpu: None,
            fusion: FusionRuntime::default(),
            camera: FlightCamera::default(),
            quality_tier: parse_f32("FUSION_VIEWER_QUALITY", 1.0).clamp(0.0, 2.0),
            draw_distance: parse_f32("FUSION_VIEWER_DRAW_DISTANCE", 140.0).clamp(20.0, 4000.0),
            planet_scale: parse_f32("FUSION_VIEWER_PLANET_SCALE", 0.0).clamp(0.0, 20.0),
            planet_time_scale: parse_f32("FUSION_VIEWER_PLANET_TIME_SCALE", 1.0).clamp(0.0, 20.0),
            gilrs: Gilrs::new().ok(),
            pad_prev: PadPrev::default(),
            held: HeldKeys::default(),
            flight_speed: 7.0,
            last_frame_at: Instant::now() - Duration::from_millis(16),
            min_frame_dt: Duration::from_millis(16),
            start_at: Instant::now(),
            last_fusion_tick_at: Instant::now(),
            last_input_at: Instant::now(),
            avg_frame_ms: 0.0,
            perf_window_start: Instant::now(),
            perf_accum: Duration::ZERO,
            perf_frames: 0,
        }
    }

    #[inline]
    fn button_pressed(now: bool, prev: &mut bool) -> bool {
        let fired = now && !*prev;
        *prev = now;
        fired
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

    fn safe_resize(&mut self, w: u32, h: u32) {
        let Some(gpu) = self.gpu.as_mut() else { return };
        let resized = panic::catch_unwind(AssertUnwindSafe(|| gpu.resize(w, h)));
        if resized.is_err() {
            log::warn!("surface resize panicked (likely device lost)");
        }
    }

    fn update_title(&self) {
        let Some(win) = self.window.as_ref() else {
            return;
        };
        let mode = if self.fusion.running { "RUN" } else { "PAUSE" };
        let (h1_mean, he4_mean) = self.fusion.mean_h1_he4();
        let p_mean = self.fusion.mean_power();
        let d = self.fusion.diag;
        let cam_r = (self.camera.pos[0] * self.camera.pos[0]
            + self.camera.pos[1] * self.camera.pos[1]
            + self.camera.pos[2] * self.camera.pos[2])
            .sqrt();
        win.set_title(&format!(
            "Fusion Viewer | {} {:.1}ms | q{} burn x{:.2} zones:{} T9[{:.3}→{:.3}] D={:.1e} | P={:.2e} H={:.4} He={:.4} dH/dt={:+.2e} dHe/dt={:+.2e} | ΣX={:.6} dev={:.2e} minX={:.2e} | Ppp={:.2e} Pcno={:.2e} P3a={:.2e} | cam ({:+.2},{:+.2},{:+.2}) yaw {:+.0}° pitch {:+.0}° roll {:+.0}° v={:.1} | scale r={:.1}R* draw={:.1}R* planets x{:.2} t{:.2}",
            mode,
            self.avg_frame_ms,
            self.quality_tier.round() as i32,
            self.fusion.speed,
            self.fusion.zone_count,
            self.fusion.core_t9(),
            self.fusion.outer_t9(),
            self.fusion.diffusion_coeff,
            p_mean,
            h1_mean,
            he4_mean,
            d.dh_dt,
            d.dhe4_dt,
            d.sum_x_mean,
            d.sum_x_max_dev,
            d.min_x,
            d.pp_power,
            d.cno_power,
            d.triple_alpha_power,
            self.camera.pos[0],
            self.camera.pos[1],
            self.camera.pos[2],
            self.camera.yaw.to_degrees(),
            self.camera.pitch.to_degrees(),
            self.camera.roll.to_degrees(),
            self.flight_speed,
            cam_r / STAR_RADIUS_UNITS,
            self.draw_distance / STAR_RADIUS_UNITS,
            self.planet_scale,
            self.planet_time_scale,
        ));
    }

    fn apply_gamepad(&mut self, el: &ActiveEventLoop, dt: f32) {
        let Some(gilrs) = self.gilrs.as_mut() else {
            return;
        };
        while gilrs.next_event().is_some() {}
        let Some((_, gp)) = gilrs.gamepads().find(|(_, g)| g.is_connected()) else {
            return;
        };
        let mut title_dirty = false;

        if Self::button_pressed(
            gp.is_pressed(Button::RightThumb),
            &mut self.pad_prev.right_thumb,
        ) {
            self.fusion.running = !self.fusion.running;
            title_dirty = true;
        }
        if Self::button_pressed(gp.is_pressed(Button::West), &mut self.pad_prev.west) {
            self.fusion.reset();
            self.fusion.running = true;
            title_dirty = true;
        }
        if Self::button_pressed(gp.is_pressed(Button::Start), &mut self.pad_prev.start) {
            self.fusion.running = !self.fusion.running;
            title_dirty = true;
        }

        let dead = 0.12;
        let lx = Self::smooth_deadzone(gp.value(Axis::LeftStickX), dead);
        let ly = Self::smooth_deadzone(gp.value(Axis::LeftStickY), dead);
        let rx = Self::smooth_deadzone(gp.value(Axis::RightStickX), dead);
        let ry = Self::smooth_deadzone(gp.value(Axis::RightStickY), dead);
        let lt = gp.value(Axis::LeftZ).max(0.0);
        let rt = gp.value(Axis::RightZ).max(0.0);
        let turbo = if gp.is_pressed(Button::South) {
            3.0
        } else {
            1.0
        };
        let speed = self.flight_speed * turbo;

        if lx != 0.0 || ly != 0.0 || (rt - lt).abs() > 0.01 {
            self.camera
                .move_local(-ly, lx, rt - lt, speed, dt.max(1e-4));
            title_dirty = true;
        }
        if rx != 0.0 || ry != 0.0 {
            self.camera.look(rx * dt * 2.3, -ry * dt * 1.8);
            title_dirty = true;
        }
        if gp.is_pressed(Button::LeftTrigger) {
            self.camera.roll -= dt * 1.8;
            title_dirty = true;
        }
        if gp.is_pressed(Button::RightTrigger) {
            self.camera.roll += dt * 1.8;
            title_dirty = true;
        }

        if gp.is_pressed(Button::DPadUp) {
            self.fusion.speed = (self.fusion.speed * 1.04).clamp(0.05, 100.0);
            title_dirty = true;
        }
        if gp.is_pressed(Button::DPadDown) {
            self.fusion.speed = (self.fusion.speed / 1.04).clamp(0.05, 100.0);
            title_dirty = true;
        }

        if gp.is_pressed(Button::DPadRight) {
            self.fusion.t9_scale = (self.fusion.t9_scale * 1.01).clamp(0.05, 8.0);
            self.fusion.rebuild_profile();
            title_dirty = true;
        }
        if gp.is_pressed(Button::DPadLeft) {
            self.fusion.t9_scale = (self.fusion.t9_scale / 1.01).clamp(0.05, 8.0);
            self.fusion.rebuild_profile();
            title_dirty = true;
        }
        if gp.is_pressed(Button::North) {
            self.draw_distance = (self.draw_distance * 1.015).clamp(20.0, 4000.0);
            title_dirty = true;
        }
        if gp.is_pressed(Button::East) {
            self.draw_distance = (self.draw_distance / 1.015).clamp(20.0, 4000.0);
            title_dirty = true;
        }
        if lt > 0.85
            && rt > 0.85
            && Self::button_pressed(gp.is_pressed(Button::Mode), &mut self.pad_prev.mode)
        {
            el.exit();
        }
        if title_dirty {
            self.update_title();
        }
    }

    fn push_frame(&mut self) {
        let now = Instant::now();
        if now.duration_since(self.last_frame_at) < self.min_frame_dt {
            return;
        }
        self.last_frame_at = now;

        let fusion_dt = now.saturating_duration_since(self.last_fusion_tick_at);
        self.last_fusion_tick_at = now;
        self.fusion.step_for_wall_dt(fusion_dt);

        let Some(win) = self.window.as_ref() else {
            return;
        };
        let sz = win.inner_size();
        let params = self.fusion.params(
            sz.width.max(1) as f32,
            sz.height.max(1) as f32,
            self.start_at.elapsed().as_secs_f32(),
            &self.camera,
            self.quality_tier,
            self.draw_distance,
            self.planet_scale,
            self.planet_time_scale,
        );

        let frame_start = Instant::now();
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
                log::error!("surface out of memory");
            }
            Err(e) => log::error!("render error: {e}"),
        }

        self.perf_accum += frame_start.elapsed();
        self.perf_frames += 1;
        if self.perf_frames >= 45 || self.perf_window_start.elapsed() >= Duration::from_secs(2) {
            self.avg_frame_ms =
                (self.perf_accum.as_secs_f64() * 1000.0 / (self.perf_frames.max(1) as f64)) as f32;
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
                    .with_title("GUTOE Fusion Viewer")
                    .with_inner_size(winit::dpi::LogicalSize::new(1000u32, 720u32)),
            )
            .unwrap(),
        );
        let gpu = Gpu::new(Arc::clone(&win));
        self.window = Some(win);
        self.gpu = Some(gpu);
        self.last_fusion_tick_at = Instant::now();
        self.last_input_at = Instant::now();
        self.update_title();
        self.push_frame();
    }

    fn window_event(&mut self, el: &ActiveEventLoop, _id: WindowId, event: WindowEvent) {
        match event {
            WindowEvent::CloseRequested => el.exit(),
            WindowEvent::Resized(size) => {
                self.safe_resize(size.width, size.height);
                self.push_frame();
            }
            WindowEvent::RedrawRequested => self.push_frame(),
            WindowEvent::KeyboardInput {
                event: KeyEvent {
                    logical_key, state, ..
                },
                ..
            } => {
                let pressed = state == ElementState::Pressed;
                match logical_key {
                    Key::Named(NamedKey::Escape) if pressed => el.exit(),
                    Key::Named(NamedKey::Shift) => {
                        self.held.turbo = pressed;
                    }
                    Key::Named(NamedKey::ArrowLeft) => self.held.yaw_left = pressed,
                    Key::Named(NamedKey::ArrowRight) => self.held.yaw_right = pressed,
                    Key::Named(NamedKey::ArrowUp) => self.held.pitch_up = pressed,
                    Key::Named(NamedKey::ArrowDown) => self.held.pitch_down = pressed,
                    Key::Named(NamedKey::Space) if pressed => {
                        self.fusion.running = !self.fusion.running;
                        self.update_title();
                        self.push_frame();
                    }
                    Key::Character(ref s) => match s.as_str() {
                        "q" | "Q" if pressed => el.exit(),
                        "w" | "W" => self.held.fwd = pressed,
                        "s" | "S" => self.held.back = pressed,
                        "a" | "A" => self.held.left = pressed,
                        "d" | "D" => self.held.right = pressed,
                        "e" | "E" => self.held.up = pressed,
                        "c" | "C" => self.held.down = pressed,
                        "z" | "Z" => self.held.roll_left = pressed,
                        "x" | "X" => self.held.roll_right = pressed,
                        "r" | "R" if pressed => {
                            self.fusion.reset();
                            self.fusion.running = true;
                            self.camera = FlightCamera::default();
                            self.update_title();
                            self.push_frame();
                        }
                        "[" if pressed => {
                            self.flight_speed = (self.flight_speed / 1.15).clamp(0.5, 80.0);
                            self.update_title();
                        }
                        "]" if pressed => {
                            self.flight_speed = (self.flight_speed * 1.15).clamp(0.5, 80.0);
                            self.update_title();
                        }
                        "j" | "J" if pressed => {
                            self.fusion.speed = (self.fusion.speed / 1.25).clamp(0.05, 100.0);
                            self.update_title();
                            self.push_frame();
                        }
                        "l" | "L" if pressed => {
                            self.fusion.speed = (self.fusion.speed * 1.25).clamp(0.05, 100.0);
                            self.update_title();
                            self.push_frame();
                        }
                        "u" | "U" if pressed => {
                            self.fusion.t9_scale = (self.fusion.t9_scale / 1.1).clamp(0.05, 8.0);
                            self.fusion.rebuild_profile();
                            self.update_title();
                            self.push_frame();
                        }
                        "i" | "I" if pressed => {
                            self.fusion.t9_scale = (self.fusion.t9_scale * 1.1).clamp(0.05, 8.0);
                            self.fusion.rebuild_profile();
                            self.update_title();
                            self.push_frame();
                        }
                        "o" | "O" if pressed => {
                            self.fusion.diffusion_coeff =
                                (self.fusion.diffusion_coeff / 1.25).max(1.0e-16);
                            self.fusion.rebuild_profile();
                            self.update_title();
                            self.push_frame();
                        }
                        "p" | "P" if pressed => {
                            self.fusion.diffusion_coeff =
                                (self.fusion.diffusion_coeff * 1.25).min(1.0e-4);
                            self.fusion.rebuild_profile();
                            self.update_title();
                            self.push_frame();
                        }
                        "-" if pressed => {
                            self.camera.fov_deg = (self.camera.fov_deg * 1.08).clamp(35.0, 120.0);
                            self.update_title();
                        }
                        "=" | "+" if pressed => {
                            self.camera.fov_deg = (self.camera.fov_deg / 1.08).clamp(35.0, 120.0);
                            self.update_title();
                        }
                        "," | "<" if pressed => {
                            self.draw_distance = (self.draw_distance / 1.2).clamp(20.0, 4000.0);
                            self.update_title();
                        }
                        "." | ">" if pressed => {
                            self.draw_distance = (self.draw_distance * 1.2).clamp(20.0, 4000.0);
                            self.update_title();
                        }
                        "n" | "N" if pressed => {
                            self.planet_scale = (self.planet_scale / 1.1).clamp(0.0, 20.0);
                            self.update_title();
                        }
                        "m" | "M" if pressed => {
                            self.planet_scale = (self.planet_scale * 1.1).clamp(0.1, 20.0);
                            self.update_title();
                        }
                        "k" | "K" if pressed => {
                            self.planet_time_scale =
                                (self.planet_time_scale / 1.2).clamp(0.0, 20.0);
                            self.update_title();
                        }
                        ";" | ":" if pressed => {
                            self.planet_time_scale =
                                (self.planet_time_scale * 1.2).clamp(0.0, 20.0);
                            self.update_title();
                        }
                        "1" if pressed => {
                            self.quality_tier = 0.0;
                            self.update_title();
                        }
                        "2" if pressed => {
                            self.quality_tier = 1.0;
                            self.update_title();
                        }
                        "3" if pressed => {
                            self.quality_tier = 2.0;
                            self.update_title();
                        }
                        _ => {}
                    },
                    _ => {}
                }
            }
            _ => {}
        }
    }

    fn about_to_wait(&mut self, el: &ActiveEventLoop) {
        let now = Instant::now();
        let dt = now
            .saturating_duration_since(self.last_input_at)
            .as_secs_f32()
            .clamp(1.0 / 1000.0, 0.1);
        self.last_input_at = now;

        self.apply_gamepad(el, dt);

        let turbo = if self.held.turbo { 3.0 } else { 1.0 };
        let speed = self.flight_speed * turbo;
        let forward = (self.held.fwd as i32 - self.held.back as i32) as f32;
        let right = (self.held.right as i32 - self.held.left as i32) as f32;
        let up = (self.held.up as i32 - self.held.down as i32) as f32;
        if forward != 0.0 || right != 0.0 || up != 0.0 {
            self.camera.move_local(forward, right, up, speed, dt);
            self.update_title();
        }
        let yaw = (self.held.yaw_right as i32 - self.held.yaw_left as i32) as f32;
        let pitch = (self.held.pitch_up as i32 - self.held.pitch_down as i32) as f32;
        if yaw != 0.0 || pitch != 0.0 {
            self.camera.look(yaw * dt * 1.8, pitch * dt * 1.6);
            self.update_title();
        }
        let roll = (self.held.roll_right as i32 - self.held.roll_left as i32) as f32;
        if roll != 0.0 {
            self.camera.roll += roll * dt * 1.8;
            self.update_title();
        }

        if let Some(win) = self.window.as_ref() {
            win.request_redraw();
        }
    }
}

fn main() {
    env_logger::init();

    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  GUTOE Fusion Viewer");
    println!("  Real-time multi-zone radial-shell burn + live 3D wgpu flight");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  Space              — run / pause");
    println!("  R                  — reset composition (solar-like seed)");
    println!("  WASD               — 3D flight (forward/strafe)");
    println!("  E / C              — move up / down");
    println!("  Arrows             — look yaw / pitch");
    println!("  Z / X              — roll left / right");
    println!("  [ / ]              — flight speed down / up");
    println!("  - / +              — FOV out / in");
    println!("  1 / 2 / 3          — quality tier q0/q1/q2");
    println!("  , / .              — draw distance down / up");
    println!("  N / M              — planet scale down / up (0 disables)");
    println!("  K / ;              — planet orbit speed down / up");
    println!("  J / L              — speed down / up");
    println!("  U / I              — temperature profile scale down / up");
    println!("  O / P              — radial diffusion down / up");
    println!("  Shift              — flight turbo");
    println!("  Q / Escape         — quit");
    println!("  Gamepad:");
    println!("    LS               — 3D flight strafe/forward");
    println!("    RS               — look yaw/pitch");
    println!("    LT/RT            — move down/up");
    println!("    LB/RB            — roll left/right");
    println!("    A hold           — flight turbo");
    println!("    R3 or Start      — run/pause");
    println!("    X                — reset composition");
    println!("    DPad Up/Down     — burn speed up/down");
    println!("    DPad Left/Right  — profile temperature down/up");
    println!("    Y / B            — draw distance up/down");
    println!("    LT+RT+Guide      — quit");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let ev = match EventLoop::new() {
        Ok(ev) => ev,
        Err(e) => {
            eprintln!("failed to create event loop: {e}");
            return;
        }
    };
    if let Err(e) = ev.run_app(&mut App::new()) {
        eprintln!("fusion_viewer runtime error: {e}");
    }
}
