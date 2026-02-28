use image::{Rgb, RgbImage};
use serde_json::Value;
use std::env;
use std::f32::consts::PI;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug)]
struct V3 {
    x: f32,
    y: f32,
    z: f32,
}

impl V3 {
    fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }
    fn dot(self, o: Self) -> f32 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }
    fn len(self) -> f32 {
        self.dot(self).sqrt()
    }
    fn norm(self) -> Self {
        let l = self.len().max(1e-8);
        Self::new(self.x / l, self.y / l, self.z / l)
    }
    fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }
    fn clamp01(self) -> Self {
        Self::new(
            self.x.clamp(0.0, 1.0),
            self.y.clamp(0.0, 1.0),
            self.z.clamp(0.0, 1.0),
        )
    }
}

use std::ops::{Add, AddAssign, Div, Mul, Sub};
impl Add for V3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        V3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}
impl AddAssign for V3 {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}
impl Sub for V3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        V3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}
impl Mul<f32> for V3 {
    type Output = Self;
    fn mul(self, rhs: f32) -> Self::Output {
        V3::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}
impl Mul for V3 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        V3::new(self.x * rhs.x, self.y * rhs.y, self.z * rhs.z)
    }
}
impl Div<f32> for V3 {
    type Output = Self;
    fn div(self, rhs: f32) -> Self::Output {
        V3::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

fn mix(a: V3, b: V3, t: f32) -> V3 {
    a * (1.0 - t) + b * t
}

#[derive(Clone, Debug)]
struct Candidate {
    rank: usize,
    z: u16,
    symbol: String,
    family: String,
    stable_like: usize,
}

fn load_candidates(path: &Path) -> Result<Vec<Candidate>, String> {
    let txt = fs::read_to_string(path)
        .map_err(|e| format!("failed to read witness json {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&txt)
        .map_err(|e| format!("failed to parse witness json {}: {e}", path.display()))?;
    let arr = v
        .get("candidates")
        .and_then(|x| x.as_array())
        .ok_or_else(|| "missing candidates[] in witness json".to_string())?;
    let mut out = Vec::new();
    for item in arr {
        let rank = item.get("rank").and_then(Value::as_u64).unwrap_or(0) as usize;
        let z = item.get("z").and_then(Value::as_u64).unwrap_or(0) as u16;
        let symbol = item
            .get("symbol")
            .and_then(Value::as_str)
            .unwrap_or("?")
            .to_string();
        let family = item
            .get("family")
            .and_then(Value::as_str)
            .unwrap_or("unknown")
            .to_string();
        let stable_like = item
            .get("stable_like")
            .and_then(Value::as_u64)
            .unwrap_or(0) as usize;
        out.push(Candidate {
            rank,
            z,
            symbol,
            family,
            stable_like,
        });
    }
    out.sort_by_key(|c| c.rank);
    Ok(out)
}

fn hash11(x: f32) -> f32 {
    (x.sin() * 43_758.547).fract().abs()
}

fn hash31(p: V3) -> f32 {
    hash11(p.x * 127.1 + p.y * 311.7 + p.z * 74.7)
}

fn noise3(p: V3) -> f32 {
    let i = V3::new(p.x.floor(), p.y.floor(), p.z.floor());
    let f = V3::new(p.x - i.x, p.y - i.y, p.z - i.z);
    let u = V3::new(
        f.x * f.x * (3.0 - 2.0 * f.x),
        f.y * f.y * (3.0 - 2.0 * f.y),
        f.z * f.z * (3.0 - 2.0 * f.z),
    );

    let c000 = hash31(i + V3::new(0.0, 0.0, 0.0));
    let c100 = hash31(i + V3::new(1.0, 0.0, 0.0));
    let c010 = hash31(i + V3::new(0.0, 1.0, 0.0));
    let c110 = hash31(i + V3::new(1.0, 1.0, 0.0));
    let c001 = hash31(i + V3::new(0.0, 0.0, 1.0));
    let c101 = hash31(i + V3::new(1.0, 0.0, 1.0));
    let c011 = hash31(i + V3::new(0.0, 1.0, 1.0));
    let c111 = hash31(i + V3::new(1.0, 1.0, 1.0));

    let x00 = c000 * (1.0 - u.x) + c100 * u.x;
    let x10 = c010 * (1.0 - u.x) + c110 * u.x;
    let x01 = c001 * (1.0 - u.x) + c101 * u.x;
    let x11 = c011 * (1.0 - u.x) + c111 * u.x;
    let y0 = x00 * (1.0 - u.y) + x10 * u.y;
    let y1 = x01 * (1.0 - u.y) + x11 * u.y;
    y0 * (1.0 - u.z) + y1 * u.z
}

fn fbm(mut p: V3) -> f32 {
    let mut amp = 0.5;
    let mut sum = 0.0;
    for _ in 0..5 {
        sum += amp * noise3(p);
        p = p * 2.03 + V3::new(7.1, 9.2, 3.7);
        amp *= 0.5;
    }
    sum
}

fn nearest_delta(v: f32, center: f32) -> f32 {
    let x = v - center;
    x - x.round()
}

// BCC(100)-like projected motif:
// - corner sites at integer lattice points
// - body-centered sites at half-shifted points
fn bcc_surface_peak(u: f32, v: f32) -> f32 {
    let dx_c = nearest_delta(u, 0.0);
    let dz_c = nearest_delta(v, 0.0);
    let d2_c = dx_c * dx_c + dz_c * dz_c;

    let dx_b = nearest_delta(u, 0.5);
    let dz_b = nearest_delta(v, 0.5);
    let d2_b = dx_b * dx_b + dz_b * dz_b;

    let corner = (-d2_c / (2.0 * 0.17 * 0.17)).exp();
    let body = 0.72 * (-d2_b / (2.0 * 0.20 * 0.20)).exp();
    corner + body
}

fn rotate_y(p: V3, a: f32) -> V3 {
    let c = a.cos();
    let s = a.sin();
    V3::new(c * p.x + s * p.z, p.y, -s * p.x + c * p.z)
}

fn rotate_x(p: V3, a: f32) -> V3 {
    let c = a.cos();
    let s = a.sin();
    V3::new(p.x, c * p.y - s * p.z, s * p.y + c * p.z)
}

fn to_local(p: V3, center: V3, yaw: f32, pitch: f32) -> V3 {
    rotate_x(rotate_y(p - center, -yaw), -pitch)
}

fn to_world_dir(p: V3, yaw: f32, pitch: f32) -> V3 {
    rotate_y(rotate_x(p, pitch), yaw)
}

fn intersect_obb(ro: V3, rd: V3, center: V3, half: V3, yaw: f32, pitch: f32) -> Option<(f32, V3, V3)> {
    let ro_l = to_local(ro, center, yaw, pitch);
    let rd_l = rotate_x(rotate_y(rd, -yaw), -pitch);

    let mut tmin = -1.0e20_f32;
    let mut tmax = 1.0e20_f32;

    for (ro_i, rd_i, h_i) in [
        (ro_l.x, rd_l.x, half.x),
        (ro_l.y, rd_l.y, half.y),
        (ro_l.z, rd_l.z, half.z),
    ] {
        if rd_i.abs() < 1e-7 {
            if ro_i < -h_i || ro_i > h_i {
                return None;
            }
        } else {
            let t1 = (-h_i - ro_i) / rd_i;
            let t2 = (h_i - ro_i) / rd_i;
            let tn = t1.min(t2);
            let tf = t1.max(t2);
            tmin = tmin.max(tn);
            tmax = tmax.min(tf);
            if tmin > tmax {
                return None;
            }
        }
    }

    let t_hit = if tmin > 1e-4 { tmin } else { tmax };
    if t_hit <= 1e-4 {
        return None;
    }

    let hp_l = ro_l + rd_l * t_hit;
    let dx = (hp_l.x.abs() - half.x).abs();
    let dy = (hp_l.y.abs() - half.y).abs();
    let dz = (hp_l.z.abs() - half.z).abs();
    let n_l = if dx <= dy && dx <= dz {
        V3::new(hp_l.x.signum(), 0.0, 0.0)
    } else if dy <= dz {
        V3::new(0.0, hp_l.y.signum(), 0.0)
    } else {
        V3::new(0.0, 0.0, hp_l.z.signum())
    };

    let n_w = to_world_dir(n_l, yaw, pitch).norm();
    let hp_w = ro + rd * t_hit;
    Some((t_hit, hp_w, n_w))
}

fn env_color(rd: V3) -> V3 {
    let t = (0.5 * (rd.y + 1.0)).clamp(0.0, 1.0);
    let top = V3::new(0.03, 0.05, 0.09);
    let bot = V3::new(0.10, 0.10, 0.11);
    mix(bot, top, t)
}

fn element_base_color(symbol: &str) -> V3 {
    match symbol {
        "Pt" => V3::new(0.90, 0.90, 0.92),
        "Hf" => V3::new(0.76, 0.78, 0.82),
        "Mo" => V3::new(0.66, 0.69, 0.73),
        "Zn" => V3::new(0.75, 0.80, 0.90),
        "Cd" => V3::new(0.80, 0.84, 0.93),
        "Cr" => V3::new(0.63, 0.68, 0.72),
        _ => V3::new(0.75, 0.77, 0.82),
    }
}

fn finish_profile(symbol: &str) -> (f32, f32) {
    match symbol {
        "Pt" => (0.07, 0.99), // mirror-ish precious metal
        "Cr" => (0.08, 0.98), // chrome-like
        "Hf" => (0.12, 0.97),
        "Mo" => (0.15, 0.96),
        "Zn" => (0.19, 0.95), // more diffuse galvanized look
        "Cd" => (0.16, 0.95),
        _ => (0.14, 0.96),
    }
}

fn fresnel_schlick(cos_theta: f32, f0: V3) -> V3 {
    f0 + (V3::new(1.0, 1.0, 1.0) - f0) * (1.0 - cos_theta).powf(5.0)
}

fn shade_material(
    c: &Candidate,
    p: V3,
    n_geom: V3,
    rd: V3,
    local_seed: f32,
    tangent_hint: V3,
) -> V3 {
    let mut n = n_geom;
    let t = tangent_hint.norm();
    let b = n.cross(t).norm();

    let a_tex = 0.24 + 0.010 * (c.stable_like as f32); // lattice spacing proxy in world units
    let u = p.dot(t) / a_tex;
    let vlat = p.dot(b) / a_tex;
    let peak = bcc_surface_peak(u, vlat);
    let du = 0.02;
    let dv = 0.02;
    let grad_u = (bcc_surface_peak(u + du, vlat) - bcc_surface_peak(u - du, vlat)) / (2.0 * du);
    let grad_v = (bcc_surface_peak(u, vlat + dv) - bcc_surface_peak(u, vlat - dv)) / (2.0 * dv);

    // Add small thermal/dislocation disorder over the ideal motif.
    let defect = fbm(p * 18.0 + V3::new(local_seed, 2.0, 7.0)) - 0.5;
    n = (n + t * (0.20 * grad_u + 0.08 * defect) + b * (0.20 * grad_v + 0.08 * defect)).norm();

    let mut albedo = element_base_color(&c.symbol);
    let grain_tint = 0.88 + 0.24 * peak.clamp(0.0, 1.0);
    albedo = (albedo * grain_tint).clamp01();

    let (rough0, metallic) = finish_profile(&c.symbol);
    let rough_base =
        (rough0 + 0.015 * ((8usize.saturating_sub(c.stable_like)) as f32)).clamp(0.06, 0.24);
    let roughness = (rough_base + 0.14 * (1.0 - peak.clamp(0.0, 1.0)) + 0.06 * defect.abs())
        .clamp(0.08, 0.42);

    let v = (rd * -1.0).norm();
    let lights = [
        (V3::new(0.62, 0.90, 0.35).norm(), V3::new(4.1, 3.9, 3.6)),
        (V3::new(-0.65, 0.48, -0.45).norm(), V3::new(1.7, 1.8, 2.1)),
        (V3::new(0.10, 0.20, -1.0).norm(), V3::new(0.8, 0.85, 1.1)),
    ];

    let f0_d = 0.04;
    let f0 = mix(V3::new(f0_d, f0_d, f0_d), albedo, metallic);
    let mut direct = V3::new(0.0, 0.0, 0.0);

    for (l, lcol) in lights {
        let h = (v + l).norm();
        let ndotl = n.dot(l).clamp(0.0, 1.0);
        if ndotl <= 0.0 {
            continue;
        }
        let ndotv = n.dot(v).clamp(0.0, 1.0);
        let ndoth = n.dot(h).clamp(0.0, 1.0);
        let vdoth = v.dot(h).clamp(0.0, 1.0);

        let a = (roughness * roughness).max(0.02);
        let a2 = a * a;
        let denom = ndoth * ndoth * (a2 - 1.0) + 1.0;
        let d = a2 / (PI * denom * denom).max(1e-7);
        let k = (roughness + 1.0).powi(2) / 8.0;
        let gv = ndotv / (ndotv * (1.0 - k) + k).max(1e-7);
        let gl = ndotl / (ndotl * (1.0 - k) + k).max(1e-7);
        let g = gv * gl;
        let f = fresnel_schlick(vdoth, f0);

        let spec = f * (d * g / (4.0 * ndotv * ndotl).max(1e-6));
        let diff = albedo * ((1.0 - metallic) / PI);
        direct += (diff + spec) * ndotl * lcol;
    }

    let refl = (v - n * (2.0 * v.dot(n))).norm();
    let env = env_color(refl);
    let rim = (1.0 - n.dot(v).clamp(0.0, 1.0)).powf(2.4);
    let ambient = albedo * 0.04 + env * 0.22 + V3::new(0.28, 0.30, 0.36) * (0.16 * rim);

    (ambient + direct).clamp01()
}

fn render_panel(c: &Candidate, w: usize, h: usize) -> Vec<V3> {
    let mut out = vec![V3::new(0.0, 0.0, 0.0); w * h];
    let aspect = w as f32 / h as f32;
    let ro = V3::new(0.0, 0.18, 1.55);

    let yaw = (12.0 + (c.z % 11) as f32).to_radians();
    let pitch = (8.0 + (c.z % 7) as f32).to_radians();
    // Macro coupon: a thin slab filling the frame, not a staged hero object.
    let box_center = V3::new(0.0, -0.02, 0.0);
    let rank_scale = 1.0 + 0.03 * (((c.rank - 1) % 3) as f32 - 1.0);
    let box_half = V3::new(1.25 * rank_scale, 0.07, 1.05);

    let seed = c.z as f32 * 0.173;
    let tangent_hint = rotate_y(V3::new(1.0, 0.0, 0.0), yaw * 0.6 + pitch * 0.4);

    for py in 0..h {
        for px in 0..w {
            let u = ((px as f32 + 0.5) / w as f32) * 2.0 - 1.0;
            let v = ((py as f32 + 0.5) / h as f32) * 2.0 - 1.0;
            let rd = V3::new(u * aspect * 0.62, -v * 0.52 - 0.14, -1.15).norm();

            let mut best_t = 1.0e9f32;
            let mut col = env_color(rd);

            if let Some((t, hp, n)) = intersect_obb(ro, rd, box_center, box_half, yaw, pitch) {
                if t < best_t {
                    best_t = t;
                    col = shade_material(c, hp, n, rd, seed, tangent_hint);
                }
            }

            // Very light atmospheric rolloff to avoid synthetic haze.
            let fog = (-(0.006 * best_t.min(6.0))).exp();
            col = col * fog + env_color(rd) * (1.0 - fog);
            out[py * w + px] = col;
        }
    }

    out
}

fn linear_to_srgb_u8(c: V3) -> [u8; 3] {
    let mapped = V3::new(
        c.x / (1.0 + c.x),
        c.y / (1.0 + c.y),
        c.z / (1.0 + c.z),
    );
    let g = V3::new(
        mapped.x.clamp(0.0, 1.0).powf(1.0 / 2.2),
        mapped.y.clamp(0.0, 1.0).powf(1.0 / 2.2),
        mapped.z.clamp(0.0, 1.0).powf(1.0 / 2.2),
    );
    [
        (g.x * 255.0 + 0.5) as u8,
        (g.y * 255.0 + 0.5) as u8,
        (g.z * 255.0 + 0.5) as u8,
    ]
}

fn main() -> Result<(), String> {
    let mut out_dir = PathBuf::from(
        env::var("GUTOE_RTSC_RENDER_OUT")
            .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_witness_macro".to_string()),
    );
    let witness_path = PathBuf::from(
        env::var("GUTOE_RTSC_WITNESS_JSON").unwrap_or_else(|_| {
            "/tmp/bh_renders/rtsc_witness_candidates/rtsc_forced_witnesses.json".to_string()
        }),
    );

    if let Some(arg) = env::args().skip(1).next() {
        if arg == "--help" || arg == "-h" {
            println!(
                "Usage: rtsc_witness_macro_render [OUT_DIR]\n\
                 Env overrides:\n\
                   GUTOE_RTSC_RENDER_OUT   (default /tmp/bh_renders/rtsc_witness_macro)\n\
                   GUTOE_RTSC_WITNESS_JSON (default /tmp/bh_renders/rtsc_witness_candidates/rtsc_forced_witnesses.json)\n"
            );
            return Ok(());
        }
        out_dir = PathBuf::from(arg);
    }
    fs::create_dir_all(&out_dir)
        .map_err(|e| format!("failed to create out dir {}: {e}", out_dir.display()))?;

    let candidates = load_candidates(&witness_path)?;
    if candidates.is_empty() {
        return Err("no RTSC witness candidates found".to_string());
    }

    // 4K UHD layout: 3 columns x 2 rows.
    let cols = 3usize;
    let rows = 2usize;
    let panel_w = 1280usize;
    let panel_h = 1080usize;
    let width = panel_w * cols;
    let height = panel_h * rows;

    let mut img = RgbImage::new(width as u32, height as u32);

    for y in 0..height {
        for x in 0..width {
            img.put_pixel(x as u32, y as u32, Rgb([8, 10, 13]));
        }
    }

    for (i, c) in candidates.iter().take(cols * rows).enumerate() {
        let row = i / cols;
        let col = i % cols;
        let x0 = col * panel_w;
        let y0 = row * panel_h;
        let panel = render_panel(c, panel_w, panel_h);
        for py in 0..panel_h {
            for px in 0..panel_w {
                let mut rgb = linear_to_srgb_u8(panel[py * panel_w + px]);
                // subtle vignette per panel for macro-lens look
                let fu = (px as f32 / panel_w as f32) * 2.0 - 1.0;
                let fv = (py as f32 / panel_h as f32) * 2.0 - 1.0;
                let vig = (1.0 - 0.12 * (fu * fu + fv * fv)).clamp(0.76, 1.0);
                rgb[0] = (rgb[0] as f32 * vig) as u8;
                rgb[1] = (rgb[1] as f32 * vig) as u8;
                rgb[2] = (rgb[2] as f32 * vig) as u8;
                img.put_pixel((x0 + px) as u32, (y0 + py) as u32, Rgb(rgb));
            }
        }

        // panel separators
        if col > 0 {
            for py in 0..panel_h {
                img.put_pixel(x0 as u32, (y0 + py) as u32, Rgb([36, 40, 48]));
            }
        }
        if row > 0 {
            for px in 0..panel_w {
                img.put_pixel((x0 + px) as u32, y0 as u32, Rgb([36, 40, 48]));
            }
        }
    }

    let png_path = out_dir.join("rtsc_witness_macro_4k.png");
    img.save(&png_path)
        .map_err(|e| format!("failed to save png {}: {e}", png_path.display()))?;

    let mut meta = String::new();
    meta.push_str("{\n");
    meta.push_str("  \"mode\": \"macro_scale_material_render\",\n");
    meta.push_str("  \"resolution\": {\"width\": 3840, \"height\": 2160},\n");
    meta.push_str(&format!(
        "  \"witness_source\": \"{}\",\n",
        witness_path.display()
    ));
    meta.push_str("  \"panels\": [\n");
    for (i, c) in candidates.iter().take(cols * rows).enumerate() {
        meta.push_str(&format!(
            "    {{\"panel\": {}, \"rank\": {}, \"z\": {}, \"symbol\": \"{}\", \"family\": \"{}\", \"stable_like\": {}}}{}\n",
            i,
            c.rank,
            c.z,
            c.symbol,
            c.family,
            c.stable_like,
            if i + 1 == candidates.len().min(cols * rows) { "" } else { "," }
        ));
    }
    meta.push_str("  ]\n}\n");

    let json_path = out_dir.join("rtsc_witness_macro_4k.json");
    fs::write(&json_path, meta)
        .map_err(|e| format!("failed to write metadata {}: {e}", json_path.display()))?;

    println!("wrote {}", png_path.display());
    println!("wrote {}", json_path.display());
    println!("rendered {} macro panels", candidates.len().min(cols * rows));

    Ok(())
}
