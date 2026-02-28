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
    fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }
    fn len(self) -> f32 {
        self.dot(self).sqrt()
    }
    fn norm(self) -> Self {
        let l = self.len().max(1.0e-8);
        self / l
    }
    fn clamp01(self) -> Self {
        Self::new(
            self.x.clamp(0.0, 1.0),
            self.y.clamp(0.0, 1.0),
            self.z.clamp(0.0, 1.0),
        )
    }
}

use std::ops::{Add, Div, Mul, Sub};
impl Add for V3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        V3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
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

fn smoothstep(edge0: f32, edge1: f32, x: f32) -> f32 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

#[derive(Clone, Debug)]
struct ReactorParams {
    material: String,
    r_major_m: f32,
    a_minor_m: f32,
    mesh_quality: f32,
    b_operating_t: f32,
    t_kev: f32,
    p_net_w: f32,
    q_engineering: f32,
}

#[derive(Clone, Copy, Debug)]
struct Camera {
    eye: V3,
    target: V3,
    up: V3,
    fov_deg: f32,
}

fn load_best(path: &Path) -> Result<ReactorParams, String> {
    let txt = fs::read_to_string(path)
        .map_err(|e| format!("failed reading reactor json {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&txt)
        .map_err(|e| format!("failed parsing reactor json {}: {e}", path.display()))?;
    let b = v
        .get("best_overall")
        .and_then(|x| x.as_object())
        .ok_or_else(|| "missing best_overall".to_string())?;

    Ok(ReactorParams {
        material: b
            .get("material")
            .and_then(|x| x.as_str())
            .unwrap_or("Cr")
            .to_string(),
        r_major_m: b.get("r_major_m").and_then(|x| x.as_f64()).unwrap_or(1.5) as f32,
        a_minor_m: b.get("a_minor_m").and_then(|x| x.as_f64()).unwrap_or(0.6) as f32,
        mesh_quality: b
            .get("mesh_quality")
            .and_then(|x| x.as_f64())
            .unwrap_or(1.0) as f32,
        b_operating_t: b
            .get("b_operating_t")
            .and_then(|x| x.as_f64())
            .unwrap_or(25.0) as f32,
        t_kev: b.get("t_kev").and_then(|x| x.as_f64()).unwrap_or(150.0) as f32,
        p_net_w: b.get("p_net_w").and_then(|x| x.as_f64()).unwrap_or(2.0e7) as f32,
        q_engineering: b
            .get("q_engineering")
            .and_then(|x| x.as_f64())
            .unwrap_or(1.8) as f32,
    })
}

fn material_color(symbol: &str) -> V3 {
    match symbol {
        "Cr" => V3::new(0.74, 0.77, 0.80),
        "Mo" => V3::new(0.70, 0.72, 0.76),
        "Pt" => V3::new(0.90, 0.90, 0.92),
        "Hf" => V3::new(0.77, 0.79, 0.83),
        "Zn" => V3::new(0.76, 0.80, 0.90),
        "Cd" => V3::new(0.80, 0.84, 0.94),
        _ => V3::new(0.75, 0.78, 0.82),
    }
}

fn hash2(x: i32, y: i32) -> f32 {
    let f = ((x as f32) * 12.9898 + (y as f32) * 78.233).sin() * 43_758.547;
    f.fract().abs()
}

fn background(x: usize, y: usize, w: usize, h: usize) -> V3 {
    let v = y as f32 / h as f32;
    let top = V3::new(0.020, 0.030, 0.055);
    let bot = V3::new(0.005, 0.007, 0.012);
    let mut c = mix(top, bot, v.powf(0.65));

    // sparse star field
    let sx = (x / 3) as i32;
    let sy = (y / 3) as i32;
    let n = hash2(sx, sy);
    if n > 0.9965 {
        let b = ((n - 0.9965) / 0.0035).clamp(0.0, 1.0);
        c = c + V3::new(0.45, 0.50, 0.60) * b * 0.5;
    }

    c
}

fn project(cam: Camera, p: V3, w: usize, h: usize) -> Option<(f32, f32, f32)> {
    let fwd = (cam.target - cam.eye).norm();
    let right = fwd.cross(cam.up).norm();
    let up = right.cross(fwd).norm();
    let rel = p - cam.eye;

    let z = rel.dot(fwd);
    if z <= 0.02 {
        return None;
    }

    let tan_half = (cam.fov_deg.to_radians() * 0.5).tan();
    let x_ndc = rel.dot(right) / (z * tan_half);
    let y_ndc = rel.dot(up) / (z * tan_half);

    let aspect = w as f32 / h as f32;
    let sx = (x_ndc / aspect * 0.5 + 0.5) * w as f32;
    let sy = (-y_ndc * 0.5 + 0.5) * h as f32;

    if sx < -8.0 || sy < -8.0 || sx > w as f32 + 8.0 || sy > h as f32 + 8.0 {
        return None;
    }

    Some((sx, sy, z))
}

fn draw_splat(
    color_buf: &mut [V3],
    z_buf: &mut [f32],
    w: usize,
    h: usize,
    sx: f32,
    sy: f32,
    z: f32,
    color: V3,
    radius_px: i32,
) {
    let cx = sx as i32;
    let cy = sy as i32;
    for dy in -radius_px..=radius_px {
        for dx in -radius_px..=radius_px {
            let x = cx + dx;
            let y = cy + dy;
            if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 {
                continue;
            }
            let d2 = (dx * dx + dy * dy) as f32;
            let r2 = (radius_px * radius_px) as f32;
            if d2 > r2 {
                continue;
            }
            let idx = y as usize * w + x as usize;
            if z > z_buf[idx] + 0.003 {
                continue;
            }
            let a = (1.0 - (d2 / r2).sqrt()).powf(0.9) * 0.95;
            if z < z_buf[idx] {
                z_buf[idx] = z;
            }
            color_buf[idx] = color_buf[idx] * (1.0 - a) + color * a;
        }
    }
}

fn shade_metal(base: V3, normal: V3, view_dir: V3, line_mask: f32) -> V3 {
    let n = normal.norm();
    let v = view_dir.norm();

    let lights = [
        (V3::new(0.70, 0.82, 0.42).norm(), V3::new(1.2, 1.15, 1.05)),
        (V3::new(-0.45, 0.32, -0.85).norm(), V3::new(0.55, 0.62, 0.75)),
    ];

    let mut out = base * (0.05 + 0.12 * (1.0 - line_mask));
    for (l, c) in lights {
        let ndotl = n.dot(l).clamp(0.0, 1.0);
        let h = (v + l).norm();
        let ndoth = n.dot(h).clamp(0.0, 1.0);
        let spec = ndoth.powf(45.0 + 20.0 * line_mask);
        out = out + base * ndotl * c + V3::new(1.0, 1.0, 1.0) * spec * c * 0.85;
    }
    let rim = (1.0 - n.dot(v).clamp(0.0, 1.0)).powf(2.1);
    out + V3::new(0.18, 0.22, 0.30) * rim * 0.30
}

fn torus_shell_line_mask(u: f32, v: f32, mesh_q: f32) -> f32 {
    // Hex-like tri-directional mesh in toroidal parameter space.
    let p = V3::new(u / (2.0 * PI) * 28.0, v / (2.0 * PI) * 18.0, 0.0);
    let g1 = (PI * p.x).sin().abs();
    let g2 = (PI * (0.5 * p.x + 0.866 * p.y)).sin().abs();
    let g3 = (PI * (0.5 * p.x - 0.866 * p.y)).sin().abs();
    let g = g1.min(g2).min(g3);
    let thickness = 0.09 - 0.03 * mesh_q.clamp(0.0, 1.0);
    1.0 - smoothstep(thickness, thickness + 0.10, g)
}

fn render_scene(params: &ReactorParams, cam: Camera, w: usize, h: usize) -> RgbImage {
    let mut color_buf = vec![V3::new(0.0, 0.0, 0.0); w * h];
    let mut z_buf = vec![1.0e9_f32; w * h];

    for y in 0..h {
        for x in 0..w {
            color_buf[y * w + x] = background(x, y, w, h);
        }
    }

    let base = material_color(&params.material);
    let r = params.r_major_m;
    let a = params.a_minor_m;

    // Outer vessel torus
    let shell_r = a * 1.22;
    let u_steps = if w >= 3000 { 540 } else { 360 };
    let v_steps = if w >= 3000 { 260 } else { 180 };
    for iu in 0..u_steps {
        let u = 2.0 * PI * (iu as f32 / u_steps as f32);
        let cu = u.cos();
        let su = u.sin();
        for iv in 0..v_steps {
            let v = 2.0 * PI * (iv as f32 / v_steps as f32);
            let cv = v.cos();
            let sv = v.sin();

            let rr = r + shell_r * cv;
            let p = V3::new(rr * cu, shell_r * sv, rr * su);
            let n = V3::new(cv * cu, sv, cv * su).norm();

            if let Some((sx, sy, z)) = project(cam, p, w, h) {
                let line = torus_shell_line_mask(u, v, params.mesh_quality);
                let albedo = mix(base * 0.55, base * 0.95, 1.0 - line * 0.65);
                let shaded = shade_metal(albedo, n, cam.eye - p, line).clamp01();
                draw_splat(
                    &mut color_buf,
                    &mut z_buf,
                    w,
                    h,
                    sx,
                    sy,
                    z,
                    shaded,
                    if w >= 3000 { 2 } else { 1 },
                );
            }
        }
    }

    // Inner coil shell
    let coil_r = a * 0.92;
    let coil_col = mix(base, V3::new(0.30, 0.55, 0.95), 0.14);
    for iu in 0..(u_steps / 2) {
        let u = 2.0 * PI * (iu as f32 / (u_steps / 2) as f32);
        let cu = u.cos();
        let su = u.sin();
        for iv in 0..(v_steps / 2) {
            let v = 2.0 * PI * (iv as f32 / (v_steps / 2) as f32);
            let cv = v.cos();
            let sv = v.sin();
            let rr = r + coil_r * cv;
            let p = V3::new(rr * cu, coil_r * sv, rr * su);
            let n = V3::new(cv * cu, sv, cv * su).norm();

            if let Some((sx, sy, z)) = project(cam, p, w, h) {
                let shade = shade_metal(coil_col, n, cam.eye - p, 0.2).clamp01();
                draw_splat(&mut color_buf, &mut z_buf, w, h, sx, sy, z, shade, 1);
            }
        }
    }

    // Plasma torus (emissive core)
    let plasma_r = a * 0.46;
    let plasma_intensity = (0.45 + 0.55 * (params.t_kev / 180.0).clamp(0.0, 1.5))
        * (0.60 + 0.40 * (params.q_engineering / 2.0).clamp(0.0, 1.0));
    for iu in 0..(u_steps / 2) {
        let u = 2.0 * PI * (iu as f32 / (u_steps / 2) as f32);
        let cu = u.cos();
        let su = u.sin();
        for iv in 0..(v_steps / 2) {
            let v = 2.0 * PI * (iv as f32 / (v_steps / 2) as f32);
            let cv = v.cos();
            let sv = v.sin();
            let rr = r + plasma_r * cv;
            let p = V3::new(rr * cu, plasma_r * sv, rr * su);
            let n = V3::new(cv * cu, sv, cv * su).norm();
            if let Some((sx, sy, z)) = project(cam, p, w, h) {
                let ndotv = n.dot((cam.eye - p).norm()).abs();
                let emissive = mix(
                    V3::new(0.95, 0.48, 0.16),
                    V3::new(1.00, 0.76, 0.24),
                    ndotv,
                ) * (0.45 + 0.75 * plasma_intensity);
                draw_splat(
                    &mut color_buf,
                    &mut z_buf,
                    w,
                    h,
                    sx,
                    sy,
                    z + 0.0005,
                    emissive,
                    if w >= 3000 { 2 } else { 1 },
                );
            }
        }
    }

    // Final post (tone-map + mild bloom from bright plasma)
    let mut img = RgbImage::new(w as u32, h as u32);
    for y in 0..h {
        for x in 0..w {
            let i = y * w + x;
            let c = color_buf[i].clamp01();
            let mapped = V3::new(
                c.x / (1.0 + 0.65 * c.x),
                c.y / (1.0 + 0.65 * c.y),
                c.z / (1.0 + 0.65 * c.z),
            );
            let srgb = V3::new(
                mapped.x.powf(1.0 / 2.2),
                mapped.y.powf(1.0 / 2.2),
                mapped.z.powf(1.0 / 2.2),
            )
            .clamp01();
            img.put_pixel(
                x as u32,
                y as u32,
                Rgb([
                    (srgb.x * 255.0 + 0.5) as u8,
                    (srgb.y * 255.0 + 0.5) as u8,
                    (srgb.z * 255.0 + 0.5) as u8,
                ]),
            );
        }
    }

    img
}

fn main() -> Result<(), String> {
    let out_dir = PathBuf::from(
        env::var("GUTOE_RTSC_FUSION_RENDER_OUT")
            .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_honeycomb_fusion_render".to_string()),
    );
    fs::create_dir_all(&out_dir)
        .map_err(|e| format!("failed creating output dir {}: {e}", out_dir.display()))?;

    let reactor_json = PathBuf::from(env::var("GUTOE_FUSION_REACTOR_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/rtsc_honeycomb_fusion_reactor/rtsc_honeycomb_fusion_reactor.json"
            .to_string()
    }));
    let frames_dir = out_dir.join("frames");
    fs::create_dir_all(&frames_dir)
        .map_err(|e| format!("failed creating frames dir {}: {e}", frames_dir.display()))?;

    let params = load_best(&reactor_json)?;

    // High-res still
    let still_cam = Camera {
        eye: V3::new(params.r_major_m * 2.7, params.a_minor_m * 1.4, params.r_major_m * 2.1),
        target: V3::new(0.0, 0.0, 0.0),
        up: V3::new(0.0, 1.0, 0.0),
        fov_deg: 38.0,
    };
    let still = render_scene(&params, still_cam, 3840, 2160);
    let png_path = out_dir.join("rtsc_honeycomb_fusion_physical_4k.png");
    still
        .save(&png_path)
        .map_err(|e| format!("failed writing PNG {}: {e}", png_path.display()))?;

    // Rotating sequence (used by ffmpeg to make MP4)
    let width = 1920usize;
    let height = 1080usize;
    let n_frames = env::var("GUTOE_RTSC_ROT_FRAMES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(240);
    let orbit_radius = params.r_major_m * 3.5;

    for i in 0..n_frames {
        let t = i as f32 / n_frames as f32;
        let ang = 2.0 * PI * t;
        let cam = Camera {
            eye: V3::new(orbit_radius * ang.cos(), params.a_minor_m * 1.2, orbit_radius * ang.sin()),
            target: V3::new(0.0, 0.0, 0.0),
            up: V3::new(0.0, 1.0, 0.0),
            fov_deg: 42.0,
        };
        let frame = render_scene(&params, cam, width, height);
        let frame_path = frames_dir.join(format!("frame_{:04}.png", i));
        frame
            .save(&frame_path)
            .map_err(|e| format!("failed writing frame {}: {e}", frame_path.display()))?;
    }

    let meta = format!(
        concat!(
            "{{\n",
            "  \"mode\": \"physical_reactor_render\",\n",
            "  \"reactor_source\": \"{}\",\n",
            "  \"best_material\": \"{}\",\n",
            "  \"geometry\": {{\"r_major_m\": {:.6}, \"a_minor_m\": {:.6}, \"mesh_quality\": {:.3}}},\n",
            "  \"operating\": {{\"b_operating_t\": {:.3}, \"t_kev\": {:.3}, \"p_net_w\": {:.6e}, \"q_engineering\": {:.6}}},\n",
            "  \"png\": \"{}\",\n",
            "  \"frames_dir\": \"{}\",\n",
            "  \"frames\": {}\n",
            "}}\n"
        ),
        reactor_json.display(),
        params.material,
        params.r_major_m,
        params.a_minor_m,
        params.mesh_quality,
        params.b_operating_t,
        params.t_kev,
        params.p_net_w,
        params.q_engineering,
        png_path.display(),
        frames_dir.display(),
        n_frames
    );

    let meta_path = out_dir.join("rtsc_honeycomb_fusion_render.json");
    fs::write(&meta_path, meta)
        .map_err(|e| format!("failed writing metadata {}: {e}", meta_path.display()))?;

    println!("wrote {}", png_path.display());
    println!("wrote {}", meta_path.display());
    println!("wrote {} frames to {}", n_frames, frames_dir.display());

    Ok(())
}
