use gutoe_physics::{scan_nuclear_chart, NucleusRecord, ScanConfig};
use image::{Rgb, RgbImage};
use std::collections::BTreeMap;
use std::env;
use std::f32::consts::PI;
use std::fs;
use std::path::PathBuf;

const ELECTRON_MASS_MEV_OBS: f64 = 0.510_998_950;
const PROTON_MASS_MEV_OBS: f64 = 938.272_088_16;
const NEUTRON_MASS_MEV_OBS: f64 = 939.565_420_52;
const BETA_MASS_COEFF_Z_MEV: f64 =
    (PROTON_MASS_MEV_OBS + ELECTRON_MASS_MEV_OBS) - NEUTRON_MASS_MEV_OBS;
const MN_MINUS_MP_MINUS_ME_MEV: f64 = 0.782_333;
const MN_MINUS_MP_MEV: f64 = 1.293_332;
const ODD_A_PAIR_RELAX_COEFF: f64 = 1.0 / 12.0;
const ODD_Z_GAP_WEAK_MARGIN_MEV: f64 = 0.85;

#[derive(Clone, Copy, Debug, Default)]
struct BetaDecayQ {
    q_beta_minus_mev: Option<f64>,
    q_ec_mev: Option<f64>,
}

#[derive(Clone, Copy, Debug, Default)]
struct BetaLocalState {
    is_local_min: bool,
    delta_to_isobar_min_mev: f64,
}

fn classify_long_lived(r: &NucleusRecord, beta_local: BetaLocalState, beta_q: BetaDecayQ) -> bool {
    let z50_corridor = r.z == 50 && (108..=126).contains(&r.a);
    let beta_q_rescue = z50_corridor && r.beta_optimal_for_a;
    let quasi_stable_even_even = z50_corridor
        && r.z % 2 == 0
        && r.n % 2 == 0
        && beta_q.q_beta_minus_mev.map(|q| q < 1.5).unwrap_or(false)
        && beta_q.q_ec_mev.map(|q| q < 0.0).unwrap_or(false);
    let odd_a_pairing_relax = (r.a % 2 == 1)
        && (beta_local.delta_to_isobar_min_mev
            <= ODD_A_PAIR_RELAX_COEFF * (12.0 / (r.a as f64).sqrt()));
    let beta_ok =
        beta_local.is_local_min || beta_q_rescue || quasi_stable_even_even || odd_a_pairing_relax;
    let weak_q_margin_mev = {
        let mut m = f64::INFINITY;
        if let Some(q) = beta_q.q_beta_minus_mev {
            m = m.min(-q);
        }
        if let Some(q) = beta_q.q_ec_mev {
            m = m.min(-q);
        }
        m
    };
    let tc_pm_weak_gap = (r.z == 43 || r.z == 61)
        && r.n >= 46
        && weak_q_margin_mev.is_finite()
        && weak_q_margin_mev < ODD_Z_GAP_WEAK_MARGIN_MEV;
    let fail_beta_optimal = !beta_ok || tc_pm_weak_gap;
    let fail_fissility = r.fissility > 1.0;
    let fail_s2n = if r.n <= 2 {
        false
    } else {
        !r.s2n_mev.map(|v| v > 0.0).unwrap_or(false)
    };
    let fail_s2p = if r.z <= 2 {
        false
    } else {
        !r.s2p_mev.map(|v| v > 0.0).unwrap_or(false)
    };
    let fail_sf = r.z > 82 && r.sf_log10_half_life_s < 20.0;
    !(fail_beta_optimal || fail_fissility || fail_s2n || fail_s2p || fail_sf)
}

fn build_beta_local_state_map(records: &[NucleusRecord]) -> BTreeMap<(u16, u16), BetaLocalState> {
    let mut mass_proxy_by_az: BTreeMap<(u16, u16), f64> = BTreeMap::new();
    for r in records {
        let mass_proxy = BETA_MASS_COEFF_Z_MEV * r.z as f64 - r.binding_mev;
        mass_proxy_by_az.insert((r.a, r.z), mass_proxy);
    }
    let mut min_proxy_by_a: BTreeMap<u16, f64> = BTreeMap::new();
    for (&(a, _), &m) in &mass_proxy_by_az {
        min_proxy_by_a
            .entry(a)
            .and_modify(|cur| {
                if m < *cur {
                    *cur = m;
                }
            })
            .or_insert(m);
    }
    let mut out = BTreeMap::new();
    for r in records {
        let Some(&m0) = mass_proxy_by_az.get(&(r.a, r.z)) else {
            out.insert((r.z, r.n), BetaLocalState::default());
            continue;
        };
        let left_ok = if r.z > 1 {
            mass_proxy_by_az
                .get(&(r.a, r.z - 1))
                .map(|&ml| m0 <= ml + 1e-9)
                .unwrap_or(true)
        } else {
            true
        };
        let right_ok = mass_proxy_by_az
            .get(&(r.a, r.z + 1))
            .map(|&mr| m0 <= mr + 1e-9)
            .unwrap_or(true);
        let min_proxy = min_proxy_by_a.get(&r.a).copied().unwrap_or(m0);
        out.insert(
            (r.z, r.n),
            BetaLocalState {
                is_local_min: left_ok && right_ok,
                delta_to_isobar_min_mev: (m0 - min_proxy).max(0.0),
            },
        );
    }
    out
}

fn build_beta_q_map(records: &[NucleusRecord]) -> BTreeMap<(u16, u16), BetaDecayQ> {
    let mut binding_by_zn: BTreeMap<(u16, u16), f64> = BTreeMap::new();
    for r in records {
        binding_by_zn.insert((r.z, r.n), r.binding_mev);
    }
    let mut out: BTreeMap<(u16, u16), BetaDecayQ> = BTreeMap::new();
    for r in records {
        let q_beta_minus_mev = if r.n > 0 {
            binding_by_zn
                .get(&(r.z + 1, r.n - 1))
                .map(|&b_d| (b_d - r.binding_mev) + MN_MINUS_MP_MINUS_ME_MEV)
        } else {
            None
        };
        let q_ec_mev = binding_by_zn
            .get(&(r.z.saturating_sub(1), r.n + 1))
            .map(|&b_d| (b_d - r.binding_mev) - MN_MINUS_MP_MEV);
        out.insert(
            (r.z, r.n),
            BetaDecayQ {
                q_beta_minus_mev,
                q_ec_mev,
            },
        );
    }
    out
}

fn env_u32(name: &str, default: u32) -> u32 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .unwrap_or(default)
}

#[derive(Clone, Copy, Debug, Default)]
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
    fn abs(self) -> Self {
        Self::new(self.x.abs(), self.y.abs(), self.z.abs())
    }
    fn max(self, o: Self) -> Self {
        Self::new(self.x.max(o.x), self.y.max(o.y), self.z.max(o.z))
    }
    fn clamp01(self) -> Self {
        Self::new(
            self.x.clamp(0.0, 1.0),
            self.y.clamp(0.0, 1.0),
            self.z.clamp(0.0, 1.0),
        )
    }
}

use std::ops::{Add, Mul, Sub};
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

fn mix(a: V3, b: V3, t: f32) -> V3 {
    a * (1.0 - t) + b * t
}

#[derive(Clone, Copy, Debug)]
enum Crystal {
    Bcc,
    Fcc,
    Hcp,
    Diamond,
    Molecular,
}

#[derive(Clone, Copy, Debug)]
struct Material {
    z: u16,
    symbol: &'static str,
    stable_like: usize,
    crystal: Crystal,
    family: &'static str,
    albedo: V3,
    metallic: f32,
    roughness: f32,
    ior: f32,
}

fn symbol_of(z: u16) -> &'static str {
    const S: [&str; 95] = [
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P",
        "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
        "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re",
        "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu",
    ];
    if (z as usize) < S.len() {
        S[z as usize]
    } else {
        "E?"
    }
}

fn family_of_z(z: u16) -> &'static str {
    match z {
        3 | 11 | 19 | 37 | 55 | 87 => "alkali",
        4 | 12 | 20 | 38 | 56 | 88 => "alkaline",
        9 | 17 | 35 | 53 | 85 => "halogen",
        2 | 10 | 18 | 36 | 54 | 86 => "noble",
        1 | 6 | 7 | 8 | 15 | 16 | 34 => "nonmetal",
        5 | 14 | 32 | 33 | 51 | 52 => "metalloid",
        57..=71 => "lanthanide",
        89..=94 => "actinide",
        13 | 31 | 49 | 50 | 81 | 82 | 83 => "post-transition",
        _ => "transition",
    }
}

fn crystal_proxy(family: &str, z: u16, stable_like: usize) -> Crystal {
    match family {
        "alkali" | "transition" => match z % 3 {
            0 => Crystal::Bcc,
            1 => Crystal::Fcc,
            _ => Crystal::Hcp,
        },
        "alkaline" | "post-transition" => match stable_like % 3 {
            0 => Crystal::Hcp,
            1 => Crystal::Fcc,
            _ => Crystal::Bcc,
        },
        "metalloid" => Crystal::Diamond,
        "nonmetal" | "halogen" | "noble" => Crystal::Molecular,
        "lanthanide" | "actinide" => Crystal::Hcp,
        _ => Crystal::Bcc,
    }
}

fn material_params(z: u16, family: &str, stable_like: usize) -> (V3, f32, f32, f32) {
    let mut albedo = match z {
        6 => V3::new(0.35, 0.37, 0.40),
        13 => V3::new(0.78, 0.81, 0.84),
        26 => V3::new(0.56, 0.58, 0.60),
        28 => V3::new(0.64, 0.67, 0.70),
        29 => V3::new(0.78, 0.49, 0.32),
        47 => V3::new(0.83, 0.87, 0.90),
        74 => V3::new(0.54, 0.56, 0.60),
        78 => V3::new(0.76, 0.79, 0.83),
        79 => V3::new(0.82, 0.70, 0.25),
        82 => V3::new(0.52, 0.54, 0.58),
        _ => match family {
            "alkali" => V3::new(0.83, 0.80, 0.75),
            "alkaline" => V3::new(0.84, 0.82, 0.78),
            "transition" => V3::new(0.68, 0.72, 0.77),
            "post-transition" => V3::new(0.72, 0.70, 0.67),
            "metalloid" => V3::new(0.50, 0.58, 0.66),
            "nonmetal" => V3::new(0.38, 0.47, 0.57),
            "halogen" => V3::new(0.47, 0.54, 0.66),
            "noble" => V3::new(0.58, 0.65, 0.76),
            "lanthanide" => V3::new(0.74, 0.76, 0.79),
            "actinide" => V3::new(0.71, 0.67, 0.60),
            _ => V3::new(0.68, 0.72, 0.77),
        },
    };
    let drift = 0.02 * ((z as f32 * 0.37).sin());
    albedo = (albedo + V3::new(drift, drift, drift)).clamp01();
    if family == "transition"
        || family == "post-transition"
        || family == "lanthanide"
        || family == "actinide"
        || family == "alkali"
        || family == "alkaline"
    {
        (albedo, 0.93, 0.12 + 0.045 * ((stable_like % 5) as f32), 2.0)
    } else if family == "metalloid" {
        (albedo, 0.20, 0.17 + 0.03 * ((z % 4) as f32), 2.6)
    } else {
        (
            albedo,
            0.04,
            0.09 + 0.03 * (((z as usize + stable_like) % 5) as f32),
            1.45 + 0.08 * (((z / 10) % 4) as f32),
        )
    }
}

fn rotate_x(p: V3, a: f32) -> V3 {
    let c = a.cos();
    let s = a.sin();
    V3::new(p.x, c * p.y - s * p.z, s * p.y + c * p.z)
}

fn rotate_y(p: V3, a: f32) -> V3 {
    let c = a.cos();
    let s = a.sin();
    V3::new(c * p.x + s * p.z, p.y, -s * p.x + c * p.z)
}

fn sdf_sphere(p: V3, r: f32) -> f32 {
    p.len() - r
}

fn sdf_box(p: V3, b: V3) -> f32 {
    let q = p.abs() - b;
    let outside = q.max(V3::new(0.0, 0.0, 0.0)).len();
    let inside = q.x.max(q.y.max(q.z)).min(0.0);
    outside + inside
}

fn sdf_octahedron(p: V3, s: f32) -> f32 {
    (p.x.abs() + p.y.abs() + p.z.abs() - s) * 0.57735026
}

fn sdf_hex_prism(p: V3, h: V3) -> f32 {
    let px = p.x.abs();
    let py = p.y.abs();
    let pz = p.z.abs();
    let d1 = (px * -0.8660254 + py * 0.5).max(py) - h.x;
    let d2 = pz - h.y;
    d1.max(d2)
}

fn sdf_crystal(p: V3, crystal: Crystal) -> f32 {
    match crystal {
        Crystal::Bcc => sdf_octahedron(p, 0.55).max(sdf_box(p, V3::new(0.44, 0.44, 0.44))),
        Crystal::Fcc => sdf_box(p, V3::new(0.48, 0.38, 0.48)).max(sdf_octahedron(p, 0.67)),
        Crystal::Hcp => sdf_hex_prism(p, V3::new(0.36, 0.40, 0.0)),
        Crystal::Diamond => {
            let p2 = rotate_y(p, 35.0_f32.to_radians());
            sdf_octahedron(p2, 0.48).max(sdf_box(p2, V3::new(0.34, 0.50, 0.34)))
        }
        Crystal::Molecular => {
            let centers = [
                V3::new(-0.18, 0.02, -0.12),
                V3::new(0.18, 0.02, -0.12),
                V3::new(-0.07, 0.16, 0.17),
                V3::new(0.10, -0.15, 0.13),
            ];
            let mut d = 1e9_f32;
            for c in centers {
                d = d.min(sdf_sphere(p - c, 0.20));
            }
            d
        }
    }
}

fn scene_sdf(p: V3, crystal: Crystal, ax: f32, ay: f32) -> (f32, u8) {
    let q = rotate_y(rotate_x(p, ax), ay);
    let d_obj = sdf_crystal(q, crystal);
    let d_floor = p.y + 0.54;
    if d_obj < d_floor {
        (d_obj, 1)
    } else {
        (d_floor, 2)
    }
}

fn estimate_normal(p: V3, crystal: Crystal, ax: f32, ay: f32) -> V3 {
    let e = 1.2e-3;
    let dx = scene_sdf(V3::new(p.x + e, p.y, p.z), crystal, ax, ay).0
        - scene_sdf(V3::new(p.x - e, p.y, p.z), crystal, ax, ay).0;
    let dy = scene_sdf(V3::new(p.x, p.y + e, p.z), crystal, ax, ay).0
        - scene_sdf(V3::new(p.x, p.y - e, p.z), crystal, ax, ay).0;
    let dz = scene_sdf(V3::new(p.x, p.y, p.z + e), crystal, ax, ay).0
        - scene_sdf(V3::new(p.x, p.y, p.z - e), crystal, ax, ay).0;
    V3::new(dx, dy, dz).norm()
}

fn env_color(rd: V3) -> V3 {
    let t = (0.5 * (rd.y + 1.0)).clamp(0.0, 1.0);
    let sky = V3::new(0.12, 0.17, 0.26);
    let horizon = V3::new(0.28, 0.26, 0.30);
    mix(horizon, sky, t)
}

fn fresnel_schlick(cos_theta: f32, f0: V3) -> V3 {
    f0 + (V3::new(1.0, 1.0, 1.0) - f0) * (1.0 - cos_theta).powf(5.0)
}

fn raytrace_pixel(material: Material, u: f32, v: f32, max_steps: u32) -> V3 {
    let ro = V3::new(0.0, 0.07, 2.35);
    let rd = V3::new(u * 0.95, -v * 0.95, -1.0).norm();
    let mut t = 0.0f32;
    let ax = (11.0 + (material.z % 13) as f32).to_radians();
    let ay = (22.0 + (material.z % 29) as f32).to_radians();
    let mut mat_id = 0u8;
    let mut hit = false;
    for _ in 0..max_steps {
        let p = ro + rd * t;
        let (d, m) = scene_sdf(p, material.crystal, ax, ay);
        if d < 9e-4 {
            mat_id = m;
            hit = true;
            break;
        }
        t += d.clamp(1e-4, 0.5);
        if t > 8.0 {
            break;
        }
    }
    if !hit {
        return env_color(rd);
    }

    let p = ro + rd * t;
    let n = estimate_normal(p, material.crystal, ax, ay);
    let vdir = (rd * -1.0).norm();
    let l = V3::new(0.42, 0.86, 0.33).norm();
    let h = (vdir + l).norm();

    let ndotl = n.dot(l).clamp(0.0, 1.0);
    let ndotv = n.dot(vdir).clamp(0.0, 1.0);
    let ndoth = n.dot(h).clamp(0.0, 1.0);
    let vdoth = vdir.dot(h).clamp(0.0, 1.0);

    if mat_id == 2 {
        let checker = (((((p.x + 4.0) * 4.0).floor() + (((p.z + 4.0) * 4.0).floor())) as i32) & 1)
            as f32;
        let base = 0.15 + checker * 0.06;
        let c = V3::new(base, base, base + 0.01);
        return (c * 0.9 + c * ndotl * 0.6).clamp01();
    }

    let metallic = material.metallic;
    let roughness = material.roughness;
    let a = (roughness * roughness).max(0.035);
    let a2 = a * a;
    let denom = ndoth * ndoth * (a2 - 1.0) + 1.0;
    let d = a2 / (PI * denom * denom).max(1e-7);
    let k = (roughness + 1.0).powi(2) / 8.0;
    let gv = ndotv / (ndotv * (1.0 - k) + k).max(1e-7);
    let gl = ndotl / (ndotl * (1.0 - k) + k).max(1e-7);
    let g = gv * gl;
    let f0_dielectric = ((material.ior - 1.0) / (material.ior + 1.0)).powi(2);
    let f0 = mix(
        V3::new(f0_dielectric, f0_dielectric, f0_dielectric),
        material.albedo,
        metallic,
    );
    let f = fresnel_schlick(vdoth, f0);
    let spec = f * (d * g / (4.0 * ndotv * ndotl).max(1e-6));
    let diffuse = material.albedo * ((1.0 - metallic) / PI);
    let refl = (vdir - n * (2.0 * vdir.dot(n))).norm();
    let env_refl = env_color(refl);
    let ambient = material.albedo * 0.18 + env_refl * (0.22 * (1.0 - roughness));
    let direct = (diffuse + spec) * ndotl * V3::new(1.5, 1.45, 1.35);
    let mut col = (ambient + direct).clamp01();
    let fog = (-0.06 * (p - ro).len()).exp();
    col = col * fog + env_color(rd) * (1.0 - fog);
    col
}

fn linear_to_srgb_u8(c: V3) -> [u8; 3] {
    let g = V3::new(
        c.x.clamp(0.0, 1.0).powf(1.0 / 2.2),
        c.y.clamp(0.0, 1.0).powf(1.0 / 2.2),
        c.z.clamp(0.0, 1.0).powf(1.0 / 2.2),
    );
    [
        (g.x * 255.0 + 0.5) as u8,
        (g.y * 255.0 + 0.5) as u8,
        (g.z * 255.0 + 0.5) as u8,
    ]
}

fn materials_from_physics() -> Vec<Material> {
    let records = scan_nuclear_chart(ScanConfig::default());
    let beta_local = build_beta_local_state_map(&records);
    let beta_q = build_beta_q_map(&records);
    let mut counts: BTreeMap<u16, usize> = BTreeMap::new();
    for r in &records {
        let bl = beta_local.get(&(r.z, r.n)).copied().unwrap_or_default();
        let bq = beta_q.get(&(r.z, r.n)).copied().unwrap_or_default();
        if classify_long_lived(r, bl, bq) {
            *counts.entry(r.z).or_insert(0) += 1;
        }
    }
    let mut mats = Vec::new();
    for (&z, &stable_like) in &counts {
        if !(1..=94).contains(&z) {
            continue;
        }
        let family = family_of_z(z);
        let crystal = crystal_proxy(family, z, stable_like);
        let (albedo, metallic, roughness, ior) = material_params(z, family, stable_like);
        mats.push(Material {
            z,
            symbol: symbol_of(z),
            stable_like,
            crystal,
            family,
            albedo,
            metallic,
            roughness,
            ior,
        });
    }
    mats.sort_by_key(|m| m.z);
    mats
}

fn main() {
    let mut out_dir = PathBuf::from(
        env::var("GUTOE_MATERIAL_RAYMARCH_OUT")
            .unwrap_or_else(|_| "/tmp/bh_renders/material_atlas_rust".to_string()),
    );
    let tile = env_u32("GUTOE_MATERIAL_TILE", 120).max(80) as usize;
    let cols = env_u32("GUTOE_MATERIAL_COLS", 10).max(1) as usize;
    let max_steps = env_u32("GUTOE_MATERIAL_STEPS", 64).max(24);
    if let Some(arg) = env::args().skip(1).next() {
        if arg == "--help" || arg == "-h" {
            println!(
                "Usage: material_raymarcher_physics [OUT_DIR]\n\
                 Env overrides:\n\
                   GUTOE_MATERIAL_RAYMARCH_OUT (default /tmp/bh_renders/material_atlas_rust)\n\
                   GUTOE_MATERIAL_TILE         (default 120)\n\
                   GUTOE_MATERIAL_COLS         (default 10)\n\
                   GUTOE_MATERIAL_STEPS        (default 64)\n"
            );
            return;
        }
        out_dir = PathBuf::from(arg);
    }
    fs::create_dir_all(&out_dir).expect("create out dir");

    let mats = materials_from_physics();
    let rows = mats.len().div_ceil(cols);
    let pad = 4usize;
    let header_h = 28usize;
    let width = cols * tile + (cols + 1) * pad;
    let height = rows * tile + (rows + 1) * pad + header_h;
    let mut img = RgbImage::new(width as u32, height as u32);
    for y in 0..height {
        for x in 0..width {
            img.put_pixel(x as u32, y as u32, Rgb([9, 11, 14]));
        }
    }

    for (i, m) in mats.iter().enumerate() {
        let r = i / cols;
        let c = i % cols;
        let x0 = pad + c * (tile + pad);
        let y0 = header_h + pad + r * (tile + pad);
        for py in 0..tile {
            for px in 0..tile {
                let u = ((px as f32 + 0.5) / tile as f32) * 2.0 - 1.0;
                let v = ((py as f32 + 0.5) / tile as f32) * 2.0 - 1.0;
                let col = linear_to_srgb_u8(raytrace_pixel(*m, u, v, max_steps));
                img.put_pixel((x0 + px) as u32, (y0 + py) as u32, Rgb(col));
            }
        }
        // thin border
        for bx in 0..tile {
            img.put_pixel((x0 + bx) as u32, y0 as u32, Rgb([44, 50, 60]));
            img.put_pixel((x0 + bx) as u32, (y0 + tile - 1) as u32, Rgb([44, 50, 60]));
        }
        for by in 0..tile {
            img.put_pixel(x0 as u32, (y0 + by) as u32, Rgb([44, 50, 60]));
            img.put_pixel((x0 + tile - 1) as u32, (y0 + by) as u32, Rgb([44, 50, 60]));
        }
    }

    let png_path = out_dir.join("material_raymarch_atlas_physics.png");
    img.save(&png_path).expect("save png");

    let mut json = String::new();
    json.push_str("{\n");
    json.push_str(&format!("  \"elements_rendered\": {},\n", mats.len()));
    json.push_str(&format!("  \"tile_size\": {},\n", tile));
    json.push_str(&format!("  \"cols\": {},\n", cols));
    json.push_str(&format!("  \"rows\": {},\n", rows));
    json.push_str(&format!("  \"raymarch_steps\": {},\n", max_steps));
    json.push_str("  \"physics_source\": \"scan_nuclear_chart + classify_long_lived (mass_periodic lane)\",\n");
    json.push_str("  \"theorem_refs\": [\n");
    json.push_str("    \"Gutoe.NuclearFirstPrinciples.nuclear_structural_bundle\",\n");
    json.push_str("    \"Gutoe.DarkMatterSector.visible_dark_state_count_split\",\n");
    json.push_str("    \"Gutoe.FineStructure.alpha_inverse_d4\"\n");
    json.push_str("  ],\n");
    json.push_str("  \"elements\": [\n");
    for (i, m) in mats.iter().enumerate() {
        let cstr = match m.crystal {
            Crystal::Bcc => "bcc",
            Crystal::Fcc => "fcc",
            Crystal::Hcp => "hcp",
            Crystal::Diamond => "diamond",
            Crystal::Molecular => "molecular",
        };
        json.push_str(&format!(
            "    {{\"z\":{},\"symbol\":\"{}\",\"stable_like\":{},\"family\":\"{}\",\"crystal\":\"{}\",\"metallic\":{:.4},\"roughness\":{:.4},\"ior\":{:.4}}}{}\n",
            m.z,
            m.symbol,
            m.stable_like,
            m.family,
            cstr,
            m.metallic,
            m.roughness,
            m.ior,
            if i + 1 == mats.len() { "" } else { "," }
        ));
    }
    json.push_str("  ]\n}\n");
    let json_path = out_dir.join("material_raymarch_atlas_physics.json");
    fs::write(&json_path, json).expect("write json");

    println!("wrote {}", png_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "rendered {} elements from physics lane; atlas {}x{}",
        mats.len(),
        width,
        height
    );
}
