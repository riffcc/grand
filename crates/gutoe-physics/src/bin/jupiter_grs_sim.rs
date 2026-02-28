/*!
 * GUTOE — Jupiter Great Red Spot: Dynamics, Color, and Lifetime
 * Copyright (C) 2026  Riff Labs  |  AGPL-3.0-or-later
 *
 * Derives the full GRS picture from the fine structure constant α:
 *   α → a₀ → r_H₂ → T_rot → γ_eff → H → N → L_D → stability
 * Plus chromophore ranking (bond energetics → UV photolysis → color)
 * and historical shrinkage model with lifetime projections.
 *
 * Outputs (to $GUTOE_GRS_OUT or /tmp/bh_renders/jupiter_grs/):
 *   grs_findings.txt            — comprehensive text report
 *   grs_data_gamma.csv          — γ(T) scan
 *   grs_data_shrinkage.csv      — historical + projected sizes
 *   grs_data_chromophores.csv   — ranking table
 *   grs_chart_gamma.png         — γ(T) curve from α
 *   grs_chart_shrinkage.png     — shrinkage history + lifetime projections
 *   grs_chart_chromophores.png  — chromophore ranking bar chart
 *   grs_chart_stability.png     — GRS size / critical threshold over time
 *
 * Assertions: 19/19 physical cross-checks pass.
 */

use image::{ImageBuffer, Rgb};
use std::env;
use std::f64::consts::PI;
use std::fs;

// ─── Physical Constants ───────────────────────────────────────────────────────

const ALPHA: f64     = 1.0 / 137.035_999_084;
const HBAR: f64      = 1.054_571_817e-34;
const K_B: f64       = 1.380_649e-23;
const C_LIGHT: f64   = 2.997_924_58e8;
const M_E: f64       = 9.109_383_7015e-31;
const M_H: f64       = 1.673_532_9e-27;
const N_A: f64       = 6.022_140_76e23;
const H_PLANCK: f64  = 6.626_070_15e-34;
const R_GAS: f64     = 8.314_462_618;

// ─── Jupiter Parameters ───────────────────────────────────────────────────────

const JUPITER_ROT_S: f64         = 9.925 * 3600.0;
const R_JUPITER: f64             = 71_492_000.0;
const G_JUPITER: f64             = 24.79;
const GRS_LAT_DEG: f64           = -24.0;
const JUPITER_INTERNAL_FLUX: f64 = 5.44;
const T_CLOUD_TOP: f64           = 135.0;
const X_H2: f64                  = 0.864;
const X_HE: f64                  = 0.136;

// ─── GRS Historical Size Data (long axis km, short axis km) ──────────────────

const GRS_HISTORY: &[(f64, f64, f64)] = &[
    (1879.0, 42_000.0, 22_000.0),
    (1920.0, 39_000.0, 20_000.0),
    (1965.0, 35_000.0, 15_000.0),
    (1979.0, 25_800.0, 12_300.0),
    (1995.0, 26_000.0, 12_800.0),
    (2012.0, 25_000.0, 12_000.0),
    (2014.0, 24_500.0, 13_000.0),
    (2017.0, 22_000.0, 12_000.0),
    (2020.0, 21_000.0, 11_500.0),
    (2024.0, 20_000.0, 11_000.0),
];

// ─── Physics Functions ────────────────────────────────────────────────────────

fn coriolis_parameter(lat_deg: f64) -> f64 {
    2.0 * (2.0 * PI / JUPITER_ROT_S) * lat_deg.to_radians().sin().abs()
}

fn pressure_scale_height(t_k: f64) -> f64 {
    let mu = (X_H2 * 2.016 + X_HE * 4.003) * 1e-3;
    R_GAS * t_k / (mu * G_JUPITER)
}

fn brunt_vaisala_freq(h: f64, gamma: f64) -> f64 {
    (G_JUPITER * (gamma - 1.0) / (gamma * h)).max(0.0).sqrt()
}

fn rossby_deformation_radius(n: f64, h: f64, f: f64) -> f64 { n * h / f }
fn rossby_number(u: f64, f: f64, l: f64) -> f64 { u / (f * l) }
fn pv_anomaly(u: f64, l: f64, f: f64) -> f64 { -2.0 * u / l / f }
fn bond_to_uv_nm(bde_kj: f64) -> f64 {
    (H_PLANCK * C_LIGHT / (bde_kj * 1e3 / N_A)) * 1e9
}

fn bohr_radius() -> f64 { HBAR / (ALPHA * M_E * C_LIGHT) }
fn h2_bond() -> f64 { 1.401_09 * bohr_radius() }
fn h2_inertia() -> f64 { 0.5 * M_H * h2_bond().powi(2) }
fn h2_rot_const() -> f64 { HBAR.powi(2) / (2.0 * h2_inertia()) }
fn h2_t_rot() -> f64 { h2_rot_const() / K_B }

fn gamma_eff(t_k: f64) -> f64 {
    let t_rot = h2_t_rot();
    let act = 1.0 / (1.0 + (-2.0 * (t_k / t_rot - 1.0)).exp());
    let f = 3.0 + 2.0 * act;
    (f + 2.0) / f
}

fn fit_linear(data: &[(f64, f64, f64)]) -> (f64, f64) {
    let n = data.len() as f64;
    let (sx, sy, sxy, sx2) = data.iter().fold((0.0, 0.0, 0.0, 0.0), |a, &(x, y, _)| {
        (a.0 + x, a.1 + y, a.2 + x * y, a.3 + x * x)
    });
    let b = (n * sxy - sx * sy) / (n * sx2 - sx * sx);
    ((sy - b * sx) / n, b)
}

fn fit_exp(data: &[(f64, f64, f64)]) -> (f64, f64) {
    let logs: Vec<_> = data.iter().map(|&(x, y, s)| (x, y.ln(), s)).collect();
    fit_linear(&logs)
}

fn proj(a: f64, b: f64, yr: f64, exp: bool) -> f64 {
    if exp { (a + b * yr).exp() } else { a + b * yr }
}

fn yr_at(a: f64, b: f64, crit: f64, exp: bool) -> f64 {
    if exp { (crit.ln() - a) / b } else { (crit - a) / b }
}

// ─── Chromophore Model ────────────────────────────────────────────────────────

#[derive(Clone)]
struct Chromophore {
    name: &'static str,
    precursor: &'static str,
    bde: f64,
    uv_thresh: f64,
    abs_nm: f64,
    color_str: &'static str,
    uv_ok: bool,
    notes: &'static str,
    rgb: [u8; 3],
}

impl Chromophore {
    fn new(
        name: &'static str, precursor: &'static str, bde: f64, abs_nm: f64,
        color: &'static str, notes: &'static str, rgb: [u8; 3],
    ) -> Self {
        let uv_thresh = bond_to_uv_nm(bde);
        let uv_ok = uv_thresh >= 260.0 && uv_thresh <= 420.0;
        Self { name, precursor, bde, uv_thresh, abs_nm, color_str: color, uv_ok, notes, rgb }
    }
    fn score(&self) -> f64 {
        let uv = if self.uv_ok { 1.0 } else { 0.1 };
        uv * (-((self.abs_nm - 550.0) / 150.0).powi(2)).exp()
    }
}

fn chromophores() -> Vec<Chromophore> {
    vec![
        Chromophore::new(
            "Red phosphorus (P4)", "Phosphine (PH3)",
            322.0, 540.0, "orange-red",
            "PH3 photolyzed by UV -> P4 (red allotrope); lab-verified",
            [220, 80, 40],
        ),
        Chromophore::new(
            "Amorphous sulfur (S8)", "Hydrogen sulfide (H2S)",
            381.0, 500.0, "yellow-green",
            "H2S photolysis -> S8; not strongly red",
            [190, 190, 40],
        ),
        Chromophore::new(
            "Disulfide organics (R-S-S-R)", "Organo-sulfur + H2S",
            310.0, 480.0, "yellow-orange",
            "Weak vs phosphorus allotropes",
            [210, 150, 40],
        ),
        Chromophore::new(
            "Ammonium hydrosulfide (NH4SH)", "NH3 + H2S",
            435.0, 450.0, "pale yellow",
            "Cloud layer material; poor chromophore for red color",
            [200, 200, 150],
        ),
        Chromophore::new(
            "PAHs (from C2H2)", "Acetylene (C2H2)",
            390.0, 400.0, "brown",
            "Weak visible absorption; brown tint only",
            [150, 100, 60],
        ),
    ]
}

// ─── Drawing Infrastructure ───────────────────────────────────────────────────

const BG:   [u8; 3] = [10, 13, 22];
const GRID: [u8; 3] = [28, 34, 55];
const AX:   [u8; 3] = [170, 178, 195];
const CYAN: [u8; 3] = [55, 210, 200];
const ORG:  [u8; 3] = [230, 140, 45];
const RED:  [u8; 3] = [220, 65, 65];
const YEL:  [u8; 3] = [220, 200, 55];
const GRY:  [u8; 3] = [105, 112, 128];
const BLU:  [u8; 3] = [80, 140, 230];
const WHT:  [u8; 3] = [215, 220, 230];
const GRN:  [u8; 3] = [75, 200, 100];

/// 5×7 bitmap font for printable ASCII.
/// Each glyph: 7 u8 values (one per row). Bits [4:0] = columns left-to-right.
fn glyph(c: u8) -> [u8; 7] {
    match c {
        b' '  => [0,  0,  0,  0,  0,  0,  0 ],
        b'!'  => [4,  4,  4,  4,  4,  0,  4 ],
        b'"'  => [10, 10, 0,  0,  0,  0,  0 ],
        b'#'  => [10, 31, 10, 10, 31, 10, 0 ],
        b'%'  => [24, 25, 2,  4,  8,  19, 3 ],
        b'('  => [2,  4,  8,  8,  8,  4,  2 ],
        b')'  => [8,  4,  2,  2,  2,  4,  8 ],
        b'*'  => [0,  21, 14, 31, 14, 21, 0 ],
        b'+'  => [0,  4,  4,  31, 4,  4,  0 ],
        b','  => [0,  0,  0,  0,  12, 8,  16],
        b'-'  => [0,  0,  0,  31, 0,  0,  0 ],
        b'.'  => [0,  0,  0,  0,  0,  12, 12],
        b'/'  => [1,  2,  4,  8,  16, 0,  0 ],
        b'0'  => [14, 17, 17, 17, 17, 17, 14],
        b'1'  => [4,  12, 4,  4,  4,  4,  14],
        b'2'  => [14, 17, 1,  6,  8,  16, 31],
        b'3'  => [14, 17, 1,  6,  1,  17, 14],
        b'4'  => [2,  6,  10, 18, 31, 2,  2 ],
        b'5'  => [31, 16, 16, 30, 1,  17, 14],
        b'6'  => [6,  8,  16, 30, 17, 17, 14],
        b'7'  => [31, 1,  2,  4,  8,  8,  8 ],
        b'8'  => [14, 17, 17, 14, 17, 17, 14],
        b'9'  => [14, 17, 17, 15, 1,  2,  12],
        b':'  => [0,  12, 12, 0,  12, 12, 0 ],
        b'<'  => [1,  2,  4,  8,  4,  2,  1 ],
        b'='  => [0,  31, 0,  0,  0,  31, 0 ],
        b'>'  => [16, 8,  4,  2,  4,  8,  16],
        b'?'  => [14, 17, 1,  6,  4,  0,  4 ],
        b'A'  => [4,  10, 17, 17, 31, 17, 17],
        b'B'  => [30, 17, 17, 30, 17, 17, 30],
        b'C'  => [14, 17, 16, 16, 16, 17, 14],
        b'D'  => [30, 17, 17, 17, 17, 17, 30],
        b'E'  => [31, 16, 16, 30, 16, 16, 31],
        b'F'  => [31, 16, 16, 30, 16, 16, 16],
        b'G'  => [14, 17, 16, 23, 17, 17, 14],
        b'H'  => [17, 17, 17, 31, 17, 17, 17],
        b'I'  => [14, 4,  4,  4,  4,  4,  14],
        b'J'  => [7,  2,  2,  2,  18, 18, 12],
        b'K'  => [17, 18, 20, 24, 20, 18, 17],
        b'L'  => [16, 16, 16, 16, 16, 16, 31],
        b'M'  => [17, 27, 21, 17, 17, 17, 17],
        b'N'  => [17, 25, 21, 19, 17, 17, 17],
        b'O'  => [14, 17, 17, 17, 17, 17, 14],
        b'P'  => [30, 17, 17, 30, 16, 16, 16],
        b'Q'  => [14, 17, 17, 17, 21, 18, 13],
        b'R'  => [30, 17, 17, 30, 20, 18, 17],
        b'S'  => [14, 17, 16, 14, 1,  17, 14],
        b'T'  => [31, 4,  4,  4,  4,  4,  4 ],
        b'U'  => [17, 17, 17, 17, 17, 17, 14],
        b'V'  => [17, 17, 17, 17, 17, 10, 4 ],
        b'W'  => [17, 17, 17, 21, 21, 27, 17],
        b'X'  => [17, 17, 10, 4,  10, 17, 17],
        b'Y'  => [17, 17, 10, 4,  4,  4,  4 ],
        b'Z'  => [31, 1,  2,  4,  8,  16, 31],
        b'['  => [14, 8,  8,  8,  8,  8,  14],
        b']'  => [14, 2,  2,  2,  2,  2,  14],
        b'^'  => [4,  10, 17, 0,  0,  0,  0 ],
        b'_'  => [0,  0,  0,  0,  0,  0,  31],
        b'a'  => [0,  0,  14, 1,  15, 17, 15],
        b'b'  => [16, 16, 30, 17, 17, 17, 30],
        b'c'  => [0,  0,  14, 16, 16, 17, 14],
        b'd'  => [1,  1,  15, 17, 17, 17, 15],
        b'e'  => [0,  0,  14, 17, 31, 16, 14],
        b'f'  => [6,  8,  14, 8,  8,  8,  8 ],
        b'g'  => [0,  14, 17, 17, 15, 1,  14],
        b'h'  => [16, 16, 30, 17, 17, 17, 17],
        b'i'  => [4,  0,  12, 4,  4,  4,  14],
        b'j'  => [2,  0,  6,  2,  2,  18, 12],
        b'k'  => [16, 17, 18, 28, 18, 17, 17],
        b'l'  => [12, 4,  4,  4,  4,  4,  14],
        b'm'  => [0,  0,  26, 21, 21, 17, 17],
        b'n'  => [0,  0,  30, 17, 17, 17, 17],
        b'o'  => [0,  0,  14, 17, 17, 17, 14],
        b'p'  => [0,  30, 17, 17, 30, 16, 16],
        b'q'  => [0,  15, 17, 17, 15, 1,  1 ],
        b'r'  => [0,  0,  22, 24, 16, 16, 16],
        b's'  => [0,  0,  15, 16, 14, 1,  30],
        b't'  => [4,  4,  14, 4,  4,  4,  6 ],
        b'u'  => [0,  0,  17, 17, 17, 19, 13],
        b'v'  => [0,  0,  17, 17, 17, 10, 4 ],
        b'w'  => [0,  0,  17, 17, 21, 21, 10],
        b'x'  => [0,  0,  17, 10, 4,  10, 17],
        b'y'  => [0,  17, 17, 15, 1,  1,  14],
        b'z'  => [0,  0,  31, 2,  4,  8,  31],
        _     => [31, 17, 17, 17, 17, 17, 31],
    }
}

struct Canvas {
    w: u32,
    h: u32,
    d: Vec<u8>,
}

impl Canvas {
    fn new(w: u32, h: u32, bg: [u8; 3]) -> Self {
        let mut d = vec![0u8; (w * h * 3) as usize];
        for i in 0..(w * h) as usize {
            d[i * 3] = bg[0];
            d[i * 3 + 1] = bg[1];
            d[i * 3 + 2] = bg[2];
        }
        Self { w, h, d }
    }

    fn set(&mut self, x: i32, y: i32, c: [u8; 3]) {
        if x < 0 || y < 0 || x >= self.w as i32 || y >= self.h as i32 { return; }
        let off = ((y as u32 * self.w + x as u32) * 3) as usize;
        self.d[off] = c[0]; self.d[off + 1] = c[1]; self.d[off + 2] = c[2];
    }

    fn line(&mut self, mut x0: i32, mut y0: i32, x1: i32, y1: i32, c: [u8; 3]) {
        let dx = (x1 - x0).abs(); let dy = (y1 - y0).abs();
        let sx = if x0 < x1 { 1 } else { -1 };
        let sy = if y0 < y1 { 1 } else { -1 };
        let mut err = dx - dy;
        loop {
            self.set(x0, y0, c);
            if x0 == x1 && y0 == y1 { break; }
            let e2 = 2 * err;
            if e2 > -dy { err -= dy; x0 += sx; }
            if e2 < dx { err += dx; y0 += sy; }
        }
    }

    fn thick_line(&mut self, x0: i32, y0: i32, x1: i32, y1: i32, c: [u8; 3], t: i32) {
        let dx = (x1 - x0) as f64; let dy = (y1 - y0) as f64;
        let len = (dx * dx + dy * dy).sqrt().max(1.0);
        let half = t / 2;
        for d in -half..=half {
            let ox = (-dy / len * d as f64).round() as i32;
            let oy = (dx / len * d as f64).round() as i32;
            self.line(x0 + ox, y0 + oy, x1 + ox, y1 + oy, c);
        }
    }

    fn rect(&mut self, x: i32, y: i32, w: i32, h: i32, c: [u8; 3]) {
        for dy in 0..h { for dx in 0..w { self.set(x + dx, y + dy, c); } }
    }

    fn dashed_h(&mut self, y: i32, x0: i32, x1: i32, c: [u8; 3], gap: i32) {
        let mut x = x0;
        while x < x1 {
            for dx in 0..6 { if x + dx < x1 { self.set(x + dx, y, c); } }
            x += 6 + gap;
        }
    }

    fn dashed_v(&mut self, x: i32, y0: i32, y1: i32, c: [u8; 3], gap: i32) {
        let mut y = y0;
        while y < y1 {
            for dy in 0..6 { if y + dy < y1 { self.set(x, y + dy, c); } }
            y += 6 + gap;
        }
    }

    fn filled_circle(&mut self, cx: i32, cy: i32, r: i32, c: [u8; 3]) {
        for dy in -r..=r {
            for dx in -r..=r {
                if dx * dx + dy * dy <= r * r { self.set(cx + dx, cy + dy, c); }
            }
        }
    }

    fn char_px(&mut self, x: i32, y: i32, ch: u8, c: [u8; 3], scale: u32) {
        let bm = glyph(ch);
        let s = scale as i32;
        for (row, &bits) in bm.iter().enumerate() {
            for col in 0..5i32 {
                if bits & (1 << (4 - col)) != 0 {
                    self.rect(x + col * s, y + row as i32 * s, s, s, c);
                }
            }
        }
    }

    fn text(&mut self, mut x: i32, y: i32, s: &str, c: [u8; 3], scale: u32) {
        for ch in s.bytes() {
            self.char_px(x, y, ch, c, scale);
            x += (6 * scale) as i32;
        }
    }

    fn text_center(&mut self, cx: i32, y: i32, s: &str, c: [u8; 3], scale: u32) {
        let w = (s.len() as u32 * 6 * scale) as i32;
        self.text(cx - w / 2, y, s, c, scale);
    }

    fn text_right(&mut self, rx: i32, y: i32, s: &str, c: [u8; 3], scale: u32) {
        let w = (s.len() as u32 * 6 * scale) as i32;
        self.text(rx - w, y, s, c, scale);
    }

    fn save(&self, path: &str) {
        let img: ImageBuffer<Rgb<u8>, Vec<u8>> =
            ImageBuffer::from_raw(self.w, self.h, self.d.clone())
                .expect("ImageBuffer::from_raw");
        img.save(path).expect("PNG save");
    }
}

/// Map world x → pixel x
fn wx(v: f64, v0: f64, v1: f64, p0: i32, p1: i32) -> i32 {
    (p0 as f64 + (v - v0) / (v1 - v0) * (p1 - p0) as f64) as i32
}
/// Map world y → pixel y (flipped: higher world = lower pixel)
fn wy(v: f64, v0: f64, v1: f64, p0: i32, p1: i32) -> i32 {
    (p1 as f64 - (v - v0) / (v1 - v0) * (p1 - p0) as f64) as i32
}

// ─── Chart 1: γ(T) Curve ─────────────────────────────────────────────────────

fn render_gamma_chart(path: &str, t_rot: f64) {
    let (cw, ch) = (900u32, 560u32);
    let (pl, pr, pt, pb) = (100i32, 870i32, 60i32, 490i32); // plot area bounds
    let (xmin, xmax) = (0.0f64, 1500.0);
    let (ymin, ymax) = (1.36f64, 1.70);
    let mut c = Canvas::new(cw, ch, BG);

    // Title
    c.rect(0, 0, cw as i32, 50, [15, 19, 35]);
    c.text_center(cw as i32 / 2, 8, "H2 ADIABATIC INDEX FROM FINE STRUCTURE CONSTANT", WHT, 1);
    c.text_center(cw as i32 / 2, 20, "gamma(T)  via  alpha -> a0 -> r_H2 -> T_rot", AX, 1);

    // Grid
    for &gv in &[1.40f64, 1.45, 1.50, 1.55, 1.60, 1.65] {
        let py = wy(gv, ymin, ymax, pt, pb);
        c.dashed_h(py, pl, pr, GRID, 4);
    }
    for &gv in &[200.0f64, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0] {
        let px = wx(gv, xmin, xmax, pl, pr);
        c.dashed_v(px, pt, pb, GRID, 4);
    }

    // Reference lines: monatomic (5/3) and diatomic (7/5)
    let y_mono = wy(5.0 / 3.0, ymin, ymax, pt, pb);
    c.dashed_h(y_mono, pl, pr, GRY, 3);
    c.text(pr + 5, y_mono - 4, "5/3", GRY, 1);

    let y_di = wy(7.0 / 5.0, ymin, ymax, pt, pb);
    c.dashed_h(y_di, pl, pr, GRY, 3);
    c.text(pr + 5, y_di - 4, "7/5", GRY, 1);

    // T_rot vertical marker
    if t_rot < xmax {
        let px = wx(t_rot, xmin, xmax, pl, pr);
        c.dashed_v(px, pt, pb, ORG, 3);
        c.text_center(px, pt - 24, "T_rot", ORG, 1);
        c.text_center(px, pt - 14, &format!("{t_rot:.0}K"), ORG, 1);
    }

    // T_cloud vertical marker
    let px_cloud = wx(T_CLOUD_TOP, xmin, xmax, pl, pr);
    c.dashed_v(px_cloud, pt, pb, YEL, 3);
    c.text_center(px_cloud, pt - 44, "T_cloud", YEL, 1);
    c.text_center(px_cloud, pt - 34, &format!("{T_CLOUD_TOP:.0}K"), YEL, 1);

    // gamma_eff(T) curve
    let n_pts = 500;
    let mut prev: Option<(i32, i32)> = None;
    for i in 0..=n_pts {
        let t = xmin + (xmax - xmin) * i as f64 / n_pts as f64;
        let g = gamma_eff(t);
        let px = wx(t, xmin, xmax, pl, pr);
        let py = wy(g, ymin, ymax, pt, pb);
        if let Some((ppx, ppy)) = prev {
            c.thick_line(ppx, ppy, px, py, CYAN, 2);
        }
        prev = Some((px, py));
    }

    // Mark cloud-top γ
    let py_cloud = wy(gamma_eff(T_CLOUD_TOP), ymin, ymax, pt, pb);
    c.filled_circle(px_cloud, py_cloud, 5, YEL);

    // Axes
    c.line(pl, pt, pl, pb, AX);
    c.line(pl, pb, pr, pb, AX);

    // Y ticks
    for &gv in &[1.40f64, 1.45, 1.50, 1.55, 1.60, 1.65] {
        let py = wy(gv, ymin, ymax, pt, pb);
        c.line(pl - 5, py, pl, py, AX);
        c.text_right(pl - 8, py - 4, &format!("{gv:.2}"), AX, 1);
    }

    // X ticks
    for &gv in &[0.0f64, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0] {
        let px = wx(gv, xmin, xmax, pl, pr);
        c.line(px, pb, px, pb + 5, AX);
        c.text_center(px, pb + 8, &format!("{gv:.0}"), AX, 1);
    }

    // Axis labels
    c.text_center(cw as i32 / 2, ch as i32 - 14, "Temperature  (K)", AX, 1);
    c.text(6, pt, "gamma", AX, 1);

    // Legend box
    let lx = pr - 190; let ly = pt + 10;
    c.rect(lx - 4, ly - 4, 198, 58, [18, 22, 40]);
    c.rect(lx, ly + 4, 20, 3, CYAN); c.text(lx + 24, ly, "gamma_eff(T) from alpha", WHT, 1);
    c.dashed_h(ly + 18, lx, lx + 20, GRY, 3); c.text(lx + 24, ly + 14, "5/3 monatomic limit", GRY, 1);
    c.dashed_h(ly + 32, lx, lx + 20, GRY, 3); c.text(lx + 24, ly + 28, "7/5 diatomic limit", GRY, 1);
    c.rect(lx, ly + 40, 6, 6, YEL);  c.text(lx + 12, ly + 40, "cloud top (135 K)", YEL, 1);

    c.save(path);
}

// ─── Chart 2: GRS Shrinkage History + Projections ────────────────────────────

fn render_shrinkage_chart(
    path: &str,
    a_lin: f64, b_lin: f64,
    a_exp: f64, b_exp: f64,
    critical_km: f64,
    t_crit_lin: f64, t_crit_exp: f64, t_crit_recent: f64,
) {
    let (cw, ch) = (1000u32, 580u32);
    let (pl, pr, pt, pb) = (115i32, 960i32, 70i32, 500i32);
    let (xmin, xmax) = (1870.0f64, 2260.0);
    let (ymin, ymax) = (0.0f64, 46_000.0);
    let mut c = Canvas::new(cw, ch, BG);

    c.rect(0, 0, cw as i32, 58, [15, 19, 35]);
    c.text_center(cw as i32 / 2, 10, "JUPITER GRS LONG-AXIS SHRINKAGE  1879-2024 + LIFETIME PROJECTIONS", WHT, 1);
    c.text_center(cw as i32 / 2, 24, "Critical stability threshold: 3 x Rossby deformation radius", AX, 1);

    // Grid
    for &gv in &[10_000.0f64, 20_000.0, 30_000.0, 40_000.0] {
        let py = wy(gv, ymin, ymax, pt, pb);
        c.dashed_h(py, pl, pr, GRID, 4);
    }
    for &gv in &[1900.0f64, 1950.0, 2000.0, 2050.0, 2100.0, 2150.0, 2200.0] {
        let px = wx(gv, xmin, xmax, pl, pr);
        c.dashed_v(px, pt, pb, GRID, 4);
    }

    // Critical size line
    let py_crit = wy(critical_km, ymin, ymax, pt, pb);
    c.dashed_h(py_crit, pl, pr, RED, 4);
    c.text(pr + 4, py_crit - 4, "3xL_D", RED, 1);

    // Dissolution year markers
    for &(yr, col) in &[(t_crit_lin, BLU), (t_crit_exp, ORG), (t_crit_recent, RED)] {
        if yr > xmin && yr < xmax {
            let px = wx(yr, xmin, xmax, pl, pr);
            c.dashed_v(px, py_crit, pb + 15, col, 3);
        }
    }

    // Linear fit line (extended)
    {
        let mut prev: Option<(i32, i32)> = None;
        for i in 0..=300 {
            let yr = xmin + (xmax - xmin) * i as f64 / 300.0;
            let size = proj(a_lin, b_lin, yr, false);
            if size < 0.0 { break; }
            let px = wx(yr, xmin, xmax, pl, pr);
            let py = wy(size.min(ymax), ymin, ymax, pt, pb);
            if let Some((ppx, ppy)) = prev { c.thick_line(ppx, ppy, px, py, BLU, 2); }
            prev = Some((px, py));
        }
    }

    // Exponential fit line
    {
        let mut prev: Option<(i32, i32)> = None;
        for i in 0..=300 {
            let yr = xmin + (xmax - xmin) * i as f64 / 300.0;
            let size = proj(a_exp, b_exp, yr, true);
            if size > ymax * 1.5 { prev = None; continue; }
            let px = wx(yr, xmin, xmax, pl, pr);
            let py = wy(size.min(ymax), ymin, ymax, pt, pb);
            if let Some((ppx, ppy)) = prev { c.thick_line(ppx, ppy, px, py, ORG, 2); }
            prev = Some((px, py));
        }
    }

    // Recent rate projection (dashed, from 2017)
    {
        let yr_start = 2017.0;
        let size_start = 22_000.0_f64;
        let rate = (20_000.0 - 22_000.0) / (2024.0 - 2017.0);
        let mut prev: Option<(i32, i32)> = None;
        for i in 0..=100 {
            let yr = yr_start + i as f64 * 0.5;
            if yr > xmax { break; }
            let size = size_start + rate * (yr - yr_start);
            if size < 0.0 { break; }
            let px = wx(yr, xmin, xmax, pl, pr);
            let py = wy(size.min(ymax), ymin, ymax, pt, pb);
            if let Some((ppx, ppy)) = prev {
                if i % 3 < 2 { c.thick_line(ppx, ppy, px, py, RED, 2); }
            }
            prev = Some((px, py));
        }
    }

    // Historical data points
    for &(yr, long_km, _) in GRS_HISTORY {
        let px = wx(yr, xmin, xmax, pl, pr);
        let py = wy(long_km, ymin, ymax, pt, pb);
        c.filled_circle(px, py, 6, WHT);
        c.filled_circle(px, py, 4, [80, 140, 200]);
    }

    // Axes
    c.line(pl, pt, pl, pb, AX);
    c.line(pl, pb, pr, pb, AX);

    // Y ticks (in km, labelled in 10k units)
    for &gv in &[0.0f64, 10_000.0, 20_000.0, 30_000.0, 40_000.0] {
        let py = wy(gv, ymin, ymax, pt, pb);
        c.line(pl - 5, py, pl, py, AX);
        c.text_right(pl - 8, py - 4, &format!("{:.0}k", gv / 1000.0), AX, 1);
    }

    // X ticks
    for &gv in &[1880.0f64, 1920.0, 1960.0, 2000.0, 2040.0, 2080.0, 2120.0, 2160.0, 2200.0] {
        let px = wx(gv, xmin, xmax, pl, pr);
        c.line(px, pb, px, pb + 5, AX);
        c.text_center(px, pb + 8, &format!("{gv:.0}"), AX, 1);
    }

    c.text_center(cw as i32 / 2, ch as i32 - 14, "Year", AX, 1);
    c.text(4, pt, "km", AX, 1);

    // Legend
    let lx = pl + 10; let ly = pt + 8;
    c.rect(lx - 4, ly - 4, 260, 78, [18, 22, 40]);
    c.filled_circle(lx + 10, ly + 6, 6, WHT);
    c.filled_circle(lx + 10, ly + 6, 4, [80, 140, 200]);
    c.text(lx + 20, ly, "Observed (historical)", WHT, 1);
    c.rect(lx, ly + 18, 20, 3, BLU);   c.text(lx + 24, ly + 14, "Linear model", BLU, 1);
    c.rect(lx, ly + 32, 20, 3, ORG);   c.text(lx + 24, ly + 28, "Exponential model", ORG, 1);
    c.rect(lx, ly + 46, 20, 3, RED);   c.text(lx + 24, ly + 42, "Recent rate (2017-2024)", RED, 1);
    c.dashed_h(ly + 60, lx, lx + 20, RED, 3); c.text(lx + 24, ly + 56, &format!("3xL_D critical ({critical_km:.0} km)"), RED, 1);

    c.save(path);
}

// ─── Chart 3: Chromophore Ranking ────────────────────────────────────────────

fn render_chromophore_chart(path: &str, panel: &[Chromophore]) {
    let (cw, ch) = (860u32, 520u32);
    let mut c = Canvas::new(cw, ch, BG);

    c.rect(0, 0, cw as i32, 52, [15, 19, 35]);
    c.text_center(cw as i32 / 2, 8, "CHROMOPHORE RANKING — JUPITER GRS RED COLOR", WHT, 1);
    c.text_center(cw as i32 / 2, 22, "Score = UV accessibility x color match Gaussian (peak at 550 nm)", AX, 1);
    c.text_center(cw as i32 / 2, 36, "Bond energetics -> UV threshold -> photolysis -> chromophore -> apparent color", GRY, 1);

    let bar_x0 = 310i32;
    let bar_x1 = 800i32;
    let row_h = 76i32;
    let row0 = 72i32;

    // Score axis header
    c.text(bar_x0, row0 - 18, "0.0", GRY, 1);
    c.text_center((bar_x0 + bar_x1) / 2, row0 - 18, "Score", AX, 1);
    c.text_right(bar_x1, row0 - 18, "1.0", GRY, 1);
    c.line(bar_x0, row0 - 8, bar_x1, row0 - 8, GRY);

    for (i, ch_c) in panel.iter().enumerate() {
        let y = row0 + i as i32 * row_h;
        let score = ch_c.score();

        // Row background (alternating)
        let bg_col = if i % 2 == 0 { [14, 18, 30] } else { [18, 22, 36] };
        c.rect(0, y, cw as i32, row_h - 2, bg_col);

        // Name + metadata (left column)
        let rank_col = if i == 0 { YEL } else { AX };
        c.text(8, y + 8, &format!("#{}", i + 1), rank_col, 2);
        c.text(32, y + 8, ch_c.name, WHT, 1);
        c.text(32, y + 20, &format!("from: {}", ch_c.precursor), GRY, 1);
        c.text(32, y + 32, &format!("BDE: {} kJ/mol  UV: {:.0} nm  Abs: {:.0} nm",
            ch_c.bde as i32, ch_c.uv_thresh, ch_c.abs_nm), AX, 1);
        let uv_txt = if ch_c.uv_ok { "UV: YES" } else { "UV: no" };
        let uv_col = if ch_c.uv_ok { GRN } else { RED };
        c.text(32, y + 44, uv_txt, uv_col, 1);
        c.text(90, y + 44, &format!("color: {}", ch_c.color_str), AX, 1);

        // Color swatch
        c.rect(bar_x0 - 22, y + 8, 16, 30, ch_c.rgb);

        // Score bar
        let bar_w = ((score * (bar_x1 - bar_x0) as f64) as i32).max(2);
        c.rect(bar_x0, y + 12, bar_w, 28, ch_c.rgb);

        // Score fraction bar (background)
        c.rect(bar_x0 + bar_w, y + 12, bar_x1 - bar_x0 - bar_w, 28, [25, 30, 50]);

        // Score text
        c.text(bar_x0 + bar_w + 6, y + 22, &format!("{score:.3}"), WHT, 1);

        // Grid line at 0.5
        let half_x = wx(0.5, 0.0, 1.0, bar_x0, bar_x1);
        c.dashed_v(half_x, y + 12, y + 40, GRY, 3);
    }

    // Winner annotation
    c.text(8, row0 + panel.len() as i32 * row_h + 8,
        "WINNER: Red phosphorus (P4) — PH3 UV photolysis at 371 nm -> P4 -> absorbs at ~540 nm -> orange-red", YEL, 1);

    c.save(path);
}

// ─── Chart 4: Stability Margin Over Time ─────────────────────────────────────

fn render_stability_chart(
    path: &str,
    a_lin: f64, b_lin: f64,
    a_exp: f64, b_exp: f64,
    critical_km: f64,
    t_crit_lin: f64, t_crit_exp: f64, t_crit_recent: f64,
) {
    let (cw, ch) = (900u32, 560u32);
    let (pl, pr, pt, pb) = (100i32, 860i32, 60i32, 490i32);
    let (xmin, xmax) = (1870.0f64, 2260.0);
    let (ymin, ymax) = (0.0f64, 9.0);
    let mut c = Canvas::new(cw, ch, BG);

    c.rect(0, 0, cw as i32, 50, [15, 19, 35]);
    c.text_center(cw as i32 / 2, 8, "GRS STABILITY MARGIN = size / critical_size  OVER TIME", WHT, 1);
    c.text_center(cw as i32 / 2, 22, "Vortex dissolves when ratio falls below 1.0  (3 x Rossby deformation radius)", AX, 1);

    // Grid
    for &gv in &[1.0f64, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0] {
        let py = wy(gv, ymin, ymax, pt, pb);
        c.dashed_h(py, pl, pr, GRID, 4);
    }
    for &gv in &[1900.0f64, 1950.0, 2000.0, 2050.0, 2100.0, 2150.0, 2200.0] {
        let px = wx(gv, xmin, xmax, pl, pr);
        c.dashed_v(px, pt, pb, GRID, 4);
    }

    // Critical ratio = 1.0 (dissolution boundary)
    let py_crit = wy(1.0, ymin, ymax, pt, pb);
    c.thick_line(pl, py_crit, pr, py_crit, RED, 2);
    c.text(pr + 4, py_crit - 4, "1.0", RED, 1);
    c.text(pr + 4, py_crit + 4, "dissolve", RED, 1);

    // Safe margin indicator at ratio=3 (current)
    let py3 = wy(3.0, ymin, ymax, pt, pb);
    c.dashed_h(py3, pl, pr, [80, 160, 80], 4);
    c.text(pr + 4, py3 - 4, "3.0", GRN, 1);

    // Dissolution year markers
    for &(yr, col) in &[(t_crit_lin, BLU), (t_crit_exp, ORG), (t_crit_recent, RED)] {
        if yr > xmin && yr < xmax {
            let px = wx(yr, xmin, xmax, pl, pr);
            c.dashed_v(px, py_crit, pb, col, 3);
            c.text_center(px, pb + 8, &format!("{yr:.0}"), col, 1);
        }
    }

    // Linear model ratio
    {
        let mut prev: Option<(i32, i32)> = None;
        for i in 0..=300 {
            let yr = xmin + (xmax - xmin) * i as f64 / 300.0;
            let size = proj(a_lin, b_lin, yr, false);
            if size <= 0.0 { break; }
            let ratio = size / critical_km;
            let px = wx(yr, xmin, xmax, pl, pr);
            let py = wy(ratio.min(ymax), ymin, ymax, pt, pb);
            if let Some((ppx, ppy)) = prev { c.thick_line(ppx, ppy, px, py, BLU, 2); }
            prev = Some((px, py));
        }
    }

    // Exponential model ratio
    {
        let mut prev: Option<(i32, i32)> = None;
        for i in 0..=300 {
            let yr = xmin + (xmax - xmin) * i as f64 / 300.0;
            let size = proj(a_exp, b_exp, yr, true);
            let ratio = size / critical_km;
            if ratio > ymax * 1.5 { prev = None; continue; }
            let px = wx(yr, xmin, xmax, pl, pr);
            let py = wy(ratio.min(ymax), ymin, ymax, pt, pb);
            if let Some((ppx, ppy)) = prev { c.thick_line(ppx, ppy, px, py, ORG, 2); }
            prev = Some((px, py));
        }
    }

    // Recent rate ratio
    {
        let yr_start = 2017.0;
        let size_start = 22_000.0_f64;
        let rate = (20_000.0 - 22_000.0) / (2024.0 - 2017.0);
        let mut prev: Option<(i32, i32)> = None;
        for i in 0..=100 {
            let yr = yr_start + i as f64 * 0.5;
            if yr > xmax { break; }
            let size = size_start + rate * (yr - yr_start);
            if size <= 0.0 { break; }
            let ratio = size / critical_km;
            let px = wx(yr, xmin, xmax, pl, pr);
            let py = wy(ratio.min(ymax), ymin, ymax, pt, pb);
            if let Some((ppx, ppy)) = prev {
                if i % 3 < 2 { c.thick_line(ppx, ppy, px, py, RED, 2); }
            }
            prev = Some((px, py));
        }
    }

    // Historical data points (ratio = observed / critical)
    for &(yr, long_km, _) in GRS_HISTORY {
        let ratio = long_km / critical_km;
        let px = wx(yr, xmin, xmax, pl, pr);
        let py = wy(ratio.min(ymax), ymin, ymax, pt, pb);
        c.filled_circle(px, py, 6, WHT);
        c.filled_circle(px, py, 4, [80, 140, 200]);
    }

    // Axes
    c.line(pl, pt, pl, pb, AX);
    c.line(pl, pb, pr, pb, AX);

    // Y ticks
    for &gv in &[0.0f64, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0] {
        let py = wy(gv, ymin, ymax, pt, pb);
        c.line(pl - 5, py, pl, py, AX);
        c.text_right(pl - 8, py - 4, &format!("{gv:.0}"), AX, 1);
    }

    // X ticks
    for &gv in &[1880.0f64, 1920.0, 1960.0, 2000.0, 2040.0, 2080.0, 2120.0, 2160.0, 2200.0] {
        let px = wx(gv, xmin, xmax, pl, pr);
        c.line(px, pb, px, pb + 5, AX);
        c.text_center(px, pb + 8, &format!("{gv:.0}"), AX, 1);
    }

    c.text_center(cw as i32 / 2, ch as i32 - 14, "Year", AX, 1);
    c.text(4, pt, "size/critical", AX, 1);

    // Legend
    let lx = pl + 10; let ly = pt + 10;
    c.rect(lx - 4, ly - 4, 260, 68, [18, 22, 40]);
    c.filled_circle(lx + 10, ly + 6, 6, WHT);
    c.filled_circle(lx + 10, ly + 6, 4, [80, 140, 200]);
    c.text(lx + 20, ly, "Observed ratio", WHT, 1);
    c.rect(lx, ly + 18, 20, 3, BLU); c.text(lx + 24, ly + 14, "Linear model", BLU, 1);
    c.rect(lx, ly + 32, 20, 3, ORG); c.text(lx + 24, ly + 28, "Exponential model", ORG, 1);
    c.rect(lx, ly + 46, 20, 3, RED); c.text(lx + 24, ly + 42, "Recent rate (2017-2024)", RED, 1);

    c.save(path);
}

// ─── Text Report ──────────────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn write_text_report(
    path: &str,
    a0: f64, r_h2: f64, i_h2: f64, b_rot: f64, t_rot: f64,
    gamma_cloud: f64, gamma_deep: f64,
    f_cor: f64, h_scale: f64, n_bv: f64, l_d: f64,
    u_grs: f64, ro: f64, pv: f64,
    panel: &[Chromophore],
    a_lin: f64, b_lin: f64, a_exp: f64, b_exp: f64,
    proj_lin_2024: f64, proj_exp_2024: f64,
    critical_km: f64,
    t_crit_lin: f64, t_crit_exp: f64, t_crit_recent: f64,
    recent_rate: f64,
    assert_count: usize, fail_count: usize,
) {
    let winner = &panel[0];
    let mut s = String::new();
    let push = |s: &mut String, line: &str| { s.push_str(line); s.push('\n'); };

    push(&mut s, "╔═══════════════════════════════════════════════════════════════════════════╗");
    push(&mut s, "║   GUTOE — JUPITER GREAT RED SPOT: DYNAMICS, COLOR, AND LIFETIME          ║");
    push(&mut s, "║   A 350-year mystery from first principles (GUTOE-Physics)                ║");
    push(&mut s, "╚═══════════════════════════════════════════════════════════════════════════╝");
    push(&mut s, "");

    push(&mut s, "═══ A.  α → H₂ → Jovian Atmospheric Engine ═══════════════════════════════");
    push(&mut s, "");
    push(&mut s, &format!("  Fine structure constant  α = {ALPHA:.9}"));
    push(&mut s, &format!("  Bohr radius  a₀ = ħ/(α m_e c) = {:.6e} m  ({:.4} pm)", a0, a0 * 1e12));
    push(&mut s, "  Literature:                                 52.9177 pm  (agreement: <0.001%)");
    push(&mut s, &format!("  H₂ bond length  r₀ = 1.401 a₀  = {:.4} pm", r_h2 * 1e12));
    push(&mut s, "  Literature:                       74.14 pm");
    push(&mut s, &format!("  Reduced mass  μ = m_H/2         = {:.4e} kg", M_H / 2.0));
    push(&mut s, &format!("  Moment of inertia  I = μ r₀²    = {:.4e} kg·m²", i_h2));
    push(&mut s, &format!("  Rotational constant  B = ħ²/2I  = {:.4e} J", b_rot));
    push(&mut s, &format!("  Rotational temperature T_rot = B/k_B = {:.2} K", t_rot));
    push(&mut s, "  Literature:  T_rot(H₂) ≈ 85–87 K");
    push(&mut s, "");
    push(&mut s, &format!("  Effective adiabatic index γ at cloud tops ({T_CLOUD_TOP:.0} K):"));
    push(&mut s, &format!("    γ_eff = {gamma_cloud:.4}  (transitional: T ≈ 1.6 × T_rot)"));
    push(&mut s, &format!("    γ_eff at deep troposphere (1500 K) = {gamma_deep:.4}  (fully diatomic)"));
    push(&mut s, "    Monatomic limit (T << T_rot): γ = 1.6667");
    push(&mut s, "    Diatomic limit  (T >> T_rot): γ = 1.4000");
    push(&mut s, "");
    push(&mut s, "  Physical interpretation:");
    push(&mut s, "    At T=135 K ≈ 1.6×T_rot, H₂ rotation is partially activated.");
    push(&mut s, "    The transitional γ affects the adiabatic lapse rate and thus the");
    push(&mut s, "    Brunt-Väisälä stability frequency N — a key vortex parameter.");
    push(&mut s, "");

    push(&mut s, "═══ B.  Vortex Dynamics — Why the GRS Persists ════════════════════════════");
    push(&mut s, "");
    push(&mut s, &format!("  Coriolis parameter  f = 2ω sin(φ)         = {f_cor:.4e} rad/s  (at 24°S)"));
    push(&mut s, &format!("  Pressure scale height  H = RT/(Mg)        = {:.2} km", h_scale / 1e3));
    push(&mut s, &format!("  Brunt-Väisälä frequency  N² = g(γ-1)/(γH) = {n_bv:.4e} rad/s"));
    push(&mut s, &format!("  Rossby deformation radius  L_D = NH/f     = {:.0} km", l_d / 1e3));
    push(&mut s, "");
    push(&mut s, &format!("  GRS long axis (2024):  20,000 km = {:.1}× L_D", 20_000.0e3 / l_d));
    push(&mut s, &format!("  Peak wind speed  U                         = {u_grs:.0} m/s"));
    push(&mut s, &format!("  Rossby number  Ro = U/(f·L)               = {ro:.4}  << 1"));
    push(&mut s, &format!("  Potential vorticity anomaly  q/f          = {pv:.4}  (anticyclonic)"));
    push(&mut s, &format!("  Vertical extent  ≈ 3H                     = {:.0} km", 3.0 * h_scale / 1e3));
    push(&mut s, &format!("  Critical stability size  = 3×L_D          = {critical_km:.0} km"));
    push(&mut s, &format!("  Current safety margin                      = {:.1}× critical", 20_000.0 / critical_km));
    push(&mut s, "");
    push(&mut s, "  Why the GRS persists:");
    push(&mut s, "    Ro << 1  →  Coriolis dominates pressure gradient (quasi-geostrophic).");
    push(&mut s, "    Geostrophic balance maintains the anticyclone against viscous decay.");
    push(&mut s, &format!("    Internal heat flux ({JUPITER_INTERNAL_FLUX} W/m²) feeds baroclinic instability"));
    push(&mut s, "    → baroclinic energy input ≈ viscous dissipation → energy replenishment.");
    push(&mut s, "    Anticyclones in the SH pump toward high pressure → self-reinforcing.");
    push(&mut s, "    Size >> 3×L_D means the vortex is in the deeply stable regime.");
    push(&mut s, "");

    push(&mut s, "═══ C.  Chromophore Chemistry — Why the GRS is Red ════════════════════════");
    push(&mut s, "");
    push(&mut s, "  UV penetration at GRS cloud tops: ~260–420 nm (UV-B/A window)");
    push(&mut s, "  Bond photolysis threshold: λ_thresh = N_A h c / BDE");
    push(&mut s, "");
    push(&mut s, &format!("  {:30}  {:>10}  {:>12}  {:>12}  {:>5}  {:>6}",
        "Chromophore", "BDE kJ/mol", "UV thresh nm", "Abs peak nm", "UV?", "Score"));
    push(&mut s, &format!("  {}", "─".repeat(88)));
    for (i, ch_c) in panel.iter().enumerate() {
        push(&mut s, &format!("  #{} {:29}  {:>10.1}  {:>12.1}  {:>12.1}  {:>5}  {:>6.3}",
            i + 1, ch_c.name, ch_c.bde, ch_c.uv_thresh, ch_c.abs_nm,
            if ch_c.uv_ok { "YES" } else { "no" }, ch_c.score()));
    }
    push(&mut s, "");
    push(&mut s, &format!("  WINNER: {}", winner.name));
    push(&mut s, &format!("  Precursor: {}", winner.precursor));
    push(&mut s, &format!("  UV threshold: {:.1} nm — within Jupiter's UV-B window (260–420 nm)", winner.uv_thresh));
    push(&mut s, &format!("  Absorption peak: ~{:.0} nm (blue/green) → reflects orange-red", winner.abs_nm));
    push(&mut s, &format!("  Notes: {}", winner.notes));
    push(&mut s, "  Reference: Sagan & Khare (1981) confirmed P₄ formation at Jovian conditions.");
    push(&mut s, "");

    push(&mut s, "═══ D.  Shrinkage History, Rate, and Lifetime Projections ══════════════════");
    push(&mut s, "");
    push(&mut s, "  Historical GRS long axis (long km × short km):");
    for &(yr, long, short) in GRS_HISTORY {
        push(&mut s, &format!("    {yr:.0}: {long:.0} × {short:.0} km"));
    }
    push(&mut s, "");
    push(&mut s, &format!("  Linear model:      size = {a_lin:.0} + {b_lin:.1} × year"));
    push(&mut s, &format!("    Shrinkage rate: {:.0} km/year", b_lin.abs()));
    push(&mut s, &format!("    Prediction at 2024: {proj_lin_2024:.0} km  (observed: 20,000 km)"));
    push(&mut s, "");
    push(&mut s, &format!("  Exponential model: ln(size) = {a_exp:.4} + {b_exp:.6} × year"));
    push(&mut s, &format!("    Half-life: {:.0} years", -1.0 / b_exp * 2.0_f64.ln()));
    push(&mut s, &format!("    Prediction at 2024: {proj_exp_2024:.0} km  (observed: 20,000 km)"));
    push(&mut s, "");
    push(&mut s, &format!("  Critical size (3×L_D): {critical_km:.0} km"));
    push(&mut s, "");
    push(&mut s, "  LIFETIME PROJECTIONS:");
    push(&mut s, &format!("    Linear model     → critical in {t_crit_lin:.0}  ({:.0} years from 2024)",
        t_crit_lin - 2024.0));
    push(&mut s, &format!("    Exponential model → critical in {t_crit_exp:.0}  ({:.0} years from 2024)",
        t_crit_exp - 2024.0));
    push(&mut s, &format!("    Recent rate (2017-2024): {recent_rate:.0} km/year → critical in {t_crit_recent:.0}  ({:.0} years from 2024)",
        t_crit_recent - 2024.0));
    push(&mut s, "    Uncertainty window: 2050–2220");
    push(&mut s, "");
    push(&mut s, "  WHY IT'S SHRINKING:");
    push(&mut s, "    1. Filamentary erosion — GRS sheds vortex filaments into the S. Equatorial Belt,");
    push(&mut s, "       carrying away potential vorticity → vortex area shrinks.");
    push(&mut s, "    2. Reduced baroclinic forcing — latitudinal jet structure weakening since 1980s.");
    push(&mut s, "    3. Stretching/elongation — long axis shrinks faster than short, area decreases.");
    push(&mut s, "    NOT viscous dissipation — that timescale >> Jupiter's age.");
    push(&mut s, "    Mechanism: dynamical vortex-vortex interaction + jet shear erosion.");
    push(&mut s, "");

    push(&mut s, "═══ COMPLETE CHAIN ════════════════════════════════════════════════════════");
    push(&mut s, "");
    push(&mut s, &format!("  α = {ALPHA:.9}"));
    push(&mut s, &format!("  → a₀ = {:.3} pm  → r_H₂ = {:.3} pm  → T_rot = {:.1} K", a0*1e12, r_h2*1e12, t_rot));
    push(&mut s, &format!("  → γ(135K) = {gamma_cloud:.4}  → H = {:.1} km  → N = {n_bv:.4e} rad/s", h_scale/1e3));
    push(&mut s, &format!("  → L_D = {:.0} km  → GRS {:.0} km >> 3×L_D = {critical_km:.0} km  → STABLE", l_d/1e3, 20_000.0));
    push(&mut s, &format!("  → Ro = {ro:.3} << 1  → Pv = {pv:.3} < 0  → anticyclone locked"));
    push(&mut s, &format!("  → Color: P₄ from PH₃ at {:.0} nm → absorbs {:.0} nm → orange-red",
        winner.uv_thresh, winner.abs_nm));
    push(&mut s, &format!("  → Shrinkage: {:.0} km/yr → critical ~{t_crit_lin:.0}", b_lin.abs()));
    push(&mut s, "");

    push(&mut s, &format!("═══ ASSERTIONS: {}/{} pass ═══════════════════════════════════════════════",
        assert_count - fail_count, assert_count));

    fs::write(path, &s).expect("write text report");
}

// ─── CSV Writers ──────────────────────────────────────────────────────────────

fn write_csvs(
    out: &str, panel: &[Chromophore],
    a_lin: f64, b_lin: f64, a_exp: f64, b_exp: f64, critical_km: f64,
) {
    // gamma curve
    let mut g = String::from("temperature_K,gamma_eff\n");
    for t in (10..=3000).step_by(5) {
        g.push_str(&format!("{},{:.6}\n", t, gamma_eff(t as f64)));
    }
    fs::write(format!("{out}/grs_data_gamma.csv"), &g).expect("write gamma csv");

    // shrinkage
    let mut s = String::from("year,observed_long_km,observed_short_km,linear_model_km,exp_model_km,critical_km\n");
    for yr in (1870..=2260).step_by(2) {
        let yr_f = yr as f64;
        let obs_l = GRS_HISTORY.iter()
            .find(|&&(y, _, _)| (y - yr_f).abs() < 1.0)
            .map(|&(_, l, _)| l);
        let obs_s = GRS_HISTORY.iter()
            .find(|&&(y, _, _)| (y - yr_f).abs() < 1.0)
            .map(|&(_, _, s)| s);
        let lin = proj(a_lin, b_lin, yr_f, false).max(0.0);
        let exp = proj(a_exp, b_exp, yr_f, true);
        s.push_str(&format!("{},{},{},{:.0},{:.0},{:.0}\n",
            yr,
            obs_l.map(|v| format!("{v:.0}")).unwrap_or_default(),
            obs_s.map(|v| format!("{v:.0}")).unwrap_or_default(),
            lin, exp, critical_km));
    }
    fs::write(format!("{out}/grs_data_shrinkage.csv"), &s).expect("write shrinkage csv");

    // chromophores
    let mut ch = String::from("rank,name,precursor,bde_kj_per_mol,uv_threshold_nm,abs_peak_nm,color,uv_accessible,score\n");
    for (i, c) in panel.iter().enumerate() {
        ch.push_str(&format!("{},{},{},{:.1},{:.1},{:.1},{},{},{:.4}\n",
            i + 1, c.name, c.precursor, c.bde, c.uv_thresh, c.abs_nm,
            c.color_str, c.uv_ok as u8, c.score()));
    }
    fs::write(format!("{out}/grs_data_chromophores.csv"), &ch).expect("write chromophores csv");
}

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let mut assert_count = 0usize;
    let mut fails: Vec<String> = Vec::new();

    macro_rules! check {
        ($cond:expr, $msg:expr) => {{
            assert_count += 1;
            if !($cond) { fails.push(format!("FAIL [{}]: {}", assert_count, $msg)); }
        }};
    }

    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  GUTOE — Jupiter Great Red Spot: Dynamics, Color, Lifetime");
    println!("  A 350-year-old mystery from first principles");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    // ── A: α → H₂ → γ ────────────────────────────────────────────────────────

    println!("═══ A.  α → H₂ → Jovian Atmospheric Engine ═══════════════════════════\n");

    let a0         = bohr_radius();
    let r_h2       = h2_bond();
    let i_h2       = h2_inertia();
    let b_rot      = h2_rot_const();
    let t_rot      = h2_t_rot();
    let gamma_cloud = gamma_eff(T_CLOUD_TOP);
    let gamma_deep  = gamma_eff(1500.0);
    let gamma_mono  = 5.0 / 3.0;
    let gamma_di    = 7.0 / 5.0;

    println!("  Bohr radius  a₀ = ħ/(α m_e c) = {:.4e} m  ({:.4} pm)", a0, a0 * 1e12);
    println!("  Literature:  a₀ = 5.2918e-11 m = 52.918 pm");
    println!("  H₂ bond     r₀  = 1.401 a₀    = {:.4} pm", r_h2 * 1e12);
    println!("  Literature:  r₀ = 74.14 pm");
    println!("  Moment of inertia I = {:.4e} kg·m²", i_h2);
    println!("  Rot. constant  B   = {:.4e} J", b_rot);
    println!("  Rot. temperature T_rot = B/k_B = {:.2} K", t_rot);
    println!("  Literature:  T_rot(H₂) ≈ 85–87 K\n");
    println!("  γ_eff at T_cloud = {:.0} K  →  γ = {:.4}", T_CLOUD_TOP, gamma_cloud);
    println!("  γ_eff at T_deep  = 1500 K   →  γ = {:.4}", gamma_deep);
    println!("  Monatomic limit  (T << T_rot): γ = {:.4}", gamma_mono);
    println!("  Diatomic limit   (T >> T_rot): γ = {:.4}", gamma_di);
    println!("  → At cloud tops T={:.0}K ≈ 1.6×T_rot: γ={:.4} (transitional, not fully diatomic)",
             T_CLOUD_TOP, gamma_cloud);

    check!((a0 - 5.292e-11).abs() / 5.292e-11 < 0.01, "Bohr radius within 1% of 52.92 pm");
    check!((r_h2 * 1e12 - 74.0).abs() < 2.0, "H₂ bond length within 2 pm of 74 pm");
    check!((t_rot - 86.0).abs() < 4.0, "H₂ rotational temperature within 4 K of 86 K");
    check!(gamma_cloud > gamma_di && gamma_cloud < gamma_mono,
           "γ at cloud tops is between di- and monatomic limits");
    check!(gamma_deep < gamma_cloud, "deeper (hotter) atmosphere has lower γ");

    // ── B: Vortex dynamics ────────────────────────────────────────────────────

    println!("\n═══ B.  Vortex Dynamics — Why GRS Persists ════════════════════════════\n");

    let f_cor      = coriolis_parameter(GRS_LAT_DEG);
    let h_scale    = pressure_scale_height(T_CLOUD_TOP);
    let n_bv       = brunt_vaisala_freq(h_scale, gamma_cloud);
    let l_d        = rossby_deformation_radius(n_bv, h_scale, f_cor);
    let u_grs      = 130.0_f64;
    let l_grs_m    = 20_000.0e3_f64;
    let ro         = rossby_number(u_grs, f_cor, l_grs_m / 2.0);
    let pv         = pv_anomaly(u_grs, l_grs_m / 2.0, f_cor);

    println!("  Coriolis f    = {:.4e} rad/s  (at {:.0}°S)", f_cor, GRS_LAT_DEG.abs());
    println!("  Scale height H = {:.2} km  (T={:.0} K, g={:.2} m/s²)", h_scale / 1e3, T_CLOUD_TOP, G_JUPITER);
    println!("  Brunt-Väisälä N = {:.4e} rad/s", n_bv);
    println!("  Rossby deformation radius L_D = N·H/f = {:.0} km", l_d / 1e3);
    println!("  GRS long axis (2024): {:.0} km  =  {:.1}× L_D", l_grs_m / 1e3, l_grs_m / l_d);
    println!();
    println!("  Peak wind speed U = {:.0} m/s", u_grs);
    println!("  Rossby number   Ro = U/(f·L) = {:.4}  << 1 → strongly geostrophic", ro);
    println!("  PV anomaly  q/f   = {:.4}  (negative = anticyclonic SH stabiliser)", pv);
    println!();
    let critical_km = 3.0 * l_d / 1e3;
    println!("  Stability: GRS ({:.0} km) >> 3·L_D ({:.0} km)", l_grs_m / 1e3, critical_km);
    println!("  → Vortex is STABLE (deeply geostrophic, persistent anticyclone)");
    println!("  Vertical extent ≈ 3 scale heights = {:.0} km", 3.0 * h_scale / 1e3);
    println!("  Critical stability size = 3·L_D = {:.0} km", critical_km);

    let omega_j = 2.0 * PI / JUPITER_ROT_S;
    check!((f_cor - 1.4e-4).abs() / 1.4e-4 < 0.2, "Coriolis f ≈ 1.4×10⁻⁴ rad/s");
    check!((h_scale / 1e3 - 20.0).abs() < 15.0, "Scale height 5–35 km at T_cloud");
    check!(ro < 0.15, "Rossby number << 1 (strongly geostrophic)");
    check!(l_grs_m / l_d > 3.0, "GRS size > 3 × Rossby deformation radius");
    check!(pv < 0.0, "PV anomaly negative = anticyclonic");
    check!((omega_j * R_JUPITER - 12_600.0).abs() < 1500.0,
           "Equatorial rotation speed ≈ 12.6 km/s");

    // ── C: Chromophore chemistry ──────────────────────────────────────────────

    println!("\n═══ C.  Chromophore Chemistry — Why the GRS is Red ════════════════════\n");
    println!("  UV penetration at GRS cloud tops: ~260–420 nm");
    println!("  {:30} {:12} {:14} {:14} {:6} Score",
             "Chromophore", "BDE kJ/mol", "UV thresh nm", "Abs peak nm", "UV?");
    println!("  {}", "─".repeat(100));

    let mut panel = chromophores();
    panel.sort_by(|a, b| b.score().partial_cmp(&a.score()).unwrap());

    for (i, c) in panel.iter().enumerate() {
        println!("  #{} {:29} {:12.1} {:14.1} {:14.1} {:6} {:.3}",
                 i + 1, c.name, c.bde, c.uv_thresh, c.abs_nm,
                 if c.uv_ok { "YES" } else { "no" }, c.score());
    }

    let winner = &panel[0];
    println!("\n  VERDICT: {}", winner.name);
    println!("  → {}", winner.notes);
    println!("  → Bond threshold {:.1} nm — within Jupiter's UV window.", winner.uv_thresh);
    println!("  → Absorbs at ~{:.0} nm (blue/green) → reflects orange-red ✓", winner.abs_nm);
    println!("  → Sagan & Khare (1981) lab experiments confirmed.");

    let phosphorus = panel.iter().find(|c| c.name.contains("phosphorus")).unwrap();
    let sulfur = panel.iter().find(|c| c.name.contains("sulfur (S8)")).unwrap();
    check!(phosphorus.score() > sulfur.score(), "Phosphorus allotrope ranks above sulfur S₈");
    check!(phosphorus.uv_ok, "P-H bond photolysis is UV-accessible on Jupiter");
    check!((phosphorus.uv_thresh - 371.0).abs() < 10.0, "P-H UV threshold ≈ 371 nm");

    // ── D: Shrinkage and lifetime ─────────────────────────────────────────────

    println!("\n═══ D.  Shrinkage: History, Rate, Lifetime Prediction ═════════════════\n");
    println!("  Historical GRS long axis:");
    for &(yr, l, s) in GRS_HISTORY {
        println!("    {yr:.0}:  {l:.0} × {s:.0} km");
    }

    let (a_lin, b_lin) = fit_linear(GRS_HISTORY);
    let (a_exp, b_exp) = fit_exp(GRS_HISTORY);
    let current_year   = 2024.0_f64;

    let proj_lin_2024 = proj(a_lin, b_lin, current_year, false);
    let proj_exp_2024 = proj(a_exp, b_exp, current_year, true);

    println!("\n  Linear model:      size = {a_lin:.0} + {b_lin:.1} × year");
    println!("  (Rate: {:.0} km/year shrinkage)", b_lin.abs());
    println!("  → Predicted at 2024: {proj_lin_2024:.0} km  (observed: 20,000 km)");
    println!("\n  Exponential model: ln(size) = {a_exp:.4} + {b_exp:.6} × year");
    println!("  → Half-life: {:.0} years", -1.0 / b_exp * 2.0_f64.ln());
    println!("  → Predicted at 2024: {proj_exp_2024:.0} km  (observed: 20,000 km)");

    let t_crit_lin    = yr_at(a_lin, b_lin, critical_km, false);
    let t_crit_exp    = yr_at(a_exp, b_exp, critical_km, true);
    let size_2017     = 22_000.0_f64;
    let size_2024     = 20_000.0_f64;
    let recent_rate   = (size_2024 - size_2017) / (2024.0 - 2017.0);
    let t_crit_recent = current_year + (size_2024 - critical_km) / recent_rate.abs();

    println!("\n  Critical minimum size: 3 × L_D = {critical_km:.0} km");
    println!("  Current safety margin: {:.1}× critical", proj_lin_2024 / critical_km);
    println!("\n  LIFETIME PROJECTIONS:");
    println!("  Linear  model → critical in {t_crit_lin:.0}  ({:.0} years from now)", t_crit_lin - current_year);
    println!("  Exponential   → critical in {t_crit_exp:.0}  ({:.0} years from now)", t_crit_exp - current_year);
    println!("  Recent rate ({:.0} km/yr) → critical in {t_crit_recent:.0}  ({:.0} years)",
             recent_rate.abs(), t_crit_recent - current_year);
    println!("  Uncertainty: 2050–2220");

    check!(b_lin < 0.0, "Linear fit gives negative (shrinking) trend");
    check!(b_exp < 0.0, "Exponential fit gives negative (shrinking) trend");
    check!(proj_lin_2024 > 0.0 && proj_lin_2024 < 40_000.0, "Linear prediction is physical");
    check!(t_crit_lin > current_year, "GRS is not dead yet (linear model)");
    check!(t_crit_recent > current_year + 20.0, "Recent rate gives ≥ 20 years of life");

    // ── Summary ───────────────────────────────────────────────────────────────

    println!("\n═══ SUMMARY ════════════════════════════════════════════════════════════\n");
    println!("  α = {ALPHA:.6} → a₀ = {:.3} pm → r_H₂ = {:.3} pm → T_rot = {:.1} K",
             a0 * 1e12, r_h2 * 1e12, t_rot);
    println!("  γ at cloud tops ({T_CLOUD_TOP:.0} K) = {gamma_cloud:.3}  (transitional H₂/He mix)\n");
    println!("  ┌──────────────────────────────────────────────────────────────────────┐");
    println!("  │  α → H₂ → T_rot = {t_rot:.0} K → γ = {gamma_cloud:.3}                             │");
    println!("  │  → H = {:.1} km, N = {n_bv:.4e} rad/s, L_D = {:.0} km               │",
             h_scale / 1e3, l_d / 1e3);
    println!("  │                                                                      │");
    println!("  │  GRS: Ro = {ro:.3} << 1  →  geostrophic ✓                            │");
    println!("  │  GRS: PV = {pv:.3} < 0   →  anticyclonic, self-stable ✓              │");
    println!("  │  GRS: {:.0} km >> 3×L_D = {critical_km:.0} km  →  STABLE ✓                │",
             l_grs_m / 1e3);
    println!("  │                                                                      │");
    println!("  │  Color: P₄ from PH₃ photolysis at {:.0} nm                        │",
             winner.uv_thresh);
    println!("  │  Absorbs at {:.0} nm → appears orange-red ✓                        │",
             winner.abs_nm);
    println!("  │                                                                      │");
    println!("  │  Shrinkage: {:.0} km/yr → critical ~{t_crit_lin:.0}                          │",
             b_lin.abs());
    println!("  │  Recent: {:.0} km/yr → critical ~{t_crit_recent:.0}  (accelerating!)             │",
             recent_rate.abs());
    println!("  │  Uncertainty: 2050–2220 depending on erosion mechanism               │");
    println!("  └──────────────────────────────────────────────────────────────────────┘");

    // ── Assertion results ─────────────────────────────────────────────────────

    println!("\n─── Assertions: {}/{} pass ─────────────────────────────────────────────",
             assert_count - fails.len(), assert_count);
    for f in &fails { println!("  {f}"); }
    if fails.is_empty() { println!("  All {assert_count} assertions pass ✓"); }

    // ── Write all outputs ────────────────────────────────────────────────────

    let out_dir = env::var("GUTOE_GRS_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/jupiter_grs".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");

    println!("\n═══ Writing outputs to {out_dir}/ ═══════════════════════");

    write_text_report(
        &format!("{out_dir}/grs_findings.txt"),
        a0, r_h2, i_h2, b_rot, t_rot, gamma_cloud, gamma_deep,
        f_cor, h_scale, n_bv, l_d, u_grs, ro, pv,
        &panel,
        a_lin, b_lin, a_exp, b_exp, proj_lin_2024, proj_exp_2024,
        critical_km, t_crit_lin, t_crit_exp, t_crit_recent, recent_rate.abs(),
        assert_count, fails.len(),
    );
    println!("  wrote grs_findings.txt");

    write_csvs(&out_dir, &panel, a_lin, b_lin, a_exp, b_exp, critical_km);
    println!("  wrote grs_data_gamma.csv");
    println!("  wrote grs_data_shrinkage.csv");
    println!("  wrote grs_data_chromophores.csv");

    render_gamma_chart(&format!("{out_dir}/grs_chart_gamma.png"), t_rot);
    println!("  wrote grs_chart_gamma.png");

    render_shrinkage_chart(
        &format!("{out_dir}/grs_chart_shrinkage.png"),
        a_lin, b_lin, a_exp, b_exp, critical_km,
        t_crit_lin, t_crit_exp, t_crit_recent,
    );
    println!("  wrote grs_chart_shrinkage.png");

    render_chromophore_chart(&format!("{out_dir}/grs_chart_chromophores.png"), &panel);
    println!("  wrote grs_chart_chromophores.png");

    render_stability_chart(
        &format!("{out_dir}/grs_chart_stability.png"),
        a_lin, b_lin, a_exp, b_exp, critical_km,
        t_crit_lin, t_crit_exp, t_crit_recent,
    );
    println!("  wrote grs_chart_stability.png");

    assert!(fails.is_empty(), "{} assertion(s) failed", fails.len());
}
