use std::f64::consts::PI;
use std::path::PathBuf;

use gutoe_physics::synth_population;
use image::{Rgb, RgbImage};

fn main() {
    let out_dir = std::env::var("SKY_TILE_OUT")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/tmp/bh_renders/sky_tiles"));
    let width = std::env::var("SKY_TILE_W")
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(2048);
    let height = std::env::var("SKY_TILE_H")
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(1024);
    let ticks = std::env::var("SKY_TILE_TICKS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(8);
    let stars_per_tick = std::env::var("SKY_TILE_STARS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(100_000);
    std::fs::create_dir_all(&out_dir).expect("create sky tile out dir");

    let mut manifest = String::from("tick,path,width,height,stars\n");
    for tick in 0..ticks {
        let stars = synth_population(stars_per_tick, 1337 + tick as u64);
        let mut img = RgbImage::from_pixel(width, height, Rgb([1, 2, 8]));
        for s in &stars {
            let r = (s.x * s.x + s.y * s.y + s.z * s.z).sqrt().max(1e-9);
            let dx = s.x / r;
            let dy = s.y / r;
            let dz = s.z / r;
            let lon = dy.atan2(dx);
            let lat = dz.asin();
            let uf = ((lon + PI) / (2.0 * PI)).clamp(0.0, 1.0);
            let vf = ((PI / 2.0 - lat) / PI).clamp(0.0, 1.0);
            let x = (uf * (width as f64 - 1.0)).round() as u32;
            let y = (vf * (height as f64 - 1.0)).round() as u32;

            let mass_boost = (s.mass_solar.log10().max(-2.0) + 2.0) / 3.5;
            let age_fade = (1.0 - s.age_gyr / 13.0).clamp(0.2, 1.0);
            let temp = (2500.0 + 12_000.0 * mass_boost * age_fade).clamp(2000.0, 20_000.0);
            let (cr, cg, cb) = temp_to_rgb(temp);
            let pix = img.get_pixel_mut(x, y);
            let old = pix.0;
            pix.0 = [
                old[0].saturating_add(cr),
                old[1].saturating_add(cg),
                old[2].saturating_add(cb),
            ];
        }
        let path = out_dir.join(format!("tile_{tick:04}.png"));
        img.save(&path).expect("write tile image");
        manifest.push_str(&format!(
            "{},{},{},{},{}\n",
            tick,
            path.display(),
            width,
            height,
            stars.len()
        ));
    }
    let manifest_path = out_dir.join("manifest.csv");
    std::fs::write(&manifest_path, manifest).expect("write manifest");
    println!("wrote {}", manifest_path.display());
}

fn temp_to_rgb(temp_k: f64) -> (u8, u8, u8) {
    let t = (temp_k / 100.0).clamp(10.0, 400.0);
    let r = if t <= 66.0 {
        255.0
    } else {
        329.698727446 * (t - 60.0).powf(-0.1332047592)
    };
    let g = if t <= 66.0 {
        99.4708025861 * t.ln() - 161.1195681661
    } else {
        288.1221695283 * (t - 60.0).powf(-0.0755148492)
    };
    let b = if t >= 66.0 {
        255.0
    } else if t <= 19.0 {
        0.0
    } else {
        138.5177312231 * (t - 10.0).ln() - 305.0447927307
    };
    (
        r.clamp(0.0, 255.0) as u8,
        g.clamp(0.0, 255.0) as u8,
        b.clamp(0.0, 255.0) as u8,
    )
}
