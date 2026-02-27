use std::path::PathBuf;

use image::{imageops::FilterType, Rgb, RgbImage};

fn main() {
    let tile_dir = std::env::var("SKY_TILE_OUT")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/tmp/bh_renders/sky_tiles"));
    let out = std::env::var("SKY_STREAM_OUT")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/tmp/bh_renders/sky_stream_probe.png"));
    let tick_a = std::env::var("SKY_STREAM_A")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let tick_b = std::env::var("SKY_STREAM_B")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(1);
    let alpha = std::env::var("SKY_STREAM_ALPHA")
        .ok()
        .and_then(|s| s.parse::<f32>().ok())
        .unwrap_or(0.5)
        .clamp(0.0, 1.0);

    let pa = tile_dir.join(format!("tile_{tick_a:04}.png"));
    let pb = tile_dir.join(format!("tile_{tick_b:04}.png"));
    let a = image::open(&pa).expect("open tile A").to_rgb8();
    let mut b = image::open(&pb).expect("open tile B").to_rgb8();
    if b.dimensions() != a.dimensions() {
        b = image::imageops::resize(&b, a.width(), a.height(), FilterType::Triangle);
    }

    let mut out_img = RgbImage::new(a.width(), a.height());
    for y in 0..a.height() {
        for x in 0..a.width() {
            let p0 = a.get_pixel(x, y).0;
            let p1 = b.get_pixel(x, y).0;
            let lerp = |u: u8, v: u8| -> u8 {
                ((1.0 - alpha) * u as f32 + alpha * v as f32).round() as u8
            };
            out_img.put_pixel(
                x,
                y,
                Rgb([lerp(p0[0], p1[0]), lerp(p0[1], p1[1]), lerp(p0[2], p1[2])]),
            );
        }
    }
    out_img.save(&out).expect("save stream probe");
    println!("wrote {}", out.display());
}
