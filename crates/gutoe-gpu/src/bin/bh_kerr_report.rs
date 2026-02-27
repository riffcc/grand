//! Kerr report utility for GRAND-159/164 groundwork.
//!
//! Prints physically checkable Kerr quantities for a given spin.

use std::f64::consts::FRAC_PI_2;

use gutoe_gpu::kerr::KerrMetric;

fn main() {
    // CLI:
    //   bh_kerr_report              -> defaults (r_s=2, a*=0.9)
    //   bh_kerr_report 2.0 0.94     -> custom r_s and dimensionless spin a*
    let args: Vec<String> = std::env::args().collect();
    let r_s = args
        .get(1)
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(2.0);
    let a_star = args
        .get(2)
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.9);

    let Some(k) = KerrMetric::new(r_s, a_star) else {
        eprintln!("invalid parameters: require r_s > 0 and |a*| <= 1");
        std::process::exit(2);
    };

    let m = k.mass();
    let a = k.a();
    let (r_plus, r_minus) = k.horizons();
    let r_ergo_eq = k.ergosphere_radius(FRAC_PI_2);
    let r_ph_pro = k.equatorial_photon_orbit_radius(true);
    let r_ph_ret = k.equatorial_photon_orbit_radius(false);
    let omega_h = k.horizon_angular_velocity();
    let omega_eq_2m = k.frame_dragging_omega(2.0 * m, FRAC_PI_2);

    println!("Kerr baseline report (GR reference for GUTOE-Kerr extension)");
    println!("inputs:");
    println!("  r_s    = {r_s:.6}");
    println!("  M      = {m:.6}");
    println!("  a*     = {:.6}", k.a_over_m());
    println!("  a      = {a:.6}");
    println!("derived:");
    println!("  r_+    = {r_plus:.6}   (outer horizon)");
    println!("  r_-    = {r_minus:.6}   (inner horizon)");
    println!("  r_erg,eq = {r_ergo_eq:.6}   (equatorial ergosphere)");
    println!("  r_ph,pro = {r_ph_pro:.6}    (equatorial prograde photon orbit)");
    println!("  r_ph,ret = {r_ph_ret:.6}    (equatorial retrograde photon orbit)");
    println!("  Ω_H      = {omega_h:.6}");
    println!("  ω(r=2M, θ=π/2) = {omega_eq_2m:.6}");
}
