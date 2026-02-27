use gutoe_physics::synth_population;

fn main() {
    let out = std::env::var("STELLAR_SANITY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/stellar_sanity_hr.csv".to_string());
    let n = std::env::var("STELLAR_SANITY_N")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(50_000);
    let stars = synth_population(n, 2026);

    let mut csv =
        String::from("id,mass,age_gyr,metallicity,temperature_k,luminosity_lsun,log_teff,log_l\n");
    for s in &stars {
        let t_eff = effective_temp_proxy(s.mass_solar, s.age_gyr, s.metallicity);
        let l = luminosity_proxy(s.mass_solar, s.age_gyr);
        csv.push_str(&format!(
            "{},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n",
            s.id,
            s.mass_solar,
            s.age_gyr,
            s.metallicity,
            t_eff,
            l,
            t_eff.log10(),
            l.log10()
        ));
    }
    std::fs::write(&out, csv).expect("write sanity csv");
    println!("wrote {out}");
}

fn luminosity_proxy(mass_solar: f64, age_gyr: f64) -> f64 {
    let m = mass_solar.clamp(0.08, 60.0);
    let main = m.powf(3.5);
    let age_term = (1.0 - 0.03 * age_gyr).clamp(0.2, 1.5);
    (main * age_term).max(1e-6)
}

fn effective_temp_proxy(mass_solar: f64, age_gyr: f64, metallicity: f64) -> f64 {
    let m = mass_solar.clamp(0.08, 60.0);
    let z = metallicity.clamp(0.0001, 0.03);
    let z_term = (0.02 / z).powf(0.06).clamp(0.8, 1.2);
    let age_term = (1.0 - 0.01 * age_gyr).clamp(0.7, 1.2);
    (5800.0 * m.powf(0.55) * z_term * age_term).clamp(2000.0, 60_000.0)
}
