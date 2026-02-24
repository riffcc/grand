use gutoe_physics::{seed_to_reactor_state, synth_population, Species};

fn main() {
    let out = std::env::var("STAR_CATALOG_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/star_catalog_probe.csv".to_string());
    let n = std::env::var("STAR_CATALOG_N")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(1024);
    let seed = std::env::var("STAR_CATALOG_SEED")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(1337);
    let stars = synth_population(n, seed);

    let mut csv = String::from("id,mass_solar,age_gyr,metallicity,h1,he4,x,y,z\n");
    for s in &stars {
        let st = seed_to_reactor_state(s);
        csv.push_str(&format!(
            "{},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n",
            s.id,
            s.mass_solar,
            s.age_gyr,
            s.metallicity,
            st.get(Species::P1),
            st.get(Species::He4),
            s.x,
            s.y,
            s.z
        ));
    }
    std::fs::write(&out, csv).expect("write star catalog probe");
    println!("wrote {out}");
}
