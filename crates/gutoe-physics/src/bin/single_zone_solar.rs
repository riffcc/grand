use gutoe_physics::{SingleZoneBurn, Species, ZoneState};

fn main() {
    let out = std::env::var("SINGLE_ZONE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/single_zone_solar.csv".to_string());
    let mut s = ZoneState::solar_like_seed();
    let burn = SingleZoneBurn::baseline();
    let steps = std::env::var("SINGLE_ZONE_STEPS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(200);

    let mut csv = String::from("step,h1,he4,c12,power\n");
    for i in 0..steps {
        burn.step(&mut s, 0.02, 1.0e6);
        csv.push_str(&format!(
            "{},{:.8e},{:.8e},{:.8e},{:.8e}\n",
            i,
            s.get(Species::P1),
            s.get(Species::He4),
            s.get(Species::C12),
            s.thermal_power
        ));
    }

    std::fs::write(&out, csv).expect("write single zone csv");
    println!("wrote {out}");
}
