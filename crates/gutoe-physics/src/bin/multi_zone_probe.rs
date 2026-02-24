use gutoe_physics::{MultiZoneBurn, Species};

fn main() {
    let out = std::env::var("MULTI_ZONE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/multi_zone_probe.csv".to_string());
    let model = MultiZoneBurn::baseline();
    let mut zones = model.seed_zones(3);
    let steps = std::env::var("MULTI_ZONE_STEPS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(200);

    let mut csv = String::from("step,zone,h1,he4,power\n");
    for step in 0..steps {
        model.step(&mut zones, 1.0e5);
        for (zi, z) in zones.iter().enumerate() {
            csv.push_str(&format!(
                "{},{},{:.8e},{:.8e},{:.8e}\n",
                step,
                zi,
                z.get(Species::P1),
                z.get(Species::He4),
                z.thermal_power
            ));
        }
    }
    std::fs::write(&out, csv).expect("write multi-zone probe");
    println!("wrote {out}");
}
