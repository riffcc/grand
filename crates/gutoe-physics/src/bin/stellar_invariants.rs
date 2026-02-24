use gutoe_physics::{MultiZoneBurn, SingleZoneBurn, Species, ZoneState};

fn main() {
    let out = std::env::var("STELLAR_INVARIANTS_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/stellar_invariants.csv".to_string());
    let steps = std::env::var("STELLAR_INVARIANTS_STEPS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(200);

    let burn = SingleZoneBurn::baseline();
    let mut single = ZoneState::solar_like_seed();
    let mut multi = MultiZoneBurn::baseline().seed_zones(3);
    let multi_burn = MultiZoneBurn::baseline();

    let mut csv = String::from("step,single_sum,single_min,single_h,single_power,multi_min,multi_max_dev\n");
    for step in 0..steps {
        burn.step(&mut single, 0.02, 1.0e6);
        multi_burn.step(&mut multi, 1.0e5);

        let single_vals: Vec<f64> = single.abund.values().copied().collect();
        let single_sum: f64 = single_vals.iter().sum();
        let single_min = single_vals.iter().copied().fold(f64::INFINITY, f64::min);
        let single_h = single.get(Species::P1);

        let mut multi_min = f64::INFINITY;
        let mut multi_max_dev: f64 = 0.0;
        for z in &multi {
            let vals: Vec<f64> = z.abund.values().copied().collect();
            let s = vals.iter().sum::<f64>();
            let mn = vals.iter().copied().fold(f64::INFINITY, f64::min);
            multi_min = multi_min.min(mn);
            multi_max_dev = multi_max_dev.max((1.0 - s).abs());
        }

        csv.push_str(&format!(
            "{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}\n",
            step, single_sum, single_min, single_h, single.thermal_power, multi_min, multi_max_dev
        ));
    }

    std::fs::write(&out, csv).expect("write invariants csv");
    println!("wrote {out}");
}
