use gutoe_physics::{compute_power_budget, RateEngine, ReactionNetwork};

fn main() {
    let out = std::env::var("STELLAR_POWER_PROBE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/stellar_power_probe.csv".to_string());
    let t9 = std::env::var("STELLAR_POWER_T9")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.02);

    let net = ReactionNetwork::baseline_p1();
    let rates = RateEngine::baseline_p1();
    let b = compute_power_budget(&net, &rates, t9);

    let mut csv =
        String::from("reaction_id,channel,rate,q_mev,nu_frac,thermal_power,neutrino_power\n");
    for r in &b.rows {
        csv.push_str(&format!(
            "{},{},{:.6e},{:.6e},{:.4},{:.6e},{:.6e}\n",
            r.reaction_id,
            r.channel,
            r.rate,
            r.q_mev,
            r.neutrino_loss_fraction,
            r.thermal_power,
            r.neutrino_power
        ));
    }
    csv.push_str(&format!(
        "TOTAL,TOTAL,0,0,0,{:.6e},{:.6e}\n",
        b.total_thermal_power, b.total_neutrino_power
    ));

    std::fs::write(&out, csv).expect("write stellar power probe csv");
    println!("wrote {out}");
}
