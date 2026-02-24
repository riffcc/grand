use gutoe_physics::{RateEngine, ReactionNetwork};

fn main() {
    let network = ReactionNetwork::baseline_p1();
    let rates = RateEngine::baseline_p1();
    let out = std::env::var("STELLAR_RATE_PROBE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/stellar_rate_probe.csv".to_string());
    let t9 = std::env::var("STELLAR_RATE_T9")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.02);

    let mut csv = String::from("reaction_id,channel,t9,rate\n");
    for r in &network.reactions {
        let rate = rates.rate_for(r.id, t9).unwrap_or(0.0);
        csv.push_str(&format!("{},{},{:.6e},{:.6e}\n", r.id, r.channel, t9, rate));
    }
    std::fs::write(&out, csv).expect("write rate probe csv");
    println!("wrote {out}");
}
