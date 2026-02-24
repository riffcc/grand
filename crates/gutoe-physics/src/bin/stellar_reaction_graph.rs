use gutoe_physics::ReactionNetwork;

fn main() {
    let network = ReactionNetwork::baseline_p1();
    let out = std::env::var("STELLAR_REACTION_GRAPH_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/stellar_reaction_graph.json".to_string());

    let mut json = String::from("{\n  \"reactions\": [\n");
    for (i, r) in network.reactions.iter().enumerate() {
        if i > 0 {
            json.push_str(",\n");
        }
        json.push_str("    {\n");
        json.push_str(&format!("      \"id\": \"{}\",\n", r.id));
        json.push_str(&format!("      \"channel\": \"{}\",\n", r.channel));
        json.push_str("      \"reactants\": [");
        for (j, s) in r.reactants.iter().enumerate() {
            if j > 0 {
                json.push_str(", ");
            }
            json.push_str(&format!(
                "{{\"species\": \"{:?}\", \"coeff\": {}}}",
                s.species, s.coeff
            ));
        }
        json.push_str("],\n      \"products\": [");
        for (j, s) in r.products.iter().enumerate() {
            if j > 0 {
                json.push_str(", ");
            }
            json.push_str(&format!(
                "{{\"species\": \"{:?}\", \"coeff\": {}}}",
                s.species, s.coeff
            ));
        }
        json.push_str("],\n");
        json.push_str(&format!(
            "      \"branching_weight\": {:.6},\n",
            r.branching_weight
        ));
        json.push_str(&format!("      \"q_mev\": {:.6}\n", r.q_mev));
        json.push_str("    }");
    }
    json.push_str("\n  ],\n");
    let stoich = network.stoichiometric_matrix();
    json.push_str("  \"stoichiometric_shape\": {");
    json.push_str(&format!(
        "\"rows\": {}, \"cols\": {}",
        stoich.len(),
        stoich.first().map_or(0, |r| r.len())
    ));
    json.push_str("},\n");
    json.push_str(&format!(
        "  \"all_conserved\": {}\n",
        if network.all_conserved() { "true" } else { "false" }
    ));
    json.push_str("}\n");

    std::fs::write(&out, json).expect("write reaction graph json");
    println!("wrote {out}");
}
