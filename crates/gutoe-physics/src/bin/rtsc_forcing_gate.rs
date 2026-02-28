use gutoe_physics::StandardModelDynamicsMap;
use std::env;
use std::fs;

fn env_f64(name: &str, default: f64) -> f64 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let out_dir = env::var("GUTOE_RTSC_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_forcing_gate".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");

    let m = StandardModelDynamicsMap::from_clifford_z3();

    // Forced lattice/filling primitives (shared with Lean theorem lane).
    let coordination = 2_u32 * m.magnetic_triplet_card; // should be 6
    let forced_lattice = if coordination == 6 {
        "simple_cubic"
    } else {
        "rejected"
    };
    let filling_multiplicity = m.z3_order; // should be 3

    // Pairing-sign gate from shared finite counts:
    // darkFraction = 5/16, repulsion = 1/(2*|grade2|) = 1/12.
    let dark_fraction = env_f64("GUTOE_RTSC_DARK_FRACTION", 5.0 / 16.0);
    let repulsion = env_f64("GUTOE_RTSC_REPULSION", 1.0 / 12.0);
    let pairing_kernel = dark_fraction - repulsion;
    let pairing_sign = if pairing_kernel > 0.0 {
        "attractive"
    } else if pairing_kernel < 0.0 {
        "repulsive"
    } else {
        "neutral"
    };

    // Structural Tc proxy aligned with Lean `tcStructuralQ = 300*(1 + kernel)`.
    let tc_proxy_k = 300.0 * (1.0 + pairing_kernel);

    let admissible = forced_lattice == "simple_cubic"
        && filling_multiplicity == 3
        && pairing_kernel > 0.0
        && tc_proxy_k >= 300.0;
    let verdict = if admissible { "ADMISSIBLE" } else { "NO_GO" };

    let mut txt = String::new();
    txt.push_str("[rtsc_forcing_gate]\n");
    txt.push_str(&format!("verdict = {}\n", verdict));
    txt.push_str("mode = forced_gate_single_shot\n\n");

    txt.push_str("[forced_lattice]\n");
    txt.push_str(&format!("coordination = {}\n", coordination));
    txt.push_str(&format!("family = {}\n\n", forced_lattice));

    txt.push_str("[forced_filling]\n");
    txt.push_str(&format!("z3_order = {}\n", filling_multiplicity));
    txt.push_str(&format!(
        "triplet_required = {}\n\n",
        filling_multiplicity == 3
    ));

    txt.push_str("[pairing_kernel]\n");
    txt.push_str(&format!("dark_fraction = {:.12e}\n", dark_fraction));
    txt.push_str(&format!("repulsion = {:.12e}\n", repulsion));
    txt.push_str(&format!("kernel = {:.12e}\n", pairing_kernel));
    txt.push_str(&format!("sign = {}\n\n", pairing_sign));

    txt.push_str("[tc_gate]\n");
    txt.push_str(&format!("tc_proxy_k = {:.9}\n", tc_proxy_k));
    txt.push_str("threshold_k = 300.000000000\n");
    txt.push_str(&format!("passes = {}\n\n", tc_proxy_k >= 300.0));

    txt.push_str("[proof_links]\n");
    txt.push_str("lean_module = Gutoe.RTSCAdmissibility\n");
    txt.push_str("theorems = forced_lattice_family_simple_cubic, forced_filling_triplet, pairing_kernel_attractive, tc_structural_ge_300, rtsc_gate_admissible\n");

    let json = format!(
        concat!(
            "{{\n",
            "  \"verdict\": \"{}\",\n",
            "  \"forced_lattice\": {{\"coordination\": {}, \"family\": \"{}\"}},\n",
            "  \"forced_filling\": {{\"z3_order\": {}, \"triplet_required\": {}}},\n",
            "  \"pairing_kernel\": {{\"dark_fraction\": {:.12e}, \"repulsion\": {:.12e}, \"kernel\": {:.12e}, \"sign\": \"{}\"}},\n",
            "  \"tc_gate\": {{\"tc_proxy_k\": {:.12e}, \"threshold_k\": 3.000000000000e2, \"passes\": {}}},\n",
            "  \"proof_links\": {{\"lean_module\": \"Gutoe.RTSCAdmissibility\"}}\n",
            "}}\n"
        ),
        verdict,
        coordination,
        forced_lattice,
        filling_multiplicity,
        if filling_multiplicity == 3 { "true" } else { "false" },
        dark_fraction,
        repulsion,
        pairing_kernel,
        pairing_sign,
        tc_proxy_k,
        if tc_proxy_k >= 300.0 { "true" } else { "false" },
    );

    let txt_path = format!("{out_dir}/rtsc_forcing_gate.txt");
    let json_path = format!("{out_dir}/rtsc_forcing_gate.json");
    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, json).expect("write json");
    println!("wrote {}", txt_path);
    println!("wrote {}", json_path);
}

