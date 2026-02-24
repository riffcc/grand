use gutoe_physics::StandardModelDynamicsMap;

fn main() {
    let m = StandardModelDynamicsMap::from_clifford_z3();
    let out = std::env::var("SM_MAP_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/sm_dynamics_map.json".to_string());
    let json = format!(
        concat!(
            "{{\n",
            "  \"clifford_dim\": {},\n",
            "  \"z3_order\": {},\n",
            "  \"magnetic_triplet_card\": {},\n",
            "  \"generations\": {},\n",
            "  \"sin2_theta_w\": {:.15},\n",
            "  \"cos2_theta_w\": {:.15},\n",
            "  \"mz_over_mw_sq\": {:.15},\n",
            "  \"alpha_leading_order\": {:.15},\n",
            "  \"lambda_qg\": {:.15},\n",
            "  \"beta0\": {:.15},\n",
            "  \"su3_generators\": {},\n",
            "  \"su2_generators\": {},\n",
            "  \"u1_generators\": {},\n",
            "  \"total_gauge_generators\": {},\n",
            "  \"valid_internal_constraints\": {}\n",
            "}}\n"
        ),
        m.clifford_dim,
        m.z3_order,
        m.magnetic_triplet_card,
        m.generations,
        m.sin2_theta_w,
        m.cos2_theta_w,
        m.mz_over_mw_sq,
        m.alpha_leading_order,
        m.lambda_qg,
        m.beta0,
        m.su3_generators,
        m.su2_generators,
        m.u1_generators,
        m.total_gauge_generators,
        if m.validate_internal_constraints() {
            "true"
        } else {
            "false"
        }
    );
    std::fs::write(&out, json).expect("write map json");
    println!("wrote {out}");
}
