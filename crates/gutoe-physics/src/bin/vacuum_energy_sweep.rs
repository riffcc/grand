use gutoe_physics::constants::{C, HBAR, LAMBDA_COSMOLOGICAL_OBSERVED};
use gutoe_physics::dark_sector::{
    vacuum_energy_density_from_lambda, vacuum_energy_density_structural,
};

#[derive(Clone, Copy, Debug)]
struct CasimirRow {
    gap_m: f64,
    pressure_pa: f64,
    energy_density_j_m3: f64,
    one_shot_work_j: f64,
    net_cycle_j: f64,
}

fn casimir_pressure(gap_m: f64) -> f64 {
    // Ideal parallel-plate Casimir pressure magnitude:
    // P = π² ħ c / (240 a⁴)
    std::f64::consts::PI.powi(2) * HBAR * C / (240.0 * gap_m.powi(4))
}

fn casimir_energy_density(gap_m: f64) -> f64 {
    // Energy density magnitude between plates (idealized):
    // |u| = π² ħ c / (720 a⁴)
    std::f64::consts::PI.powi(2) * HBAR * C / (720.0 * gap_m.powi(4))
}

fn casimir_one_shot_work(area_m2: f64, a_min: f64, a_max: f64) -> f64 {
    // Ideal reversible extracted work when moving from a_max -> a_min:
    // E(a) = -K A / a^3, K = π² ħ c / 720
    // W_out = K A (1/a_min^3 - 1/a_max^3)
    let k = std::f64::consts::PI.powi(2) * HBAR * C / 720.0;
    k * area_m2 * (1.0 / a_min.powi(3) - 1.0 / a_max.powi(3)).max(0.0)
}

fn cycle_net_energy(one_shot_work_j: f64, harvest_eff: f64, reset_eff: f64) -> f64 {
    // Close stroke: harvest usable energy = η_h * W.
    // Re-open stroke: actuator energy cost >= W / η_reset.
    // Net per cycle:
    //   E_net = η_h * W - W / η_reset
    one_shot_work_j * harvest_eff - one_shot_work_j / reset_eff
}

fn main() {
    let out_dir = std::env::var("GUTOE_VACUUM_SWEEP_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/vacuum_energy_sweep".to_string());
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    let rho_struct = vacuum_energy_density_structural();
    let rho_obs = vacuum_energy_density_from_lambda(LAMBDA_COSMOLOGICAL_OBSERVED).unwrap_or(0.0);

    let u_struct = rho_struct * C * C;
    let u_obs = rho_obs * C * C;

    // Sweep settings: intentionally optimistic lab geometry.
    let area_m2 = 1.0e-4; // 1 cm^2
    let gap_min_m = 1.0e-9; // 1 nm
    let gap_max_m = 1.0e-8; // 10 nm stroke upper anchor
    let harvest_eff = 0.90;
    let reset_eff = 0.90;

    let gaps_m = [
        1.0e-9, 2.0e-9, 5.0e-9, 1.0e-8, 2.0e-8, 5.0e-8, 1.0e-7, 2.0e-7, 5.0e-7, 1.0e-6, 2.0e-6,
        5.0e-6, 1.0e-5,
    ];

    let one_shot = casimir_one_shot_work(area_m2, gap_min_m, gap_max_m);
    let net_cycle = cycle_net_energy(one_shot, harvest_eff, reset_eff);

    let mut rows = Vec::new();
    for &a in &gaps_m {
        let p = casimir_pressure(a);
        let u = casimir_energy_density(a);
        rows.push(CasimirRow {
            gap_m: a,
            pressure_pa: p,
            energy_density_j_m3: u,
            one_shot_work_j: one_shot,
            net_cycle_j: net_cycle,
        });
    }

    let mut csv = String::from(
        "gap_m,pressure_pa,casimir_energy_density_j_m3,one_shot_work_j_for_1cm2_10to1nm,net_cycle_j_eta_h0p9_eta_reset0p9\n",
    );
    for r in &rows {
        csv.push_str(&format!(
            "{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n",
            r.gap_m, r.pressure_pa, r.energy_density_j_m3, r.one_shot_work_j, r.net_cycle_j
        ));
    }

    let mut txt = String::new();
    txt.push_str("[vacuum_energy_density]\n");
    txt.push_str(&format!("rho_lambda_struct_kg_m3 = {:.12e}\n", rho_struct));
    txt.push_str(&format!("rho_lambda_obs_kg_m3 = {:.12e}\n", rho_obs));
    txt.push_str(&format!("u_lambda_struct_j_m3 = {:.12e}\n", u_struct));
    txt.push_str(&format!("u_lambda_obs_j_m3 = {:.12e}\n\n", u_obs));

    txt.push_str("[casimir_cycle_assumptions]\n");
    txt.push_str(&format!("plate_area_m2 = {:.6e}\n", area_m2));
    txt.push_str(&format!("stroke_from_m = {:.6e}\n", gap_max_m));
    txt.push_str(&format!("stroke_to_m = {:.6e}\n", gap_min_m));
    txt.push_str(&format!("harvest_efficiency = {:.3}\n", harvest_eff));
    txt.push_str(&format!("reset_efficiency = {:.3}\n", reset_eff));
    txt.push_str(&format!("one_shot_work_j = {:.12e}\n", one_shot));
    txt.push_str(&format!("net_cycle_j = {:.12e}\n", net_cycle));
    txt.push_str("net_cycle_sign = ");
    txt.push_str(if net_cycle > 0.0 {
        "positive\n\n"
    } else if net_cycle < 0.0 {
        "negative\n\n"
    } else {
        "zero\n\n"
    });

    txt.push_str("[verdict]\n");
    txt.push_str("uniform_vacuum_harvestable = false (no gradient/work cycle in equilibrium)\n");
    txt.push_str(
        "casimir_one_shot_exists = true (geometry-dependent differential vacuum energy)\n",
    );
    txt.push_str("casimir_cycle_free_energy = false (reset work cancels/exceeds harvest)\n");

    let csv_path = format!("{out_dir}/vacuum_energy_sweep.csv");
    let txt_path = format!("{out_dir}/vacuum_energy_sweep.txt");
    std::fs::write(&csv_path, csv).expect("write csv");
    std::fs::write(&txt_path, txt).expect("write txt");

    println!("wrote {}", csv_path);
    println!("wrote {}", txt_path);
    println!(
        "u_lambda_struct={:.6e} J/m^3, one_shot={:.6e} J, net_cycle={:.6e} J",
        u_struct, one_shot, net_cycle
    );
}
