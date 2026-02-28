//! Proton spin decomposition report from shared Cl(1,3) primitives.
//!
//! Structural channel split:
//!   quark  = 5 channels (dark-sector finite split count)
//!   gluon  = 8 channels (SU(3) generators)
//!   orbital = 4 channels (grade-1 state count)
//! total = 17 = 16 + 1 (Clifford basis plus identity lane).

use gutoe_physics::constants::{
    CLIFFORD_STATE_COUNT_STRUCTURAL, DARK_STATE_COUNT_STRUCTURAL, GRADE1_STATE_COUNT_STRUCTURAL,
};
use gutoe_physics::StandardModelDynamicsMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let out_dir = std::env::var("GUTOE_PROTON_SPIN_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/proton_spin_report".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let sm = StandardModelDynamicsMap::from_clifford_z3();

    let n_quark_channels = DARK_STATE_COUNT_STRUCTURAL; // 5
    let n_gluon_channels = sm.su3_generators as f64; // 8
    let n_orbital_channels = GRADE1_STATE_COUNT_STRUCTURAL; // 4
    let channel_total = n_quark_channels + n_gluon_channels + n_orbital_channels; // 17

    let clifford_plus_identity = CLIFFORD_STATE_COUNT_STRUCTURAL + 1.0; // 16 + 1
    let channel_total_matches_clifford_plus_identity =
        (channel_total - clifford_plus_identity).abs() < 1.0e-12;

    let quark_fraction = n_quark_channels / channel_total;
    let gluon_fraction = n_gluon_channels / channel_total;
    let orbital_fraction = n_orbital_channels / channel_total;
    let fraction_sum = quark_fraction + gluon_fraction + orbital_fraction;

    let j_total = 0.5_f64;
    let j_quark = j_total * quark_fraction;
    let j_gluon = j_total * gluon_fraction;
    let j_orbital = j_total * orbital_fraction;
    let j_sum = j_quark + j_gluon + j_orbital;

    // Broad phenomenology windows for a first-pass spin-crisis lane.
    let quark_window = (0.25_f64, 0.35_f64);
    let gluon_window = (0.35_f64, 0.55_f64);
    let orbital_window = (0.15_f64, 0.35_f64);

    let quark_in_window = quark_fraction >= quark_window.0 && quark_fraction <= quark_window.1;
    let gluon_in_window = gluon_fraction >= gluon_window.0 && gluon_fraction <= gluon_window.1;
    let orbital_in_window =
        orbital_fraction >= orbital_window.0 && orbital_fraction <= orbital_window.1;

    let txt_path = out.join("proton_spin_report.txt");
    let json_path = out.join("proton_spin_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[channels]").expect("write");
    writeln!(txt, "quark_channels = {:.0}", n_quark_channels).expect("write");
    writeln!(txt, "gluon_channels = {:.0}", n_gluon_channels).expect("write");
    writeln!(txt, "orbital_channels = {:.0}", n_orbital_channels).expect("write");
    writeln!(txt, "total_channels = {:.0}", channel_total).expect("write");
    writeln!(txt, "clifford_plus_identity = {:.0}", clifford_plus_identity).expect("write");
    writeln!(
        txt,
        "total_matches_clifford_plus_identity = {}",
        channel_total_matches_clifford_plus_identity
    )
    .expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[fractions]").expect("write");
    writeln!(txt, "quark_fraction = {:.12}", quark_fraction).expect("write");
    writeln!(txt, "gluon_fraction = {:.12}", gluon_fraction).expect("write");
    writeln!(txt, "orbital_fraction = {:.12}", orbital_fraction).expect("write");
    writeln!(txt, "fraction_sum = {:.12}", fraction_sum).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[angular_momentum_hbar]").expect("write");
    writeln!(txt, "J_total = {:.12}", j_total).expect("write");
    writeln!(txt, "J_quark = {:.12}", j_quark).expect("write");
    writeln!(txt, "J_gluon = {:.12}", j_gluon).expect("write");
    writeln!(txt, "J_orbital = {:.12}", j_orbital).expect("write");
    writeln!(txt, "J_sum = {:.12}", j_sum).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[benchmark_windows]").expect("write");
    writeln!(
        txt,
        "quark_window = [{:.3}, {:.3}] | pass = {}",
        quark_window.0, quark_window.1, quark_in_window
    )
    .expect("write");
    writeln!(
        txt,
        "gluon_window = [{:.3}, {:.3}] | pass = {}",
        gluon_window.0, gluon_window.1, gluon_in_window
    )
    .expect("write");
    writeln!(
        txt,
        "orbital_window = [{:.3}, {:.3}] | pass = {}",
        orbital_window.0, orbital_window.1, orbital_in_window
    )
    .expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"channels\": {{\"quark\": {:.0}, \"gluon\": {:.0}, \"orbital\": {:.0}, \"total\": {:.0}, \"clifford_plus_identity\": {:.0}, \"total_matches_clifford_plus_identity\": {}}},\n  \"fractions\": {{\"quark\": {:.12}, \"gluon\": {:.12}, \"orbital\": {:.12}, \"sum\": {:.12}}},\n  \"angular_momentum_hbar\": {{\"total\": {:.12}, \"quark\": {:.12}, \"gluon\": {:.12}, \"orbital\": {:.12}, \"sum\": {:.12}}},\n  \"benchmark_windows\": {{\"quark\": {{\"min\": {:.3}, \"max\": {:.3}, \"pass\": {}}}, \"gluon\": {{\"min\": {:.3}, \"max\": {:.3}, \"pass\": {}}}, \"orbital\": {{\"min\": {:.3}, \"max\": {:.3}, \"pass\": {}}}}}\n}}",
        n_quark_channels,
        n_gluon_channels,
        n_orbital_channels,
        channel_total,
        clifford_plus_identity,
        channel_total_matches_clifford_plus_identity,
        quark_fraction,
        gluon_fraction,
        orbital_fraction,
        fraction_sum,
        j_total,
        j_quark,
        j_gluon,
        j_orbital,
        j_sum,
        quark_window.0,
        quark_window.1,
        quark_in_window,
        gluon_window.0,
        gluon_window.1,
        gluon_in_window,
        orbital_window.0,
        orbital_window.1,
        orbital_in_window
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "spin fractions q={:.4}, g={:.4}, L={:.4} (sum {:.4})",
        quark_fraction, gluon_fraction, orbital_fraction, fraction_sum
    );
}

