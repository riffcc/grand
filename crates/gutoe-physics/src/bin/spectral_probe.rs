use std::collections::HashMap;

use gutoe_physics::{synthesize_spectrum, Species};

fn main() {
    let out = std::env::var("SPECTRAL_PROBE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/spectral_probe.csv".to_string());
    let mut abund = HashMap::new();
    abund.insert(Species::P1, 0.70);
    abund.insert(Species::He4, 0.28);
    abund.insert(Species::C12, 0.01);
    abund.insert(Species::N14, 0.005);
    abund.insert(Species::O15, 0.005);
    let s = synthesize_spectrum(&abund, 5800.0, 256);

    let mut csv = String::from("type,wavelength_nm,value,name\n");
    for p in &s.continuum {
        csv.push_str(&format!(
            "continuum,{:.6},{:.12e},-\n",
            p.wavelength_nm, p.intensity
        ));
    }
    for l in &s.lines {
        csv.push_str(&format!(
            "line,{:.6},{:.12e},{}\n",
            l.wavelength_nm, l.strength, l.name
        ));
    }
    std::fs::write(&out, csv).expect("write spectral probe csv");
    println!("wrote {out}");
}
