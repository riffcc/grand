use gutoe_em::{neutrino_hierarchy_prediction, neutrino_texture_eigenvalues};
use std::fs::{self, File};
use std::io::Write;

fn main() {
    let out_dir = std::env::var("GUTOE_NEUTRINO_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/neutrino_hierarchy_report".to_string());
    let _ = fs::create_dir_all(&out_dir);

    let m = neutrino_texture_eigenvalues();
    let dm21 = m[1] - m[0];
    let dm31 = m[2] - m[0];
    let hierarchy = neutrino_hierarchy_prediction();

    let txt_path = format!("{out_dir}/neutrino_hierarchy_report.txt");
    let json_path = format!("{out_dir}/neutrino_hierarchy_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[neutrino_texture_eigenvalues]").expect("write");
    writeln!(txt, "m1 = {:.12e}", m[0]).expect("write");
    writeln!(txt, "m2 = {:.12e}", m[1]).expect("write");
    writeln!(txt, "m3 = {:.12e}", m[2]).expect("write");
    writeln!(txt, "delta_m21 = {:.12e}", dm21).expect("write");
    writeln!(txt, "delta_m31 = {:.12e}", dm31).expect("write");
    writeln!(txt, "hierarchy_prediction = {hierarchy}").expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"m1\": {:.12e},\n  \"m2\": {:.12e},\n  \"m3\": {:.12e},\n  \"delta_m21\": {:.12e},\n  \"delta_m31\": {:.12e},\n  \"hierarchy_prediction\": \"{}\"\n}}",
        m[0], m[1], m[2], dm21, dm31, hierarchy
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "hierarchy_prediction={} (m1={:.3e}, m2={:.3e}, m3={:.3e})",
        hierarchy, m[0], m[1], m[2]
    );
}
