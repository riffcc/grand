/*!
 * Weinberg-angle M_Z target report.
 *
 * Quantifies the correction needed to move from the structural value
 * sin²(theta_W)=3/13 to the observed M_Z-scale value.
 */

use gutoe_em::weak::sin2_weinberg;
use std::fs::{self, File};
use std::io::Write;

const SIN2_THETA_W_MZ_OBSERVED: f64 = 0.23122;

fn main() {
    let structural = sin2_weinberg();
    let observed = SIN2_THETA_W_MZ_OBSERVED;
    let delta = observed - structural;
    let rel_pct = (delta / observed).abs() * 100.0;

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let csv_path = format!("{out_dir}/weinberg_mz_report.csv");
    let txt_path = format!("{out_dir}/weinberg_mz_summary.txt");

    let mut csv = File::create(&csv_path).expect("create weinberg csv");
    writeln!(csv, "quantity,value").expect("write csv header");
    writeln!(csv, "sin2_theta_w_structural,{structural:.12}").expect("write structural");
    writeln!(csv, "sin2_theta_w_observed_mz,{observed:.12}").expect("write observed");
    writeln!(csv, "delta_observed_minus_structural,{delta:.12}").expect("write delta");
    writeln!(csv, "relative_percent_vs_observed,{rel_pct:.9}").expect("write rel pct");

    let mut txt = File::create(&txt_path).expect("create weinberg summary");
    writeln!(txt, "sin²θ_W structural (3/13): {structural:.12}").expect("summary structural");
    writeln!(txt, "sin²θ_W observed at M_Z:   {observed:.12}").expect("summary observed");
    writeln!(txt, "delta (obs - structural):  {delta:.12}").expect("summary delta");
    writeln!(txt, "relative deviation:        {rel_pct:.9}%").expect("summary rel");
    writeln!(
        txt,
        "Interpretation: this is the RG/loop correction budget that GRAND-61 must explain."
    )
    .expect("summary interpretation");

    println!("wrote {csv_path}");
    println!("wrote {txt_path}");
    println!("sin²θ_W structural = {structural:.9}");
    println!("sin²θ_W observed   = {observed:.9}");
    println!("delta             = {delta:.9}");
    println!("relative error    = {rel_pct:.6}%");
}
