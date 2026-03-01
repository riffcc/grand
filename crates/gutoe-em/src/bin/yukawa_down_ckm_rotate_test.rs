//! Step-1 test: rotate down-sector by structural CKM dagger and recheck closure.
//!
//! Requested check:
//!   d_mass = V_CKM * d_weak  =>  d_weak = V_CKM^† * d_mass
//!
//! We run this in amplitude space (sqrt masses), then map back to masses by |amp|^2.
//! No optimizer, no free fitted parameter.

use gutoe_em::alpha::z3_extract_params;
use gutoe_em::ckm_from_clifford;
use num_complex::Complex64;
use serde::Serialize;
use std::fs;

const D: f64 = 4.67;
const S: f64 = 93.0;
const B: f64 = 4180.0;
const C: f64 = 1270.0;
const T: f64 = 172_760.0;

#[derive(Debug, Clone, Serialize)]
struct FitSummary {
    masses_mev: [f64; 3],
    m_scale: f64,
    s: f64,
    s2: f64,
    delta_deg: f64,
    ratio_ms_over_md: f64,
    ratio_mb_over_ms: f64,
    down_closure_rms_log: f64,
    cross_ratio_mc_over_ms: f64,
    cross_ratio_mt_over_mb: f64,
}

#[derive(Debug, Clone, Serialize)]
struct Report {
    ckm_structural: CkmDump,
    down_mass_basis: FitSummary,
    down_weak_basis_rotated: FitSummary,
    target_ms_over_md: f64,
    target_mb_over_ms: f64,
    target_mc_over_ms: f64,
    target_mt_over_mb: f64,
    closure_improved: bool,
    summary: String,
}

#[derive(Debug, Clone, Serialize)]
struct CkmDump {
    s12: f64,
    s23: f64,
    s13: f64,
    delta_deg: f64,
    theta12_deg: f64,
    theta23_deg: f64,
    theta13_deg: f64,
}

fn ckm_matrix_pdg(
    s12: f64,
    s23: f64,
    s13: f64,
    delta_rad: f64,
) -> [[Complex64; 3]; 3] {
    let c12 = (1.0 - s12 * s12).sqrt();
    let c23 = (1.0 - s23 * s23).sqrt();
    let c13 = (1.0 - s13 * s13).sqrt();

    let e_pos = Complex64::from_polar(1.0, delta_rad);
    let e_neg = Complex64::from_polar(1.0, -delta_rad);

    [
        [
            Complex64::new(c12 * c13, 0.0),
            Complex64::new(s12 * c13, 0.0),
            Complex64::new(s13, 0.0) * e_neg,
        ],
        [
            Complex64::new(-s12 * c23, 0.0) - Complex64::new(c12 * s23 * s13, 0.0) * e_pos,
            Complex64::new(c12 * c23, 0.0) - Complex64::new(s12 * s23 * s13, 0.0) * e_pos,
            Complex64::new(s23 * c13, 0.0),
        ],
        [
            Complex64::new(s12 * s23, 0.0) - Complex64::new(c12 * c23 * s13, 0.0) * e_pos,
            Complex64::new(-c12 * s23, 0.0) - Complex64::new(s12 * c23 * s13, 0.0) * e_pos,
            Complex64::new(c23 * c13, 0.0),
        ],
    ]
}

fn apply_vdag_to_vec(v: [[Complex64; 3]; 3], x: [Complex64; 3]) -> [Complex64; 3] {
    // y_i = sum_j conj(V_{j i}) x_j
    let mut y = [Complex64::new(0.0, 0.0); 3];
    for i in 0..3 {
        let mut s = Complex64::new(0.0, 0.0);
        for j in 0..3 {
            s += v[j][i].conj() * x[j];
        }
        y[i] = s;
    }
    y
}

fn rms_log(a: f64, b: f64, ta: f64, tb: f64) -> f64 {
    let x = (a / ta).ln();
    let y = (b / tb).ln();
    ((x * x + y * y) / 2.0).sqrt()
}

fn summarize(masses: [f64; 3], target_msmd: f64, target_mbms: f64) -> FitSummary {
    let (m_scale, s, delta_rad) = z3_extract_params(masses);
    let ms_md = masses[1] / masses[0];
    let mb_ms = masses[2] / masses[1];
    FitSummary {
        masses_mev: masses,
        m_scale,
        s,
        s2: s * s,
        delta_deg: delta_rad.to_degrees(),
        ratio_ms_over_md: ms_md,
        ratio_mb_over_ms: mb_ms,
        down_closure_rms_log: rms_log(ms_md, mb_ms, target_msmd, target_mbms),
        cross_ratio_mc_over_ms: C / masses[1],
        cross_ratio_mt_over_mb: T / masses[2],
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_YUKAWA_DOWN_CKM_ROTATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/yukawa_down_ckm_rotate_test.txt");
    let json_path = format!("{out_dir}/yukawa_down_ckm_rotate_test.json");

    let ckm = ckm_from_clifford();
    let v = ckm_matrix_pdg(ckm.s12, ckm.s23, ckm.s13, ckm.delta_rad);

    // Structural targets from existing quark-ratio lane.
    let target_msmd = 19.0;
    let target_mbms = (8.0 / 3.0) * 19.0 * (67.0 / 66.0);
    let target_mcms = (13.0 / 21.0) * 19.0 * (67.0 / 66.0);
    let target_mtmb = (13.0 / 6.0) * 19.0 * (67.0 / 66.0);

    let down_mass_basis = [D, S, B];
    let a_mass = [
        Complex64::new(down_mass_basis[0].sqrt(), 0.0),
        Complex64::new(down_mass_basis[1].sqrt(), 0.0),
        Complex64::new(down_mass_basis[2].sqrt(), 0.0),
    ];
    let a_weak = apply_vdag_to_vec(v, a_mass);
    let down_weak_basis = [a_weak[0].norm_sqr(), a_weak[1].norm_sqr(), a_weak[2].norm_sqr()];

    let mass_summary = summarize(down_mass_basis, target_msmd, target_mbms);
    let weak_summary = summarize(down_weak_basis, target_msmd, target_mbms);

    let closure_improved = weak_summary.down_closure_rms_log < mass_summary.down_closure_rms_log;

    let summary = format!(
        "down closure RMS-log: mass-basis {:.6} -> weak-basis {:.6}; \
cross mc/ms: {:.6} -> {:.6}; cross mt/mb: {:.6} -> {:.6}; improved={}.",
        mass_summary.down_closure_rms_log,
        weak_summary.down_closure_rms_log,
        mass_summary.cross_ratio_mc_over_ms,
        weak_summary.cross_ratio_mc_over_ms,
        mass_summary.cross_ratio_mt_over_mb,
        weak_summary.cross_ratio_mt_over_mb,
        closure_improved
    );

    let report = Report {
        ckm_structural: CkmDump {
            s12: ckm.s12,
            s23: ckm.s23,
            s13: ckm.s13,
            delta_deg: ckm.delta_deg,
            theta12_deg: ckm.theta12_deg,
            theta23_deg: ckm.theta23_deg,
            theta13_deg: ckm.theta13_deg,
        },
        down_mass_basis: mass_summary,
        down_weak_basis_rotated: weak_summary,
        target_ms_over_md: target_msmd,
        target_mb_over_ms: target_mbms,
        target_mc_over_ms: target_mcms,
        target_mt_over_mb: target_mtmb,
        closure_improved,
        summary,
    };

    let mut txt = String::new();
    txt.push_str("[yukawa_down_ckm_rotate_test]\n");
    txt.push_str(&format!(
        "theta12_deg={:.9} theta23_deg={:.9} theta13_deg={:.9} delta_deg={:.9}\n\n",
        report.ckm_structural.theta12_deg,
        report.ckm_structural.theta23_deg,
        report.ckm_structural.theta13_deg,
        report.ckm_structural.delta_deg
    ));
    txt.push_str("[down_mass_basis]\n");
    txt.push_str(&format!(
        "masses={:?}\nms/md={:.9} mb/ms={:.9} closure_rms_log={:.9}\nmc/ms={:.9} mt/mb={:.9}\n\n",
        report.down_mass_basis.masses_mev,
        report.down_mass_basis.ratio_ms_over_md,
        report.down_mass_basis.ratio_mb_over_ms,
        report.down_mass_basis.down_closure_rms_log,
        report.down_mass_basis.cross_ratio_mc_over_ms,
        report.down_mass_basis.cross_ratio_mt_over_mb,
    ));
    txt.push_str("[down_weak_basis_rotated]\n");
    txt.push_str(&format!(
        "masses={:?}\nms/md={:.9} mb/ms={:.9} closure_rms_log={:.9}\nmc/ms={:.9} mt/mb={:.9}\n\n",
        report.down_weak_basis_rotated.masses_mev,
        report.down_weak_basis_rotated.ratio_ms_over_md,
        report.down_weak_basis_rotated.ratio_mb_over_ms,
        report.down_weak_basis_rotated.down_closure_rms_log,
        report.down_weak_basis_rotated.cross_ratio_mc_over_ms,
        report.down_weak_basis_rotated.cross_ratio_mt_over_mb,
    ));
    txt.push_str(&format!(
        "[targets]\nms/md={:.9} mb/ms={:.9} mc/ms={:.9} mt/mb={:.9}\n\n",
        report.target_ms_over_md,
        report.target_mb_over_ms,
        report.target_mc_over_ms,
        report.target_mt_over_mb
    ));
    txt.push_str(&format!(
        "closure_improved={}\nsummary={}\n",
        report.closure_improved, report.summary
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&report).expect("serialize"),
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!("{summary}", summary = report.summary);
}
