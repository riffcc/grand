//! Structural ratio probe for the CTC origin-gain lane.
//!
//! Goal: test whether the required product `(eta * infra_gain)` can be
//! expressed by existing Cl(1,3) structural counts/rationals, with exact
//! rational arithmetic (no float fitting).

use serde_json::json;
use std::cmp::Ordering;
use std::fs;
use std::path::PathBuf;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Rat {
    n: i128,
    d: i128,
}

impl Rat {
    fn new(n: i128, d: i128) -> Self {
        assert!(d != 0);
        let mut n = n;
        let mut d = d;
        if d < 0 {
            n = -n;
            d = -d;
        }
        let g = gcd_i128(n.abs(), d.abs());
        Self { n: n / g, d: d / g }
    }

    fn mul(self, other: Self) -> Self {
        Self::new(self.n * other.n, self.d * other.d)
    }

    fn div(self, other: Self) -> Self {
        Self::new(self.n * other.d, self.d * other.n)
    }

    fn sub(self, other: Self) -> Self {
        Self::new(self.n * other.d - other.n * self.d, self.d * other.d)
    }

    fn abs(self) -> Self {
        Self::new(self.n.abs(), self.d)
    }

    fn to_f64(self) -> f64 {
        self.n as f64 / self.d as f64
    }

    fn as_string(self) -> String {
        format!("{}/{}", self.n, self.d)
    }
}

fn gcd_i128(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a.max(1)
}

#[derive(Clone)]
struct Candidate {
    name: &'static str,
    value: Rat,
    note: &'static str,
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_GAIN_RATIO_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_gain_ratio_structural_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Structural constants from Cl(1,3) rails.
    let z3 = Rat::new(3, 1);
    let grade1 = Rat::new(4, 1);
    let grade2 = Rat::new(6, 1);
    let basis = Rat::new(16, 1);
    let ew_sum = Rat::new(10, 1); // grade1 + grade2
    let z3_fixed_grade1 = Rat::new(1, 1);
    let sin2w = Rat::new(3, 13);
    let lambda_h = Rat::new(13, 100);

    // In the closure lane, with branching=3 and merge=3/16:
    // 1 = (3)*(3/16)*(eta*infra) => eta*infra = 16/9.
    let branching = z3;
    let merge = Rat::new(3, 16);
    let required_eta_infra = Rat::new(1, 1).div(branching.mul(merge));
    let target = Rat::new(16, 9);

    let candidates = vec![
        Candidate {
            name: "basis_16 / z3^2",
            value: basis.div(z3.mul(z3)),
            note: "clean: microstate count over generation-orbit square",
        },
        Candidate {
            name: "grade1^2 / z3^2",
            value: grade1.mul(grade1).div(z3.mul(z3)),
            note: "clean: spacetime-vector square over generation-orbit square",
        },
        Candidate {
            name: "basis_16 / (ew_sum - z3_fixed_grade1)",
            value: basis.div(Rat::new(ew_sum.n - z3_fixed_grade1.n, ew_sum.d)),
            note: "clean: 16 over (4+6-1)=9",
        },
        Candidate {
            name: "sin2w / lambda_h",
            value: sin2w.div(lambda_h),
            note: "near: coupled EW/Higgs ratio",
        },
        Candidate {
            name: "1 + cos2w",
            value: Rat::new(1, 1).add(Rat::new(10, 13)),
            note: "near: electroweak complement sum",
        },
        Candidate {
            name: "24/13 (su3* sin2w)",
            value: Rat::new(24, 13),
            note: "near-high: color-weighted weak angle",
        },
    ];

    let mut rows = Vec::new();
    for c in candidates {
        let delta = c.value.sub(target);
        rows.push(json!({
            "name": c.name,
            "value_q": c.value.as_string(),
            "value_f64": c.value.to_f64(),
            "delta_q": delta.as_string(),
            "delta_abs_f64": delta.abs().to_f64(),
            "exact_target_match": c.value == target,
            "note": c.note
        }));
    }
    rows.sort_by(|a, b| {
        let da = a["delta_abs_f64"].as_f64().unwrap_or(f64::INFINITY);
        let db = b["delta_abs_f64"].as_f64().unwrap_or(f64::INFINITY);
        da.partial_cmp(&db).unwrap_or(Ordering::Equal)
    });

    let payload = json!({
        "structural_constants": {
            "z3": z3.as_string(),
            "grade1": grade1.as_string(),
            "grade2": grade2.as_string(),
            "basis_16": basis.as_string(),
            "ew_sum": ew_sum.as_string(),
            "z3_fixed_grade1": z3_fixed_grade1.as_string(),
            "sin2w": sin2w.as_string(),
            "lambda_h": lambda_h.as_string(),
            "merge_void": merge.as_string()
        },
        "closure_constraint": {
            "equation": "1 = branching * merge * (eta*infra)",
            "branching": branching.as_string(),
            "merge": merge.as_string(),
            "required_eta_infra_q": required_eta_infra.as_string(),
            "required_eta_infra_f64": required_eta_infra.to_f64(),
            "target_16_over_9_q": target.as_string()
        },
        "candidates_sorted_by_abs_delta": rows
    });

    let txt_path = out.join("ctc_gain_ratio_structural_probe.txt");
    let json_path = out.join("ctc_gain_ratio_structural_probe.json");

    let mut txt = String::new();
    txt.push_str("[ctc_gain_ratio_structural_probe]\n");
    txt.push_str(&format!(
        "required_eta_infra = {} ({:.12})\n",
        required_eta_infra.as_string(),
        required_eta_infra.to_f64()
    ));
    txt.push_str("target = 16/9\n\n");
    txt.push_str("[candidates]\n");
    for row in payload["candidates_sorted_by_abs_delta"]
        .as_array()
        .expect("array")
    {
        txt.push_str(&format!(
            "{} = {} ({:.12}), exact={}, |delta|={:.6e}\n",
            row["name"].as_str().unwrap_or(""),
            row["value_q"].as_str().unwrap_or(""),
            row["value_f64"].as_f64().unwrap_or(0.0),
            row["exact_target_match"].as_bool().unwrap_or(false),
            row["delta_abs_f64"].as_f64().unwrap_or(0.0),
        ));
    }

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}

trait RatAdd {
    fn add(self, other: Self) -> Self;
}

impl RatAdd for Rat {
    fn add(self, other: Self) -> Self {
        Rat::new(self.n * other.d + other.n * self.d, self.d * other.d)
    }
}
