use gutoe_physics::{IntegratorConfig, Rosenbrock1};

fn main() {
    let out = std::env::var("STIFF_PROBE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/stiff_probe.csv".to_string());
    let rb = Rosenbrock1::new(IntegratorConfig::default());
    let mut y = vec![1.0_f64];
    let mut t = 0.0_f64;
    let mut dt = 0.05_f64;
    let t_end = 1.0_f64;

    let f = |yy: &[f64]| vec![-100.0 * yy[0]];
    let j = |_yy: &[f64]| vec![vec![-100.0]];

    let mut csv = String::from("t,y,dt,accepted,err\n");
    while t < t_end {
        let s = rb.step(&y, dt, f, j);
        csv.push_str(&format!(
            "{:.8},{:.12e},{:.6e},{},{}\n",
            t, y[0], s.dt_used, s.accepted, s.err_norm
        ));
        if s.accepted {
            y = s.y_next;
            t += s.dt_used;
        }
        dt = s.dt_suggested;
    }
    std::fs::write(&out, csv).expect("write stiff probe");
    println!("wrote {out}");
}
