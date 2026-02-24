//! Minimal stiff integrator scaffold (Rosenbrock-like, adaptive dt).

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntegratorConfig {
    pub dt_min: f64,
    pub dt_max: f64,
    pub tol: f64,
    pub gamma: f64,
}

impl Default for IntegratorConfig {
    fn default() -> Self {
        Self {
            dt_min: 1.0e-8,
            dt_max: 1.0,
            tol: 1.0e-6,
            gamma: 1.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StepResult {
    pub y_next: Vec<f64>,
    pub dt_used: f64,
    pub dt_suggested: f64,
    pub accepted: bool,
    pub err_norm: f64,
}

pub struct Rosenbrock1 {
    pub cfg: IntegratorConfig,
}

impl Rosenbrock1 {
    pub fn new(cfg: IntegratorConfig) -> Self {
        Self { cfg }
    }

    pub fn step<F, J>(&self, y: &[f64], dt: f64, f: F, jac: J) -> StepResult
    where
        F: Fn(&[f64]) -> Vec<f64>,
        J: Fn(&[f64]) -> Vec<Vec<f64>>,
    {
        let n = y.len();
        let dt = dt.clamp(self.cfg.dt_min, self.cfg.dt_max);
        let fy = f(y);
        let jy = jac(y);

        // Solve (I - gamma*dt*J) k = f(y)
        let mut a = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                a[i][j] = if i == j { 1.0 } else { 0.0 } - self.cfg.gamma * dt * jy[i][j];
            }
        }
        let k = solve_linear(&a, &fy).unwrap_or_else(|| vec![0.0; n]);
        let y_ros = y
            .iter()
            .zip(k.iter())
            .map(|(yi, ki)| yi + dt * ki)
            .collect::<Vec<_>>();

        // Embedded explicit Euler estimate for error control.
        let y_euler = y
            .iter()
            .zip(fy.iter())
            .map(|(yi, fi)| yi + dt * fi)
            .collect::<Vec<_>>();

        let err_norm = rms_diff(&y_ros, &y_euler);
        let accepted = err_norm <= self.cfg.tol || dt <= self.cfg.dt_min * 1.01;
        let scale = if err_norm > 0.0 {
            (self.cfg.tol / err_norm).powf(0.5)
        } else {
            2.0
        };
        let dt_suggested = (0.9 * dt * scale).clamp(self.cfg.dt_min, self.cfg.dt_max);

        StepResult {
            y_next: if accepted { y_ros } else { y.to_vec() },
            dt_used: dt,
            dt_suggested,
            accepted,
            err_norm,
        }
    }
}

fn rms_diff(a: &[f64], b: &[f64]) -> f64 {
    if a.is_empty() || a.len() != b.len() {
        return f64::INFINITY;
    }
    let m = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        / a.len() as f64;
    m.sqrt()
}

fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = a.len();
    if n == 0 || b.len() != n || a.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut aug = vec![vec![0.0_f64; n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j];
        }
        aug[i][n] = b[i];
    }

    for col in 0..n {
        let mut pivot = col;
        let mut best = aug[col][col].abs();
        for (r, row) in aug.iter().enumerate().skip(col + 1).take(n - (col + 1)) {
            let v = row[col].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best < 1e-14 {
            return None;
        }
        if pivot != col {
            aug.swap(pivot, col);
        }
        let diag = aug[col][col];
        for j in col..=n {
            aug[col][j] /= diag;
        }
        for r in 0..n {
            if r == col {
                continue;
            }
            let factor = aug[r][col];
            if factor.abs() < 1e-20 {
                continue;
            }
            for j in col..=n {
                aug[r][j] -= factor * aug[col][j];
            }
        }
    }

    Some((0..n).map(|i| aug[i][n]).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stiff_decay_step_is_stable() {
        // y' = -100 y ; true solution at dt=0.1 from y0=1 is exp(-10) ~= 4.5e-5.
        let rb = Rosenbrock1::new(IntegratorConfig::default());
        let mut y = vec![1.0];
        let mut dt = 0.1;
        let f = |y: &[f64]| vec![-100.0 * y[0]];
        let j = |_y: &[f64]| vec![vec![-100.0]];
        for _ in 0..8 {
            let s = rb.step(&y, dt, f, j);
            assert!(s.y_next[0].is_finite());
            if s.accepted {
                y = s.y_next;
                break;
            }
            dt = s.dt_suggested;
        }
        assert!(y[0] >= 0.0);
    }

    #[test]
    fn adaptive_rejects_when_error_is_large() {
        let rb = Rosenbrock1::new(IntegratorConfig {
            tol: 1e-12,
            ..IntegratorConfig::default()
        });
        let y0 = vec![1.0];
        let f = |y: &[f64]| vec![-10.0 * y[0]];
        let j = |_y: &[f64]| vec![vec![-10.0]];
        let s = rb.step(&y0, 1.0, f, j);
        assert!(!s.accepted);
        assert!(s.dt_suggested < s.dt_used);
    }
}
