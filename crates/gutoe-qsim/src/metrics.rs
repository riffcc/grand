/// Per-round snapshot of the content-addressed simulation.
pub struct RoundMetrics {
    pub round: usize,
    pub gate_name: String,
    pub actors_targeted: usize,
    pub unique_states_before: usize,
    pub unique_states_after: usize,
    pub new_states_created: usize,
    pub compression_factor: f64, // actors / unique_states
}

impl std::fmt::Display for RoundMetrics {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "round {}: {} -> {} unique ({} new), compression {:.1}x",
            self.round,
            self.gate_name,
            self.unique_states_after,
            self.new_states_created,
            self.compression_factor,
        )
    }
}

/// Final results of a content-addressed simulation experiment.
pub struct ExperimentResult {
    pub total_actors: usize,
    pub total_rounds: usize,
    pub final_unique_states: usize,
    pub peak_compression: f64,
    pub growth_curve: Vec<usize>, // unique_states[round]
    pub growth_exponent: f64,     // polynomial fit: S ~ t^b
}

impl ExperimentResult {
    /// Fit `unique_states ~ t^b` via log-log linear regression.
    ///
    /// Uses least-squares on `ln(unique_states) = b * ln(t) + c`.
    /// Skips round 0 (ln(0) undefined) and entries where unique_states <= 1.
    pub fn fit_growth_exponent(growth_curve: &[usize]) -> f64 {
        let mut sum_x = 0.0_f64;
        let mut sum_y = 0.0_f64;
        let mut sum_xx = 0.0_f64;
        let mut sum_xy = 0.0_f64;
        let mut n = 0u64;

        for (t, &s) in growth_curve.iter().enumerate() {
            if t == 0 || s <= 1 {
                continue;
            }
            let x = (t as f64).ln();
            let y = (s as f64).ln();
            sum_x += x;
            sum_y += y;
            sum_xx += x * x;
            sum_xy += x * y;
            n += 1;
        }

        if n < 2 {
            return 0.0;
        }

        let nf = n as f64;
        let denom = nf * sum_xx - sum_x * sum_x;
        if denom.abs() < 1e-15 {
            return 0.0;
        }

        (nf * sum_xy - sum_x * sum_y) / denom
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fit_linear_growth() {
        // S = t^1.0 → exponent should be ~1.0
        let curve: Vec<usize> = (0..=100).map(|t| t.max(1)).collect();
        let b = ExperimentResult::fit_growth_exponent(&curve);
        assert!(
            (b - 1.0).abs() < 0.1,
            "linear growth should give exponent ~1.0, got {b:.3}"
        );
    }

    #[test]
    fn fit_sqrt_growth() {
        // S = t^0.5 → exponent should be ~0.5
        let curve: Vec<usize> = (0..=100)
            .map(|t| ((t as f64).sqrt() * 10.0).max(1.0) as usize)
            .collect();
        let b = ExperimentResult::fit_growth_exponent(&curve);
        assert!(
            (b - 0.5).abs() < 0.15,
            "sqrt growth should give exponent ~0.5, got {b:.3}"
        );
    }

    #[test]
    fn fit_constant() {
        // S = constant → exponent should be ~0
        let curve = vec![5usize; 101];
        let b = ExperimentResult::fit_growth_exponent(&curve);
        assert!(
            b.abs() < 0.1,
            "constant should give exponent ~0, got {b:.3}"
        );
    }

    #[test]
    fn round_metrics_display() {
        let m = RoundMetrics {
            round: 3,
            gate_name: "Z3Cycle".to_string(),
            actors_targeted: 4000,
            unique_states_before: 2,
            unique_states_after: 3,
            new_states_created: 1,
            compression_factor: 1333.3,
        };
        let s = format!("{m}");
        assert!(s.contains("round 3"), "display should show round number");
        assert!(s.contains("Z3Cycle"), "display should show gate name");
        assert!(s.contains("1333.3x"), "display should show compression");
    }
}
