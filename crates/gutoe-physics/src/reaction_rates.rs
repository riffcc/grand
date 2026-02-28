//! Temperature-dependent reaction-rate tables with interpolation.
//!
//! This is the P1 baseline rate engine scaffold for stellar burning.
//! Rates are represented in table form and interpolated in log-log space.

use crate::stellar_reactions::ReactionNetwork;

#[derive(Debug, Clone, PartialEq)]
pub struct RatePoint {
    pub t9: f64,
    pub rate: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RateTable {
    pub reaction_id: &'static str,
    pub points: Vec<RatePoint>,
}

impl RateTable {
    pub fn interpolate_loglog(&self, t9: f64) -> Option<f64> {
        if self.points.len() < 2 || t9 <= 0.0 {
            return None;
        }
        let mut pts = self.points.clone();
        pts.sort_by(|a, b| a.t9.partial_cmp(&b.t9).unwrap_or(std::cmp::Ordering::Equal));
        let first = pts.first()?;
        let last = pts.last()?;
        if t9 <= first.t9 {
            return Some(first.rate);
        }
        if t9 >= last.t9 {
            return Some(last.rate);
        }
        for w in pts.windows(2) {
            let a = &w[0];
            let b = &w[1];
            if (a.t9..=b.t9).contains(&t9) {
                let xa = a.t9.ln();
                let xb = b.t9.ln();
                let ya = a.rate.ln();
                let yb = b.rate.ln();
                let x = t9.ln();
                let u = (x - xa) / (xb - xa);
                return Some((ya + u * (yb - ya)).exp());
            }
        }
        None
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RateEngine {
    pub tables: Vec<RateTable>,
}

impl RateEngine {
    pub fn baseline_p1() -> Self {
        // Baseline anchors for first pass. Units are normalized proxy rates.
        // Real REACLIB/NACRE ingestion can replace these tables in P1.1.
        let tables = vec![
            table("pp_1", &[(0.01, 1.0e-22), (0.02, 5.0e-21), (0.05, 8.0e-19)]),
            table("pp_2", &[(0.01, 2.0e-18), (0.02, 1.0e-16), (0.05, 2.0e-14)]),
            table("pp_3", &[(0.01, 5.0e-15), (0.02, 1.0e-12), (0.05, 2.0e-10)]),
            table(
                "cno_1",
                &[(0.01, 1.0e-25), (0.02, 1.0e-21), (0.05, 8.0e-16)],
            ),
            table(
                "cno_2",
                &[(0.01, 5.0e-11), (0.02, 5.5e-11), (0.05, 6.0e-11)],
            ),
            table(
                "cno_3",
                &[(0.01, 2.0e-24), (0.02, 6.0e-21), (0.05, 2.0e-15)],
            ),
            table(
                "cno_4",
                &[(0.01, 1.0e-28), (0.02, 5.0e-24), (0.05, 6.0e-18)],
            ),
            table("cno_5", &[(0.01, 2.0e-9), (0.02, 2.2e-9), (0.05, 2.5e-9)]),
            table(
                "cno_6",
                &[(0.01, 1.0e-23), (0.02, 3.0e-20), (0.05, 9.0e-15)],
            ),
            table(
                "triple_alpha",
                &[(0.05, 5.0e-31), (0.1, 3.0e-24), (0.2, 4.0e-16)],
            ),
            // Proxy Pop-II lithium burn anchors around the 2.5e6 K threshold.
            // Units are normalized per-year effective rates for envelope depletion.
            table(
                "li7_burn",
                &[
                    (0.0015, 2.0e-10),
                    (0.0020, 6.0e-10),
                    (0.0025, 1.7e-9),
                    (0.0030, 4.5e-9),
                    (0.0035, 1.1e-8),
                    (0.0040, 2.5e-8),
                ],
            ),
        ];
        Self { tables }
    }

    pub fn rate_for(&self, reaction_id: &str, t9: f64) -> Option<f64> {
        self.tables
            .iter()
            .find(|t| t.reaction_id == reaction_id)
            .and_then(|t| t.interpolate_loglog(t9))
    }

    pub fn covers_network(&self, network: &ReactionNetwork) -> bool {
        network
            .reactions
            .iter()
            .all(|r| self.tables.iter().any(|t| t.reaction_id == r.id))
    }
}

fn table(id: &'static str, pts: &[(f64, f64)]) -> RateTable {
    RateTable {
        reaction_id: id,
        points: pts
            .iter()
            .map(|(t9, rate)| RatePoint {
                t9: *t9,
                rate: *rate,
            })
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stellar_reactions::ReactionNetwork;

    #[test]
    fn interpolation_is_monotone_for_monotone_table() {
        let t = table("x", &[(0.01, 1.0e-5), (0.02, 1.0e-3), (0.05, 1.0e-1)]);
        let r1 = t.interpolate_loglog(0.015).expect("r1");
        let r2 = t.interpolate_loglog(0.03).expect("r2");
        assert!(r2 > r1);
    }

    #[test]
    fn baseline_rate_engine_covers_reaction_graph() {
        let net = ReactionNetwork::baseline_p1();
        let eng = RateEngine::baseline_p1();
        assert!(eng.covers_network(&net));
        assert!(eng.rate_for("pp_1", 0.02).unwrap_or(0.0) > 0.0);
        assert!(eng.rate_for("triple_alpha", 0.1).unwrap_or(0.0) > 0.0);
    }
}
