//! Multi-zone stellar burn model with simple radial mixing hooks.

use std::collections::HashMap;

use crate::{SingleZoneBurn, Species, ZoneState};

#[derive(Debug, Clone)]
pub struct MultiZoneConfig {
    pub diffusion_coeff: f64,
    pub zone_temperatures_t9: Vec<f64>,
    /// Relative reaction-rate scale per zone (e.g., density^2 weighting).
    pub zone_rate_scales: Vec<f64>,
}

impl Default for MultiZoneConfig {
    fn default() -> Self {
        Self {
            diffusion_coeff: 1.0e-10,
            zone_temperatures_t9: vec![0.04, 0.03, 0.02],
            zone_rate_scales: vec![1.0, 1.0, 1.0],
        }
    }
}

#[derive(Debug, Clone)]
pub struct MultiZoneBurn {
    pub core: SingleZoneBurn,
    pub cfg: MultiZoneConfig,
}

impl MultiZoneBurn {
    pub fn baseline() -> Self {
        Self {
            core: SingleZoneBurn::baseline(),
            cfg: MultiZoneConfig::default(),
        }
    }

    pub fn seed_zones(&self, n: usize) -> Vec<ZoneState> {
        (0..n).map(|_| ZoneState::solar_like_seed()).collect()
    }

    pub fn step(&self, zones: &mut [ZoneState], dt: f64) {
        if zones.is_empty() {
            return;
        }
        for (i, z) in zones.iter_mut().enumerate() {
            let t9 = self
                .cfg
                .zone_temperatures_t9
                .get(i)
                .copied()
                .or_else(|| self.cfg.zone_temperatures_t9.last().copied())
                .unwrap_or(0.02);
            let rate_scale = self
                .cfg
                .zone_rate_scales
                .get(i)
                .copied()
                .or_else(|| self.cfg.zone_rate_scales.last().copied())
                .unwrap_or(1.0)
                .max(0.0);
            self.core.step(z, t9, dt * rate_scale);
        }
        self.mix_adjacent(zones, dt);
    }

    fn mix_adjacent(&self, zones: &mut [ZoneState], dt: f64) {
        let k = (self.cfg.diffusion_coeff * dt).max(0.0);
        if k <= 0.0 {
            return;
        }
        for i in 0..zones.len().saturating_sub(1) {
            let left = zones[i].abund.clone();
            let right = zones[i + 1].abund.clone();
            let species_union = left
                .keys()
                .chain(right.keys())
                .copied()
                .collect::<std::collections::HashSet<_>>();
            for s in species_union {
                let l = *left.get(&s).unwrap_or(&0.0);
                let r = *right.get(&s).unwrap_or(&0.0);
                let mut flux = k * (r - l);
                flux = flux.max(-l).min(r);
                set_abund(&mut zones[i].abund, s, (l + flux).max(0.0));
                set_abund(&mut zones[i + 1].abund, s, (r - flux).max(0.0));
            }
        }
        for z in zones {
            renorm(&mut z.abund);
        }
    }
}

fn set_abund(map: &mut HashMap<Species, f64>, s: Species, v: f64) {
    if v > 0.0 {
        map.insert(s, v);
    } else {
        map.remove(&s);
    }
}

fn renorm(abund: &mut HashMap<Species, f64>) {
    let sum = abund.values().sum::<f64>();
    if sum > 1e-12 {
        for v in abund.values_mut() {
            *v /= sum;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multi_zone_runs_and_preserves_norm() {
        let m = MultiZoneBurn::baseline();
        let mut z = m.seed_zones(3);
        for _ in 0..50 {
            m.step(&mut z, 1.0e5);
        }
        for zone in &z {
            let s = zone.abund.values().sum::<f64>();
            assert!((s - 1.0).abs() < 1e-8);
        }
    }
}
