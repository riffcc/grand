//! Catalog population pipeline: sampled star population -> reactor seeds.

use crate::{Species, ZoneState};

#[derive(Debug, Clone, PartialEq)]
pub struct StarSeed {
    pub id: u64,
    pub mass_solar: f64,
    pub age_gyr: f64,
    pub metallicity: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

pub fn synth_population(n: usize, seed: u64) -> Vec<StarSeed> {
    let mut rng = Lcg::new(seed);
    (0..n)
        .map(|i| {
            let mass_solar = sample_imf_powerlaw(&mut rng, 0.08, 60.0, 2.35);
            let age_gyr = rng.next_f64() * 13.0;
            let metallicity = 0.0001 + rng.next_f64() * 0.03;
            let r = 20_000.0 * rng.next_f64().sqrt();
            let theta = 2.0 * std::f64::consts::PI * rng.next_f64();
            let z = (rng.next_f64() - 0.5) * 600.0;
            StarSeed {
                id: i as u64,
                mass_solar,
                age_gyr,
                metallicity,
                x: r * theta.cos(),
                y: r * theta.sin(),
                z,
            }
        })
        .collect()
}

pub fn seed_to_reactor_state(s: &StarSeed) -> ZoneState {
    let z = s.metallicity.clamp(0.0, 0.2);
    let y = (0.24 + 2.0 * z).clamp(0.20, 0.60); // helium mass fraction proxy
    let x = (1.0 - y - z).max(0.01); // hydrogen
    let cno = (z * 0.7).max(0.0);
    let mut st = ZoneState::solar_like_seed();
    st.abund.insert(Species::P1, x);
    st.abund.insert(Species::He4, y);
    st.abund.insert(Species::C12, cno * 0.4);
    st.abund.insert(Species::N14, cno * 0.3);
    st.abund.insert(Species::O15, cno * 0.3);
    let sum = st.abund.values().sum::<f64>();
    if sum > 1e-12 {
        for v in st.abund.values_mut() {
            *v /= sum;
        }
    }
    st
}

fn sample_imf_powerlaw(rng: &mut Lcg, m_min: f64, m_max: f64, alpha: f64) -> f64 {
    let u = rng.next_f64().clamp(1e-12, 1.0 - 1e-12);
    let p = 1.0 - alpha;
    let a = m_min.powf(p);
    let b = m_max.powf(p);
    (a + u * (b - a)).powf(1.0 / p)
}

#[derive(Debug, Clone)]
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self { state: seed | 1 }
    }
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.state
    }
    fn next_f64(&mut self) -> f64 {
        let x = self.next_u64() >> 11;
        (x as f64) / ((1u64 << 53) as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn population_is_deterministic_for_seed() {
        let a = synth_population(8, 42);
        let b = synth_population(8, 42);
        assert_eq!(a, b);
    }

    #[test]
    fn reactor_seed_abundances_normalize() {
        let s = StarSeed {
            id: 1,
            mass_solar: 1.0,
            age_gyr: 4.5,
            metallicity: 0.013,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let st = seed_to_reactor_state(&s);
        let sum = st.abund.values().sum::<f64>();
        assert!((sum - 1.0).abs() < 1e-9);
    }
}
