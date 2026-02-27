//! Runtime Standard-Model emulator scaffold.
//!
//! Deterministic coupled gauge+matter evolution loop driven by the
//! theorem-linked [`StandardModelDynamicsMap`](crate::dynamics_map::StandardModelDynamicsMap).

use crate::dynamics_map::StandardModelDynamicsMap;

#[derive(Debug, Clone, PartialEq)]
pub struct EmulatorState {
    pub tick: u64,
    pub gauge: [f64; 3],
    pub matter: [f64; 4],
}

impl Default for EmulatorState {
    fn default() -> Self {
        Self {
            tick: 0,
            gauge: [0.0, 0.0, 0.0],
            matter: [1.0, 0.0, 0.0, 0.0],
        }
    }
}

#[derive(Debug, Clone)]
pub struct EmulatorConfig {
    pub dt: f64,
    pub damping: f64,
}

impl Default for EmulatorConfig {
    fn default() -> Self {
        Self {
            dt: 0.01,
            damping: 0.02,
        }
    }
}

#[derive(Debug, Clone)]
pub struct RuntimeSmEmulator {
    pub map: StandardModelDynamicsMap,
    pub cfg: EmulatorConfig,
}

pub fn encode_emulator_state_payload(s: &EmulatorState) -> Vec<u8> {
    let mut out = Vec::with_capacity(8 + (3 + 4) * 8);
    out.extend_from_slice(&s.tick.to_le_bytes());
    for &g in &s.gauge {
        out.extend_from_slice(&g.to_le_bytes());
    }
    for &m in &s.matter {
        out.extend_from_slice(&m.to_le_bytes());
    }
    out
}

pub fn decode_emulator_state_payload(bytes: &[u8]) -> Option<EmulatorState> {
    const N: usize = 8 + (3 + 4) * 8;
    if bytes.len() != N {
        return None;
    }
    let mut o = 0usize;
    let read_u64 = |b: &[u8], off: &mut usize| {
        let v = u64::from_le_bytes(b[*off..*off + 8].try_into().ok()?);
        *off += 8;
        Some(v)
    };
    let read_f64 = |b: &[u8], off: &mut usize| {
        let v = f64::from_le_bytes(b[*off..*off + 8].try_into().ok()?);
        *off += 8;
        Some(v)
    };
    let tick = read_u64(bytes, &mut o)?;
    let mut gauge = [0.0_f64; 3];
    for g in &mut gauge {
        *g = read_f64(bytes, &mut o)?;
    }
    let mut matter = [0.0_f64; 4];
    for m in &mut matter {
        *m = read_f64(bytes, &mut o)?;
    }
    Some(EmulatorState {
        tick,
        gauge,
        matter,
    })
}

impl RuntimeSmEmulator {
    pub fn new(map: StandardModelDynamicsMap, cfg: EmulatorConfig) -> Self {
        Self { map, cfg }
    }

    pub fn step(&self, s: &mut EmulatorState) {
        let dt = self.cfg.dt.max(0.0);
        let damp = self.cfg.damping.max(0.0);

        // Gauge coupling proxy from theorem-linked electroweak structure.
        let ew_scale = self.map.sin2_theta_w; // 3/13
        let qcd_scale = self.map.beta0 / 20.0; // ~0.97 from 58/3

        let g0 = s.gauge[0];
        let g1 = s.gauge[1];
        let g2 = s.gauge[2];
        let m0 = s.matter[0];
        let m1 = s.matter[1];
        let m2 = s.matter[2];
        let m3 = s.matter[3];

        // Coupled, deterministic, bounded update.
        s.gauge[0] = g0 + dt * (ew_scale * (m0 - m1) - damp * g0);
        s.gauge[1] = g1 + dt * (ew_scale * (m2 - m3) - damp * g1);
        s.gauge[2] = g2 + dt * (qcd_scale * (m0 + m1 + m2 + m3 - 1.0) - damp * g2);

        let drift = self.map.lambda_qg * 0.1; // theorem-linked UV correction scale
        s.matter[0] = (m0 + dt * (-g0 + drift * (m1 - m0))).max(0.0);
        s.matter[1] = (m1 + dt * (g0 + drift * (m0 - m1))).max(0.0);
        s.matter[2] = (m2 + dt * (-g1 + drift * (m3 - m2))).max(0.0);
        s.matter[3] = (m3 + dt * (g1 + drift * (m2 - m3))).max(0.0);

        // Renormalize matter occupancy to a stable simplex.
        let sum = s.matter.iter().sum::<f64>();
        if sum > 1e-12 {
            for v in &mut s.matter {
                *v /= sum;
            }
        } else {
            s.matter = [1.0, 0.0, 0.0, 0.0];
        }

        s.tick = s.tick.saturating_add(1);
    }

    pub fn run(&self, steps: usize, mut state: EmulatorState) -> EmulatorState {
        for _ in 0..steps {
            self.step(&mut state);
        }
        state
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_same_input_same_output() {
        let emu = RuntimeSmEmulator::new(
            StandardModelDynamicsMap::from_clifford_z3(),
            EmulatorConfig::default(),
        );
        let s0 = EmulatorState::default();
        let a = emu.run(250, s0.clone());
        let b = emu.run(250, s0);
        assert_eq!(a, b);
    }

    #[test]
    fn matter_distribution_stays_normalized() {
        let emu = RuntimeSmEmulator::new(
            StandardModelDynamicsMap::from_clifford_z3(),
            EmulatorConfig::default(),
        );
        let s = emu.run(500, EmulatorState::default());
        let sum = s.matter.iter().sum::<f64>();
        assert!((sum - 1.0).abs() < 1e-9, "matter sum drifted: {sum}");
        assert!(s.matter.iter().all(|v| *v >= 0.0));
    }

    #[test]
    fn state_payload_roundtrip() {
        let s = EmulatorState {
            tick: 123,
            gauge: [0.1, -0.2, 0.3],
            matter: [0.2, 0.3, 0.1, 0.4],
        };
        let p = encode_emulator_state_payload(&s);
        let s2 = decode_emulator_state_payload(&p).expect("decode payload");
        assert_eq!(s, s2);
    }
}
