use crate::snapshot::read_snapshot_file;
use crate::snapshot::UniverseSnapshot;
use std::fs;
use std::io;
use std::path::Path;

#[derive(Debug, Clone, PartialEq)]
pub struct SnapshotTrack {
    frames: Vec<UniverseSnapshot>,
}

impl SnapshotTrack {
    pub fn load_from_dir(dir: &Path) -> io::Result<Self> {
        let mut frames = Vec::new();
        for ent in fs::read_dir(dir)? {
            let p = ent?.path();
            if p.extension().and_then(|s| s.to_str()) != Some("gts") {
                continue;
            }
            if let Ok(snap) = read_snapshot_file(&p) {
                frames.push(snap);
            }
        }
        Ok(Self::new(frames))
    }

    pub fn new(mut frames: Vec<UniverseSnapshot>) -> Self {
        frames.sort_by_key(|f| f.tick);
        Self { frames }
    }

    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    pub fn len(&self) -> usize {
        self.frames.len()
    }

    pub fn frames(&self) -> &[UniverseSnapshot] {
        &self.frames
    }

    pub fn seek_tick(&self, tick: u64) -> Option<(&UniverseSnapshot, &UniverseSnapshot, f64)> {
        if self.frames.is_empty() {
            return None;
        }
        if self.frames.len() == 1 {
            let f = &self.frames[0];
            return Some((f, f, 0.0));
        }
        if tick <= self.frames[0].tick {
            let f = &self.frames[0];
            return Some((f, f, 0.0));
        }
        let last = self.frames.len() - 1;
        if tick >= self.frames[last].tick {
            let f = &self.frames[last];
            return Some((f, f, 0.0));
        }
        let hi = self.frames.partition_point(|f| f.tick < tick);
        let lo = hi - 1;
        let a = &self.frames[lo];
        let b = &self.frames[hi];
        let dt = (b.tick - a.tick) as f64;
        let alpha = if dt <= 0.0 {
            0.0
        } else {
            ((tick - a.tick) as f64 / dt).clamp(0.0, 1.0)
        };
        Some((a, b, alpha))
    }

    pub fn seek_time(&self, sim_time: f64) -> Option<(&UniverseSnapshot, &UniverseSnapshot, f64)> {
        if self.frames.is_empty() {
            return None;
        }
        if self.frames.len() == 1 {
            let f = &self.frames[0];
            return Some((f, f, 0.0));
        }
        if sim_time <= self.frames[0].sim_time {
            let f = &self.frames[0];
            return Some((f, f, 0.0));
        }
        let last = self.frames.len() - 1;
        if sim_time >= self.frames[last].sim_time {
            let f = &self.frames[last];
            return Some((f, f, 0.0));
        }
        let hi = self.frames.partition_point(|f| f.sim_time < sim_time);
        let lo = hi - 1;
        let a = &self.frames[lo];
        let b = &self.frames[hi];
        let dt = b.sim_time - a.sim_time;
        let alpha = if dt <= 1e-12 {
            0.0
        } else {
            ((sim_time - a.sim_time) / dt).clamp(0.0, 1.0)
        };
        Some((a, b, alpha))
    }

    pub fn sample_tick_with<T, F>(&self, tick: u64, interp: F) -> Option<T>
    where
        F: Fn(&UniverseSnapshot, &UniverseSnapshot, f64) -> T,
    {
        let (a, b, alpha) = self.seek_tick(tick)?;
        Some(interp(a, b, alpha))
    }

    pub fn sample_time_with<T, F>(&self, sim_time: f64, interp: F) -> Option<T>
    where
        F: Fn(&UniverseSnapshot, &UniverseSnapshot, f64) -> T,
    {
        let (a, b, alpha) = self.seek_time(sim_time)?;
        Some(interp(a, b, alpha))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame(tick: u64, sim_time: f64, byte: u8) -> UniverseSnapshot {
        UniverseSnapshot {
            tick,
            seed: 7,
            sim_time,
            payload: vec![byte],
        }
    }

    #[test]
    fn seek_tick_interpolates_between_neighbors() {
        let t = SnapshotTrack::new(vec![
            frame(10, 1.0, 1),
            frame(20, 2.0, 2),
            frame(30, 3.0, 3),
        ]);
        let (a, b, alpha) = t.seek_tick(15).expect("seek");
        assert_eq!(a.tick, 10);
        assert_eq!(b.tick, 20);
        assert!((alpha - 0.5).abs() < 1e-9);
    }

    #[test]
    fn seek_time_interpolates_between_neighbors() {
        let t = SnapshotTrack::new(vec![
            frame(10, 1.0, 1),
            frame(20, 2.0, 2),
            frame(30, 3.0, 3),
        ]);
        let (a, b, alpha) = t.seek_time(2.5).expect("seek");
        assert_eq!(a.tick, 20);
        assert_eq!(b.tick, 30);
        assert!((alpha - 0.5).abs() < 1e-9);
    }

    #[test]
    fn sample_with_closure_works() {
        let t = SnapshotTrack::new(vec![frame(0, 0.0, 10), frame(10, 1.0, 20)]);
        let v = t
            .sample_tick_with(5, |a, b, alpha| {
                (a.payload[0] as f64 * (1.0 - alpha) + b.payload[0] as f64 * alpha).round() as u8
            })
            .expect("sample");
        assert_eq!(v, 15);
    }
}
