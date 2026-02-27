use std::fs;
use std::path::PathBuf;

use gutoe_gpu::playback::SnapshotTrack;
use gutoe_gpu::snapshot::{write_snapshot_file, UniverseSnapshot};
use gutoe_physics::{
    decode_emulator_state_payload, encode_emulator_state_payload, EmulatorConfig, EmulatorState,
    RuntimeSmEmulator, StandardModelDynamicsMap,
};

fn parse_env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn parse_env_u64(name: &str, default: u64) -> u64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(default)
}

fn parse_env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let out_dir = std::env::var("SM_SNAPSHOT_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/tmp/bh_renders/sm_snapshots"));
    let total_steps = parse_env_usize("SM_TICK_STEPS", 1000);
    let stride = parse_env_usize("SM_SNAPSHOT_STRIDE", 50).max(1);
    let seed = parse_env_u64("SM_TICK_SEED", 1337);
    let dt = parse_env_f64("SM_DT", 0.01);
    let damping = parse_env_f64("SM_DAMPING", 0.02);

    fs::create_dir_all(&out_dir).expect("create snapshot dir");

    let map = StandardModelDynamicsMap::from_clifford_z3();
    let emu = RuntimeSmEmulator::new(map, EmulatorConfig { dt, damping });
    let mut state = EmulatorState::default();

    let mut written = 0usize;
    for step in 0..=total_steps {
        if step % stride == 0 {
            let payload = encode_emulator_state_payload(&state);
            let snap = UniverseSnapshot {
                tick: state.tick,
                seed,
                sim_time: state.tick as f64 * dt,
                payload,
            };
            let p = out_dir.join(format!("tick_{:08}.gts", state.tick));
            write_snapshot_file(&p, &snap).expect("write snapshot");
            written += 1;
        }
        if step < total_steps {
            emu.step(&mut state);
        }
    }

    let track = SnapshotTrack::load_from_dir(&out_dir).expect("load snapshot track");
    let probe_tick = (total_steps as u64 / 2).max(1);
    let interpolated_q0 = track.sample_tick_with(probe_tick, |a, b, alpha| {
        let sa = decode_emulator_state_payload(&a.payload).expect("decode a");
        let sb = decode_emulator_state_payload(&b.payload).expect("decode b");
        sa.matter[0] * (1.0 - alpha) + sb.matter[0] * alpha
    });

    println!(
        "sm_tick complete: steps={} stride={} snapshots={} dir={}",
        total_steps,
        stride,
        written,
        out_dir.display()
    );
    if let Some(q0) = interpolated_q0 {
        println!(
            "playback probe: tick={} interpolated matter[0]={:.9}",
            probe_tick, q0
        );
    }
}
