# Finding 008: Seeded Universe Tick Scheduler + Deterministic Replay

Date: 2026-02-24  
Related issues: GRAND-192, GRAND-196, GRAND-197, GRAND-198

## What was implemented

- Added deterministic runtime tick driver:
  - `crates/gutoe-gpu/src/bin/sm_tick.rs`
- Uses theorem-linked runtime map + emulator:
  - `gutoe-physics::StandardModelDynamicsMap`
  - `gutoe-physics::RuntimeSmEmulator`
- Emits versioned snapshots every `SM_SNAPSHOT_STRIDE` ticks via:
  - `gutoe-gpu::snapshot::{write_snapshot_file, read_snapshot_file}`
- Replays snapshots with interpolation via:
  - `gutoe-gpu::playback::SnapshotTrack`

## Run and artifact

- Command:
  - `cargo run -p gutoe-gpu --bin sm_tick`
- Default output:
  - `/tmp/bh_renders/sm_snapshots`
- Default run produced:
  - `steps=1000`, `stride=50`, `snapshots=21`
  - Replay probe at tick `500` with interpolated `matter[0]=0.053233037`

## Why this matters

- Confirms seeded, deterministic tick progression is now executable.
- Confirms deterministic checkpoint serialization/deserialization path.
- Confirms fast seek/interpolation over checkpoint timelines.
- Provides direct scaffolding for renderer-time streaming (`GRAND-191`) and full runtime bridge completion.

