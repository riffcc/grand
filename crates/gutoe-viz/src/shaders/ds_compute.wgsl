// GUTOE Double-Slit Quantum Compute Shader
// Each invocation simulates one particle.
// Dead particles (state==2) are respawned at the source.
// Particles pass through slits with quantum-mechanical diffraction kicks.
// Interference pattern accumulates in `pattern` via atomicAdd.

struct Particle {
    pos: vec2<f32>,   // pixel coords
    vel: vec2<f32>,   // pixels/tick
    state: u32,       // 0=pre-barrier, 1=post-slit, 2=dead
    seed: u32,        // per-particle rng seed
}

struct Uniforms {
    width:          f32,
    height:         f32,
    barrier_x:      f32,
    slit1_y:        f32,
    slit2_y:        f32,
    slit_half_w:    f32,
    screen_x:       f32,
    source_x:       f32,
    source_y:       f32,
    particle_speed: f32,
    frame:          u32,
    running:        u32,
    n_particles:    u32,
    n_buckets:      u32,
    _pad0:          u32,
    _pad1:          u32,
}

@group(0) @binding(0) var<storage, read_write> particles: array<Particle>;
@group(0) @binding(1) var<storage, read_write> pattern:   array<atomic<u32>>;
@group(0) @binding(2) var<uniform>             u:         Uniforms;

// --- Pseudorandom (Wang hash) ---
fn wang_hash(v: u32) -> u32 {
    var x = v;
    x  = (x ^ 61u) ^ (x >> 16u);
    x *= 9u;
    x ^= x >> 4u;
    x *= 0x27d4eb2du;
    x ^= x >> 15u;
    return x;
}

// Returns uniform float in [0, 1)
fn frand(seed: u32) -> f32 {
    return f32(wang_hash(seed) & 0x00FFFFFFu) / f32(0x01000000u);
}

// Returns uniform float in [lo, hi)
fn frange(seed: u32, lo: f32, hi: f32) -> f32 {
    return lo + (hi - lo) * frand(seed);
}

// Spawn a fresh particle at the source
fn spawn_particle(idx: u32, frame_seed: u32) -> Particle {
    let s0 = wang_hash(frame_seed ^ (idx * 2654435761u));
    let s1 = wang_hash(s0 + 1u);
    let s2 = wang_hash(s1 + 2u);
    let jitter_y = frange(s2, -4.0, 4.0);
    let vy       = frange(wang_hash(s2 + 3u), -0.4, 0.4);
    return Particle(
        vec2<f32>(u.source_x, u.source_y + jitter_y),
        vec2<f32>(u.particle_speed, vy),
        0u,
        wang_hash(s0 ^ s2),
    );
}

// Sample deflection vy from the double-slit interference distribution
// P(y_screen) ∝ cos²(π·d·y / λ·L)
// Uses rejection sampling with up to `max_tries` attempts.
fn sample_diffraction_vy(seed: u32) -> f32 {
    let slit_sep      = abs(u.slit2_y - u.slit1_y);
    let L             = u.screen_x - u.barrier_x;
    let lambda        = 22.0;  // effective de Broglie λ in screen pixels — tunable
    let k             = 3.14159265358979 * slit_sep / (lambda * L);
    let max_offset    = u.height * 0.44;

    var s = seed;
    for (var i = 0u; i < 40u; i++) {
        s = wang_hash(s + i * 1013u);
        let y_off = frange(s, -max_offset, max_offset);
        s = wang_hash(s + 7u);
        let cos_val  = cos(k * y_off);
        let prob     = cos_val * cos_val;         // cos² ∈ [0, 1]
        if frand(s) < prob {
            // Convert y offset at screen to vy needed to reach it
            return (y_off / L) * u.particle_speed;
        }
    }
    // Rejection ran out — return small random deflection
    return frange(wang_hash(seed + 9999u), -1.0, 1.0);
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let idx = gid.x;
    if idx >= u.n_particles { return; }

    var p = particles[idx];

    // ── Respawn dead particles ──────────────────────────────────────────────
    if p.state == 2u {
        particles[idx] = spawn_particle(idx, u.frame * 6271u);
        return;
    }

    // ── Simulation paused ───────────────────────────────────────────────────
    if u.running == 0u {
        particles[idx] = p;
        return;
    }

    // ── Advance position ────────────────────────────────────────────────────
    p.pos += p.vel;
    p.seed = wang_hash(p.seed ^ (u.frame * 31337u));

    // ── Barrier collision (pre-slit only) ───────────────────────────────────
    if p.state == 0u {
        let dx = p.pos.x - u.barrier_x;
        if dx >= -2.0 && dx <= 4.0 {
            let in_s1 = abs(p.pos.y - u.slit1_y) <= u.slit_half_w;
            let in_s2 = abs(p.pos.y - u.slit2_y) <= u.slit_half_w;

            if in_s1 || in_s2 {
                p.state    = 1u;
                p.vel.y    = sample_diffraction_vy(p.seed);
                p.pos.x    = u.barrier_x + 5.0;  // snap past barrier face
            } else {
                p.state = 2u;  // blocked — respawn next tick
            }
        }
    }

    // ── Detection screen ────────────────────────────────────────────────────
    if p.state == 1u && p.pos.x >= u.screen_x {
        let y_frac = p.pos.y / u.height;
        let bidx   = i32(y_frac * f32(u.n_buckets));
        if bidx >= 0 && bidx < i32(u.n_buckets) {
            atomicAdd(&pattern[u32(bidx)], 1u);
        }
        p.state = 2u;
    }

    // ── Out of bounds ───────────────────────────────────────────────────────
    if p.pos.y < -10.0 || p.pos.y > u.height + 10.0 { p.state = 2u; }
    if p.pos.x < -20.0 || p.pos.x > u.width  + 20.0 { p.state = 2u; }

    particles[idx] = p;
}
