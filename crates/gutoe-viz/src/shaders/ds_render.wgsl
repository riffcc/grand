// GUTOE Double-Slit Render Shader
// Three draw calls share this file via entry-point selection:
//   vs_barrier / fs_barrier   — opaque barrier geometry
//   vs_particle / fs_particle — additive glow particles
//   vs_pattern  / fs_pattern  — additive interference pattern bar

struct Particle {
    pos:   vec2<f32>,
    vel:   vec2<f32>,
    state: u32,
    seed:  u32,
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
    // pattern normalization: expected peak bucket count this frame
    expected_peak:  f32,
    _pad1:          u32,
}

@group(0) @binding(0) var<storage, read> particles: array<Particle>;
@group(0) @binding(1) var<storage, read> pattern:   array<u32>;
@group(0) @binding(2) var<uniform>       u:         Uniforms;

// ── Helpers ──────────────────────────────────────────────────────────────────

// Pixel → NDC (wgpu: y=0 is top in pixel space, y=+1 is top in NDC)
fn to_ndc(px: vec2<f32>) -> vec2<f32> {
    return vec2<f32>(
         px.x / u.width  * 2.0 - 1.0,
        -(px.y / u.height * 2.0 - 1.0),
    );
}

// ── Barrier ──────────────────────────────────────────────────────────────────
// Fullscreen quad; fragment discards everything except the barrier strip and
// the two slit openings.

struct BarrierVOut {
    @builtin(position) pos: vec4<f32>,
}

const QUAD: array<vec2<f32>, 6> = array<vec2<f32>, 6>(
    vec2<f32>(-1.0, -1.0),
    vec2<f32>( 1.0, -1.0),
    vec2<f32>(-1.0,  1.0),
    vec2<f32>(-1.0,  1.0),
    vec2<f32>( 1.0, -1.0),
    vec2<f32>( 1.0,  1.0),
);

@vertex
fn vs_barrier(@builtin(vertex_index) vi: u32) -> BarrierVOut {
    return BarrierVOut(vec4<f32>(QUAD[vi], 0.0, 1.0));
}

@fragment
fn fs_barrier(@builtin(position) frag: vec4<f32>) -> @location(0) vec4<f32> {
    let px = frag.x;
    let py = frag.y;

    // Only render the barrier column
    if abs(px - u.barrier_x) > 4.0 { discard; }

    // Slit openings: let light through
    if abs(py - u.slit1_y) <= u.slit_half_w { discard; }
    if abs(py - u.slit2_y) <= u.slit_half_w { discard; }

    // Barrier body: dark steel with a slight blue tinge
    return vec4<f32>(0.28, 0.30, 0.38, 1.0);
}

// ── Slit glow overlay (separate from barrier, draws over it) ─────────────────
// Drawn as fullscreen quad after barrier, shows faint cyan glow at slit edges.
@vertex
fn vs_slit_glow(@builtin(vertex_index) vi: u32) -> BarrierVOut {
    return BarrierVOut(vec4<f32>(QUAD[vi], 0.0, 1.0));
}

@fragment
fn fs_slit_glow(@builtin(position) frag: vec4<f32>) -> @location(0) vec4<f32> {
    let px = frag.x;
    let py = frag.y;

    // Only near barrier
    let dist_to_barrier = abs(px - u.barrier_x);
    if dist_to_barrier > 12.0 { discard; }

    let d1 = abs(py - u.slit1_y) - u.slit_half_w;  // >0 outside slit
    let d2 = abs(py - u.slit2_y) - u.slit_half_w;
    let edge_dist = min(max(d1, 0.0), max(d2, 0.0));

    // Glow on slit edges
    let glow = exp(-edge_dist * 0.15) * exp(-dist_to_barrier * 0.25) * 0.6;
    if glow < 0.01 { discard; }

    return vec4<f32>(0.0, 0.85, 1.0, glow);
}

// ── Particles ─────────────────────────────────────────────────────────────────
// Instanced quads: instance_index → particle index.

struct ParticleVOut {
    @builtin(position) pos:   vec4<f32>,
    @location(0)       uv:    vec2<f32>,   // [-1,1]² within the quad
    @location(1)       color: vec4<f32>,
    @location(2)       alive: f32,
}

// Unit quad offsets for 2-triangle strip (6 verts)
const PQUAD: array<vec2<f32>, 6> = array<vec2<f32>, 6>(
    vec2<f32>(-1.0, -1.0),
    vec2<f32>( 1.0, -1.0),
    vec2<f32>(-1.0,  1.0),
    vec2<f32>(-1.0,  1.0),
    vec2<f32>( 1.0, -1.0),
    vec2<f32>( 1.0,  1.0),
);

@vertex
fn vs_particle(
    @builtin(vertex_index)   vi: u32,
    @builtin(instance_index) ii: u32,
) -> ParticleVOut {
    let p      = particles[ii];
    let offset = PQUAD[vi];

    let alive = f32(p.state != 2u);

    // Particle radius in pixels
    let r = 5.0;
    let ndc_pos    = to_ndc(p.pos);
    let ndc_offset = vec2<f32>(offset.x * r / u.width * 2.0,
                               offset.y * r / u.height * 2.0);

    // Color by phase:
    //   pre-barrier  → warm amber
    //   post-slit    → cool cyan
    var col = vec4<f32>(0.0);
    if p.state == 0u {
        col = vec4<f32>(1.0, 0.55, 0.05, 1.0);   // amber
    } else if p.state == 1u {
        // Slight speed-based tint
        let speed_t = clamp(abs(p.vel.y) / 3.0, 0.0, 1.0);
        col = mix(vec4<f32>(0.0, 1.0, 0.92, 1.0),   // cyan
                  vec4<f32>(0.4, 0.0, 1.0, 1.0),     // violet
                  speed_t);
    }

    var out: ParticleVOut;
    out.pos   = vec4<f32>(ndc_pos + ndc_offset, 0.0, 1.0);
    out.uv    = offset;
    out.color = col;
    out.alive = alive;
    return out;
}

@fragment
fn fs_particle(in: ParticleVOut) -> @location(0) vec4<f32> {
    if in.alive < 0.5 { discard; }

    let dist = length(in.uv);
    if dist > 1.0 { discard; }

    // Soft circular glow: bright centre, fades at edge
    let glow = (1.0 - dist * dist) * (1.0 - dist * 0.3);
    return in.color * glow;
}

// ── Interference pattern bar ──────────────────────────────────────────────────
// Drawn as N_BUCKETS instanced horizontal bars on the right side of the screen.
// instance_index → bucket index.

struct PatternVOut {
    @builtin(position) pos:       vec4<f32>,
    @location(0)       intensity: f32,
}

@vertex
fn vs_pattern(
    @builtin(vertex_index)   vi: u32,
    @builtin(instance_index) ii: u32,
) -> PatternVOut {
    let count     = f32(pattern[ii]);
    let peak      = max(u.expected_peak, 1.0);
    let intensity = clamp(count / peak, 0.0, 1.0);

    // Each bucket is a horizontal bar
    let n         = f32(u.n_buckets);
    let bucket_h  = u.height / n;
    let bar_max_w = u.width * 0.18;  // max bar width = 18% of screen

    let y_top = f32(ii)       * bucket_h;
    let y_bot = f32(ii + 1u)  * bucket_h;
    let x_right = u.screen_x + 6.0;          // start right of screen line
    let x_left  = x_right + intensity * bar_max_w;

    // Two triangles forming the bar
    let corners = array<vec2<f32>, 6>(
        vec2<f32>(x_right, y_top),
        vec2<f32>(x_left,  y_top),
        vec2<f32>(x_right, y_bot),
        vec2<f32>(x_right, y_bot),
        vec2<f32>(x_left,  y_top),
        vec2<f32>(x_left,  y_bot),
    );

    var out: PatternVOut;
    out.pos       = vec4<f32>(to_ndc(corners[vi]), 0.0, 1.0);
    out.intensity = intensity;
    return out;
}

@fragment
fn fs_pattern(in: PatternVOut) -> @location(0) vec4<f32> {
    if in.intensity < 0.001 { discard; }

    // Colour: interference peaks are bright cyan-white, troughs fade out
    let t   = in.intensity;
    let col = mix(vec4<f32>(0.0, 0.5, 0.8, 1.0),    // dark blue trough
                  vec4<f32>(0.1, 1.0, 0.95, 1.0),   // bright cyan peak
                  t);
    return col * (t * 0.9 + 0.1);
}

// ── Screen line ───────────────────────────────────────────────────────────────
// A thin vertical line at screen_x showing the detection screen.

@vertex
fn vs_screen_line(@builtin(vertex_index) vi: u32) -> BarrierVOut {
    return BarrierVOut(vec4<f32>(QUAD[vi], 0.0, 1.0));
}

@fragment
fn fs_screen_line(@builtin(position) frag: vec4<f32>) -> @location(0) vec4<f32> {
    if abs(frag.x - u.screen_x) > 1.0 { discard; }
    return vec4<f32>(0.15, 0.45, 0.6, 0.8);
}
