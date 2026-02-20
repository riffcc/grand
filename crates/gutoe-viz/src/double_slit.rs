/*!
 * GUTOE Viz - Double-Slit Experiment
 * Copyright (C) 2026  Wings
 *
 * AGPL-3.0-or-later
 */

//! GPU-accelerated double-slit quantum simulation.
//!
//! Architecture:
//!   1. Compute pass: 65 536 particles simulated in parallel on the GPU.
//!      Each particle either travels free-flight, diffracts through a slit,
//!      or dies and respawns. The interference pattern accumulates in an
//!      atomic storage buffer.
//!   2. Render pass (4 sub-draws):
//!      a. Barrier (opaque) — wall with two slit openings
//!      b. Slit glow (additive) — subtle edge illumination
//!      c. Particles (additive) — glowing dots coloured by phase
//!      d. Pattern bar (additive) — right-side histogram

use bytemuck::{Pod, Zeroable};
use wgpu::{util::DeviceExt, *};

/// Number of simulated particles (must be a multiple of workgroup_size=64).
const N_PARTICLES: u32 = 65_536;

/// Number of histogram buckets for the interference pattern.
const N_BUCKETS: u32 = 512;

// ─── GPU-side structs (must match WGSL layouts exactly) ─────────────────────

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct Particle {
    pos:   [f32; 2],
    vel:   [f32; 2],
    state: u32,   // 0=pre-barrier, 1=post-slit, 2=dead
    seed:  u32,
}

/// Uniform block shared between compute and render shaders.
/// **Must** be 16-byte aligned (pad to 64 bytes total = 16 u32s).
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct DsUniforms {
    width:          f32,
    height:         f32,
    barrier_x:      f32,
    slit1_y:        f32,   // pixel y of centre of slit 1
    slit2_y:        f32,   // pixel y of centre of slit 2
    slit_half_w:    f32,   // half-width of each slit in pixels
    screen_x:       f32,   // x coordinate of detection screen
    source_x:       f32,   // x coordinate of particle source
    source_y:       f32,   // y coordinate of particle source (= height/2)
    particle_speed: f32,   // horizontal velocity in pixels/tick
    frame:          u32,
    running:        u32,   // 1 = simulating, 0 = paused
    n_particles:    u32,
    n_buckets:      u32,
    expected_peak:  f32,   // used by pattern shader to normalise bar widths
    _pad1:          u32,
}

// ─── Blend state helpers ─────────────────────────────────────────────────────

fn additive_blend() -> BlendState {
    BlendState {
        color: BlendComponent {
            src_factor: BlendFactor::One,
            dst_factor: BlendFactor::One,
            operation:  BlendOperation::Add,
        },
        alpha: BlendComponent {
            src_factor: BlendFactor::One,
            dst_factor: BlendFactor::One,
            operation:  BlendOperation::Add,
        },
    }
}

fn alpha_blend() -> BlendState {
    BlendState {
        color: BlendComponent {
            src_factor: BlendFactor::SrcAlpha,
            dst_factor: BlendFactor::OneMinusSrcAlpha,
            operation:  BlendOperation::Add,
        },
        alpha: BlendComponent {
            src_factor: BlendFactor::One,
            dst_factor: BlendFactor::Zero,
            operation:  BlendOperation::Add,
        },
    }
}

// ─── Main state object ───────────────────────────────────────────────────────

pub struct DoubleSlitState {
    // Simulation parameters (CPU side)
    width:          u32,
    height:         u32,
    running:        bool,
    frame:          u32,
    slit_gap:       f32,  // distance between slit centres in pixels

    // GPU buffers
    particle_buf:   Buffer,   // Particle array (compute rw / render ro)
    pattern_buf:    Buffer,   // u32 histogram (compute atomic rw / render ro)
    uniform_buf:    Buffer,   // DsUniforms (UNIFORM | COPY_DST)
    _clear_buf:     Buffer,   // zero-filled source for resetting pattern (reserved)

    // Bind groups
    compute_bg:     BindGroup,
    render_bg:      BindGroup,

    // Pipelines
    compute_pipeline:      ComputePipeline,
    barrier_pipeline:      RenderPipeline,
    slit_glow_pipeline:    RenderPipeline,
    particle_pipeline:     RenderPipeline,
    pattern_pipeline:      RenderPipeline,
    screen_line_pipeline:  RenderPipeline,
}

impl DoubleSlitState {
    pub fn new(device: &Device, _queue: &Queue, format: TextureFormat, w: u32, h: u32) -> Self {
        // ── Geometry ─────────────────────────────────────────────────────────
        let (barrier_x, screen_x, slit_gap) = geometry(w, h);
        let slit1_y = h as f32 / 2.0 - slit_gap / 2.0;
        let slit2_y = h as f32 / 2.0 + slit_gap / 2.0;

        // ── Uniform buffer ───────────────────────────────────────────────────
        let uniforms = build_uniforms(w, h, barrier_x, slit1_y, slit2_y, slit_gap, screen_x, true, 0, 1.0);
        let uniform_buf = device.create_buffer_init(&util::BufferInitDescriptor {
            label:    Some("ds-uniforms"),
            contents: bytemuck::bytes_of(&uniforms),
            usage:    BufferUsages::UNIFORM | BufferUsages::COPY_DST,
        });

        // ── Particle buffer ──────────────────────────────────────────────────
        let particles = init_particles(w, h, N_PARTICLES);
        let particle_buf = device.create_buffer_init(&util::BufferInitDescriptor {
            label:    Some("ds-particles"),
            contents: bytemuck::cast_slice(&particles),
            usage:    BufferUsages::STORAGE | BufferUsages::COPY_DST,
        });

        // ── Pattern buffer ───────────────────────────────────────────────────
        let pattern_bytes = vec![0u8; (N_BUCKETS * 4) as usize];
        let pattern_buf = device.create_buffer_init(&util::BufferInitDescriptor {
            label:    Some("ds-pattern"),
            contents: &pattern_bytes,
            usage:    BufferUsages::STORAGE | BufferUsages::COPY_DST | BufferUsages::COPY_SRC,
        });

        // ── Clear buffer (zeroes, used to reset pattern each experiment) ─────
        let clear_buf = device.create_buffer_init(&util::BufferInitDescriptor {
            label:    Some("ds-pattern-clear"),
            contents: &pattern_bytes,
            usage:    BufferUsages::COPY_SRC,
        });

        // ── Shader modules ───────────────────────────────────────────────────
        let compute_shader = device.create_shader_module(ShaderModuleDescriptor {
            label:  Some("ds-compute"),
            source: ShaderSource::Wgsl(
                include_str!("shaders/ds_compute.wgsl").into(),
            ),
        });
        let render_shader = device.create_shader_module(ShaderModuleDescriptor {
            label:  Some("ds-render"),
            source: ShaderSource::Wgsl(
                include_str!("shaders/ds_render.wgsl").into(),
            ),
        });

        // ── Compute bind group layout ────────────────────────────────────────
        let compute_bgl = device.create_bind_group_layout(&BindGroupLayoutDescriptor {
            label:   Some("ds-compute-bgl"),
            entries: &[
                BindGroupLayoutEntry {
                    binding:    0,
                    visibility: ShaderStages::COMPUTE,
                    ty: BindingType::Buffer {
                        ty:                 BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size:   None,
                    },
                    count: None,
                },
                BindGroupLayoutEntry {
                    binding:    1,
                    visibility: ShaderStages::COMPUTE,
                    ty: BindingType::Buffer {
                        ty:                 BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size:   None,
                    },
                    count: None,
                },
                BindGroupLayoutEntry {
                    binding:    2,
                    visibility: ShaderStages::COMPUTE,
                    ty: BindingType::Buffer {
                        ty:                 BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size:   None,
                    },
                    count: None,
                },
            ],
        });

        let compute_bg = device.create_bind_group(&BindGroupDescriptor {
            label:   Some("ds-compute-bg"),
            layout:  &compute_bgl,
            entries: &[
                BindGroupEntry { binding: 0, resource: particle_buf.as_entire_binding() },
                BindGroupEntry { binding: 1, resource: pattern_buf.as_entire_binding() },
                BindGroupEntry { binding: 2, resource: uniform_buf.as_entire_binding() },
            ],
        });

        // ── Render bind group layout ─────────────────────────────────────────
        let render_bgl = device.create_bind_group_layout(&BindGroupLayoutDescriptor {
            label:   Some("ds-render-bgl"),
            entries: &[
                BindGroupLayoutEntry {
                    binding:    0,
                    visibility: ShaderStages::VERTEX | ShaderStages::FRAGMENT,
                    ty: BindingType::Buffer {
                        ty:                 BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size:   None,
                    },
                    count: None,
                },
                BindGroupLayoutEntry {
                    binding:    1,
                    visibility: ShaderStages::VERTEX | ShaderStages::FRAGMENT,
                    ty: BindingType::Buffer {
                        ty:                 BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size:   None,
                    },
                    count: None,
                },
                BindGroupLayoutEntry {
                    binding:    2,
                    visibility: ShaderStages::VERTEX | ShaderStages::FRAGMENT,
                    ty: BindingType::Buffer {
                        ty:                 BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size:   None,
                    },
                    count: None,
                },
            ],
        });

        let render_bg = device.create_bind_group(&BindGroupDescriptor {
            label:   Some("ds-render-bg"),
            layout:  &render_bgl,
            entries: &[
                BindGroupEntry { binding: 0, resource: particle_buf.as_entire_binding() },
                BindGroupEntry { binding: 1, resource: pattern_buf.as_entire_binding() },
                BindGroupEntry { binding: 2, resource: uniform_buf.as_entire_binding() },
            ],
        });

        // ── Compute pipeline ─────────────────────────────────────────────────
        let compute_layout = device.create_pipeline_layout(&PipelineLayoutDescriptor {
            label:                Some("ds-compute-layout"),
            bind_group_layouts:   &[&compute_bgl],
            push_constant_ranges: &[],
        });

        let compute_pipeline = device.create_compute_pipeline(&ComputePipelineDescriptor {
            label:       Some("ds-compute"),
            layout:      Some(&compute_layout),
            module:      &compute_shader,
            entry_point: "main",
            compilation_options: PipelineCompilationOptions::default(),
            cache:       None,
        });

        // ── Render pipeline layout (shared) ──────────────────────────────────
        let render_layout = device.create_pipeline_layout(&PipelineLayoutDescriptor {
            label:                Some("ds-render-layout"),
            bind_group_layouts:   &[&render_bgl],
            push_constant_ranges: &[],
        });

        // Helper: build a render pipeline with a given blend state and entry points.
        let make_pipeline = |label: &str,
                             vs_entry: &str,
                             fs_entry: &str,
                             blend: Option<BlendState>| {
            device.create_render_pipeline(&RenderPipelineDescriptor {
                label: Some(label),
                layout: Some(&render_layout),
                vertex: VertexState {
                    module:             &render_shader,
                    entry_point:        vs_entry,
                    buffers:            &[],
                    compilation_options: PipelineCompilationOptions::default(),
                },
                fragment: Some(FragmentState {
                    module:             &render_shader,
                    entry_point:        fs_entry,
                    compilation_options: PipelineCompilationOptions::default(),
                    targets: &[Some(ColorTargetState {
                        format,
                        blend:      Some(blend.unwrap_or(BlendState::REPLACE)),
                        write_mask: ColorWrites::ALL,
                    })],
                }),
                primitive: PrimitiveState {
                    topology:           PrimitiveTopology::TriangleList,
                    ..Default::default()
                },
                depth_stencil: None,
                multisample:   MultisampleState::default(),
                multiview:     None,
                cache:         None,
            })
        };

        let barrier_pipeline     = make_pipeline("barrier",     "vs_barrier",     "fs_barrier",     None);
        let slit_glow_pipeline   = make_pipeline("slit-glow",   "vs_slit_glow",   "fs_slit_glow",   Some(additive_blend()));
        let particle_pipeline    = make_pipeline("particles",   "vs_particle",    "fs_particle",    Some(additive_blend()));
        let pattern_pipeline     = make_pipeline("pattern",     "vs_pattern",     "fs_pattern",     Some(additive_blend()));
        let screen_line_pipeline = make_pipeline("screen-line", "vs_screen_line", "fs_screen_line", Some(alpha_blend()));

        Self {
            width:  w,
            height: h,
            running: true,
            frame: 0,
            slit_gap,

            particle_buf,
            pattern_buf,
            uniform_buf,
            _clear_buf: clear_buf,

            compute_bg,
            render_bg,

            compute_pipeline,
            barrier_pipeline,
            slit_glow_pipeline,
            particle_pipeline,
            pattern_pipeline,
            screen_line_pipeline,
        }
    }

    // ── Public controls ──────────────────────────────────────────────────────

    pub fn toggle_running(&mut self) {
        self.running = !self.running;
    }

    pub fn is_running(&self) -> bool { self.running }

    pub fn reset(&mut self, _device: &Device, queue: &Queue) {
        let (barrier_x, screen_x, slit_gap) = geometry(self.width, self.height);
        let slit1_y = self.height as f32 / 2.0 - slit_gap / 2.0;
        let slit2_y = self.height as f32 / 2.0 + slit_gap / 2.0;

        // Re-init particles on CPU, upload
        let particles = init_particles(self.width, self.height, N_PARTICLES);
        queue.write_buffer(&self.particle_buf, 0, bytemuck::cast_slice(&particles));

        // Zero the pattern buffer
        let zeros = vec![0u8; (N_BUCKETS * 4) as usize];
        queue.write_buffer(&self.pattern_buf, 0, &zeros);

        self.frame   = 0;
        self.running = true;
        self.slit_gap = slit_gap;
        self.update_uniforms(queue, barrier_x, slit1_y, slit2_y, screen_x, 1.0);
    }

    pub fn resize(&mut self, w: u32, h: u32, device: &Device, queue: &Queue) {
        self.width  = w;
        self.height = h;
        self.reset(device, queue);
    }

    // ── Per-frame update ─────────────────────────────────────────────────────

    pub fn update(&mut self, queue: &Queue) {
        if !self.running { return; }
        self.frame = self.frame.wrapping_add(1);

        let (barrier_x, screen_x, _) = geometry(self.width, self.height);
        let slit1_y = self.height as f32 / 2.0 - self.slit_gap / 2.0;
        let slit2_y = self.height as f32 / 2.0 + self.slit_gap / 2.0;

        // Rough expected peak for normalization:
        // particles/frame reaching screen ≈ N_PARTICLES * slit_coverage
        // slit_coverage ≈ 2*slit_half_w / height
        let slit_half_w  = self.slit_half_w();
        let slit_coverage = (2.0 * slit_half_w) / self.height as f32;
        let hits_per_frame = N_PARTICLES as f32 * slit_coverage;
        // Peak bucket gets ~5x the average (interference)
        let peak = (hits_per_frame / N_BUCKETS as f32) * 5.0 * self.frame as f32;

        self.update_uniforms(queue, barrier_x, slit1_y, slit2_y, screen_x, peak.max(1.0));
    }

    fn slit_half_w(&self) -> f32 { (self.height as f32 * 0.030).max(12.0) }

    fn update_uniforms(&self, queue: &Queue, barrier_x: f32, slit1_y: f32, slit2_y: f32, screen_x: f32, peak: f32) {
        let u = build_uniforms(
            self.width, self.height,
            barrier_x, slit1_y, slit2_y,
            self.slit_gap,
            screen_x,
            self.running,
            self.frame,
            peak,
        );
        queue.write_buffer(&self.uniform_buf, 0, bytemuck::bytes_of(&u));
    }

    // ── Compute pass ─────────────────────────────────────────────────────────

    pub fn compute(&self, enc: &mut CommandEncoder) {
        if !self.running { return; }

        let mut pass = enc.begin_compute_pass(&ComputePassDescriptor {
            label:              Some("ds-compute"),
            timestamp_writes:   None,
        });
        pass.set_pipeline(&self.compute_pipeline);
        pass.set_bind_group(0, &self.compute_bg, &[]);
        // Ceiling division: (N_PARTICLES + 63) / 64
        let groups = N_PARTICLES.div_ceil(64);
        pass.dispatch_workgroups(groups, 1, 1);
    }

    // ── Render pass ──────────────────────────────────────────────────────────

    pub fn render(&self, enc: &mut CommandEncoder, view: &TextureView) {
        let mut pass = enc.begin_render_pass(&RenderPassDescriptor {
            label: Some("ds-render"),
            color_attachments: &[Some(RenderPassColorAttachment {
                view,
                resolve_target: None,
                ops: Operations {
                    load:  LoadOp::Clear(Color { r: 0.02, g: 0.02, b: 0.06, a: 1.0 }),
                    store: StoreOp::Store,
                },
            })],
            depth_stencil_attachment: None,
            timestamp_writes:         None,
            occlusion_query_set:      None,
        });

        pass.set_bind_group(0, &self.render_bg, &[]);

        // 1. Barrier (opaque) — draw a fullscreen quad, fragment discards non-barrier pixels
        pass.set_pipeline(&self.barrier_pipeline);
        pass.draw(0..6, 0..1);

        // 2. Slit glow (additive)
        pass.set_pipeline(&self.slit_glow_pipeline);
        pass.draw(0..6, 0..1);

        // 3. Detection screen line (alpha blend)
        pass.set_pipeline(&self.screen_line_pipeline);
        pass.draw(0..6, 0..1);

        // 4. Particles (additive, instanced)
        pass.set_pipeline(&self.particle_pipeline);
        pass.draw(0..6, 0..N_PARTICLES);

        // 5. Interference pattern bars (additive, instanced)
        pass.set_pipeline(&self.pattern_pipeline);
        pass.draw(0..6, 0..N_BUCKETS);
    }
}

// ─── Pure helper functions ───────────────────────────────────────────────────

/// Derive experiment geometry from window dimensions.
fn geometry(w: u32, h: u32) -> (f32, f32, f32) {
    let fw = w as f32;
    let fh = h as f32;
    let barrier_x = fw * 0.38;
    let screen_x  = fw * 0.78;
    // Slit gap proportional to height — chosen so the interference pattern fills the screen
    let slit_gap  = fh * 0.12;
    (barrier_x, screen_x, slit_gap)
}

fn build_uniforms(
    w:          u32,
    h:          u32,
    barrier_x:  f32,
    slit1_y:    f32,
    slit2_y:    f32,
    slit_gap:   f32,
    screen_x:   f32,
    running:    bool,
    frame:      u32,
    peak:       f32,
) -> DsUniforms {
    let fh = h as f32;
    let slit_half_w = (fh * 0.030).max(12.0);
    let _ = slit_gap; // used implicitly via slit1_y / slit2_y

    DsUniforms {
        width:          w as f32,
        height:         fh,
        barrier_x,
        slit1_y,
        slit2_y,
        slit_half_w,
        screen_x,
        source_x:       w as f32 * 0.05,
        source_y:       fh / 2.0,
        particle_speed: 4.5,
        frame,
        running:        running as u32,
        n_particles:    N_PARTICLES,
        n_buckets:      N_BUCKETS,
        expected_peak:  peak,
        _pad1:          0,
    }
}

/// Initialise N particles at the source with staggered velocities.
fn init_particles(w: u32, h: u32, n: u32) -> Vec<Particle> {
    let fw = w as f32;
    let fh = h as f32;
    let source_x = fw * 0.05;
    let source_y = fh / 2.0;

    // Simple LCG for deterministic init
    let mut rng = 0x12345678u32;
    let mut next = move || -> f32 {
        rng = rng.wrapping_mul(1664525).wrapping_add(1013904223);
        (rng >> 8) as f32 / (1u32 << 24) as f32
    };

    (0..n).map(|i| {
        let jitter = (next() * 2.0 - 1.0) * 5.0;
        let vy     = (next() * 2.0 - 1.0) * 0.35;
        Particle {
            pos:   [source_x, source_y + jitter],
            vel:   [4.5, vy],
            state: 0,
            seed:  i.wrapping_mul(2654435761).wrapping_add(rng),
        }
    }).collect()
}
