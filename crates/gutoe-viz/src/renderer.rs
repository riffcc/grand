//! GUTOE wgpu renderer — instanced hex grid, 4×3 layout of 12 layers.
//!
//! Vertex buffer  (static): 18 verts, local pixel-space offsets for one pointy-top hex
//! Instance buffer (dynamic): 1728 × InstanceData { world_pos: [f32;2], color: [f32;4] }
//! Uniform buffer: ScreenUniforms { screen_size: [f32;2], _pad: [f32;2] }
//!
//! WGSL pixel → NDC: ndc.x = px.x/W*2-1, ndc.y = 1-px.y/H*2

use wgpu::util::DeviceExt;
use bytemuck::{Pod, Zeroable};

// ── Layout constants ───────────────────────────────────────────────────────────

pub const HEX_SIZE: f32    = 10.0;                              // pointy-top radius px
pub const HEX_W: f32       = 1.732_050_8 * HEX_SIZE;           // √3 × size
pub const HEX_H: f32       = 2.0 * HEX_SIZE;
pub const COL_STEP: f32    = HEX_W;
pub const ROW_STEP: f32    = HEX_H * 0.75;
pub const ODD_OFFSET: f32  = HEX_W * 0.5;

// Panel layout: 4 columns × 3 rows of layers, 15px gap
pub const GAP: f32         = 15.0;
pub const PANEL_W: f32     = 280.0;
pub const PANEL_H: f32     = 215.0;

// ── Vertex / instance data ─────────────────────────────────────────────────────

/// Per-vertex data: local offset in pixel space from hex center
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct HexVertex {
    pub local_pos: [f32; 2],
}

/// Per-instance data: world-space center + RGBA color
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct InstanceData {
    pub world_pos: [f32; 2],
    pub color:     [f32; 4],
}

/// Screen uniforms for NDC conversion
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct ScreenUniforms {
    pub screen_size: [f32; 2],
    pub _pad: [f32; 2],
}

// ── Hex vertex geometry (pointy-top) ──────────────────────────────────────────

/// Build the 18 local vertices (6 triangles) for one pointy-top hex.
/// Pointy-top: vertex i at angle (30 + 60*i)°.
pub fn hex_vertices() -> Vec<HexVertex> {
    let mut verts = Vec::with_capacity(18);
    let center = [0.0f32, 0.0f32];
    for i in 0..6u32 {
        let a0 = std::f32::consts::PI / 180.0 * (30.0 + 60.0 * i as f32);
        let a1 = std::f32::consts::PI / 180.0 * (30.0 + 60.0 * (i + 1) as f32);
        let v0 = [HEX_SIZE * a0.cos(), HEX_SIZE * a0.sin()];
        let v1 = [HEX_SIZE * a1.cos(), HEX_SIZE * a1.sin()];
        verts.push(HexVertex { local_pos: center });
        verts.push(HexVertex { local_pos: v0 });
        verts.push(HexVertex { local_pos: v1 });
    }
    verts
}

// ── Color scheme ───────────────────────────────────────────────────────────────

pub fn state_color(state: u8) -> [f32; 4] {
    match state {
        0       => [0.06, 0.06, 0.10, 1.0], // VOID — dark blue-black
        2       => [0.0,  1.0,  1.0,  1.0], // LEPTON γ⁰ — bright cyan
        3 | 5 | 9  => [0.95, 0.45, 0.1,  1.0], // Z3 orbit 0 — orange
        4 | 6 | 10 => [0.8,  0.2,  0.3,  1.0], // Z3 orbit 1 — crimson
        7 | 11 | 13 => [0.5, 0.2,  0.9,  1.0], // Z3 orbit 2 — purple
        8 | 12 | 14 => [0.2, 0.6,  0.9,  1.0], // Z3 orbit 3 — blue
        1       => [1.0,  1.0,  0.8,  1.0], // grade-0 scalar — cream
        16      => [1.0,  0.2,  0.2,  1.0], // grade-4 pseudoscalar — red
        s => {
            // grade-2 or grade-3
            let grade = (s - 1).count_ones();
            if grade == 2 {
                [0.3, 0.5, 0.95, 1.0] // bivector — steel blue
            } else {
                [0.7, 0.3, 0.85, 1.0] // trivector — violet
            }
        }
    }
}

/// World-space center position of hex (r, c) in layer `layer_idx`.
pub fn hex_world_pos(r: usize, c: usize, layer_idx: usize) -> [f32; 2] {
    let col_in_grid = (layer_idx % 4) as f32;
    let row_in_grid = (layer_idx / 4) as f32;
    let panel_ox = GAP + col_in_grid * (PANEL_W + GAP);
    let panel_oy = GAP + row_in_grid * (PANEL_H + GAP);

    let x = panel_ox + HEX_SIZE + c as f32 * COL_STEP
        + if r % 2 == 1 { ODD_OFFSET } else { 0.0 };
    let y = panel_oy + HEX_SIZE + r as f32 * ROW_STEP;
    [x, y]
}

// ── WGSL shader ────────────────────────────────────────────────────────────────

pub const SHADER_SRC: &str = r#"
struct ScreenUniforms {
    screen_size: vec2<f32>,
    _pad: vec2<f32>,
};

@group(0) @binding(0)
var<uniform> screen: ScreenUniforms;

struct VertexInput {
    @location(0) local_pos: vec2<f32>,
    @location(1) world_pos: vec2<f32>,
    @location(2) color:     vec4<f32>,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
};

@vertex
fn vs_main(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    let px = in.local_pos + in.world_pos;
    let ndc = vec2<f32>(
        px.x / screen.screen_size.x * 2.0 - 1.0,
        1.0 - px.y / screen.screen_size.y * 2.0
    );
    out.clip_position = vec4<f32>(ndc, 0.0, 1.0);
    out.color = in.color;
    return out;
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    return in.color;
}
"#;

// ── Renderer ───────────────────────────────────────────────────────────────────

pub struct Renderer {
    pub surface:        wgpu::Surface<'static>,
    pub device:         wgpu::Device,
    pub queue:          wgpu::Queue,
    pub config:         wgpu::SurfaceConfiguration,
    pub pipeline:       wgpu::RenderPipeline,
    pub vertex_buf:     wgpu::Buffer,
    pub instance_buf:   wgpu::Buffer,
    pub uniform_buf:    wgpu::Buffer,
    pub bind_group:     wgpu::BindGroup,
    pub n_vertices:     u32,
    pub n_instances:    u32,
    pub scale_factor:   f64,
}

impl Renderer {
    pub async fn new(
        instance: &wgpu::Instance,
        surface: wgpu::Surface<'static>,
        width: u32,
        height: u32,
        scale_factor: f64,
        n_hex: usize,
    ) -> Self {
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .expect("No suitable GPU adapter found");

        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor {
                label: Some("gutoe-device"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::default(),
                memory_hints: wgpu::MemoryHints::default(),
            }, None)
            .await
            .expect("Failed to create device");

        let surface_caps = surface.get_capabilities(&adapter);
        // Prefer non-sRGB: our colors are already in display/sRGB space.
        // Using an sRGB render target would apply gamma correction a second
        // time, washing everything out to beige/brown.
        let format = surface_caps.formats.iter()
            .find(|f| !f.is_srgb())
            .copied()
            .unwrap_or(surface_caps.formats[0]);

        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format,
            width,
            height,
            present_mode: wgpu::PresentMode::Fifo,
            alpha_mode: surface_caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        // Vertex buffer
        let verts = hex_vertices();
        let n_vertices = verts.len() as u32;
        let vertex_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("hex-verts"),
            contents: bytemuck::cast_slice(&verts),
            usage: wgpu::BufferUsages::VERTEX,
        });

        // Instance buffer (initially zeroed)
        let n_instances = n_hex as u32;
        let instance_data: Vec<InstanceData> = vec![InstanceData {
            world_pos: [0.0; 2],
            color: [0.0; 4],
        }; n_hex];
        let instance_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("hex-instances"),
            contents: bytemuck::cast_slice(&instance_data),
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        });

        // Uniform buffer — store LOGICAL size so hex world_pos (in logical px)
        // maps correctly to NDC regardless of Retina DPI.
        let logical_w = width  as f32 / scale_factor as f32;
        let logical_h = height as f32 / scale_factor as f32;
        let uniforms = ScreenUniforms {
            screen_size: [logical_w, logical_h],
            _pad: [0.0; 2],
        };
        let uniform_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("screen-uniforms"),
            contents: bytemuck::bytes_of(&uniforms),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        // Bind group layout & bind group
        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("bgl"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("bg"),
            layout: &bgl,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buf.as_entire_binding(),
            }],
        });

        // Shader
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("gutoe-shader"),
            source: wgpu::ShaderSource::Wgsl(SHADER_SRC.into()),
        });

        // Pipeline layout
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("pl"),
            bind_group_layouts: &[&bgl],
            push_constant_ranges: &[],
        });

        // Vertex buffer layouts
        let vertex_buf_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<HexVertex>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[wgpu::VertexAttribute {
                offset: 0,
                shader_location: 0,
                format: wgpu::VertexFormat::Float32x2,
            }],
        };
        let instance_buf_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<InstanceData>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &[
                wgpu::VertexAttribute {
                    offset: 0,
                    shader_location: 1,
                    format: wgpu::VertexFormat::Float32x2,
                },
                wgpu::VertexAttribute {
                    offset: 8,
                    shader_location: 2,
                    format: wgpu::VertexFormat::Float32x4,
                },
            ],
        };

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("gutoe-pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: "vs_main",
                buffers: &[vertex_buf_layout, instance_buf_layout],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: "fs_main",
                targets: &[Some(wgpu::ColorTargetState {
                    format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                ..Default::default()
            },
            depth_stencil: None,
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
            cache: None,
        });

        Renderer {
            surface, device, queue, config, pipeline,
            vertex_buf, instance_buf, uniform_buf, bind_group,
            n_vertices, n_instances, scale_factor,
        }
    }

    /// Upload new instance data to GPU
    pub fn update_instances(&self, instances: &[InstanceData]) {
        self.queue.write_buffer(
            &self.instance_buf,
            0,
            bytemuck::cast_slice(instances),
        );
    }

    pub fn render(&self) -> Result<(), wgpu::SurfaceError> {
        let output = self.surface.get_current_texture()?;
        let view = output.texture.create_view(&wgpu::TextureViewDescriptor::default());
        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("render-encoder"),
        });
        {
            let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("rpass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.02, g: 0.02, b: 0.04, a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                occlusion_query_set: None,
                timestamp_writes: None,
            });
            rpass.set_pipeline(&self.pipeline);
            rpass.set_bind_group(0, &self.bind_group, &[]);
            rpass.set_vertex_buffer(0, self.vertex_buf.slice(..));
            rpass.set_vertex_buffer(1, self.instance_buf.slice(..));
            rpass.draw(0..self.n_vertices, 0..self.n_instances);
        }
        self.queue.submit(std::iter::once(encoder.finish()));
        output.present();
        Ok(())
    }

    pub fn resize(&mut self, width: u32, height: u32, scale_factor: f64) {
        if width == 0 || height == 0 { return; }
        self.scale_factor = scale_factor;
        self.config.width = width;
        self.config.height = height;
        self.surface.configure(&self.device, &self.config);
        let logical_w = width  as f32 / scale_factor as f32;
        let logical_h = height as f32 / scale_factor as f32;
        let uniforms = ScreenUniforms {
            screen_size: [logical_w, logical_h],
            _pad: [0.0; 2],
        };
        self.queue.write_buffer(&self.uniform_buf, 0, bytemuck::bytes_of(&uniforms));
    }
}
