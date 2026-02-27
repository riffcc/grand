//! GUTOE wgpu renderer — instanced hex grid in 3D, perspective camera.
//!
//! 12 layers of 12×12 hex grids stacked along the Z axis.
//! Camera orbits around the centre of the stack.
//!
//! InstanceData:  world_pos: [f32;3], _pad: f32, color: [f32;4]  (32 bytes)
//! CameraUniforms: view_proj: mat4x4 (64 bytes, column-major)

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

// ── 3D layout constants ────────────────────────────────────────────────────────

pub const HEX_SIZE: f32 = 8.0;
pub const HEX_W: f32 = 13.856_407; // √3 × HEX_SIZE
pub const COL_STEP: f32 = HEX_W;
pub const ROW_STEP: f32 = 12.0; // HEX_SIZE × 1.5
pub const ODD_OFFSET: f32 = 6.928_203; // HEX_W × 0.5
pub const LAYER_SPACING: f32 = 22.0;

// ── GPU data types ─────────────────────────────────────────────────────────────

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct HexVertex {
    pub local_pos: [f32; 2],
}

/// Per-instance: 3D world-space centre + RGBA color.
/// _pad aligns `color` to 16-byte boundary (required by WGSL vec4 alignment).
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct InstanceData {
    pub world_pos: [f32; 3],
    pub _pad: f32,
    pub color: [f32; 4],
}

/// Column-major 4×4 matrix for view-projection uniform.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct CameraUniforms {
    pub view_proj: [[f32; 4]; 4],
}

// ── Hex geometry (pointy-top) ──────────────────────────────────────────────────

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
        0 => [0.06, 0.06, 0.10, 1.0],        // VOID
        2 => [0.0, 1.0, 1.0, 1.0],           // LEPTON γ⁰ — cyan
        3 | 5 | 9 => [0.95, 0.45, 0.1, 1.0], // Z3 orbit 0 — orange
        4 | 6 | 10 => [0.8, 0.2, 0.3, 1.0],  // Z3 orbit 1 — crimson
        7 | 11 | 13 => [0.5, 0.2, 0.9, 1.0], // Z3 orbit 2 — purple
        8 | 12 | 14 => [0.2, 0.6, 0.9, 1.0], // Z3 orbit 3 — blue
        1 => [1.0, 1.0, 0.8, 1.0],           // grade-0 scalar — cream
        16 => [1.0, 0.2, 0.2, 1.0],          // grade-4 pseudoscalar — red
        s => {
            if (s - 1).count_ones() == 2 {
                [0.3, 0.5, 0.95, 1.0] // bivector — steel blue
            } else {
                [0.7, 0.3, 0.85, 1.0] // trivector — violet
            }
        }
    }
}

/// 3D world-space centre for hex at grid position (r, c) in layer `layer`.
/// Y increases upward (row 0 = top, row hex_rows-1 = bottom).
pub fn hex_world_pos(r: usize, c: usize, layer: usize, hex_rows: usize) -> [f32; 3] {
    let x = c as f32 * COL_STEP + if r % 2 == 1 { ODD_OFFSET } else { 0.0 };
    let y = (hex_rows as f32 - 1.0 - r as f32) * ROW_STEP;
    let z = layer as f32 * LAYER_SPACING;
    [x, y, z]
}

// ── WGSL shader ────────────────────────────────────────────────────────────────

pub const SHADER_SRC: &str = r#"
struct CameraUniforms {
    view_proj: mat4x4<f32>,
};

@group(0) @binding(0)
var<uniform> camera: CameraUniforms;

struct VertexInput {
    @location(0) local_pos: vec2<f32>,
    @location(1) world_pos: vec3<f32>,
    @location(2) color:     vec4<f32>,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
};

@vertex
fn vs_main(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    // Hex faces lie in the XY plane; local_pos gives the XY offset from centre.
    let world = vec4<f32>(
        in.local_pos.x + in.world_pos.x,
        in.local_pos.y + in.world_pos.y,
        in.world_pos.z,
        1.0,
    );
    out.clip_position = camera.view_proj * world;
    out.color = in.color;
    return out;
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    return in.color;
}
"#;

// ── Depth texture helper ───────────────────────────────────────────────────────

fn make_depth_view(device: &wgpu::Device, width: u32, height: u32) -> wgpu::TextureView {
    device
        .create_texture(&wgpu::TextureDescriptor {
            label: Some("depth"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        })
        .create_view(&wgpu::TextureViewDescriptor::default())
}

// ── Renderer ───────────────────────────────────────────────────────────────────

pub struct Renderer {
    pub surface: wgpu::Surface<'static>,
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
    pub config: wgpu::SurfaceConfiguration,
    pub pipeline: wgpu::RenderPipeline,
    pub vertex_buf: wgpu::Buffer,
    pub instance_buf: wgpu::Buffer,
    pub uniform_buf: wgpu::Buffer,
    pub bind_group: wgpu::BindGroup,
    pub depth_view: wgpu::TextureView,
    pub n_vertices: u32,
    pub n_instances: u32,
}

impl Renderer {
    pub async fn new(
        instance: &wgpu::Instance,
        surface: wgpu::Surface<'static>,
        width: u32,
        height: u32,
        n_hex: usize,
    ) -> Self {
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .expect("No suitable GPU adapter");

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("gutoe"),
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::default(),
                    memory_hints: wgpu::MemoryHints::default(),
                },
                None,
            )
            .await
            .expect("Device creation failed");

        let caps = surface.get_capabilities(&adapter);
        // Non-sRGB: our colors are already in display space.
        let format = caps
            .formats
            .iter()
            .find(|f| !f.is_srgb())
            .copied()
            .unwrap_or(caps.formats[0]);

        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format,
            width,
            height,
            present_mode: wgpu::PresentMode::Fifo,
            alpha_mode: caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        // ── Vertex buffer ──────────────────────────────────────────────────────
        let verts = hex_vertices();
        let n_vertices = verts.len() as u32;
        let vertex_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("hex-verts"),
            contents: bytemuck::cast_slice(&verts),
            usage: wgpu::BufferUsages::VERTEX,
        });

        // ── Instance buffer ────────────────────────────────────────────────────
        let n_instances = n_hex as u32;
        let init_inst: Vec<InstanceData> = vec![
            InstanceData {
                world_pos: [0.0; 3],
                _pad: 0.0,
                color: [0.0; 4],
            };
            n_hex
        ];
        let instance_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("hex-instances"),
            contents: bytemuck::cast_slice(&init_inst),
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        });

        // ── Camera uniform buffer (identity initially) ─────────────────────────
        let identity = CameraUniforms {
            view_proj: identity4(),
        };
        let uniform_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("camera-uniforms"),
            contents: bytemuck::bytes_of(&identity),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        // ── Bind group ─────────────────────────────────────────────────────────
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

        // ── Shader ─────────────────────────────────────────────────────────────
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("gutoe-3d"),
            source: wgpu::ShaderSource::Wgsl(SHADER_SRC.into()),
        });

        // ── Pipeline ───────────────────────────────────────────────────────────
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("pl"),
            bind_group_layouts: &[&bgl],
            push_constant_ranges: &[],
        });

        let vtx_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<HexVertex>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[wgpu::VertexAttribute {
                offset: 0,
                shader_location: 0,
                format: wgpu::VertexFormat::Float32x2,
            }],
        };
        let inst_layout = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<InstanceData>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &[
                wgpu::VertexAttribute {
                    offset: 0,
                    shader_location: 1,
                    format: wgpu::VertexFormat::Float32x3,
                },
                wgpu::VertexAttribute {
                    offset: 16, // after world_pos[3] + _pad = 4*f32
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
                buffers: &[vtx_layout, inst_layout],
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
                cull_mode: None,
                ..Default::default()
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::Less,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
            cache: None,
        });

        let depth_view = make_depth_view(&device, width, height);

        Renderer {
            surface,
            device,
            queue,
            config,
            pipeline,
            vertex_buf,
            instance_buf,
            uniform_buf,
            bind_group,
            depth_view,
            n_vertices,
            n_instances,
        }
    }

    pub fn update_instances(&self, data: &[InstanceData]) {
        self.queue
            .write_buffer(&self.instance_buf, 0, bytemuck::cast_slice(data));
    }

    pub fn update_camera(&self, uniforms: &CameraUniforms) {
        self.queue
            .write_buffer(&self.uniform_buf, 0, bytemuck::bytes_of(uniforms));
    }

    pub fn render(&self) -> Result<(), wgpu::SurfaceError> {
        let output = self.surface.get_current_texture()?;
        let view = output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());
        let mut enc = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("enc") });
        {
            let mut rpass = enc.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("rpass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.02,
                            g: 0.02,
                            b: 0.04,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                occlusion_query_set: None,
                timestamp_writes: None,
            });
            rpass.set_pipeline(&self.pipeline);
            rpass.set_bind_group(0, &self.bind_group, &[]);
            rpass.set_vertex_buffer(0, self.vertex_buf.slice(..));
            rpass.set_vertex_buffer(1, self.instance_buf.slice(..));
            rpass.draw(0..self.n_vertices, 0..self.n_instances);
        }
        self.queue.submit(std::iter::once(enc.finish()));
        output.present();
        Ok(())
    }

    pub fn resize(&mut self, width: u32, height: u32) {
        if width == 0 || height == 0 {
            return;
        }
        self.config.width = width;
        self.config.height = height;
        self.surface.configure(&self.device, &self.config);
        self.depth_view = make_depth_view(&self.device, width, height);
    }
}

// ── Matrix helpers (column-major, matching WGSL mat4x4) ───────────────────────

pub fn identity4() -> [[f32; 4]; 4] {
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

/// Column-major 4×4 matrix multiply: C = A × B
pub fn mat4_mul(a: &[[f32; 4]; 4], b: &[[f32; 4]; 4]) -> [[f32; 4]; 4] {
    let mut c = [[0.0f32; 4]; 4];
    for col in 0..4 {
        for row in 0..4 {
            c[col][row] = (0..4).map(|k| a[k][row] * b[col][k]).sum();
        }
    }
    c
}

/// Perspective projection (right-hand, Z in [0, 1] per WebGPU convention).
/// `fov_y` in radians.
pub fn perspective(fov_y: f32, aspect: f32, near: f32, far: f32) -> [[f32; 4]; 4] {
    let f = 1.0 / (fov_y * 0.5).tan();
    let z_a = -far / (far - near);
    let z_b = -near * far / (far - near);
    [
        [f / aspect, 0.0, 0.0, 0.0],
        [0.0, f, 0.0, 0.0],
        [0.0, 0.0, z_a, -1.0],
        [0.0, 0.0, z_b, 0.0],
    ]
}

fn dot3(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm3(a: [f32; 3]) -> [f32; 3] {
    let l = dot3(a, a).sqrt();
    [a[0] / l, a[1] / l, a[2] / l]
}

/// LookAt view matrix (right-hand, column-major).
pub fn look_at(eye: [f32; 3], center: [f32; 3], up: [f32; 3]) -> [[f32; 4]; 4] {
    let f = norm3(sub3(center, eye)); // forward
    let r = norm3(cross3(f, up)); // right
    let u = cross3(r, f); // corrected up

    let tx = -dot3(r, eye);
    let ty = -dot3(u, eye);
    let tz = dot3(f, eye);
    [
        [r[0], u[0], -f[0], 0.0],
        [r[1], u[1], -f[1], 0.0],
        [r[2], u[2], -f[2], 0.0],
        [tx, ty, tz, 1.0],
    ]
}
