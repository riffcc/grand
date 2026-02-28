/*!
 * GUTOE — Gaia DR3 MVL Streaming Life-Map (wgpu)
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * 1.81 billion Gaia DR3 stars → GUTOE habitability model → GPU point cloud.
 * Stars glow green when habitable, red when filtered.
 * Additive blending makes dense regions bloom — the galaxy emerges in real time.
 *
 * Controls:
 *   Arrow keys — pan    +/- — zoom    R — reset    Q/Escape — quit
 *
 * Env vars:
 *   GAIA_MVL       — path to gaia_dr3.mvl (default /mnt/riffcastle/gaia_dr3.mvl)
 *   GAIA_CACHE_DB  — SQLite cache path (default /tmp/gaia_life_cache.sqlite)
 */

use anyhow::{Context, Result};
use blake3::Hasher as B3;
use bytemuck::{Pod, Zeroable};
use gutoe_physics::{
    habitability_score, infer_component_from_position, is_habitable, main_sequence_lifetime_gyr,
    GalacticLifeSeed,
};
use memmap2::MmapOptions;
use rusqlite::{params, Connection};
use std::{
    fs::File,
    path::PathBuf,
    sync::{
        mpsc::{self, Receiver, Sender, TryRecvError},
        Arc,
    },
    thread,
    time::{Duration, Instant},
};
use winit::{
    application::ApplicationHandler,
    dpi::PhysicalSize,
    event::{ElementState, KeyEvent, WindowEvent},
    event_loop::{ActiveEventLoop, EventLoop},
    keyboard::{Key, NamedKey},
    window::{Window, WindowId},
};

// ─── MVL column byte offsets (binary analysis of /mnt/riffcastle/gaia_dr3.mvl)
const NROWS: usize = 1_811_709_771;
// DOUBLE (f64): file_offset + 64 + row*8
const COL_PARALLAX_OFF: u64 = 235_328_265_152;
const COL_L_OFF: u64 = 1_046_921_704_640;
const COL_B_OFF: u64 = 1_061_415_382_880;
// FLOAT (f32): file_offset + 64 + row*4
const COL_TEFF_OFF: u64 = 1_152_000_878_592;
const COL_LOGG_OFF: u64 = 1_173_741_396_096;
const COL_MH_OFF: u64 = 1_195_481_913_600;
const COL_DIST_OFF: u64 = 1_217_222_431_104;

const PC_TO_LY: f64 = 3.261_56;

/// Max habitable stars in GPU buffer (green)
const MAX_HAB: usize = 3_000_000;
/// Max non-habitable stars in GPU buffer (dim red) — second section
const MAX_FIL: usize = 1_000_000;
/// Total GPU buffer size — Sol + galactic centre always at fixed indices
const MAX_STARS: usize = 2 + MAX_HAB + MAX_FIL;
/// Stars per CPU processing batch
const BATCH_SIZE: usize = 200_000;

// ─── GPU vertex type ─────────────────────────────────────────────────────────

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct StarVertex {
    x: f32, // galactic X (parsecs)
    y: f32, // galactic Y (parsecs)
    r: f32,
    g: f32,
    b: f32,
}

#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable)]
struct CameraUniform {
    center: [f32; 2], // galactic pc at screen centre
    zoom: f32,        // NDC units per parsec (multiply galactic → NDC)
    aspect: f32,      // width / height (for non-square windows)
}

// ─── WGSL shader ─────────────────────────────────────────────────────────────

const SHADER: &str = r#"
struct Camera {
    center : vec2<f32>,
    zoom   : f32,
    aspect : f32,
}
@group(0) @binding(0) var<uniform> cam: Camera;

struct VIn {
    @location(0) pos : vec2<f32>,
    @location(1) col : vec3<f32>,
}
struct VOut {
    @builtin(position) clip : vec4<f32>,
    @location(0)       col  : vec3<f32>,
}

@vertex fn vs_main(v: VIn) -> VOut {
    var out : VOut;
    let rel = v.pos - cam.center;
    let ndc = rel * cam.zoom;
    out.clip = vec4<f32>(ndc.x / cam.aspect, ndc.y, 0.0, 1.0);
    out.col  = v.col;
    return out;
}

@fragment fn fs_main(in: VOut) -> @location(0) vec4<f32> {
    return vec4<f32>(in.col, 1.0);
}
"#;

// ─── Data flowing from processor to renderer ──────────────────────────────────

#[derive(Clone)]
struct StarResult {
    x: f32,
    y: f32,
    habitable: bool,
}

#[derive(Clone, Default)]
struct Stats {
    processed: usize,
    skipped: usize,
    habitable: usize,
    cached: usize,
    rate: f64,
    elapsed: Duration,
    done: bool,
}

struct Batch {
    stars: Vec<StarResult>,
    stats: Stats,
}

enum Update {
    Data(Batch),
    Done(Stats),
}

// ─── MVL byte readers ────────────────────────────────────────────────────────

#[inline(always)]
fn rd_f64(m: &[u8], col: u64, row: usize) -> f64 {
    let off = (col + 64 + row as u64 * 8) as usize;
    f64::from_le_bytes(m[off..off + 8].try_into().unwrap())
}

#[inline(always)]
fn rd_f32(m: &[u8], col: u64, row: usize) -> f32 {
    let off = (col + 64 + row as u64 * 4) as usize;
    f32::from_le_bytes(m[off..off + 4].try_into().unwrap())
}

// ─── Physics helpers ─────────────────────────────────────────────────────────

fn estimate_mass(teff_k: f32) -> f64 {
    (teff_k as f64 / 5778.0).powi(2).clamp(0.08, 5.0)
}

fn estimate_age(mh: f32) -> f64 {
    (5.0_f64 * 10.0_f64.powf(-0.4 * mh as f64)).clamp(0.5, 13.0)
}

fn metallicity_z(mh: f32) -> f64 {
    0.014_2 * 10.0_f64.powf(mh as f64)
}

fn param_hash(l: f64, b: f64, parallax: f64, teff: f32, logg: f32, mh: f32) -> [u8; 32] {
    let mut h = B3::new();
    h.update(&l.to_le_bytes());
    h.update(&b.to_le_bytes());
    h.update(&parallax.to_le_bytes());
    h.update(&teff.to_le_bytes());
    h.update(&logg.to_le_bytes());
    h.update(&mh.to_le_bytes());
    *h.finalize().as_bytes()
}

// ─── SQLite helpers ───────────────────────────────────────────────────────────

fn db_init(c: &Connection) -> Result<()> {
    c.execute_batch(
        "PRAGMA journal_mode = WAL;
         PRAGMA synchronous  = NORMAL;
         CREATE TABLE IF NOT EXISTS star_fate (
             row_idx   INTEGER PRIMARY KEY,
             phash     BLOB    NOT NULL,
             habitable INTEGER NOT NULL,
             x_gal     REAL    NOT NULL,
             y_gal     REAL    NOT NULL
         );",
    )?;
    Ok(())
}

struct CacheEntry {
    row: usize,
    hash: [u8; 32],
    habitable: bool,
    x: f32,
    y: f32,
}

fn db_lookup(c: &Connection, row: usize, hash: &[u8]) -> Option<(bool, f32, f32)> {
    c.query_row(
        "SELECT phash, habitable, x_gal, y_gal FROM star_fate WHERE row_idx = ?1",
        params![row as i64],
        |r| {
            let stored: Vec<u8> = r.get(0)?;
            if stored.as_slice() == hash {
                Ok(Some((
                    r.get::<_, i64>(1)? != 0,
                    r.get::<_, f64>(2)? as f32,
                    r.get::<_, f64>(3)? as f32,
                )))
            } else {
                Ok(None)
            }
        },
    )
    .unwrap_or(None)
}

fn db_flush(c: &Connection, rows: &[CacheEntry]) -> Result<()> {
    if rows.is_empty() {
        return Ok(());
    }
    c.execute("BEGIN", [])?;
    for e in rows {
        c.execute(
            "INSERT OR REPLACE INTO star_fate (row_idx, phash, habitable, x_gal, y_gal)
             VALUES (?1, ?2, ?3, ?4, ?5)",
            params![
                e.row as i64,
                e.hash.to_vec(),
                e.habitable as i64,
                e.x as f64,
                e.y as f64,
            ],
        )?;
    }
    c.execute("COMMIT", [])?;
    Ok(())
}

// ─── Processor thread ────────────────────────────────────────────────────────

fn run_processor(mvl_path: PathBuf, db_path: PathBuf, tx: Sender<Update>) {
    let result: Result<()> = (|| {
        let file = File::open(&mvl_path).context("open MVL")?;
        let mmap = unsafe { MmapOptions::new().map(&file).context("mmap")? };
        let db = Connection::open(&db_path)?;
        db_init(&db)?;

        let mut stats = Stats::default();
        let t0 = Instant::now();
        let mut pending: Vec<CacheEntry> = Vec::with_capacity(4096);

        for batch_start in (0..NROWS).step_by(BATCH_SIZE) {
            let batch_end = (batch_start + BATCH_SIZE).min(NROWS);
            let mut batch_stars: Vec<StarResult> = Vec::new();

            for row in batch_start..batch_end {
                let teff = rd_f32(&mmap, COL_TEFF_OFF, row);
                if teff.is_nan() || teff < 2_800.0 || teff > 10_000.0 {
                    stats.skipped += 1;
                    continue;
                }
                let logg = rd_f32(&mmap, COL_LOGG_OFF, row);
                if logg.is_nan() || logg < 3.5 {
                    stats.skipped += 1;
                    continue;
                }

                let mh_raw = rd_f32(&mmap, COL_MH_OFF, row);
                let mh = if mh_raw.is_nan() { 0.0_f32 } else { mh_raw };
                let parallax = rd_f64(&mmap, COL_PARALLAX_OFF, row);
                let l = rd_f64(&mmap, COL_L_OFF, row);
                let b = rd_f64(&mmap, COL_B_OFF, row);
                let dist_raw = rd_f32(&mmap, COL_DIST_OFF, row);

                let hash = param_hash(l, b, parallax, teff, logg, mh);

                // Cache lookup
                if let Some((hab, x, y)) = db_lookup(&db, row, &hash) {
                    stats.cached += 1;
                    stats.processed += 1;
                    if hab {
                        stats.habitable += 1;
                    }
                    batch_stars.push(StarResult { x, y, habitable: hab });
                    continue;
                }

                // Compute
                let dist_pc = if !dist_raw.is_nan() && dist_raw > 0.0 {
                    dist_raw as f64
                } else if parallax > 0.5 {
                    1000.0 / parallax
                } else {
                    stats.skipped += 1;
                    continue;
                };

                let l_r = l.to_radians();
                let b_r = b.to_radians();
                let d_proj = dist_pc * b_r.cos();
                let x_gal = (d_proj * l_r.cos()) as f32;
                let y_gal = (d_proj * l_r.sin()) as f32;
                let z_gal = dist_pc * b_r.sin();

                let mass = estimate_mass(teff);
                let age = estimate_age(mh);
                let metal = metallicity_z(mh);
                let ms_life = main_sequence_lifetime_gyr(mass);
                let comp = infer_component_from_position(
                    x_gal as f64 * PC_TO_LY,
                    y_gal as f64 * PC_TO_LY,
                    z_gal * PC_TO_LY,
                );

                let seed = GalacticLifeSeed {
                    id: row as u64,
                    component: comp,
                    x_ly: x_gal as f64 * PC_TO_LY,
                    y_ly: y_gal as f64 * PC_TO_LY,
                    z_ly: z_gal * PC_TO_LY,
                    galactic_radius_ly: d_proj * PC_TO_LY,
                    mass_solar: mass,
                    age_gyr: age,
                    metallicity: metal,
                    main_sequence_lifetime_gyr: ms_life,
                };

                let score = habitability_score(seed);
                let hab = is_habitable(seed, score);

                pending.push(CacheEntry {
                    row,
                    hash,
                    habitable: hab,
                    x: x_gal,
                    y: y_gal,
                });

                stats.processed += 1;
                if hab {
                    stats.habitable += 1;
                }
                batch_stars.push(StarResult {
                    x: x_gal,
                    y: y_gal,
                    habitable: hab,
                });
            }

            if !pending.is_empty() {
                db_flush(&db, &pending).ok();
                pending.clear();
            }

            stats.elapsed = t0.elapsed();
            let seen = (stats.processed + stats.skipped) as f64;
            stats.rate = if stats.elapsed.as_secs_f64() > 0.0 {
                seen / stats.elapsed.as_secs_f64()
            } else {
                0.0
            };

            if tx.send(Update::Data(Batch { stars: batch_stars, stats: stats.clone() })).is_err() {
                break;
            }
        }

        stats.done = true;
        let _ = tx.send(Update::Done(stats));
        Ok(())
    })();

    if let Err(e) = result {
        eprintln!("Processor error: {e}");
    }
}

// ─── wgpu state ──────────────────────────────────────────────────────────────

struct Gpu {
    _instance: wgpu::Instance,
    surface: wgpu::Surface<'static>,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    pipeline: wgpu::RenderPipeline,
    vertex_buf: wgpu::Buffer,
    cam_buf: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    window: Arc<Window>,
}

impl Gpu {
    fn new(window: Arc<Window>) -> Self {
        pollster::block_on(Self::new_async(window))
    }

    async fn new_async(window: Arc<Window>) -> Self {
        let backends = match std::env::var("WGPU_BACKEND").as_deref() {
            Ok("vulkan") => wgpu::Backends::VULKAN,
            Ok("gl") => wgpu::Backends::GL,
            _ => wgpu::Backends::all(),
        };
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends,
            ..Default::default()
        });

        let surface = instance.create_surface(Arc::clone(&window)).unwrap();

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                force_fallback_adapter: false,
                compatible_surface: Some(&surface),
            })
            .await
            .expect("no suitable GPU adapter found");

        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor::default(), None)
            .await
            .expect("device request failed");

        let size = window.inner_size();
        let caps = surface.get_capabilities(&adapter);
        let format = caps
            .formats
            .iter()
            .copied()
            .find(|f| f.is_srgb())
            .unwrap_or(caps.formats[0]);

        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format,
            width: size.width.max(1),
            height: size.height.max(1),
            present_mode: wgpu::PresentMode::Fifo,
            alpha_mode: caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("gaia_shader"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });

        let cam_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("cam_buf"),
            size: std::mem::size_of::<CameraUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

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
                resource: cam_buf.as_entire_binding(),
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("pl"),
            bind_group_layouts: &[&bgl],
            push_constant_ranges: &[],
        });

        // Pre-allocate vertex buffer for MAX_STARS (Sol + GC + habitable + filtered)
        let vertex_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("star_buf"),
            size: (MAX_STARS * std::mem::size_of::<StarVertex>()) as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // Sol at index 0 (yellow), galactic centre at index 1 (magenta)
        let landmark_verts = [
            StarVertex { x: 0.0, y: 0.0, r: 1.0, g: 1.0, b: 0.0 },      // ☀ Sol
            StarVertex { x: 8122.0, y: 0.0, r: 0.9, g: 0.3, b: 1.0 },   // ⦿ GC
        ];
        queue.write_buffer(&vertex_buf, 0, bytemuck::cast_slice(&landmark_verts));

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("gaia_pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: "vs_main",
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<StarVertex>() as u64,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &[
                        wgpu::VertexAttribute {
                            format: wgpu::VertexFormat::Float32x2,
                            offset: 0,
                            shader_location: 0,
                        },
                        wgpu::VertexAttribute {
                            format: wgpu::VertexFormat::Float32x3,
                            offset: 8,
                            shader_location: 1,
                        },
                    ],
                }],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: "fs_main",
                targets: &[Some(wgpu::ColorTargetState {
                    format,
                    // Additive blending: dense regions bloom bright
                    blend: Some(wgpu::BlendState {
                        color: wgpu::BlendComponent {
                            src_factor: wgpu::BlendFactor::SrcAlpha,
                            dst_factor: wgpu::BlendFactor::One,
                            operation: wgpu::BlendOperation::Add,
                        },
                        alpha: wgpu::BlendComponent::OVER,
                    }),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::PointList,
                ..Default::default()
            },
            depth_stencil: None,
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
            cache: None,
        });

        Self {
            _instance: instance,
            surface,
            device,
            queue,
            config,
            pipeline,
            vertex_buf,
            cam_buf,
            bind_group,
            window,
        }
    }

    fn resize(&mut self, size: PhysicalSize<u32>) {
        if size.width > 0 && size.height > 0 {
            self.config.width = size.width;
            self.config.height = size.height;
            self.surface.configure(&self.device, &self.config);
        }
    }

    fn upload_stars(&self, stars: &[StarVertex], byte_offset: u64) {
        self.queue.write_buffer(&self.vertex_buf, byte_offset, bytemuck::cast_slice(stars));
    }

    fn render(&self, n_total: u32, cam: CameraUniform) {
        self.queue.write_buffer(&self.cam_buf, 0, bytemuck::bytes_of(&cam));

        let frame = match self.surface.get_current_texture() {
            Ok(f) => f,
            Err(_) => return,
        };
        let view = frame.texture.create_view(&Default::default());
        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("gaia_pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.0,
                            g: 0.0,
                            b: 0.015,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                ..Default::default()
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &self.bind_group, &[]);
            pass.set_vertex_buffer(0, self.vertex_buf.slice(..));
            pass.draw(0..n_total, 0..1);
        }
        self.queue.submit(std::iter::once(enc.finish()));
        frame.present();
        self.window.request_redraw();
    }
}

// ─── Application state ────────────────────────────────────────────────────────

struct AppState {
    // GPU-side star counts (2 fixed landmarks + accumulated stars)
    n_hab: usize,  // habitable stars uploaded (placed after the 2 landmarks)
    n_fil: usize,  // filtered stars uploaded  (placed after habitable section)
    // Stats
    stats: Stats,
    // Camera
    cam_center: [f32; 2],
    cam_zoom: f32,  // NDC per parsec
    done: bool,
    // Buffer upload tracking
    last_uploaded: usize, // total bytes written past the 2-landmark preamble
}

impl AppState {
    fn new() -> Self {
        Self {
            n_hab: 0,
            n_fil: 0,
            stats: Stats::default(),
            cam_center: [0.0, 0.0],
            cam_zoom: 1.0 / 25_000.0, // ±25 kpc fills the screen
            done: false,
            last_uploaded: 0,
        }
    }

    fn camera_uniform(&self, width: u32, height: u32) -> CameraUniform {
        CameraUniform {
            center: self.cam_center,
            zoom: self.cam_zoom,
            aspect: width as f32 / height.max(1) as f32,
        }
    }

    fn n_total(&self) -> u32 {
        // 2 landmarks + habitable + filtered
        (2 + self.n_hab + self.n_fil) as u32
    }
}

// ─── App (winit application handler) ─────────────────────────────────────────

struct App {
    gpu: Option<Gpu>,
    rx: Receiver<Update>,
    state: AppState,
    // Staging buffer: CPU-side vertices for next GPU upload
    upload_buf: Vec<StarVertex>,
}

impl App {
    fn new(rx: Receiver<Update>) -> Self {
        Self {
            gpu: None,
            rx,
            state: AppState::new(),
            upload_buf: Vec::new(),
        }
    }

    fn drain_channel(&mut self) {
        loop {
            match self.rx.try_recv() {
                Ok(Update::Data(batch)) => {
                    self.state.stats = batch.stats;
                    for s in batch.stars {
                        let v = if s.habitable {
                            if self.state.n_hab < MAX_HAB {
                                let v = StarVertex {
                                    x: s.x,
                                    y: s.y,
                                    r: 0.0,
                                    g: 1.0,
                                    b: 0.2,
                                };
                                self.state.n_hab += 1;
                                Some((v, false))
                            } else {
                                None
                            }
                        } else if self.state.n_fil < MAX_FIL {
                            let v = StarVertex {
                                x: s.x,
                                y: s.y,
                                r: 0.7,
                                g: 0.05,
                                b: 0.05,
                            };
                            self.state.n_fil += 1;
                            Some((v, true))
                        } else {
                            None
                        };
                        if let Some((vertex, is_fil)) = v {
                            let _ = (is_fil, vertex); // stored via staging below
                            self.upload_buf.push(vertex);
                        }
                    }
                }
                Ok(Update::Done(stats)) => {
                    self.state.stats = stats;
                    self.state.done = true;
                }
                Err(TryRecvError::Empty) => break,
                Err(TryRecvError::Disconnected) => {
                    self.state.done = true;
                    break;
                }
            }
        }
    }

    fn flush_upload(&mut self) {
        if self.upload_buf.is_empty() {
            return;
        }
        if let Some(gpu) = &self.gpu {
            // Write after the 2 landmark vertices
            let landmark_bytes = 2 * std::mem::size_of::<StarVertex>() as u64;
            let offset = landmark_bytes + self.state.last_uploaded as u64;
            gpu.upload_stars(&self.upload_buf, offset);
            self.state.last_uploaded +=
                self.upload_buf.len() * std::mem::size_of::<StarVertex>();
            self.upload_buf.clear();
        }
    }

    fn update_title(&self) {
        if let Some(gpu) = &self.gpu {
            let s = &self.state.stats;
            let hab_pct = if s.processed > 0 {
                s.habitable as f64 / s.processed as f64 * 100.0
            } else {
                0.0
            };
            let scan_pct = (s.processed + s.skipped) as f64 / NROWS as f64 * 100.0;
            let status = if self.state.done { "DONE" } else { "SCANNING" };
            gpu.window.set_title(&format!(
                "GUTOE Galactic Life Map | {status} {scan_pct:.1}% | \
                 Habitable: {} ({hab_pct:.1}%) | {:.0} stars/s | ☀=Sol  ⦿=GC  q=quit",
                s.habitable, s.rate,
            ));
        }
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, el: &ActiveEventLoop) {
        let win = Arc::new(
            el.create_window(
                winit::window::WindowAttributes::default()
                    .with_title("GUTOE Galactic Life Map — Loading…")
                    .with_inner_size(winit::dpi::LogicalSize::new(1400u32, 900u32)),
            )
            .unwrap(),
        );
        self.gpu = Some(Gpu::new(Arc::clone(&win)));
        win.request_redraw();
    }

    fn window_event(&mut self, el: &ActiveEventLoop, _id: WindowId, event: WindowEvent) {
        match event {
            WindowEvent::CloseRequested => el.exit(),

            WindowEvent::KeyboardInput {
                event: KeyEvent { logical_key, state: ElementState::Pressed, .. },
                ..
            } => {
                let pan_step = 1_000.0_f32 / (self.state.cam_zoom * 25_000.0);
                match logical_key {
                    Key::Named(NamedKey::Escape) => el.exit(),
                    Key::Character(c) if c.as_str() == "q" || c.as_str() == "Q" => el.exit(),
                    Key::Character(c) if c.as_str() == "r" || c.as_str() == "R" => {
                        self.state.cam_center = [0.0, 0.0];
                        self.state.cam_zoom = 1.0 / 25_000.0;
                    }
                    Key::Character(c) if c.as_str() == "+" || c.as_str() == "=" => {
                        self.state.cam_zoom *= 1.5;
                    }
                    Key::Character(c) if c.as_str() == "-" => {
                        self.state.cam_zoom /= 1.5;
                    }
                    Key::Named(NamedKey::ArrowLeft) => {
                        self.state.cam_center[0] -= pan_step;
                    }
                    Key::Named(NamedKey::ArrowRight) => {
                        self.state.cam_center[0] += pan_step;
                    }
                    Key::Named(NamedKey::ArrowUp) => {
                        self.state.cam_center[1] += pan_step;
                    }
                    Key::Named(NamedKey::ArrowDown) => {
                        self.state.cam_center[1] -= pan_step;
                    }
                    _ => {}
                }
                if let Some(gpu) = &self.gpu {
                    gpu.window.request_redraw();
                }
            }

            WindowEvent::Resized(size) => {
                if let Some(gpu) = &mut self.gpu {
                    gpu.resize(size);
                }
            }

            WindowEvent::RedrawRequested => {
                self.drain_channel();
                self.flush_upload();
                self.update_title();

                if let Some(gpu) = &self.gpu {
                    let cam =
                        self.state.camera_uniform(gpu.config.width, gpu.config.height);
                    gpu.render(self.state.n_total(), cam);
                }
            }

            _ => {}
        }
    }

    fn about_to_wait(&mut self, _el: &ActiveEventLoop) {
        if let Some(gpu) = &self.gpu {
            gpu.window.request_redraw();
        }
    }
}

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() -> Result<()> {
    env_logger::init();

    let mvl_path = PathBuf::from(
        std::env::var("GAIA_MVL").unwrap_or_else(|_| "/mnt/riffcastle/gaia_dr3.mvl".into()),
    );
    let db_path = PathBuf::from(
        std::env::var("GAIA_CACHE_DB").unwrap_or_else(|_| "/tmp/gaia_life_cache.sqlite".into()),
    );

    if !mvl_path.exists() {
        anyhow::bail!(
            "Gaia DR3 MVL not found at {}\nSet GAIA_MVL env var to override.",
            mvl_path.display()
        );
    }

    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  GUTOE Galactic Life Map — Gaia DR3 × 1.81 billion stars");
    println!("  Green = habitable   Red = filtered   Yellow ☀ = Sol");
    println!("  Purple ⦿ = Galactic centre (~8.1 kpc)");
    println!("  Arrow keys: pan   +/-: zoom   R: reset   Q: quit");
    println!("  GAIA_MVL = {}", mvl_path.display());
    println!("  GAIA_CACHE_DB = {}", db_path.display());
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let (tx, rx) = mpsc::channel();
    {
        let mvl = mvl_path.clone();
        let db = db_path.clone();
        thread::spawn(move || run_processor(mvl, db, tx));
    }

    let ev = EventLoop::new().context("create event loop")?;
    let mut app = App::new(rx);
    ev.run_app(&mut app).context("event loop")?;
    Ok(())
}
