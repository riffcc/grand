/*!
 * GUTOE Viz - GPU Quantum Experiment Visualizer
 * Copyright (C) 2026  Wings
 *
 * AGPL-3.0-or-later
 */

//! Entry point.
//!
//! Usage:
//!   gutoe-viz                          # interactive window
//!   gutoe-viz --headless [--frames N] [--output DIR]

mod double_slit;
mod gpu;

use std::sync::Arc;
use winit::{
    application::ApplicationHandler,
    dpi::LogicalSize,
    event::{ElementState, KeyEvent, WindowEvent},
    event_loop::{ActiveEventLoop, EventLoop},
    keyboard::{KeyCode, PhysicalKey},
    window::{Window, WindowId},
};

fn main() {
    env_logger::init();

    let args: Vec<String> = std::env::args().collect();
    if args.iter().any(|a| a == "--headless") {
        headless::run(&args);
    } else {
        interactive();
    }
}

// ─── Interactive mode ────────────────────────────────────────────────────────

fn interactive() {
    let event_loop = EventLoop::new().expect("event loop");
    let mut app = GutoeApp::default();
    event_loop.run_app(&mut app).expect("event loop run");
}

#[derive(Default)]
struct GutoeApp {
    window:  Option<Arc<Window>>,
    gpu:     Option<gpu::GpuCtx>,
    ds:      Option<double_slit::DoubleSlitState>,
}

impl ApplicationHandler for GutoeApp {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.window.is_some() { return; } // already initialised

        let window = Arc::new(
            event_loop
                .create_window(
                    Window::default_attributes()
                        .with_title("GUTOE  ·  Double Slit  ·  [Space] pause  [R] reset  [Q] quit")
                        .with_inner_size(LogicalSize::new(1280u32, 720u32))
                        .with_resizable(true),
                )
                .expect("window creation"),
        );

        let gpu = pollster::block_on(gpu::GpuCtx::new(Arc::clone(&window)));
        let (w, h) = gpu.size();
        let ds = double_slit::DoubleSlitState::new(&gpu.device, &gpu.queue, gpu.format(), w, h);

        self.window = Some(window);
        self.gpu    = Some(gpu);
        self.ds     = Some(ds);
    }

    fn window_event(
        &mut self,
        event_loop: &ActiveEventLoop,
        _id: WindowId,
        event: WindowEvent,
    ) {
        match event {
            WindowEvent::CloseRequested => event_loop.exit(),

            WindowEvent::KeyboardInput {
                event: KeyEvent {
                    physical_key: PhysicalKey::Code(code),
                    state: ElementState::Pressed,
                    ..
                },
                ..
            } => match code {
                KeyCode::Escape | KeyCode::KeyQ => event_loop.exit(),

                KeyCode::Space => {
                    if let Some(ds) = &mut self.ds {
                        ds.toggle_running();
                        let status = if ds.is_running() { "running" } else { "paused" };
                        log::info!("Simulation {status}");
                    }
                }

                KeyCode::KeyR => {
                    if let (Some(ds), Some(gpu)) = (self.ds.as_mut(), self.gpu.as_ref()) {
                        ds.reset(&gpu.device, &gpu.queue);
                        log::info!("Simulation reset");
                    }
                }

                _ => {}
            },

            WindowEvent::Resized(size) => {
                if let (Some(gpu), Some(ds)) = (self.gpu.as_mut(), self.ds.as_mut()) {
                    gpu.resize(size.width, size.height);
                    ds.resize(size.width, size.height, &gpu.device, &gpu.queue);
                }
            }

            WindowEvent::RedrawRequested => {
                self.render();
            }

            _ => {}
        }
    }

    fn about_to_wait(&mut self, _: &ActiveEventLoop) {
        if let Some(w) = &self.window {
            w.request_redraw();
        }
    }
}

impl GutoeApp {
    fn render(&mut self) {
        let (Some(gpu), Some(ds)) = (self.gpu.as_mut(), self.ds.as_mut()) else { return };

        let output = match gpu.surface.get_current_texture() {
            Ok(t) => t,
            Err(SurfaceError::Outdated | SurfaceError::Lost) => {
                gpu.reconfigure();
                return;
            }
            Err(e) => {
                log::warn!("surface error: {e}");
                return;
            }
        };

        let view = output.texture.create_view(&Default::default());
        let mut enc = gpu.device.create_command_encoder(&Default::default());

        ds.update(&gpu.queue);
        ds.compute(&mut enc);
        ds.render(&mut enc, &view);

        gpu.queue.submit([enc.finish()]);
        output.present();
    }
}

use wgpu::SurfaceError;

// ─── Headless mode ───────────────────────────────────────────────────────────

mod headless {
    use super::*;
    use image::{ImageBuffer, Rgba};
    use std::path::PathBuf;
    use wgpu::*;

    pub fn run(args: &[String]) {
        let n_frames: u32 = args
            .windows(2)
            .find(|w| w[0] == "--frames")
            .and_then(|w| w[1].parse().ok())
            .unwrap_or(300);

        let output_dir: PathBuf = args
            .windows(2)
            .find(|w| w[0] == "--output")
            .map(|w| PathBuf::from(&w[1]))
            .unwrap_or_else(|| PathBuf::from("gutoe-frames"));

        std::fs::create_dir_all(&output_dir).expect("create output dir");

        let ctx = pollster::block_on(gpu::HeadlessCtx::new());
        let w = 1280u32;
        let h = 720u32;
        let format = TextureFormat::Rgba8UnormSrgb;

        // Offscreen render texture
        let tex = ctx.device.create_texture(&TextureDescriptor {
            label:           Some("hl-target"),
            size:            Extent3d { width: w, height: h, depth_or_array_layers: 1 },
            mip_level_count: 1,
            sample_count:    1,
            dimension:       TextureDimension::D2,
            format,
            usage:           TextureUsages::RENDER_ATTACHMENT | TextureUsages::COPY_SRC,
            view_formats:    &[],
        });
        let view = tex.create_view(&Default::default());

        // Readback buffer: each pixel = 4 bytes (RGBA8)
        let bytes_per_row = align256(w * 4);
        let readback = ctx.device.create_buffer(&BufferDescriptor {
            label:              Some("hl-readback"),
            size:               (bytes_per_row * h) as u64,
            usage:              BufferUsages::COPY_DST | BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        let mut ds = double_slit::DoubleSlitState::new(&ctx.device, &ctx.queue, format, w, h);

        let save_every = (n_frames / 10).max(1);

        for frame in 0..n_frames {
            let mut enc = ctx.device.create_command_encoder(&Default::default());
            ds.update(&ctx.queue);
            ds.compute(&mut enc);
            ds.render(&mut enc, &view);

            if frame % save_every == save_every - 1 {
                enc.copy_texture_to_buffer(
                    tex.as_image_copy(),
                    ImageCopyBuffer {
                        buffer:  &readback,
                        layout:  ImageDataLayout {
                            offset:         0,
                            bytes_per_row:  Some(bytes_per_row),
                            rows_per_image: Some(h),
                        },
                    },
                    Extent3d { width: w, height: h, depth_or_array_layers: 1 },
                );
            }

            ctx.queue.submit([enc.finish()]);

            if frame % save_every == save_every - 1 {
                // Block until GPU is done, then save
                ctx.device.poll(wgpu::Maintain::Wait);

                let (tx, rx) = std::sync::mpsc::channel();
                readback.slice(..).map_async(MapMode::Read, move |r| { tx.send(r).ok(); });
                ctx.device.poll(wgpu::Maintain::Wait);
                rx.recv().unwrap().expect("map_async failed");

                let data = readback.slice(..).get_mapped_range();
                let img: ImageBuffer<Rgba<u8>, Vec<u8>> = ImageBuffer::from_fn(w, h, |x, y| {
                    let row_start = (y * bytes_per_row) as usize;
                    let p = row_start + x as usize * 4;
                    Rgba([data[p], data[p+1], data[p+2], data[p+3]])
                });
                drop(data);
                readback.unmap();

                let path = output_dir.join(format!("frame_{frame:04}.png"));
                img.save(&path).expect("save png");
                log::info!("Saved {}", path.display());
            }
        }

        log::info!("Headless render complete ({n_frames} frames)");
    }

    /// Round up to the nearest multiple of 256 (wgpu requirement for copy row pitch).
    fn align256(n: u32) -> u32 {
        (n + 255) & !255
    }
}
