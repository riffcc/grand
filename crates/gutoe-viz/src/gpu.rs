/*!
 * GUTOE Viz - GPU Context
 * Copyright (C) 2026  Wings
 *
 * AGPL-3.0-or-later
 */

//! wgpu device / queue / surface initialization.

use std::sync::Arc;
use winit::window::Window;
use wgpu::*;

pub struct GpuCtx {
    pub device:  Device,
    pub queue:   Queue,
    pub surface: Surface<'static>,
    config:      SurfaceConfiguration,
}

impl GpuCtx {
    pub async fn new(window: Arc<Window>) -> Self {
        let instance = Instance::new(InstanceDescriptor {
            backends: Backends::all(),
            ..Default::default()
        });

        let surface = instance.create_surface(Arc::clone(&window)).unwrap();

        let adapter = instance
            .request_adapter(&RequestAdapterOptions {
                power_preference:       PowerPreference::HighPerformance,
                compatible_surface:     Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .expect("no suitable wgpu adapter");

        log::info!("Adapter: {:?}", adapter.get_info());

        let (device, queue) = adapter
            .request_device(
                &DeviceDescriptor {
                    label:             Some("gutoe-device"),
                    required_features: Features::empty(),
                    required_limits:   Limits::default(),
                    memory_hints:      MemoryHints::default(),
                },
                None,
            )
            .await
            .expect("device creation failed");

        let size = window.inner_size();
        let caps = surface.get_capabilities(&adapter);

        // Prefer sRGB, fall back to whatever the adapter gives us.
        let format = caps
            .formats
            .iter()
            .copied()
            .find(|f| f.is_srgb())
            .unwrap_or(caps.formats[0]);

        let config = SurfaceConfiguration {
            usage:                        TextureUsages::RENDER_ATTACHMENT,
            format,
            width:                        size.width.max(1),
            height:                       size.height.max(1),
            present_mode:                 PresentMode::AutoVsync,
            alpha_mode:                   caps.alpha_modes[0],
            view_formats:                 vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        Self { device, queue, surface, config }
    }

    pub fn format(&self) -> TextureFormat { self.config.format }
    pub fn size(&self) -> (u32, u32) { (self.config.width, self.config.height) }

    pub fn resize(&mut self, width: u32, height: u32) {
        if width == 0 || height == 0 { return; }
        self.config.width  = width;
        self.config.height = height;
        self.surface.configure(&self.device, &self.config);
    }

    pub fn reconfigure(&mut self) {
        self.surface.configure(&self.device, &self.config);
    }
}

/// Headless wgpu context: no window, no surface.
pub struct HeadlessCtx {
    pub device: Device,
    pub queue:  Queue,
}

impl HeadlessCtx {
    pub async fn new() -> Self {
        let instance = Instance::new(InstanceDescriptor {
            backends: Backends::all(),
            ..Default::default()
        });

        let adapter = instance
            .request_adapter(&RequestAdapterOptions {
                power_preference:       PowerPreference::HighPerformance,
                compatible_surface:     None,
                force_fallback_adapter: false,
            })
            .await
            .expect("no headless adapter");

        let (device, queue) = adapter
            .request_device(
                &DeviceDescriptor {
                    label:             Some("gutoe-headless"),
                    required_features: Features::empty(),
                    required_limits:   Limits::default(),
                    memory_hints:      MemoryHints::default(),
                },
                None,
            )
            .await
            .expect("headless device creation failed");

        Self { device, queue }
    }
}
