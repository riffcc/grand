//! GUTOE Live Visualizer — winit 0.30 + wgpu 22, full 3D perspective
//!
//! Controls:
//!   Left drag   — orbit camera
//!   Scroll      — zoom
//!   Space       — pause / resume
//!   Up / Down   — faster / slower simulation
//!   Escape / Q  — quit

mod gauge;
mod renderer;
mod sim;

use std::collections::HashSet;
use std::sync::Arc;

use rand::rngs::StdRng;
use rand::SeedableRng;

use winit::application::ApplicationHandler;
use winit::dpi::PhysicalPosition;
use winit::event::{ElementState, KeyEvent, MouseButton, MouseScrollDelta, WindowEvent};
use winit::event_loop::{ActiveEventLoop, EventLoop};
use winit::keyboard::{Key, NamedKey};
use winit::window::{Window, WindowId};

use gauge::{update_gauge, GaugeFields, NbrCache};
use renderer::{
    hex_world_pos, look_at, mat4_mul, perspective, state_color, CameraUniforms, InstanceData,
    Renderer, COL_STEP, LAYER_SPACING, ROW_STEP,
};
use sim::{
    detect_quarks, find_proton_triplets, inject_leptons, proton_shell, site_coords, GutoeState,
    LatticeConfig,
};

// ── Phase constants ────────────────────────────────────────────────────────────

const PHASE1_STEPS: u64 = 150;
const N_INJECT: usize = 20;
const GAUGE_EVERY: u64 = 5;

// ── Camera ─────────────────────────────────────────────────────────────────────

struct Camera {
    azimuth: f32,   // horizontal orbit angle (radians)
    elevation: f32, // vertical orbit angle (radians)
    distance: f32,  // distance from target
    target: [f32; 3],
}

impl Camera {
    fn default_for(cfg: &LatticeConfig) -> Self {
        let cx = (cfg.hex_cols as f32 - 1.0) * COL_STEP * 0.5;
        let cy = (cfg.hex_rows as f32 - 1.0) * ROW_STEP * 0.5;
        let cz = (cfg.layers as f32 - 1.0) * LAYER_SPACING * 0.5;
        Self {
            azimuth: 0.6,
            elevation: 0.65,
            distance: 280.0,
            target: [cx, cy, cz],
        }
    }

    fn eye(&self) -> [f32; 3] {
        let az = self.azimuth;
        let el = self.elevation;
        let d = self.distance;
        [
            self.target[0] + d * el.cos() * az.sin(),
            self.target[1] + d * el.sin(),
            self.target[2] - d * el.cos() * az.cos(),
        ]
    }

    fn view_proj(&self, width: u32, height: u32) -> [[f32; 4]; 4] {
        let aspect = width as f32 / height.max(1) as f32;
        let view = look_at(self.eye(), self.target, [0.0, 1.0, 0.0]);
        let proj = perspective(45_f32.to_radians(), aspect, 0.5, 3000.0);
        mat4_mul(&proj, &view)
    }
}

// ── App ────────────────────────────────────────────────────────────────────────

struct App {
    window: Option<Arc<Window>>,
    renderer: Option<Renderer>,
    cfg: LatticeConfig,
    state: GutoeState,
    gauge: GaugeFields,
    nbr_cache: NbrCache,
    rng: StdRng,
    proton_sites: HashSet<usize>,
    proton_shell_sites: HashSet<usize>,
    world_positions: Vec<[f32; 3]>,
    camera: Camera,
    mouse_pressed: bool,
    last_mouse: PhysicalPosition<f64>,
    paused: bool,
    steps_per_frame: u32,
    n_protons: usize,
    n_leptons: usize,
    n_hydrogen: usize,
    enrichment: f32,
    injected: bool,
}

impl App {
    fn new() -> Self {
        let cfg = LatticeConfig {
            void_votes: 4,
            ..Default::default()
        };
        let state = GutoeState::new(&cfg);
        let gauge = GaugeFields::new(&cfg);
        let nbr_cache = NbrCache::build(&cfg);

        let n = cfg.hex_rows * cfg.hex_cols * cfg.layers;
        let world_positions: Vec<[f32; 3]> = (0..n)
            .map(|site| {
                let (r, c, z) = site_coords(site, &cfg);
                hex_world_pos(r, c, z, cfg.hex_rows)
            })
            .collect();

        App {
            window: None,
            renderer: None,
            camera: Camera::default_for(&cfg),
            cfg,
            state,
            gauge,
            nbr_cache,
            rng: StdRng::seed_from_u64(42),
            proton_sites: HashSet::new(),
            proton_shell_sites: HashSet::new(),
            world_positions,
            mouse_pressed: false,
            last_mouse: PhysicalPosition { x: 0.0, y: 0.0 },
            paused: false,
            steps_per_frame: 1,
            n_protons: 0,
            n_leptons: 0,
            n_hydrogen: 0,
            enrichment: 0.0,
            injected: false,
        }
    }

    fn simulation_step(&mut self) {
        if self.state.step_count == PHASE1_STEPS && !self.injected {
            let quarks = detect_quarks(&self.state, &self.cfg);
            let triplets = find_proton_triplets(&quarks, &self.cfg);
            inject_leptons(
                &mut self.state,
                &triplets,
                N_INJECT,
                &mut self.rng,
                &self.cfg,
            );
            self.state.phase = 2;
            self.injected = true;
            log::info!(
                "Phase 2: injected {} leptons, {} proton triplets",
                N_INJECT,
                triplets.len()
            );
        }

        if self.state.phase == 2 && self.state.step_count % GAUGE_EVERY == 0 {
            let quarks = detect_quarks(&self.state, &self.cfg);
            let triplets = find_proton_triplets(&quarks, &self.cfg);
            self.proton_sites = triplets
                .iter()
                .flat_map(|&(d, u1, u2)| [d, u1, u2])
                .collect();
            self.proton_shell_sites = proton_shell(&triplets, &self.cfg);

            let q_map: std::collections::HashMap<usize, sim::QuarkType> = quarks
                .iter()
                .map(|q| (q.site, q.quark_type.clone()))
                .collect();

            update_gauge(
                &mut self.gauge,
                &self.state.lattice,
                &q_map,
                &self.proton_sites,
                &self.cfg,
                &self.nbr_cache,
            );

            self.n_protons = triplets.len();
            self.n_leptons = self
                .state
                .lattice
                .iter()
                .filter(|&&s| s == sim::LEPTON_SEED)
                .count();

            let shell_lep = self
                .proton_shell_sites
                .iter()
                .filter(|&&s| self.state.lattice[s] == sim::LEPTON_SEED)
                .count();
            let layer_stride = self.cfg.hex_rows * self.cfg.hex_cols;
            let proton_layers: HashSet<usize> = self
                .proton_sites
                .iter()
                .map(|&s| s / layer_stride)
                .collect();
            let bg_sites: Vec<usize> = (0..self.state.n(&self.cfg))
                .filter(|s| {
                    proton_layers.contains(&(s / layer_stride))
                        && !self.proton_sites.contains(s)
                        && !self.proton_shell_sites.contains(s)
                })
                .collect();
            let bg_lep = bg_sites
                .iter()
                .filter(|&&s| self.state.lattice[s] == sim::LEPTON_SEED)
                .count();
            let rs = shell_lep as f32 / self.proton_shell_sites.len().max(1) as f32;
            let rb = bg_lep as f32 / bg_sites.len().max(1) as f32;
            self.enrichment = if rb > 1e-6 {
                (rs / rb).min(20.0)
            } else if rs > 0.0 {
                20.0
            } else {
                0.0
            };
            self.n_hydrogen = shell_lep;
        }

        self.state.step(
            &mut self.rng,
            &self.cfg,
            &self.gauge.phi,
            &self.proton_sites,
        );
    }

    fn build_instances(&self) -> Vec<InstanceData> {
        let n = self.state.n(&self.cfg);
        let phi_max = self.gauge.phi.iter().cloned().fold(0.01_f64, f64::max);

        (0..n)
            .map(|site| {
                let s = self.state.lattice[site];
                let wpos = self.world_positions[site];
                let phi_n = (self.gauge.phi[site] / phi_max).clamp(0.0, 1.0) as f32;
                let mut col = state_color(s);
                col[1] = (col[1] + phi_n * 0.4).min(1.0); // φ → green tint
                if self.proton_sites.contains(&site) {
                    col = [1.0, 0.92, 0.0, 1.0]; // proton quark → yellow
                }
                if s == sim::LEPTON_SEED && self.proton_shell_sites.contains(&site) {
                    col = [0.2, 1.0, 0.5, 1.0]; // H-atom lepton → green
                }
                InstanceData {
                    world_pos: wpos,
                    _pad: 0.0,
                    color: col,
                }
            })
            .collect()
    }

    fn window_title(&self) -> String {
        format!(
            "GUTOE | Phase {} | t={} | p={} | γ⁰={} | H={} | enrich={:.2}×  \
             [drag=orbit  scroll=zoom  space=pause  ↑↓=speed]",
            self.state.phase,
            self.state.step_count,
            self.n_protons,
            self.n_leptons,
            self.n_hydrogen,
            self.enrichment,
        )
    }
}

// ── ApplicationHandler ────────────────────────────────────────────────────────

impl ApplicationHandler for App {
    fn resumed(&mut self, el: &ActiveEventLoop) {
        let attrs = Window::default_attributes()
            .with_title("GUTOE Visualizer")
            .with_inner_size(winit::dpi::LogicalSize::new(1200u32, 900u32));
        let window = Arc::new(el.create_window(attrs).expect("window creation failed"));

        let wgpu_instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });
        let surface = wgpu_instance
            .create_surface(Arc::clone(&window))
            .expect("surface creation failed");

        let size = window.inner_size();
        let n = self.cfg.hex_rows * self.cfg.hex_cols * self.cfg.layers;

        let renderer = pollster::block_on(Renderer::new(
            &wgpu_instance,
            surface,
            size.width,
            size.height,
            n,
        ));

        let instances = self.build_instances();
        renderer.update_instances(&instances);

        let cam_uni = CameraUniforms {
            view_proj: self.camera.view_proj(size.width, size.height),
        };
        renderer.update_camera(&cam_uni);

        self.renderer = Some(renderer);
        window.request_redraw();
        self.window = Some(window);
    }

    fn window_event(&mut self, el: &ActiveEventLoop, _id: WindowId, event: WindowEvent) {
        match event {
            WindowEvent::CloseRequested => el.exit(),

            WindowEvent::KeyboardInput {
                event:
                    KeyEvent {
                        logical_key,
                        state: ElementState::Pressed,
                        ..
                    },
                ..
            } => match logical_key {
                Key::Named(NamedKey::Escape) => el.exit(),
                Key::Named(NamedKey::Space) => {
                    self.paused = !self.paused;
                }
                Key::Named(NamedKey::ArrowUp) => {
                    self.steps_per_frame = (self.steps_per_frame * 2).min(64);
                }
                Key::Named(NamedKey::ArrowDown) => {
                    self.steps_per_frame = (self.steps_per_frame / 2).max(1);
                }
                Key::Character(ref s) if s == "q" || s == "Q" => el.exit(),
                _ => {}
            },

            // ── Mouse controls ─────────────────────────────────────────────────
            WindowEvent::MouseInput {
                button: MouseButton::Left,
                state,
                ..
            } => {
                self.mouse_pressed = state == ElementState::Pressed;
            }

            WindowEvent::CursorMoved { position, .. } => {
                if self.mouse_pressed {
                    let dx = (position.x - self.last_mouse.x) as f32;
                    let dy = (position.y - self.last_mouse.y) as f32;
                    self.camera.azimuth += dx * 0.005;
                    self.camera.elevation = (self.camera.elevation - dy * 0.004).clamp(0.05, 1.50);
                }
                self.last_mouse = position;
            }

            WindowEvent::MouseWheel { delta, .. } => {
                let scroll = match delta {
                    MouseScrollDelta::LineDelta(_, y) => y,
                    MouseScrollDelta::PixelDelta(p) => p.y as f32 * 0.01,
                };
                self.camera.distance =
                    (self.camera.distance * 0.9_f32.powf(scroll)).clamp(40.0, 1200.0);
            }

            WindowEvent::Resized(size) => {
                if let Some(r) = &mut self.renderer {
                    r.resize(size.width, size.height);
                }
            }

            WindowEvent::RedrawRequested => {
                if !self.paused {
                    for _ in 0..self.steps_per_frame {
                        self.simulation_step();
                    }
                }

                if let Some(r) = &self.renderer {
                    let instances = self.build_instances();
                    r.update_instances(&instances);

                    let cam_uni = CameraUniforms {
                        view_proj: self.camera.view_proj(r.config.width, r.config.height),
                    };
                    r.update_camera(&cam_uni);

                    match r.render() {
                        Ok(()) => {}
                        Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                            if let Some(w) = &self.window {
                                if let Some(r) = &mut self.renderer {
                                    let s = w.inner_size();
                                    r.resize(s.width, s.height);
                                }
                            }
                        }
                        Err(wgpu::SurfaceError::OutOfMemory) => el.exit(),
                        Err(e) => log::error!("render error: {e:?}"),
                    }
                }

                if let Some(w) = &self.window {
                    w.set_title(&self.window_title());
                    w.request_redraw();
                }
            }

            _ => {}
        }
    }
}

// ── Entry point ────────────────────────────────────────────────────────────────

fn main() {
    env_logger::init();
    let event_loop = EventLoop::new().expect("event loop failed");
    let mut app = App::new();
    event_loop.run_app(&mut app).expect("run failed");
}
