//! GUTOE Live Visualizer — winit 0.30 + wgpu 22
//!
//! Controls:
//!   Space     — pause / resume
//!   Up        — faster (more steps per frame)
//!   Down      — slower
//!   Escape/Q  — quit

mod sim;
mod gauge;
mod renderer;

use std::collections::HashSet;
use std::sync::Arc;

use rand::rngs::StdRng;
use rand::SeedableRng;

use winit::application::ApplicationHandler;
use winit::event::{ElementState, KeyEvent, WindowEvent};
use winit::event_loop::{ActiveEventLoop, EventLoop};
use winit::keyboard::{Key, NamedKey};
use winit::window::{Window, WindowId};

use sim::{
    LatticeConfig, GutoeState,
    detect_quarks, find_proton_triplets, inject_leptons, proton_shell, site_coords,
};
use gauge::{GaugeFields, NbrCache, update_gauge};
use renderer::{Renderer, InstanceData, hex_world_pos, state_color};

// ── Phase constants ────────────────────────────────────────────────────────────

const PHASE1_STEPS: u64 = 150;
const N_INJECT: usize   = 20;
const GAUGE_EVERY: u64  = 5;

// ── App state ──────────────────────────────────────────────────────────────────

struct App {
    // winit / wgpu
    window:          Option<Arc<Window>>,
    renderer:        Option<Renderer>,
    // simulation
    cfg:             LatticeConfig,
    state:           GutoeState,
    gauge:           GaugeFields,
    nbr_cache:       NbrCache,
    rng:             StdRng,
    proton_sites:    HashSet<usize>,
    proton_shell_sites: HashSet<usize>,
    // precomputed world positions (one per site)
    world_positions: Vec<[f32; 2]>,
    // control
    paused:          bool,
    steps_per_frame: u32,
    // stats
    n_protons:       usize,
    n_leptons:       usize,
    n_hydrogen:      usize,
    enrichment:      f32,
    injected:        bool,
}

impl App {
    fn new() -> Self {
        let cfg = LatticeConfig { void_votes: 4, ..Default::default() };
        let state  = GutoeState::new(&cfg);
        let gauge  = GaugeFields::new(&cfg);
        let nbr_cache = NbrCache::build(&cfg);

        // Pre-compute world positions for all sites
        let n = cfg.hex_rows * cfg.hex_cols * cfg.layers;
        let mut world_positions = Vec::with_capacity(n);
        for site in 0..n {
            let (r, c, z) = site_coords(site, &cfg);
            world_positions.push(hex_world_pos(r, c, z));
        }

        App {
            window: None,
            renderer: None,
            cfg,
            state,
            gauge,
            nbr_cache,
            rng: StdRng::seed_from_u64(42),
            proton_sites: HashSet::new(),
            proton_shell_sites: HashSet::new(),
            world_positions,
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
        // Phase 1 → 2 transition
        if self.state.step_count == PHASE1_STEPS && !self.injected {
            let quarks   = detect_quarks(&self.state, &self.cfg);
            let triplets = find_proton_triplets(&quarks, &self.cfg);
            inject_leptons(&mut self.state, &triplets, N_INJECT, &mut self.rng, &self.cfg);
            self.state.phase = 2;
            self.injected = true;
            log::info!("Phase 2: injected {} leptons. Proton triplets: {}", N_INJECT, triplets.len());
        }

        // Gauge update every GAUGE_EVERY steps in phase 2
        if self.state.phase == 2 && self.state.step_count % GAUGE_EVERY == 0 {
            let quarks   = detect_quarks(&self.state, &self.cfg);
            let triplets = find_proton_triplets(&quarks, &self.cfg);
            self.proton_sites = triplets.iter()
                .flat_map(|&(d, u1, u2)| [d, u1, u2])
                .collect();
            self.proton_shell_sites = proton_shell(&triplets, &self.cfg);

            let q_map: std::collections::HashMap<usize, sim::QuarkType> = quarks.iter()
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

            // Update stats
            self.n_protons = triplets.len();
            self.n_leptons = self.state.lattice.iter().filter(|&&s| s == sim::LEPTON_SEED).count();
            let shell_lep = self.proton_shell_sites.iter()
                .filter(|&&s| self.state.lattice[s] == sim::LEPTON_SEED)
                .count();
            let layer_stride = self.cfg.hex_rows * self.cfg.hex_cols;
            let proton_layers: HashSet<usize> = self.proton_sites.iter()
                .map(|&s| s / layer_stride)
                .collect();
            let bg_sites: Vec<usize> = (0..self.state.n(&self.cfg))
                .filter(|s| {
                    (s / layer_stride) < self.cfg.layers
                        && proton_layers.contains(&(s / layer_stride))
                        && !self.proton_sites.contains(s)
                        && !self.proton_shell_sites.contains(s)
                })
                .collect();
            let bg_lep = bg_sites.iter()
                .filter(|&&s| self.state.lattice[s] == sim::LEPTON_SEED)
                .count();
            let shell_sz = self.proton_shell_sites.len().max(1);
            let bg_sz    = bg_sites.len().max(1);
            let rs = shell_lep as f32 / shell_sz as f32;
            let rb = bg_lep   as f32 / bg_sz   as f32;
            self.enrichment = if rb > 1e-6 { (rs / rb).min(20.0) } else if rs > 0.0 { 20.0 } else { 0.0 };

            // H atom count: leptons in proton shell
            self.n_hydrogen = shell_lep;
        }

        self.state.step(&mut self.rng, &self.cfg, &self.gauge.phi, &self.proton_sites);
    }

    fn build_instances(&self) -> Vec<InstanceData> {
        let n = self.state.n(&self.cfg);
        let phi_max = self.gauge.phi.iter().cloned().fold(0.01_f64, f64::max);
        let mut instances = Vec::with_capacity(n);

        for site in 0..n {
            let s     = self.state.lattice[site];
            let wpos  = self.world_positions[site];
            let phi_n = (self.gauge.phi[site] / phi_max).clamp(0.0, 1.0) as f32;

            let mut color = state_color(s);

            // φ overlay: green tint
            color[1] = (color[1] + phi_n * 0.4).min(1.0);

            // Proton quark override → bright yellow
            if self.proton_sites.contains(&site) {
                color = [1.0, 0.92, 0.0, 1.0];
            }

            // Lepton in shell (H atom) → bright green
            if s == sim::LEPTON_SEED && self.proton_shell_sites.contains(&site) {
                color = [0.2, 1.0, 0.5, 1.0];
            }

            instances.push(InstanceData { world_pos: wpos, color });
        }
        instances
    }

    fn window_title(&self) -> String {
        format!(
            "GUTOE | Phase {} | t={} | p={} | γ⁰={} | H={} | enrich={:.2}×",
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
        let win_attrs = Window::default_attributes()
            .with_title("GUTOE Visualizer")
            .with_inner_size(winit::dpi::LogicalSize::new(1200u32, 900u32));
        let window = Arc::new(el.create_window(win_attrs).expect("Failed to create window"));

        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        let surface = instance.create_surface(Arc::clone(&window))
            .expect("Failed to create surface");

        let size = window.inner_size();
        let scale_factor = window.scale_factor();
        let n = self.cfg.hex_rows * self.cfg.hex_cols * self.cfg.layers;

        let renderer = pollster::block_on(Renderer::new(
            &instance,
            surface,
            size.width,
            size.height,
            scale_factor,
            n,
        ));

        // Upload initial (all-VOID) instance data
        let instances = self.build_instances();
        renderer.update_instances(&instances);

        self.renderer = Some(renderer);
        self.window = Some(window);
    }

    fn window_event(&mut self, el: &ActiveEventLoop, _id: WindowId, event: WindowEvent) {
        match event {
            WindowEvent::CloseRequested => {
                el.exit();
            }

            WindowEvent::KeyboardInput {
                event: KeyEvent { logical_key, state: ElementState::Pressed, .. },
                ..
            } => match logical_key {
                Key::Named(NamedKey::Escape) => el.exit(),
                Key::Named(NamedKey::Space) => {
                    self.paused = !self.paused;
                    if let Some(w) = &self.window {
                        let base = self.window_title();
                        w.set_title(&if self.paused { format!("{base} [PAUSED]") } else { base });
                    }
                }
                Key::Named(NamedKey::ArrowUp) => {
                    self.steps_per_frame = (self.steps_per_frame * 2).min(64);
                    log::info!("steps_per_frame → {}", self.steps_per_frame);
                }
                Key::Named(NamedKey::ArrowDown) => {
                    self.steps_per_frame = (self.steps_per_frame / 2).max(1);
                    log::info!("steps_per_frame → {}", self.steps_per_frame);
                }
                Key::Character(ref s) if s == "q" || s == "Q" => el.exit(),
                _ => {}
            },

            WindowEvent::Resized(size) => {
                let sf = self.window.as_ref().map_or(1.0, |w| w.scale_factor());
                if let Some(r) = &mut self.renderer {
                    r.resize(size.width, size.height, sf);
                }
            }

            WindowEvent::ScaleFactorChanged { scale_factor, inner_size_writer: _ } => {
                let size = self.window.as_ref().map(|w| w.inner_size()).unwrap_or_default();
                if let Some(r) = &mut self.renderer {
                    r.resize(size.width, size.height, scale_factor);
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
                    match r.render() {
                        Ok(()) => {}
                        Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                            if let Some(w) = &self.window {
                                let size = w.inner_size();
                                let sf   = w.scale_factor();
                                if let Some(r) = &mut self.renderer {
                                    r.resize(size.width, size.height, sf);
                                }
                            }
                        }
                        Err(wgpu::SurfaceError::OutOfMemory) => el.exit(),
                        Err(e) => log::error!("Surface error: {e:?}"),
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
    let event_loop = EventLoop::new().expect("Failed to create event loop");
    let mut app = App::new();
    event_loop.run_app(&mut app).expect("Event loop failed");
}
