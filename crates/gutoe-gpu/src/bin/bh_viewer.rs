//! GUTOE Black Hole Viewer — real-time interactive geodesic ray tracer
//!
//! Runs the GUTOE-corrected Schwarzschild metric null-geodesic integrator
//! entirely on the GPU in a WGSL fragment shader.
//!
//! Controls:
//!   Left drag (vertical)   — inclination: 0° face-on ↔ 90° edge-on
//!   Left drag (horizontal) — disk plane rotation (azimuth)
//!   Scroll                 — zoom (field of view in units of r_s)
//!   +/-                    — disk outer radius
//!   G                      — toggle GUTOE lattice core (r_c on/off)
//!   R                      — reset camera
//!   Q / Escape             — quit

use std::sync::Arc;

use winit::{
    application::ApplicationHandler,
    dpi::PhysicalPosition,
    event::{ElementState, KeyEvent, MouseButton, MouseScrollDelta, WindowEvent},
    event_loop::{ActiveEventLoop, EventLoop},
    keyboard::{Key, NamedKey},
    window::{Window, WindowId},
};

// ── WGSL shader ───────────────────────────────────────────────────────────────

const SHADER: &str = r#"
// GUTOE Black Hole Shader — null geodesic tracer with stars, Doppler, gravitational redshift
//
// Physics:
//   d²r/dφ² = r(2r²+r_c²)/b² − r + r_s·r(r²+2r_c²)/(2r_eff³)   [GUTOE Schwarzschild]
//   Disk: Novikov-Thorne T ∝ (r_ISCO/r)^(3/4)
//   Doppler: D³ factor for approaching / receding disk gas (Keplerian orbit)
//   Gravitational redshift: sqrt(1 - r_s/r_eff)
//   Stars: procedural hash field; gravitationally lensed by tracing exit angle

struct Params {
    r_s      : f32,   // Schwarzschild radius (internal units = 1)
    r_c      : f32,   // GUTOE lattice core radius = sqrt(C_inf) × l_P
    disk_in  : f32,   // disk inner areal radius (× r_s)
    disk_out : f32,   // disk outer areal radius (× r_s)
    sin_inc  : f32,   // sin(observer inclination): 1=edge-on, 0=face-on
    fov      : f32,   // half-width of image in r_s
    width    : f32,   // viewport width  (pixels)
    height   : f32,   // viewport height (pixels)
    max_phi  : f32,   // max integration angle (radians)
    dphi     : f32,   // RK4 step size (radians)
    az       : f32,   // azimuth offset (rotates disk orientation on screen)
    _pad     : f32,
}
@group(0) @binding(0) var<uniform> P : Params;

// ── Hash / Noise ──────────────────────────────────────────────────────────────

fn hash21(p: vec2<f32>) -> f32 {
    var q = fract(p * vec2(127.1, 311.7));
    q += dot(q, q.yx + 19.19);
    return fract(q.x * q.y);
}

// ── Star field ────────────────────────────────────────────────────────────────
// Sample a procedural star field at sky direction 'sky' (in r_s units).
// ~8% of cells contain a star; magnitude and color vary per cell.

fn starfield(sky: vec2<f32>) -> vec3<f32> {
    let cell_size = P.r_s * 0.50;          // angular cell size in r_s units
    let cell      = floor(sky / cell_size);
    let local     = fract(sky / cell_size) - 0.5;  // local coords in [-0.5, 0.5]²

    let h = hash21(cell);
    if h > 0.92 { return vec3(0.0); }      // only ~8% of cells have a star

    let magnitude = 1.0 - h / 0.92;        // [0,1]: h≈0 → bright, h≈0.92 → faint

    // Star position within cell (sub-cell random offset)
    let jitter    = vec2(hash21(cell + 0.31), hash21(cell + 0.73)) - 0.5;
    let d2        = dot(local - jitter * 0.7, local - jitter * 0.7);

    // PSF: sharp Airy-like core + faint diffraction halo
    let sigma     = 0.003 * (0.15 + magnitude);
    let core      = exp(-d2 / sigma);
    let halo      = exp(-d2 / (sigma * 8.0)) * 0.12;

    // Color temperature: blue-white (hot) ↔ orange-red (cool)
    let temp      = hash21(cell + 0.99);
    let star_col  = mix(vec3(1.0, 0.50, 0.20), vec3(0.70, 0.85, 1.0), temp);

    return star_col * (core + halo) * magnitude * 2.5;
}

// ── Orbit integrator ──────────────────────────────────────────────────────────

fn accel(r: f32, b: f32, r_s: f32, r_c: f32) -> f32 {
    let re2 = r * r + r_c * r_c;
    let re3 = re2 * sqrt(re2);
    return r * (2.0 * r * r + r_c * r_c) / (b * b)
         - r
         + r_s * r * (r * r + 2.0 * r_c * r_c) / (2.0 * re3);
}

fn orbit_vr_sq(r: f32, b: f32, r_s: f32, r_c: f32) -> f32 {
    let re2 = r * r + r_c * r_c;
    let re  = sqrt(re2);
    let f   = 1.0 - r_s / re;
    return r * r * re2 / (b * b) - r * r * f;
}

fn rk4(r: f32, p: f32, b: f32, r_s: f32, r_c: f32, h: f32) -> vec2<f32> {
    let k1r = p;                               let k1p = accel(r,                     b, r_s, r_c);
    let k2r = p + 0.5 * h * k1p;              let k2p = accel(r + 0.5 * h * k1r,    b, r_s, r_c);
    let k3r = p + 0.5 * h * k2p;              let k3p = accel(r + 0.5 * h * k2r,    b, r_s, r_c);
    let k4r = p + h * k3p;                    let k4p = accel(r + h * k3r,           b, r_s, r_c);
    return vec2(
        r + h * (k1r + 2.0*k2r + 2.0*k3r + k4r) / 6.0,
        p + h * (k1p + 2.0*k2p + 2.0*k3p + k4p) / 6.0,
    );
}

// ── Tracer ────────────────────────────────────────────────────────────────────
// Returns vec4(r_eff_hit, phi_total, f32(n_cross), kind)
//   kind 0 = captured, 1 = disk hit, 2 = escaped
// phi_total is the total orbital angle (used for gravitational lensing of stars).

fn trace(bx_in: f32, by_in: f32) -> vec4<f32> {
    // Rotate impact parameters by azimuth (spins the disk orientation on screen)
    let ca = cos(P.az); let sa = sin(P.az);
    let bx = ca * bx_in - sa * by_in;
    let by = sa * bx_in + ca * by_in;

    let b = sqrt(bx * bx + by * by);
    if b < 0.0001 { return vec4(0.0, 0.0, 0.0, 0.0); }

    let r_s   = P.r_s;
    let r_c   = P.r_c;

    // Deep-shadow shortcut: b < b_crit/2 → definitely captured.
    // b_crit = (3√3/2) r_s ≈ 2.598 r_s; half is (3√3/4) r_s ≈ 1.299 r_s.
    let b_crit_half = 0.75 * sqrt(3.0) * r_s;
    if b < b_crit_half { return vec4(0.0, 0.0, 0.0, 0.0); }

    let sin_i = by / b;
    let is_eq = abs(sin_i) < 1e-5;

    // r_start = 3b: orbit_accel ≈ 2r³/b², p ≈ r²/b → change per step ≈ 6·dphi = 6%.
    // Scale-independent and stable for all b > 0.5·b_crit (deep shadow already handled).
    let r_start = 3.0 * b;

    // Initial radial velocity: ingoing from r_start
    let re0   = sqrt(r_start * r_start + r_c * r_c);
    let v0sq  = r_start * r_start * re0 * re0 / (b * b)
              - r_start * r_start * (1.0 - r_s / re0);
    var p_init = -r_start * r_start / b;
    if v0sq > 0.0 { p_init = -sqrt(v0sq); }

    // Capture threshold: areal radius r_eff = sqrt(r²+r_c²) drops below the horizon.
    // The event horizon is at areal radius r_s (not the coordinate radius sqrt(r_s²-r_c²)).
    let r_cap = r_s * 0.99;

    var r          = r_start;
    var p          = p_init;
    var phi        = 0.0f;
    var n_cross    = 0u;
    var turned     = false;
    var in_disk_eq = false;

    let max_steps = i32(P.max_phi / P.dphi) + 1;

    for (var i = 0; i < max_steps; i++) {
        let s      = rk4(r, p, b, r_s, r_c, P.dphi);
        let rn     = s.x;
        let pn_rk4 = s.y;
        // Enforce orbital constraint p² = orbit_vr_sq(r). Prevents centrifugal blowup;
        // preserves direction sign from RK4 → correctly triggers turned / capture.
        let vr2n   = max(orbit_vr_sq(rn, b, r_s, r_c), 0.0);
        let pn     = select(-sqrt(vr2n), sqrt(vr2n), pn_rk4 >= 0.0);
        let phin = phi + P.dphi;
        let ren  = sqrt(rn * rn + r_c * r_c);

        // Capture: inside horizon, or coordinate r went below core / negative
        if rn < r_c * 0.01 || ren < r_cap { return vec4(0.0, 0.0, 0.0, 0.0); }

        // Turning point (ingoing → outgoing)
        if !turned && p < 0.0 && pn >= 0.0 { turned = true; }

        // Escape: returned to launch radius
        if turned && rn >= r_start * 0.99 { return vec4(0.0, phin, 0.0, 2.0); }

        // Disk detection
        if is_eq {
            let re_cur = sqrt(r * r + r_c * r_c);
            let in_d   = re_cur >= P.disk_in && re_cur <= P.disk_out;
            if !in_disk_eq && in_d && p < 0.0 {
                return vec4(re_cur, phi, 1.0, 1.0);
            }
            in_disk_eq = in_d;
        } else {
            let phi_cross = f32(n_cross + 1u) * 3.14159265359;
            if phi < phi_cross && phin >= phi_cross {
                let t    = (phi_cross - phi) / P.dphi;
                let r_x  = r + t * (rn - r);
                let re_x = sqrt(r_x * r_x + r_c * r_c);
                n_cross += 1u;
                if re_x >= P.disk_in && re_x <= P.disk_out {
                    return vec4(re_x, phi_cross, f32(n_cross), 1.0);
                }
            }
        }

        r   = rn;
        p   = pn;
        phi = phin;
    }

    if r >= r_start * 0.5 { return vec4(0.0, phi, 0.0, 2.0); }
    return vec4(0.0, 0.0, 0.0, 0.0);
}

// ── Disk colour: Novikov-Thorne + Doppler + gravitational redshift ─────────────

fn disk_color(r_eff: f32, n_cross: u32, sx: f32) -> vec4<f32> {
    let r_isco = 3.0 * P.r_s;

    // Novikov-Thorne temperature profile: T ∝ (r_ISCO/r)^(3/4)
    let t_nt  = pow(clamp(r_isco / r_eff, 0.01, 1.0), 0.75);

    // Gravitational redshift: photons climbing out of the gravitational well
    let grav  = sqrt(max(1.0 - P.r_s / r_eff, 0.01));

    // Keplerian orbital velocity at r_eff (Schwarzschild, units c = 1)
    let v_k   = sqrt(max(0.5 * P.r_s / r_eff, 0.0));

    // Line-of-sight component: sx ≈ r_eff × cos(azimuthal angle)
    // Approaching side (sx < 0 for prograde): blueshifted → brighter
    let cos_az = clamp(sx / r_eff, -1.0, 1.0);
    let beta   = v_k * (-cos_az) * P.sin_inc;   // + = approaching (blueshift)

    // Relativistic Doppler + beaming: I ∝ D³ where D = sqrt((1+β)/(1-β))
    let D3     = pow(clamp((1.0 + beta) / max(1.0 - beta, 0.01), 0.1, 10.0), 3.0);

    // Higher-order images (photon ring, secondary ring...) are progressively fainter
    let fade   = pow(0.55, f32(n_cross) - 1.0);

    // Combined specific intensity (not clamped — allow HDR to show bright crescent)
    let bright = t_nt * grav * D3 * fade;

    // Spectral mapping: dim = deep amber, bright = yellow, very bright = blue-white
    let b1 = clamp(bright, 0.0, 1.0);
    let bx = clamp(bright - 1.0, 0.0, 2.0);   // HDR overflow → white/blue glow
    let r_ch  = clamp(pow(b1, 0.28) + bx * 0.55, 0.0, 1.0);
    let g_ch  = clamp(pow(b1, 0.62) * 0.88 + bx * 0.50, 0.0, 1.0);
    let bl_ch = clamp(pow(b1, 2.1)  * 0.42 + bx * 0.45, 0.0, 1.0);
    return vec4(r_ch, g_ch, bl_ch, 1.0);
}

// ── Vertex / Fragment ─────────────────────────────────────────────────────────

struct VOut { @builtin(position) pos: vec4<f32> }

@vertex
fn vs(@builtin(vertex_index) vi: u32) -> VOut {
    var xy = array<vec2<f32>, 6>(
        vec2(-1.0, -1.0), vec2(1.0, -1.0), vec2(-1.0, 1.0),
        vec2(-1.0,  1.0), vec2(1.0, -1.0), vec2(1.0,  1.0),
    );
    return VOut(vec4(xy[vi], 0.0, 1.0));
}

@fragment
fn fs(@builtin(position) frag: vec4<f32>) -> @location(0) vec4<f32> {
    let scale = 2.0 * P.fov * P.r_s / P.width;
    let sx = (frag.x - P.width  * 0.5) * scale;
    let sy = (P.height * 0.5 - frag.y) * scale;

    let bx = sx;
    let by = sy * P.sin_inc;

    let hit  = trace(bx, by);
    let kind = hit.w;

    // Shadow — pure black
    if kind < 0.5 { return vec4(0.0, 0.0, 0.0, 1.0); }

    // Disk hit — Novikov-Thorne + Doppler + gravitational redshift
    // For Doppler, use the x-component in the disk frame (rotated by azimuth),
    // so the bright crescent stays on the correct side when azimuth is changed.
    if kind < 1.5 {
        let sx_disk = cos(P.az) * sx - sin(P.az) * sy;
        return disk_color(hit.x, u32(hit.z), sx_disk);
    }

    // Escaped photon — background with gravitationally lensed star field
    //
    // The photon swept phi_total radians in its orbital plane.
    // A straight photon (no BH) sweeps exactly π. The excess delta = phi_total − π
    // is the deflection angle. We reverse-rotate the sky direction by delta so
    // stars appear at their true (source) positions, giving lensing arcs near
    // the photon sphere automatically.
    let phi_total  = hit.y;
    let deflection = phi_total - 3.14159265359;
    let cd = cos(deflection);
    let sd = sin(deflection);
    let sky = vec2(cd * sx - sd * sy,
                   sd * sx + cd * sy);

    let stars = starfield(sky);

    // Very faint blue-purple nebula haze, brighter near image edges
    let r_dist = length(vec2(sx, sy)) / max(P.r_s * P.fov, 0.01);
    let nebula = vec3(
        0.002 + 0.003 * r_dist * r_dist,
        0.001 + 0.001 * r_dist,
        0.010 + 0.012 * r_dist,
    );

    return vec4(nebula + stars, 1.0);
}
"#;

// ── CPU-side uniform struct (must match WGSL layout, std140 / 16-byte align) ──

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    r_s:      f32,
    r_c:      f32,
    disk_in:  f32,
    disk_out: f32,
    sin_inc:  f32,
    fov:      f32,
    width:    f32,
    height:   f32,
    max_phi:  f32,
    dphi:     f32,
    az:       f32,
    _pad:     f32,
}

// ── Camera state ──────────────────────────────────────────────────────────────

struct Camera {
    inclination: f32, // degrees, 0 = face-on, 90 = edge-on
    azimuth:     f32, // radians, rotates disk orientation on screen
    fov_rs:      f32, // half-width of image in r_s
    disk_outer:  f32, // disk outer radius in r_s
    gutoe_core:  bool, // toggle GUTOE lattice correction
}

impl Default for Camera {
    fn default() -> Self {
        Self {
            // EHT M87 geometry: ~17° from face-on, bright crescent at bottom.
            // azimuth = -π/2 rotates the disk so the approaching (bright) side is
            // at screen-bottom, matching the classic EHT image orientation.
            inclination: 17.0,
            azimuth:     -std::f32::consts::FRAC_PI_2,
            fov_rs:      7.0,
            disk_outer:  20.0,
            gutoe_core:  true,
        }
    }
}

impl Camera {
    fn params(&self, width: f32, height: f32) -> Params {
        // r_s = 1 in internal units; r_core = sqrt(C_inf) * l_P
        const C_INF: f32 = 0.5466;
        let r_c = if self.gutoe_core { C_INF.sqrt() } else { 0.0 };
        Params {
            r_s:      1.0,
            r_c,
            disk_in:  3.0,  // r_ISCO = 3 r_s
            disk_out: self.disk_outer,
            sin_inc:  self.inclination.to_radians().sin(),
            fov:      self.fov_rs,
            width,
            height,
            max_phi:  20.0 * std::f32::consts::PI,
            dphi:     0.02,
            az:       self.azimuth,
            _pad:     0.0,
        }
    }
}

// ── wgpu resources ────────────────────────────────────────────────────────────

struct Gpu {
    surface:    wgpu::Surface<'static>,
    device:     wgpu::Device,
    queue:      wgpu::Queue,
    config:     wgpu::SurfaceConfiguration,
    pipeline:   wgpu::RenderPipeline,
    bind_group: wgpu::BindGroup,
    uni_buf:    wgpu::Buffer,
}

impl Gpu {
    fn new(window: Arc<Window>) -> Self {
        pollster::block_on(Self::new_async(window))
    }

    async fn new_async(window: Arc<Window>) -> Self {
        let size = window.inner_size();

        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::PRIMARY,
            ..Default::default()
        });
        let surface = instance.create_surface(window).unwrap();
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference:       wgpu::PowerPreference::HighPerformance,
                compatible_surface:     Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .expect("no GPU adapter found — is a display connected?");

        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor::default(), None)
            .await
            .unwrap();

        let caps   = surface.get_capabilities(&adapter);
        let format = caps.formats.iter()
            .copied()
            .find(|f| f.is_srgb())
            .unwrap_or(caps.formats[0]);

        let config = wgpu::SurfaceConfiguration {
            usage:        wgpu::TextureUsages::RENDER_ATTACHMENT,
            format,
            width:        size.width.max(1),
            height:       size.height.max(1),
            present_mode: wgpu::PresentMode::AutoNoVsync, // no latency cap — render as fast as GPU allows
            alpha_mode:   caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 1, // minimal latency: submit → present immediately
        };
        surface.configure(&device, &config);

        // Uniform buffer
        let uni_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label:              Some("params"),
            size:               std::mem::size_of::<Params>() as u64,
            usage:              wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // Bind group layout
        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label:   Some("bgl"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding:    0,
                visibility: wgpu::ShaderStages::FRAGMENT,
                ty:         wgpu::BindingType::Buffer {
                    ty:                 wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size:   None,
                },
                count: None,
            }],
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label:  Some("bg"),
            layout: &bgl,
            entries: &[wgpu::BindGroupEntry {
                binding:  0,
                resource: uni_buf.as_entire_binding(),
            }],
        });

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label:  Some("bh"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label:                Some("pl"),
            bind_group_layouts:   &[&bgl],
            push_constant_ranges: &[],
        });
        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label:  Some("bh_pipe"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module:              &shader,
                entry_point:        "vs",
                compilation_options: Default::default(),
                buffers:            &[],
            },
            fragment: Some(wgpu::FragmentState {
                module:              &shader,
                entry_point:        "fs",
                compilation_options: Default::default(),
                targets: &[Some(wgpu::ColorTargetState {
                    format,
                    blend:      Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive:    wgpu::PrimitiveState::default(),
            depth_stencil: None,
            multisample:  wgpu::MultisampleState::default(),
            multiview:    None,
            cache:        None,
        });

        Self { surface, device, queue, config, pipeline, bind_group, uni_buf }
    }

    fn resize(&mut self, w: u32, h: u32) {
        self.config.width  = w.max(1);
        self.config.height = h.max(1);
        self.surface.configure(&self.device, &self.config);
    }

    fn upload_params(&self, params: &Params) {
        self.queue.write_buffer(&self.uni_buf, 0, bytemuck::bytes_of(params));
    }

    fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        let frame = self.surface.get_current_texture()?;
        let view  = frame.texture.create_view(&Default::default());
        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut rp = enc.begin_render_pass(&wgpu::RenderPassDescriptor {
                label:                    Some("frame"),
                color_attachments:        &[Some(wgpu::RenderPassColorAttachment {
                    view:          &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load:  wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes:         None,
                occlusion_query_set:      None,
            });
            rp.set_pipeline(&self.pipeline);
            rp.set_bind_group(0, &self.bind_group, &[]);
            rp.draw(0..6, 0..1); // fullscreen quad (2 triangles = 6 vertices)
        }
        self.queue.submit([enc.finish()]);
        frame.present();
        Ok(())
    }
}

// ── Application ───────────────────────────────────────────────────────────────

struct App {
    window:       Option<Arc<Window>>,
    gpu:          Option<Gpu>,
    camera:       Camera,
    mouse_down:   bool,
    last_mouse:   PhysicalPosition<f64>,
    win_size:     (f32, f32), // (width, height) in physical pixels
}

impl App {
    fn new() -> Self {
        Self {
            window:     None,
            gpu:        None,
            camera:     Camera::default(),
            mouse_down: false,
            last_mouse: PhysicalPosition::new(0.0, 0.0),
            win_size:   (800.0, 800.0),
        }
    }

    fn update_title(&self) {
        let Some(win) = self.window.as_ref() else { return };
        let core = if self.camera.gutoe_core { "GUTOE r_c" } else { "GR" };
        let az_deg = self.camera.azimuth.to_degrees().rem_euclid(360.0);
        win.set_title(&format!(
            "GUTOE BH  |  inc {:.0}°  az {:.0}°  fov {:.1} r_s  disk {:.0} r_s  [{}]",
            self.camera.inclination, az_deg, self.camera.fov_rs, self.camera.disk_outer, core,
        ));
    }

    fn push_frame(&mut self) {
        let Some(gpu) = self.gpu.as_mut() else { return };
        let Some(win) = self.window.as_ref() else { return };
        let sz = win.inner_size();
        let params = self.camera.params(sz.width as f32, sz.height as f32);
        gpu.upload_params(&params);
        match gpu.render() {
            Ok(_) => {}
            Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                gpu.resize(sz.width, sz.height);
            }
            Err(e) => log::error!("render error: {e}"),
        }
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, el: &ActiveEventLoop) {
        let win = Arc::new(
            el.create_window(
                winit::window::WindowAttributes::default()
                    .with_title("GUTOE Black Hole — GRAND-151")
                    .with_inner_size(winit::dpi::LogicalSize::new(900u32, 900u32)),
            )
            .unwrap(),
        );
        let sz = win.inner_size();
        self.win_size = (sz.width.max(1) as f32, sz.height.max(1) as f32);
        let gpu = Gpu::new(Arc::clone(&win));
        self.window = Some(win);
        self.gpu    = Some(gpu);
        self.update_title();
    }

    fn window_event(&mut self, el: &ActiveEventLoop, _id: WindowId, event: WindowEvent) {
        match event {
            // ── Quit ──────────────────────────────────────────────────────────
            WindowEvent::CloseRequested => el.exit(),

            WindowEvent::KeyboardInput { event: KeyEvent { logical_key, state: ElementState::Pressed, .. }, .. } => {
                match logical_key {
                    Key::Named(NamedKey::Escape) => el.exit(),
                    Key::Named(NamedKey::ArrowUp) => {
                        // Arrow up = tilt toward edge-on
                        self.camera.inclination = (self.camera.inclination + 5.0).min(90.0);
                        self.update_title(); self.push_frame();
                    }
                    Key::Named(NamedKey::ArrowDown) => {
                        // Arrow down = tilt toward face-on
                        self.camera.inclination = (self.camera.inclination - 5.0).max(1.0);
                        self.update_title(); self.push_frame();
                    }
                    Key::Named(NamedKey::ArrowLeft) => {
                        self.camera.azimuth -= std::f32::consts::PI / 12.0; // 15° per press
                        self.update_title(); self.push_frame();
                    }
                    Key::Named(NamedKey::ArrowRight) => {
                        self.camera.azimuth += std::f32::consts::PI / 12.0;
                        self.update_title(); self.push_frame();
                    }
                    Key::Character(ref s) => match s.as_str() {
                        "q" | "Q" => el.exit(),
                        "r" | "R" => {
                            self.camera = Camera::default();
                            self.update_title(); self.push_frame();
                        }
                        "g" | "G" => {
                            self.camera.gutoe_core = !self.camera.gutoe_core;
                            log::info!("GUTOE lattice core: {}", self.camera.gutoe_core);
                            self.update_title(); self.push_frame();
                        }
                        "=" | "+" => {
                            self.camera.disk_outer = (self.camera.disk_outer + 1.0).min(30.0);
                            self.update_title(); self.push_frame();
                        }
                        "-" => {
                            self.camera.disk_outer = (self.camera.disk_outer - 1.0).max(3.5);
                            self.update_title(); self.push_frame();
                        }
                        _ => {}
                    },
                    _ => {}
                }
            }

            // ── Resize ────────────────────────────────────────────────────────
            WindowEvent::Resized(size) => {
                self.win_size = (size.width.max(1) as f32, size.height.max(1) as f32);
                if let Some(gpu) = self.gpu.as_mut() {
                    gpu.resize(size.width, size.height);
                    self.push_frame();
                }
            }

            // ── Mouse ─────────────────────────────────────────────────────────
            WindowEvent::MouseInput { button: MouseButton::Left, state, .. } => {
                self.mouse_down = state == ElementState::Pressed;
            }

            WindowEvent::CursorMoved { position, .. } => {
                if self.mouse_down {
                    let dx = (position.x - self.last_mouse.x) as f32;
                    let dy = (position.y - self.last_mouse.y) as f32;
                    let (w, h) = self.win_size;

                    // Vertical drag → inclination: full-height drag = 90° swing
                    // (screen-relative so it feels the same on any resolution/DPI)
                    self.camera.inclination =
                        (self.camera.inclination - dy * 90.0 / h).clamp(1.0, 90.0);

                    // Horizontal drag → azimuth: full-width drag = one full rotation
                    self.camera.azimuth -= dx * std::f32::consts::TAU / w;

                    self.update_title();
                    self.push_frame();
                }
                self.last_mouse = position;
            }

            // ── Scroll / zoom ─────────────────────────────────────────────────
            WindowEvent::MouseWheel { delta, .. } => {
                // Positive scroll = zoom in (fov shrinks), negative = zoom out.
                // Line delta: 1 notch ≈ 12% fov change — snappy but not jumpy.
                // Pixel delta (trackpad): scale so ~8px movement = same as one notch.
                let scroll = match delta {
                    MouseScrollDelta::LineDelta(_, y)   => y,
                    MouseScrollDelta::PixelDelta(p)     => p.y as f32 / 8.0,
                };
                self.camera.fov_rs = (self.camera.fov_rs * (1.0 - scroll * 0.12))
                    .clamp(3.0, 40.0);
                self.update_title();
                self.push_frame();
            }

            WindowEvent::RedrawRequested => self.push_frame(),
            _ => {}
        }
    }

    fn about_to_wait(&mut self, _el: &ActiveEventLoop) {
        if let Some(win) = self.window.as_ref() {
            win.request_redraw();
        }
    }
}

// ── Entry point ───────────────────────────────────────────────────────────────

fn main() {
    env_logger::init();

    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  GUTOE Black Hole Viewer — GRAND-151");
    println!("  GUTOE Schwarzschild metric + SC lattice (Cl(1,3)) + EHT-like rendering");
    println!("  Physics: Novikov-Thorne disk · Doppler beaming · gravitational redshift");
    println!("           Lensed star field · Hawking temperature · singularity-free core");
    println!("  Default: M87-like — 17° inclination, bright crescent at bottom");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
    println!("  Drag (vertical)    — inclination: 0°=face-on  90°=edge-on [↑↓ 5°/step]");
    println!("  Drag (horizontal)  — disk azimuth / rotate bright crescent  [←→ 15°/step]");
    println!("  Scroll             — zoom (field of view in units of r_s)");
    println!("  + / -              — disk outer radius (grow / shrink accretion disk)");
    println!("  G                  — toggle GUTOE lattice core r_c  (GR ↔ GUTOE)");
    println!("  R                  — reset to M87-like defaults");
    println!("  Q / Escape         — quit");
    println!("  Title bar shows live: inclination, azimuth, fov, disk size, mode");
    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");

    let ev = EventLoop::new().unwrap();
    ev.run_app(&mut App::new()).unwrap();
}
