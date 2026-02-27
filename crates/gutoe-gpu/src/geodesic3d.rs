// GUTOE 3D geodesic camera/ray reduction helpers
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// This module adds a true 3D camera/ray layer and maps those rays to the
// conserved quantities used by the current Schwarzschild tracer.
//
// NOTE:
//   The current `trace_photon` integrator is axisymmetric (r,phi) and does not
//   carry full 3D orientation through integration yet. This module gives us the
//   physically correct 3D ray construction and conserved-impact reduction that
//   we can evolve into a full 3D integrator.

use crate::metric::GutoeMetric;
use crate::tracer::{trace_photon, TraceResult};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    pub fn cross(self, o: Self) -> Self {
        Self {
            x: self.y * o.z - self.z * o.y,
            y: self.z * o.x - self.x * o.z,
            z: self.x * o.y - self.y * o.x,
        }
    }

    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    pub fn normalized(self) -> Option<Self> {
        let n = self.norm();
        if n <= 1e-15 {
            None
        } else {
            Some(Self::new(self.x / n, self.y / n, self.z / n))
        }
    }
}

impl std::ops::Add for Vec3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl std::ops::Sub for Vec3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl std::ops::Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct CameraFrame {
    pub position: Vec3,
    pub forward: Vec3,
    pub right: Vec3,
    pub up: Vec3,
    pub fov_y_rad: f64,
    pub aspect: f64,
}

impl CameraFrame {
    pub fn new(
        position: Vec3,
        forward: Vec3,
        world_up: Vec3,
        fov_y_rad: f64,
        aspect: f64,
    ) -> Option<Self> {
        let f = forward.normalized()?;
        let mut r = f.cross(world_up).normalized()?;
        let mut u = r.cross(f).normalized()?;
        // Re-orthogonalize for numerical stability.
        r = u.cross(f).normalized()?;
        u = r.cross(f).normalized()?;
        Some(Self {
            position,
            forward: f,
            right: r,
            up: u,
            fov_y_rad,
            aspect,
        })
    }

    /// True 3D pinhole ray from pixel center in NDC coordinates [-1,1]^2.
    pub fn ray_dir_from_ndc(&self, ndc_x: f64, ndc_y: f64) -> Vec3 {
        let tan_half = (0.5 * self.fov_y_rad).tan();
        let sx = ndc_x * self.aspect * tan_half;
        let sy = ndc_y * tan_half;
        (self.forward + self.right * sx + self.up * sy)
            .normalized()
            .unwrap_or(self.forward)
    }
}

/// Conserved-quantity reduction from a 3D null ray to current axisymmetric tracer
/// parameters `(bx, by)`.
#[derive(Debug, Clone, Copy)]
pub struct ReducedImpact {
    pub b: f64,
    pub bx: f64,
    pub by: f64,
    pub sin_i: f64,
}

/// Lean parity helper (Gutoe.Geodesic3DProjection.rayVec):
/// unnormalized pinhole ray through image-plane coordinates `(alpha, beta)`.
pub fn ray_vec(alpha: f64, beta: f64) -> Vec3 {
    Vec3::new(alpha, beta, 1.0)
}

/// Lean parity helper (impactRadiusSq).
pub fn impact_radius_sq(alpha: f64, beta: f64) -> f64 {
    alpha * alpha + beta * beta
}

/// Lean parity helper (rayNormSq_eval): `rayNormSq = impactRadiusSq + 1`.
pub fn ray_norm_sq(alpha: f64, beta: f64) -> f64 {
    let v = ray_vec(alpha, beta);
    v.dot(v)
}

/// Lean parity helper (rayDir + rayDir_unit_normSq).
pub fn ray_dir(alpha: f64, beta: f64) -> Option<Vec3> {
    ray_vec(alpha, beta).normalized()
}

/// Reduce 3D observer position `r` and unit direction `n` to Schwarzschild
/// impact invariants.
///
/// Angular momentum vector (up to E scale) is `L = r × n`.
/// - `b = |L|` (impact parameter in geometric-optics normalization E=1)
/// - orbital-plane inclination from equatorial plane:
///   `sin(i) = sqrt(1 - (Lz/|L|)^2)`
///
/// The current axisymmetric tracer only needs `b` and `sin(i)` (via `by/b`).
pub fn reduce_3d_to_axisym(observer_pos: Vec3, ray_dir_unit: Vec3) -> Option<ReducedImpact> {
    let n = ray_dir_unit.normalized()?;
    let l = observer_pos.cross(n);
    let b = l.norm();
    if b <= 1e-12 {
        return None;
    }
    let cos_i = (l.z / b).clamp(-1.0, 1.0);
    let sin_i = (1.0 - cos_i * cos_i).max(0.0).sqrt();
    let by = b * sin_i;
    // Keep sign from Lz so bx carries a stable orientation branch.
    let bx_mag = (b * b - by * by).max(0.0).sqrt();
    let bx = l.z.signum() * bx_mag;
    Some(ReducedImpact { b, bx, by, sin_i })
}

/// 3D-camera entry point that traces using reduced conserved quantities.
pub fn trace_photon_from_3d(
    metric: &GutoeMetric,
    disk_inner_re: f64,
    disk_outer_re: f64,
    observer_pos: Vec3,
    ray_dir_unit: Vec3,
    max_phi: f64,
    dphi: f64,
) -> TraceResult {
    let Some(red) = reduce_3d_to_axisym(observer_pos, ray_dir_unit) else {
        return TraceResult::Captured;
    };
    trace_photon(
        metric,
        disk_inner_re,
        disk_outer_re,
        red.bx,
        red.by,
        max_phi,
        dphi,
    )
}

fn orbit_accel_3d(r: f64, b: f64, r_s: f64, r_c: f64) -> f64 {
    let re2 = r * r + r_c * r_c;
    let re3 = re2 * re2.sqrt();
    r * (2.0 * r * r + r_c * r_c) / (b * b) - r + r_s * r * (r * r + 2.0 * r_c * r_c) / (2.0 * re3)
}

fn orbit_vr_sq_3d(r: f64, b: f64, r_s: f64, r_c: f64) -> f64 {
    let re2 = r * r + r_c * r_c;
    let re = re2.sqrt();
    let f = 1.0 - r_s / re;
    r * r * re2 / (b * b) - r * r * f
}

fn rk4_step_3d(r: f64, p: f64, b: f64, r_s: f64, r_c: f64, dphi: f64) -> (f64, f64) {
    let a = |ri: f64| orbit_accel_3d(ri, b, r_s, r_c);
    let k1r = p;
    let k1p = a(r);
    let k2r = p + 0.5 * dphi * k1p;
    let k2p = a(r + 0.5 * dphi * k1r);
    let k3r = p + 0.5 * dphi * k2p;
    let k3p = a(r + 0.5 * dphi * k2r);
    let k4r = p + dphi * k3p;
    let k4p = a(r + dphi * k3r);
    (
        r + dphi * (k1r + 2.0 * k2r + 2.0 * k3r + k4r) / 6.0,
        p + dphi * (k1p + 2.0 * k2p + 2.0 * k3p + k4p) / 6.0,
    )
}

/// True-3D world-space disk intersection using the same radial geodesic dynamics.
///
/// Geodesics in spherical symmetry lie in a plane; we integrate in that orbital
/// plane and lift each step back to world space to test `z=0` disk crossings.
pub fn trace_photon_3d_schwarzschild(
    metric: &GutoeMetric,
    disk_inner_re: f64,
    disk_outer_re: f64,
    observer_pos: Vec3,
    ray_dir_unit: Vec3,
    max_phi: f64,
    dphi: f64,
) -> TraceResult {
    let Some(d) = ray_dir_unit.normalized() else {
        return TraceResult::Captured;
    };
    let l = observer_pos.cross(d);
    let b = l.norm();
    if b < 1e-12 {
        return TraceResult::Captured;
    }
    let Some(n_hat) = l.normalized() else {
        return TraceResult::Captured;
    };

    // Closest-approach vector of the incoming asymptotic line.
    let p = observer_pos - d * observer_pos.dot(d);
    let b_line = p.norm();
    if b_line < 1e-12 {
        return TraceResult::Captured;
    }

    // Start from camera-side far point at r_start = 3b (same stability rule).
    let r_s = metric.r_s;
    let r_c = metric.r_core();
    let r_start = (3.0 * b_line).max(3.0 * b);
    let s = -(r_start * r_start - b_line * b_line).max(0.0).sqrt();
    let start_world = p + d * s;
    let Some(e_x) = start_world.normalized() else {
        return TraceResult::Captured;
    };
    let Some(e_y) = n_hat.cross(e_x).normalized() else {
        return TraceResult::Captured;
    };

    let vr0_sq = orbit_vr_sq_3d(r_start, b, r_s, r_c);
    let mut p_r = if vr0_sq > 0.0 {
        -vr0_sq.sqrt()
    } else {
        -r_start * r_start / b
    };
    let mut r = r_start;
    let mut phi = 0.0_f64;
    let mut turned = false;
    let mut n_cross = 0_u32;
    let r_capture_re = r_s * 0.99;
    let max_steps = (max_phi / dphi).ceil() as usize + 1;

    let world_at =
        |rad: f64, ang: f64| -> Vec3 { e_x * (rad * ang.cos()) + e_y * (rad * ang.sin()) };
    let mut w_prev = world_at(r, phi);

    for _ in 0..max_steps {
        let (r_new, p_rk4) = rk4_step_3d(r, p_r, b, r_s, r_c, dphi);
        let vr2_new = orbit_vr_sq_3d(r_new, b, r_s, r_c).max(0.0);
        let p_new = if p_rk4 >= 0.0 {
            vr2_new.sqrt()
        } else {
            -vr2_new.sqrt()
        };
        let phi_new = phi + dphi;
        let re_new = (r_new * r_new + r_c * r_c).sqrt();
        if !re_new.is_finite() || r_new.is_nan() {
            return TraceResult::Captured;
        }
        if re_new < r_capture_re || r_new < r_c * 0.01 {
            return TraceResult::Captured;
        }
        if !turned && p_r < 0.0 && p_new >= 0.0 {
            turned = true;
        }
        if turned && r_new >= r_start * 0.99 {
            return TraceResult::Escaped { phi_total: phi_new };
        }

        let w_new = world_at(r_new, phi_new);
        if w_prev.z * w_new.z <= 0.0 {
            let t = if (w_new.z - w_prev.z).abs() > 1e-15 {
                (-w_prev.z / (w_new.z - w_prev.z)).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let r_cross = r + t * (r_new - r);
            let re_cross = (r_cross * r_cross + r_c * r_c).sqrt();
            n_cross += 1;
            if re_cross >= disk_inner_re && re_cross <= disk_outer_re {
                return TraceResult::DiskHit {
                    r_eff: re_cross,
                    phi_orb: phi_new,
                    n_cross,
                };
            }
        }

        r = r_new;
        p_r = p_new;
        phi = phi_new;
        w_prev = w_new;
    }

    if r >= r_start * 0.5 {
        TraceResult::Escaped { phi_total: phi }
    } else {
        TraceResult::Captured
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_shift_keeps_impact_parameter() {
        let r0 = Vec3::new(0.0, 0.0, 100.0);
        let d = Vec3::new(0.3, -0.2, -1.0).normalized().unwrap();
        let r1 = r0 + d * 15.0; // same geometric line
        let a = reduce_3d_to_axisym(r0, d).unwrap();
        let b = reduce_3d_to_axisym(r1, d).unwrap();
        assert!((a.b - b.b).abs() < 1e-9, "b0={} b1={}", a.b, b.b);
    }

    #[test]
    fn camera_frame_emits_unit_rays() {
        let cam = CameraFrame::new(
            Vec3::new(0.0, 0.0, 30.0),
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, 1.0, 0.0),
            std::f64::consts::FRAC_PI_3,
            16.0 / 9.0,
        )
        .unwrap();
        for (x, y) in [
            (-1.0, -1.0),
            (-0.2, 0.6),
            (0.0, 0.0),
            (0.7, -0.4),
            (1.0, 1.0),
        ] {
            let d = cam.ray_dir_from_ndc(x, y);
            assert!((d.norm() - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn trace_3d_head_on_captures() {
        let m = GutoeMetric::new(1.0, 1.0);
        let r = trace_photon_3d_schwarzschild(
            &m,
            3.0,
            12.0,
            Vec3::new(0.0, 0.0, 60.0),
            Vec3::new(0.0, 0.0, -1.0),
            60.0,
            0.01,
        );
        assert!(matches!(r, TraceResult::Captured));
    }

    #[test]
    fn trace_3d_offaxis_nontrivial() {
        let m = GutoeMetric::new(1.0, 1.0);
        let r = trace_photon_3d_schwarzschild(
            &m,
            3.0,
            12.0,
            Vec3::new(0.0, 0.0, 60.0),
            Vec3::new(0.25, 0.08, -1.0),
            80.0,
            0.01,
        );
        assert!(
            matches!(r, TraceResult::Escaped { .. } | TraceResult::DiskHit { .. }),
            "expected off-axis ray to avoid direct capture, got {r:?}"
        );
    }

    #[test]
    fn lean_ray_norm_sq_eval_parity() {
        // Parity with theorem: Gutoe.Geodesic3DProjection.rayNormSq_eval
        for (a, b) in [
            (-2.0, -1.5),
            (-0.4, 0.7),
            (0.0, 0.0),
            (1.2, -0.8),
            (3.5, 4.0),
        ] {
            let lhs = ray_norm_sq(a, b);
            let rhs = impact_radius_sq(a, b) + 1.0;
            assert!((lhs - rhs).abs() < 1e-12, "a={a} b={b} lhs={lhs} rhs={rhs}");
        }
    }

    #[test]
    fn lean_impact_radius_even_beta_parity() {
        // Parity with theorem: Gutoe.Geodesic3DProjection.impactRadius_even_beta
        for (a, b) in [(-1.0, -2.0), (0.25, 0.9), (2.5, -0.3)] {
            let left = impact_radius_sq(a, -b).sqrt();
            let right = impact_radius_sq(a, b).sqrt();
            assert!(
                (left - right).abs() < 1e-12,
                "a={a} b={b} left={left} right={right}"
            );
        }
    }

    #[test]
    fn lean_ray_dir_unit_norm_parity() {
        // Parity with theorem: Gutoe.Geodesic3DProjection.rayDir_unit_normSq
        for (a, b) in [(-1.2, 0.4), (0.0, 0.0), (1.1, 2.2)] {
            let d = ray_dir(a, b).expect("ray_dir should normalize for finite inputs");
            assert!(
                (d.dot(d) - 1.0).abs() < 1e-12,
                "a={a} b={b} norm2={}",
                d.dot(d)
            );
            assert!(d.z > 0.0, "expected positive forward z component");
        }
    }
}
