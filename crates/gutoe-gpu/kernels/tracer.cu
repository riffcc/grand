/*
 * GUTOE Black Hole Ray Tracer — CUDA/HIP kernel
 * Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
 *
 * One GPU thread per pixel.  Ports trace_photon + pixel_color + star_field_color
 * from crates/gutoe-gpu/src/tracer.rs and src/bin/bh_render.rs.
 *
 * Compile (NVIDIA sm_89 = RTX 4070 Ti / Ada):
 *   nvcc -O3 -arch=sm_89 -Xcompiler -fPIC -c tracer.cu -o tracer_gpu.o
 *
 * For AMD ROCm (same source via HIP):
 *   hipcc -O3 -fPIC --offload-arch=gfx1151 -c tracer.cu -o tracer_gpu.o
 */

/* ── Portability ──────────────────────────────────────────────────────────── */
#ifdef __HIP_PLATFORM_AMD__
#  include <hip/hip_runtime.h>
#  define cudaError_t            hipError_t
#  define cudaSuccess            hipSuccess
#  define cudaGetErrorString     hipGetErrorString
#  define cudaMalloc             hipMalloc
#  define cudaFree               hipFree
#  define cudaMemcpy             hipMemcpy
#  define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#  define cudaDeviceSynchronize  hipDeviceSynchronize
#else
#  include <cuda_runtime.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TRACER_BLOCK 256
#define TRACER_PI    3.141592653589793238462643383279502884

#define TRACER_CUDA_CHECK(x) do { \
    cudaError_t _e = (x); \
    if (_e != cudaSuccess) { \
        fprintf(stderr, "CUDA error %s at %s:%d\n", \
                cudaGetErrorString(_e), __FILE__, __LINE__); \
        exit(1); \
    } \
} while(0)

/* ── Orbit equations (exact port of tracer.rs) ────────────────────────────── */

/* d²r/dφ² for a null geodesic in the GUTOE Schwarzschild metric.
 * Derivation: H = (dr/dφ)² = r²·r_eff²·(1/b² − V), d²r/dφ² = ½ dH/dr */
static __device__ __forceinline__
double bh_orbit_accel(double r, double b, double r_s, double r_c)
{
    double re2 = r*r + r_c*r_c;
    double re3 = re2 * sqrt(re2);
    return r*(2.0*r*r + r_c*r_c)/(b*b)
           - r
           + r_s*r*(r*r + 2.0*r_c*r_c)/(2.0*re3);
}

/* (dr/dφ)² = r²·r_eff²·(1/b² − V), V = (1 − r_s/r_eff)/r_eff² */
static __device__ __forceinline__
double bh_orbit_vr_sq(double r, double b, double r_s, double r_c)
{
    double re2 = r*r + r_c*r_c;
    double re  = sqrt(re2);
    double f   = 1.0 - r_s/re;
    return r*r*re2/(b*b) - r*r*f;
}

/* 4th-order Runge–Kutta step for (r, p = dr/dφ) */
static __device__ __forceinline__
void bh_rk4_step(double r, double p, double b, double r_s, double r_c, double dphi,
                  double* r_out, double* p_out)
{
    double k1r = p;
    double k1p = bh_orbit_accel(r, b, r_s, r_c);
    double k2r = p + 0.5*dphi*k1p;
    double k2p = bh_orbit_accel(r + 0.5*dphi*k1r, b, r_s, r_c);
    double k3r = p + 0.5*dphi*k2p;
    double k3p = bh_orbit_accel(r + 0.5*dphi*k2r, b, r_s, r_c);
    double k4r = p + dphi*k3p;
    double k4p = bh_orbit_accel(r + dphi*k3r, b, r_s, r_c);
    *r_out = r + dphi*(k1r + 2.0*k2r + 2.0*k3r + k4r)/6.0;
    *p_out = p + dphi*(k1p + 2.0*k2p + 2.0*k3p + k4p)/6.0;
}

/* Return codes for trace_photon_gpu */
#define TRACE_CAPTURED 0
#define TRACE_DISKHIT  1
#define TRACE_ESCAPED  2

/* Trace a null geodesic through the GUTOE Schwarzschild metric.
 *
 * Outputs:
 *   out_r_eff    — areal radius at disk hit (DISKHIT only)
 *   out_n_cross  — crossing number (DISKHIT only; 1 = direct, 2 = secondary, …)
 *   out_phi_total — total orbital angle (ESCAPED only)
 */
static __device__
int trace_photon_gpu(
    double r_s, double r_c,
    double disk_inner, double disk_outer,
    double bx, double by,
    double max_phi, double dphi,
    double* out_r_eff, unsigned int* out_n_cross, double* out_phi_total)
{
    double b = sqrt(bx*bx + by*by);
    if (b < 1e-12) return TRACE_CAPTURED;

    /* Deep-shadow shortcut: b_crit = (3√3/2)·r_s */
    double b_crit = 1.5 * 1.7320508075688772935 * r_s;
    if (b < 0.5 * b_crit) return TRACE_CAPTURED;

    double sin_i = by / b;
    int is_equatorial = (fabs(sin_i) < 1e-6);

    double r_start  = 3.0 * b;
    double vr0_sq   = bh_orbit_vr_sq(r_start, b, r_s, r_c);
    double p_start  = (vr0_sq > 0.0) ? -sqrt(vr0_sq) : -(r_start*r_start/b);

    double r_capture_re = r_s * 0.99;

    double r   = r_start;
    double p   = p_start;
    double phi = 0.0;
    unsigned int n_cross = 0;
    int turned = 0;

    int max_steps = (int)(max_phi / dphi) + 1;

    for (int step = 0; step < max_steps; step++) {
        double r_new, p_rk4;
        bh_rk4_step(r, p, b, r_s, r_c, dphi, &r_new, &p_rk4);

        /* Enforce orbital constraint: keep integration on the constraint surface */
        double vr2_new = bh_orbit_vr_sq(r_new, b, r_s, r_c);
        if (vr2_new < 0.0) vr2_new = 0.0;
        double p_new = (p_rk4 >= 0.0) ? sqrt(vr2_new) : -sqrt(vr2_new);
        double phi_new = phi + dphi;

        double re_new = sqrt(r_new*r_new + r_c*r_c);

        /* Capture: areal radius fell inside the horizon */
        if (re_new < r_capture_re || r_new < r_c * 0.01)
            return TRACE_CAPTURED;

        /* Detect turning point (ingoing → outgoing) */
        if (!turned && p < 0.0 && p_new >= 0.0)
            turned = 1;

        /* Escape: photon turned and returned to starting radius */
        if (turned && r_new >= r_start * 0.99) {
            *out_phi_total = phi_new;
            return TRACE_ESCAPED;
        }

        /* ── Disk hit detection ─────────────────────────────────────────── */
        if (is_equatorial) {
            /* Equatorial orbit: first inward crossing of the disk zone */
            double re_cur = sqrt(r*r + r_c*r_c);
            if (re_cur >= disk_inner && re_cur <= disk_outer && p < 0.0) {
                *out_r_eff   = re_cur;
                *out_n_cross = 1;
                return TRACE_DISKHIT;
            }
        } else {
            /* Tilted orbit: disk crossings at φ = n×π */
            double target = ((double)(n_cross) + 1.0) * TRACER_PI;
            if (phi < target && phi_new >= target) {
                double t       = (target - phi) / dphi;
                double r_cross = r + t * (r_new - r);
                double re_cross = sqrt(r_cross*r_cross + r_c*r_c);
                n_cross++;
                if (re_cross >= disk_inner && re_cross <= disk_outer) {
                    *out_r_eff   = re_cross;
                    *out_n_cross = n_cross;
                    return TRACE_DISKHIT;
                }
            }
        }

        r   = r_new;
        p   = p_new;
        phi = phi_new;
    }

    /* Ran out of steps — classify by current radius */
    if (r >= r_start * 0.5) {
        *out_phi_total = phi;
        return TRACE_ESCAPED;
    }
    return TRACE_CAPTURED;
}

/* ── Interior-camera trace: camera at r_cam < r_horizon, photons fire outward ─ */

/* Inside the horizon V(r) < 0, so orbit_vr_sq > 0 for ALL b.
 * Every photon can start moving outward.
 *   b < b_crit  → clears the photon-sphere barrier → ESCAPED (outside universe)
 *   b > b_crit  → turns around before the photon sphere → falls to GUTOE core
 * After turning, the photon falls back; we capture it once it returns
 * below the camera's areal radius (r_eff < r_eff_cam * 1.1). */
static __device__
int trace_photon_gpu_interior(
    double r_s, double r_c,
    double disk_inner, double disk_outer,
    double r_cam,
    double bx, double by,
    double max_phi, double dphi,
    double* out_r_eff, unsigned int* out_n_cross, double* out_phi_total)
{
    double b = sqrt(bx*bx + by*by);

    /* Pure-radial photon: b=0, no angular momentum → always escapes */
    if (b < 1e-12) {
        *out_phi_total = max_phi;
        return TRACE_ESCAPED;
    }

    /* From inside the horizon, orbit_vr_sq(r_cam) ≥ 0 always (V < 0) */
    double vr0_sq = bh_orbit_vr_sq(r_cam, b, r_s, r_c);
    if (vr0_sq < 0.0) vr0_sq = 0.0;
    double p_start = sqrt(vr0_sq);   /* outward */

    /* Camera areal radius: after turning, capture when re < this */
    double re_cam = sqrt(r_cam*r_cam + r_c*r_c);
    double re_cap = re_cam * 1.05;   /* 5 % buffer */

    /* Escape once the photon reaches the exterior region (≥ 3b from BH or
     * well past the disk outer edge — whichever is larger). */
    double r_escape = 3.0 * b;
    if (r_escape < disk_outer * 1.5) r_escape = disk_outer * 1.5;
    if (r_escape < 20.0 * r_s)       r_escape = 20.0 * r_s;

    double sin_i = by / b;
    int is_equatorial = (fabs(sin_i) < 1e-6);

    double r   = r_cam;
    double p   = p_start;
    double phi = 0.0;
    unsigned int n_cross = 0;
    int turned = 0;

    int max_steps = (int)(max_phi / dphi) + 1;

    for (int step = 0; step < max_steps; step++) {
        double r_new, p_rk4;
        bh_rk4_step(r, p, b, r_s, r_c, dphi, &r_new, &p_rk4);

        double vr2_new = bh_orbit_vr_sq(r_new, b, r_s, r_c);
        if (vr2_new < 0.0) vr2_new = 0.0;
        double p_new   = (p_rk4 >= 0.0) ? sqrt(vr2_new) : -sqrt(vr2_new);
        double phi_new = phi + dphi;
        double re_new  = sqrt(r_new*r_new + r_c*r_c);

        /* Detect turning point: p was positive (outward), now negative (inward) */
        if (!turned && p > 0.0 && p_new <= 0.0)
            turned = 1;

        /* After turning: capture once back below the camera areal radius */
        if (turned && re_new < re_cap)
            return TRACE_CAPTURED;

        /* Escape: photon has risen past the outer disk / strong-gravity region */
        if (!turned && r_new >= r_escape) {
            *out_phi_total = phi_new;
            return TRACE_ESCAPED;
        }

        /* Disk hit detection (same equatorial-crossing logic as exterior) */
        if (is_equatorial) {
            double re_cur = sqrt(r*r + r_c*r_c);
            if (!turned && re_cur >= disk_inner && re_cur <= disk_outer && p > 0.0) {
                *out_r_eff   = re_cur;
                *out_n_cross = 1;
                return TRACE_DISKHIT;
            }
        } else {
            double target = ((double)(n_cross) + 1.0) * TRACER_PI;
            if (phi < target && phi_new >= target) {
                double t        = (target - phi) / dphi;
                double r_cross  = r + t * (r_new - r);
                double re_cross = sqrt(r_cross*r_cross + r_c*r_c);
                n_cross++;
                if (re_cross >= disk_inner && re_cross <= disk_outer) {
                    *out_r_eff   = re_cross;
                    *out_n_cross = n_cross;
                    return TRACE_DISKHIT;
                }
            }
        }

        r   = r_new;
        p   = p_new;
        phi = phi_new;
    }

    /* Timed out */
    if (r >= r_escape * 0.5 && !turned) {
        *out_phi_total = phi;
        return TRACE_ESCAPED;
    }
    return TRACE_CAPTURED;
}

/* ── Interior-core trace: camera at r_cam < r_horizon, rays plunge toward core ─ */

static __device__
int trace_photon_gpu_interior_core(
    double r_s, double r_c,
    double r_cam,
    double bx, double by,
    double max_phi, double dphi,
    double* out_r_eff, double* out_phi_total)
{
    double b = sqrt(bx*bx + by*by);
    if (b < 1e-12) {
        *out_r_eff = r_c;
        *out_phi_total = 0.0;
        return TRACE_DISKHIT;
    }

    double vr0_sq = bh_orbit_vr_sq(r_cam, b, r_s, r_c);
    if (vr0_sq < 0.0) vr0_sq = 0.0;
    double p_start = -sqrt(vr0_sq);  /* inward */

    double re_core_cap = fmax(1.02 * r_c, r_c + 1e-9);

    double r = r_cam;
    double p = p_start;
    double phi = 0.0;
    int max_steps = (int)(max_phi / dphi) + 1;

    for (int step = 0; step < max_steps; step++) {
        double r_new, p_rk4;
        bh_rk4_step(r, p, b, r_s, r_c, dphi, &r_new, &p_rk4);

        double vr2_new = bh_orbit_vr_sq(r_new, b, r_s, r_c);
        if (vr2_new < 0.0) vr2_new = 0.0;
        double p_new   = (p_rk4 >= 0.0) ? sqrt(vr2_new) : -sqrt(vr2_new);
        double phi_new = phi + dphi;
        double re_new  = sqrt(r_new*r_new + r_c*r_c);

        if (re_new <= re_core_cap || r_new <= r_c * 0.01) {
            *out_r_eff = re_new;
            *out_phi_total = phi_new;
            return TRACE_DISKHIT;
        }

        /* Defensive: should not happen for inward branch, but keep sane output if it does. */
        if (p < 0.0 && p_new >= 0.0 && r_new > r_cam * 1.01) {
            *out_phi_total = phi_new;
            return TRACE_ESCAPED;
        }

        r = r_new;
        p = p_new;
        phi = phi_new;
    }

    *out_r_eff = re_core_cap;
    *out_phi_total = phi;
    return TRACE_DISKHIT;
}

/* ── Colour functions (port of bh_render.rs) ─────────────────────────────── */

/* False-colour the shadow interior by b/b_crit (interior_mode).
 * n ≈ −ln(1 − b/b_crit) / π   half-orbits before capture.
 * Cycles: orange→cyan→violet→warm-white per 4 half-orbits. */
static __device__
void bh_shadow_interior_color(double bx, double by, double r_s,
                                unsigned char* r_o, unsigned char* g_o, unsigned char* b_o)
{
    double b      = sqrt(bx*bx + by*by);
    double b_crit = 1.5 * 1.7320508075688772935 * r_s;
    double ratio  = b / b_crit;
    if (ratio >= 1.0) ratio = 1.0 - 1e-9;
    double n_float = -log(1.0 - ratio) / TRACER_PI;
    unsigned int n  = (unsigned int)n_float;
    float frac      = (float)(n_float - (double)n);

    float brightness = tanhf((float)(n + 1) * 0.5f) * (0.3f + 0.7f * frac);
    float rr, gg, bb;
    switch (n % 4) {
        case 0:  rr=1.00f; gg=0.55f; bb=0.10f; break;
        case 1:  rr=0.15f; gg=0.90f; bb=0.80f; break;
        case 2:  rr=0.75f; gg=0.20f; bb=1.00f; break;
        default: rr=1.00f; gg=0.95f; bb=0.60f; break;
    }
    *r_o = (unsigned char)fmin(fmax(rr * brightness * 255.0f, 0.0f), 255.0f);
    *g_o = (unsigned char)fmin(fmax(gg * brightness * 255.0f, 0.0f), 255.0f);
    *b_o = (unsigned char)fmin(fmax(bb * brightness * 255.0f, 0.0f), 255.0f);
}

/* GUTOE lattice-floor glow for the camera-inside view (b > b_crit photons).
 *
 * These are outward photons that turned around before the photon sphere.
 * The closer to b_crit (the photon ring edge) the brighter — those photons
 * orbited near r_ph many times and carry energy from the ring.
 * Deep inside the core (large excess) the glow fades to a dark amber. */
static __device__
void bh_gutoe_core_color(double b, double b_crit,
                          unsigned char* r_o, unsigned char* g_o, unsigned char* b_o)
{
    double excess = fmax(b - b_crit, 0.0) / b_crit;  /* 0 at ring edge */
    double glow   = exp(-excess * excess * 3.0);      /* Gaussian falloff */
    double L      = glow * 4.0;                       /* peak > 1 → Reinhard maps to ~0.8 */
    double bv     = L / (1.0 + L);                    /* Reinhard */

    /* Hot amber-orange palette: bright near b_crit, dim and red-orange further out */
    *r_o = (unsigned char)fmin(fmax(255.0 * pow(bv, 0.35), 0.0), 255.0);
    *g_o = (unsigned char)fmin(fmax(160.0 * pow(bv, 0.65), 0.0), 255.0);
    *b_o = (unsigned char)fmin(fmax( 30.0 * pow(bv, 2.0),  0.0), 255.0);
}

/* Core-facing interior palette derived only from traced geodesic invariants. */
static __device__
void bh_gutoe_core_physics_color(double b, double b_crit, double r_eff_hit, double phi_orb,
                                  double r_cam, double r_core,
                                  unsigned char* r_o, unsigned char* g_o, unsigned char* b_o)
{
    double eta = fmin(fmax(b / fmax(b_crit, 1e-9), 0.0), 2.0);
    double n_half = fmax(phi_orb / TRACER_PI, 0.0);
    double re_span = fmax(r_cam - r_core, 1e-9);
    double depth = fmin(fmax((r_eff_hit - r_core) / re_span, 0.0), 1.0);

    double d = eta - 1.0;
    double near_crit = exp(-(d * d) / 0.08);
    double winding = fmin(n_half / 12.0, 1.0);
    double plunge = 1.0 - depth;
    double luminance = fmax(0.10 + 0.75 * near_crit + 0.45 * winding + 0.35 * plunge, 0.0);
    double tone = luminance / (1.0 + luminance);

    *r_o = (unsigned char)fmin(fmax(255.0 * pow(tone, 0.36), 0.0), 255.0);
    *g_o = (unsigned char)fmin(fmax(170.0 * pow(tone, 0.62), 0.0), 255.0);
    *b_o = (unsigned char)fmin(fmax( 45.0 * pow(tone, 1.50), 0.0), 255.0);
}

/* splitmix64-style hash of two integer sky-grid coordinates */
static __device__ __forceinline__
unsigned long long bh_star_hash(long long x, long long y)
{
    unsigned long long s =
        ((unsigned long long)x * 0x9e3779b97f4a7c15ULL)
      + ((unsigned long long)y * 0x6c62272e07bb0142ULL);
    s ^= (s >> 30);
    s *= 0xbf58476d1ce4e5b9ULL;
    s ^= (s >> 27);
    s *= 0x94d049bb133111ebULL;
    s ^= (s >> 31);
    return s;
}

/* RGB for escaped (background) ray — hashed on sky coordinates */
static __device__
void bh_star_field_color(double bx, double by, double phi_total,
                          unsigned char* r_o, unsigned char* g_o, unsigned char* b_o)
{
    /* Rotate impact vector (bx, by) by phi_total to get sky direction.
     * phi_total ≈ π   → direct ray (weakly deflected background star)
     * phi_total ≈ 2π  → Einstein ring copy (same stars)
     * phi_total ≫ π   → photon ring regime (rapidly changing backgrounds) */
    double sky_x = bx * cos(phi_total) - by * sin(phi_total);
    double sky_y = bx * sin(phi_total) + by * cos(phi_total);

    const double SKY_SCALE = 50.0;
    long long hx = (long long)round(sky_x * SKY_SCALE);
    long long hy = (long long)round(sky_y * SKY_SCALE);
    unsigned long long h = bh_star_hash(hx, hy);

    /* ~1.5 % star density (1000/65536) */
    if ((h & 0xFFFFULL) >= 1000ULL) {
        unsigned char v = (unsigned char)((h >> 40) & 7);
        *r_o = v >> 2;
        *g_o = v >> 2;
        *b_o = (v >> 1) + 8;
        return;
    }

    unsigned int bright_raw = (unsigned int)((h >> 16) & 0xFF);
    unsigned char bright    = (unsigned char)(80u + bright_raw * 175u / 255u);

    unsigned int spec = (unsigned int)((h >> 24) & 0xF);
    if (spec <= 1) {
        /* Deep red — cool M-dwarf */
        *r_o = bright/2; *g_o = bright/5; *b_o = bright/10;
    } else if (spec <= 4) {
        /* Orange — K-star */
        *r_o = bright;
        *g_o = (unsigned char)((unsigned int)bright*68u/100u);
        *b_o = (unsigned char)((unsigned int)bright*38u/100u);
    } else if (spec <= 9) {
        /* Warm white — F/G sun-like */
        *r_o = bright; *g_o = bright;
        *b_o = (unsigned char)((unsigned int)bright*88u/100u);
    } else if (spec <= 12) {
        /* Pure white — A-star */
        *r_o = bright; *g_o = bright; *b_o = bright;
    } else {
        /* Blue-white — hot B/O star */
        unsigned char rg = (unsigned char)((unsigned int)bright*82u/100u);
        *r_o = rg; *g_o = rg; *b_o = bright;
    }
}

/* Colour by photon ring order (ring_mode) */
static __device__ __forceinline__
void bh_ring_order_color(unsigned int n_cross,
                          unsigned char* r_o, unsigned char* g_o, unsigned char* b_o)
{
    switch (n_cross) {
        case 1:  *r_o=255; *g_o=165; *b_o= 40; return;   /* orange — direct */
        case 2:  *r_o= 40; *g_o=220; *b_o=200; return;   /* cyan — secondary */
        case 3:  *r_o=190; *g_o= 60; *b_o=255; return;   /* purple — tertiary */
        default: *r_o=100; *g_o=100; *b_o=100; return;   /* grey — higher order */
    }
}

/* Novikov–Thorne temperature colour with Reinhard tone mapping */
static __device__
void bh_pixel_color(double r_eff, double r_isco, double r_outer, double r_s,
                     double bx_raw, double sin_inc,
                     unsigned int n_cross, int doppler, int ring_mode,
                     unsigned char* r_o, unsigned char* g_o, unsigned char* b_o)
{
    if (ring_mode) { bh_ring_order_color(n_cross, r_o, g_o, b_o); return; }

    /* Novikov–Thorne temperature: T ∝ (r_ISCO/r_eff)^{3/4} */
    double t_rel = pow(r_isco / r_eff, 0.75);

    /* Smooth outer disk taper: avoids hard edge at disk_outer */
    double excess = fmax(r_eff - r_outer, 0.0) / (0.5 * r_outer);
    double outer_taper = exp(-excess*excess);

    /* Higher-order images dimmer */
    double fade = pow(0.65, (double)((int)n_cross - 1));

    /* Relativistic Keplerian Doppler D⁴ */
    double doppler_d4 = 1.0;
    if (doppler) {
        double r_safe  = fmax(r_eff, 1e-12);
        double beta    = fmin(sqrt(r_s / (2.0 * r_safe)), 0.5);
        double beta_obs = beta * sin_inc * fmax(fmin(bx_raw / r_safe, 1.0), -1.0);
        double D = 1.0 / (1.0 - beta_obs);
        doppler_d4 = fmax(fmin(D*D*D*D, 200.0), 0.01);
    }

    /* Reinhard tone mapping */
    double luminance = fmax(t_rel * fade * doppler_d4 * outer_taper, 0.0);
    double bv = luminance / (1.0 + luminance);

    /* Orange-white thermal palette */
    *r_o = (unsigned char)fmin(fmax(255.0 * pow(bv, 0.35), 0.0), 255.0);
    *g_o = (unsigned char)fmin(fmax(210.0 * pow(bv, 0.60), 0.0), 255.0);
    *b_o = (unsigned char)fmin(fmax(130.0 * pow(bv, 1.60), 0.0), 255.0);
}

/* ── Render kernel: one thread per pixel ─────────────────────────────────── */

__global__
void bh_render_kernel(
    int width, int height,
    double fov_rs, double sin_inc,
    double r_s, double r_c,
    double disk_inner, double disk_outer,
    double max_phi, double dphi,
    double az_cos, double az_sin,
    int doppler, int ring_mode, int interior_mode, int core_look_mode,
    double r_cam,                   /* 0.0 = exterior; >0 = interior camera */
    unsigned char* pixels)          /* output: width*height*3 RGB bytes */
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= width * height) return;

    int iy = idx / width;
    int ix = idx % width;

    /* Camera coordinates — orthographic, y-axis upward */
    double scale  = 2.0 * fov_rs * r_s / (double)width;
    double sx     = ((double)ix - 0.5 * ((double)width  - 1.0)) * scale;
    double sy     = (0.5 * ((double)height - 1.0) - (double)iy) * scale;
    double bx_raw = sx;
    double by_raw = sy * sin_inc;

    /* Apply azimuth rotation in the screen plane */
    double bx = az_cos * bx_raw - az_sin * by_raw;
    double by = az_sin * bx_raw + az_cos * by_raw;

    double r_eff = 0.0, phi_total = 0.0;
    unsigned int n_cross = 0;
    int result;

    if (r_cam > 0.0) {
        if (core_look_mode) {
            /* Interior camera looking down toward the core. */
            result = trace_photon_gpu_interior_core(
                r_s, r_c, r_cam, bx, by, max_phi, dphi, &r_eff, &phi_total);
        } else {
            /* Interior camera: fire photons outward from r_cam */
            result = trace_photon_gpu_interior(r_s, r_c, disk_inner, disk_outer,
                                                r_cam, bx, by, max_phi, dphi,
                                                &r_eff, &n_cross, &phi_total);
        }
    } else {
        result = trace_photon_gpu(r_s, r_c, disk_inner, disk_outer, bx, by,
                                   max_phi, dphi, &r_eff, &n_cross, &phi_total);
    }

    unsigned char r_px, g_px, b_px;
    double r_isco = 3.0 * r_s;
    double b_crit = 1.5 * 1.7320508075688772935 * r_s;

    switch (result) {
        case TRACE_CAPTURED:
            if (r_cam > 0.0) {
                /* Photon turned around before photon sphere → GUTOE core glow */
                double b_mag = sqrt(bx*bx + by*by);
                bh_gutoe_core_color(b_mag, b_crit, &r_px, &g_px, &b_px);
            } else if (interior_mode) {
                /* False-colour shadow: colour by half-orbit count before capture */
                bh_shadow_interior_color(bx, by, r_s, &r_px, &g_px, &b_px);
            } else {
                r_px = 0; g_px = 0; b_px = 0;
            }
            break;
        case TRACE_DISKHIT:
            if (r_cam > 0.0 && core_look_mode) {
                double b_mag = sqrt(bx*bx + by*by);
                bh_gutoe_core_physics_color(b_mag, b_crit, r_eff, phi_total, r_cam, r_c, &r_px, &g_px, &b_px);
            } else {
                bh_pixel_color(r_eff, r_isco, disk_outer, r_s, bx_raw, sin_inc,
                                n_cross, doppler, ring_mode, &r_px, &g_px, &b_px);
            }
            break;
        default: /* TRACE_ESCAPED */
            bh_star_field_color(bx, by, phi_total, &r_px, &g_px, &b_px);
            break;
    }

    pixels[idx * 3 + 0] = r_px;
    pixels[idx * 3 + 1] = g_px;
    pixels[idx * 3 + 2] = b_px;
}

/* ── Host driver ─────────────────────────────────────────────────────────── */

extern "C"
void gutoe_render_bh(
    int width, int height,
    double fov_rs, double inclination_deg,
    double r_s, double r_c,
    double disk_inner_rs, double disk_outer_rs,
    double max_phi, double dphi,
    double az_deg,
    int doppler, int ring_mode, int interior_mode, int core_look_mode,
    double r_cam_rs,                /* 0.0 = exterior; >0 = interior camera at r_cam_rs * r_s */
    unsigned char* out_pixels)      /* host output buffer: width*height*3 */
{
    double sin_inc    = sin(inclination_deg * TRACER_PI / 180.0);
    double az_rad     = az_deg * TRACER_PI / 180.0;
    double az_cos     = cos(az_rad);
    double az_sin     = sin(az_rad);
    double disk_inner = disk_inner_rs * r_s;
    double disk_outer = disk_outer_rs * r_s;
    double r_cam      = r_cam_rs * r_s;

    int n = width * height;
    unsigned char* d_pixels;
    TRACER_CUDA_CHECK(cudaMalloc(&d_pixels, (size_t)n * 3));

    int n_blk = (n + TRACER_BLOCK - 1) / TRACER_BLOCK;
    bh_render_kernel<<<n_blk, TRACER_BLOCK>>>(
        width, height, fov_rs, sin_inc,
        r_s, r_c, disk_inner, disk_outer,
        max_phi, dphi, az_cos, az_sin,
        doppler, ring_mode, interior_mode, core_look_mode, r_cam,
        d_pixels);

    TRACER_CUDA_CHECK(cudaDeviceSynchronize());
    TRACER_CUDA_CHECK(cudaMemcpy(out_pixels, d_pixels, (size_t)n * 3,
                                  cudaMemcpyDeviceToHost));
    cudaFree(d_pixels);
}
