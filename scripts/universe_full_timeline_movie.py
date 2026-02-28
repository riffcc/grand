#!/usr/bin/env python3
"""
Render a 1080p full cosmic timeline movie from ~t=0 to today,
including a star-formation visual lane.

Outputs (default /tmp/bh_renders):
- universe_full_timeline_1080p.mp4
- universe_full_timeline_1080p.gif
- universe_full_timeline_summary.json
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import shutil
import subprocess
import tempfile

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

# GUTOE constants / ratios mirrored from crates/gutoe-physics/src/constants.rs
C = 299_792_458.0
METER_PER_MPC = 3.085_677_581_491_367e22
SEC_PER_YEAR = 31_557_600.0
SEC_PER_GYR = 1.0e9 * SEC_PER_YEAR
PLANCK_TIME = 5.391247e-44
T_CMB0_K = 2.7255
OMEGA_B0 = 0.0493
DARK_TO_VISIBLE_GEOMETRIC_RATIO = 60.0 / 11.0
OMEGA_DM0 = OMEGA_B0 * DARK_TO_VISIBLE_GEOMETRIC_RATIO
OMEGA_M0 = OMEGA_B0 + OMEGA_DM0
OMEGA_R0 = 9.0e-5
OMEGA_K0 = 0.0
LAMBDA_FULL = 1.1033233175245378e-52

# Reference milestones (seconds)
T_BBN = 180.0
T_RECOMB = 380_000.0 * SEC_PER_YEAR
T_FIRST_STARS = 180_000_000.0 * SEC_PER_YEAR
T_SFR_PEAK = 3_500_000_000.0 * SEC_PER_YEAR


def h0_from_lambda_and_omega_lambda(lambda_full: float, omega_lambda0: float) -> float:
    h0_s_inv = C * math.sqrt(lambda_full / (3.0 * omega_lambda0))
    return h0_s_inv * METER_PER_MPC / 1000.0


def h_of_a_sinv(a: np.ndarray, h0_km_s_mpc: float, omega_r0: float, omega_m0: float, omega_k0: float, omega_lambda0: float) -> np.ndarray:
    h0_s_inv = h0_km_s_mpc * 1000.0 / METER_PER_MPC
    e2 = omega_r0 / (a**4) + omega_m0 / (a**3) + omega_k0 / (a**2) + omega_lambda0
    return h0_s_inv * np.sqrt(np.maximum(e2, 1e-300))


def omega_components_at_a(a: np.ndarray, omega_r0: float, omega_m0: float, omega_k0: float, omega_lambda0: float):
    er = omega_r0 / (a**4)
    em = omega_m0 / (a**3)
    ek = omega_k0 / (a**2)
    el = np.full_like(a, omega_lambda0)
    et = er + em + ek + el
    return er / et, em / et, el / et


def build_t_a_table(h0_km_s_mpc: float, omega_r0: float, omega_m0: float, omega_k0: float, omega_lambda0: float):
    a = np.geomspace(1e-22, 1.0, 450_000)
    h = h_of_a_sinv(a, h0_km_s_mpc, omega_r0, omega_m0, omega_k0, omega_lambda0)
    dt_da = 1.0 / np.maximum(a * h, 1e-300)
    da = np.diff(a)
    dt_mid = 0.5 * (dt_da[:-1] + dt_da[1:]) * da
    t = np.empty_like(a)
    t[0] = 0.0
    t[1:] = np.cumsum(dt_mid)
    return t, a


def sfr_proxy_madau_dickinson(z: np.ndarray) -> np.ndarray:
    # shape-only proxy; normalized later
    zp1 = 1.0 + np.maximum(z, 0.0)
    num = zp1**2.7
    den = 1.0 + (zp1 / 2.9) ** 5.6
    return num / np.maximum(den, 1e-30)


def gaussian_blobs(h: int, w: int, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    y = np.linspace(-1.0, 1.0, h)
    x = np.linspace(-1.0, 1.0, w)
    X, Y = np.meshgrid(x, y)
    img = np.zeros((h, w), dtype=float)
    for _ in range(n):
        cx, cy = rng.uniform(-1.0, 1.0, size=2)
        sx, sy = rng.uniform(0.08, 0.35, size=2)
        amp = rng.uniform(0.3, 1.0)
        img += amp * np.exp(-(((X - cx) / sx) ** 2 + ((Y - cy) / sy) ** 2))
    img -= img.min()
    img /= max(img.max(), 1e-9)
    return img


def stage_weights(t_cur: float):
    logt = math.log10(max(t_cur, PLANCK_TIME))

    # Smooth sigmoids in log-time
    def s(x, c, w):
        return 1.0 / (1.0 + math.exp(-(x - c) / w))

    plasma = 1.0 - s(logt, math.log10(T_RECOMB), 0.25)
    cmb = s(logt, math.log10(T_RECOMB), 0.20) * (1.0 - s(logt, math.log10(T_FIRST_STARS), 0.25))
    stars = s(logt, math.log10(T_FIRST_STARS), 0.20)

    norm = plasma + cmb + stars
    if norm <= 0.0:
        return 0.0, 0.0, 1.0
    return plasma / norm, cmb / norm, stars / norm


def render_frame(
    frame_path: pathlib.Path,
    i: int,
    nframes: int,
    t_start: float,
    t_end: float,
    t_table: np.ndarray,
    a_table: np.ndarray,
    h0_km_s_mpc: float,
    omega_r0: float,
    omega_m0: float,
    omega_k0: float,
    omega_lambda0: float,
    z_hist: np.ndarray,
    sfr_hist: np.ndarray,
    sfr_cum: np.ndarray,
    stars_x: np.ndarray,
    stars_y: np.ndarray,
    stars_size: np.ndarray,
    stars_birth: np.ndarray,
    plasma_tex: np.ndarray,
    cmb_tex: np.ndarray,
    nebula_tex: np.ndarray,
):
    u = 0.0 if nframes <= 1 else i / (nframes - 1)
    logt = math.log10(t_start) + u * (math.log10(t_end) - math.log10(t_start))
    t_cur = 10.0**logt

    a_cur = float(np.interp(t_cur, t_table, a_table))
    z_cur = 1.0 / max(a_cur, 1e-300) - 1.0
    temp_cur = T_CMB0_K / max(a_cur, 1e-300)
    h_cur = float(h_of_a_sinv(np.array([a_cur]), h0_km_s_mpc, omega_r0, omega_m0, omega_k0, omega_lambda0)[0])
    h_cur_km_s_mpc = h_cur * METER_PER_MPC / 1000.0
    omega_r_cur, omega_m_cur, omega_l_cur = [
        float(x[0])
        for x in omega_components_at_a(np.array([a_cur]), omega_r0, omega_m0, omega_k0, omega_lambda0)
    ]

    # local stage blending
    w_plasma, w_cmb, w_stars = stage_weights(t_cur)
    sf = float(np.interp(z_cur, z_hist[::-1], sfr_hist[::-1]))
    sf_cum = float(np.interp(z_cur, z_hist[::-1], sfr_cum[::-1]))

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    gs = GridSpec(2, 3, width_ratios=[1.5, 1.0, 1.0], height_ratios=[1.0, 1.0], figure=fig)

    # Left cinematic panel
    ax0 = fig.add_subplot(gs[:, 0])
    ax0.set_facecolor("black")
    ax0.imshow(plasma_tex, cmap="inferno", alpha=np.clip(0.95 * w_plasma, 0.0, 1.0), origin="lower")
    ax0.imshow(cmb_tex, cmap="coolwarm", alpha=np.clip(0.75 * w_cmb, 0.0, 1.0), origin="lower")

    # stars + nebula become visible in late universe
    n_vis = int(sf_cum * len(stars_x))
    if n_vis > 0:
        boost = 0.5 + 1.2 * np.clip(sf / (sfr_hist.max() + 1e-9), 0.0, 1.0)
        ax0.imshow(nebula_tex, cmap="magma", alpha=np.clip(0.20 * w_stars * boost, 0.0, 0.45), origin="lower")
        ax0.scatter(
            stars_x[:n_vis],
            stars_y[:n_vis],
            s=stars_size[:n_vis] * (0.7 + 0.6 * w_stars),
            c="white",
            alpha=np.clip(0.25 + 0.7 * w_stars, 0.0, 0.95),
            linewidths=0,
        )

    ax0.set_axis_off()
    ax0.set_title("Cosmic Visual Lane: Plasma -> CMB -> Dark Ages -> Stars", fontsize=15, pad=10)

    # Top-right: scale factor and temperature history
    ax1 = fig.add_subplot(gs[0, 1:])
    t_plot = np.geomspace(max(t_start, 1e-40), t_end, 1600)
    a_plot = np.interp(t_plot, t_table, a_table)
    z_plot = 1.0 / np.maximum(a_plot, 1e-300) - 1.0
    temp_plot = T_CMB0_K / np.maximum(a_plot, 1e-300)

    ax1.plot(t_plot, temp_plot, color="#ff7f11", lw=2.0, label="T(t)")
    ax1.scatter([t_cur], [temp_cur], color="#ffd166", s=45, zorder=6)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("time [s]")
    ax1.set_ylabel("temperature [K]")
    ax1.grid(alpha=0.25)

    ax1b = ax1.twinx()
    ax1b.plot(t_plot, z_plot, color="#48bfe3", lw=1.5, alpha=0.9, label="z(t)")
    ax1b.scatter([t_cur], [z_cur], color="#90e0ef", s=30, zorder=6)
    ax1b.set_yscale("log")
    ax1b.set_ylabel("redshift z")

    # Bottom-middle: density fractions
    ax2 = fig.add_subplot(gs[1, 1])
    omr_plot, omm_plot, oml_plot = omega_components_at_a(a_plot, omega_r0, omega_m0, omega_k0, omega_lambda0)
    ax2.plot(t_plot, omr_plot, color="#f94144", lw=2.0, label=r"$\Omega_r$")
    ax2.plot(t_plot, omm_plot, color="#43aa8b", lw=2.0, label=r"$\Omega_m$")
    ax2.plot(t_plot, oml_plot, color="#577590", lw=2.0, label=r"$\Omega_\Lambda$")
    ax2.scatter([t_cur], [omega_r_cur], color="#f94144", s=26)
    ax2.scatter([t_cur], [omega_m_cur], color="#43aa8b", s=26)
    ax2.scatter([t_cur], [omega_l_cur], color="#577590", s=26)
    ax2.set_xscale("log")
    ax2.set_xlabel("time [s]")
    ax2.set_ylabel("density fraction")
    ax2.set_ylim(-0.02, 1.02)
    ax2.grid(alpha=0.25)
    ax2.legend(loc="upper right", fontsize=10)

    # Bottom-right: star formation lane + telemetry
    ax3 = fig.add_subplot(gs[1, 2])
    t_from_z = np.interp(z_hist[::-1], z_plot[::-1], t_plot[::-1], left=t_plot[0], right=t_plot[-1])
    ax3.plot(t_from_z, sfr_hist[::-1] / (sfr_hist.max() + 1e-9), color="#bc6c25", lw=2.0, label="SFR proxy (norm)")
    ax3.plot(t_from_z, sfr_cum[::-1], color="#e9c46a", lw=2.0, label="cumulative stars")
    ax3.scatter([t_cur], [sf / (sfr_hist.max() + 1e-9)], color="#f4a261", s=30)
    ax3.scatter([t_cur], [sf_cum], color="#e9c46a", s=30)
    ax3.set_xscale("log")
    ax3.set_ylim(-0.02, 1.02)
    ax3.set_xlabel("time [s]")
    ax3.set_ylabel("normalized")
    ax3.grid(alpha=0.25)
    ax3.legend(loc="upper left", fontsize=9)

    # Text overlay
    if t_cur < T_BBN:
        era = "Primordial plasma"
    elif t_cur < T_RECOMB:
        era = "Nucleosynthesis / radiation era"
    elif t_cur < T_FIRST_STARS:
        era = "Post-recombination dark ages"
    elif t_cur < T_SFR_PEAK:
        era = "First stars and galaxy assembly"
    else:
        era = "Mature structure / late universe"

    text = (
        f"{era}\n"
        f"t = {t_cur: .3e} s ({t_cur/SEC_PER_YEAR: .3e} yr)\n"
        f"z = {z_cur: .3e}\n"
        f"T = {temp_cur: .3e} K\n"
        f"H = {h_cur_km_s_mpc: .3e} km/s/Mpc\n"
        f"star_fraction ~= {sf_cum: .3f}"
    )
    ax3.text(0.02, 0.98, text, transform=ax3.transAxes, va="top", ha="left", fontsize=10, family="monospace")

    fig.suptitle("GUTOE Full Timeline Simulation — from t~0 to today", fontsize=18, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(frame_path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/tmp/bh_renders")
    ap.add_argument("--clip-seconds", type=float, default=30.0)
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--skip-gif", action="store_true")
    args = ap.parse_args()

    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    omega_lambda0 = 1.0 - OMEGA_M0 - OMEGA_R0 - OMEGA_K0
    h0_km_s_mpc = h0_from_lambda_and_omega_lambda(LAMBDA_FULL, omega_lambda0)

    t_table, a_table = build_t_a_table(h0_km_s_mpc, OMEGA_R0, OMEGA_M0, OMEGA_K0, omega_lambda0)
    t_start = max(PLANCK_TIME, 1e-36)
    t_end = float(t_table[-1])

    # star-formation proxy in z-space
    z_hist = np.geomspace(1.0, 81.0, 2500) - 1.0
    sfr_hist = sfr_proxy_madau_dickinson(z_hist)
    x = np.log1p(z_hist)
    dx = np.diff(x)
    trap = 0.5 * (sfr_hist[:-1] + sfr_hist[1:]) * dx
    sfr_cum = np.empty_like(sfr_hist)
    sfr_cum[0] = 0.0
    sfr_cum[1:] = np.cumsum(trap)
    sfr_cum -= sfr_cum.min()
    sfr_cum /= max(sfr_cum.max(), 1e-12)

    rng = np.random.default_rng(7)
    nstars = 3500
    stars_x = rng.uniform(0, 1, size=nstars)
    stars_y = rng.uniform(0, 1, size=nstars)
    stars_size = rng.uniform(1.0, 9.0, size=nstars)
    stars_birth = np.sort(rng.uniform(0.0, 1.0, size=nstars))

    h, w = 700, 1050
    plasma_tex = gaussian_blobs(h, w, n=55, seed=3)
    cmb_tex = gaussian_blobs(h, w, n=26, seed=11)
    nebula_tex = gaussian_blobs(h, w, n=18, seed=19)

    nframes = max(2, int(round(args.clip_seconds * args.fps)))
    frame_dir = pathlib.Path(tempfile.mkdtemp(prefix="universe_full_frames_"))

    try:
        for i in range(nframes):
            render_frame(
                frame_dir / f"frame_{i:05d}.png",
                i,
                nframes,
                t_start,
                t_end,
                t_table,
                a_table,
                h0_km_s_mpc,
                OMEGA_R0,
                OMEGA_M0,
                OMEGA_K0,
                omega_lambda0,
                z_hist,
                sfr_hist,
                sfr_cum,
                stars_x,
                stars_y,
                stars_size,
                stars_birth,
                plasma_tex,
                cmb_tex,
                nebula_tex,
            )

        ffmpeg = shutil.which("ffmpeg")
        if ffmpeg is None:
            raise RuntimeError("ffmpeg not found in PATH")

        mp4 = out_dir / "universe_full_timeline_1080p.mp4"
        subprocess.run(
            [
                ffmpeg,
                "-y",
                "-framerate",
                str(args.fps),
                "-i",
                str(frame_dir / "frame_%05d.png"),
                "-c:v",
                "libx264",
                "-preset",
                "slow",
                "-crf",
                "18",
                "-pix_fmt",
                "yuv420p",
                str(mp4),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        gif = out_dir / "universe_full_timeline_1080p.gif"
        if not args.skip_gif:
            palette = frame_dir / "palette.png"
            subprocess.run(
                [ffmpeg, "-y", "-i", str(frame_dir / "frame_%05d.png"), "-vf", "palettegen", str(palette)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            subprocess.run(
                [
                    ffmpeg,
                    "-y",
                    "-framerate",
                    str(args.fps),
                    "-i",
                    str(frame_dir / "frame_%05d.png"),
                    "-i",
                    str(palette),
                    "-lavfi",
                    "paletteuse",
                    str(gif),
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

        summary = {
            "clip_seconds": args.clip_seconds,
            "fps": args.fps,
            "frames": nframes,
            "time_start_s": t_start,
            "time_end_s": t_end,
            "time_end_gyr": t_end / SEC_PER_GYR,
            "h0_km_s_mpc": h0_km_s_mpc,
            "omega_b0": OMEGA_B0,
            "omega_dm0": OMEGA_DM0,
            "omega_m0": OMEGA_M0,
            "omega_r0": OMEGA_R0,
            "omega_lambda0": omega_lambda0,
            "lambda_full": LAMBDA_FULL,
            "artifacts": {
                "mp4": str(mp4),
                "gif": None if args.skip_gif else str(gif),
            },
        }
        summary_path = out_dir / "universe_full_timeline_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

        print("wrote", mp4)
        if not args.skip_gif:
            print("wrote", gif)
        print("wrote", summary_path)
    finally:
        shutil.rmtree(frame_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
