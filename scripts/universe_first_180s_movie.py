#!/usr/bin/env python3
"""
Render a 1080p movie of the first 180 seconds of cosmic evolution
from the GUTOE FRW lane (derived Λ + matter fractions).

Default outputs:
- /tmp/bh_renders/universe_first_180s_1080p.mp4
- /tmp/bh_renders/universe_first_180s_1080p.gif
- /tmp/bh_renders/universe_first_180s_summary.json
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
T_CMB0_K = 2.7255
OMEGA_B0 = 0.0493
DARK_TO_VISIBLE_GEOMETRIC_RATIO = 60.0 / 11.0
OMEGA_DM0 = OMEGA_B0 * DARK_TO_VISIBLE_GEOMETRIC_RATIO
OMEGA_M0 = OMEGA_B0 + OMEGA_DM0
OMEGA_R0 = 9.0e-5
OMEGA_K0 = 0.0
LAMBDA_FULL = 1.1033233175245378e-52


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


def build_time_to_scale_table(t_max: float, h0_km_s_mpc: float, omega_r0: float, omega_m0: float, omega_k0: float, omega_lambda0: float):
    # Radiation-era estimate for a(t_max), padded up for safety.
    h0_s_inv = h0_km_s_mpc * 1000.0 / METER_PER_MPC
    a_est = math.sqrt(max(1e-40, 2.0 * h0_s_inv * math.sqrt(omega_r0) * t_max))

    a_min = 1e-16
    a_max = max(1e-6, 4.0 * a_est)
    n = 250_000
    a = np.geomspace(a_min, a_max, n)

    h = h_of_a_sinv(a, h0_km_s_mpc, omega_r0, omega_m0, omega_k0, omega_lambda0)
    dt_da = 1.0 / np.maximum(a * h, 1e-300)

    # cumulative trapezoid integration
    da = np.diff(a)
    dt_mid = 0.5 * (dt_da[:-1] + dt_da[1:]) * da
    t = np.empty_like(a)
    t[0] = 0.0
    t[1:] = np.cumsum(dt_mid)

    # Ensure table covers target window.
    if t[-1] < t_max:
        # Expand range once more if required.
        a2 = np.geomspace(a_max, a_max * 20.0, n // 5)
        h2 = h_of_a_sinv(a2, h0_km_s_mpc, omega_r0, omega_m0, omega_k0, omega_lambda0)
        dt_da2 = 1.0 / np.maximum(a2 * h2, 1e-300)
        da2 = np.diff(a2)
        dt_mid2 = 0.5 * (dt_da2[:-1] + dt_da2[1:]) * da2
        t2 = np.empty_like(a2)
        t2[0] = t[-1]
        t2[1:] = t2[0] + np.cumsum(dt_mid2)
        a = np.concatenate([a, a2[1:]])
        t = np.concatenate([t, t2[1:]])

    return t, a


def plasma_image(size: int, temp_norm: float, phase: float) -> np.ndarray:
    x = np.linspace(-1.0, 1.0, size)
    y = np.linspace(-1.0, 1.0, size)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X * X + Y * Y)

    core = np.exp(-4.5 * R**2)
    shell = np.exp(-18.0 * (R - 0.55) ** 2)
    ripple = 0.5 + 0.5 * np.sin(14.0 * R - phase)
    swirl = 0.5 + 0.5 * np.sin(6.0 * X + 5.0 * Y + 0.7 * phase)

    img = core * (1.1 + 0.8 * temp_norm) + 0.35 * shell * ripple + 0.12 * swirl
    img *= (R <= 1.0)
    img = np.clip(img, 0.0, None)
    return img


def render_frame(frame_path: pathlib.Path, i: int, nframes: int, t_now: float, t_series: np.ndarray, a_series: np.ndarray, h0_km_s_mpc: float, omega_r0: float, omega_m0: float, omega_k0: float, omega_lambda0: float):
    u = 0.0 if nframes <= 1 else i / (nframes - 1)
    t_cur = t_now * u

    a_cur = float(np.interp(t_cur, t_series, a_series))
    z_cur = 1.0 / max(a_cur, 1e-300) - 1.0
    temp_cur = T_CMB0_K / max(a_cur, 1e-300)
    h_cur = float(h_of_a_sinv(np.array([a_cur]), h0_km_s_mpc, omega_r0, omega_m0, omega_k0, omega_lambda0)[0])
    h_cur_km_s_mpc = h_cur * METER_PER_MPC / 1000.0
    or_cur, om_cur, ol_cur = [float(x[0]) for x in omega_components_at_a(np.array([a_cur]), omega_r0, omega_m0, omega_k0, omega_lambda0)]

    # Precompute trajectories for right-side plots.
    t_plot = np.linspace(1e-6, t_now, 1200)
    a_plot = np.interp(t_plot, t_series, a_series)
    z_plot = 1.0 / np.maximum(a_plot, 1e-300) - 1.0
    temp_plot = T_CMB0_K / np.maximum(a_plot, 1e-300)
    or_plot, om_plot, ol_plot = omega_components_at_a(a_plot, omega_r0, omega_m0, omega_k0, omega_lambda0)

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    gs = GridSpec(2, 3, width_ratios=[1.5, 1.0, 1.0], height_ratios=[1.0, 1.0], figure=fig)

    # Left cinematic plasma panel.
    ax0 = fig.add_subplot(gs[:, 0])
    temp_norm = np.clip((math.log10(max(temp_cur, 1.0)) - 9.0) / 3.0, 0.0, 1.0)
    img = plasma_image(620, temp_norm=temp_norm, phase=10.0 * u)
    ax0.imshow(img, cmap="inferno", origin="lower")
    ax0.set_axis_off()
    ax0.set_title("Early-Universe Plasma Proxy (GUTOE FRW lane)", fontsize=15, pad=10)

    # Top-right: temperature and redshift trajectories.
    ax1 = fig.add_subplot(gs[0, 1:])
    ax1.plot(t_plot, temp_plot, color="#ff6b00", lw=2.0, label="T(t) [K]")
    ax1.scatter([max(t_cur, 1e-6)], [temp_cur], color="#ffd166", s=45, zorder=5)
    ax1.set_yscale("log")
    ax1.set_xlabel("time [s]")
    ax1.set_ylabel("temperature [K]")
    ax1.grid(alpha=0.25)
    ax1_t = ax1.twinx()
    ax1_t.plot(t_plot, z_plot, color="#4ea8de", lw=1.6, alpha=0.85, label="z(t)")
    ax1_t.scatter([max(t_cur, 1e-6)], [z_cur], color="#90e0ef", s=30, zorder=5)
    ax1_t.set_yscale("log")
    ax1_t.set_ylabel("redshift z")

    # Bottom-middle: density fractions.
    ax2 = fig.add_subplot(gs[1, 1])
    ax2.plot(t_plot, or_plot, color="#f94144", lw=2.0, label=r"$\Omega_r$")
    ax2.plot(t_plot, om_plot, color="#43aa8b", lw=2.0, label=r"$\Omega_m$")
    ax2.plot(t_plot, ol_plot, color="#577590", lw=2.0, label=r"$\Omega_\Lambda$")
    ax2.scatter([max(t_cur, 1e-6)], [or_cur], color="#f94144", s=28)
    ax2.scatter([max(t_cur, 1e-6)], [om_cur], color="#43aa8b", s=28)
    ax2.scatter([max(t_cur, 1e-6)], [ol_cur], color="#577590", s=28)
    ax2.set_xscale("log")
    ax2.set_xlabel("time [s]")
    ax2.set_ylabel("density fraction")
    ax2.set_ylim(-0.02, 1.02)
    ax2.grid(alpha=0.25)
    ax2.legend(loc="upper right", fontsize=10)

    # Bottom-right: telemetry text.
    ax3 = fig.add_subplot(gs[1, 2])
    ax3.axis("off")
    text = (
        f"First 180 Seconds\n\n"
        f"t = {t_cur:10.6f} s\n"
        f"z = {z_cur: .3e}\n"
        f"T = {temp_cur: .3e} K\n"
        f"H = {h_cur_km_s_mpc: .3e} km/s/Mpc\n\n"
        f"Omega_r = {or_cur: .6f}\n"
        f"Omega_m = {om_cur: .6f}\n"
        f"Omega_L = {ol_cur: .6f}\n\n"
        f"H0 = {h0_km_s_mpc:.4f} km/s/Mpc\n"
        f"Omega_b0 = {OMEGA_B0:.4f}\n"
        f"Omega_dm0 = {OMEGA_DM0:.4f}\n"
        f"Lambda = {LAMBDA_FULL:.3e} 1/m^2"
    )
    ax3.text(0.02, 0.98, text, va="top", ha="left", family="monospace", fontsize=12)

    fig.suptitle("GUTOE Universe Simulation — t≈0+ to 180 s", fontsize=18, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(frame_path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/tmp/bh_renders")
    ap.add_argument("--duration-s", type=float, default=180.0)
    ap.add_argument("--clip-seconds", type=float, default=30.0, help="movie duration in wall-clock seconds")
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--skip-gif", action="store_true")
    args = ap.parse_args()

    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    omega_lambda0 = 1.0 - OMEGA_M0 - OMEGA_R0 - OMEGA_K0
    h0_km_s_mpc = h0_from_lambda_and_omega_lambda(LAMBDA_FULL, omega_lambda0)

    t_series, a_series = build_time_to_scale_table(
        t_max=args.duration_s,
        h0_km_s_mpc=h0_km_s_mpc,
        omega_r0=OMEGA_R0,
        omega_m0=OMEGA_M0,
        omega_k0=OMEGA_K0,
        omega_lambda0=omega_lambda0,
    )

    nframes = max(2, int(round(args.clip_seconds * args.fps)))
    frame_dir = pathlib.Path(tempfile.mkdtemp(prefix="universe_180s_frames_"))

    try:
        for i in range(nframes):
            frame_path = frame_dir / f"frame_{i:05d}.png"
            render_frame(
                frame_path,
                i=i,
                nframes=nframes,
                t_now=args.duration_s,
                t_series=t_series,
                a_series=a_series,
                h0_km_s_mpc=h0_km_s_mpc,
                omega_r0=OMEGA_R0,
                omega_m0=OMEGA_M0,
                omega_k0=OMEGA_K0,
                omega_lambda0=omega_lambda0,
            )

        ffmpeg = shutil.which("ffmpeg")
        if ffmpeg is None:
            raise RuntimeError("ffmpeg not found in PATH")

        mp4 = out_dir / "universe_first_180s_1080p.mp4"
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

        gif = out_dir / "universe_first_180s_1080p.gif"
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
            "duration_sim_seconds": args.duration_s,
            "clip_seconds": args.clip_seconds,
            "fps": args.fps,
            "frames": nframes,
            "h0_km_s_mpc": h0_km_s_mpc,
            "omega_b0": OMEGA_B0,
            "omega_dm0": OMEGA_DM0,
            "omega_m0": OMEGA_M0,
            "omega_r0": OMEGA_R0,
            "omega_lambda0": omega_lambda0,
            "lambda_full": LAMBDA_FULL,
            "earliest_table_time_seconds": float(t_series[0]),
            "latest_table_time_seconds": float(t_series[-1]),
            "artifacts": {
                "mp4": str(mp4),
                "gif": None if args.skip_gif else str(gif),
            },
        }
        summary_path = out_dir / "universe_first_180s_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

        print("wrote", mp4)
        if not args.skip_gif:
            print("wrote", gif)
        print("wrote", summary_path)
    finally:
        shutil.rmtree(frame_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
