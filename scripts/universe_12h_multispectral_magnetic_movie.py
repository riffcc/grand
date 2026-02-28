#!/usr/bin/env python3
"""
Render first-12h universe movie in multiple spectral views + magnetic overlay.

Outputs (default /tmp/bh_renders/universe_12h):
- universe_12h_multispectral_magnetic_1080p.mp4
- universe_12h_multispectral_magnetic_1080p.gif
- universe_12h_multispectral_magnetic_summary.json
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

# Constants mirrored from the physics lane.
C = 299_792_458.0
METER_PER_MPC = 3.085_677_581_491_367e22
SEC_PER_YEAR = 31_557_600.0
PLANCK_TIME = 5.391247e-44
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


def build_t_a_table(h0_km_s_mpc: float, omega_r0: float, omega_m0: float, omega_k0: float, omega_lambda0: float):
    a = np.geomspace(1e-22, 1e-6, 300_000)
    h = h_of_a_sinv(a, h0_km_s_mpc, omega_r0, omega_m0, omega_k0, omega_lambda0)
    dt_da = 1.0 / np.maximum(a * h, 1e-300)
    da = np.diff(a)
    dt_mid = 0.5 * (dt_da[:-1] + dt_da[1:]) * da
    t = np.empty_like(a)
    t[0] = 0.0
    t[1:] = np.cumsum(dt_mid)
    return t, a


def plasma_field(size: int, phase: float, temp_norm: float):
    x = np.linspace(-1.0, 1.0, size)
    y = np.linspace(-1.0, 1.0, size)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X * X + Y * Y)

    core = np.exp(-3.8 * R**2)
    turbulence = (
        0.45 * np.sin(8.0 * X + 6.0 * Y + phase)
        + 0.35 * np.sin(14.0 * R - 1.3 * phase)
        + 0.20 * np.sin(17.0 * X - 11.0 * Y + 0.5 * phase)
    )
    shell = np.exp(-14.0 * (R - 0.62) ** 2)

    field = core * (1.2 + 1.4 * temp_norm) + 0.18 * shell + 0.16 * turbulence
    field *= (R <= 1.0)
    field -= field.min()
    field /= max(field.max(), 1e-9)
    return field


def spectrum_maps(field: np.ndarray):
    # Microwave: broad smooth anisotropy
    mw = np.sqrt(np.clip(field, 0.0, 1.0))
    # Infrared-ish: warm dense structures
    ir = np.clip(field ** 1.4, 0.0, 1.0)
    # X-ray/UV-ish: highlight hottest peaks
    xr = np.clip((field - 0.55) / 0.45, 0.0, 1.0) ** 1.1
    return mw, ir, xr


def magnetic_vectors(field: np.ndarray):
    gy, gx = np.gradient(field)
    # 2D magnetic proxy from perpendicular gradient
    bx = -gy
    by = gx
    mag = np.sqrt(bx * bx + by * by)
    m = max(mag.max(), 1e-9)
    bx /= m
    by /= m
    return bx, by, mag / m


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
):
    u = 0.0 if nframes <= 1 else i / (nframes - 1)
    t_cur = 10.0 ** (math.log10(t_start) + u * (math.log10(t_end) - math.log10(t_start)))

    a_cur = float(np.interp(t_cur, t_table, a_table))
    z_cur = 1.0 / max(a_cur, 1e-300) - 1.0
    temp_cur = T_CMB0_K / max(a_cur, 1e-300)
    h_cur = float(h_of_a_sinv(np.array([a_cur]), h0_km_s_mpc, omega_r0, omega_m0, omega_k0, omega_lambda0)[0])
    h_cur_km_s_mpc = h_cur * METER_PER_MPC / 1000.0

    temp_norm = np.clip((math.log10(max(temp_cur, 1.0)) - 8.0) / 5.0, 0.0, 1.0)
    field = plasma_field(size=480, phase=12.0 * u, temp_norm=temp_norm)
    mw, ir, xr = spectrum_maps(field)
    bx, by, bmag = magnetic_vectors(field)

    fig, axes = plt.subplots(2, 2, figsize=(19.2, 10.8), dpi=100)

    panels = [
        (axes[0, 0], field, "inferno", "Plasma / Visible Proxy"),
        (axes[0, 1], mw, "coolwarm", "Microwave Proxy"),
        (axes[1, 0], ir, "magma", "Infrared Proxy"),
        (axes[1, 1], xr, "viridis", "X-ray / UV Proxy"),
    ]

    for ax, img, cmap, title in panels:
        ax.imshow(img, cmap=cmap, origin="lower")
        ax.set_title(title, fontsize=12)
        ax.set_axis_off()

    # Magnetic overlay on X-ray panel
    axm = axes[1, 1]
    step = 18
    y = np.arange(0, field.shape[0], step)
    x = np.arange(0, field.shape[1], step)
    X, Y = np.meshgrid(x, y)
    U = bx[::step, ::step]
    V = by[::step, ::step]
    Cc = bmag[::step, ::step]
    axm.quiver(X, Y, U, V, Cc, cmap="cividis", alpha=0.75, scale=28, width=0.0025)

    fig.suptitle("GUTOE First 12 Hours — Multispectral + Magnetic View", fontsize=18, y=0.98)

    text = (
        f"t = {t_cur:.3e} s ({t_cur/SEC_PER_YEAR:.3e} yr)\n"
        f"z = {z_cur:.3e}\n"
        f"T = {temp_cur:.3e} K\n"
        f"H = {h_cur_km_s_mpc:.3e} km/s/Mpc\n"
        f"Omega_b0={OMEGA_B0:.4f}, Omega_dm0={OMEGA_DM0:.4f}\n"
        f"Lambda={LAMBDA_FULL:.3e} 1/m^2"
    )
    fig.text(0.012, 0.03, text, family="monospace", fontsize=11, color="white", bbox=dict(facecolor="black", alpha=0.6, boxstyle="round,pad=0.4"))

    fig.tight_layout(rect=[0, 0.06, 1, 0.95])
    fig.savefig(frame_path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/tmp/bh_renders/universe_12h")
    ap.add_argument("--duration-s", type=float, default=43_200.0)
    ap.add_argument("--clip-seconds", type=float, default=20.0)
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--skip-gif", action="store_true")
    args = ap.parse_args()

    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    omega_lambda0 = 1.0 - OMEGA_M0 - OMEGA_R0 - OMEGA_K0
    h0_km_s_mpc = h0_from_lambda_and_omega_lambda(LAMBDA_FULL, omega_lambda0)

    t_table, a_table = build_t_a_table(h0_km_s_mpc, OMEGA_R0, OMEGA_M0, OMEGA_K0, omega_lambda0)
    t_start = max(PLANCK_TIME, 1e-36)
    t_end = min(args.duration_s, float(t_table[-1]))

    nframes = max(2, int(round(args.clip_seconds * args.fps)))
    frame_dir = pathlib.Path(tempfile.mkdtemp(prefix="universe_12h_multi_frames_"))

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
            )

        ffmpeg = shutil.which("ffmpeg")
        if ffmpeg is None:
            raise RuntimeError("ffmpeg not found in PATH")

        mp4 = out_dir / "universe_12h_multispectral_magnetic_1080p.mp4"
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

        gif = out_dir / "universe_12h_multispectral_magnetic_1080p.gif"
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
            "duration_sim_seconds": t_end,
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
            "artifacts": {
                "mp4": str(mp4),
                "gif": None if args.skip_gif else str(gif),
            },
        }
        summary_path = out_dir / "universe_12h_multispectral_magnetic_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

        print("wrote", mp4)
        if not args.skip_gif:
            print("wrote", gif)
        print("wrote", summary_path)
    finally:
        shutil.rmtree(frame_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
