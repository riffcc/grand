#!/usr/bin/env python3
"""
Extract X-ray/UV lane from the 12h simulation, render a standalone movie,
and quantify rotational geometry over time.

Outputs (default /tmp/bh_renders/universe_12h):
- universe_12h_xray_only_1080p.mp4
- universe_12h_xray_only_1080p.gif
- universe_12h_xray_rotation_trace.png
- universe_12h_xray_rotation_summary.json
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

# Constants mirrored from physics lane
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


def xray_channel(field: np.ndarray) -> np.ndarray:
    return np.clip((field - 0.55) / 0.45, 0.0, 1.0) ** 1.1


def orientation_from_second_moment(img: np.ndarray) -> float:
    """Return principal orientation angle (radians), pi-periodic."""
    h, w = img.shape
    yy, xx = np.mgrid[0:h, 0:w]
    x = (xx - 0.5 * (w - 1)) / max(w - 1, 1)
    y = (yy - 0.5 * (h - 1)) / max(h - 1, 1)

    weights = np.clip(img, 0.0, None) ** 2.0
    s = float(weights.sum())
    if s <= 1e-12:
        return 0.0

    mx = float((weights * x).sum() / s)
    my = float((weights * y).sum() / s)
    dx = x - mx
    dy = y - my

    cxx = float((weights * dx * dx).sum() / s)
    cyy = float((weights * dy * dy).sum() / s)
    cxy = float((weights * dx * dy).sum() / s)

    return 0.5 * math.atan2(2.0 * cxy, cxx - cyy)


def best_mod90_residual_deg(theta_deg: np.ndarray) -> tuple[float, float]:
    """
    Fit an offset phi0 such that theta is closest to (phi0 + 90*k).
    Returns (phi0_deg, rms_residual_deg).
    """
    phi_grid = np.linspace(0.0, 90.0, 3601)  # 0.025° resolution
    best_phi = 0.0
    best_rms = float("inf")
    for phi in phi_grid:
        r = (theta_deg - phi + 45.0) % 90.0 - 45.0
        rms = float(np.sqrt(np.mean(r * r)))
        if rms < best_rms:
            best_rms = rms
            best_phi = float(phi)
    return best_phi, best_rms


def jump_stats(theta_deg: np.ndarray) -> dict:
    d = np.diff(theta_deg)
    absd = np.abs(d)
    med = float(np.median(absd))
    mad = float(np.median(np.abs(absd - med)))
    thr = med + 6.0 * max(mad, 1e-6)
    jump_idx = np.where(absd > thr)[0]
    jump_sizes = d[jump_idx]
    return {
        "median_abs_step_deg": med,
        "mad_abs_step_deg": mad,
        "jump_threshold_deg": float(thr),
        "jump_count": int(jump_idx.size),
        "jump_indices": jump_idx.astype(int).tolist(),
        "jump_sizes_deg": [float(x) for x in jump_sizes],
        "jump_median_abs_deg": float(np.median(np.abs(jump_sizes))) if jump_idx.size > 0 else 0.0,
    }


def render_xray_frame(frame_path: pathlib.Path, xr: np.ndarray, t_cur: float, z_cur: float, temp_cur: float, h_cur_km_s_mpc: float):
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    ax = fig.add_subplot(111)

    im = ax.imshow(xr, cmap="viridis", origin="lower")
    levels = np.linspace(0.35, 0.95, 7)
    ax.contour(xr, levels=levels, colors="white", linewidths=0.9, alpha=0.75)

    ax.set_title("X-ray / UV Proxy — First 12h Universe Lane", fontsize=18, pad=12)
    ax.set_axis_off()

    text = (
        f"t = {t_cur:.3e} s ({t_cur/SEC_PER_YEAR:.3e} yr)\n"
        f"z = {z_cur:.3e}\n"
        f"T = {temp_cur:.3e} K\n"
        f"H = {h_cur_km_s_mpc:.3e} km/s/Mpc"
    )
    ax.text(0.02, 0.98, text, transform=ax.transAxes, va="top", ha="left", fontsize=12,
            family="monospace", color="white", bbox=dict(facecolor="black", alpha=0.55, boxstyle="round,pad=0.35"))

    cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.02)
    cbar.set_label("X-ray proxy intensity")

    fig.tight_layout()
    fig.savefig(frame_path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/tmp/bh_renders/universe_12h")
    ap.add_argument("--duration-s", type=float, default=43_200.0)
    ap.add_argument("--clip-seconds", type=float, default=20.0)
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--size", type=int, default=520)
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
    frame_dir = pathlib.Path(tempfile.mkdtemp(prefix="universe_12h_xray_frames_"))

    orient = []
    sim_t = []

    try:
        for i in range(nframes):
            u = 0.0 if nframes <= 1 else i / (nframes - 1)
            logt = math.log10(t_start) + u * (math.log10(t_end) - math.log10(t_start))
            t_cur = 10.0**logt

            a_cur = float(np.interp(t_cur, t_table, a_table))
            z_cur = 1.0 / max(a_cur, 1e-300) - 1.0
            temp_cur = T_CMB0_K / max(a_cur, 1e-300)
            h_cur = float(h_of_a_sinv(np.array([a_cur]), h0_km_s_mpc, OMEGA_R0, OMEGA_M0, OMEGA_K0, omega_lambda0)[0])
            h_cur_km_s_mpc = h_cur * METER_PER_MPC / 1000.0

            temp_norm = np.clip((math.log10(max(temp_cur, 1.0)) - 8.0) / 5.0, 0.0, 1.0)
            field = plasma_field(size=args.size, phase=12.0 * u, temp_norm=temp_norm)
            xr = xray_channel(field)

            theta = orientation_from_second_moment(xr)
            orient.append(theta)
            sim_t.append(t_cur)

            render_xray_frame(
                frame_dir / f"frame_{i:05d}.png",
                xr,
                t_cur=t_cur,
                z_cur=z_cur,
                temp_cur=temp_cur,
                h_cur_km_s_mpc=h_cur_km_s_mpc,
            )

        orient = np.array(orient)
        sim_t = np.array(sim_t)

        # Unwrap pi-periodic orientation by doubling, unwrap, then halve.
        orient_unwrapped = 0.5 * np.unwrap(2.0 * orient)
        frame_idx = np.arange(nframes, dtype=float)

        # Linear drift in movie frame index.
        coeff = np.polyfit(frame_idx, orient_unwrapped, deg=1)
        slope = float(coeff[0])
        intercept = float(coeff[1])
        fit = slope * frame_idx + intercept

        # R^2
        ss_res = float(np.sum((orient_unwrapped - fit) ** 2))
        ss_tot = float(np.sum((orient_unwrapped - orient_unwrapped.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-15 else 1.0

        total_rotation_rad = float(orient_unwrapped[-1] - orient_unwrapped[0])
        total_rotation_deg = total_rotation_rad * 180.0 / math.pi
        turns = total_rotation_deg / 360.0

        # Quantization diagnostics
        theta_deg = orient_unwrapped * 180.0 / math.pi
        phi0_deg, mod90_rms_deg = best_mod90_residual_deg(theta_deg)
        jstats = jump_stats(theta_deg)

        # trace plot
        fig = plt.figure(figsize=(12, 5), dpi=140)
        ax = fig.add_subplot(111)
        ax.plot(frame_idx, theta_deg, label="orientation (deg)", color="#2a9d8f", lw=2)
        ax.plot(frame_idx, fit * 180.0 / math.pi, label="linear drift fit", color="#e76f51", lw=1.6, ls="--")
        ax.set_xlabel("frame")
        ax.set_ylabel("orientation [deg]")
        ax.set_title("X-ray geometric orientation trace")
        ax.grid(alpha=0.25)
        ax.legend(loc="best")
        trace_path = out_dir / "universe_12h_xray_rotation_trace.png"
        fig.tight_layout()
        fig.savefig(trace_path)
        plt.close(fig)

        # Save raw orientation trace for independent checks.
        csv_path = out_dir / "universe_12h_xray_orientation_trace.csv"
        with csv_path.open("w", encoding="utf-8") as f:
            f.write("frame,time_seconds,orientation_deg,fit_deg,residual_deg\n")
            fit_deg = fit * 180.0 / math.pi
            residual = theta_deg - fit_deg
            for fi in range(nframes):
                f.write(
                    f"{fi},{sim_t[fi]:.12e},{theta_deg[fi]:.9f},{fit_deg[fi]:.9f},{residual[fi]:.9f}\n"
                )

        ffmpeg = shutil.which("ffmpeg")
        if ffmpeg is None:
            raise RuntimeError("ffmpeg not found in PATH")

        mp4 = out_dir / "universe_12h_xray_only_1080p.mp4"
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

        gif = out_dir / "universe_12h_xray_only_1080p.gif"
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
            "rotation": {
                "total_rotation_rad": total_rotation_rad,
                "total_rotation_deg": total_rotation_deg,
                "turns": turns,
                "linear_slope_rad_per_frame": slope,
                "linear_slope_deg_per_frame": slope * 180.0 / math.pi,
                "linear_r2": r2,
                "mod90_best_offset_deg": phi0_deg,
                "mod90_rms_residual_deg": mod90_rms_deg,
                "jump_stats": jstats,
            },
            "artifacts": {
                "mp4": str(mp4),
                "gif": None if args.skip_gif else str(gif),
                "trace": str(trace_path),
                "trace_csv": str(csv_path),
            },
        }

        summary_path = out_dir / "universe_12h_xray_rotation_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

        print("wrote", mp4)
        if not args.skip_gif:
            print("wrote", gif)
        print("wrote", trace_path)
        print("wrote", csv_path)
        print("wrote", summary_path)
        print(
            f"rotation: total={total_rotation_deg:.3f} deg ({turns:.4f} turns), slope={slope*180.0/math.pi:.5f} deg/frame, R2={r2:.4f}"
        )
        print(
            f"quantization: mod90_rms={mod90_rms_deg:.3f} deg @ offset={phi0_deg:.3f} deg, jumps={jstats['jump_count']} (median |jump|={jstats['jump_median_abs_deg']:.3f} deg)"
        )
    finally:
        shutil.rmtree(frame_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
