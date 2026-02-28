#!/usr/bin/env python3
"""
Render first-second universe as a multispectral 1080p movie (no petal recolor overlays).

Outputs (default /tmp/bh_renders/first_second_multispectral):
- universe_first_second_multispectral_1080p.mp4
- universe_first_second_multispectral_1080p.gif
- universe_first_second_multispectral_summary.json
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
from scipy import ndimage

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PLANCK_TIME = 5.391247e-44
T_END = 1.0

# Milestones (seconds)
T_INFLATION = 1e-36
T_EW_BREAK = 1e-12
T_QCD = 1e-6
T_NEUTRINO = 1e-2


def smoothstep(x: float, c: float, w: float) -> float:
    return 1.0 / (1.0 + math.exp(-(x - c) / max(w, 1e-9)))


def temp_radiation_era_k(t: float) -> float:
    # Approx: T ~ 1.16e10 K / sqrt(t/s)
    return 1.16045e10 / math.sqrt(max(t, 1e-40))


def phase_fractions(t: float):
    logt = math.log10(max(t, PLANCK_TIME))
    s_infl = smoothstep(logt, math.log10(T_INFLATION), 0.35)
    s_ew = smoothstep(logt, math.log10(T_EW_BREAK), 0.25)
    s_qcd = smoothstep(logt, math.log10(T_QCD), 0.30)

    foam = (1.0 - s_infl)
    inflation = s_infl * (1.0 - s_ew)
    plasma = s_ew * (1.0 - s_qcd)
    hadronic = s_qcd

    norm = foam + inflation + plasma + hadronic
    if norm <= 0.0:
        return 1.0, 0.0, 0.0, 0.0
    return foam / norm, inflation / norm, plasma / norm, hadronic / norm


def make_grid(size: int):
    x = np.linspace(-1.0, 1.0, size)
    y = np.linspace(-1.0, 1.0, size)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X * X + Y * Y)
    A = np.arctan2(Y, X)
    return X, Y, R, A


def foam_layer(X: np.ndarray, Y: np.ndarray, phase: float):
    z = (
        np.sin(22.0 * X + 17.0 * Y + phase)
        + np.sin(31.0 * X - 12.0 * Y + 1.7 * phase)
        + np.sin(19.0 * (X + Y) - 0.9 * phase)
    )
    z = (z - z.min()) / max(z.max() - z.min(), 1e-9)
    return z


def inflation_layer(R: np.ndarray, A: np.ndarray, phase: float, inflation_progress: float):
    ring = np.exp(-24.0 * (R - (0.10 + 0.55 * inflation_progress)) ** 2)
    rays = 0.5 + 0.5 * np.sin(9.0 * A + 3.0 * phase)
    core = np.exp(-20.0 * R * R)
    z = 0.7 * ring * rays + 0.6 * core
    z = (z - z.min()) / max(z.max() - z.min(), 1e-9)
    return z


def plasma_layer(X: np.ndarray, Y: np.ndarray, R: np.ndarray, phase: float, temp_norm: float):
    core = np.exp(-4.0 * R * R)
    turb = (
        0.45 * np.sin(8.0 * X + 6.5 * Y + phase)
        + 0.35 * np.sin(14.0 * R - 1.4 * phase)
        + 0.20 * np.sin(11.0 * X - 13.0 * Y + 0.7 * phase)
    )
    shell = np.exp(-16.0 * (R - 0.58) ** 2)
    z = core * (1.1 + 1.4 * temp_norm) + 0.20 * shell + 0.18 * turb
    z *= (R <= 1.0)
    z = (z - z.min()) / max(z.max() - z.min(), 1e-9)
    return z


def spectral_maps(base: np.ndarray):
    # Derived pseudo-spectral channels from the same state field
    radio = ndimage.gaussian_filter(base, sigma=6.0)
    microwave = np.sqrt(np.clip(base, 0.0, 1.0))
    infrared = np.clip(base ** 1.35, 0.0, 1.0)
    visible = np.clip(base, 0.0, 1.0)
    uv = np.clip((base - 0.45) / 0.55, 0.0, 1.0)
    xray = np.clip((base - 0.62) / 0.38, 0.0, 1.0) ** 1.15
    gamma = np.clip((base - 0.80) / 0.20, 0.0, 1.0) ** 1.35

    gy, gx = np.gradient(base)
    bmag = np.sqrt(gx * gx + gy * gy)
    bmag = bmag / max(float(bmag.max()), 1e-9)

    entropy = -np.clip(base, 1e-9, 1.0) * np.log(np.clip(base, 1e-9, 1.0))
    entropy = entropy / max(float(entropy.max()), 1e-9)

    return {
        "Radio": (radio, "cividis"),
        "Microwave": (microwave, "coolwarm"),
        "Infrared": (infrared, "magma"),
        "Visible": (visible, "inferno"),
        "Ultraviolet": (uv, "plasma"),
        "X-ray": (xray, "viridis"),
        "Gamma": (gamma, "turbo"),
        "Magnetic |B|": (bmag, "cubehelix"),
        "Entropy Proxy": (entropy, "Spectral_r"),
    }


def stage_name(t: float) -> str:
    if t < T_INFLATION:
        return "Quantum foam"
    if t < T_EW_BREAK:
        return "Inflation / expansion"
    if t < T_QCD:
        return "Electroweak broken phase"
    if t < T_NEUTRINO:
        return "QCD confinement onset"
    return "Approaching neutrino decoupling"


def render_frame(frame_path: pathlib.Path, i: int, nframes: int, X, Y, R, A):
    u = 0.0 if nframes <= 1 else i / (nframes - 1)
    logt = math.log10(PLANCK_TIME) + u * (math.log10(T_END) - math.log10(PLANCK_TIME))
    t = 10.0**logt

    temp_k = temp_radiation_era_k(t)
    temp_norm = np.clip((math.log10(max(temp_k, 1.0)) - 10.0) / 18.0, 0.0, 1.0)
    foam_w, infl_w, plasma_w, had_w = phase_fractions(t)

    infl_prog = np.clip(
        (logt - math.log10(PLANCK_TIME)) / (math.log10(T_EW_BREAK) - math.log10(PLANCK_TIME)),
        0.0,
        1.0,
    )

    phase = 16.0 * u
    foam = foam_layer(X, Y, phase)
    infl = inflation_layer(R, A, phase, infl_prog)
    plasma = plasma_layer(X, Y, R, phase, temp_norm)

    base = (
        0.95 * foam_w * foam
        + 0.95 * infl_w * infl
        + 1.05 * plasma_w * plasma
        + 0.30 * had_w * ndimage.gaussian_filter(plasma, sigma=2.2)
    )
    base = np.clip(base, 0.0, None)
    base = (base - base.min()) / max(float(base.max() - base.min()), 1e-9)

    maps = spectral_maps(base)

    fig, axes = plt.subplots(3, 3, figsize=(19.2, 10.8), dpi=100)
    for ax, (name, (img, cmap)) in zip(axes.ravel(), maps.items()):
        ax.imshow(img, cmap=cmap, origin="lower")
        ax.set_title(name, fontsize=11)
        ax.set_axis_off()

    stage = stage_name(t)
    telemetry = (
        f"t = {t:.3e} s | T = {temp_k:.3e} K | stage = {stage}\n"
        f"weights [foam,infl,plasma,had] = [{foam_w:.2f}, {infl_w:.2f}, {plasma_w:.2f}, {had_w:.2f}]"
    )

    fig.suptitle("First Second of Universe — Multispectral View", fontsize=18, y=0.985)
    fig.text(
        0.5,
        0.02,
        telemetry,
        ha="center",
        va="bottom",
        fontsize=11,
        family="monospace",
        color="white",
        bbox=dict(facecolor="black", alpha=0.60, boxstyle="round,pad=0.35"),
    )

    fig.tight_layout(rect=[0, 0.06, 1, 0.96])
    fig.savefig(frame_path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/tmp/bh_renders/first_second_multispectral")
    ap.add_argument("--clip-seconds", type=float, default=120.0)
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--size", type=int, default=420)
    ap.add_argument("--skip-gif", action="store_true")
    args = ap.parse_args()

    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    nframes = max(2, int(round(args.clip_seconds * args.fps)))
    X, Y, R, A = make_grid(args.size)

    frame_dir = pathlib.Path(tempfile.mkdtemp(prefix="first_second_multi_frames_"))
    try:
        for i in range(nframes):
            render_frame(frame_dir / f"frame_{i:05d}.png", i, nframes, X, Y, R, A)

        ffmpeg = shutil.which("ffmpeg")
        if ffmpeg is None:
            raise RuntimeError("ffmpeg not found in PATH")

        mp4 = out_dir / "universe_first_second_multispectral_1080p.mp4"
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

        gif = out_dir / "universe_first_second_multispectral_1080p.gif"
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
            "sim_start_s": PLANCK_TIME,
            "sim_end_s": 1.0,
            "clip_seconds": args.clip_seconds,
            "fps": args.fps,
            "frames": nframes,
            "spectra": [
                "Radio",
                "Microwave",
                "Infrared",
                "Visible",
                "Ultraviolet",
                "X-ray",
                "Gamma",
                "Magnetic |B|",
                "Entropy Proxy",
            ],
            "artifacts": {
                "mp4": str(mp4),
                "gif": None if args.skip_gif else str(gif),
            },
        }
        summary_path = out_dir / "universe_first_second_multispectral_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

        print("wrote", mp4)
        if not args.skip_gif:
            print("wrote", gif)
        print("wrote", summary_path)
    finally:
        shutil.rmtree(frame_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
