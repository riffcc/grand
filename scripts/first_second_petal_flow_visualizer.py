#!/usr/bin/env python3
"""
Color-code petal sectors and overlay rotational flow from first-second sampled frames.

Inputs:
- frames from tile export (default: /tmp/bh_renders/first_second/ai_tiles/frames)

Outputs:
- /tmp/bh_renders/first_second/petal_flow/petal_flow_colored.mp4
- /tmp/bh_renders/first_second/petal_flow/petal_flow_colored.gif
- /tmp/bh_renders/first_second/petal_flow/petal_flow_summary.json
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
from PIL import Image
from scipy import ndimage

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def principal_orientation_deg(img: np.ndarray) -> float:
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
    theta = 0.5 * math.atan2(2.0 * cxy, cxx - cyy)
    return theta * 180.0 / math.pi


def build_colored_overlay(gray: np.ndarray, petals: int = 9):
    h, w = gray.shape
    yy, xx = np.mgrid[0:h, 0:w]
    cx, cy = (w - 1) / 2.0, (h - 1) / 2.0
    X = xx - cx
    Y = yy - cy
    R = np.sqrt(X * X + Y * Y)
    A = (np.arctan2(Y, X) + 2 * np.pi) % (2 * np.pi)

    rmin = 0.10 * min(h, w)
    rmax = 0.33 * min(h, w)
    ann = (R >= rmin) & (R <= rmax)

    g = gray.astype(float)
    g = (g - np.percentile(g, 5)) / max(np.percentile(g, 99) - np.percentile(g, 5), 1e-9)
    g = np.clip(g, 0.0, 1.0)

    # Detect petal blobs from bright annulus structure, then color per blob
    # (not by fixed angular sector) so each physical petal carries one color.
    ann_vals = g[ann]
    thr = np.percentile(ann_vals, 82.0) if ann_vals.size > 0 else 0.8
    binary = ann & (g >= thr)
    binary = ndimage.binary_opening(binary, structure=np.ones((3, 3)))
    binary = ndimage.binary_closing(binary, structure=np.ones((3, 3)))

    labels, nlab = ndimage.label(binary)
    comp = []
    for lab in range(1, nlab + 1):
        m = labels == lab
        area = int(m.sum())
        if area < 18:
            continue
        wmean = float(g[m].mean())
        # Use angular centroid for stable ordering around the ring.
        theta = float(np.angle(np.mean(np.exp(1j * A[m]))))
        if theta < 0:
            theta += 2 * np.pi
        score = area * wmean
        comp.append((score, theta, m))

    # Keep strongest components and sort by angle.
    comp.sort(key=lambda x: x[0], reverse=True)
    comp = comp[: max(1, petals)]
    comp.sort(key=lambda x: x[1])

    weight = np.zeros_like(g)
    color_img = np.zeros((h, w, 3), dtype=float)
    cmap = plt.get_cmap("tab10")
    for i, (_, _, m) in enumerate(comp):
        c = np.array(cmap((i % 10) / 9.0)[:3])
        local_w = np.clip((g[m] - thr) / max(1e-6, 1.0 - thr), 0.0, 1.0)
        weight[m] = np.maximum(weight[m], local_w)
        color_img[m] = c

    base = np.dstack([g, g, g])
    out = base * (1.0 - 0.78 * weight[..., None]) + color_img * (0.78 * weight[..., None])

    # soft ring guide
    ring = np.exp(-((R - 0.21 * min(h, w)) / (0.018 * min(h, w))) ** 2)
    out[..., 1] += 0.08 * ring
    out[..., 2] += 0.05 * ring
    out = np.clip(out, 0.0, 1.0)

    # rotational flow field (tangent vectors)
    mag = np.clip(weight * 1.2, 0.0, 1.0)
    nx = np.zeros_like(X, dtype=float)
    ny = np.zeros_like(Y, dtype=float)
    rr = np.maximum(R, 1e-9)
    nx = -Y / rr
    ny = X / rr

    return out, ann, nx * mag, ny * mag


def render_frame(src_path: pathlib.Path, dst_path: pathlib.Path, frame_i: int, nframes: int, petals: int):
    rgb = np.array(Image.open(src_path).convert("RGB"), dtype=float) / 255.0

    h, w, _ = rgb.shape
    # Crop central square where morphology lives (avoid HUD/text-heavy borders)
    s = int(min(h, w) * 0.74)
    y0 = (h - s) // 2
    x0 = (w - s) // 2
    crop = rgb[y0 : y0 + s, x0 : x0 + s]
    gray = 0.2126 * crop[..., 0] + 0.7152 * crop[..., 1] + 0.0722 * crop[..., 2]

    overlay, ann, u, v = build_colored_overlay(gray, petals=petals)
    theta = principal_orientation_deg(gray)

    fig = plt.figure(figsize=(12.8, 7.2), dpi=150)
    ax = fig.add_subplot(111)
    ax.imshow(overlay, origin="lower")

    # quiver flow on sparse grid inside annulus
    step = max(10, s // 48)
    yy, xx = np.mgrid[0:s:step, 0:s:step]
    ann_s = ann[::step, ::step]
    uu = u[::step, ::step]
    vv = v[::step, ::step]

    # mask vectors outside annulus
    uu = np.where(ann_s, uu, np.nan)
    vv = np.where(ann_s, vv, np.nan)

    ax.quiver(xx, yy, uu, vv, color="#00e5ff", alpha=0.70, scale=26, width=0.003)

    ax.set_axis_off()
    ax.set_title(f"Petal Flow Map ({petals} sectors) — frame {frame_i+1}/{nframes}", fontsize=14)
    ax.text(
        0.015,
        0.96,
        f"orientation ~ {theta:.2f} deg",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=10,
        color="white",
        family="monospace",
        bbox=dict(facecolor="black", alpha=0.55, boxstyle="round,pad=0.35"),
    )

    fig.tight_layout()
    fig.savefig(dst_path)
    plt.close(fig)
    return theta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames-dir", default="/tmp/bh_renders/first_second/ai_tiles/frames")
    ap.add_argument("--out-dir", default="/tmp/bh_renders/first_second/petal_flow")
    ap.add_argument("--petals", type=int, default=9)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--skip-gif", action="store_true")
    args = ap.parse_args()

    src_dir = pathlib.Path(args.frames_dir)
    frames = sorted(src_dir.glob("frame_*.png"))
    if not frames:
        raise FileNotFoundError(f"no frames found in {src_dir}")

    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    render_dir = pathlib.Path(tempfile.mkdtemp(prefix="petal_flow_frames_"))
    thetas = []

    try:
        for i, f in enumerate(frames):
            t = render_dir / f"frame_{i:05d}.png"
            theta = render_frame(f, t, i, len(frames), args.petals)
            thetas.append(theta)

        ffmpeg = shutil.which("ffmpeg")
        if ffmpeg is None:
            raise RuntimeError("ffmpeg not found in PATH")

        mp4 = out_dir / "petal_flow_colored.mp4"
        subprocess.run(
            [
                ffmpeg,
                "-y",
                "-framerate",
                str(args.fps),
                "-i",
                str(render_dir / "frame_%05d.png"),
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

        gif = out_dir / "petal_flow_colored.gif"
        if not args.skip_gif:
            palette = render_dir / "palette.png"
            subprocess.run(
                [ffmpeg, "-y", "-i", str(render_dir / "frame_%05d.png"), "-vf", "palettegen", str(palette)],
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
                    str(render_dir / "frame_%05d.png"),
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
            "source_frames_dir": str(src_dir),
            "frames": len(frames),
            "petals": args.petals,
            "fps": args.fps,
            "orientation_deg": {
                "min": float(np.min(thetas)),
                "max": float(np.max(thetas)),
                "mean": float(np.mean(thetas)),
            },
            "artifacts": {
                "mp4": str(mp4),
                "gif": None if args.skip_gif else str(gif),
            },
        }
        summary_path = out_dir / "petal_flow_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

        print("wrote", mp4)
        if not args.skip_gif:
            print("wrote", gif)
        print("wrote", summary_path)
    finally:
        shutil.rmtree(render_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
