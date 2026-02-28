#!/usr/bin/env python3
"""
Grid invariance falsification gate for multispectral renderer.

Runs square vs hex renders with identical physics params, extracts one panel
from frame 0, computes dominant angular mode + orientation, and fails when
drift exceeds thresholds.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import shutil
import subprocess
import tempfile

import numpy as np
from PIL import Image


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def extract_panel_png(mp4: pathlib.Path, out_png: pathlib.Path, panel_index: int) -> None:
    # 3x3 layout panel mapping row-major.
    tx = panel_index % 3
    ty = panel_index // 3
    # Renderer currently writes 1920x1080 by default unless overridden.
    w, h = 1920, 1080
    tw, th = w // 3, h // 3
    x, y = tx * tw, ty * th
    vf = f"select=eq(n\\,0),crop={tw}:{th}:{x}:{y}"
    run(
        [
            "ffmpeg",
            "-y",
            "-i",
            str(mp4),
            "-vf",
            vf,
            "-frames:v",
            "1",
            str(out_png),
        ]
    )


def angular_mode_metrics(img_path: pathlib.Path, max_mode: int = 24) -> tuple[int, float, float]:
    img = np.array(Image.open(img_path).convert("L"), dtype=np.float64) / 255.0
    h, w = img.shape
    cy, cx = h / 2.0, w / 2.0
    y, x = np.indices((h, w))
    xr = (x - cx) / (w / 2.0)
    yr = (y - cy) / (h / 2.0)
    r = np.sqrt(xr * xr + yr * yr)
    th = np.arctan2(yr, xr)
    mask = (r > 0.25) & (r < 0.7)
    vals = img[mask] - img[mask].mean()
    ang = th[mask]

    coeffs: list[complex] = []
    for m in range(1, max_mode + 1):
        coeffs.append((vals * np.exp(-1j * m * ang)).sum())
    powers = [abs(c) for c in coeffs]
    mode = int(np.argmax(powers)) + 1
    c = coeffs[mode - 1]
    orient_deg = (np.angle(c) / mode) * 180.0 / math.pi
    return mode, orient_deg, float(abs(c))


def wrap_angle_delta_deg(a: float, b: float) -> float:
    d = (a - b) % 360.0
    if d > 180.0:
        d -= 360.0
    return abs(d)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--renderer-bin",
        default="./target/release/universe_first_second_multispectral_movie",
    )
    ap.add_argument("--out-dir", default="/tmp/bh_renders/grid_invariance")
    ap.add_argument("--work-dir", default="/tmp/bh_renders/grid_invariance_work")
    ap.add_argument("--clip-seconds", type=float, default=1.0)
    ap.add_argument("--fps", type=int, default=1)
    ap.add_argument("--first-seconds", type=float, default=1e-43)
    ap.add_argument("--size", type=int, default=420)
    ap.add_argument("--panel-index", type=int, default=3, help="0..8, default Visible panel")
    ap.add_argument("--max-mode", type=int, default=24)
    ap.add_argument("--orientation-threshold-deg", type=float, default=5.0)
    ap.add_argument("--mode-must-match", action="store_true")
    args = ap.parse_args()

    out_root = pathlib.Path(args.out_dir)
    work_root = pathlib.Path(args.work_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)

    results: dict[str, dict[str, float | int | str]] = {}
    with tempfile.TemporaryDirectory(prefix="grid_inv_") as td:
        td_path = pathlib.Path(td)
        for grid in ("square", "hex"):
            od = out_root / grid
            od.mkdir(parents=True, exist_ok=True)
            cmd = [
                args.renderer_bin,
                "--out-dir",
                str(od),
                "--work-dir",
                str(work_root),
                "--clip-seconds",
                str(args.clip_seconds),
                "--fps",
                str(args.fps),
                "--first-seconds",
                str(args.first_seconds),
                "--size",
                str(args.size),
                "--grid",
                grid,
                "--skip-gif",
            ]
            run(cmd)
            mp4 = od / "universe_first_second_multispectral_1080p.mp4"
            png = td_path / f"{grid}_panel.png"
            extract_panel_png(mp4, png, args.panel_index)
            mode, orient_deg, power = angular_mode_metrics(png, args.max_mode)
            results[grid] = {
                "mode": mode,
                "orientation_deg": orient_deg,
                "power": power,
                "mp4": str(mp4),
            }

    square = results["square"]
    hex_ = results["hex"]
    mode_match = int(square["mode"]) == int(hex_["mode"])
    orient_delta = wrap_angle_delta_deg(
        float(square["orientation_deg"]), float(hex_["orientation_deg"])
    )
    pass_orientation = orient_delta <= args.orientation_threshold_deg
    passes = pass_orientation and (mode_match if args.mode_must_match else True)

    report = {
        "panel_index": args.panel_index,
        "params": {
            "clip_seconds": args.clip_seconds,
            "fps": args.fps,
            "first_seconds": args.first_seconds,
            "size": args.size,
            "max_mode": args.max_mode,
            "orientation_threshold_deg": args.orientation_threshold_deg,
            "mode_must_match": args.mode_must_match,
        },
        "square": square,
        "hex": hex_,
        "mode_match": mode_match,
        "orientation_delta_deg": orient_delta,
        "pass_orientation": pass_orientation,
        "passes": passes,
    }
    report_path = out_root / "grid_invariance_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {report_path}")
    print(
        f"mode square={square['mode']} hex={hex_['mode']} | "
        f"orientation_delta={orient_delta:.3f} deg | passes={passes}"
    )
    return 0 if passes else 2


if __name__ == "__main__":
    raise SystemExit(main())

